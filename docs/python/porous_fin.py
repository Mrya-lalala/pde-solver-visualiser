"""Hermite solver for Eq. (13)-(14) of Pavan Kumar et al. (2026).

  (A theta_X)_X - sqrt(A) F(theta) = A theta_tau, A=1-X**2,
  F(theta)=Nc*theta**2+Nr*theta+n1_squared*theta**(p+1)+wet_term*theta**p.

wet_term means n0**2*n3, not an independently inferred material property.
The initial interior value is zero and the base steps to one at tau=0+.
Uses piecewise shifted physicists' Hermite polynomials with C1 interface
constraints, base Dirichlet and tip Neumann constraints. A weighted Galerkin
form avoids dividing by zero at X=1. Time integration is implicit BDF.
An independent conservative finite-volume calculation is an optional numerical
comparison, never an exact solution or an input to the Hermite solve.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import cholesky, null_space, solve_triangular
from scipy.sparse import diags

from rosenau_hyman import hermite_integrals

DEFAULTS = dict(nc=0.5, nr=0.1, n1_squared=0.5, wet_term=0.1,
                power=1, duration=0.8, cells=12, count=6, reference=True)
PARAMETER_NOTE = ('Illustrative dimensionless coefficients for Eq. (13); '
                  'not a reconstruction of either ceramic material baseline.')
DATA_TOLERANCE = 1e-10
COMPARISON_TOLERANCE = 2e-3
WEAK_RESIDUAL_TOLERANCE = 1e-5
STRONG_RESIDUAL_TOLERANCE = 2e-2


def validate_parameters(p):
    if not isinstance(p, dict) or set(p) != set(DEFAULTS):
        raise ValueError('Provide exactly the documented porous-fin numeric parameters.')
    if type(p['reference']) is not bool:
        raise ValueError('reference must be a boolean (independent numerical comparison).')
    for key in set(p)-{'reference'}:
        if type(p[key]) not in (int, float) or not np.isfinite(p[key]):
            raise ValueError(f'{key} must be a finite number.')
    for key in ('power', 'cells', 'count'):
        if type(p[key]) is not int:
            raise ValueError(f'{key} must be an integer.')
    if p['power'] not in (1, 2) or p['count'] not in (4, 6) or p['cells'] not in (8, 12, 16):
        raise ValueError('Use p=1 or 2, N=4 or 6, and 8, 12 or 16 cells.')
    if not 0.2 <= p['duration'] <= 1:
        raise ValueError('Time horizon must be between 0.2 and 1.')
    if any(not 0 <= p[key] <= 2 for key in ('nc', 'nr', 'n1_squared')) or not -0.5 <= p['wet_term'] <= 2:
        raise ValueError('Use Nc, Nr and n1^2 in [0,2], and n0^2*n3 in [-0.5,2].')
    return dict(p)


def reaction(u, p):
    return (p['nc']*u**2 + p['nr']*u + p['n1_squared']*u**(p['power']+1)
            + p['wet_term']*u**p['power'])


def reaction_derivative(u, p):
    return (2*p['nc']*u+p['nr']+(p['power']+1)*p['n1_squared']*u**p['power']
            + p['power']*p['wet_term']*u**(p['power']-1))


class HermiteFin:
    def __init__(self, cells, count):
        self.cells, self.count = cells, count
        self.edges = (1-np.cos(np.linspace(0, np.pi, cells+1)))/2
        self.polynomials = hermite_integrals(count, 1, 0)[0]
        size = cells*count
        constraints = []
        def row(cell, where, order):
            out = np.zeros(size)
            width = self.edges[cell+1]-self.edges[cell]
            out[cell*count:(cell+1)*count] = [f.deriv(order)(where)/width**order
                                                        for f in self.polynomials]
            return out
        constraints.append(row(0, 0, 0))
        constraints.append(row(cells-1, 1, 1))
        for cell in range(cells-1):
            for order in (0, 1):
                constraints.append(row(cell, 1, order)-row(cell+1, 0, order))
        constraint = np.asarray(constraints)
        constraint /= np.linalg.norm(constraint, axis=1)[:, None]
        z = null_space(constraint)
        # Gauss nodes never sample the degenerate endpoint.
        nodes, weights = np.polynomial.legendre.leggauss(32)
        qx, qw = [], []
        for left, right in zip(self.edges[:-1], self.edges[1:]):
            qx.extend(left+(nodes+1)*(right-left)/2)
            qw.extend(weights*(right-left)/2)
        self.qx, self.qw = np.asarray(qx), np.asarray(qw)
        raw, derivative, _ = self.raw(self.qx)
        v = raw @ z
        mass = v.T @ ((self.qw*(1-self.qx**2))[:, None]*v)
        self.mass_condition_number = float(np.linalg.cond(mass))
        chol = cholesky(mass, lower=True)
        self.transform = z @ solve_triangular(chol.T, np.eye(len(chol)), lower=False)
        self.v = raw @ self.transform
        self.d = derivative @ self.transform
        self.weights_a = self.qw*(1-self.qx**2)
        self.weights_s = self.qw*np.sqrt(1-self.qx**2)
        self.stiffness = self.d.T @ (self.weights_a[:, None]*self.d)
        # Weighted L2 projection of the zero interior initial state, using lift 1.
        self.initial = self.v.T @ (-self.weights_a)

    def raw(self, points):
        points = np.asarray(points)
        cell = np.clip(np.searchsorted(self.edges, points, side='right')-1, 0, self.cells-1)
        width = self.edges[cell+1]-self.edges[cell]
        local = (points-self.edges[cell])/width
        maps = []
        for order in (0, 1, 2):
            values = np.zeros((len(points), self.cells*self.count))
            for j, polynomial in enumerate(self.polynomials):
                values[np.arange(len(points)), cell*self.count+j] = polynomial.deriv(order)(local)/width**order
            maps.append(values)
        return maps

    def operators(self, points):
        return [v @ self.transform for v in self.raw(points)]

    def rhs(self, time, coefficients, p, source=None):
        temperature = 1+self.v @ coefficients
        result = -self.stiffness @ coefficients-self.v.T @ (self.weights_s*reaction(temperature, p))
        if source is not None:
            result += self.v.T @ (self.weights_a*source(self.qx, time))
        return result

    def integrate(self, p, initial=None, source=None):
        def jac(time, c):
            u = 1+self.v @ c
            return -self.stiffness-self.v.T @ ((self.weights_s*reaction_derivative(u, p))[:, None]*self.v)
        solution = solve_ivp(lambda t,c:self.rhs(t,c,p,source), (0,p['duration']),
                             self.initial if initial is None else initial,
                             method='BDF', jac=jac, rtol=2e-8, atol=2e-10, dense_output=True)
        if not solution.success or not np.isfinite(solution.y).all():
            raise RuntimeError('Hermite time integration failed: '+solution.message)
        return solution


def finite_volume(p, grid=400):
    """Independent cell-centred conservative solve; zero-area tip flux is zero."""
    dx = 1/grid
    centres = (np.arange(grid)+0.5)*dx
    faces = np.arange(grid+1)*dx
    # Exact cell integrals of A and sqrt(A).
    a_primitive = lambda x:x-x**3/3
    s_primitive = lambda x:0.5*(x*np.sqrt(np.maximum(0,1-x*x))+np.arcsin(x))
    mass = np.diff(a_primitive(faces))
    surface = np.diff(s_primitive(faces))
    conductance = (1-faces**2)/dx
    conductance[0] *= 2
    diag = -(conductance[:-1]+conductance[1:])/mass
    upper = conductance[1:-1]/mass[:-1]
    lower = conductance[1:-1]/mass[1:]
    matrix = diags([lower,diag,upper],[-1,0,1],format='csc')
    forcing = np.zeros(grid);forcing[0] = conductance[0]/mass[0]
    def rhs(t,u):return matrix @ u+forcing-surface/mass*reaction(u,p)
    def jac(t,u):return matrix-diags(surface/mass*reaction_derivative(u,p),format='csc')
    result=solve_ivp(rhs,(0,p['duration']),np.zeros(grid),method='BDF',jac=jac,
                     rtol=1e-8,atol=1e-10,dense_output=True)
    if not result.success:raise RuntimeError('Independent finite-volume solve failed.')
    return centres,result


def time_derivative(solution, time, duration):
    """Differentiate the completed BDF interpolant, independently of the RHS."""
    step=1e-5
    if time+step>duration:
        return (3*solution.sol(time)-4*solution.sol(time-step)+solution.sol(time-2*step))/(2*step)
    return (solution.sol(time+step)-solution.sol(time-step))/(2*step)


def compute(parameters=None):
    p=validate_parameters(DEFAULTS if parameters is None else parameters)
    model=HermiteFin(p['cells'],p['count'])
    solution=model.integrate(p)
    xs=np.linspace(0,1,161)
    ts=np.linspace(0,p['duration'],101)
    v,d,d2=model.operators(xs)
    coeff=solution.sol(ts)
    u=1+v @ coeff
    # The startup corner is a deliberate boundary step, not compatible smooth data.
    u[:,0]=0;u[0,0]=1
    boundary_error=float(np.max(np.abs(u[0]-1)))
    tip_error=float(np.max(np.abs(d[-1] @ coeff)))
    minimum=float(u[:,1:].min());maximum=float(u[:,1:].max())
    # Strong residual: report it as a diagnostic; Galerkin enforces weak equations.
    # Omit t=0 and exact interfaces where second derivatives are discontinuous.
    probe_x=np.concatenate([np.linspace(a,b,11)[1:-1] for a,b in zip(model.edges[:-1],model.edges[1:])])
    probe_t=np.linspace(0.05,p['duration'],21)
    pv,pd,pd2=model.operators(probe_x)
    cc=solution.sol(probe_t)
    uu=1+pv @ cc
    ut=pv @ np.column_stack([time_derivative(solution,t,p['duration']) for t in probe_t])
    residual=(1-probe_x[:,None]**2)*(ut-pd2 @ cc)+2*probe_x[:,None]*(pd @ cc)+np.sqrt(1-probe_x[:,None]**2)*reaction(uu,p)
    # Independent quadrature evaluates the conservative weak equation.
    nodes,weights=np.polynomial.legendre.leggauss(48)
    qx=np.concatenate([a+(nodes+1)*(b-a)/2 for a,b in zip(model.edges[:-1],model.edges[1:])])
    qw=np.concatenate([weights*(b-a)/2 for a,b in zip(model.edges[:-1],model.edges[1:])])
    qv,qd,_=model.operators(qx)
    weak=[]
    for time in probe_t:
        c=solution.sol(time);ct=time_derivative(solution,time,p['duration']);theta=1+qv @ c
        balance=qv.T @ (qw*(1-qx*qx)*(qv @ ct))+qd.T @ (qw*(1-qx*qx)*(qd @ c))+qv.T @ (qw*np.sqrt(1-qx*qx)*reaction(theta,p))
        weak.append(np.max(np.abs(balance)))
    report=dict(equation='(1-X^2)*theta_tau = ((1-X^2)*theta_X)_X - sqrt(1-X^2)*F(theta)',
                source='Pavan Kumar, Gireesha, Venkatesh (2026), DOI 10.1007/s10999-026-09893-6, Eq. (13)-(14)',
                method='C1 piecewise shifted Hermite weighted Galerkin; implicit BDF in time',
                parameters=p,parameter_note=PARAMETER_NOTE,
                basis_count_per_cell=p['count'],cells=p['cells'],unknown_coefficients=len(model.initial),
                space_interval=[0,1],time_interval=[0,p['duration']],
                initial_profile='theta(X,0)=0 for X>0; theta(0,0+)=1 (base step)',
                boundary_conditions='theta(0,tau)=1; theta_X(1,tau)=0',
                startup_note='Zero interior initial state is projected in weighted L2. The initial corner is discontinuous; polynomial startup traces are not pointwise exact.',
                initial_weighted_L2_projection_error=float(np.sqrt(np.dot(model.weights_a,(1+model.v @ model.initial)**2))),
                nonlinear_solver_message=solution.message,
                time_integrator_rtol=2e-8,time_integrator_atol=2e-10,
                mass_condition_number_before_orthonormalization=model.mass_condition_number,
                base_max_abs_error=boundary_error,tip_max_abs_derivative=tip_error,
                data_error_tolerance=DATA_TOLERANCE,
                independent_weak_residual=float(max(weak)),
                weak_residual_tolerance=WEAK_RESIDUAL_TOLERANCE,
                independent_weighted_strong_pde_residual=float(np.max(np.abs(residual))),
                weighted_strong_residual_tolerance=STRONG_RESIDUAL_TOLERANCE,
                residual_check_time_interval=[0.05,p['duration']],
                temporal_residual_check='Second-order finite differences of the completed BDF interpolant, step 1e-5',
                minimum_temperature_after_start=minimum,maximum_temperature_after_start=maximum,
                reference_used_during_solve=False,
                comparison_tolerance=COMPARISON_TOLERANCE,
                limitation='Weak-form approximation of dimensionless Eq. (13); no exact solution supplied by the paper. Finite-volume agreement and refinement are sampled checks. No full paper/material reproduction claimed.')
    arrays=dict(x=xs,t=ts,u=u,coefficients=coeff)
    ideal=(p['nc']+p['nr']+p['n1_squared']+p['wet_term'])*np.pi/4
    if ideal>0:
        efficiency=(model.weights_s @ reaction(1+model.v @ coeff,p))/ideal
        efficiency[0]=0
        arrays['efficiency']=efficiency
        report['final_fin_efficiency_eq16']=float(efficiency[-1])
    else:
        report['fin_efficiency_note']='Eq. (16) is undefined because its denominator is zero.'
    comparison_error=None
    if p['reference']:
        centres,reference=finite_volume(p)
        reference_values=reference.sol(ts)
        comparison=np.column_stack([np.interp(xs,np.r_[0,centres,1],np.r_[1,row,row[-1]]) for row in reference_values.T])
        comparison[:,0]=0;comparison[0,0]=1
        arrays['comparison']=comparison
        report['comparison_label']='Independent conservative finite-volume solution (400 cells), not exact'
        comparison_error=float(np.max(np.abs(u[:,ts>=0.05]-comparison[:,ts>=0.05])))
        report['max_abs_comparison_difference_after_start']=comparison_error
        report['comparison_time_interval']=[0.05,p['duration']]
    report['validation_passed']=bool(np.isfinite(u).all() and boundary_error<=DATA_TOLERANCE
        and tip_error<=DATA_TOLERANCE and min(u[:,ts>=0.05].ravel())>=-2e-3
        and max(u[:,ts>=0.05].ravel())<=1.002 and max(weak)<=WEAK_RESIDUAL_TOLERANCE
        and np.max(np.abs(residual))<=STRONG_RESIDUAL_TOLERANCE
        and (comparison_error is None or comparison_error<=COMPARISON_TOLERANCE))
    return arrays,report


def run(parameters=None, output=Path('results-fin')):
    from build_viewer import render_viewer,viewer_payload
    arrays,report=compute(parameters)
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output/'solution.npz',**arrays)
    (output/'results.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    (output/'interactive.html').write_text(render_viewer(viewer_payload(arrays,report)),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figure,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for j in (0,50,100):
        line,=axes[0].plot(arrays['x'],arrays['u'][:,j],label=f"tau={arrays['t'][j]:.2f}")
        if 'comparison' in arrays:
            axes[0].plot(arrays['x'][::8],arrays['comparison'][::8,j],'o',
                         color=line.get_color(),fillstyle='none',markersize=4)
    axes[0].set(xlabel='Dimensionless position X',ylabel='Dimensionless temperature theta',
                title='Hermite solution (circles: numerical comparison)' if 'comparison' in arrays else 'Hermite solution')
    axes[0].legend()
    if 'efficiency' in arrays:
        axes[1].plot(arrays['t'],arrays['efficiency'])
        axes[1].set(xlabel='Dimensionless time tau',ylabel='Fin efficiency eta',title='Eq. (16)')
    else:
        axes[1].text(.5,.5,'Efficiency undefined for zero losses',ha='center',va='center',transform=axes[1].transAxes)
    for axis in axes:axis.grid(alpha=.2)
    figure.suptitle('Semi-spherical porous fin - illustrative dimensionless coefficients')
    figure.savefig(output/'solution.png',dpi=160)
    plt.close(figure)
    print(json.dumps(report,indent=2))
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in DEFAULTS.keys()-{'reference'}:
        parser.add_argument('--'+key.replace('_','-'),type=int if key in ('cells','count','power') else float,default=DEFAULTS[key])
    parser.add_argument('--no-reference',action='store_true')
    parser.add_argument('--output',type=Path,default=Path('results-fin'))
    args=parser.parse_args()
    p={k:getattr(args,k) for k in DEFAULTS if k!='reference'};p['reference']=not args.no_reference
    try:report=run(p,args.output)
    except (ValueError,RuntimeError) as exc:parser.exit(1,str(exc)+'\n')
    if not report['validation_passed']:parser.exit(2,'Validation failed; diagnostic output retained.\n')

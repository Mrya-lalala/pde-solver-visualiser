"""Legendre-wavelet solver for Pavan Kumar et al. (2024), Eq. (14)-(15).

Inclined dovetail fin, DOI 10.1140/epjp/s13360-023-04848-8.
Localized orthonormal Legendre wavelets on dyadic cells are coupled by C1
constraints; a weighted Galerkin system is integrated by implicit BDF.
The printed dimensionless PDE is retained, including its powers and signs.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import Legendre
from scipy.integrate import solve_ivp
from scipy.linalg import cholesky, null_space, solve_triangular
from scipy.sparse import diags

DEFAULTS = dict(nc=1, nr=1, m_squared=1, taper=-0.1, q=0.1, epsilon_g=0.2,
                theta_a=0.2, wet_parameter=1, power=1, angle_degrees=45,
                phi1=0.02, phi2=0.02, phi3=0.02, count=6, cells=8,
                duration=0.8, reference=True)
PRESETS = {
    'paper_baseline': {k:v for k,v in DEFAULTS.items() if k not in ('count','cells','duration','reference')},
}
PRESETS['wider_dovetail'] = dict(PRESETS['paper_baseline'], taper=-0.3)
PRESETS['no_generation'] = dict(PRESETS['paper_baseline'], q=0)
PRESETS['linear_benchmark'] = dict(PRESETS['paper_baseline'], nc=0,nr=0,m_squared=1,
    taper=0,q=0,epsilon_g=0,theta_a=0,wet_parameter=0,power=0,angle_degrees=90,
    phi1=0,phi2=0,phi3=0)
STRONG_TOL = 2e-3
WEAK_TOL = 1e-5
COMPARISON_TOL = 2e-4
DATA_TOL = 1e-10


def validate_parameters(parameters):
    if not isinstance(parameters,dict) or set(parameters)!=set(DEFAULTS):
        raise ValueError('Provide exactly the documented dovetail-fin parameters.')
    p=dict(parameters)
    if type(p['reference']) is not bool:
        raise ValueError('reference must be a boolean.')
    for key in p.keys()-{'reference'}:
        if type(p[key]) not in (int,float) or not np.isfinite(p[key]):
            raise ValueError(f'{key} must be a finite number.')
    for key in ('cells','count','power'):
        if type(p[key]) is not int:
            raise ValueError(f'{key} must be an integer.')
    if p['cells'] not in (4,8,16) or p['count'] not in (4,6) or p['power'] not in (0,1):
        raise ValueError('Use 4/8/16 dyadic cells, 4/6 Legendre functions, and p=0/1.')
    if not -0.3<=p['taper']<=0 or not 0<=p['theta_a']<=0.3 or not 0.1<=p['duration']<=5:
        raise ValueError('Use C in [-0.3,0], ambient ratio in [0,0.3], and duration in [0.1,5].')
    if not 0<=p['angle_degrees']<=90 or any(not 0<=p[k]<=0.02 for k in ('phi1','phi2','phi3')):
        raise ValueError('Use inclination in [0,90] degrees and nanoparticle fractions in [0,0.02].')
    if any(not 0<=p[k]<=2 for k in ('nc','nr','m_squared','wet_parameter')) or not 0<=p['q']<=0.2 or not 0<=p['epsilon_g']<=0.2:
        raise ValueError('Unsupported coefficient range.')
    return p


def material_ratios(p):
    """Table 1 data and dimensionless ratios on the RHS of Table 2.

    Base fluid: 50/50 ethylene glycol-water. Particles: MWCNT, Ag, Cu.
    Table 2's density/heat/expansion RHSs are ratios despite dimensional LHS labels.
    """
    rho=np.array([1063.8,1600,10500,8930.])
    cp=np.array([3630,796,235,385.])
    beta=np.array([58,44,1.89,1.65])*1e-5
    phi1,phi2,phi3=(p[k] for k in ('phi1','phi2','phi3'))
    def nested(values):
        r=values/values[0]
        return float(phi1*r[1]+(1-phi1)*((1-phi2)*((1-phi3)+phi3*r[3])+phi2*r[2]))
    conductivity=0.37
    for fraction,particle in zip((phi1,phi2,phi3),(3000.,429.,400.)):
        conductivity*=((particle+2*conductivity-2*fraction*(conductivity-particle)) /
                       (particle+2*conductivity+fraction*(conductivity-particle)))
    return dict(rho=nested(rho),heat=nested(rho*cp),beta=nested(rho*beta),
                mu=float(((1-phi1)*(1-phi2)*(1-phi3))**-2.5),k=conductivity/0.37)


def coefficients(p):
    r=material_ratios(p);h=1/r['heat']
    return dict(diffusion=h*r['k'],generation=h*p['q'],
                convection=r['beta']/r['mu']*p['nc']*np.sin(np.deg2rad(p['angle_degrees'])),
                radiation=h*p['nr'],
                wet=h*p['m_squared']*(1+p['wet_parameter']*h*r['rho'])/(1-p['theta_a'])**p['power'])


def losses(theta,p):
    b=coefficients(p);s=theta-p['theta_a']
    return b['convection']*s*s+b['radiation']*s+b['wet']*s**(p['power']+1)


def loss_derivative(theta,p):
    b=coefficients(p);s=theta-p['theta_a']
    return 2*b['convection']*s+b['radiation']+(p['power']+1)*b['wet']*s**p['power']


class LegendreWavelets:
    def __init__(self,cells,count,taper):
        self.cells,self.count,self.taper=cells,count,taper
        self.edges=np.linspace(0,1,cells+1)
        self.polynomials=[Legendre.basis(m) for m in range(count)]
        self.normalization=np.sqrt((2*np.arange(count)+1)*cells)
        size=cells*count
        def row(cell,local,order):
            out=np.zeros(size)
            out[cell*count:(cell+1)*count]=[self.normalization[m]*f.deriv(order)(2*local-1)*(2*cells)**order
                                          for m,f in enumerate(self.polynomials)]
            return out
        constraints=[row(0,0,1),row(cells-1,1,0)]
        for cell in range(cells-1):
            for order in (0,1):constraints.append(row(cell,1,order)-row(cell+1,0,order))
        constraint=np.asarray(constraints)
        constraint/=np.linalg.norm(constraint,axis=1)[:,None]
        z=null_space(constraint)
        self.qx,self.qw=self.quadrature(24)
        raw,raw_d,_=self.raw(self.qx)
        v=raw @ z
        self.area=1-taper*self.qx
        mass=v.T @ ((self.qw*self.area)[:,None]*v)
        chol=cholesky(mass,lower=True)
        self.transform=z @ solve_triangular(chol.T,np.eye(len(chol)),lower=False)
        self.v=raw @ self.transform;self.d=raw_d @ self.transform
        self.stiffness=self.d.T @ ((self.qw*self.area)[:,None]*self.d)
        self.initial=-self.v.T @ (self.qw*self.area)

    def quadrature(self,order):
        nodes,weights=np.polynomial.legendre.leggauss(order)
        return (np.concatenate([a+(nodes+1)*(b-a)/2 for a,b in zip(self.edges[:-1],self.edges[1:])]),
                np.tile(weights/(2*self.cells),self.cells))

    def raw(self,points):
        points=np.asarray(points)
        cell=np.clip((points*self.cells).astype(int),0,self.cells-1)
        local=(points-self.edges[cell])*self.cells
        operators=[]
        for order in (0,1,2):
            values=np.zeros((len(points),self.cells*self.count))
            for m,f in enumerate(self.polynomials):
                values[np.arange(len(points)),cell*self.count+m]=self.normalization[m]*f.deriv(order)(2*local-1)*(2*self.cells)**order
            operators.append(values)
        return operators

    def operators(self,points):return [a @ self.transform for a in self.raw(points)]

    def integrate(self,p,initial=None,source=None):
        b=coefficients(p)
        def rhs(t,c):
            theta=1+self.v @ c
            reaction=b['generation']*self.area*(1+p['epsilon_g']*(theta-p['theta_a']))-losses(theta,p)
            if source is not None:reaction+=self.area*source(self.qx,t)
            return -b['diffusion']*self.stiffness @ c+self.v.T @ (self.qw*reaction)
        def jac(t,c):
            theta=1+self.v @ c
            derivative=b['generation']*self.area*p['epsilon_g']-loss_derivative(theta,p)
            return -b['diffusion']*self.stiffness+self.v.T @ ((self.qw*derivative)[:,None]*self.v)
        solution=solve_ivp(rhs,(0,p['duration']),self.initial if initial is None else initial,
                           method='BDF',jac=jac,rtol=1e-9,atol=1e-11,dense_output=True)
        if not solution.success or not np.isfinite(solution.y).all():raise RuntimeError('Legendre-wavelet integration failed: '+solution.message)
        return solution


def finite_volume(p,grid=400):
    """Independent conservative finite-volume calculation, insulated LEFT tip."""
    b=coefficients(p);dx=1/grid
    x=(np.arange(grid)+0.5)*dx;faces=np.arange(grid+1)*dx
    mass=dx*(1-p['taper']*x)
    flux=b['diffusion']*(1-p['taper']*faces)/dx
    flux[0]=0;flux[-1]*=2
    matrix=diags([flux[1:-1]/mass[1:],-(flux[:-1]+flux[1:])/mass,
                  flux[1:-1]/mass[:-1]],[-1,0,1],format='csc')
    boundary=np.zeros(grid);boundary[-1]=flux[-1]/mass[-1]
    area=1-p['taper']*x
    def rhs(t,u):return matrix @ u+boundary+b['generation']*(1+p['epsilon_g']*(u-p['theta_a']))-losses(u,p)/area
    def jac(t,u):return matrix+diags(b['generation']*p['epsilon_g']-loss_derivative(u,p)/area,format='csc')
    solution=solve_ivp(rhs,(0,p['duration']),np.zeros(grid),method='BDF',jac=jac,
                       rtol=2e-10,atol=2e-12,dense_output=True)
    if not solution.success:raise RuntimeError('Independent finite-volume integration failed.')
    return x,solution


def analytic_benchmark(x,t):
    """Post-solve reference for theta_t=theta_XX-theta, NOT the nonlinear case."""
    x=np.asarray(x);t=np.asarray(t)
    n=np.arange(160);lam=(n+0.5)*np.pi
    weights=-2*lam*(-1.)**n/(lam*lam+1)
    u=np.cosh(x)[:,None]/np.cosh(1)+(np.cos(x[:,None]*lam)*weights) @ np.exp(-(lam*lam+1)[:,None]*t)
    u[:,t==0]=0;u[np.ix_(x==1,t==0)]=1
    return u


def compute(parameters=None):
    p=validate_parameters(DEFAULTS if parameters is None else parameters)
    model=LegendreWavelets(p['cells'],p['count'],p['taper']);b=coefficients(p)
    solution=model.integrate(p)
    x=np.linspace(0,1,161);t=np.linspace(0,p['duration'],101)
    v,d,_=model.operators(x);cc=solution.sol(t);u=1+v @ cc
    u[:,0]=0;u[-1,0]=1
    # Differentiate the completed interpolant, not the integrator's RHS.
    def ct(t):
        step=1e-5
        return ((3*solution.sol(t)-4*solution.sol(t-step)+solution.sol(t-2*step))/(2*step)
                if t+step>p['duration'] else (solution.sol(t+step)-solution.sol(t-step))/(2*step))
    px=np.concatenate([np.linspace(a,b,13)[1:-1] for a,b in zip(model.edges[:-1],model.edges[1:])])
    pt=np.linspace(0.05,p['duration'],25)
    pv,pd,pd2=model.operators(px);a=1-p['taper']*px[:,None]
    c=solution.sol(pt);theta=1+pv @ c;derivative=pv @ np.column_stack([ct(s) for s in pt])
    residual=a*derivative-b['diffusion']*(a*(pd2 @ c)-p['taper']*(pd @ c))-b['generation']*a*(1+p['epsilon_g']*(theta-p['theta_a']))+losses(theta,p)
    qx,qw=model.quadrature(36);qv,qd,_=model.operators(qx);qa=1-p['taper']*qx
    weak=[]
    for s in pt:
        c=solution.sol(s);theta=1+qv @ c
        reaction=losses(theta,p)-b['generation']*qa*(1+p['epsilon_g']*(theta-p['theta_a']))
        balance=qv.T @ (qw*qa*(qv @ ct(s)))+b['diffusion']*qd.T @ (qw*qa*(qd @ c))+qv.T @ (qw*reaction)
        weak.append(float(np.max(abs(balance))))
    base_error=float(np.max(abs(u[-1]-1)));tip_error=float(np.max(abs(d[0] @ cc)))
    report=dict(problem='dovetail_fin',basis='legendre',parameters=p,
        source='Pavan Kumar, Gireesha, Venkatesh (2024), DOI 10.1140/epjp/s13360-023-04848-8, Eq. (14)-(15)',
        equation='a*theta_tau = D*(a*theta_X)_X + g*a*(1+epsilon_g*s) - B*s^2 - R*s - W*s^(p+1); a=1-C*X; s=theta-theta_a',
        method='Localized orthonormal Legendre wavelets on dyadic cells; C1-constrained weighted Galerkin; BDF time integration',
        material_ratios=material_ratios(p),effective_coefficients=b,wavelet_level=int(np.log2(p['cells']))+1,
        unknown_coefficients=len(model.initial),space_interval=[0,1],time_interval=[0,p['duration']],
        initial_profile='theta(X,0)=0 for X<1; theta(1,0+)=1',
        boundary_conditions='theta_X(0,tau)=0 (insulated tip); theta(1,tau)=1 (base)',
        startup_note='The zero interior initial state is projected in weighted L2; the base steps to one at X=1. Saved time zero shows the prescribed data. Startup before tau=0.05 is not certified.',
        initial_weighted_L2_projection_error=float(np.sqrt(np.dot(model.qw*model.area,(1+model.v @ model.initial)**2))),
        nonlinear_solver_message=solution.message,time_integrator_rtol=1e-9,time_integrator_atol=1e-11,
        independent_weighted_strong_pde_residual=float(np.max(abs(residual))),weighted_strong_residual_tolerance=STRONG_TOL,
        independent_weak_residual=max(weak),weak_residual_tolerance=WEAK_TOL,
        base_max_abs_error=base_error,tip_max_abs_derivative=tip_error,data_error_tolerance=DATA_TOL,
        residual_check_time_interval=[0.05,p['duration']],reference_used_during_solve=False,
        interpretation='Printed dimensionless Eq. (14); Table 2 RHSs interpreted as property ratios; baseline G=0.1 interpreted as Q=0.1. Printed dimensional-time scaling is unresolved.',
        limitation='Numerical solution of the stated dimensionless problem; sampled checks are not rigorous bounds or full reproduction of the material curves.')
    arrays=dict(x=x,t=t,u=u,coefficients=cc,wavelet_coefficients=model.transform @ cc,cell_edges=model.edges)
    if p['reference']:
        centres,reference=finite_volume(p)
        values=reference.sol(t)
        comparison=np.column_stack([np.interp(x,np.r_[0,centres,1],np.r_[row[0],row,1]) for row in values.T])
        comparison[:,0]=0;comparison[-1,0]=1
        arrays['comparison']=comparison
        report['comparison_label']='Independent conservative finite-volume solution (400 cells), not exact'
        report['max_abs_comparison_difference_after_start']=float(np.max(abs(u[:,t>=0.05]-comparison[:,t>=0.05])))
        report['comparison_tolerance']=COMPARISON_TOL
    report['validation_passed']=bool(np.isfinite(u).all() and base_error<=DATA_TOL and tip_error<=DATA_TOL
        and max(weak)<=WEAK_TOL and np.max(abs(residual))<=STRONG_TOL
        and report.get('max_abs_comparison_difference_after_start',0)<=COMPARISON_TOL)
    return arrays,report


def run(parameters=None,output=Path('results-dovetail')):
    from build_viewer import render_viewer,viewer_payload
    arrays,report=compute(parameters);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output/'solution.npz',**arrays)
    (output/'results.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (output/'interactive.html').write_text(render_viewer(viewer_payload(arrays,report)),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(8,4.8),layout='constrained')
    for j in (0,25,50,100):
        line,=ax.plot(arrays['x'],arrays['u'][:,j],label=f"tau={arrays['t'][j]:.2f}")
        if 'comparison' in arrays:ax.plot(arrays['x'][::10],arrays['comparison'][::10,j],'o',color=line.get_color(),fillstyle='none',ms=4)
    ax.set(xlabel='Dimensionless position X (tip 0, base 1)',ylabel='Dimensionless temperature theta',title='Inclined dovetail fin: Legendre-wavelet numerical solution')
    ax.grid(alpha=.2);ax.legend();fig.savefig(output/'solution.png',dpi=160);plt.close(fig)
    print(json.dumps(report,indent=2))
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preset',choices=PRESETS,default='paper_baseline')
    parser.add_argument('--cells',type=int,default=8);parser.add_argument('--count',type=int,default=6)
    parser.add_argument('--time',type=float,default=.8);parser.add_argument('--no-reference',action='store_true')
    parser.add_argument('--output',type=Path,default=Path('results-dovetail'))
    args=parser.parse_args()
    p=dict(DEFAULTS,**PRESETS[args.preset]);p.update(cells=args.cells,count=args.count,duration=args.time,reference=not args.no_reference)
    try:report=run(p,args.output)
    except (ValueError,RuntimeError) as error:parser.exit(1,str(error)+'\n')
    if not report['validation_passed']:parser.exit(2,'Independent validation failed; diagnostics saved.\n')

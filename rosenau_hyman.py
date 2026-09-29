#!/usr/bin/env python3
"""Rosenau-Hyman benchmark using single-cell (k=1) Hermite collocation.

Solve u_t = u*u_xxx + u*u_x + 3*u_x*u_xx.
Run: python rosenau_hyman.py --output results
Dependencies: numpy, scipy, sympy, matplotlib.

This is a benchmark implementation, not an arbitrary-PDE solver. The prescribed
u, u_x and u_xx traces at x=0 supplement the initial profile. Their values are
explicit inputs, NOT extracted from EXACT_SOLUTION by the solver.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp
from numpy.polynomial import Polynomial
from scipy.optimize import root

# ========================= EDITABLE INPUTS =============================
x, t = sp.symbols('x t', real=True)
SPACE_END = 1.0                  # Domain: 0 <= x <= SPACE_END
TIME_END = 1.0                   # Time:   0 <= t <= TIME_END
BASIS_COUNT = 6                  # Functions per variable; total unknowns = N*N
ROOT_TOLERANCE = 1e-10
COLLOCATION_TOLERANCE = 1e-8
VALIDATION_RESIDUAL_TOLERANCE = 1e-5

# Try SPEED = 0.5, 1, or 2.
SPEED = sp.Rational(1)
SHIFT = sp.Rational(0)

INITIAL_PROFILE = (
    -8*SPEED/3 * sp.cos((x - SHIFT)/4)**2
)

BOUNDARY_U = (
    -8*SPEED/3 * sp.cos((SPEED*t + SHIFT)/4)**2
)

BOUNDARY_UX = (
    -2*SPEED/3 * sp.sin((SPEED*t + SHIFT)/2)
)

BOUNDARY_UXX = (
    SPEED/3 * sp.cos((SPEED*t + SHIFT)/2)
)

# Used only for checking numerical error.
EXACT_SOLUTION = (
    -8*SPEED/3 * sp.cos((x - SPEED*t - SHIFT)/4)**2
)
# ======================================================================


def numerical_function(expression):
    """Broadcast constant expressions as well as expressions in x and t."""
    fn = sp.lambdify((x, t), expression, modules='numpy')

    def evaluate(X, T):
        X, T = np.broadcast_arrays(np.asarray(X, float), np.asarray(T, float))
        return np.broadcast_to(np.asarray(fn(X, T), float), X.shape)
    return evaluate


def construct_base(initial, traces):
    """Integrate u_xxxt and impose data; compute all needed derivatives."""
    initial = sp.sympify(initial)
    traces = tuple(sp.sympify(g) for g in traces)
    if initial.free_symbols - {x}:
        raise ValueError('Initial profile must depend only on x.')
    for order, trace in enumerate(traces):
        if trace.free_symbols - {t}:
            raise ValueError('Boundary traces must depend only on t.')
        mismatch = sp.simplify(sp.diff(initial, x, order).subs(x, 0)
                               - trace.subs(t, 0))
        if mismatch != 0 and abs(float(mismatch)) > 1e-10:
            raise ValueError(f'Incompatible initial/boundary data at derivative {order}.')
    base = initial + sum(
        x**m / sp.factorial(m) * (g - g.subs(t, 0))
        for m, g in enumerate(traces)
    )
    spatial = [numerical_function(sp.diff(base, x, i)) for i in range(4)]
    temporal = numerical_function(sp.diff(base, t))
    return spatial, temporal


def hermite_integrals(count, length, max_order):
    """Shifted physicists' Hermite basis and exact polynomial antiderivatives.

    k=1 uses one interval. Scaling each basis function improves numerical
    conditioning and changes coefficients, not the approximation space.
    Every repeated integral has lower endpoint 0; interval scaling is automatic.
    """
    z = Polynomial([-1.0, 2.0 / length])
    basis = [Polynomial([1.0]), 2 * z]
    for m in range(2, count):
        basis.append(2 * z * basis[-1] - 2 * (m - 1) * basis[-2])
    sample = np.linspace(0, length, 501)
    basis = [p / np.max(np.abs(p(sample))) for p in basis[:count]]
    return [[p if order == 0 else p.integ(order, lbnd=0)
             for p in basis] for order in range(max_order + 1)]


class HermiteProblem:
    def __init__(self, count, length, duration, initial, traces):
        if not 2 <= count <= 10:
            raise ValueError('This small dense implementation supports 2 <= N <= 10.')
        if not np.isfinite([length, duration]).all() or min(length, duration) <= 0:
            raise ValueError('Spatial and time intervals must have positive finite lengths.')
        self.count, self.length, self.duration = count, length, duration
        self.spatial, self.temporal = construct_base(initial, traces)
        self.x_integrals = hermite_integrals(count, length, 3)
        self.t_integrals = hermite_integrals(count, duration, 1)

    def operators(self, X, T):
        X, T = np.broadcast_arrays(np.asarray(X, float), np.asarray(T, float))
        ax = [np.stack([p(X) for p in row], axis=-1) for row in self.x_integrals]
        at = [np.stack([p(T) for p in row], axis=-1) for row in self.t_integrals]
        def tensor(a, b):
            return np.einsum('...i,...j->...ij', a, b).reshape(X.shape + (-1,))
        # Maps flattened coefficients to u, u_x, u_xx, u_xxx, u_t.
        maps = [tensor(ax[3-i], at[1]) for i in range(4)]
        maps.append(tensor(ax[3], at[0]))
        bases = [fn(X, T) for fn in self.spatial] + [self.temporal(X, T)]
        return maps, bases

    @staticmethod
    def evaluate(coefficients, operators):
        maps, bases = operators
        return [b + a @ coefficients for a, b in zip(maps, bases)]

    @staticmethod
    def residual(values):
        u, ux, uxx, uxxx, ut = values
        return ut - u * uxxx - u * ux - 3 * ux * uxx

    def solve(self):
        q = (np.arange(self.count) + 0.5) / self.count
        X, T = np.meshgrid(q*self.length, q*self.duration, indexing='ij')
        ops = self.operators(X.ravel(), T.ravel())
        maps, _ = ops
        def fun(c):
            return self.residual(self.evaluate(c, ops))
        def jac(c):
            u, ux, uxx, uxxx, _ = self.evaluate(c, ops)
            a0, a1, a2, a3, at = maps
            return (at - (uxxx+ux)[:, None]*a0 - u[:, None]*(a3+a1)
                    - 3*uxx[:, None]*a1 - 3*ux[:, None]*a2)
        result = root(fun, np.zeros(self.count**2), jac=jac,
                      method='hybr', options={'xtol': ROOT_TOLERANCE, 'maxfev': 4000})
        residual = float(np.max(np.abs(fun(result.x))))
        if not np.isfinite(result.x).all() or not np.isfinite(residual) or residual > COLLOCATION_TOLERANCE:
            raise RuntimeError(f'Nonlinear solve failed: {result.message}; residual={residual:.3e}. '
                               'Try fewer basis functions or a shorter time interval.')
        message = str(result.message)
        if not result.success:
            message = 'Residual tolerance satisfied despite optimizer status: ' + message
        return result.x, residual, message


def compute(count=BASIS_COUNT, length=SPACE_END, duration=TIME_END,
            initial=INITIAL_PROFILE, traces=None, exact=EXACT_SOLUTION):
    """Compute and validate in memory; reference is used only after solving."""
    # The solver never receives the optional exact solution.
    traces = (BOUNDARY_U, BOUNDARY_UX, BOUNDARY_UXX) if traces is None else traces
    problem = HermiteProblem(count, length, duration, initial, traces)
    coefficients, coll_residual, message = problem.solve()
    xs, ts = np.linspace(0, length, 121), np.linspace(0, duration, 101)
    X, T = np.meshgrid(xs, ts, indexing='ij')
    values = problem.evaluate(coefficients, problem.operators(X, T))
    solution = values[0]
    validation_residual = float(np.max(np.abs(problem.residual(values))))
    if not all(np.isfinite(v).all() for v in values) or not np.isfinite(validation_residual):
        raise RuntimeError('Non-finite values in independent validation.')
    initial_error = float(np.max(np.abs(solution[:, 0] - numerical_function(initial)(xs, 0))))
    boundary_errors = [float(np.max(np.abs(values[i][0, :] - numerical_function(g)(0, ts))))
                       for i, g in enumerate(traces)]
    report = {
        'equation': 'u_t = u*u_xxx + u*u_x + 3*u_x*u_xx',
        'method': 'k=1 shifted Hermite, mixed-derivative integration, tensor midpoint collocation',
        'basis_count_per_variable': count, 'unknown_coefficients': count**2,
        'space_interval': [0, length], 'time_interval': [0, duration],
        'initial_profile': str(initial),
        'boundary_traces_u_ux_uxx': [str(g) for g in traces],
        'nonlinear_solver_message': message,
        'collocation_max_abs_residual': coll_residual,
        'collocation_residual_tolerance': COLLOCATION_TOLERANCE,
        'validation_grid_shape': list(X.shape),
        'validation_max_abs_pde_residual': validation_residual,
        'validation_residual_tolerance': VALIDATION_RESIDUAL_TOLERANCE,
        'validation_passed': bool(validation_residual <= VALIDATION_RESIDUAL_TOLERANCE
                                  and initial_error <= 1e-12 and max(boundary_errors) <= 1e-12),
        'data_error_tolerance': 1e-12,
        'initial_max_abs_error': initial_error,
        'boundary_max_abs_errors_u_ux_uxx': boundary_errors,
        'reference_used_during_solve': False,
        'limitation': 'Benchmark with prescribed u, u_x, u_xx at x=0; no general well-posedness guarantee.'
    }
    reference = None if exact is None else numerical_function(exact)(X, T)
    if reference is not None:
        report['max_abs_solution_error'] = float(np.max(np.abs(solution-reference)))
        report['rms_solution_error'] = float(np.sqrt(np.mean((solution-reference)**2)))
    arrays = dict(x=xs, t=ts, u=solution, coefficients=coefficients.reshape(count, count))
    if reference is not None:
        arrays['exact'] = reference
    return arrays, report, problem.residual(values)


def travelling_wave(speed, shift):
    """Generate compatible benchmark data, independently of reference evaluation."""
    c, a = sp.Rational(str(speed)), sp.Rational(str(shift))
    initial = -8*c/3 * sp.cos((x-a)/4)**2
    traces = (-8*c/3 * sp.cos((c*t+a)/4)**2,
              -2*c/3 * sp.sin((c*t+a)/2), c/3 * sp.cos((c*t+a)/2))
    reference = -8*c/3 * sp.cos((x-c*t-a)/4)**2
    return initial, traces, reference


def run(count=BASIS_COUNT, length=SPACE_END, duration=TIME_END,
        output=Path('results-original'), exact=EXACT_SOLUTION):
    arrays, report, residual_grid = compute(count, length, duration, exact=exact)
    xs, ts, solution = arrays['x'], arrays['t'], arrays['u']
    reference = arrays.get('exact')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/'solution.npz', **arrays)
    (output/'results.json').write_text(json.dumps(report, indent=2) + '\n')

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    colors = ['#2563eb', '#b45309', '#15803d']
    for j, color in zip([0, 50, 100], colors):
        axes[0].plot(xs, solution[:, j], color=color, lw=2, label=f't = {ts[j]:.2f}')
        if reference is not None:
            axes[0].plot(xs[::6], reference[::6, j], 'o', ms=4, mfc='none', color=color)
            axes[1].semilogy(xs, np.maximum(np.abs(solution[:, j]-reference[:, j]), 1e-16),
                             color=color, label=f't = {ts[j]:.2f}')
        else:
            axes[1].semilogy(xs, np.maximum(np.abs(residual_grid[:, j]), 1e-16),
                             color=color, label=f't = {ts[j]:.2f}')
    axes[0].set(title='Hermite solution (circles: exact)' if reference is not None else 'Hermite solution',
                xlabel='Position x', ylabel='u(x,t)')
    axes[1].set(title='Absolute solution error' if reference is not None else 'Absolute PDE residual',
                xlabel='Position x', ylabel='Absolute error (floor 1e-16)' if reference is not None else 'Absolute residual (floor 1e-16)')
    for ax in axes:
        ax.grid(alpha=.2); ax.legend()
    fig.suptitle(f'Rosenau-Hyman benchmark | {count} x {count} Hermite coefficients')
    fig.savefig(output/'solution.png', dpi=180)
    plt.close(fig)
    print(json.dumps(report, indent=2))
    print('\nSample numerical values:')
    print(' x       t=0          t=T/2        t=T')
    for i in [0, 30, 60, 90, 120]:
        print(f'{xs[i]:.2f}  ' + '  '.join(f'{solution[i,j]: .9f}' for j in [0, 50, 100]))
    print(f'\nOutputs: {output.resolve()}')
    if not report['validation_passed']:
        print('WARNING: independent PDE residual exceeds tolerance; inspect/refine before trusting the result.')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--basis', type=int, default=BASIS_COUNT)
    parser.add_argument('--length', type=float, default=SPACE_END)
    parser.add_argument('--time', type=float, default=TIME_END)
    parser.add_argument('--output', type=Path, default=Path('results-original'))
    parser.add_argument('--no-reference', action='store_true', help='Disable exact-solution comparison.')
    args = parser.parse_args()
    try:
        result = run(args.basis, args.length, args.time, args.output,
                     None if args.no_reference else EXACT_SOLUTION)
    except (ValueError, RuntimeError) as exc:
        parser.exit(1, f'Error: {exc}\n')
    if not result['validation_passed']:
        parser.exit(2, 'Validation failed; outputs retained for diagnosis.\n')

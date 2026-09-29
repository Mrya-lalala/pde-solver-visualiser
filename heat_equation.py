#!/usr/bin/env python3
"""Heat-equation benchmark using single-cell (k=1) Hermite collocation.

Solve
    u_t = alpha * u_xx            on  0 <= x <= L,  0 <= t <= T,
with
    u(x, 0) = sin(pi*x/L),        u(0, t) = u(L, t) = 0,
and the reference
    u(x, t) = exp(-alpha*(pi/L)**2 * t) * sin(pi*x/L).

Run:  python heat_equation.py --output results-heat

Reconstruction (every integral starts at the lower endpoint 0)
    q(x, t) = u_xxt(x, t) ~ phi_x(x)^T C phi_t(t)
    A_i(x)  = I_x^2(phi_i)(x) - (x/L) * I_x^2(phi_i)(L)
    B_j(t)  = I_t(phi_j)(t)
    u(x, t) = f(x) + A(x)^T C B(t)

Why the endpoint values are built in:
  * B_j(0) = 0 for every j (the antiderivative has zero constant term), so
    u(x, 0) = f(x) exactly: the initial condition needs no extra equation.
  * A_i(0) = I_x^2(phi_i)(0) - 0 = 0, and
    A_i(L) = I_x^2(phi_i)(L) - (L/L) I_x^2(phi_i)(L) = 0,
    so u(0, t) = f(0) = 0 and u(L, t) = f(L) = sin(pi) = 0 for every t.
    Both homogeneous endpoints are satisfied without prescribing boundary
    slopes (unlike the Rosenau-Hyman benchmark, which needs u_x and u_xx).

Differentiating the reconstruction gives u_t = A^T C phi_t and
u_xx = f'' + phi_x^T C B (the linear endpoint term is annihilated by two
x-derivatives). Collocating u_t - alpha*u_xx = 0 at the tensor product of the
N midpoints of each interval therefore yields the linear system

    [ A(x) (x) phi_t(t) - alpha * phi_x(x) (x) B(t) ] c = alpha * f''(x),

where (x) is the Kronecker/outer product over the two one-dimensional bases.
It is solved directly with numpy.linalg.solve; no nonlinear optimizer is used
and the reference solution never enters assembly or the solve.

This is a fixed one-dimensional benchmark, not a general PDE solver. The
single-cell shifted physicists' Hermite basis (and its rescaling) is reused
unchanged from rosenau_hyman.hermite_integrals.

Dependencies: numpy, matplotlib (the current requirements.txt also carries
scipy and sympy, which this module does not need).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial

from rosenau_hyman import hermite_integrals

# ========================= EDITABLE INPUTS =============================
DEFAULT_ALPHA = 1.0              # Diffusivity, must be finite and > 0
DEFAULT_LENGTH = 1.0             # Domain: 0 <= x <= L
DEFAULT_TIME = 1.0               # Time:   0 <= t <= T
DEFAULT_BASIS = 16               # Functions per variable; unknowns = N*N

MIN_BASIS = 2
MAX_BASIS = 24                   # Above ~18 conditioning dominates (see README)
SOLUTION_ERROR_TOLERANCE = 1e-6
PDE_RESIDUAL_TOLERANCE = 1e-5
DATA_ERROR_TOLERANCE = 1e-12
REFINEMENT_SIZES = (8, 10, 12, 14, 16, 18, 20)
# ======================================================================


def initial_profile(x_values, length):
    """u(x, 0) = sin(pi*x/L); also the t=0 base f(x) of the reconstruction."""
    return np.sin(np.pi * np.asarray(x_values, float) / length)


def initial_second_derivative(x_values, length):
    """f''(x) for f(x) = sin(pi*x/L); the right-hand side source term."""
    return -((np.pi / length) ** 2) * initial_profile(x_values, length)


def exact_solution(x_values, t_values, alpha, length):
    """Reference solution, used only for post-solve diagnostics."""
    return (np.exp(-alpha * (np.pi / length) ** 2 * np.asarray(t_values, float))
            * np.sin(np.pi * np.asarray(x_values, float) / length))


def validate_parameters(count, length, duration, alpha):
    """Reject out-of-range or non-finite inputs before any assembly."""
    if not isinstance(count, (int, np.integer)) or isinstance(count, bool):
        raise ValueError('basis count must be an integer')
    if not MIN_BASIS <= count <= MAX_BASIS:
        raise ValueError(f'basis count must satisfy {MIN_BASIS} <= N <= {MAX_BASIS}')
    for name, value in (('alpha', alpha), ('length', length), ('duration', duration)):
        if not np.isfinite(value):
            raise ValueError(f'{name} must be finite')
        if value <= 0:
            raise ValueError(f'{name} must be strictly positive')


def build_operators(count, length, duration):
    """Basis functions and the endpoint-constrained x antiderivatives A, B.

    Reuses rosenau_hyman.hermite_integrals so the shifted physicists' Hermite
    basis and its per-interval rescaling are identical to the other benchmark.
    """
    spatial = hermite_integrals(count, length, 2)   # [phi, I_x phi, I_x^2 phi]
    temporal = hermite_integrals(count, duration, 1)  # [phi_t, I_t phi_t]
    phi_x, ix2 = spatial[0], spatial[2]
    phi_t, b_t = temporal[0], temporal[1]
    # A_i(x) = I_x^2(phi_i)(x) - (x/L) * I_x^2(phi_i)(L) vanishes at x=0 and x=L.
    a_x = [ix2[i] - Polynomial([0.0, ix2[i](length) / length]) for i in range(count)]
    return phi_x, phi_t, a_x, b_t


def collocation_nodes(count, length, duration):
    """Tensor product of the N interval midpoints (N*N collocation points)."""
    q = (np.arange(count) + 0.5) / count
    grid_x, grid_t = np.meshgrid(q * length, q * duration, indexing='ij')
    return grid_x.ravel(), grid_t.ravel()


def independent_grid(end, node_count, basis_count):
    """Uniform diagnostic grid with every collocation node removed.

    Sampling the PDE residual at collocation points would report a value that
    is zero by construction, so the validation grid is kept disjoint from them.
    Endpoints of the interval are always retained.
    """
    nodes = np.linspace(0.0, end, node_count)
    midpoints = (np.arange(basis_count) + 0.5) / basis_count * end
    nodes = nodes[np.min(np.abs(nodes[:, None] - midpoints[None, :]), axis=1) > 1e-9]
    return np.unique(np.concatenate(([0.0], nodes, [end])))


def _tensor(rows_a, rows_b):
    """Outer product per row: (P, n), (P, n) -> (P, n*n)."""
    return np.einsum('pi,pj->pij', rows_a, rows_b).reshape(rows_a.shape[0], -1)


def assemble(count, length, duration, alpha):
    """Build the linear collocation system M c = rhs for the coefficients."""
    phi_x, phi_t, a_x, b_t = build_operators(count, length, duration)
    nodes_x, nodes_t = collocation_nodes(count, length, duration)
    ax = np.stack([p(nodes_x) for p in a_x], axis=-1)
    bt = np.stack([p(nodes_t) for p in b_t], axis=-1)
    px = np.stack([p(nodes_x) for p in phi_x], axis=-1)
    pt = np.stack([p(nodes_t) for p in phi_t], axis=-1)
    matrix = _tensor(ax, pt) - alpha * _tensor(px, bt)
    rhs = alpha * initial_second_derivative(nodes_x, length)
    return matrix, rhs, (phi_x, phi_t, a_x, b_t)


def solve_linear_system(matrix, rhs):
    """Direct solve with explicit handling of singular and non-finite systems."""
    if not np.isfinite(matrix).all() or not np.isfinite(rhs).all():
        raise RuntimeError('Collocation system contains non-finite entries.')
    condition_number = float(np.linalg.cond(matrix))
    try:
        coefficients = np.linalg.solve(matrix, rhs)
    except np.linalg.LinAlgError as exc:
        raise RuntimeError(f'Collocation matrix is singular; no solution: {exc}') from exc
    if not np.isfinite(coefficients).all():
        raise RuntimeError('Linear solve produced non-finite coefficients.')
    algebraic_residual = float(np.max(np.abs(matrix @ coefficients - rhs)))
    return coefficients, condition_number, algebraic_residual


def evaluate_fields(operators, coefficients, x_values, t_values, length, alpha):
    """Return u, u_t and u_xx from the reconstructed coefficient matrix."""
    phi_x, phi_t, a_x, b_t = operators
    count = len(a_x)
    matrix_c = coefficients.reshape(count, count)
    ax = np.stack([p(x_values) for p in a_x], axis=-1)
    bt = np.stack([p(t_values) for p in b_t], axis=-1)
    px = np.stack([p(x_values) for p in phi_x], axis=-1)
    pt = np.stack([p(t_values) for p in phi_t], axis=-1)
    u = initial_profile(x_values, length) + np.einsum('...i,ij,...j->...', ax, matrix_c, bt)
    u_t = np.einsum('...i,ij,...j->...', ax, matrix_c, pt)
    u_xx = (initial_second_derivative(x_values, length)
            + np.einsum('...i,ij,...j->...', px, matrix_c, bt))
    return u, u_t, u_xx


def refinement_sweep(length=DEFAULT_LENGTH, duration=DEFAULT_TIME, alpha=DEFAULT_ALPHA,
                     sizes=REFINEMENT_SIZES):
    """Solve at several basis sizes to show convergence before round-off."""
    rows = []
    for count in sizes:
        matrix, rhs, operators = assemble(count, length, duration, alpha)
        coefficients, condition_number, algebraic = solve_linear_system(matrix, rhs)
        xs = independent_grid(length, 121, count)
        ts = independent_grid(duration, 101, count)
        grid_x, grid_t = np.meshgrid(xs, ts, indexing='ij')
        u, u_t, u_xx = evaluate_fields(operators, coefficients, grid_x, grid_t, length, alpha)
        residual = float(np.max(np.abs(u_t - alpha * u_xx)))
        reference = exact_solution(grid_x, grid_t, alpha, length)
        rows.append({
            'basis_count': count,
            'unknown_coefficients': count * count,
            'max_abs_solution_error': float(np.max(np.abs(u - reference))),
            'validation_max_abs_pde_residual': residual,
            'collocation_max_abs_residual': algebraic,
            'condition_number_2norm': condition_number,
        })
    return rows


def run(count=DEFAULT_BASIS, length=DEFAULT_LENGTH, duration=DEFAULT_TIME,
        alpha=DEFAULT_ALPHA, output=Path('results-heat'), compare_reference=True):
    """Solve, validate on an independent grid, and write NPZ/JSON/PNG outputs."""
    validate_parameters(count, length, duration, alpha)
    matrix, rhs, operators = assemble(count, length, duration, alpha)
    coefficients, condition_number, algebraic_residual = solve_linear_system(matrix, rhs)

    xs = independent_grid(length, 121, count)
    ts = independent_grid(duration, 101, count)
    grid_x, grid_t = np.meshgrid(xs, ts, indexing='ij')
    solution, u_t, u_xx = evaluate_fields(operators, coefficients, grid_x, grid_t, length, alpha)
    residual = u_t - alpha * u_xx
    validation_residual = float(np.max(np.abs(residual)))
    if not np.isfinite(solution).all() or not np.isfinite(validation_residual):
        raise RuntimeError('Non-finite values in independent validation.')

    initial_error = float(np.max(np.abs(solution[:, 0] - initial_profile(xs, length))))
    endpoint_errors = [float(np.max(np.abs(solution[0, :]))), float(np.max(np.abs(solution[-1, :])))]
    data_error = max(initial_error, *endpoint_errors)

    reference = None
    if compare_reference:
        reference = exact_solution(grid_x, grid_t, alpha, length)
        solution_error = float(np.max(np.abs(solution - reference)))
        rms_error = float(np.sqrt(np.mean((solution - reference) ** 2)))
    else:
        solution_error = rms_error = None

    # Conditioning and accuracy are reported separately: the system is very
    # ill-conditioned, yet the sampled accuracy is checked independently.
    ill_conditioned = condition_number * float(np.finfo(float).eps) > 1.0
    passed = (validation_residual <= PDE_RESIDUAL_TOLERANCE
              and data_error <= DATA_ERROR_TOLERANCE)
    if solution_error is not None:
        passed = passed and solution_error <= SOLUTION_ERROR_TOLERANCE

    report = {
        'equation': 'u_t = alpha*u_xx',
        'method': 'k=1 shifted Hermite, endpoint-constrained integration, tensor midpoint collocation',
        'basis_count_per_variable': count,
        'unknown_coefficients': count * count,
        'parameters': {'alpha': alpha, 'length': length, 'duration': duration},
        'initial_profile': 'sin(pi*x/L)',
        'boundary_conditions': 'u(0,t) = u(L,t) = 0',
        'exact_solution': 'exp(-alpha*(pi/L)**2*t)*sin(pi*x/L)',
        'linear_solve': 'numpy.linalg.solve (direct)',
        'collocation_max_abs_residual': algebraic_residual,
        'condition_number_2norm': condition_number,
        'condition_number_times_eps': condition_number * float(np.finfo(float).eps),
        'conditioning_warning': bool(ill_conditioned),
        'validation_grid_shape': list(grid_x.shape),
        'validation_max_abs_pde_residual': validation_residual,
        'validation_residual_tolerance': PDE_RESIDUAL_TOLERANCE,
        'solution_error_tolerance': SOLUTION_ERROR_TOLERANCE,
        'data_error_tolerance': DATA_ERROR_TOLERANCE,
        'initial_max_abs_error': initial_error,
        'endpoint_max_abs_errors_u0_uL': endpoint_errors,
        'reference_used_during_solve': False,
        'reference_compared_after_solve': bool(compare_reference),
        'validation_passed': bool(passed),
        'limitation': ('Fixed 1D benchmark. Accuracy degrades when alpha*(pi/L)**2*T '
                       'is large because the fast exponential decay is hard for a '
                       'low-degree time polynomial. High condition numbers do not by '
                       'themselves imply the sampled errors above.'),
    }
    if solution_error is not None:
        report['max_abs_solution_error'] = solution_error
        report['rms_solution_error'] = rms_error

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    arrays = dict(x=xs, t=ts, u=solution, coefficients=coefficients.reshape(count, count))
    if reference is not None:
        arrays['exact'] = reference
    np.savez_compressed(output / 'solution.npz', **arrays)
    (output / 'results.json').write_text(json.dumps(report, indent=2) + '\n')

    _plot(xs, ts, solution, residual, reference, operators, coefficients, length, alpha,
          count, output / 'solution.png')

    print(json.dumps(report, indent=2))
    print('\nSample numerical values:')
    print(' x       t=0          t=T/2        t=T')
    for i in [0, 30, 60, 90, len(xs) - 1]:
        if i >= len(xs):
            continue
        row = [solution[i, int(round(f * (len(ts) - 1)))] for f in (0.0, 0.5, 1.0)]
        print(f'{xs[i]:.2f}  ' + '  '.join(f'{v: .9f}' for v in row))
    print(f'\nOutputs: {output.resolve()}')
    if report['conditioning_warning']:
        print('NOTE: condition number * eps > 1; the system is ill-conditioned. '
              'Accuracy claims rest only on the sampled diagnostics above.')
    return report


def _plot(xs, ts, solution, residual, reference, operators, coefficients, length, alpha,
          count, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    columns = [0, (len(ts) - 1) // 2, len(ts) - 1]
    colors = ['#2563eb', '#b45309', '#15803d']
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    for j, color in zip(columns, colors):
        axes[0].plot(xs, solution[:, j], color=color, lw=2, label=f't = {ts[j]:.2f}')
        if reference is not None:
            axes[0].plot(xs[::6], reference[::6, j], 'o', ms=4, mfc='none', color=color)
            axes[1].semilogy(xs, np.maximum(np.abs(solution[:, j] - reference[:, j]), 1e-16),
                             color=color, label=f't = {ts[j]:.2f}')
        else:
            axes[1].semilogy(xs, np.maximum(np.abs(residual[:, j]), 1e-16),
                             color=color, label=f't = {ts[j]:.2f}')
    axes[0].set(title='Hermite solution (circles: exact)' if reference is not None else 'Hermite solution',
                xlabel='Position x', ylabel='u(x,t)')
    axes[1].set(title='Absolute solution error' if reference is not None else 'Absolute PDE residual',
                xlabel='Position x',
                ylabel='Absolute error (floor 1e-16)' if reference is not None else 'Absolute residual (floor 1e-16)')
    for ax in axes:
        ax.grid(alpha=.2)
        ax.legend()
    fig.suptitle(f'Heat equation | {count} x {count} Hermite coefficients | alpha={alpha:g}')
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _print_refinement(rows):
    print('basis  unknowns  solution_error  independent_pde_residual  algebraic_residual  cond(2)')
    for row in rows:
        print(f"{row['basis_count']:5d}  {row['unknown_coefficients']:8d}  "
              f"{row['max_abs_solution_error']:.3e}         "
              f"{row['validation_max_abs_pde_residual']:.3e}              "
              f"{row['collocation_max_abs_residual']:.2e}           "
              f"{row['condition_number_2norm']:.1e}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--basis', type=int, default=DEFAULT_BASIS,
                        help=f'Hermite functions per variable ({MIN_BASIS}-{MAX_BASIS})')
    parser.add_argument('--alpha', type=float, default=DEFAULT_ALPHA, help='Diffusivity (> 0)')
    parser.add_argument('--length', type=float, default=DEFAULT_LENGTH, help='Spatial length L (> 0)')
    parser.add_argument('--time', type=float, default=DEFAULT_TIME, help='Time horizon T (> 0)')
    parser.add_argument('--output', type=Path, default=Path('results-heat'))
    parser.add_argument('--no-reference', action='store_true',
                        help='Solve and validate without comparing to the exact solution.')
    parser.add_argument('--refine', action='store_true',
                        help='Print a basis-refinement table instead of writing outputs.')
    args = parser.parse_args()
    try:
        if args.refine:
            validate_parameters(args.basis, args.length, args.time, args.alpha)
            _print_refinement(refinement_sweep(args.length, args.time, args.alpha))
            raise SystemExit(0)
        result = run(args.basis, args.length, args.time, args.alpha, args.output,
                     compare_reference=not args.no_reference)
    except ValueError as exc:
        parser.exit(1, f'Error: {exc}\n')
    except RuntimeError as exc:
        parser.exit(1, f'Solve failed: {exc}\n')
    if not result['validation_passed']:
        parser.exit(2, 'Validation failed; outputs retained for diagnosis.\n')

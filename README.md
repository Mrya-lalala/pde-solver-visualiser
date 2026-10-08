# PDE Solver / Visualiser

A local Hermite-basis numerical solver and interactive visualiser for supported
PDE benchmarks. The interface supports Rosenau-Hyman and semi-spherical porous-fin
problems; the heat-equation benchmark is available as a separate script.

The local app now also includes the semi-spherical porous-fin PDE from Pavan
Kumar et al. (2026), implemented with Hermite polynomials. See [POROUS_FIN.md](POROUS_FIN.md)
for its equation, startup interpretation, tested dimensionless presets and
independent numerical validation. Select it in the problem menu and press Solve.
The heat solver remains separate from the local interface.

For academic review, see the [porous-fin equation comparison and implementation
note](output/pdf/Hermite_Porous_Fin_Comparison_Note.pdf). A saved porous-fin run
and offline viewer are included in [results-fin/](results-fin/).

The runnable programs use shifted physicists' Hermite polynomials and a shared viewer:

- `rosenau_hyman.py` - nonlinear Rosenau-Hyman benchmark.
- `heat_equation.py` - linear heat-equation benchmark (see its own section).
- `porous_fin.py` - transient semi-spherical porous-fin PDE from the supplied PDF.
- `build_viewer.py` - offline interactive HTML exporter for solver archives.

These solvers do not reproduce every table in the papers or accept arbitrary PDE
expressions. The Rosenau-Hyman benchmark solves

    u_t = u*u_xxx + u*u_x + 3*u_x*u_xx

on 0 <= x <= 1, 0 <= t <= 1. Current defaults: `SPEED = 1`, `SHIFT = 0`, `N = 6`
Hermite functions per variable. It implements the mixed-derivative
integration/collocation approach.

## Local solver interface

Run from this directory with the existing environment:

```powershell
.\.venv\Scripts\python.exe local_server.py
```

For a fresh download, create the environment and install both dependency lists first:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt -r requirements-visual.txt
```

Open **http://127.0.0.1:8765**. Stop the service with Ctrl+C in its terminal.
No new runtime dependencies, remote access, deployment or GPU are required.
The existing visual dependencies (`requirements-visual.txt`) supply Plotly.
The heat solver is unchanged and is not offered in the local interface.

The original preset is c=1, a=0, N=6, L=T=1. Press **Solve** to run the numerical
collocation solver. Editing inputs does not solve; time scrubbing, playback,
reference overlay and optional heatmap/error views use the completed arrays.
**Reset to original benchmark** restores the input preset and closes optional
views; press Solve to apply it. **Reset zoom** only resets the plots.

Supported discrete choices: c in {0.5, 1, 2}, a in {0, 0.5}, N in {5, 6},
T in {0.5, 1}; L is fixed at 1. Of these 24 combinations, 23 passed local
validation. The combination c=2, a=0.5, N=5, T=1 is rejected because its PDE
residual exceeds 1e-5. See `validated_combinations.json` for measurements.
These are tested combinations, not a guarantee for unrestricted numeric inputs.
The Rosenau-Hyman solver uses the single-cell Hermite basis. No arbitrary initial-function
entry or separate amplitude control is provided. Initial data and all three
boundary traces are generated together from speed and shift and shown in the UI.

Every Rosenau-Hyman solve checks finite fields, algebraic residual <=1e-8, independently sampled
PDE residual <=1e-5, and initial/boundary errors <=1e-12. Exact comparison is
optional and runs only after the numerical solve; it never supplies coefficients.
Sampled errors are not rigorous error bounds or evidence of general well-posedness.

One computation runs at a time. The displayed parameter label belongs to the
completed result, even when controls have been edited. Solver or validation
failure preserves the previous successful result. The latest successful result
remains available across page reloads for the lifetime of the service.
Each successful run saves NPZ, JSON diagnostics and an offline HTML snapshot in
`local-runs/<run-id>/`. Offline snapshots cannot start new solves.
Existing `results/`, `results-heat/` and `interactive_preview.html` are preserved.
The CLI's restored original benchmark now defaults to `results-original/`.

Verification commands:

```powershell
.\.venv\Scripts\python.exe -B -m unittest test_local_solver -v
node test_viewer.js
```

The optional Node test uses the already installed Node runtime and checks UI
state logic with test doubles; it does not test real-browser layout or Plotly
rendering. Automated checks cover all enabled combinations, reference independence,
original/speed-2 regression, request rejection, concurrent solves, failure retention
and offline data identity. Real-browser visual verification remains outstanding
because no connected browser was available during implementation.

## Setup and run

Open a terminal in this folder and install the dependencies once.

Windows (PowerShell):

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe rosenau_hyman.py
```

macOS/Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python rosenau_hyman.py
```

The script prints diagnostic results and sample solution values. Open
`results-original/solution.png` for the plotted curves and errors. No manual
precalculation or input is needed to run the included example.

## Input data

The labelled EDITABLE INPUTS section at the top of the script supplies:

- Initial profile: u(x,0) = -(8/3) cos(x/4)^2.
- Boundary value: u(0,t) = -(8/3) cos(t/4)^2.
- Boundary slope: u_x(0,t) = -(2/3) sin(t/2).
- Boundary curvature: u_xx(0,t) = (1/3) cos(t/2).
- Space/time interval lengths, number of basis functions and tolerances.
- An optional exact reference: -(8/3) cos((x-t)/4)^2.

These values restore the original `SPEED = 1`, `SHIFT = 0` benchmark.
The existing `results/` files retain the previous speed-2 run.
`interactive_preview.html` remains an older speed-1 viewer snapshot.

These derivative traces are additional explicit benchmark inputs, consistent
with the known exact solution. The paper's displayed example gives the initial
profile and boundary value; its reconstruction also needs derivative traces.
This implementation makes that extra information visible. It does not infer
spatial derivatives from the time-dependent boundary value, or establish
well-posedness for arbitrary choices of these traces.

The solve never receives the exact reference. Run without comparison using:

```bash
python rosenau_hyman.py --no-reference --output results-no-reference
```

Other options:

```bash
python rosenau_hyman.py --basis 6 --length 1 --time 0.5 --output results-short
python rosenau_hyman.py --help
```

N is the number of Hermite functions per variable, not the maximum degree:
N=6 uses degrees 0 through 5 and 36 coefficients. The interval starts at zero.
Wavelet subdivision is fixed at k=1; changing N changes polynomial resolution,
not the number of cells. Higher N is not automatically more reliable.

If changing the input functions, use SymPy expressions such as `sp.sin(x)`
and `sp.Rational(1, 3)`. All three boundary traces are required by this version.
If the exact reference no longer applies, replace it or set it to `None`.
The equation itself is fixed in `HermiteProblem.residual`; changing the PDE
also requires changing its algebraic Jacobian and checking the formulation.

## What is automated

1. Check initial/boundary compatibility at x=t=0.
2. Build shifted physicists' Hermite polynomials H_m(2x/L-1) and their time
   counterparts. Each is rescaled for conditioning; this changes coefficients
   but not the polynomial space.
3. Compute exact polynomial antiderivatives with lower endpoint zero.
4. Symbolically differentiate the supplied initial and boundary functions.
5. Approximate u_xxxt with a tensor-product Hermite expansion.
6. Reconstruct u and its derivatives by integration, incorporating the input
   functions as integration constants.
7. Enforce the PDE on an N-by-N tensor midpoint collocation grid.
8. Solve the nonlinear equations with SciPy's hybrid root method and an
   analytic Jacobian. No fitted coefficients are stored in the source.
9. Check the PDE on a separate dense grid, along with initial/boundary errors.
10. Plot and export numerical results; compare to the exact solution only if
    supplied.

The reconstruction is

    B(x,t) = f(x) + sum_{m=0}^2 x^m/m! * [g_m(t) - g_m(0)]
    u(x,t) = B(x,t) + I_x^3(phi_x)^T C I_t(phi_t)

where f is the initial profile and g_0, g_1, g_2 are the prescribed boundary
traces. The integrated polynomials are computed directly instead of manually
transcribing the paper's printed integration matrices.

## Output files

- `solution.png`: profiles at t=0, T/2 and T; exact error or PDE residual.
- `results.json`: input description, convergence message and error diagnostics.
- `solution.npz`: x, t, u, coefficient matrix, and optional exact values.
  `u[i,j]` is the solution at `x[i], t[j]`.

Read numerical arrays with:

```python
import numpy as np
r = np.load('results/solution.npz')
print(r['u'].shape)
print(r['u'][:, -1])  # final-time profile
```

## Validation performed

For the preserved previous run (SPEED = 2, SHIFT = 0, N = 6, unit square):

- Maximum sampled solution error: 3.42e-8.
- Maximum sampled PDE residual: 2.06e-6.
- Initial and supplied boundary trace errors: at most 8.9e-16.
- Disabling reference comparison produced identical numerical solution arrays.
- Refining N = 4 to N = 6 reduced the sampled solution error from 1.09e-5
  to 3.42e-8. At these defaults N = 4 fails validation (residual 7.1e-4);
  N = 5 (residual 9.4e-6) and N = 6 pass.
- Incompatible initial/boundary data was rejected.

Restored original benchmark, reproduced locally: at SPEED = 1, N = 6
the sampled solution error was about 1.99e-10 and the PDE residual about
1.98e-8. Those values apply to the restored speed-1 defaults.

These are sampled checks, not rigorous error bounds. If the independent PDE
residual exceeds its configured tolerance, the script saves diagnostic output
but returns exit code 2. An invalid input or failed algebraic solve returns 1.
A small algebraic residual can be accepted even if the optimizer stops for lack
of progress; that optimizer message is retained, and independent validation
still runs. More coefficients or different inputs can expose conditioning,
nonlinear convergence and well-posedness limitations.

## Heat-equation benchmark (`heat_equation.py`)

A separate, linear benchmark on `0 <= x <= L` and `0 <= t <= T`:

    u_t = alpha*u_xx,   u(x,0) = sin(pi*x/L),   u(0,t) = u(L,t) = 0
    u(x,t) = exp(-alpha*(pi/L)**2*t) * sin(pi*x/L)   (reference)

Defaults: `alpha = 1`, `L = 1`, `T = 1`, `N = 16` Hermite functions per
variable. Outputs go to `results-heat/`. This is deliberately separate from the
Rosenau-Hyman benchmark, whose original defaults are restored and now write to `results-original/`.

### What it supports

Only this fixed heat problem with homogeneous endpoints and the sine initial
profile. The diffusivity must be finite and strictly positive, as must L and T.
No reaction-diffusion, Burgers, general PDE parsing, other bases or 3D output
are included.

### How both endpoints are enforced

Every integral starts at the lower endpoint 0. With

    q(x,t) = u_xxt ~ phi_x(x)^T C phi_t(t)
    A_i(x) = I_x^2(phi_i)(x) - (x/L) * I_x^2(phi_i)(L)
    B_j(t) = I_t(phi_j)(t)
    u(x,t) = f(x) + A(x)^T C B(t)

`B_j(0) = 0`, so `u(x,0) = f(x)` exactly, and `A_i(0) = A_i(L) = 0`, so
`u(0,t) = f(0) = 0` and `u(L,t) = f(L) = sin(pi) = 0` for every t. Both
homogeneous endpoints are therefore satisfied without prescribing boundary
slopes (unlike the Rosenau-Hyman benchmark, which needs u_x and u_xx at x=0).
Differentiating gives `u_t = A^T C phi_t` and `u_xx = f'' + phi_x^T C B`.
Collocating `u_t - alpha*u_xx = 0` at the tensor product of the N interval
midpoints gives one linear system, solved directly with `numpy.linalg.solve`:

    [ A(x) (x) phi_t(t) - alpha * phi_x(x) (x) B(t) ] c = alpha * f''(x)

The reference never enters assembly or the solve.

### Run, validate and export the viewer

Windows PowerShell:

```powershell
.\.venv\Scripts\python.exe heat_equation.py --output results-heat
.\.venv\Scripts\python.exe heat_equation.py --refine
.\.venv\Scripts\python.exe build_viewer.py results-heat/solution.npz --output results-heat/interactive.html
```

Options: `--basis`, `--alpha`, `--length`, `--time`, `--output`, `--no-reference`,
`--refine`. A solve or validation failure returns exit code 1 or 2; a failing
run still writes its diagnostics so they can be inspected.

### What was validated

- Default unit square (N = 16): sampled solution error 4.9e-8 (<= 1e-6),
  independent PDE residual 5.7e-6 (<= 1e-5), initial and endpoint errors
  <= 1.2e-16 (<= 1e-12).
- Basis refinement (N = 8, 10, 12, 14, 16, 18, 20): solution error and PDE
  residual fall monotonically from 2.3e-3 / 1.3e-1 to 2.2e-9 / 3.1e-7 at
  N = 18, then stop improving at N = 20 as round-off begins to dominate.
- A scaling case with L != T and alpha != 1 (alpha = 0.4, L = 0.7, T = 1.3)
  gives solution error 1.0e-7 and PDE residual 9.2e-6.
- Enabling and disabling the reference comparison produced byte-identical
  `x`, `t`, `u` and coefficient arrays.
- Invalid parameters (non-positive or non-finite alpha, L, T, or an
  out-of-range basis count) are rejected with exit code 1.

### Why this is not an exact solution

`u` is an explicit polynomial approximation, not the exact solution; it only
matches the reference within the sampled tolerances above. Those are sampled
diagnostics, not rigorous error bounds, and small sampled errors do not imply
support for arbitrary PDEs, initial data or boundary conditions.

Numerical limitation: the time factor is a polynomial of degree N, so accuracy
degrades when `alpha*(pi/L)**2*T` is large because the fast exponential decay is
hard to represent. For example alpha = 2, L = 0.5, T = 1.2 gives a sampled
solution error near 1e-1 at N = 16, and alpha = 4, L = T = 1 gives about 1e-2.
The reported 2-norm condition number is very large (about 1e19 at N = 16) and
grows with N; it is reported separately from the observed accuracy, which is
always checked on an independent grid and should not be trusted on
conditioning alone.

## Interactive viewer (`build_viewer.py`)

`build_viewer.py` turns any saved solver archive into a self-contained offline
HTML file. It never runs a solver and never modifies the archive.

- `results/interactive.html` matches the preserved previous speed-2 Rosenau-Hyman archive.
- `results-heat/interactive.html` matches the heat-equation archive.
- `interactive_preview.html` is a historical preview generated from an older
  **speed-1** Rosenau-Hyman archive; it uses the original parameters but the older interface.

Regenerate a viewer after any new solver run:

```powershell
.\.venv\Scripts\python.exe build_viewer.py <archive.npz> --output <output.html>
```

## Reference

S. Kumbinarasaiah and Waleed Adel, "Hermite wavelet method for solving nonlinear
Rosenau-Hyman equation", Partial Differential Equations in Applied Mathematics
4 (2021), 100062. DOI: 10.1016/j.padiff.2021.100062.

Only documented problem formulations are implemented. The porous-fin solver uses
multiple spatial cells with Hermite polynomials; the other benchmarks use a single
cell. No general basis selection, boundary-condition engine or PDE classifier is included.

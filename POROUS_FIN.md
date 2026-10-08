# Semi-spherical porous fin: Hermite implementation

Source: P. L. Pavan Kumar, B. J. Gireesha and P. Venkatesh,
"Ceramic-based semi-spherical porous fin under condensation condition: a transient
study on thermal performance and efficiency", Int J Mech Mater Des 22:72 (2026),
DOI 10.1007/s10999-026-09893-6. The supplied PDF was read and the displayed
equations checked visually. This implements the dimensionless Eq. (13) and (14).
The paper used finite differences; this implementation uses Hermite polynomials.

## Equation and data

With A(X)=1-X^2 and S(X)=sqrt(A),

    A theta_tau = (A theta_X)_X - S F(theta), 0<X<1
    F(theta) = Nc theta^2 + Nr theta + n1_squared theta^(p+1) + wet_term theta^p
    n1_squared = n1^2, wet_term = n0^2*n3
    theta(X,0)=0 for X>0
    theta(0,tau)=1 for tau>=0 (base step)
    theta_X(1,tau)=0

The signed moisture coefficient is retained. It is not replaced by its absolute
value. Fractional exponents are not exposed in this initial implementation.

## Formulation

Use theta=1+sum_i c_i(tau) v_i(X). Construct v_i from shifted, scaled physicists'
Hermite polynomials on cells with cosine-clustered endpoints. Linear constraints
enforce continuous values and first derivatives between cells, v_i(0)=0,
and v_i'(1)=0. Six local functions use polynomial degrees zero through five;
12 cells yield 48 free coefficients after the constraints.

Multiplying the conservative PDE by a test function and integrating gives

    M c' + K c + integral_0^1 S v F(1+V c) dX = 0
    M_ij = integral_0^1 A v_i v_j dX
    K_ij = integral_0^1 A v_i' v_j' dX

The boundary term vanishes. No division by A or evaluation of the PDE at the
degenerate tip is required. A Cholesky change of basis makes the mass matrix
identity; the Hermite approximation space is preserved. Gauss quadrature uses
32 nodes per cell. SciPy's implicit BDF integrator solves the resulting nonlinear
ODE with an analytic Jacobian. This is spatial Hermite approximation with numerical
time integration, rather than the Rosenau-Hyman tensor space-time reconstruction.

## Startup and source limitations

Eq. (14) has a discontinuous starting corner: zero initial interior data and
base value one. Interpret it as switching the base on at time zero. The initial
interior state is projected in weighted L2, which cannot impose the discontinuous
profile pointwise. The reported initial projection error is explicit. The first
saved frame shows the prescribed initial data with base value one; the coefficient
representation is the projection. Later frames are evaluated from the numerical
solution. Startup is not certified for tau<0.05.

The paper's dimensional Eq. (11), T(x,0)=0, does not transform to the zero theta
initial condition under its own definition of theta. We follow dimensionless
Eq. (14). Table 3 reports theta at X=0 below one, conflicting with its base
condition. Eq. (10)'s printed radiation term omits the Ta^3 factor used in
Eq. (12). These prevent treating the printed tables as unquestioned benchmarks.

The presets below supply explicit dimensionless coefficients. They are not
calibrated SiC or Si3N4 cases and do not claim reproduction of the paper's figures.
Material calibration needs clarification of the dimensional inputs, latent heat,
and the intended initial/boundary data. No condensation switch at the dew point
is added: this implements the printed dimensionless equation.

## Run

Use the existing project environment and dependencies:

```powershell
.\.venv\Scripts\python.exe local_server.py
```

Open http://127.0.0.1:8765, select **Semi-spherical porous fin (PDF Eq. 13)**,
choose a preset/resolution and press **Solve**. Start the server again if an old
process is still serving the previous interface. On macOS activate `.venv` and
run `python local_server.py` as before.

Standalone example:

```powershell
.\.venv\Scripts\python.exe porous_fin.py
.\.venv\Scripts\python.exe porous_fin.py --cells 16 --no-reference --output results-fin-no-comparison
```

The default is Nc=0.5, Nr=0.1, n1^2=0.5, n0^2*n3=0.1, p=1, tau_max=0.8,
12 cells and six Hermite functions per cell. It writes `results-fin/solution.npz`,
`results.json`, `solution.png`, and the offline `interactive.html` viewer.

The CLI additionally accepts `--nc`, `--nr`, `--n1-squared`, `--wet-term`,
`--power`, `--duration`, `--cells`, `--count`, `--output` and `--no-reference`.
CLI arithmetic bounds permit exploratory coefficients; these are not all tested
combinations. Solve failure returns exit 1; validation failure returns exit 2
with diagnostic files retained. The local app accepts only the registry below
and preserves the previous successful result after failure.

## Tested combinations and diagnostics

The local app exposes four dimensionless presets:

| Preset | Nc | Nr | n1^2 | n0^2*n3 | p |
|---|---:|---:|---:|---:|---:|
| Illustrative losses | 0.5 | 0.1 | 0.5 | 0.1 | 1 |
| Zero-loss diffusion check | 0 | 0 | 0 | 0 | 1 |
| Signed moisture term | 0.5 | 0.1 | 0.5 | -0.1 | 1 |
| Power p=2 | 0.5 | 0.1 | 0.5 | 0.1 | 2 |

For each preset, 8/12/16 cells, six Hermite functions per cell, and durations
0.4/0.8 passed: 24 enabled combinations. Four-function cases were also tested;
their strong residuals failed and they are not offered in the app. Measurements
are recorded in `validated_fin_combinations.json`.

Acceptance requires finite data, base/tip errors <=1e-10, weak residual <=1e-5,
weighted strong residual <=0.02, and sampled temperature in [-0.002,1.002] after
startup. The weak residual uses independent 48-node quadrature. The strong
residual samples nine interior locations per cell at 21 times between 0.05 and
the final time. Temporal derivatives are estimated from the completed BDF
interpolant with a separate second-order finite difference, rather than substituting
the ODE right-hand side. Reported strong residuals use the undivided PDE A theta_tau;
they are not residuals of a singular divided equation.

Default measurements: comparison difference about 3.6e-6, weak residual 1.1e-6,
weighted strong residual 4.6e-3. Strong residual falls from about 6.8e-3 to 4.6e-3
to 3.4e-3 for 8, 12, 16 cells. These sampled diagnostics have different scales
from the existing Rosenau-Hyman checks and are not rigorous error bounds.

Optional comparison independently solves the conservative equation with 400
finite-volume cells and BDF time integration. With comparison enabled, maximum
difference after startup must be <=0.002. Tests also compare against 800 cells
and check a manufactured solution with analytically derived forcing. Turning
comparison off gives identical Hermite arrays. The PDF supplies no exact solution;
no numerical reference is labelled exact.

Eq. (16) efficiency is computed by quadrature and saved as a time series in NPZ
and a final value in JSON. It is undefined for zero-loss cases, which are marked
explicitly. The viewer overlays the optional independent comparison and labels
its difference view accordingly. `exact` retains its existing meaning for genuine
exact references; this solver uses a separate optional `comparison` array.

```powershell
.\.venv\Scripts\python.exe -B -m unittest test_local_solver test_porous_fin -v
node test_viewer.js
```

Numerical, API, export and UI-state tests passed. No browser was connected for
visual testing of the updated interface. PDF parsing used PyMuPDF only during
development; running the solvers does not require it or access to the original PDF.

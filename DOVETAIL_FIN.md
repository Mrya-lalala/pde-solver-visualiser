# Inclined dovetail fin with Legendre wavelets

Source: P. L. Pavan Kumar, B. J. Gireesha and P. Venkatesh,
"Impact of trihybrid nanofluid on the transient thermal performance of inclined
dovetail fin with emphasis on internal heat generation", European Physical
Journal Plus 139:60 (2024), DOI
[10.1140/epjp/s13360-023-04848-8](https://doi.org/10.1140/epjp/s13360-023-04848-8).
The supplied paper's governing equation is Eq. (14), with data in Eq. (15).
Its Eq. (13) defines dimensionless quantities. This is a different problem from
the 2026 semi-spherical fin already in the app.

## Equation implemented

Define the material-property ratios

    H = (rho Cp)_thnf / (rho Cp)_f
    Rrho = rho_thnf / rho_f
    Rbeta = (rho beta)_thnf / (rho beta)_f
    Rmu = mu_thnf / mu_f
    Rk = k_thnf / k_f

With `a=1-CX` and `s=theta-theta_a`, the printed equation is

    a theta_tau = D (a theta_X)_X + g a (1+epsilon_G s)
                  - B s^2 - R s - W s^(p+1)
    D = Rk/H
    g = Q/H
    B = (Rbeta/Rmu) Nc sin(alpha)
    R = Nr/H
    W = (M^2/H) (1+n1 Rrho/H) / (1-theta_a)^p

Expanding `(a theta_X)_X` gives `a theta_XX-C theta_X`, reproducing
the two diffusion terms of Eq. (14). All printed signs and integer powers are
retained; no absolute-value substitution, clipping or condensation switch is used.

The boundary orientation differs from the earlier semi-spherical problem:

    theta_X(0,tau)=0        insulated tip
    theta(1,tau)=1          fixed-temperature base
    theta(X,0)=0 for X<1   initial interior

Material ratios use Table 1's 50/50 ethylene glycol-water base fluid and MWCNT,
Ag and Cu properties. The mixture laws follow the normalized right-hand sides
of Table 2. The baseline has phi1=phi2=phi3=0.02. Linear verification uses zero
particle fractions, reducing all property ratios to one.

## Legendre-wavelet method

For `J=2^(k-1)` dyadic cells of width `h=1/J`, use the conventional localized basis

    psi_(n,m)(X) = sqrt((2m+1)/h) P_m(2(X-nh)/h-1)
                   on nh <= X < (n+1)h; zero elsewhere

Here `n=0,...,J-1` and `m=0,...,M-1`. Six functions mean Legendre degrees 0-5.
These functions are orthonormal before constraints and have compact cell support.
This is a localized Legendre-wavelet expansion, rather than one global Legendre
polynomial. See, for example, the basis definition in
[Legendre wavelets approaches](https://www.sciencedirect.com/science/article/pii/S0307904X19300319).

Write `theta=1+V c(tau)`. Linear constraints on combinations of the wavelets
enforce continuous values/slopes at cell interfaces, zero slope at the tip and
zero basis value at the base. A mass-matrix change of coordinates improves
conditioning without changing the approximation space. The resulting weak form is

    Mass c' + D Stiffness c + integral v [losses(theta)
        - g a (1+epsilon_G (theta-theta_a))] dX = 0
    Mass_ij = integral a v_i v_j dX
    Stiffness_ij = integral a v_i' v_j' dX

Gauss quadrature uses 24 nodes per cell. SciPy's implicit BDF method advances time
with an analytic Jacobian, rtol=1e-9 and atol=1e-11. Only space uses wavelets;
this is not a tensor space-time wavelet discretization. Eight cells and six local
functions produce 32 free coefficients after constraints.

The prescribed zero interior and base-one step are incompatible at the starting
corner. Initial coefficients come from a weighted L2 projection. Its default
error is about 0.0636. The saved time-zero display frame explicitly shows the
prescribed data; it is not the polynomial reconstruction of the initial
coefficients. Later frames are numerical reconstructions. Startup before tau=0.05
is not certified.

## Interpretations needed for the paper

- Eq. (13)'s printed time scaling omits the factor `1/L^2` normally needed to
  make diffusivity times dimensional time dimensionless. This code solves
  Eq. (14) in its stated time coordinate; it does not assign that coordinate a
  physical value in seconds or silently repair the dimensional derivation.
- Table 2 labels several left-hand sides as dimensional properties, while their
  right-hand sides are dimensionless ratios. The implementation explicitly uses
  those right-hand sides as ratios. Density/heat/expansion use the printed nesting;
  conductivity follows the printed successive mixture laws.
- The results section gives `G=0.1`, while the PDE uses `Q`. The baseline preset
  explicitly interprets this as `Q=0.1`. It uses M^2=1, Nc=Nr=n1=1, p=1,
  theta_a=0.2, epsilon_G=0.2, C=-0.1 and alpha=45 degrees.
- The early dimensional velocity equation uses sin(tau), but later governing
  equations use the inclination sin(alpha). The solver follows Eq. (14).

These choices define a computable dimensionless problem. They do not establish
that all dimensional equations and plotted material cases in the paper are
mutually consistent. The supplied initial temperature is retained as printed,
rather than replaced by an ambient-temperature initial state.

## Supported interface controls

Select **Inclined dovetail fin (2024, Eq. 14)**. Its implemented basis is
**Legendre wavelets**. The earlier two interface problems continue to use Hermite.
There is no unsupported basis switch or arbitrary PDE/function entry.

Four coefficient presets are offered:

| Preset | Change from the baseline |
| --- | --- |
| Paper baseline | Values stated above, with the explicit Q interpretation |
| Wider dovetail | C=-0.3 |
| No internal heat generation | Q=0 |
| Linear verification | Pure fluid; C=Nc=Nr=Q=n1=theta_a=epsilon_G=0; M^2=1; p=0 |

For the three nonlinear presets, 4/8/16 cells with six functions and durations
0.4/0.8 passed: 18 combinations. The linear benchmark uses 8/16 cells, six
functions and duration 5: two combinations. Other combinations are rejected.
Four-function candidates and the coarsest linear case failed the same thresholds
and are not offered as validated options. All 42 measured candidates, including
failures, are recorded in `validated_dovetail_combinations.json`.

Press **Solve** to run a new numerical calculation. Editing controls does not
solve. Time scrubbing uses completed arrays. The optional comparison is an
independent 400-cell finite-volume calculation, never an exact-formula substitute
for the wavelet solve. On the deployed site the same Python runs in a browser
worker; failed solves retain the last successful displayed result.

## Validation

Default nonlinear case, after tau=0.05:

| Diagnostic | Value |
| --- | --- |
| Maximum difference from independent 400-cell finite-volume result | 2.84e-6 |
| Independent weak residual | 5.00e-8 |
| Sampled weighted strong PDE residual | 1.47e-4 |
| Base-value error | 0 (desktop run) |
| Insulated-tip derivative error | 1.93e-14 |

Acceptance thresholds are weak residual <=1e-5, sampled weighted strong residual
<=2e-3, comparison difference <=2e-4 and boundary errors <=1e-10. The sampled
strong residual decreases from 1.77e-3 to 1.47e-4 to 9.90e-6 for 4/8/16 cells.
An additional 800-cell reference also agrees within 2e-6.

Tests verify raw-wavelet orthonormality/local support, interface and boundary
constraints, a manufactured nonlinear solution with independently derived forcing,
reference independence, invalid inputs, API solving and exported arrays.

For the linear limit `theta_tau=theta_XX-theta`, an independent analytic transient
reference is

    theta = cosh(X)/cosh(1)
            - sum_(n>=0) [2 lambda_n (-1)^n/(lambda_n^2+1)]
              cos(lambda_n X) exp(-(lambda_n^2+1)tau)
    lambda_n=(n+1/2) pi

The default eight-cell solution differs from it by at most 2.33e-8 after startup.
At tau=5, its samples agree with Table 3's "Present work" column within 8e-7.
This is a limiting verification case; it is not an exact solution of the nonlinear
trihybrid problem. No exact array is attached to that problem in the interface.

## Run and output

```powershell
.\.venv\Scripts\python.exe local_server.py
.\.venv\Scripts\python.exe dovetail_fin.py
.\.venv\Scripts\python.exe dovetail_fin.py --preset wider_dovetail --cells 16 --time 0.4
.\.venv\Scripts\python.exe -m unittest test_dovetail_fin
```

No new runtime dependencies are needed. The CLI also accepts `--no-reference` and
`--output`. Its arithmetic bounds allow exploration; the service additionally
requires an exact match to the validated registry. Validation failure saves
diagnostics and returns exit code 2. A computation/input failure returns 1.

`results-dovetail/` contains `solution.npz`, `results.json`, `solution.png` and the
self-contained offline viewer. The existing `x,t,u` array contract is preserved.
`coefficients` are mass-normalized constrained coordinates; `wavelet_coefficients`
are coefficients of the raw localized basis in the formula above, and `cell_edges`
specify its cells. Reconstruction is `1+psi @ wavelet_coefficients`, except for
the explicitly prescribed saved time-zero frame.

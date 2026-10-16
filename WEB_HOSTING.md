# Shareable browser solver

Published URL: https://mrya-lalala.github.io/pde-solver-visualiser/

The public GitHub Pages site serves `docs/` from `main`. Visitors do not need
Python, a GitHub account, or access to the developer's computer. The initial
porous-fin plot is an explicitly labelled saved numerical sample. Pressing
**Solve** executes a new numerical computation on the visitor's device.

## How it works

- `build_web.py` uses the existing viewer template and Plotly exporter.
- `browser_client.js` maps the existing job/status/result contract to a Web Worker.
  Jobs and results belong to that browser tab, independently of other visitors.
- `browser_worker.js` loads pinned Pyodide 0.28.3 and its NumPy, SciPy and SymPy
  packages from jsDelivr. It imports copies of the original Python modules and
  calls the shared `solver_service.solve_request` entry point.
- `solver_service.py` contains the same structured-input validation, supported
  combinations and solver dispatch used by `local_server.py`. The Rosenau-Hyman
  and porous-fin numerical formulations and tolerances are unchanged.
- Exact/reference comparisons remain optional and separate from the numerical
  solve. Porous-fin comparison uses the independent finite-volume solver, not an
  exact solution. All results run the existing independent validation checks.

No Python backend runs on GitHub Pages. Python executes as WebAssembly in a
worker, keeping the interface responsive. No solver inputs are posted to a shared
service. Invalid inputs are rejected; failed or timed-out jobs retain the last
successful result. Reloading the page returns to the included sample.

The first Solve needs an internet connection to download the runtime and numerical
libraries. Later solves reuse the runtime while the tab stays open. Performance
depends on the device and connection; use a current desktop browser for review.
The included numerical sample, optional heatmap, time controls and comparison PDF
are available without starting a new solve.

## Update the deployment

After modifying the solver, template, supported registries or browser adapter:

```powershell
.\.venv\Scripts\python.exe -m unittest test_local_solver test_porous_fin
node test_viewer.js
node test_browser_client.js
.\.venv\Scripts\python.exe build_web.py
```

Commit the source changes and rebuilt `docs/`, then push to `main`. GitHub Pages
automatically publishes the `docs/` folder. On macOS/Linux use the activated
environment's `python` in place of the Windows Python path.

For a local browser preview, serve `docs/` over HTTP:

```powershell
.\.venv\Scripts\python.exe -m http.server 8890 --bind 127.0.0.1 --directory docs
```

Open http://127.0.0.1:8890. Opening `docs/index.html` directly as a `file:` URL does
not support the numerical worker. Existing self-contained offline exports and
the loopback-only Python application continue to work as before.

## Verification performed on 16 October 2026

- The 10 existing native numerical/API tests pass after extracting shared dispatch.
- Viewer state checks and browser-adapter checks pass, including invalid inputs,
  visitor isolation, stale messages, failure retention and timeout recovery.
- Headless Chrome executes real browser solves for all 47 exposed combinations:
  24 porous-fin and 23 Rosenau-Hyman cases. All pass their existing numerical
  validation thresholds under Pyodide 0.28.3.
- The default browser porous-fin result has a maximum difference of about
  `3.58e-6` from the independent 400-cell finite-volume calculation, consistent
  with the desktop result. The default Rosenau-Hyman independent PDE residual is
  about `1.98e-8`.
- Browser checks exercise both problem families, parameter changes and immediate
  time scrubbing; desktop rendering is visually inspected.

These checks cover the documented supported combinations, not arbitrary PDEs or
unrestricted inputs. The porous-fin coefficients remain illustrative dimensionless
presets; the paper's unresolved formulation questions are documented in the linked
comparison PDF.

References: [Pyodide Web Workers](https://pyodide.org/en/0.28.3/usage/webworker.html)
and [GitHub Pages publishing sources](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site).

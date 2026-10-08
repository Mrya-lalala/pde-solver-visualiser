# Local solving and offline export

For the local Hermite PDE interface, run:

```powershell
.\.venv\Scripts\python.exe local_server.py
```

Open http://127.0.0.1:8765, select Rosenau-Hyman or the semi-spherical porous-fin
problem and press Solve. Each uses its tested combinations; see POROUS_FIN.md
for the PDF equation, coefficient presets and numerical comparison. The heat solver remains separate.
The service binds only to loopback, validates numeric JSON and request origin,
and runs a single numerical job at a time. No solver runs on slider movement.
Successful runs are saved separately under `local-runs/<run-id>/`, including an
offline HTML export and the parameter/validation report. Failed jobs retain the
last successful displayed result. Stop with Ctrl+C.

The exporter below still accepts the existing NPZ contract. Its updated template
puts the solution curve first, with optional heatmap/error panels. Saved exports
have no live Python connection. Existing historical HTML files are unchanged.

# Interactive viewer add-on (solver unchanged)

Copy these three files into your existing solver folder, alongside
`rosenau_hyman.py`:

- `build_viewer.py`
- `viewer_template.html` (keep beside build_viewer.py)
- `requirements-visual.txt`

Use the Python environment you already use for the solver. From that folder:

```bash
python -m pip install -r requirements-visual.txt
python build_viewer.py results/solution.npz --output results/interactive.html
```

Open `results/interactive.html` in Edge, Chrome or Firefox. On Windows you can
simply double-click it. In WSL, open the project folder through Windows File
Explorer, then double-click the HTML. WSL or a special GPU setup is not required.

No changes to the existing solver are required. No need to rerun it if
`results/solution.npz` already exists. This exporter never imports the solver,
never solves a PDE and never modifies the input archive. After a future solver
run produces new results, run the exporter again. Existing HTML is a snapshot,
not a live connection to Python.

## Features

- Navy/violet/cyan/yellow gradient with fine contours.
- Position x horizontally; time t vertically.
- Fixed colour scale and profile y range across time for honest comparison.
- Selected-time line linked to the solution profile.
- Play/pause, speed control, time slider, and click-to-select saved samples.
- Native plot zoom, pan, reset and PNG export in the plot toolbar.
- Absolute-error landscape if the archive includes an exact reference.
- Fully embedded Plotly and data: the generated HTML works offline.
- No server, Node.js, CUDA or SwapTube installation required.

The included `interactive_preview.html` uses the previously generated benchmark
output, without recomputing it. Open that immediately to see the design.

## Optional automatic integration

For now, running the exporter as a second command is simplest. If you later
want one command to produce both numerical and interactive results, add this
inside the existing solver's `run` function AFTER it writes `solution.npz`:

```python
from build_viewer import build_viewer
build_viewer(output / 'solution.npz', output / 'interactive.html')
```

This assumes `output` is a pathlib.Path, as in the delivered solver. Keep both
viewer source files beside the solver. This hook is optional; the add-on does
not edit your source automatically.

## Data contract and interpretation

Input must contain finite numeric arrays:

- x: one-dimensional, strictly increasing, at least two entries.
- t: one-dimensional, strictly increasing, at least two entries.
- u: shape `(len(x), len(t))`.
- exact: optional, same shape as u.

Colours represent the saved computed values. Contours interpolate the samples;
zooming does not add precision. Slider playback visits saved time samples; it
does not recompute intermediate solutions. Clicking selects the nearest saved
sample; hover values can reflect the plotting library's contour interpolation.
The lower plot always shows the solution, even while the upper plot shows error.

Absolute error is `abs(u - exact)`. Its colours use log10, with a labelled
positive display floor to include zeros. Reported point errors remain the
unfloored values. If exact is absent, error mode is disabled. A small error in
this benchmark is not evidence that arbitrary PDEs or conditions are supported.

This viewer deliberately does not invent particle paths, physical flow or
fractal structure. The benchmark produces smooth bands. It is a 2D interactive
viewer; a rotatable 3D surface is not part of this add-on.

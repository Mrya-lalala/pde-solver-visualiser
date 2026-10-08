"""Export an existing solver NPZ as a self-contained interactive HTML viewer.
No solver import or PDE computation. Usage: python build_viewer.py results/solution.npz
"""
import argparse
import json
from pathlib import Path
import numpy as np
from plotly.offline import get_plotlyjs


def render_viewer(payload):
    template = Path(__file__).with_name('viewer_template.html').read_text(encoding='utf-8')
    return template.replace('/*__PLOTLY__*/', get_plotlyjs()).replace(
        '/*__DATA__*/', json.dumps(payload, allow_nan=False).replace('<', '\\u003c'))


def viewer_payload(arrays, report=None, **extra):
    return dict(x=arrays['x'].tolist(), t=arrays['t'].tolist(), u=arrays['u'].T.tolist(),
                exact=None if 'exact' not in arrays else arrays['exact'].T.tolist(),
                comparison=None if 'comparison' not in arrays else arrays['comparison'].T.tolist(),
                report=report, **extra)


def build_viewer(source, output):
    source, output = Path(source), Path(output)
    if source.resolve() == output.resolve():
        raise ValueError('Output must not overwrite input.')
    with np.load(source, allow_pickle=False) as data:
        x, t, u = (np.asarray(data[k], dtype=float) for k in ('x', 't', 'u'))
        exact = np.asarray(data['exact'], dtype=float) if 'exact' in data else None
        comparison = np.asarray(data['comparison'], dtype=float) if 'comparison' in data else None
    for name, a in [('x', x), ('t', t)]:
        if a.ndim != 1 or len(a) < 2 or not np.isfinite(a).all() or not (np.diff(a) > 0).all():
            raise ValueError(f'{name} must be a finite, strictly increasing 1D array with at least 2 values.')
    if u.shape != (len(x), len(t)) or not np.isfinite(u).all():
        raise ValueError('u must be finite and have shape (len(x), len(t)).')
    if exact is not None and (exact.shape != u.shape or not np.isfinite(exact).all()):
        raise ValueError('exact must be finite and have the same shape as u.')
    if comparison is not None and (comparison.shape != u.shape or not np.isfinite(comparison).all()):
        raise ValueError('comparison must be finite and have the same shape as u.')
    # Rows of browser field correspond to times, columns to positions.
    payload = dict(x=x.tolist(), t=t.tolist(), u=u.T.tolist(),
                   exact=None if exact is None else exact.T.tolist(),
                   comparison=None if comparison is None else comparison.T.tolist(), source=source.name)
    report_path = source.with_name('results.json')
    payload['report'] = json.loads(report_path.read_text()) if report_path.exists() else None
    text = render_viewer(payload)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(text, encoding='utf-8')
    print(f'Viewer created: {output.resolve()}\nOpen this HTML in your browser. No server or internet required.')
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', nargs='?', default='results/solution.npz')
    parser.add_argument('--output', default='results/interactive.html')
    args = parser.parse_args()
    try:
        build_viewer(args.input, args.output)
    except (OSError, ValueError, KeyError) as e:
        parser.exit(1, f'Error: {e}\n')

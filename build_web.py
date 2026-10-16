"""Build the shareable browser solver for static hosting (no backend account)."""
import argparse
import json
from pathlib import Path
import shutil

import numpy as np
from build_viewer import render_viewer, viewer_payload
from solver_service import DEFAULTS, SUPPORTED, FIN_SUPPORTED, FIN_PRESETS

ROOT = Path(__file__).resolve().parent


def build(output=ROOT / 'docs'):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    # Saved numerical run is available immediately, before any runtime download.
    with np.load(ROOT/'results-fin/solution.npz', allow_pickle=False) as saved:
        arrays = {key:saved[key] for key in saved.files}
    report = json.loads((ROOT/'results-fin/results.json').read_text())
    report['run_id'] = 'included-porous-fin-sample'
    report['execution_environment'] = 'Included numerical run; press Solve for a new browser computation'
    initial = dict(live=True, browser=True, token='browser-worker', defaults=DEFAULTS,
                   supported=SUPPORTED, fin_supported=FIN_SUPPORTED, fin_presets=FIN_PRESETS,
                   saved_result=viewer_payload(arrays, report))
    html = render_viewer(initial).replace('</main><script>',
        '</main><script src="browser_client.js"></script><script>', 1)
    # Keep the professor-facing mathematical note one click away.
    html = html.replace('<section id="editor" hidden>',
        '<p class="hint">Solves run on your device. The first Solve downloads Python and numerical libraries; '
        'later solves reuse them while this tab stays open. '
        '<a href="comparison-note.pdf" style="color:#82dddd">Equation comparison PDF</a> · '
        '<a href="https://github.com/Mrya-lalala/pde-solver-visualiser" style="color:#82dddd">Source and documentation</a></p>'
        '<section id="editor" hidden>', 1)
    (output/'index.html').write_text(html, encoding='utf-8')
    (output/'.nojekyll').touch()
    for name in ('browser_client.js', 'browser_worker.js'):
        shutil.copyfile(ROOT/name, output/name)
    python = output/'python'
    python.mkdir(exist_ok=True)
    for name in ('rosenau_hyman.py', 'porous_fin.py', 'solver_service.py',
                 'validated_combinations.json', 'validated_fin_combinations.json'):
        shutil.copyfile(ROOT/name, python/name)
    shutil.copyfile(ROOT/'output/pdf/Hermite_Porous_Fin_Comparison_Note.pdf', output/'comparison-note.pdf')
    print(f'Browser solver built: {output.resolve()}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'docs')
    build(parser.parse_args().output)

"""Build the shareable browser solver for static hosting (no backend account)."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from build_viewer import render_viewer, viewer_payload
from solver_service import (DEFAULTS,SUPPORTED,FIN_SUPPORTED,FIN_PRESETS,
                            DOVETAIL_SUPPORTED,DOVETAIL_PRESETS,DOVETAIL_MODELS)

ROOT = Path(__file__).resolve().parent


def build(output=ROOT / 'docs'):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    # The newest numerical sample is immediately available, before runtime download.
    with np.load(ROOT/'results-dovetail/solution.npz', allow_pickle=False) as saved:
        arrays = {key:saved[key] for key in saved.files}
    report = json.loads((ROOT/'results-dovetail/results.json').read_text())
    report['run_id'] = 'included-dovetail-fin-sample'
    report['execution_environment'] = 'Included numerical run; press Solve for a new browser computation'
    assets=('rosenau_hyman.py','porous_fin.py','dovetail_fin.py','solver_service.py',
            'validated_combinations.json','validated_fin_combinations.json','validated_dovetail_combinations.json',
            'browser_client.js','browser_worker.js','viewer_template.html')
    version=hashlib.sha256(b''.join((ROOT/name).read_bytes() for name in assets)).hexdigest()[:12]
    initial = dict(live=True, browser=True, token='browser-worker',site_version=version,defaults=DEFAULTS,
                   supported=SUPPORTED, fin_supported=FIN_SUPPORTED, fin_presets=FIN_PRESETS,
                   dovetail_supported=DOVETAIL_SUPPORTED,dovetail_presets=DOVETAIL_PRESETS,
                   dovetail_models=DOVETAIL_MODELS,
                   saved_result=viewer_payload(arrays, report))
    html = render_viewer(initial).replace('</main><script>',
        f'</main><script src="browser_client.js?v={version}"></script><script>', 1)
    # Keep the professor-facing mathematical note one click away.
    html = html.replace('<section id="editor" hidden>',
        '<p class="hint">Solves run on your device. The first Solve downloads Python and numerical libraries; '
        'later solves reuse them while this tab stays open. '
        '<a href="comparison-note.pdf" style="color:#82dddd">Earlier semi-spherical fin comparison PDF</a> · '
        '<a href="https://github.com/Mrya-lalala/pde-solver-visualiser" style="color:#82dddd">Source and documentation</a></p>'
        '<section id="editor" hidden>', 1)
    (output/'index.html').write_text(html, encoding='utf-8')
    (output/'.nojekyll').touch()
    for name in ('browser_client.js', 'browser_worker.js'):
        shutil.copyfile(ROOT/name, output/name)
    python = output/'python'
    python.mkdir(exist_ok=True)
    for name in ('rosenau_hyman.py','porous_fin.py','dovetail_fin.py','solver_service.py',
                 'validated_combinations.json','validated_fin_combinations.json','validated_dovetail_combinations.json'):
        shutil.copyfile(ROOT/name, python/name)
    shutil.copyfile(ROOT/'output/pdf/Hermite_Porous_Fin_Comparison_Note.pdf', output/'comparison-note.pdf')
    print(f'Browser solver built: {output.resolve()}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'docs')
    build(parser.parse_args().output)

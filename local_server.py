"""Loopback-only Hermite PDE application. Run with the existing venv Python."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import math
from pathlib import Path
import secrets
import threading
import uuid

import numpy as np
from build_viewer import render_viewer, viewer_payload
from rosenau_hyman import compute, travelling_wave
import porous_fin

ROOT = Path(__file__).resolve().parent
DEFAULTS = dict(speed=1, shift=0, count=6, length=1, duration=1, reference=True)
COMBINATIONS = json.loads((ROOT / 'validated_combinations.json').read_text())
SUPPORTED = [row['parameters'] for row in COMBINATIONS if row['passed']]
FIN_COMBINATIONS = json.loads((ROOT / 'validated_fin_combinations.json').read_text())
FIN_SUPPORTED = [row['parameters'] for row in FIN_COMBINATIONS if row['passed']]
FIN_PRESETS = {row['preset']: {k:v for k,v in row['parameters'].items()
               if k in ('nc','nr','n1_squared','wet_term','power')} for row in FIN_COMBINATIONS if row['passed']}


def validate_request(data):
    if isinstance(data, dict) and data.get('problem') == 'porous_fin':
        params = {k:v for k,v in data.items() if k != 'problem'}
        porous_fin.validate_parameters(params)
        if {k:v for k,v in params.items() if k!='reference'} not in FIN_SUPPORTED:
            raise ValueError('Unsupported porous-fin combination. Choose a tested coefficient preset and resolution.')
        return dict(data)
    if isinstance(data, dict) and 'problem' in data:
        if data['problem'] != 'rosenau_hyman':
            raise ValueError('Unsupported problem family.')
        validate_request({k:v for k,v in data.items() if k!='problem'})
        return dict(data)
    if not isinstance(data, dict) or set(data) != set(DEFAULTS):
        raise ValueError('Provide exactly speed, shift, count, length, duration and reference.')
    if type(data['reference']) is not bool:
        raise ValueError('reference must be a boolean.')
    for key in set(DEFAULTS) - {'reference'}:
        if type(data[key]) not in (int, float) or not math.isfinite(data[key]):
            raise ValueError(f'{key} must be a finite number.')
    if type(data['count']) is not int:
        raise ValueError('count must be an integer.')
    if {k: v for k, v in data.items() if k != 'reference'} not in SUPPORTED:
        raise ValueError('Unsupported combination. Choose a tested combination; use N=6 for c=2, a=0.5, T=1.')
    return dict(data)


def solve_request(parameters):
    p = validate_request(parameters)
    if p.get('problem') == 'porous_fin':
        arrays, report = porous_fin.compute({k:v for k,v in p.items() if k!='problem'})
        report['parameters'] = p
        report['run_id'] = uuid.uuid4().hex
        return arrays, report
    initial, traces, reference = travelling_wave(p['speed'], p['shift'])
    arrays, report, _ = compute(p['count'], p['length'], p['duration'], initial, traces,
                               reference if p['reference'] else None)
    report['parameters'] = p
    report['run_id'] = uuid.uuid4().hex
    report['reference_compared_after_solve'] = p['reference']
    return arrays, report


class Application:
    def __init__(self, output):
        self.output = Path(output)
        self.lock = threading.Lock()
        self.worker = ThreadPoolExecutor(max_workers=1)
        self.token = secrets.token_urlsafe(32)
        self.state = dict(status='idle')
        self.latest = None

    def submit(self, parameters):
        p = validate_request(parameters)
        with self.lock:
            if self.state['status'] == 'running':
                raise RuntimeError('A solve is already running.')
            job = uuid.uuid4().hex
            self.state = dict(status='running', job_id=job, parameters=p)
            self.worker.submit(self._work, job, p)
        return job

    def _work(self, job, parameters):
        try:
            arrays, report = solve_request(parameters)
            if not report['validation_passed']:
                with self.lock:
                    self.state = dict(status='failed', job_id=job, parameters=parameters,
                                      error='Independent validation failed. Previous result retained.', report=report)
                return
            folder = self.output / report['run_id']
            folder.mkdir(parents=True, exist_ok=False)
            np.savez_compressed(folder / 'solution.npz', **arrays)
            (folder / 'results.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
            payload = viewer_payload(arrays, report)
            (folder / 'interactive.html').write_text(render_viewer(payload), encoding='utf-8')
            with self.lock:
                self.latest = payload
                self.state = dict(status='complete', job_id=job, parameters=parameters,
                                  run_id=report['run_id'], report=report)
        except Exception as exc:
            with self.lock:
                self.state = dict(status='failed', job_id=job, parameters=parameters, error=str(exc))


def make_server(port=8765, output=ROOT / 'local-runs'):
    app = Application(output)

    class Handler(BaseHTTPRequestHandler):
        def respond(self, status, data, content_type='application/json'):
            body = (json.dumps(data, allow_nan=False) if content_type == 'application/json' else data).encode('utf-8')
            self.send_response(status)
            self.send_header('Content-Type', content_type+'; charset=utf-8')
            self.send_header('Content-Length', str(len(body)))
            self.send_header('Cache-Control', 'no-store')
            self.send_header('X-Content-Type-Options', 'nosniff')
            self.end_headers()
            self.wfile.write(body)

        def local_host(self):
            return self.headers.get('Host') == f'127.0.0.1:{self.server.server_port}'

        def do_GET(self):
            if not self.local_host():
                return self.respond(403, {'error': 'Loopback host required.'})
            if self.path == '/':
                return self.respond(200, render_viewer(dict(live=True, token=app.token,
                    defaults=DEFAULTS, supported=SUPPORTED, fin_defaults=porous_fin.DEFAULTS,
                    fin_supported=FIN_SUPPORTED, fin_presets=FIN_PRESETS)), 'text/html')
            with app.lock:
                if self.path == '/api/status':
                    return self.respond(200, app.state)
                if self.path == '/api/result' and app.latest is not None:
                    return self.respond(200, app.latest)
            return self.respond(404, {'error': 'Not found.'})

        def do_POST(self):
            origin = f'http://127.0.0.1:{self.server.server_port}'
            if (not self.local_host() or self.headers.get('Origin') != origin
                    or self.headers.get('X-App-Token') != app.token):
                return self.respond(403, {'error': 'Local application request required.'})
            if self.path != '/api/solve':
                return self.respond(404, {'error': 'Not found.'})
            try:
                size = int(self.headers.get('Content-Length', '0'))
                if not 0 < size <= 4096 or self.headers.get('Content-Type') != 'application/json':
                    raise ValueError('Expected a small JSON request.')
                data = json.loads(self.rfile.read(size))
                job = app.submit(data)
                self.respond(202, {'job_id': job})
            except (ValueError, OverflowError) as exc:
                self.respond(400, {'error': str(exc)})
            except RuntimeError as exc:
                self.respond(409, {'error': str(exc)})

        def log_message(self, *args):
            pass

    server = ThreadingHTTPServer(('127.0.0.1', port), Handler)
    server.app = app
    return server


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--port', type=int, default=8765)
    args = parser.parse_args()
    server = make_server(args.port)
    print(f'Open http://127.0.0.1:{server.server_port} - press Ctrl+C to stop.', flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
        server.app.worker.shutdown(wait=True)

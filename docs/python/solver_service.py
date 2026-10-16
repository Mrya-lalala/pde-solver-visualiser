"""Shared numeric validation and real solver entry points for local and web use."""
import json
import math
from pathlib import Path
import uuid

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


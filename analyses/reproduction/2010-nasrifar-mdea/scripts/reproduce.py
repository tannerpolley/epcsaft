"""Declare neutral isobars; native Engine owns all thermodynamic calculations."""
import csv
from datetime import datetime, timezone
from decimal import Decimal
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time
import uuid

LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[2]
SOURCE = LEAF / 'inputs/source'
RESULTS = LEAF / 'results'
QA = ROOT / 'build/issue-734'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def serial(value):
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(v) for v in value]
    if type(value).__name__ in ('Stop', 'Category'):
        return value.name
    return serial({k: getattr(value, k) for k in dir(value)
                   if not k.startswith('_') and not callable(getattr(value, k))})


def dump(path, value):
    data = (json.dumps(serial(value), indent=2, allow_nan=False) + '\n').encode()
    path.write_bytes(gzip.compress(data, mtime=0) if path.suffix == '.gz' else data)


def declarations():
    from epcsaft import Mixture, Parameters

    return Mixture(Parameters.from_json(LEAF / 'inputs/parameters/parameters.json', components=['MDEA', 'water']))


def calculate(targets):
    import epcsaft
    from epcsaft import _core, trace
    from epcsaft.equilibrium import Free, MoleFractions, Phase, Pinned, Problem, solve_equilibrium

    assert Path(sys.prefix).resolve() == ROOT / '.venv', 'Use this worktree’s environment'
    started = time.perf_counter()
    parameter = LEAF / 'inputs/parameters/parameters.json'
    materialization = json.loads((parameter.parent / 'provenance.json').read_text())
    assert sha(parameter) == materialization['sha256'] == sha(ROOT / materialization['packet_path'] / materialization['packet_member'])
    model = declarations()
    identity = dict(run_id=str(uuid.uuid4()), started_utc=datetime.now(timezone.utc).isoformat(),
                    checkout_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                    checkout_status=subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True),
                    target_freeze_commit=subprocess.check_output(['git', 'log', '-1', '--format=%H', '--', str(SOURCE / 'figure-targets.csv')], cwd=ROOT, text=True).strip(),
                    engine_provenance=serial(epcsaft.provenance()), extension_path=_core.__file__,
                    extension_sha256=sha(Path(_core.__file__)), interpreter=sys.executable,
                    producer_sha256=sha(Path(__file__)), materialization=materialization,
                    inputs_sha256={str(p.relative_to(LEAF)): sha(p) for p in sorted(SOURCE.iterdir()) if p.is_file()})
    run = dict(execution_identity=identity, points=[], wall_seconds=None)
    traces = []
    for row in targets:
        if row['panel'] == '15':
            continue
        for offset in [-1, 0, 1]:
            x = float(row['x']) + offset * float(row['x_half_width'])
            branch = 'bubble' if row['series'] == 'measured' else row['series']
            # Each endpoint is an independent native declaration, not a Python root search.
            tolerances = [1e-10, 1e-11] if offset == 0 and row['row_id'].endswith('/1') and row['series'] != 'measured' else [1e-10]
            for tolerance in tolerances:
                point = dict(row_id=row['row_id'], offset=offset, feed_mdea=x,
                             requested_tolerance=tolerance, run_id=identity['run_id'],
                             pressure_Pa=float(row['pressure_Pa']), temperature_K=None)
                phases = [Phase('L', kind='liquid', amount=Pinned(0.0) if branch == 'dew' else Free()),
                          Phase('V', kind='vapor', amount=Pinned(0.0) if branch == 'bubble' else Free())]
                problem = Problem(T=Free(float(row['y'])), P=point['pressure_Pa'],
                                  feed=MoleFractions({'MDEA': x, 'water': 1-x}),
                                  phases=phases, row_tolerance=tolerance)
                clock = time.perf_counter()
                with trace(diagnostics=True) as capture:
                    try:
                        result = solve_equilibrium(model, problem)
                        point.update({k: serial(getattr(result, k)) for k in
                                      ['success', 'message', 'requested_tolerance_met', 'mole_fractions',
                                       'molar_densities', 'packing_fractions', 'amounts', 'residuals',
                                       'classified_residuals', 'convergence', 'solver_options']})
                        point['temperature_K'] = result.temperature
                    except Exception as error:
                        point.update(failure_type=type(error).__name__, failure=str(error))
                point['wall_seconds'] = time.perf_counter() - clock
                point['trace_index'] = len(traces)
                traces.append(serial(capture.trace))
                spans = traces[-1]['spans']
                point['inner_residual_maxima'] = {name: max((s['convergence']['residual'] for s in spans
                    if s['name'] == name and s['convergence'] is not None), default=None)
                    for name in ['association_newton', 'solve_density']}
                run['points'].append(point)
                print(row['row_id'], offset, tolerance, point['temperature_K'], point.get('message', point.get('failure')), flush=True)
    run['wall_seconds'] = time.perf_counter() - started
    QA.mkdir(parents=True, exist_ok=True)
    (QA / 'native-traces.json.gz').write_bytes(gzip.compress(json.dumps(traces, allow_nan=False).encode(), mtime=0))
    dump(RESULTS / 'native-results.json.gz', run)
    return run


def compare(targets, run):
    records = []
    for target in targets:
        rows = [p for p in run['points'] if p['row_id'] == target['row_id'] and p['requested_tolerance'] == 1e-10]
        values = [p['temperature_K'] for p in rows if p.get('success') and p.get('requested_tolerance_met') and p['temperature_K'] is not None]
        row = dict(row_id=target['row_id'], panel=target['panel'], series=target['series'], evidence=target['evidence'],
                   source_x=target['x'], source_y=target['y'], x_half_width=target['x_half_width'],
                   y_half_width=target['y_half_width'], y_scale=target['y_scale'], pressure_Pa=target['pressure_Pa'],
                   prediction='', delta='', prediction_low='', prediction_high='', outcome='', cause='', run_id='')
        row['extrapolation'] = ''
        if len(values) == 3:
            prediction = rows[1]['temperature_K']
            if prediction > 459 or prediction < 350:
                row['extrapolation'] = 'outside Table 3 MDEA-water k_ij fit interval 350-459 K; Table 2 pure fit intervals respected'
            low, high = min(values), max(values)
            reference, half = float(target['y']), float(target['y_half_width'])
            row.update(prediction=prediction, delta=prediction-reference, prediction_low=low, prediction_high=high,
                       outcome='reproduced' if high >= reference-half and low <= reference+half else 'not reproduced', run_id=rows[1]['run_id'])
            if row['outcome'] == 'not reproduced':
                row['cause'] = 'converged disagreement; cause undetermined'
            if target['condition_note']:
                row.update(outcome='condition unresolved', cause=target['condition_note'])
        elif target['panel'] == '15':
            row.update(outcome='unavailable', cause='Engine #776: ideal mole-fraction chemistry and native effective-water conversion')
        else:
            row.update(outcome='undetermined', cause='native result absent or requested residual unmet; see native record')
        records.append(row)
    for filename, field, scale in [('kamps-2001-tables-2-3.csv', 'ten_times_total_pressure_MPa', Decimal('100000')),
                                   ('rho-1997-tables-1-3-selected.csv', 'CO2_partial_pressure_kPa', Decimal('1000'))]:
        for i, source in enumerate(read(SOURCE / filename), 1):
            printed = Decimal(source[field])
            row = dict.fromkeys(records[0], '')
            row.update(row_id=f'{filename}/{i}', panel=('16a' if source['gas'] == 'CO2' else '16b') if filename.startswith('kamps') else ('17a' if source['temperature_K'] == '323.15' else '17b'),
                       series='measured', evidence='physical comparison', source_y=str(printed*scale),
                       y_half_width=str(Decimal(10)**printed.as_tuple().exponent*scale/2), y_scale='linear Pa',
                       outcome='unavailable', cause='Engine #776: ideal mole-fraction chemistry and native effective-water conversion')
            records.append(row)
    with (RESULTS / 'comparisons.csv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=records[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    return records


if __name__ == '__main__':
    targets = read(SOURCE / 'figure-targets.csv')
    run = json.loads(gzip.decompress((RESULTS / 'native-results.json.gz').read_bytes())) if '--compare-only' in sys.argv else calculate(targets)
    compare(targets, run)

"""Declare independent source states; the native Engine owns every thermodynamic solve."""
import csv
from datetime import datetime, timezone
from decimal import Decimal
from enum import Enum
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
SOURCE, RESULTS = LEAF / 'inputs/source', LEAF / 'results'
QA = ROOT / 'build/issue-736'
CRITERIA = json.loads((SOURCE / 'numerical-criteria.json').read_text())
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
    if isinstance(value, (list, tuple)):
        return [serial(v) for v in value]
    if isinstance(value, Enum):
        return value.name
    return {k: serial(getattr(value, k)) for k in dir(value)
            if not k.startswith('_') and not callable(getattr(value, k))}

def dump(path, value):
    data = (json.dumps(serial(value), indent=2, allow_nan=False) + '\n').encode()
    path.write_bytes(gzip.compress(data, mtime=0) if path.suffix == '.gz' else data)

def targets():
    rows = []
    for figure in [1, 4]:
        for source in read(SOURCE / f'figure-{figure}-model-points.csv'):
            pure = figure == 1
            rows.append(dict(row_id=source['row_id'], figure=str(figure), branch='bubble' if pure else source['series'],
                mdea='MDEA_' + ('2B' if source['mdea_scheme'].startswith('2') else source['mdea_scheme'][0]+'C'),
                water=None if pure else 'water_' + ('2B' if source['water_scheme'].startswith('2') else '4C'),
                pressure_Pa=float(source['pressure_Pa']), source_x=float(source['pressure_Pa'] if pure else source['x']),
                x_half_width=float(source['pressure_half_width_Pa'] if pure else source['x_half_width']),
                source_y=float(source['density_kg_per_m3'] if pure else source['temperature_K']),
                y_half_width=float(source['density_half_width_kg_per_m3'] if pure else source['temperature_half_width_K']),
                unit='kg/m3' if pure else 'K'))
    return rows

def declarations():
    from epcsaft import Mixture, Pure, Parameters
    keys = {(r['mdea'], r['water']) for r in targets()}
    return {key:(Pure if key[1] is None else Mixture)(Parameters.from_json(
        LEAF/'inputs/parameters/parameters.json', components=[c for c in key if c is not None])) for key in keys}

def calculate(rows):
    import epcsaft
    from epcsaft import equilibrium as eq
    from epcsaft import _core, trace
    from epcsaft.equilibrium import Free, MoleFractions, Phase, Pinned, Problem, solve_equilibrium
    from epcsaft.regression import observation, predict

    assert Path(sys.prefix).resolve() == ROOT / '.venv', 'Use this worktree environment'
    started = time.perf_counter()
    parameter = LEAF / 'inputs/parameters/parameters.json'
    materialization = json.loads((parameter.parent / 'provenance.json').read_text())
    assert sha(parameter) == materialization['sha256'] == sha(ROOT / materialization['packet_path'] / materialization['packet_member'])
    identity = dict(run_id=str(uuid.uuid4()), started_utc=datetime.now(timezone.utc).isoformat(),
        checkout_commit=subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip(),
        checkout_status=subprocess.check_output(['git','status','--porcelain'], cwd=ROOT, text=True),
        target_freeze_commit=subprocess.check_output(['git','log','-1','--format=%H','--',str(SOURCE)], cwd=ROOT, text=True).strip(),
        engine_provenance=serial(epcsaft.provenance()), extension_path=_core.__file__, extension_sha256=sha(Path(_core.__file__)),
        interpreter=sys.executable, producer_sha256=sha(Path(__file__)), materialization=materialization,
        input_sha256={str(p.relative_to(LEAF)):sha(p) for p in sorted(SOURCE.iterdir()) if p.is_file()})
    run = dict(execution_identity=identity, equation='Engine sigma-cubed association; raw printed kappa; Figure 4 k_ij=0', points=[])
    traces = []
    models = declarations()
    for row in rows:
        pure = row['figure'] == '1'
        names = [row['mdea']] + ([] if pure else [row['water']])
        model = models[(row['mdea'],row['water'])]
        parameters = model.parameters
        for offset in [-1, 0, 1]:
            coordinate = row['source_x'] + offset*row['x_half_width']
            tolerances = [CRITERIA['row_tolerance']]
            if offset == 0 and (pure or row['pressure_Pa'] == 40000):
                tolerances.append(CRITERIA['refinement_tolerance'])
            for tolerance in tolerances:
                phases = [Phase('L', kind='liquid', amount=Pinned(0) if row['branch']=='dew' else Free()),
                          Phase('V', kind='vapor', amount=Pinned(0) if row['branch']=='bubble' else Free())]
                problem = Problem(T=Free(420 if pure else row['source_y']), P=coordinate if pure else row['pressure_Pa'],
                    feed=MoleFractions({names[0]:1} if pure else {names[0]:coordinate,names[1]:1-coordinate}),
                    phases=phases, row_tolerance=tolerance)
                point = dict(row_id=row['row_id'], offset=offset, coordinate=coordinate, requested_tolerance=tolerance,
                             run_id=identity['run_id'], prediction=None)
                clock = time.perf_counter()
                with trace(diagnostics=True) as capture:
                    try:
                        result = solve_equilibrium(model, problem)
                        point.update({k:serial(getattr(result,k)) for k in ['success','message','solver_status','requested_tolerance_met',
                            'temperature','pressure','mole_fractions','molar_densities','packing_fractions','amounts',
                            'residuals','classified_residuals','convergence','solver_options']})
                        point['prediction'] = result.temperature
                        if result.success:
                            observable = eq.SolvedStateObservable(eq.SolvedStateObservableKind.PhaseMolarDensity if pure else eq.SolvedStateObservableKind.Temperature)
                            direction = eq.SolvedStateActionDirection(0, 1 if pure else 0, [0] if pure else [1,-1], [])
                            action = eq.solved_state_actions(eq.compile_problem(model,problem),result,eq.SolvedStateActionRequest([observable],[direction]))
                            point['branch_action'] = serial(action)
                        if pure:
                            obs = observation(parameters,'phase_density',problem,phase='L',observed=row['source_y'],
                                form='difference',mass_basis=True,scale=row['y_half_width'])
                            density = predict(model, obs)
                            point.update(prediction=density.value, density_observation=serial(density))
                    except Exception as error:
                        point.update(failure_type=type(error).__name__, failure=str(error))
                point['wall_seconds'] = time.perf_counter()-clock
                point['trace_index'] = len(traces)
                traces.append(serial(capture.trace))
                point['inner_residual_maxima'] = {name:max((s['convergence']['residual'] for s in traces[-1]['spans']
                    if s['name']==name and s['convergence'] is not None),default=None) for name in ['association_newton','solve_density']}
                run['points'].append(point)
                print(row['row_id'],offset,tolerance,point['prediction'],point.get('message',point.get('failure')),flush=True)
    run['wall_seconds'] = time.perf_counter()-started
    QA.mkdir(parents=True,exist_ok=True)
    dump(QA/'native-traces.json.gz',traces)
    dump(RESULTS/'native-results.json.gz',run)
    return run

def compare(rows, run):
    records = []
    for row in rows:
        points = [p for p in run['points'] if p['row_id']==row['row_id'] and p['requested_tolerance']==CRITERIA['row_tolerance']]
        valid = [p for p in points if p.get('success') and p.get('requested_tolerance_met') and p['prediction'] is not None
                 and (row['figure']!='1' or p['density_observation']['requested_tolerance_met'])]
        values = [p['prediction'] for p in valid]
        monotone = len(values)==3 and (values[0]<=values[1]<=values[2] or values[0]>=values[1]>=values[2])
        slopes = [p.get('branch_action',{}).get('results',[{}])[0].get('action') for p in valid]
        monotone = monotone and len(slopes)==3 and all(v is not None for v in slopes) and (all(v>0 for v in slopes) or all(v<0 for v in slopes))
        branches = len(valid)==3 and all(p['molar_densities'][0]>p['molar_densities'][1] for p in valid)
        refined = [p for p in run['points'] if p['row_id']==row['row_id'] and p['requested_tolerance']==CRITERIA['refinement_tolerance']]
        refinement = abs(refined[0]['prediction']-values[1]) if refined and refined[0]['prediction'] is not None and len(valid)==3 else None
        okay = monotone and branches and (not refined or (refined[0].get('requested_tolerance_met') and refinement is not None and refinement<row['y_half_width']*CRITERIA['refinement_fraction_of_source_half_width']))
        interval = row['figure']!='1' or all(CRITERIA['pure_temperature_interval_K'][0]<=p['temperature']<=CRITERIA['pure_temperature_interval_K'][1] for p in valid)
        overlap = bool(values) and max(values)>=row['source_y']-row['y_half_width'] and min(values)<=row['source_y']+row['y_half_width']
        outcome = ('reproduced' if overlap and interval else 'not reproduced') if okay else 'undetermined'
        records.append(row | dict(evidence='numerical verification',equation=run['equation'],run_id=run['execution_identity']['run_id'],
            prediction=values[1] if len(valid)==3 else '',prediction_low=min(values) if values else '',prediction_high=max(values) if values else '',
            delta=values[1]-row['source_y'] if len(valid)==3 else '',refinement_change=refinement,valid_solves=len(valid),
            monotone_samples=monotone,separated_phases=branches,pure_temperature_interval_met=interval,outcome=outcome,
            solve_outcome='admitted' if len(valid)==3 else '; '.join(sorted({p.get('message',p.get('failure','')) for p in points})),
            cause='' if outcome=='reproduced' else ('raw σ³ disagreement; see retained paired-convention diagnosis' if okay else 'failed solve, branch change, unresolved extremum or refinement; see native evidence')))
    for filename in ['kuranov-1996-table-2.csv','kamps-2001-table-2-co2.csv']:
        for i,source in enumerate(read(SOURCE/filename),1):
            printed = Decimal(source['ten_times_total_pressure_MPa'])
            record = dict.fromkeys(records[0],'')
            record.update(row_id=f'{filename}/{i}',figure='8' if filename.startswith('kuranov') else '9',branch='measured total pressure',
                source_y=str(printed*100000), y_half_width=str(Decimal(10)**printed.as_tuple().exponent*50000),unit='Pa',
                evidence='physical comparison',outcome='unavailable',cause='Engine #776 ideal chemistry; printed Eq. (14)/Table 5 source conflicts and unspecified molecular-feed mapping')
            records.append(record)
    with (RESULTS/'comparisons.csv').open('w') as handle:
        writer = csv.DictWriter(handle,fieldnames=records[0],lineterminator='\n');writer.writeheader();writer.writerows(records)

if __name__ == '__main__':
    rows = targets()
    run = json.loads(gzip.decompress((RESULTS/'native-results.json.gz').read_bytes())) if '--compare-only' in sys.argv else calculate(rows)
    compare(rows,run)

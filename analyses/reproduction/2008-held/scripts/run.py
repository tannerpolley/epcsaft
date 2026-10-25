"""Held 2008 declarations, native results and fixed source-interval comparisons."""
import csv
from datetime import datetime, timezone
from functools import lru_cache
import gzip
import hashlib
import json
import math
from multiprocessing import Pool
import os
from pathlib import Path
import subprocess
import sys
import time
import uuid

import epcsaft
from epcsaft import Mixture, Parameters, _core, regression as reg
from epcsaft import equilibrium as q

LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[2]
SOURCE, OUT = LEAF/'inputs/source', LEAF/'results'
PACKET = ROOT/'data/packets/held-cameretti-sadowski-2008/4/parameters/parameters.json'
SUPPLEMENT = json.loads((SOURCE/'supplements.json').read_text())
SALTS = {r['salt']: r for r in csv.DictReader((SOURCE/'salt-components.csv').open())}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serial(v):
    if v is None or isinstance(v, (str, bool, int)):
        return v
    if isinstance(v, float):
        return v if math.isfinite(v) else str(v)
    if isinstance(v, dict):
        return {str(k): serial(x) for k, x in v.items()}
    if isinstance(v, (tuple, list)):
        return [serial(x) for x in v]
    if hasattr(v, 'name') and not hasattr(v, 'category'):
        return v.name
    return {k: serial(getattr(v, k)) for k in dir(v) if not k.startswith('_') and not callable(getattr(v, k))}


def fields(v, names):
    return serial({k: getattr(v, k) for k in names.split()})


def write(name, rows):
    with (OUT/name).open('w') as f:
        w = csv.DictWriter(f, fieldnames=rows[0], lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


@lru_cache(None)
def declaration(salt, T, quantity, mode='primary', prepare_only=False):
    row = SALTS[salt]
    ids = ['water', row['cation'], row['anion']]
    nu = (int(row['nu_cation']), int(row['nu_anion']))
    doc = json.loads(PACKET.read_text())
    if T != 298.15 or quantity == 'density_difference' or mode != 'primary':
        doc['document_id'] = f'held-2008-local-{T}-{mode}'
        doc['domains'].append(dict(domain_id='figure-condition', kind='reported-conditions',
                                   temperature_min=dict(magnitude=T, unit='kelvin'), temperature_max=dict(magnitude=T, unit='kelvin')))
        if T != 298.15 and mode != 'fixed-dielectric':
            doc['correlations'][0]['constant']['value']['magnitude'] = SUPPLEMENT['dielectric_constants'][str(T)]
        if quantity == 'density_difference':
            for c in doc['components']:
                mass = SUPPLEMENT['masses'][c['component_id']]
                c['formula'] = mass['formula']
                c['fixed']['molar_mass'] = {'value': dict(magnitude=mass['molar_mass_kg_mol'], unit='kilogram / mole')}
        if mode == 'short-screening':
            for c in doc['components']:
                for coefficient in c['coefficients']:
                    if coefficient['family'] == 'debye_huckel_diameter':
                        coefficient['value']['magnitude'] *= 0.88
        params = Parameters.from_mapping(doc, components=ids)
    else:
        params = Parameters.from_json(PACKET, components=ids)
    model = Mixture(params)
    solvent = [float(cid == 'water') for cid in params.component_ids]
    ref = None if prepare_only or quantity in ('vapor_pressure', 'density_difference') else model.reference(
        T, 100000, solvent_fractions=solvent,
        **(dict(ions=ids[1:], stoichiometry=nu) if quantity == 'mean_ionic_activity_coefficient_molal' else {}))
    return params, model, ref, solvent, nu, hashlib.sha256(json.dumps(doc, sort_keys=True).encode()).hexdigest()


def declarations():
    for panel in json.loads((SOURCE/'figure-calibration.json').read_text())['panels']:
        declaration(panel['salt'], panel['temperature_K'], panel['quantity'], prepare_only=True)


@lru_cache(None)
def calculate(salt, T, quantity, m, mode='primary'):
    record = dict(salt=salt, temperature_K=T, quantity=quantity, molality=m, mode=mode, value=None, cause='')
    try:
        params, model, ref, solvent, nu, doc_hash = declaration(salt, T, quantity, mode)
        ions = tuple(params.component_ids.index(SALTS[salt][k]) for k in ('cation', 'anion'))
        x = _core.formula_composition(solvent, params._native, ions, nu, m) if m > 0 else solvent
        record.update(component_ids=params.component_ids, composition=x, declaration_sha256=doc_hash)
        if quantity == 'vapor_pressure':
            problem = q.Problem(T=T, P=q.Free(), feed=q.Molalities(
                {params.component_ids[i]: n*m for i, n in zip(ions, nu)}, solvent=['water'], solvent_mass_kg=1),
                phases=[q.Phase('L', kind='liquid'), q.Phase('V', kind='vapor', support=('water',), amount=q.Pinned(0.0))])
            result = q.solve_equilibrium(model, problem)
            record.update(fields(result, 'success message pressure convergence residuals mole_fractions amounts molar_densities iterations solver_options'))
            record['value'] = result.pressure/1e5 if result.pressure is not None else None
        else:
            state = model.state(T, P=100000, x=x, phase='liquid')
            record.update(fields(state, 'pressure pressure_density_derivative density_convergence density_branch'))
            if quantity == 'density_difference':
                obs = reg.density(params, T, 100000, dict(zip(params.component_ids, x)), phase='liquid',
                                  observed=0, form='difference', scale=1, mass_basis=True)
                pred = reg.predict(model, obs)
                record.update(fields(pred, 'value status message convergence requested_tolerance_met'))
                record['mass_density_kg_m3'] = pred.value
                record['pure_water_density_kg_m3'] = calculate(salt, T, quantity, 0, mode)['mass_density_kg_m3'] if m else pred.value
                record['value'] = pred.value-record['pure_water_density_kg_m3'] if pred.value is not None else None
            else:
                record['reference'] = fields(ref._native, 'molality convergence_error solvent_cause solvent_convergence solver_options')
                if quantity == 'water_activity_coefficient':
                    activity = ref.solvent_activity(x=x)
                    aw, xw = activity.solvent_activity, x[params.component_ids.index('water')]
                    record.update(fields(activity, 'status cause'))
                    record.update(solvent_activity=aw, water_mole_fraction=xw)
                    record['value'] = aw/xw if aw is not None and math.isfinite(aw) and xw > 0 else None
                    record['normalization_error'] = abs(xw*record['value']-aw) if record['value'] is not None else None
                elif m:
                    activity = ref.activity(m)
                    record.update(value=activity.mean_ionic_molality_scale, osmotic_coefficient=serial(activity.osmotic_coefficient))
                else:
                    record.update(value=1.0, cause='ionic infinite-dilution limit; no activity(0) call')
    except Exception as error:
        record.update(cause=f'{type(error).__name__}: {error}')
    convergence = record.get('density_convergence', record.get('convergence'))
    record['numerically_qualified'] = bool(convergence and convergence['residual'] <= convergence['criterion']['value']
        and record.get('pressure_density_derivative', 1) > 0 and isinstance(record['value'], (float, int))
        and math.isfinite(record['value']) and not record.get('reference', {}).get('solvent_cause'))
    return serial(record)


def series(rows, mode='primary'):
    comparisons, records = [], {}
    for index, row in rows:
        m, dx = float(row['molality']), float(row['molality_half_width'])
        lo, hi = max(0, m-dx), max(0, m+dx)
        samples = [lo+(hi-lo)*f for f in SUPPLEMENT['native_interval_samples']]
        samples[2] = max(0, m)
        values = []
        for abscissa in samples:
            key = (row['salt'], float(row['temperature_K']), row['quantity'], abscissa, mode)
            result = records.setdefault(key, calculate(*key))
            values.append(result['value'])
        finite = all(isinstance(v, (float, int)) and math.isfinite(v) for v in values)
        qualified = all(records[(row['salt'], float(row['temperature_K']), row['quantity'], a, mode)]['numerically_qualified'] for a in samples)
        y, dy = float(row['value']), float(row['value_half_width'])
        overlap = finite and max(values) >= y-dy and min(values) <= y+dy
        monotone = finite and (all(b >= a for a, b in zip(values, values[1:])) or all(b <= a for a, b in zip(values, values[1:])))
        outcome = 'reproduced' if overlap else 'not reproduced' if monotone else 'undetermined'
        cause = '' if overlap else 'directional native range disjoint from source ordinate interval' if monotone else 'unresolved interior extremum or unavailable native value'
        if row['qa'] != 'pass':
            outcome, cause = 'undetermined', row['cause']
        if not qualified:
            outcome, cause = 'undetermined', 'native residual, mechanical slope or reference not admitted'
        row.update(target_id=index, mode=mode, native_value=values[2], signed_error=values[2]-y if finite else None,
                   abscissae=json.dumps(samples), native_ordinates=json.dumps(values), outcome=outcome, comparison_cause=cause,
                   condition='unchanged packet 4' if float(row['temperature_K']) == 298.15 else 'conditional: Ascani-Held 2021 dielectric; fixed P assumed 1 bar where applicable')
        comparisons.append(row)
    print(rows[0][1]['figure'], rows[0][1]['panel'], rows[0][1]['salt'], len(comparisons), flush=True)
    return comparisons, list(records.values())


def main():
    started = time.perf_counter()
    OUT.mkdir(exist_ok=True)
    for line in (SOURCE/'frozen-targets.sha256').read_text().splitlines():
        expected, path = line.split('  ', 1)
        assert digest(LEAF/path) == expected, path
    if '--compare-only' in sys.argv:
        # Reclassify source ownership using only the original native evaluations.
        global calculate
        original = list(csv.DictReader(gzip.decompress((OUT/'figure-comparisons.csv.gz').read_bytes()).decode().splitlines()))
        targets = list(csv.DictReader((SOURCE/'figure-targets.csv').open()))
        native = [json.loads(line) for line in gzip.decompress((OUT/'native-results.jsonl.gz').read_bytes()).splitlines()]
        retained = {(r['salt'], r['temperature_K'], r['quantity'], r['molality'], r['mode']): r for r in native}
        calculate = lambda *key: retained[key]  # A missing retained evaluation stops; no Engine fallback.
        comparisons, probes = [], []
        for index, (old, target) in enumerate(zip(original, targets, strict=True)):
            assert all(old[k] == v for k, v in target.items() if k not in ('qa', 'cause')), index
            if old['cause'].startswith('same source ink'):
                target.update(qa='review', cause=old['cause'])
            updated = series([(index, target.copy())])[0][0]
            updated['execution_run_id'] = old['execution_run_id']
            assert all(str(updated[k]) == v for k, v in old.items() if k not in ('qa', 'cause', 'outcome', 'comparison_cause')), index
            comparisons.append(updated)
            if target['source_kind'] == 'line' and 'activity_coefficient' in target['quantity']:
                probe = series([(index, target.copy())], 'short-screening')[0][0]
                probe['execution_run_id'] = old['execution_run_id']
                probes.extend((updated, probe))
        write('figure-comparisons.csv', comparisons)
        path = OUT/'figure-comparisons.csv'
        path.with_suffix('.csv.gz').write_bytes(gzip.compress(path.read_bytes(), mtime=0))
        path.unlink()
        columns = next(csv.DictReader((OUT/'convention-probes.csv').open())).keys()
        probes.sort(key=lambda r: (r['mode'] != 'primary', r['target_id']))
        write('convention-probes.csv', [{k:r[k] for k in columns} for r in probes])
        summary = json.loads((OUT/'run-summary.json').read_text())
        summary['outcomes'] = {o: sum(r['outcome'] == o for r in comparisons) for o in ('reproduced', 'not reproduced', 'undetermined')}
        (OUT/'run-summary.json').write_text(json.dumps(summary, indent=2)+'\n')
        return
    identity = dict(run_id=str(uuid.uuid4()), started_utc=datetime.now(timezone.utc).isoformat(),
                    source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
                    dirty_paths=subprocess.check_output(['git', 'status', '--porcelain'], text=True).splitlines(),
                    engine_provenance=serial(epcsaft.provenance()), core_path=_core.__file__, core_sha256=digest(Path(_core.__file__)),
                    runner_sha256=digest(Path(__file__)), packet_sha256=digest(PACKET), supplements_sha256=digest(SOURCE/'supplements.json'),
                    frozen_targets_sha256=digest(SOURCE/'frozen-targets.sha256'), python=sys.executable)
    (OUT/'run-identity.json').write_text(json.dumps(identity, indent=2)+'\n')
    grouped = {}
    targets = list(csv.DictReader((SOURCE/'figure-targets.csv').open()))
    ink = {}
    for r in targets:
        if r['source_kind'] == 'line':
            ink.setdefault((r['figure'], r['pixel_x'], r['pixel_y']), set()).add((r['panel'], r['salt']))
    for i, r in enumerate(targets):
        if r['source_kind'] == 'line' and len(ink[(r['figure'], r['pixel_x'], r['pixel_y'])]) > 1:
            r.update(qa='review', cause='same source ink assigned to overlapping model series; separate line identity unresolved')
        grouped.setdefault((r['figure'], r['panel'], r['salt']), []).append((i, r))
    with Pool(os.cpu_count()) as pool:
        batches = pool.map(series, grouped.values())
    comparisons = [r for batch, _ in batches for r in batch]
    records = [r for _, batch in batches for r in batch]
    activity = [r.copy() for r in comparisons if r['source_kind']=='line' and 'activity_coefficient' in r['quantity']]
    with Pool(os.cpu_count()) as pool:
        probe_batches = pool.starmap(series, [([(r['target_id'], r.copy()) for r in activity if (r['figure'], r['panel'], r['salt'])==key], 'short-screening')
                                               for key in dict.fromkeys((r['figure'], r['panel'], r['salt']) for r in activity)])
    probes = activity+[r for batch, _ in probe_batches for r in batch]
    records += [r for _, batch in probe_batches for r in batch]
    # Fixed convention probes, declared before examining results; never replace the primary curve.
    fixed_probes = [('CsBr', 298.15, 'mean_ionic_activity_coefficient_molal', 1),
              ('CsBr', 293.15, 'density_difference', 4), ('CsCl', 333.15, 'vapor_pressure', 6),
              ('NaI', 303.15, 'water_activity_coefficient', 3), ('Cs2SO4', 303.15, 'water_activity_coefficient', 1)]
    diagnostics = [calculate(*p, mode) for p in fixed_probes for mode in ('primary', 'short-screening', 'fixed-dielectric')]
    shared = calculate('LiCl', 298.15, 'mean_ionic_activity_coefficient_molal', .001)
    shared.update(gamma_error=abs(shared['value']-0.9645487485200079), osmotic_error=abs(shared['osmotic_coefficient']-0.988174383011451))
    assert shared['gamma_error'] <= 1e-12 and shared['osmotic_error'] <= 1e-12, shared
    for r in comparisons+probes:
        r['execution_run_id'] = identity['run_id']
    write('convention-probes.csv', [{k:r[k] for k in ('target_id','figure','panel','salt','temperature_K','quantity','molality','molality_half_width','value','value_half_width','qa','cause','mode','native_value','signed_error','abscissae','native_ordinates','outcome','comparison_cause','execution_run_id')} for r in probes])
    write('figure-comparisons.csv', comparisons)
    csv_path = OUT/'figure-comparisons.csv'
    csv_path.with_suffix('.csv.gz').write_bytes(gzip.compress(csv_path.read_bytes(), mtime=0))
    csv_path.unlink()
    for r in records+diagnostics+[shared]:
        r['execution_run_id'] = identity['run_id']
    raw = ''.join(json.dumps(r, allow_nan=False)+'\n' for r in records).encode()
    (OUT/'native-results.jsonl.gz').write_bytes(gzip.compress(raw, mtime=0))
    (OUT/'diagnostics.json').write_text(json.dumps(dict(shared_point=shared, fixed_probes=diagnostics), indent=2)+'\n')
    summary = dict(wall_seconds=time.perf_counter()-started, source_rows=len(comparisons), native_rows=len(records),
                   execution_run_id=identity['run_id'], outcomes={o: sum(r['outcome'] == o for r in comparisons) for o in ('reproduced', 'not reproduced', 'undetermined')})
    (OUT/'run-summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()

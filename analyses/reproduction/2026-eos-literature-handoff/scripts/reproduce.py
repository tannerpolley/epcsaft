"""Declare the selected source comparisons; native Engine owns all thermodynamics."""
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys
import time
import uuid

import epcsaft
from epcsaft import Mixture, Parameters, _core, regression as reg

LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[2]
PACKET = ROOT / 'data/packets/figiel-2025-reference-electrolytes'
TARGETS = LEAF / 'inputs/source/literature-targets.csv'
SOLVENT = LEAF / 'inputs/source/solvent-compositions.csv'
MASSES = LEAF / 'inputs/source/density-molar-masses.csv'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fields(value, names):
    return {name: getattr(value, name) for name in names.split()}


def serial(value):
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else str(value)
    if isinstance(value, dict):
        return {str(k): serial(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [serial(v) for v in value]
    if hasattr(value, 'name') and not hasattr(value, 'category'):
        return value.name
    return serial({k: getattr(value, k) for k in dir(value)
                   if not k.startswith('_') and not callable(getattr(value, k))})


def parameters(version, ids, density=False):
    document = json.loads((PACKET / version / 'parameters/parameters.json').read_text())
    if density:
        masses = {r['component_id']: r for r in csv.DictReader(MASSES.open())}
        for component in document['components']:
            if component['component_id'] in masses:
                row = masses[component['component_id']]
                component['fixed']['molar_mass'] = {'value': {'magnitude': float(row['molar_mass_kg_mol']), 'unit': 'kilogram / mole'}}
    return Parameters.from_mapping(document, components=ids)


def declarations():
    solvent = next(csv.DictReader(SOLVENT.open()))
    cases = [('hamer-nacl-1m', 'activity', ('water', 'sodium-cation', 'chloride-anion'), (1, 0, 0), 1.0),
             ('esteso-nacl-60wt-0.5m', 'activity', ('water', 'ethanol', 'sodium-cation', 'chloride-anion'),
              (float(solvent['water_mole_fraction']), float(solvent['ethanol_mole_fraction']), 0, 0), 0.5),
             ('held-nabr-ethanol-0.1422m', 'density', ('ethanol', 'sodium-cation', 'bromide-anion'), (1, 0, 0), 0.1422)]
    for ion in ('lithium-cation', 'sodium-cation', 'potassium-cation', 'chloride-anion', 'bromide-anion'):
        pair = (ion, 'chloride-anion') if ion.endswith('cation') else ('sodium-cation', ion)
        cases.append(('figiel-' + ion, ion, ('water', *pair), (1, 0, 0), 0))
    return [(case, kind, parameters('4', ids, kind == 'density'), x, m)
            for case, kind, ids, x, m in cases]


def calculate(case, kind, params, solvent, molality, identity, version='4'):
    started = time.perf_counter()
    record = dict(case_id=case, execution_identity=identity, packet_version=version,
                  packet_sha256=digest(PACKET / version / 'parameters/parameters.json'),
                  component_ids=params.component_ids, solvent_fractions=solvent, molality=molality,
                  temperature_K=298.15, pressure_Pa=100000.0, ions=params.component_ids[-2:], stoichiometry=(1, 1))
    if kind == 'density':
        record['molar_mass_declaration_sha256'] = digest(MASSES)
    try:
        model = Mixture(params)
        if kind == 'density':
            x = _core.formula_composition(solvent, params._native, (1, 2), (1, 1), molality)
            obs = reg.density(params, 298.15, 100000.0, dict(zip(params.component_ids, x)),
                              phase='liquid', observed=795.96, form='difference', mass_basis=True, scale=0.005)
            state = model.state(298.15, P=100000.0, x=x, phase='liquid')
            record.update(composition=x, molar_masses_kg_mol=params.molar_masses,
                          density_diagnostic_recomputation=fields(state, 'molar_density pressure pressure_density_derivative density_convergence density_branch stable_root_count'))
        else:
            reference = model.reference(298.15, 100000.0, solvent_fractions=solvent,
                                        ions=params.component_ids[-2:], stoichiometry=(1, 1))
            record['reference'] = fields(reference, 'molality convergence_error log_fugacity')
            record['reference_solver_options'] = reference._native.solver_options
            record['reference_role'] = 'consumed activity reference' if kind == 'activity' else 'separate diagnostic recomputation; not predictor consumed record'
            if kind == 'activity':
                activity = reference.activity(molality)
                record.update(value=activity.mean_ionic_molality_scale, unit='1', status='Available',
                              solver_options=activity.solver_options)
            else:
                obs = reg.electrolyte(params, 'solvation_gibbs_energy', 298.15, 100000.0,
                                     {'water': 1}, ions=params.component_ids[-2:], stoichiometry=(1, 1),
                                     molality=0, species=kind, observed=0, form='difference', scale=50)
        if kind != 'activity':
            record.update(fields(reg.predict(model, obs), 'value unit status message convergence branch stable_root_count requested_tolerance_met solver_options'))
    except Exception as error:
        record.update(value=None, status=type(error).__name__, message=str(error))
    record['wall_seconds'] = time.perf_counter() - started
    return serial(record)


def main():
    output = LEAF / 'results'
    output.mkdir(exist_ok=True)
    path = output / 'native-results.json'
    if '--compare-only' in sys.argv:
        run = json.loads(path.read_text())
    else:
        identity = dict(run_id=str(uuid.uuid4()), started_utc=datetime.now(timezone.utc).isoformat(),
                        engine_provenance=serial(epcsaft.provenance()), core_sha256=digest(Path(_core.__file__)),
                        runner_sha256=digest(Path(__file__)), targets_sha256=digest(TARGETS), solvent_sha256=digest(SOLVENT),
                        density_molar_masses_sha256=digest(MASSES))
        started = time.perf_counter()
        cases = declarations()
        records = [calculate(*case, identity) for case in cases]
        # One changed physical declaration isolates the mixed-solvent association rule.
        case, kind, params, x, m = cases[1]
        diagnostic = []
        if records[1]['value'] is not None and abs(records[1]['value'] - 0.415) > 0.0005:
            previous = parameters('3', params.component_ids)
            diagnostic.append(calculate(case, kind, previous, x, m, identity, '3'))
        run = dict(execution_identity=identity, quantities=records, diagnostics=diagnostic,
                   reference_policy=[serial(t) for t in _core.tolerances() if 'REFERENCE' in t.name],
                   wall_seconds=time.perf_counter() - started)
        path.write_text(json.dumps(run, indent=2, allow_nan=False) + '\n')
    rows = list(csv.DictReader(TARGETS.open()))
    assert len(rows) == 13 and len({r['case_id'] for r in rows}) == 13
    for row in rows:
        key = row['case_id'].removesuffix('-literature-average').removesuffix('-paper-model')
        result = next(r for r in run['quantities'] if r['case_id'] == key)
        value = result['value']
        error = value - float(row['source_value']) if value is not None else None
        row.update(computed_value=value, signed_error=error, absolute_error=abs(error) if error is not None else None,
                   outcome=('reproduced' if abs(error) <= float(row['absolute_tolerance']) else 'not reproduced') if error is not None else 'unavailable: ' + (result.get('message') or result['status']),
                   execution_run_id=result['execution_identity']['run_id'], engine_commit=result['execution_identity']['engine_provenance'][0])
    with (output / 'comparisons.csv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({'wall_seconds': run['wall_seconds'], 'comparisons': len(rows), 'run_id': run['execution_identity']['run_id']}))


if __name__ == '__main__':
    main()

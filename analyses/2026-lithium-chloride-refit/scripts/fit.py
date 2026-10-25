"""#640 stage A: frozen declarations and the #575 native Ceres multistart fit."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

from epcsaft import Parameters, regression as r

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
RESULTS = HERE / 'results'
DECL = HERE / 'inputs/fit-declaration.json'
D = json.loads(DECL.read_text())
T, P = D['conditions']['temperature_K'], D['conditions']['pressure_Pa']
COMPONENTS = tuple(D['conditions']['components'])
ION = 'lithium-cation'
LABELS = ['segment_diameter/lithium-cation', 'packing_diameter/lithium-cation',
          'debye_huckel_diameter/lithium-cation', 'born_diameter/lithium-cation',
          'k_ij/lithium-cation+water', 'k_ij/chloride-anion+lithium-cation']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(line for line in stream if not line.startswith('#')))


def write(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def verify():
    assert sha(DECL) == (HERE / 'inputs/declaration.sha256').read_text().split()[0]
    assert sha(ROOT / f"data/packets/{D['base_packet']}/parameters/parameters.json") == D['base_packet_sha256']
    for path, digest in D['inputs_sha256'].items():
        assert sha(ROOT / path) == digest, path


def parameters(physical=None, components=COMPONENTS):
    """The base document with only the six tied Li physical values substituted."""
    document = json.loads((ROOT / f"data/packets/{D['base_packet']}/parameters/parameters.json").read_text())
    if physical is not None:
        values = dict(zip(LABELS, physical, strict=True))
        changed = []
        for owner, coefficients in [(c['component_id'], c.get('coefficients', [])) for c in document['components']] + [
                ('+'.join(sorted((q['component_id_a'], q['component_id_b']))), q['coefficients']) for q in document['pairs']]:
            for coefficient in coefficients:
                label = f"{coefficient['family']}/{owner}"
                if label in values:
                    coefficient['value']['magnitude'] = values[label]
                    changed.append(label)
        assert set(changed) == set(LABELS)
    return Parameters.from_mapping(document, components=components)


def composition(p, molalities):
    amounts = {'water': 1 / p.molar_masses[0], **molalities}
    total = math.fsum(amounts.values())
    return {c: amounts.get(c, 0) / total for c in p.component_ids}


def training(p):
    rows = []
    for spec in D['training']:
        data = read(ROOT / spec['path'])
        assert len(data) == spec['rows']
        for i, row in enumerate(data, 1):
            if spec['source'] == 'Hamer-Wu':
                m, phi = float(row['molality_mol_kg']), float(row['phi'])
                li, na = m, 0.0
                observation = r.electrolyte(p, 'osmotic_coefficient', T, P, {'water': 1},
                                           ions=(ION, 'chloride-anion'), stoichiometry=(1, 1), molality=m,
                                           observed=phi, form='difference', scale=spec['scale'])
                locator = row['locator']
            else:
                li, na, phi = (float(row[k]) for k in ('m_LiCl', 'm_NaCl', 'phi'))
                observation = r.mixture_osmotic_coefficient(p, T, P, composition(p, {
                    ION: li, 'sodium-cation': na, 'chloride-anion': li + na}), observed=phi, scale=spec['scale'])
                locator = D['row_lists']['El Guendouzi'][i - 1]['locator']
            rows.append({'source': spec['source'], 'salt': 'LiCl' if na == 0 else 'LiCl-NaCl', 'row': i, 'locator': locator, 'm_LiCl': li, 'm_NaCl': na,
                         'observed': phi, 'scale': spec['scale'], 'observation': observation})
    return rows


def coordinates(p):
    b = D['bounds']
    diameters = (b['diameter_floor_angstrom'], b['diameter_numerical_upper_angstrom'])
    pairs = (b['k_ij_numerical_lower'], b['k_ij_upper'])
    return [r.coordinate(p, [('segment_diameter', ION, 1), ('packing_diameter', ION, 0.88)], origin=0, scale=1, bounds=diameters),
            r.coordinate(p, [('debye_huckel_diameter', ION, 1)], origin=0, scale=1, bounds=diameters),
            r.coordinate(p, [('born_diameter', ION, 1)], origin=0, scale=1, bounds=diameters),
            r.coordinate(p, [('k_ij', (ION, 'water'), 1)], origin=0, scale=1, bounds=pairs),
            r.coordinate(p, [('k_ij', (ION, 'chloride-anion'), 1)], origin=0, scale=1, bounds=pairs)]


def declarations():
    verify()
    p = parameters()
    assert len(training(p)) == sum(s['rows'] for s in D['training'])
    assert all(len(s['values']) == len(coordinates(p)) for s in D['starts'])


def main():
    verify()
    assert 'EPCSAFT_HSL_LIBRARY' not in os.environ
    assert os.environ['OMP_NUM_THREADS'] == '1'
    assert not (RESULTS / 'fit.json').exists(), 'retain the first frozen run; no implicit reruns'
    prior = load('chloride_fit', ROOT / 'analyses/2026-chloride-brine-refit/scripts/fit.py')
    p = parameters()
    rows, cs, c = training(p), coordinates(p), D['controls']
    began = time.perf_counter()
    results = r.fit(p, cs, [row['observation'] for row in rows],
                    controls=r.FitControls(maximum_iterations=c['maximum_iterations'],
                                           maximum_elapsed_time_seconds=c['maximum_elapsed_time_seconds'], workers=c['workers']),
                    starts=[s['values'] for s in D['starts']], shared_iterations=c['shared_iterations'],
                    shared_seconds=c['shared_seconds'])
    elapsed = time.perf_counter() - began
    records = [{'start': start['id'], **prior.summary(result, LABELS, rows, cs),
                'predictions': list(result.predictions), 'optimizer': list(result.optimizer)}
               for start, result in zip(D['starts'], results, strict=True)]
    found = prior.optima(records)
    best = next((s for s in sorted(records, key=lambda s: s['final_cost']) if s['status'] == 'converged'), None)
    record = {'declaration_sha256': sha(DECL), 'engine_commit': subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip(),
              'engine_module': __import__('epcsaft')._core.__file__, 'python': sys.executable,
              'engine_dirty': bool(subprocess.check_output(['git','status','--porcelain','engine'], text=True).strip()),
              'wall_seconds': elapsed, 'labels': LABELS, 'starts': records, 'optima': found,
              'selected_start': best['start'] if best else None}
    (RESULTS / 'fit.json').write_text(json.dumps(record, indent=2) + '\n')
    for s in records:
        print(s['start'], s['status'], 'cost', s['final_cost'], 'iterations', s['iterations'], 'native failures', s['native_failures'], flush=True)
    assert best is not None, 'no converged start: stop for diagnosis'
    write(RESULTS / 'training-residuals.csv', [
        {**{k:v for k,v in row.items() if k != 'observation'}, 'model': value, 'residual': value-row['observed']}
        for row, value in zip(rows, best['predictions'], strict=True)])
    print('Selected', best['start'], 'cost', best['final_cost'], 'optimizer', best['optimizer'], 'wall seconds', elapsed, flush=True)


if __name__ == '__main__':
    main()

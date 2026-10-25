"""Copy selected prior numerical values and same-host competitor times, with lineage."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

from prepare import BUILD, MAIN

HERE = Path(__file__).resolve().parents[1]
PHYSICAL = ('temperature', 'pressure', 'amounts', 'mole_fractions', 'molar_densities')
REFERENCES, SOURCES = {}, {}


def read(path):
    SOURCES[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return list(csv.DictReader(path.open()))


def add(version, operation, point, quantity, index, value, source, kind, tolerance=1e-10):
    key = (version, operation, point, quantity, int(index))
    REFERENCES[key] = dict(version=version, operation=operation, point=point, quantity=quantity,
        index=int(index), value=float(value), source=str(source), reference_kind=kind, tolerance=tolerance)


def main():
    environments = json.loads((BUILD / 'environments.json').read_text())
    for v in environments:
        version, source = v['version'], Path(v['source'])
        prior = sorted((MAIN / 'build/wheels').glob(f'key-solves-{version}-*/states.csv'))
        if prior:
            path = prior[-1]
            for r in read(path):
                if r['wheel'] == 'B' and r['sample'] == 'warmup' and r['quantity'] in (*PHYSICAL,
                        'predictions', 'gradient', 'optimizer_gradient', 'residual_optimizer_gradient'):
                    add(version, r['case'], r['point'], r['quantity'], r['index'], r['value'], path, 'same revision key-solve return')
        else:
            # The earliest wheel predates the automatic key-solve archive. Its checkout retains the ionic result.
            path = source / 'analyses/2026-tpd-phase-discovery/results/qualification-8470b49a/ionic.csv'
            for r in read(path):
                quantity, _, index = r['quantity'].partition('/')
                if r['case'] == 'case-6/figure-2-point-7/declared' and quantity in PHYSICAL:
                    for mode in ('cold', 'accepted'):
                        op = f'case-6/figure-2-point-7/{mode}'
                        add(version, op, op, quantity, index or 0, r['value'], path, 'earlier return retained in this checkout')

        path = source / 'analyses/2026-competitor-cross-check/results/density-closure-final-values.csv'
        for r in read(path):
            if r['implementation'] != 'Engine' or r['id'] not in ('C1-300', 'C2-0.5', 'C4-298.15', 'C4-323.15'):
                continue
            n = 1 if r['id'] == 'C1-300' else 2
            phase = 0 if r['phase'] in ('liquid', 'water-rich') else 1
            quantity = {'x': 'mole_fractions', 'rho': 'molar_densities', 'pressure': 'pressure'}.get(r['quantity'])
            if quantity is None:
                continue
            index = phase * n + int(r['component']) if quantity == 'mole_fractions' else phase
            add(version, r['id'], r['id'], quantity, index, r['value'], path,
                'prior cross-check return retained in this checkout', 1e-8)

        path = source / 'analyses/2026-binary-methane-ethane-bubble-point/results/values.csv'
        r = next(r for r in read(path) if r['implementation'] == 'ePC-SAFT')
        for quantity, fields in [('pressure', ('pressure_Pa', 'pressure_Pa')),
                ('molar_densities', ('rho_liquid_mol_m3', 'rho_vapor_mol_m3')),
                ('mole_fractions', ('liquid_x_methane', 'liquid_x_ethane', 'vapor_y_methane', 'vapor_y_ethane'))]:
            for i, field in enumerate(fields):
                add(version, 'methane-bubble', 'methane-bubble', quantity, i, r[field], path,
                    'prior bubble return retained in this checkout', 1e-8)

        path = source / 'analyses/2026-binary-methane-ethane-tp-flash/results/performance-current.csv'
        r = next(r for r in read(path) if r['case'] == '0' and r['mode'] == 'cold')
        for quantity, values in [('pressure', [r['P_Pa'], r['P_Pa']]),
                ('molar_densities', [r['rho_L'], r['rho_V']]),
                ('mole_fractions', [float(r['x_methane']), 1-float(r['x_methane']), float(r['y_methane']), 1-float(r['y_methane'])]),
                ('amounts', [1-float(r['beta']), float(r['beta'])])]:
            for i, value in enumerate(values):
                add(version, 'methane-flash', 'methane-flash', quantity, i, value, path,
                    'prior flash return retained in this checkout', 1e-8)

    # The fixed F1 declaration was captured before the latest optimizations; compare older wheels to that actual return.
    path = MAIN / 'build/issue-533/source/mea-F1-historical-values.json'
    SOURCES[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    f1 = json.loads(path.read_text())
    packet = gzip.decompress((HERE.parents[1] / 'tools/fixtures/mea-f1.json.gz').read_bytes())
    assert f1['packet_sha256'] == hashlib.sha256(packet).hexdigest(), 'historical F1 declaration differs'
    for v in environments:
        for quantity in ('predictions', 'gradient', 'optimizer_gradient', 'residual_optimizer_gradient'):
            for i, value in enumerate(f1[quantity]):
                key = (v['version'], 'mea-f1/residual-jacobian', 'mea-f1/residual-jacobian', quantity, i)
                if key not in REFERENCES:
                    add(*key, value, path, 'F1 captured on 49aac619; same frozen packet')

    references = HERE / 'results/retained-values.csv.gz'
    with gzip.open(references, 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(next(iter(REFERENCES.values()))), lineterminator='\n')
        writer.writeheader()
        writer.writerows(REFERENCES.values())
    path = MAIN / 'build/issue-533/selected-row-statistics.csv'
    selected = read(path)
    identities = {}
    for library in ('Clapeyron', 'FeOs'):
        identity_path = MAIN / f'build/issue-533/identity/{library}.json'
        SOURCES[str(identity_path)] = hashlib.sha256(identity_path.read_bytes()).hexdigest()
        identities[library] = json.loads(identity_path.read_text())
    competitors = []
    for op in ('C1-300', 'C2-0.5', 'methane-bubble', 'C4-298.15', 'C4-323.15', 'methane-flash', 'A-propane-liquid'):
        r = next(r for r in selected if r['id'] == op)
        for library in ('Clapeyron', 'FeOs'):
            identity = identities[library].get('identity', identities[library])
            competitors.append(dict(operation=op, library=library, version=identity['version'],
                seconds=float(r[f'{library}_median_us'])/1e6,
                conditions=r['timing_conditions'], source=str(path), samples=int(float(r[f'{library}_n']))))
    with (HERE / 'results/competitors.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(competitors[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(competitors)
    path = HERE.parents[1] / 'analyses/2026-audit-timing-campaign/results/statistics.csv'
    audit = read(path)
    extension = json.loads((HERE / 'inputs/extension.json').read_text())
    matches = {op: op for op in ('C1-300', 'C2-0.5', 'methane-bubble')}
    matches.update({op: f"methane_ethane/{q['T_K']}-{q['P_Pa']}/discover"
                    for op, q in extension.items() if q['kind'] == 'discover'})
    retained = []
    for op, prior_op in matches.items():
        for r in audit:
            if r['operation'] == prior_op and r['tool'] in ('FeOs', 'Clapeyron'):
                retained.append(dict(operation=op, library=r['tool'], seconds=r['median_s'],
                    source=str(path), source_operation=prior_op, samples=r['samples'],
                    load_min=r['load1_min'], load_max=r['load1_max']))
    with (HERE / 'results/competitors-578.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(retained[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(retained)
    (HERE / 'results/competitor-identities.json').write_text(json.dumps(identities, indent=2) + '\n')
    (HERE / 'results/prior-source-hashes.json').write_text(json.dumps(SOURCES, indent=2) + '\n')
    print(f'Retained {len(REFERENCES)} prior values and {len(competitors)} competitor times from {len(SOURCES)} files')


if __name__ == '__main__':
    main()

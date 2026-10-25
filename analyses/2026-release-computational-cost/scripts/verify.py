"""Value agreement before timing; no thermodynamic calculations or historical answers."""
import csv
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
TOL = json.loads((HERE / 'inputs/inputs.json').read_text())['tolerances']


def flat(value, prefix=''):
    if isinstance(value, dict):
        return {key: v for k, item in value.items() for key, v in flat(item, f'{prefix}/{k}').items()}
    if isinstance(value, list):
        return {key: v for k, item in enumerate(value) for key, v in flat(item, f'{prefix}/{k}').items()}
    return {prefix: value}


def values(operation, tool, sample):
    v = sample['value']
    if operation == 'property-bundle' or '/actions-' in operation:
        return v
    if operation == 'pressure-state':
        assert len(v) == 32
        return {str(i): dict(rho=row['rho'] if tool == 'Engine' else row['rho_mol_m3']) for i, row in enumerate(v)}
    if operation == 'electrolyte-eight-feeds':
        assert len(v) == 8
        for r in v:
            assert r['success'] and r['criterion_met'] and r['phases'] == 2, r
        return {str(r['point']): {k: r[k] for k in ('amounts', 'mole_fractions', 'molar_densities', 'max_residual')} for r in v}
    if operation == 'explicit-discovery':
        if tool == 'Engine':
            r = v[0]
            assert r['success'] and r['criterion_met'] and r['phases'] == 2, r
            assert r['discovery']['attempts'] > 0
            rho = r['molar_densities']
            xs = [r['mole_fractions'][i * 2:i * 2 + 2] for i in range(2)]
            amounts = r['amounts']
            fraction = [a / sum(amounts) for a in amounts]
        elif tool == 'FeOs':
            assert v['phases'] == 2, v
            rho, xs = v['densities'], v['compositions']
            fraction = [1 - v['vapor_fraction'], v['vapor_fraction']]
        else:
            assert v['phases'] == 2, v
            rho, xs, fraction = [1 / a for a in v['volumes']], v['compositions'], v['fractions']
        order = sorted(range(2), key=lambda j: -rho[j])
        return dict(rho=[rho[j] for j in order], x=[xs[j] for j in order], fractions=[fraction[j] for j in order])
    if tool in ('Engine', 'FeOs'):
        phases = sorted(v['phases'], key=lambda r: -r['rho'])
        result = dict(P=v['P'], rho=[r['rho'] for r in phases], x=[r['x'] for r in phases])
    else:
        order = sorted(range(2), key=lambda j: -v['densities'][j])
        result = dict(P=v['P'], rho=[v['densities'][j] for j in order], x=[v['compositions'][j] for j in order])
    if operation.startswith('C4'):
        result['beta'] = v['beta']
    return result


def compare(operation, tool, reference, other):
    a, b = flat(reference), flat(other)
    assert a.keys() == b.keys(), (a.keys(), b.keys())
    rows = []
    for key, x in a.items():
        y = b[key]
        if '/actions-' in operation:
            absolute, relative = TOL['action_abs'], TOL['action_rel']
        elif '/x/' in key or key.startswith('/fractions/') or key == '/beta':
            absolute, relative = TOL['composition_abs'], 0
        else:
            relative = TOL['property_rel'] if operation == 'property-bundle' else TOL['density_rel'] if operation == 'pressure-state' else TOL['equilibrium_rel']
            absolute = relative if operation == 'property-bundle' else 0
        error = abs(x - y) if x is not None and y is not None else None
        allowed = (relative * max(abs(x), abs(y), 1.) if operation == 'property-bundle' else absolute + relative * max(abs(x), abs(y))) if error is not None else 0
        rows.append(dict(operation=operation, comparator=tool, quantity=key, package_value=x, comparator_value=y,
                         absolute_error=error, relative_error=error / max(abs(x), abs(y)) if error is not None and max(abs(x), abs(y)) else 0,
                         absolute_tolerance=absolute, relative_tolerance=relative, allowed_error=allowed,
                         passed=error is not None and math.isfinite(error) and error <= allowed, cause=''))
    assert rows
    return rows


def main():
    d = json.loads((HERE / 'results/preflight.json').read_text())
    rows, normalized = [], {}
    for operation, arms in d['operations'].items():
        normalized[operation] = {}
        for tool, reply in arms.items():
            try:
                normalized[operation][tool] = values(operation, tool, reply['samples'][0])
            except Exception as error:
                rows.append(dict(operation=operation, comparator=tool, quantity='complete numerical return', passed=False,
                                 cause=reply.get('error') or f'{type(error).__name__}: {error}'))
        if operation.endswith('/actions-fd'):
            exact_op = operation.replace('/actions-fd', '/actions-exact')
            try:
                rows += compare(operation, 'package central differences', normalized[exact_op]['Engine'], normalized[operation]['Engine'])
            except Exception as error:
                rows.append(dict(operation=operation, comparator='package central differences', passed=False, cause=str(error)))
        elif not '/actions-' in operation and operation != 'electrolyte-eight-feeds':
            for tool, value in normalized[operation].items():
                if tool != 'Engine':
                    rows += compare(operation, tool, normalized[operation]['Engine'], value)
        elif operation == 'electrolyte-eight-feeds':
            feeds = {r.pop('point'): [float(v) for v in r.values()] for r in csv.DictReader((HERE / 'inputs/electrolyte-feeds.csv').open())}
            for point, r in normalized[operation]['Engine'].items():
                n, x = r['amounts'], r['mole_fractions']
                for i, feed in enumerate(feeds[point]):
                    value = sum(n[j] * x[5 * j + i] for j in range(2))
                    rows.append(dict(operation=operation, comparator='component balance', quantity=f'feed-{point}/component-{i}',
                                     package_value=value, comparator_value=feed, absolute_error=abs(value - feed), absolute_tolerance=1e-10,
                                     passed=abs(value - feed) <= 1e-10, cause=''))
                for j in range(2):
                    for name, value in [('mole-fraction sum', sum(x[5*j:5*j+5]) - 1), ('charge', x[5*j+3] - x[5*j+4])]:
                        rows.append(dict(operation=operation, comparator=name, quantity=f'feed-{point}/phase-{j}',
                                         package_value=value, comparator_value=0, absolute_error=abs(value), absolute_tolerance=1e-10,
                                         passed=abs(value) <= 1e-10, cause=''))
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (HERE / 'results/checks.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    (HERE / 'results/preflight-values.json').write_text(json.dumps(normalized, indent=2) + '\n')
    print('value checks:', sum(r['passed'] for r in rows), 'passed of', len(rows))
    for r in rows:
        if not r['passed']:
            print(r)


if __name__ == '__main__':
    main()

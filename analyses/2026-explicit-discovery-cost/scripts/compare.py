"""Compare frozen #855 states, trial outcomes, causes and native numerical rows.

The ordinary scalar owner is #828/check-equilibrium-corpora: relative 1e-10,
absolute 1e-12. Raw charts and iteration diagnostics remain in the source records.
"""
import copy
import csv
import itertools
import json
import math
from pathlib import Path
import re
import sys

TASK, OUT = map(Path, sys.argv[1:3])
IGNORED = {'work', 'iterations', 'evaluations', 'hessian_calls', 'coordinates', 'raw_coordinates',
           'raw_residuals', 'candidate_coordinates', 'declaration_attempt_coordinates', 'split_coordinates',
           'residual_over_criterion', 'residual_over_floor', 'amplification', 'contraction', 'wall_s'}
NUMBER = re.compile(r'[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?')
rows, failures, counts = [], [], {}


def scalar(a, b, path, family):
    allowed = max(1e-12, 1e-10 * max(abs(a), abs(b)))
    scaled = 0.0 if a == b else abs(a - b) / allowed
    ok = a == b or math.isfinite(scaled) and scaled <= 1.0
    rows.append(dict(family=family, quantity=path, A=a, B=b, allowance=allowed, error_over_allowance=scaled, passed=ok))
    if not ok:
        failures.append(dict(family=family, quantity=path, A=a, B=b, error_over_allowance=scaled))


def compare(a, b, path, family):
    if isinstance(a, dict) and isinstance(b, dict):
        a, b = copy.deepcopy(a), copy.deepcopy(b)
        if a.get('type') == b.get('type') == 'EquilibriumResult' and a['phase_count'] == b['phase_count'] and a['phase_count'] > 0:
            p = a['phase_count']
            if len(a['mole_fractions']) == len(b['mole_fractions']) and a['mole_fractions']:
                n = len(a['mole_fractions']) // p
                phases = [v.get('discovery', {}).get('declaration', {}).get('phases', []) if v.get('discovery') else [] for v in (a, b)]
                def admissible(order):
                    return not all(phases) or all((phases[0][j]['kind'], phases[0][j]['support']) ==
                        (phases[1][order[j]]['kind'], phases[1][order[j]]['support']) for j in range(p))
                orders = [o for o in itertools.permutations(range(p)) if admissible(o)]
                if not orders:
                    failures.append(dict(family=family, quantity=path, reason='phase support/kind mismatch')); return
                order = min(orders, key=lambda o: max(abs(a['mole_fractions'][j*n+i] - b['mole_fractions'][o[j]*n+i]) for j in range(p) for i in range(n)))
                for key in ('amounts', 'molar_densities', 'packing_fractions', 'mole_fractions'):
                    width = n if key == 'mole_fractions' else 1
                    b[key] = [b[key][j*width+i] for j in order for i in range(width)]
                if order != tuple(range(p)):
                    failures.append(dict(family=family, quantity=path, reason='phase permutation requires defining-row remapping', permutation=order))
        keys = (a.keys() | b.keys()) - IGNORED
        for k in sorted(keys):
            if k not in a or k not in b:
                failures.append(dict(family=family, quantity=path+'/'+k, reason='missing field')); continue
            compare(a[k], b[k], path+'/'+k, family)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            failures.append(dict(family=family, quantity=path, reason='changed population', A=len(a), B=len(b))); return
        if a and all(isinstance(v, dict) and v.get('type') == 'StabilityTrial' for v in a+b):
            identity = lambda t: (t['seed'], t['phase']['kind'], tuple(t['phase']['support']))
            a, b = sorted(a, key=identity), sorted(b, key=identity)
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, path+'/'+str(i), family)
    elif type(a) in (float, int) and type(b) in (float, int):
        scalar(a, b, path, family)
    elif a != b:
        if path.endswith('/message') and isinstance(a, str) and isinstance(b, str):
            # Successful iteration counts are diagnostic; retain every other causal operand.
            if 'converged:' in a and 'converged:' in b:
                a, b = [re.sub(r'after \d+ iterations', 'after <iterations> iterations', s) for s in (a, b)]
            if NUMBER.sub('<number>', a) == NUMBER.sub('<number>', b):
                for i, (x, y) in enumerate(zip(NUMBER.findall(a), NUMBER.findall(b))):
                    scalar(float(x), float(y), path+'/operand-'+str(i), family)
                return
        failures.append(dict(family=family, quantity=path, reason='changed outcome or cause', A=a, B=b))


def records(arm):
    result = {}
    for p in (TASK / ('inventory-'+arm)).glob('*.json'):
        r = json.loads(p.read_text())
        if 'key' in r:
            assert r['key'] not in result
            result[r['key']] = r
    return result


aa, bb = records('A'), records('B')
inventory = json.loads((TASK/'discovery-inventory.json').read_text())['cases']
expected = {r['key']+'/'+role for r in inventory for role in r['roles']}
assert set(aa) == set(bb) == expected and len(inventory) == 51
for key in sorted(expected):
    assert aa[key]['declaration'] == bb[key]['declaration']
    compare(aa[key]['value'], bb[key]['value'], key, 'primary')
counts['primary'] = len(inventory)
affected = json.loads((TASK/'affected-check-inventory.json').read_text())['cases']
for case in affected:
    key, family = case['key'], case['runner']
    if family == 'native':
        def native(arm):
            lines = (TASK/('affected-'+arm)/(key+'.log')).read_text().splitlines()
            assert any('EPCSAFT_CASE schema_version=1 test_id='+key+' execution=Passed' in s for s in lines)
            rr = [json.loads(s.split(' ', 1)[1]) for s in lines if s.startswith('EPCSAFT_COMPARE_JSON ')]
            if not rr:
                assert any('actual=' in s and 'reference=' in s for s in lines)
                return [dict(text=NUMBER.sub('<number>', s), values=list(map(float, NUMBER.findall(s)))) for s in lines]
            assert len({r['comparison_id'] for r in rr}) == len(rr)
            assert all(r['verdict'] if 'verdict' in r else r['judgment'] == 'Observation only' for r in rr)
            return {r['comparison_id']:r for r in rr}
        compare(native('A'), native('B'), key, family)
    else:
        rr = [json.loads((TASK/case['output'].replace('{arm}', arm)).read_text()) for arm in ('A', 'B')]
        for r in rr:
            for k in ('engine_commit', 'engine_module', 'python', 'trace'): r.pop(k, None)
        compare(*rr, key, family)
    counts[family] = counts.get(family, 0)+1
OUT.mkdir(parents=True, exist_ok=True)
with (OUT/'comparison-values.csv').open('w') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
summary = [dict(family=k, entries=n, comparisons=sum(r['family']==k for r in rows),
    max_error_over_allowance=max((r['error_over_allowance'] for r in rows if r['family']==k), default=0),
    changed_outcomes=sum(r['family']==k for r in failures)) for k,n in counts.items()]
(OUT/'equivalence.json').write_text(json.dumps(dict(summary=summary, failures=failures), indent=2)+'\n')
print(json.dumps(summary, indent=2)); print('Failures:', len(failures))
sys.exit(bool(failures))

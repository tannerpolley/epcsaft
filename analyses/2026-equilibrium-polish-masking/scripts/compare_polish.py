"""Compare retained states in the owner-approved coordinates (#272, amendment 2)."""
import csv
import json
from qualify_polish import ROOT, OUT, save_csv


def state(label):
    values = {}
    owner = json.loads((OUT / f'{label}-wheel.json').read_text())['states_path']
    for row in csv.DictReader((ROOT / owner).open()):
        values.setdefault(row['case'], {}).setdefault(row['quantity'], []).append(row['value'])
    return values


base, candidate = state('baseline'), state('candidate')
epsilon = json.loads((OUT / 'candidate-policy.json').read_text())['EQUILIBRIUM_POLISH_CORRECTION']
inputs = json.loads((OUT / 'inputs.json').read_text())
traces = list(csv.DictReader((OUT / 'candidate-trace.csv').open()))
comparisons, work = [], []
def compare_state(case, b, c, declaration):
    result = []
    def compare(budget, quantity, index, x, y, scale):
        difference = abs(x-y) / abs(scale) if scale else (0. if x == y else float('inf'))
        result.append(dict(case=case, budget=budget, quantity=quantity, index=index,
                                baseline=x, candidate=y, scale=scale, difference=difference))
    for quantity in ('temperature', 'pressure', 'molar_densities'):
        for index, (x, y) in enumerate(zip(b[quantity], c[quantity], strict=True)):
            compare('observable-relative', quantity, index, float(x), float(y), float(x))
    feed = declaration['feed']
    total = sum(feed.values())
    offset = 0
    for phase, specification in enumerate(declaration['phases']):
        amount, other = float(b['amounts'][phase]), float(c['amounts'][phase])
        compare('phase-relative', 'phase_amount', phase, amount, other, amount)
        for species in specification['support'] or declaration['components']:
            x = amount * float(b['mole_fractions'][offset])
            y = other * float(c['mole_fractions'][offset])
            trace = feed.get(species, 0.) < 1e-6 * total
            compare('trace-relative' if trace else 'phase-relative', f'amount/{phase}/{species}',
                    offset, x, y, x if trace else amount)
            offset += 1
    return result


for case, b in base.items():
    c, declaration = candidate[case], inputs[case]
    comparisons.extend(compare_state(case, b, c, declaration))
    directions = [r for r in traces if r['case'] == case and r['event'] == 'polish_correction']
    stops = [r for r in directions if float(r['v8']) == 1.]
    chosen = stops[-1] if stops else directions[-1] if directions else None
    row = dict(case=case)
    for label, values in [('baseline', b), ('candidate', c)]:
        for key in ('iterations', 'evaluations', 'hessian_calls'):
            row[f'{label}_{key}'] = int(values[key][0])
    row.update(fired=bool(stops), controlling_component=chosen['component'] if chosen else 'none',
               scaled_correction=float(chosen['v0']) if chosen else 0.,
               feed_relative_amount=float(chosen['v2']) if chosen else 0.,
               forward_bound=float(chosen['v3']) if chosen else 0., step=int(chosen['step']) if chosen else -1)
    work.append(row)
# The timing tool also compares common accepted starts. Apply the same budgets.
if (OUT / 'timing/states.csv').exists():
    timed = {}
    for r in csv.DictReader((OUT / 'timing/states.csv').open()):
        key = (r['case'], r['sample'], r['point'])
        timed.setdefault(key, {}).setdefault(r['wheel'], {}).setdefault(r['quantity'], []).append(r['value'])
    for (case, sample, point), wheels in timed.items():
        declaration = inputs[point if point.startswith('khudaida/') else
                             'mea-nominally-unloaded' if case.startswith('mea-nominally-unloaded/') else
                             case.rsplit('/', 1)[0]]
        comparisons.extend(compare_state(f'timing/{case}/{sample}/{point}', wheels['A'], wheels['B'], declaration))
save_csv(OUT / 'physical-state-comparison.csv', comparisons)
save_csv(OUT / 'work.csv', work)
maxima = [max((r for r in comparisons if r['budget'] == budget), key=lambda r: r['difference'])
          for budget in ('observable-relative', 'trace-relative', 'phase-relative')]
save_csv(OUT / 'budget-maxima.csv', maxima)
for row in maxima: print(row)
assert all(r['difference'] <= 10 * epsilon for r in comparisons)
assert not next(r for r in work if r['case'] == 'masking')['fired']
small = [
    r for r in traces if r['case'] == 'small-positive-phase' and r['event'] == 'polish_correction']
assert any(r['component'] == 'amount/V' and float(r['v0']) > epsilon
           and float(r['v2']) <= epsilon and float(r['v8']) == 0. for r in small)
ionic = next(r for r in work if r['case'] == 'case-6/figure-2-point-7')
assert ionic['candidate_evaluations'] <= 49
mea = next(r for r in work if r['case'] == 'case-7/flash-5918.5286-Pa')
assert all(mea[f'candidate_{k}'] <= limit for k, limit in
           zip(('iterations', 'evaluations', 'hessian_calls'), (43, 57, 44)))
kh = [r for r in work if r['case'].startswith('khudaida/')]
assert len(kh) == 8
assert all(sum(r[f'candidate_{k}'] for r in kh) <= limit for k, limit in
           zip(('iterations', 'evaluations', 'hessian_calls'), (66, 119, 74)))
assert any(r['case'] == 'mea-nominally-unloaded' and r['event'] == 'polish_trial'
           and r['step'] == '3' and float(r['v7']) == 1. for r in traces)
for row in work: print(row)

"""#243 D2: quantity-based comparison of two baseline.py replays, all four sections.
usage: replay_compare.py BEFORE_DIR AFTER_DIR OUT.json
Work counts and equilibrium row residuals are reported, not compared; every result that reports
requested_tolerance_met True must have its residuals within its retained criterion (ASL central states: the
anchor's). Keys must be unique, sections nonempty, kinds equal and floats finite; every other text/int/bool
row must be equal. Fitting residuals are observation residuals and compared as values."""
import csv, json, math, re, sys
from collections import defaultdict
from pathlib import Path

WORK = {'iterations', 'evaluations', 'hessian_calls', 'residual_evaluations', 'jacobian_evaluations', 'iteration_costs',
        'state_factorizations', 'linear_solves', 'steps', 'trial_failures'}
DIAGNOSTIC = {'residuals', 'classified_residuals', 'residual', 'floor', 'amplification', 'condition_estimate',
              'maximum_linear_residual', 'central_residual', 'central_floor', 'classified_residual', 'water_gap'}


def base(quantity):
    return re.sub(r'/\d+', '', quantity)


def load(path, failures):
    rows = {}
    for r in csv.DictReader(path.open()):
        key = (r['case'], r['quantity'])
        if key in rows:
            failures.append(dict(key=key, reason=f'duplicate key in {path}'))
        rows[key] = r
    if not rows:
        failures.append(dict(key=None, reason=f'no rows in {path}'))
    return rows


def unmet(rows, label, failures):
    """Results claiming requested_tolerance_met whose residual exceeds the retained criterion."""
    classified = defaultdict(list)
    for (case, q), r in rows.items():
        if base(q) in ('classified_residuals', 'classified_residual'):
            classified[case].append(r['value'])
    for (case, q), r in rows.items():
        if q != 'requested_tolerance_met' or r['value'] != 'True':
            continue
        anchor = 'asl/anchor/declared' if case.endswith('/diagnostics') else case
        criterion = rows.get((anchor + '/convergence', 'criterion'), {}).get('value')
        residuals = classified[case] + [rows.get((case + '/convergence', 'residual'), {}).get('value', '0.')]
        try:
            limit, worst = float(criterion), max(abs(float(v)) for v in residuals)
        except (TypeError, ValueError):
            limit, worst = math.nan, math.nan
        if not worst <= limit:
            failures.append(dict(key=(case, q), reason=f'{label}: requested tolerance met but residual {worst:.3g} > {limit:.3g}'))


def main():
    before_dir, after_dir, out = map(Path, sys.argv[1:4])
    report = {}
    for section in ('declared', 'continuation', 'fitting', 'asl'):
        failures = []
        before, after = load(before_dir / f'{section}.csv', failures), load(after_dir / f'{section}.csv', failures)
        unmet(before, 'before', failures)
        unmet(after, 'after', failures)
        totals = defaultdict(float)  # per case-state: feed total from the before amounts
        for (case, q), r in before.items():
            if base(q) == 'amounts' and r['kind'] == 'float':
                totals[case] += abs(float(r['value']))
        work, maxima, identical = [], defaultdict(float), 0
        for key in sorted(set(before) | set(after)):
            a, b = before.get(key), after.get(key)
            if a is None or b is None:
                failures.append(dict(key=key, reason='row missing', before=a and a['value'], after=b and b['value']))
                continue
            if a['kind'] != b['kind']:
                failures.append(dict(key=key, reason='kind changed', before=a['kind'], after=b['kind']))
                continue
            if a['kind'] == 'float' and not (math.isfinite(float(a['value'])) and math.isfinite(float(b['value']))):
                failures.append(dict(key=key, reason='non-finite value', before=a['value'], after=b['value']))
                continue
            if a['value'] == b['value']:
                identical += 1
                continue
            case, q = key
            name = base(q).rsplit('/', 1)[-1]
            if name in WORK:
                work.append(dict(key=key, before=a['value'], after=b['value']))
                continue
            if a['kind'] != 'float':
                failures.append(dict(key=key, reason='non-numeric row changed', before=a['value'], after=b['value']))
                continue
            x, y = float(a['value']), float(b['value'])
            if name in DIAGNOSTIC and section != 'fitting':
                maxima['diagnostic:' + name] = max(maxima['diagnostic:' + name], abs(y))
                continue
            N = max(1., totals[case.split('/row/')[0]])
            if section == 'asl' and name in ('amounts', 'gradient', 'hessian'):
                error, limit = abs(x - y), 1e-11 + 1e-9 * max(abs(x), abs(y))  # run_live_settler.py transport criterion
            elif name == 'mole_fractions':
                error, limit = abs(x - y), 1e-9
            elif name == 'amounts':
                error, limit = abs(x - y) / N, 1e-12
            elif name in ('molar_densities', 'packing_fractions', 'temperature', 'pressure', 'path_values',
                          'reached_path', 'attempted_path'):
                error, limit = abs(x / y - 1.) if y else abs(x), 1e-9
            elif name == 'objective':
                error, limit = abs(x - y) / N / max(1., abs(y) / N), 1e-10
            elif name == 'reduced_hessian_minimum':
                error, limit = (0. if (x > 0) == (y > 0) else 1.), 0.5  # inertia/admission meaning, not a feed identity
            else:  # multipliers, coordinates, actions and fitting values: reference-or-one scale
                error, limit = abs(x - y) / max(1., abs(x), abs(y)), 1e-9
            maxima[name] = max(maxima[name], error)
            if not (math.isfinite(error) and error <= limit):
                failures.append(dict(key=key, reason=f'{name} {error:.3g} > {limit:.3g}', before=x, after=y))
        report[section] = dict(rows=len(before), identical=identical, failures=failures, work=work, maxima=dict(maxima))
        print(section, len(before), 'identical', identical, 'failures', len(failures), 'work rows changed', len(work))
        for f in failures[:20]:
            print('  ', f)
    Path(out).write_text(json.dumps(report, indent=2, default=str) + '\n')


if __name__ == '__main__':
    main()

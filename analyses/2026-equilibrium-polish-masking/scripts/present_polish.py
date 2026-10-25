"""Format retained polish comparisons for Quarto; no Engine imports or calculations."""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'analyses/2026-equilibrium-polish-masking/results/resumed'
lines = []


def rows(name):
    return list(csv.DictReader((OUT / name).open()))


def table(title, headers, values):
    lines.extend(['### ' + title, '', '| ' + ' | '.join(headers) + ' |',
                  '| ' + ' | '.join(['---'] * len(headers)) + ' |'])
    lines.extend('| ' + ' | '.join(str(x) for x in value) + ' |' for value in values)
    lines.append('')


policy = json.loads((OUT / 'candidate-policy.json').read_text())
metadata = json.loads((OUT / 'qualification.json').read_text())
inputs = json.loads((OUT / 'inputs.json').read_text())
work = rows('work.csv')
traces = rows('candidate-trace.csv')
small = next(r for r in work if r['case'] == 'small-positive-phase')
fourth = next(r for r in traces if r['case'] == 'mea-nominally-unloaded'
              and r['event'] == 'polish_trial' and r['step'] == '3' and float(r['v7']) == 1)
epsilon = policy['EQUILIBRIUM_POLISH_CORRECTION']
mask = inputs['masking']
historical = list(csv.DictReader((OUT.parent / 'physical-state-comparison.csv').open()))
small_input = json.loads((OUT.parent / 'amendment-2/small-phase-input.json').read_text())
table('Basis and falsifiers', ['Quantity', 'Unit or basis', 'Value'], [
    ('Local correction tolerance', 'scaled chart', epsilon), ('Each preservation budget', 'relative', 10 * epsilon),
    ('Declared qualification cases', 'count', len(work)), ('Paired timing states', 'count', metadata['timed_state_pairs']),
    ('Masking temperature', 'K', mask['T']), ('Masking pressure', 'Pa', mask['P']),
    ('Masking trace feed', 'mol', mask['feed']['s1']), ('Masking row tolerance', 'classified rows', mask['row_tolerance']),
    ('Small vapor nominal amount/feed', 'mol/mol', small_input['beta']),
    ('Topology cutoff', 'amount/feed', policy['STABILITY_PHASE_AMOUNT_FLOOR']),
    ('Small vapor correction', 'phase-relative', small['scaled_correction']),
    ('Same amount correction', 'feed-relative', small['feed_relative_amount']),
    ('Controlling coordinate', 'chart', small['controlling_component']),
    ('Unloaded fourth undamped correction', 'infinity norm', fourth['v0']),
    ('Rejected norm-band maximum state change', 'own-relative', max(float(r['relative_change']) for r in historical))])
table('Maximum preserved-state differences', ['Budget', 'Case', 'Quantity', 'Maximum difference'],
      [(r['budget'], r['case'], r['quantity'], f"{float(r['difference']):.6g}") for r in rows('budget-maxima.csv')])
def counts(r, label):
    return '/'.join(r[f'{label}_{k}'] for k in ('iterations', 'evaluations', 'hessian_calls'))
table('Native solve work', ['Case', 'Unchanged polish', 'Candidate', 'Shortcut fired'],
      [(r['case'], counts(r, 'baseline'), counts(r, 'candidate'), r['fired']) for r in work])
kh = [r for r in work if r['case'].startswith('khudaida/')]
table('Sequential solve total', ['Sequence', 'Unchanged polish', 'Candidate', 'Solved points'],
      [('Khudaida', *('/'.join(str(sum(int(r[f'{label}_{k}']) for r in kh))
          for k in ('iterations', 'evaluations', 'hessian_calls')) for label in ('baseline', 'candidate')), len(kh))])
scaling = rows('amount-scaling.csv')[0]
table('Small-phase full-coordinate scaling replay', ['Quantity', 'Value'],
      [(k.replace('_', ' '), v) for k, v in scaling.items()])
table('Native solve-admission falsifier', ['Quantity', 'Unit', 'Computed', 'Reference', 'Relative tolerance', 'Passed'],
      [(r['quantity'], r['unit'], r['model'], r['reference'], r['rtol'], r['verdict']) for r in rows('admission-falsifier.csv')])
table('Unchanged checks that fail on the candidate',
      ['Comparison', 'Candidate', 'Reference', 'Absolute error', 'Allowed error', 'Passes on main'],
      [(r['comparison_id'], r['candidate'], r['reference'], r['absolute_error'], r['allowed_error'],
        r['baseline_check_passed']) for r in rows('gate-findings.csv')])
kkt = next(r for r in rows('baseline-kkt-values.csv') if r['comparison_id'].endswith('/kkt_residual'))
stop = next(r for r in rows('candidate-kkt-trace.csv') if r['event'] == 'polish_correction' and r['v8'] == '1.0')
table('KKT correction admission at the failed point', ['Quantity', 'Value'],
      [('Baseline KKT residual', kkt['model']), ('Candidate forward bound', stop['v3']),
       ('Classified residual', stop['v4']), ('Full-rank factorization rank', stop['v7']),
       ('Current best', stop['v5']), ('Positive free amount', stop['v6']), ('Polish step', stop['step'])])
timing = json.loads((OUT / 'timing/summary.json').read_text())
table('Paired solve time', ['Case', 'Reference median (s)', 'Candidate median (s)', 'Median change (%)',
       'Observed range (%)', 'Median absolute deviation (%)'],
      [(case, f"{r['median_A_seconds']:.6g}", f"{r['median_B_seconds']:.6g}", f"{r['median_change_percent']:+.3f}",
        f"[{r['range_percent'][0]:+.3f}, {r['range_percent'][1]:+.3f}]", f"{r['mad_percent']:.3f}")
       for case, r in timing['cases'].items()])
table('Timing controls', ['Quantity', 'Value'], [('Counterbalanced pairs', timing['pairs']),
      ('Numerical threads', timing['threads']), ('CPU', timing['cpu']), ('Niceness', timing['nice'])])
table('Source and wheel identities', ['Quantity', 'Value'], [('Reference source', metadata['reference_source']),
      ('Reference states owner', metadata['reference_states_owner']),
      *[(f'{label} wheel SHA-256', json.loads((OUT / f'{label}-wheel.json').read_text())['sha256'])
        for label in ('baseline', 'candidate')]])
destination = ROOT / 'build/notebook-data/polish.md'
destination.parent.mkdir(parents=True, exist_ok=True)
destination.write_text('\n'.join(lines))

# Read the unchanged-polish diagnosis independently of the earlier candidate.
OUT = OUT.parent / 'damped-stop'
lines = []
trials = rows('trials.csv')
predicates = rows('predicates.csv')
identity = json.loads((OUT / 'provenance.json').read_text())
table('Recording verification', ['Quantity', 'Value'], [
    ('Unchanged main source', identity['main']),
    ('Recorded trials', len(trials)),
    ('Returned states compared exactly', identity['states_compared']),
    ('Exact equality of all recorded states', identity['pristine_states_identical']),
    ('Productive damped trials in the fit', sum(r['case'] == 'fit' and int(r['accepted'])
        and float(r['alpha']) < 1 for r in trials)),
    *[(f"{k} preservation budget", v) for k, v in identity['preservation_budgets'].items()]])
table('Predicate falsification', ['Predicate', 'Ionic rejected caught', 'Ionic productive caught',
      'Khudaida rejected caught', 'Khudaida productive caught', 'Protected productive caught'], [
    (name, *[sum(int(r[key]) for r in predicates if r['predicate'] == name and selected(r['case']))
     for selected, key in ((lambda c: c == 'case-6_figure-2-point-7', 'rejected_caught'),
                          (lambda c: c == 'case-6_figure-2-point-7', 'productive_caught'),
                          (lambda c: c.startswith('khudaida'), 'rejected_caught'),
                          (lambda c: c.startswith('khudaida'), 'productive_caught'),
                          (lambda c: c in ('isomerization', 'caloric-action', 'two-phase-action',
                                          'fit', 'masking', 'small-positive-phase'), 'productive_caught'))])
    for name in sorted({r['predicate'] for r in predicates})])
witnesses = [r for r in trials if r['case'] == 'case-6_figure-2-point-7'
    or (r['case'] == 'two-phase-action' and r['solve'] == '18' and float(r['alpha']) < 1)
    or (r['case'] == 'fit' and r['solve'] == '27' and r['step'] == '1' and r['alpha'] == '0.25')]
table('Trial witnesses', ['Case', 'Solve:step', 'α', 'Maximum coordinate ulps',
      'Predicted merit decrease', 'Assembly rounding estimate', 'Raw after', 'Outcome'], [
    (r['case'], r['solve'] + ':' + r['step'], r['alpha'], r['max_ulps'],
     *[f"{float(r[k]):.6g}" for k in ('predicted_decrease', 'rounding_max', 'raw_after')], r['outcome'])
    for r in witnesses])
table('Unchanged work', ['Case', 'Iterations / evaluations / Hessians', 'Raw maximum'], [
    (r['case'], '/'.join(r[k] for k in ('iterations', 'evaluations', 'hessian_calls')), r['raw_max'])
    for r in rows('work.csv')])
table('Existing comparison checks', ['Case', 'Computed', 'Reference', 'Absolute limit', 'Relative limit'], [
    (r['case'], r['model'], r['reference'], r['atol'], r['rtol']) for r in rows('checks.csv')])
(destination.parent / 'polish-damped.md').write_text('\n'.join(lines))

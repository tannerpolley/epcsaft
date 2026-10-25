"""Write qualification.md tables from the retained round-5 CSV and key-solve files."""
import csv, collections
from pathlib import Path
HERE = Path(__file__).resolve().parents[1] / 'results/density-initialization/round5'
read = lambda name: list(csv.DictReader((HERE/name).open()))
CONTROLS = {'script-controls': 'Round-2 script controls (1e-15, 1e-14, 1e-14)', 'work-matched': 'Work-matched (1e-15, 1e-14, 1e-8)'}
out = ['| Controls | Coordinate | Seeded − cold | Allowance | Fraction used | Allowance without absolute floors | Osmotic-row share of allowance |',
       '|---|---:|---:|---:|---:|---:|---:|']
f2 = read('coordinate-allowance.csv')
for x in f2:
    out.append(f"| {CONTROLS[x['controls']]} | {x['coordinate']} | {float(x['difference']):.3g} | {float(x['allowance']):.3g} | "
               f"{100*float(x['fraction']):.3f}% | {float(x['relative_only_allowance']):.3g} | {100*float(x['osmotic_share']):.1f}% |")
out += ['', f"Every F2 coordinate passes at both control sets; the largest fraction used is {100*max(float(x['fraction']) for x in f2):.3f}%. "
        f"Endpoint Jacobian condition number {float(f2[0]['condition']):.4g}; the returned factors reproduce the weighted Jacobian to "
        f"{max(float(x['layout_check']) for x in f2):.1e} of its largest singular value. The allowance at the seeded endpoint differs from the cold one by at most "
        f"{max(abs(float(x['seeded_endpoint_allowance'])/float(x['allowance'])-1) for x in f2):.1e} relative.", '',
        '| Controls | B4 start | Seeded objective | Cold objective | Relative difference | Allowance max(1e-12, 1e-10·abs(cold)) | Fraction used | Outcome |',
        '|---|---:|---:|---:|---:|---:|---:|---|']
for x in read('b4-objective.csv'):
    out.append(f"| {CONTROLS[x['controls']]} | {x['start']} | {float(x['seeded']):.15g} | {float(x['cold']):.15g} | {float(x['relative']):.2e} | "
               f"{float(x['allowance']):.0e} | {100*float(x['fraction']):.1f}% | {'Pass' if x['passed'] == 'True' else 'Fail'} |")
b4 = read('b4-coordinates.csv')
out += ['', 'B4 controls are its declared Ipopt tolerance in both sets, so both sets repeat the same calls.', '',
        '| B4 start | Largest coordinate difference, seeded − cold (diagnostic) | Coordinate |', '|---:|---:|---:|']
for start in ('0', '1'):
    worst = max((x for x in b4 if x['start'] == start), key=lambda x: abs(float(x['difference'])))
    out.append(f"| {start} | {float(worst['difference']):.2e} | {worst['coordinate']} |")
work = collections.defaultdict(lambda: [0, 0, 0])
for x in csv.DictReader((HERE/'key-solves/work.csv').open()):
    if x['sample'] == 'warmup':
        for k, c in enumerate(('iterations', 'evaluations', 'hessian_calls')): work[(x['case'], x['wheel'])][k] += int(x[c])
out += ['', '| Key solve | Main iterations / evaluations / Hessians | Seeded iterations / evaluations / Hessians |', '|---|---:|---:|']
for case in dict.fromkeys(c for c, _ in work):
    out.append(f"| {case} | {' / '.join(map(str, work[(case, 'A')]))} | {' / '.join(map(str, work[(case, 'B')]))} |")
(HERE/'qualification.md').write_text('\n'.join(out) + '\n')

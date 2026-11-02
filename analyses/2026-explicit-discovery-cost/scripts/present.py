"""Format retained comparisons and categorical work counts; never run the Engine."""
import csv
import gzip
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
R,F = HERE/'results',HERE/'figures'
F.mkdir(exist_ok=True)
def table(name, headings, rows):
    text = '| '+' | '.join(headings)+' |\n|'+'|'.join('---' for _ in headings)+'|\n'
    (F/name).write_text(text+''.join('| '+' | '.join(map(str,r))+' |\n' for r in rows))

equivalence = json.loads((R/'equivalence.json').read_text())
assert not equivalence['failures']
gradients = list(csv.DictReader((R/'terminal-gradients.csv').open()))
table('gradients.md',['Supported gradient comparisons','Maximum error / allowance'],
      [[len(gradients),max(float(r['absolute_error'])/float(r['allowance']) for r in gradients)]])
table('equivalence.md',['Inventory','Entries','Scalar comparisons','Maximum error / allowance','Changed outcomes'],
      [[{'811':'Solver parent/trial pairs'}.get(r[k],r[k]) for k in ('family','entries','comparisons','max_error_over_allowance','changed_outcomes')] for r in equivalence['summary']])
work = list(csv.DictReader((R/'work.csv').open()))
table('work.md',['Case','Trial points A/B','Trial Hessians A/B','Density roots A/B','Helmholtz evaluations A/B'],
      [[r['case']]+[r['A_'+k]+' / '+r['B_'+k] for k in ('points','hessians','roots','helmholtz')] for r in work])
stages = list(csv.DictReader((R/'stages.csv').open()))
table('stages.md',['Case / arm','Compile / µs','Declarations / µs','Seed roots / µs','Trials / µs','Other / µs','Load range'],
      [[r['case']+' / '+r['arm']]+[f'{float(r[k]):.6g}' for k in ('compile','declarations','seed_roots','trials','other')]+[r['load_range']] for r in stages])
neutral = [r for r in work if r['B_hessians']!='shared driver']
fig,ax = plt.subplots(figsize=(9,4.8),layout='constrained')
y = list(range(len(neutral)))
ax.barh([i-.18 for i in y],[float(r['A_hessians']) for r in neutral],height=.36,label='A: eager requests',color='#0072B2')
ax.barh([i+.18 for i in y],[float(r['B_hessians']) for r in neutral],height=.36,label='B: second-order attempts',color='#D55E00',hatch='//')
ax.set(yticks=y,yticklabels=[r['case'] for r in neutral],xlabel='Trial isothermal phase-Hessian requests per discovery / count')
ax.invert_yaxis(); ax.legend(); fig.savefig(F/'hessian-requests.png',dpi=170); plt.close(fig)
if (R/'ratios.csv').exists():
    ratios = list(csv.DictReader((R/'ratios.csv').open()))
    table('cost.md',['Boundary','Ratio','Median','Minimum','Maximum','Target met'],[[r['operation'],r['comparison']]+[f'{float(r[k]):.6g}' for k in ('median','minimum','maximum')]+[r['meets_target']] for r in ratios])
    statistics = list(csv.DictReader((R/'statistics.csv').open()))
    table('times.md',['Boundary','Arm','Median wall / µs','Minimum / µs','Maximum / µs','Median CPU / µs'],[[r['operation'],r['arm']]+[f'{float(r[k]):.3f}' if r[k] else '' for k in ('wall_median_us','wall_min_us','wall_max_us','cpu_median_us')] for r in statistics])
else:
    (F/'cost.md').write_text('Quiet-host cost is deferred. No saving or comparator ratio has been measured under the required admission.\n')
    (F/'times.md').write_text('Ordinary and prepared timing remain unjudged. Deferred host observations are retained in `results/deferred-host.jsonl`.\n')

# Explicit reduced discovery is a separate search population from the historical #855 comparison.
fast = R/'fast-discovery'
if fast.exists():
    verification = json.loads((fast/'verification.json').read_text())
    table('fast-verification.md', ['Population', 'Entries', 'Numerical comparisons', 'Differences'],
          [[name, verification[name]['entries'], verification[name]['comparisons'], verification[name]['differences']]
           for name in ('full', 'fast', 'native', 'pairs')])
    work = list(csv.DictReader((fast/'fast-work.csv').open()))
    with gzip.open(HERE/'inputs/fast-discovery.json.gz', 'rt') as stream:
        inputs = {r['key']: r['input']['problem'] for r in json.load(stream)['cases']}
    table('fast-inputs.md', ['Case', 'T / K', 'P / Pa', 'Feed basis / values in declared order'],
          [[r['case'], inputs[r['case']]['T'], inputs[r['case']]['P'],
            inputs[r['case']]['feed']['basis'].split('.')[-1]+' / ('+
            ', '.join(f'{x:.8g}' for x in inputs[r['case']]['feed']['values'].values())+')'] for r in work])
    table('fast-work.md', ['Case', 'Full trials / H / R', 'Fast trials / H / R', 'Fast limits: trials / H / R'],
          [[r['case'], ' / '.join(str(int(float(r[k]))) for k in ('full_trials', 'full_helmholtz_evaluations', 'full_density_solves')),
            ' / '.join(r[k] for k in ('trials', 'helmholtz_evaluations', 'density_solves')),
            ' / '.join(r[k] or 'report only' for k in ('trial_limit', 'helmholtz_limit', 'density_limit'))] for r in work])
    if (fast/'ratios.csv').exists():
        ratios = list(csv.DictReader((fast/'ratios.csv').open()))
        table('fast-ratios.md', ['Case / boundary', 'Comparison', 'Median', 'Observed range', 'Owner ceiling', 'Control rule'],
              [[r['case']+' / '+r['boundary'], r['comparison'], f"{float(r['median']):.4g}",
                f"{float(r['minimum']):.4g}–{float(r['maximum']):.4g}", r['limit'] or 'report only', r['control_pass'] or '—']
               for r in ratios])
        cases = list(dict.fromkeys(r['case'] for r in work))
        labels = ['Methane/ethane', 'Ionic Figure 2 point 7', 'Water/ethanol', 'water_np/hexane',
                  'Water/isobutanol', 'Dilute isobutanol', 'Aqueous salt']
        fig, ax = plt.subplots(figsize=(8.5, 5), layout='constrained')
        for boundary, offset, color, marker in [('ordinary', -.12, '#0072B2', 'o'), ('prepared', .12, '#D55E00', 's')]:
            rows = [next(r for r in ratios if r['case']==case and r['boundary']==boundary and r['comparison']=='fast/Engine-A') for case in cases]
            medians = [float(r['median']) for r in rows]
            ax.errorbar(medians, [i+offset for i in range(len(cases))],
                        xerr=[[m-float(r['minimum']) for m,r in zip(medians,rows)],
                              [float(r['maximum'])-m for m,r in zip(medians,rows)]],
                        fmt=marker, color=color, label=boundary, capsize=3)
            if boundary=='ordinary':
                at = [(i,float(r['limit'])) for i,r in enumerate(rows) if r['limit']]
                ax.scatter([v for _,v in at], [i for i,_ in at], marker='x', color='black', label='owner ceiling')
        ax.set(yticks=range(len(cases)), yticklabels=labels, xlabel='Fast / full wall time, dimensionless', xlim=(0,1.05))
        ax.axvline(1, color='0.6', linestyle=':', linewidth=1); ax.invert_yaxis(); ax.legend()
        fig.savefig(F/'fast-discovery-cost.png', dpi=170); plt.close(fig)

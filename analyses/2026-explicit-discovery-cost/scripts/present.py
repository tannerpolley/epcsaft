"""Format retained comparisons and categorical work counts; never run the Engine."""
import csv
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

"""Present retained native values only; no Engine import or calculation."""
import csv
import json
from pathlib import Path
import matplotlib.pyplot as plt

LEAF = Path(__file__).resolve().parents[1]
rows = list(csv.DictReader((LEAF/'results/data/comparison.csv').open()))
run = json.loads((LEAF/'results/data/run.json').read_text())
reported = {r['quantity']: float(r['value']) for r in csv.DictReader(
    (LEAF/'inputs/source/schick-2023-reported-metrics.csv').open())}
diag = json.loads((LEAF/'results/data/native-diagnostics.json').read_text())
out = LEAF/'results/rendered'
text = []


def table(title, headers, values):
    text.extend(['\n### '+title+'\n', '| '+' | '.join(headers)+' |', '| '+' | '.join(['---']*len(headers))+' |'])
    text.extend('| '+' | '.join(str(x) for x in row)+' |' for row in values)


table('Comparison outcomes', ['Quantity', 'Evidence', 'Outcome'], [
    ['Source AARD / %', 'numerical verification', f"{run['source_aard_percent']:.12f}; printed target {reported['paper_aard']:.2f}"],
    ['Printed row ARD interval overlaps', 'numerical verification', f"{sum(r['source_ard_interval_overlap']=='True' for r in rows)}/{len(rows)}"],
    ['Native regular closure', 'numerical verification', f"{sum(d['regular']['closed'] for d in diag)}/{len(rows)}"],
    ['Numerical refinement qualified', 'numerical verification', f"{sum(r['numerically_qualified']=='True' for r in rows)}/{len(rows)}"],
    ['Paper liquid x reproduced', 'numerical verification', f"{sum(r['paper_reproduced']=='True' for r in rows)}/{len(rows)}"],
    ['Experimental liquid x reproduced', 'physical comparison', f"{sum(r['measurement_reproduced']=='True' for r in rows)}/{len(rows)}"],
    ['Engine AARD / %', 'physical comparison', f"{run['engine_aard_percent']:.12f}; all {len(rows)} regular closed values; differs from paper agreement level"],
    ['Native run wall time / s', 'provenance', run['wall_seconds']]])
table('Frozen tolerances', ['Quantity', 'Limit'], run['tolerances'].items())
table('Row 1 initial check and native correction', ['Solve', 'Native row tolerance', 'Closure tolerance', 'Liquid x / mol mol⁻¹', 'Classified residual', 'Pressure relative error', 'Max log fugacity difference', 'Closed'],
    [[name, d['requested_row_tolerance'], d['closure_tolerance'], d['mole_fractions'][0],
      d['convergence']['residual'], d['pressure_relative_error'],
      max(map(abs, d['log_fugacity_difference'])), d['closed']]
     for name in ('initial_regular', 'initial_refinement', 'correction', 'regular', 'refinement')
     for d in [diag[0][name]]])
failed = [r for r in rows if r['paper_reproduced'] != 'True']
table('Signs of paper mismatches', ['Mismatches', 'Positive Engine minus paper', 'Negative Engine minus paper'],
    [[len(failed), sum(float(r['paper_difference']) > 0 for r in failed),
      sum(float(r['paper_difference']) < 0 for r in failed)]])
table('Per-isotherm physical comparison', ['T / K', 'Rows', 'Engine AARD / %', 'Paper AARD / %'],
    [[t, len(s), sum(float(r['ard_percent']) for r in s)/len(s), sum(float(r['source_ard_percent']) for r in s)/len(s)]
     for t in sorted({r['temperature_k'] for r in rows}) for s in [[r for r in rows if r['temperature_k']==t]]])
table('Every source row', ['Row', 'T / K', 'P / bar', 'Experiment x', 'Paper x', 'Engine x', 'Vapor y', 'Refinement change', 'Qualified', 'Paper reproduced', 'Signed measurement error', 'ARD / %', 'Bubble delta P / bar', 'Diagnosis'],
    [[r.get(k, '') for k in ('source_row','temperature_k','pressure_bar','experimental_liquid_co2_mole_fraction',
      'paper_model_liquid_co2_mole_fraction','engine_x','engine_y','refinement_change','numerically_qualified',
      'paper_reproduced','signed_error','ard_percent','bubble_pressure_difference_bar','diagnosis')] for r in rows])
table('Independent closure and rounding diagnostics', ['Row', 'Regular residual', 'Refined residual', 'Refined pressure relative error', 'Refined max log fugacity difference', 'dx for +0.05 K', 'dx for +pressure half-width'],
    [[d['source_row'], d['regular']['convergence']['residual'], d['refinement']['convergence']['residual'],
      d['refinement']['pressure_relative_error'], max(map(abs,d['refinement']['log_fugacity_difference'])),
      d['temperature_rounding_action']['action'], d['pressure_rounding_action']['action']] for d in diag])
(LEAF/'results/tables.md').write_text('\n'.join(text)+'\n')
fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
for ax, t in zip(axes, sorted({r['temperature_k'] for r in rows})):
    selected = [r for r in rows if r['temperature_k']==t]
    for key, label, marker in [('experimental_liquid_co2_mole_fraction','Knez via Schick','o'),
                             ('paper_model_liquid_co2_mole_fraction','Paper model','+'), ('engine_x','Engine binary VLE','x')]:
        available = [r for r in selected if r[key]]
        ax.scatter([float(r['pressure_bar']) for r in available], [float(r[key]) for r in available], label=label, marker=marker, s=28)
    ax.set(title=f'{t} K', xlabel='pressure / bar', ylim=(0, 1)); ax.grid(alpha=.2)
axes[0].set_ylabel('liquid CO₂ mole fraction / mol mol⁻¹')
fig.legend(*axes[-1].get_legend_handles_labels(), loc='outside lower center', ncol=3)
fig.suptitle('Full binary VLE at the published interaction parameter')
fig.tight_layout(rect=(0, .09, 1, .93))
fig.savefig(out/'schick-2023-table-s10-co2-thf.png', dpi=180)
plt.close(fig)

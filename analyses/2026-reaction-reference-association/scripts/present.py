"""Summarize retained numerical comparisons; no model evaluation."""
import csv
import json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
summary = defaultdict(lambda: dict(count=0, max_error=0., failures=0))
for name in ('e3', 'mea'):
    rows = list(csv.DictReader((HERE/f'results/{name}.csv').open()))
    for row in rows:
        if row['pass']=='record':
            continue
        group = summary[row['control']]
        group['count'] += 1
        group['max_error'] = max(group['max_error'],float(row['error']))
        group['failures'] += row['pass']!='1'
        group['limit'] = float(row['limit'])
(HERE/'results/e3-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
with (HERE/'results/e3-table.md').open('w') as stream:
    stream.write('| Comparison | Count | Maximum disagreement | Failed comparisons |\n|---|---:|---:|---:|\n')
    for name, group in summary.items():
        stream.write(f'| {name} | {group["count"]} | {group["max_error"]:.6g} | {group["failures"]} |\n')
e4 = json.loads((HERE/'results/e4.json').read_text())
with (HERE/'results/e4-table.md').open('w') as stream:
    stream.write('| Training / held-out pressure | Current squared error / K² | Nasrifar squared error / K² | Reduction | Required |\n|---|---:|---:|---:|---:|\n')
    for panel, values in e4['MDEA_held_out_scores'].items():
        pressures = '40 / 53.3 kPa' if panel=='6a' else '53.3 / 40 kPa'
        stream.write(f'| {pressures} | {2*values["current_held_cost"]:.6f} | {2*values["candidate_held_cost"]:.6f} | {100*values["improvement"]:.6f}% | ≥20% |\n')
with (HERE/'results/sensitivity-table.md').open('w') as stream:
    stream.write('| Control | Minimum signed value | Maximum signed value |\n|---|---:|---:|\n')
    controls = defaultdict(list)
    for row in csv.DictReader((HERE/'results/e3.csv').open()):
        if row['control'].endswith('_signed') or row['control'].startswith('c_density_derivative_'):
            controls[row['control']].append(float(row['actual']))
    for name, values in controls.items():
        stream.write(f'| {name} | {min(values):.8g} | {max(values):.8g} |\n')
with (HERE/'results/starts-table.md').open('w') as stream:
    stream.write('| Candidate / training panel | Physical starts | Signed omitted-start errors | Final coefficient spread | Held cost spread / K² |\n|---|---|---|---:|---:|\n')
    for control in e4['start_controls']:
        fits = [f for f in e4['fits'] if f['candidate']==control['candidate'] and f['training_panel']==control['panel']]
        stream.write(f'| {control["candidate"]} / {control["panel"]} | {[f["physical_start"] for f in fits]} | {control["omitted_start_defects"]} | {control["final_parameter_spread"]:.8g} | {max(f["held_cost"] for f in fits)-min(f["held_cost"] for f in fits):.8g} |\n')


def e5_present():
    import math
    target = {11, 18, 50, 51}
    names = ('A0', 'A6a', 'A6b', 'B6a', 'B6b')
    raw = {name: list(csv.DictReader((HERE/f'results/e5-reactive-{name}.csv').open())) for name in names}
    errors = {name: {int(r['row']): float(r['actual']) for r in rows if r['control'] == 'signed_error_percent'} for name, rows in raw.items()}
    assert all(len(values) == 67 for values in errors.values())
    def score(values, selected):
        v = [values[i] for i in selected]
        return dict(count=len(v), mean_absolute_percent=sum(map(abs, v))/len(v), signed_mean_percent=sum(v)/len(v), maximum_absolute_percent=max(map(abs, v)), mean_absolute_log=sum(abs(math.log1p(x/100)) for x in v)/len(v))
    inventory = json.loads((HERE/'inputs/e5-inventory.json').read_text())
    report = dict(fitted_parameters=0, baseline_mean=32.82626030221431, target_limit=26.26100824177145, declarations={})
    for name in names:
        values = errors[name]
        other = set(values)-target
        neutral = list(csv.DictReader((HERE/f'results/e5-neutral-{name}.csv').open()))
        report['declarations'][name] = dict(target=score(values,target), other=score(values,other), all=score(values,values), target_errors={i:values[i] for i in sorted(target)}, neutral_predictions_K=[float(r['actual']) for r in neutral], other_rows_worse_than_1pp=[i for i in sorted(other) if abs(values[i])>abs(errors['A0'][i])+1], isotherms=[])
        for group in inventory['isotherms']:
            ids = set(group['source_rows'])
            report['declarations'][name]['isotherms'].append(dict(wt_percent=group['nominal_wt_percent'], T_K=group['nominal_T_K'], all=score(values,ids), other=score(values,ids-target)))
        report['declarations'][name]['equation_checks'] = {control:dict(maximum_absolute_error=max(float(r['error']) for r in raw[name] if r['control']==control), failures=sum(r['pass']!='1' for r in raw[name] if r['control']==control)) for control in ('classified_rows','inventory_balance','phase_sum','bubble_fugacity','equilibrium_affinity')}
    for name, control in (('B6a','A6a'),('B6b','A6b')):
        d = report['declarations'][name]
        d['target_pass'] = d['target']['mean_absolute_percent']<=report['target_limit'] and all(abs(errors[name][i])<=abs(errors['A0'][i])+1 for i in target)
        d['matched_4C_pass'] = d['target']['mean_absolute_percent']<report['declarations'][control]['target']['mean_absolute_percent'] and all(abs(errors[name][i])<=abs(errors[control][i])+1 for i in target)
        d['other_pass'] = not d['other_rows_worse_than_1pp'] and d['other']['mean_absolute_percent']<=report['declarations']['A0']['other']['mean_absolute_percent']+1e-8 and all(g['other']['mean_absolute_percent']<=a['other']['mean_absolute_percent']+1e-8 for g,a in zip(d['isotherms'],report['declarations']['A0']['isotherms']))
    (HERE/'results/e5-summary.json').write_text(json.dumps(report,indent=2)+'\n')
    with (HERE/'results/e5-table.md').open('w') as out:
        out.write('| Row | A0 / % | A6a 4C / % | B6a 3B / % | A6b 4C / % | B6b 3B / % |\n|---|---:|---:|---:|---:|---:|\n')
        for i in sorted(target):
            out.write(f'| {i} | '+ ' | '.join(f'{errors[n][i]:.8f}' for n in ('A0','A6a','B6a','A6b','B6b'))+' |\n')
        out.write('| Mean absolute | '+' | '.join(f'{report["declarations"][n]["target"]["mean_absolute_percent"]:.8f}' for n in ('A0','A6a','B6a','A6b','B6b'))+' |\n')
        out.write('| Other 63 mean absolute | '+' | '.join(f'{report["declarations"][n]["other"]["mean_absolute_percent"]:.8f}' for n in ('A0','A6a','B6a','A6b','B6b'))+' |\n')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 3.5), layout='constrained')
    for n, color, marker in zip(('A0','A6a','B6a','A6b','B6b'), ('#333333','#0072B2','#D55E00','#009E73','#CC79A7'), ('o','s','^','D','v')):
        ax.scatter(n, report['declarations'][n]['target']['mean_absolute_percent'], color=color, marker=marker, s=65)
    ax.axhline(report['target_limit'], color='#555555', linestyle='--', label='E5 ceiling: 26.2610%')
    ax.set(xlabel='Frozen neutral declaration (A: 4C; B: 3B)', ylabel='Four-row mean absolute pressure error / %', ylim=(0,45))
    ax.legend(frameon=False)
    fig.savefig(HERE/'results/e5-pressure-error.svg')
    plt.close(fig)

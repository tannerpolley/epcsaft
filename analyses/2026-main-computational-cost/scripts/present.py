"""Tables and one cost figure from retained measurements; never import the Engine."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / 'figures'
LABELS = {'property-bundle': 'Fixed-density property bundle', 'pressure-state': 'Pressure state / density closure',
          'C1-300-prepared': 'Prepared propane saturation', 'C1-300': 'Propane saturation', 'C2-0.5': 'Bubble point', 'C4-298.15': 'Liquid–liquid split',
          'electrolyte-eight-feeds': 'Electrolyte split / eight feeds', 'explicit-discovery': 'Explicit phase discovery',
          'C2-0.5/actions-exact': 'Bubble / actions on solved state', 'C2-0.5/actions-fd': 'Bubble / central differences',
          'C4-298.15/actions-exact': 'Liquid split / actions on solved state', 'C4-298.15/actions-fd': 'Liquid split / central differences'}


def read(name):
    return list(csv.DictReader((HERE / 'results' / name).open()))


def table(name, columns, rows):
    text = '| ' + ' | '.join(columns) + ' |\n|' + '|'.join('---' for _ in columns) + '|\n'
    for row in rows:
        text += '| ' + ' | '.join(map(str, row)) + ' |\n'
    (OUT / name).write_text(text)


def main():
    OUT.mkdir(exist_ok=True)
    ratios, stats, checks, host = [read(name) for name in ('ratios.csv', 'statistics.csv', 'checks.csv', 'host.csv')]
    rows = []
    for r in ratios:
        op = r['operation'].replace('/actions-on-solved-state', '')
        label = LABELS[op] + (' / marginal actions' if '/actions-on-solved-state' in r['operation'] else '')
        rows.append([label, r['comparator'], f"{float(r['package_over_comparator']):.4g}",
                     f"{float(r['minimum_paired_ratio']):.4g}–{float(r['maximum_paired_ratio']):.4g}"])
    expected = [(op, tool) for op in ('property-bundle', 'pressure-state', 'C1-300', 'C1-300-prepared', 'C2-0.5', 'C4-298.15', 'explicit-discovery') for tool in ('FeOs', 'Clapeyron')]
    expected += [(f'{case}/actions-on-solved-state', 'package central differences') for case in ('C2-0.5', 'C4-298.15')]
    for op, tool in expected:
        if not any((r['operation'], r['comparator']) == (op, tool) for r in ratios):
            cause = '; '.join(sorted({r.get('cause') or r.get('quantity', 'value mismatch') for r in checks if r['operation'] in (op, op.replace('/actions-on-solved-state', '/actions-fd')) and r['passed'] != 'True'}))
            rows.append([LABELS[op.replace('/actions-on-solved-state', '')], tool, 'Unavailable', cause or 'Timed value mismatch or incomplete admitted attempt; see timed-value-checks.csv and blocks.csv'])
    table('ratios.md', ['Operation', 'Comparator', 'Package/comparator time', 'Range of paired block ratios'], rows)
    table('timings.md', ['Operation', 'Worker', 'Median (ms)', 'Block medians, min–max (ms)', 'Samples'], [
        [LABELS[r['operation']], r['arm'], f"{1000*float(r['median_seconds']):.6g}",
         f"{1000*float(r['minimum_block_median_seconds']):.6g}–{1000*float(r['maximum_block_median_seconds']):.6g}", r['samples']] for r in stats])
    grouped = {}
    for r in checks:
        grouped.setdefault((r['operation'], r['comparator']), []).append(r)
    def tolerance_fraction(r):
        allowed = float(r.get('allowed_error') or r.get('absolute_tolerance') or 0)
        return float(r.get('absolute_error') or 0) / allowed if allowed else 0
    table('agreement.md', ['Operation', 'Comparison', 'Checked quantities', 'Largest error/allowed error (dimensionless)', 'All pass'], [
        [LABELS[op], comparator, len(rs), f"{max(tolerance_fraction(r) for r in rs):.4g}", all(r['passed']=='True' for r in rs)]
        for (op, comparator), rs in grouped.items()])
    initial = json.loads((HERE / 'results/preflight.json').read_text())
    inputs = json.loads((HERE / 'inputs/inputs.json').read_text())
    declared = json.loads((HERE / 'inputs/declared-cases.json').read_text())
    t1, t2 = inputs['t1'][0], inputs['t2']['ordinary'][0]
    conditions = [[LABELS['property-bundle'], 'propane', t1['T'], f"density {t1['rho']:g} mol/m³", 'pure'],
                  [LABELS['pressure-state'], 'propane', t2['T'], f"{min(s['P'] for s in t2['states']):.12g}–{max(s['P'] for s in t2['states']):.12g} Pa; liquid branch", 'pure']]
    for req in declared['requests']:
        names = ' + '.join(c['id'] for c in declared['models'][req['case']]['components'])
        conditions.append([LABELS[req['id']], names, req['T'], f"{req['P']:g} Pa" if 'P' in req else 'pressure solved', str(req.get('x', [1.]))])
    electrolyte = initial['operations']['electrolyte-eight-feeds']['Engine']['samples'][0]['value'][0]
    conditions.append([LABELS['electrolyte-eight-feeds'], 'water + ethanol + isobutanol + sodium + chloride', electrolyte['temperature'][0], f"{electrolyte['pressure'][0]:g} Pa", 'eight feeds in inputs/electrolyte-feeds.csv'])
    T, P = inputs['t4']['methane_ethane']['points'][0]
    conditions.append([LABELS['explicit-discovery'], 'methane + ethane', T, f'{P:g} Pa', str(inputs['t4']['methane_ethane']['x'])])
    table('conditions.md', ['Operation', 'Components', 'Temperature (K)', 'Fixed density or pressure specification', 'Mole fractions / feed'], conditions)
    tol = inputs['tolerances']
    table('tolerances.md', ['Check', 'Numerical rule'], [
        ['Property bundle', f"absolute difference ≤ {tol['property_rel']:g} × max(absolute values, 1 in each reported unit)"],
        ['Pressure-state density', f"relative difference ≤ {tol['density_rel']:g}"],
        ['Equilibrium pressure and density', f"relative difference ≤ {tol['equilibrium_rel']:g}"],
        ['Equilibrium mole fractions and phase fractions', f"absolute difference ≤ {tol['composition_abs']:g}"],
        ['Exact actions versus central differences', f"absolute difference ≤ {tol['action_abs']:g} + {tol['action_rel']:g} × max(absolute action values)"],
        ['Repeated timed returns versus preflight', 'absolute difference / max(absolute values, 1 in each reported unit) ≤ 1e-10'],
        ['Electrolyte component balance, charge and mole-fraction sum', 'absolute difference ≤ 1e-10']])
    selected_samples = [r for r in read('samples.csv') if r['selected']=='True']
    timed_checks, all_blocks = read('timed-value-checks.csv'), read('blocks.csv')
    table('selection.md', ['Quantity', 'Count'], [
        ['Preflight numerical checks passed', sum(r['passed']=='True' for r in checks)], ['Timed returned-value checks passed', sum(r['passed']=='True' for r in timed_checks)],
        ['Selected blocks', sum(r['selected']=='True' for r in all_blocks)], ['Selected samples', len(selected_samples)],
        ['Disturbed blocks excluded', sum(r['disposition']=='disturbed; excluded' for r in all_blocks)],
        ['Incomplete or superseded blocks', sum(r['disposition']=='incomplete or superseded attempt' for r in all_blocks)]])
    cpu = json.loads(initial['host']['lscpu'])['lscpu']
    model = next(r['data'] for r in cpu if r['field'] == 'Model name:')
    loads = [float(r['load1']) for r in host if r['disposition'] == 'selected quiet block']
    table('host.md', ['Quantity', 'Retained observation'], [
        ['CPU', model], ['Logical CPU for serial workers', initial['host']['cpu']], ['Niceness', initial['host']['nice']],
        ['Numerical threads per worker', 'one'], ['Observed one-minute load range', f'{min(loads):.3g}–{max(loads):.3g}'],
        ['Accepted measurement snapshots', len(loads)], ['Compiler observations in accepted blocks', sum(r['compile_processes'] != '[]' for r in host if r['disposition']=='selected quiet block')]])
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [1.55, 1]})
    operations = ['property-bundle', 'pressure-state', 'C1-300', 'C1-300-prepared', 'C2-0.5', 'C4-298.15', 'explicit-discovery']
    for tool, offset, color, marker in [('FeOs', -.12, '#0072B2', 'o'), ('Clapeyron', .12, '#D55E00', 's')]:
        for i, op in enumerate(operations):
            r = next((r for r in ratios if (r['operation'], r['comparator']) == (op, tool)), None)
            if r is None:
                continue
            x, lo, hi = [float(r[k]) for k in ('package_over_comparator', 'minimum_paired_ratio', 'maximum_paired_ratio')]
            axes[0].errorbar(x, i+offset, xerr=[[x-lo], [hi-x]], fmt=marker, color=color, capsize=3, label=tool if i==0 else None)
    axes[0].set_yticks(range(len(operations)), [LABELS[o] for o in operations])
    axes[0].invert_yaxis()
    axes[0].legend()
    axes[0].set_title('Matched operations')
    derivatives = [r for r in ratios if r['comparator']=='package central differences']
    for i, r in enumerate(derivatives):
        x, lo, hi = [float(r[k]) for k in ('package_over_comparator', 'minimum_paired_ratio', 'maximum_paired_ratio')]
        axes[1].errorbar(x, i, xerr=[[x-lo], [hi-x]], fmt='o', color='#009E73', capsize=3)
    axes[1].set_yticks(range(len(derivatives)), ['Bubble actions', 'Liquid-split actions'][:len(derivatives)])
    axes[1].invert_yaxis()
    axes[1].set_title('Actions on an already solved state')
    for axis, label in zip(axes, ['Package / comparator time\n(dimensionless)', 'Exact / central-difference time\n(dimensionless)']):
        axis.set_xscale('log')
        axis.axvline(1, color='.4', linestyle='--', linewidth=.8)
        axis.set_xlabel(label)
        axis.grid(axis='x', alpha=.2)
    fig.suptitle('Current-main computational cost; intervals show paired block spread', fontsize=12)
    fig.tight_layout()
    fig.savefig(OUT / 'cost.svg', metadata={'Date': None})
    plt.close(fig)


if __name__ == '__main__':
    main()

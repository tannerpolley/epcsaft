"""Render figures/family.png from results/comparison.csv; never runs the Engine."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
rows = list(csv.DictReader((ROOT / 'results/comparison.csv').open()))
for r in rows:
    r['before'], r['after'] = r['before_success'] == 'True', r['after_success'] == 'True'
STATUS = [('converged before and after', lambda r: r['before'] and r['after'], dict(marker='o', color='#0072B2')),
          ('converged after only', lambda r: r['after'] and not r['before'], dict(marker='^', color='#E69F00')),
          ('fails after', lambda r: not r['after'], dict(marker='x', color='#D55E00'))]

fig = plt.figure(figsize=(10, 7.5), layout='constrained')
grid = fig.add_gridspec(2, 3, height_ratios=(1, 1.3))
temperatures = sorted({float(r['T']) for r in rows})
for column, T in enumerate(temperatures):
    ax = fig.add_subplot(grid[0, column])
    for label, test, style in STATUS:
        pick = [r for r in rows if float(r['T']) == T and test(r)]
        ax.scatter([float(r['loading']) for r in pick], [float(r['P']) / 1e3 for r in pick], s=36, label=label, **style)
    ax.set(xscale='log', title=f'T = {T:.2f} K', xlabel='CO$_2$ loading (mol/mol MEA)')
    if column == 0:
        ax.set_ylabel('P (kPa)')
fig.axes[1].legend(fontsize=8, loc='center left', bbox_to_anchor=(0.06, 0.42))

ax = fig.add_subplot(grid[1, :])
common = [r for r in rows if r['before'] and r['after']]
x, y = [int(r['before_iterations']) for r in common], [int(r['after_iterations']) for r in common]
ax.scatter(x, y, marker='o', color='#0072B2', label='common converged point')
key = rows[0]
ax.scatter([int(key['before_iterations'])], [int(key['after_iterations'])], marker='s', s=90, facecolor='none',
           edgecolor='#D55E00', linewidth=2, label='key flash (318.15 K, 5918.5 Pa)')
limit = (10, 200)
ax.plot(limit, limit, color='0.4', linestyle='--', linewidth=1, label='equal iterations')
ax.set(xscale='log', yscale='log', xlim=limit, ylim=limit, xlabel='iterations before PR #370 (7fa511b6)',
       ylabel='iterations on main 6bbaaf1e')
ax.legend(fontsize=8)
(ROOT / 'figures').mkdir(exist_ok=True)
fig.savefig(ROOT / 'figures/family.png', dpi=150)

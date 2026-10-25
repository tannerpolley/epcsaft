"""Present retained numerical identity errors; no thermodynamic calculation."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt'] = 'trace-solvent-activity'
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
rows = list(csv.DictReader((HERE / 'results/resolved-comparison.csv').open()))
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6), constrained_layout=True)
for salt, color, marker in zip(('LiCl', 'NaCl', 'MgCl2', 'CaCl2'),
                                ('#0072B2', '#D55E00', '#009E73', '#CC79A7'), ('o', 'x', '^', 's')):
    selected = [r for r in rows if r['salt'] == salt]
    x = [float(r['molality_mol_kg']) for r in selected]
    delta = [float(r['osmotic_change']) for r in selected]
    residual = [float(r['osmotic_change']) - float(r['pressure_term']) for r in selected]
    axes[0].scatter(x, delta, s=12, color=color, marker=marker, label=salt, alpha=0.75)
    axes[1].scatter(x, residual, s=12, color=color, marker=marker, alpha=0.75)
for ax, limit in zip(axes, (1e-6, 5e-9)):
    ax.axhline(limit, color='0.25', linestyle='--', linewidth=1)
    ax.axhline(-limit, color='0.25', linestyle='--', linewidth=1)
    ax.set_xscale('log')
    ax.set_xlabel('Formula-unit molality (mol/kg water)')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
axes[0].set_ylabel('Common-pressure − retained Φ\n(dimensionless)')
axes[1].set_ylabel('Pressure-corrected difference in Φ\n(dimensionless)')
axes[0].legend(frameon=False, fontsize=9)
(HERE / 'figures').mkdir(exist_ok=True)
fig.savefig(HERE / 'figures/resolved-identity.svg', metadata={'Date': None})
plt.close(fig)

"""Compare retained native Mg/Ca predictions with evaluated tables; no model calls."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
fig, axes = plt.subplots(2, 2, figsize=(9, 6), layout='constrained', sharex=True)
residual, error_axes = plt.subplots(2, 2, figsize=(9, 6), layout='constrained', sharex=True)
for column, salt in enumerate(('MgCl2', 'CaCl2')):
    rows = list(csv.DictReader((HERE/f'results/stage2-{salt}-stored-rows.csv').open()))
    for row, (quantity, label) in enumerate((('gamma', 'γ±'), ('osmotic', 'Φ'))):
        data = [v for v in rows if v['quantity'] == quantity]
        x = [float(v['molality']) for v in data]
        ax = axes[row, column]
        ax.plot(x, [float(v['prediction']) for v in data], color='#0072B2', label='Native minimax candidate')
        ax.scatter(x, [float(v['observed']) for v in data], marker='o', facecolors='none', edgecolors='#D55E00', label='Evaluated reference')
        ax.set(ylabel=f'{label} (1)', title=salt if row==0 else '')
        ax.set_xscale('log')
        error = error_axes[row, column]
        error.scatter(x, [float(v['error_decade']) for v in data], color='#0072B2', s=15)
        error.axhspan(-.02, .02, color='0.88', zorder=0, label='Owner criterion ±0.02 decade')
        error.axhline(0, color='0.4', linewidth=.6)
        error.set(ylabel=f'log₁₀(model / evaluated) {label}\n(decade)', title=salt if row==0 else '')
        error.set_xscale('log')
        if row==1:
            ax.set_xlabel('Formula-unit molality (mol/kg water)')
            error.set_xlabel('Formula-unit molality (mol/kg water)')
for figure, panel, name in ((fig, axes, 'stage2-predictions.png'), (residual, error_axes, 'stage2-residuals.png')):
    handles, labels = panel[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='outside lower center', ncol=2)
    figure.suptitle('298.15 K, 100000 Pa; evaluated-reference fitting, no physical validation', fontsize=11)
    figure.savefig(HERE/'figures'/name, dpi=160)
    plt.close(figure)

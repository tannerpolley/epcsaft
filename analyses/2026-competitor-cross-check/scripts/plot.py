"""Presentation only: ordered input contrasts and strict Helmholtz differences."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
def rows(name):
    return list(csv.DictReader((HERE / 'results' / name).open()))
values = rows('values.csv')
decomposition = rows('decomposition.csv')
adjusted = rows('source-adjusted.csv')
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout='constrained')
gamma = {(r['id'], r['implementation']): float(r['value']) for r in values if r['quantity'] == 'gamma'}
axes[0].plot(range(6), [gamma[r['id'], r['implementation']] for r in decomposition], 'o', color='#0072B2')
axes[0].set(xticks=range(6), xticklabels=['Held', '0.88σ', 'constants', 'kCl=0', 'finite ref.', 'zero ref.'], ylabel='Molality-scale γ± [1]', title='Ordered NaCl input/reference contrasts')
labels = [(case, lib) for case in ('C1','C2','C3','C4','C5') for lib in ('FeOs','Clapeyron')]
for marker, color, table, field, label in [('o','#D55E00',values,'difference','raw'), ('x','#0072B2',adjusted,'remaining_signed_difference','source-adjusted')]:
    errors = [max(abs(float(r[field])) for r in table if r['case']==case and r['implementation']==lib and r['quantity']=='ar' and r.get(field) and (r.get('target') or table is adjusted)) for case,lib in labels]
    axes[1].plot(range(len(labels)), errors, marker, color=color, label=label)
axes[1].set(yscale='log', xticks=range(len(labels)), xticklabels=[f'{c}\n{l}' for c,l in labels], ylabel='Maximum |Δaʳ| [1]', title='Common-state energy: raw versus diagnosed source error')
axes[1].legend()
for ax in axes:
    ax.tick_params(axis='x', labelsize=8)
    ax.grid(axis='y', alpha=.2)
(HERE/'figures').mkdir(exist_ok=True)
fig.savefig(HERE/'figures/comparison.svg')

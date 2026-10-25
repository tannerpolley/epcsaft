"""Plot retained F3 action errors; no thermodynamic calculation runs here."""
import csv
import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt

study = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--zero-inventory', action='store_true')
parser.add_argument('--lithium-trace', action='store_true')
args = parser.parse_args()
if args.lithium_trace:
    rows = list(csv.DictReader((study / 'results/lithium-trace-endpoints.csv').open()))
    fig, ax = plt.subplots(figsize=(8, 4.5), layout='constrained')
    for label, color, marker, name in [('historical', '#D55E00', 'o', 'Original wheel ddf63226 (reconstructed)'),
                                        ('current', '#0072B2', 'x', 'Corrected wheel e28d258d')]:
        selected = [row for row in rows if row['source_label'] == label and row['path'] == 'frozen_asl']
        if len(selected) != 10:
            raise SystemExit(f'Expected ten retained frozen-path endpoints for {label}')
        style = dict(facecolors='none', edgecolors=color) if marker == 'o' else dict(color=color)
        ax.scatter([float(row['inlet_li_mol_l']) for row in selected],
                   [float(row['organic_lithium_carrier_fraction']) for row in selected],
                   marker=marker, s=60, linewidths=1.5, label=name, **style)
    ax.set(xscale='log', yscale='log', xlabel='Declared aqueous inlet Li concentration (mol/L)',
           ylabel='Organic lithium-carrier mole fraction (1)',
           title='Historical trace-Li endpoint fractions are preserved')
    ax.legend(frameon=False)
    ax.grid(alpha=.2)
    fig.savefig(study / 'figures/lithium-trace-endpoints.png', dpi=180)
    raise SystemExit(0)
if args.zero_inventory:
    rows = list(csv.DictReader((study / 'results/zero-inventory.csv').open()))
    coordinate_names = ['k_MEACOO_water', 'k_MEAH_water', 'k_HCO3_water', 'k_MEACOO_MEAH',
                        'inverse_T_k_MEAH_water', 'R2_b', 'R2_c']
    fig, ax = plt.subplots(figsize=(8, 4.5), layout='constrained')
    for alpha, color, marker in [('0.110000', '#0072B2', 'o'), ('0.450000', '#D55E00', 's')]:
        selected = [row for row in rows if row['comparison'].startswith(f'zero-inventory/{alpha}/column/')]
        if len(selected) != 7:
            raise SystemExit(f'Expected seven retained native fitted-coordinate errors for alpha={alpha}')
        ax.scatter([coordinate_names.index(row['comparison'].split('/')[-1]) for row in selected],
                   [json.loads(row['inputs'])['normalized_error'] for row in selected],
                   label=f'α={alpha}', color=color, marker=marker, s=38)
    ax.axhline(1, color='black', linestyle='--', linewidth=1, label='F3 allowed error')
    ax.set(yscale='log', ylabel='Native action error / allowed error (1)', xlabel='Fitted coordinate',
           xticks=range(7), xticklabels=['k₄₂', 'k₃₂', 'k₅₂', 'k₄₃', '1/T slope', 'R2 b', 'R2 c'],
           title='Fresh-solvent finite-addition derivatives agree with re-solves')
    ax.legend(frameon=False)
    ax.grid(axis='y', alpha=.2)
    fig.savefig(study / 'figures/zero-inventory-action-errors.png', dpi=180)
    raise SystemExit(0)
rows = list(csv.DictReader((study / 'results/comparisons.csv').open()))
coordinates = ['kij', 'epsilon_tied', 'reaction_a', 'reaction_b', 'reaction_c', 'reaction_d']
variants = ['fixed_system', 'fixed_source', 'free_pressure', 'two_phase', 'unequal_temperature']
colors = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#8B6500']
fig, ax = plt.subplots(figsize=(8, 4.5), layout='constrained')
for offset, (variant, color, marker) in enumerate(zip(variants, colors, ['o', 'x', 's', '^', 'D'])):
    selected = {r['comparison'].split('/')[2]: r for r in rows
                if r['comparison'].startswith(f'addition-f3/{variant}/') and r['comparison'].endswith('/fine')}
    error = [float(selected[c]['absolute_error']) / float(selected[c]['atol']) for c in coordinates]
    ax.scatter([i + (offset - 2) * .1 for i in range(6)], error, color=color, marker=marker,
               label=variant.replace('_', ' '), s=38)
ax.axhline(1, color='black', linestyle='--', linewidth=1, label='F3 allowed error')
ax.set(yscale='log', ylabel='Absolute action error / allowed error (1)', xlabel='Fitted coordinate',
       xticks=range(6), xticklabels=['kᵢⱼ', 'tied ε/k', 'a', 'b', 'c', 'd'],
       title='Finite-addition derivatives agree with centered re-solves')
ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False)
ax.grid(axis='y', alpha=.2)
(study / 'figures').mkdir(exist_ok=True)
fig.savefig(study / 'figures/enthalpy-action-errors.png', dpi=180)

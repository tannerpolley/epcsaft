"""Show native derivative refinement from retained CSV; no model execution."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'results/enthalpy-shift/columns.csv'
OUTPUT = ROOT / 'figures/enthalpy-shift-refinement.png'
rows = list(csv.DictReader(SOURCE.open()))
colors = dict(zip((293.15, 313.15, 353.15, 393.15), ('#0072B2', '#D55E00', '#009E73', '#CC79A7')))
fig, axes = plt.subplots(1, 3, figsize=(11.5, 4), sharey=True)
for ax, kind in zip(axes, ('pressure', 'species', 'kim-addition')):
    for temperature, color in colors.items():
        for reaction, marker in (('R2', 'o'), ('R5', 'x')):
            selected = [r for r in rows if r['kind'] == kind and r['reaction'] == reaction
                        and float(r['temperature']) == temperature and float(r['h0']) == 0.0]
            fractions = [float(r['absolute_error']) / (1e-8 if abs(float(r['exact_weighted_residual'])) <= 1e-8 else
                         1e-5 * max(abs(float(r['exact_weighted_residual'])), abs(float(r['fd_weighted_residual'])))) for r in selected]
            ax.scatter([float(r['step']) for r in selected], fractions, color=color, marker=marker, s=27)
    ax.axhline(1, color='black', lw=.8, ls='--')
    ax.set(xlabel='Centered physical step (J/mol)', title=kind, xscale='log', yscale='symlog', ylim=(-1e-8, 2))
    ax.set_yscale('symlog', linthresh=1e-7)
    ax.set_xticks((3.125, 6.25, 12.5, 25), ('3.125', '6.25', '12.5', '25'))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.grid(axis='y', alpha=.18)
axes[0].set_ylabel('Residual-column error / MEA criterion (1)')
handles = [Line2D([], [], color=c, marker='s', ls='', label=f'{t:g} K') for t, c in colors.items()]
handles += [Line2D([], [], color='black', marker=m, ls='', label=r) for r, m in (('R2', 'o'), ('R5', 'x'))]
fig.legend(handles=handles, loc='upper center', ncol=6, frameon=False)
fig.tight_layout(rect=(0, 0, 1, .9))
fig.savefig(OUTPUT, dpi=180)
plt.close(fig)
(ROOT / 'results/enthalpy-shift/figure-identity.json').write_text(json.dumps(
    {'source': str(SOURCE.relative_to(ROOT)), 'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
     'figure': str(OUTPUT.relative_to(ROOT)), 'basis': 'zero physical shifts; native weighted residual columns; MEA 1e-5 relative and 1e-8 near-zero criterion'}, indent=2) + '\n')

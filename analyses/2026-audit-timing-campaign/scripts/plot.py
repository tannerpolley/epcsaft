"""Figure from retained statistics only: Engine time over comparator time for value-matched equal work (no Engine import)."""
import csv
from pathlib import Path

import matplotlib

matplotlib.use('svg')
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
rows = [r for r in csv.DictReader((HERE / 'results/statistics.csv').open())
        if r['tool'] in ('FeOs', 'Clapeyron') and r.get('value_matched') == 'True' and r['section'] in ('T1', 'T2', 'T3', 'T4', 'T6')]
labels = {'T1': 'T1 fixed-density properties', 'T2': 'T2 density closure', 'T3': 'T3 saturation and bubble', 'T4': 'T4 phase discovery',
          'T6': 'T6 activity curve'}
fig, ax = plt.subplots(figsize=(8, 4.8))
sections = [s for s in labels if any(r['section'] == s for r in rows)]
for tool, marker, colour in (('FeOs', 'o', '#3b6ea8'), ('Clapeyron', 's', '#c05a2b')):
    for i, s in enumerate(sections):
        pts = [r for r in rows if r['section'] == s and r['tool'] == tool]
        y = [1 / float(r['ratio_to_engine']) for r in pts]  # Engine time / comparator time
        x = [i + (-0.12 if tool == 'FeOs' else 0.12) + 0.06 * ((k % 5) - 2) / 2 for k in range(len(y))]
        ax.scatter(x, y, marker=marker, s=18, color=colour, alpha=0.75, label=tool if i == 0 else None, linewidths=0)
ax.axhline(1, color='0.3', lw=0.8)
ax.set_yscale('log')
ax.set_xticks(range(len(sections)), [labels[s] for s in sections], rotation=15, ha='right')
ax.set_ylabel('Engine time / comparator time (median, same work)')
ax.set_title('Above 1: Engine slower; each point one value-matched operation')
ax.legend(frameon=False)
fig.tight_layout()
(HERE / 'figures').mkdir(exist_ok=True)
fig.savefig(HERE / 'figures/engine-over-comparator.svg')

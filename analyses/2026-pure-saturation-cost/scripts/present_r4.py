"""Render the retained 30-pair prepared-call comparison without solving a state."""
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / 'results/dedicated-integrated-r4'
FIGURES = HERE / 'figures'
r = json.loads((DATA / 'saving.json').read_text())[0]
c = json.loads((DATA / 'campaign.json').read_text())
pairs = [x for x in csv.DictReader((DATA / 'paired-ratios.csv').open()) if x['comparison'] == 'B/A']
assert len(pairs) == c['pairs_completed'] == 30
rows = [(label, f'{r[key]:.12g}') for label, key in (
    ('Paired median ratio R / 1', 'R'), ('−ln R / 1', 'negative_log_R'),
    ('Same-image envelope C / 1', 'C'), ('Median paired saving / µs', 'saving_median_us'))]
rows += [('Paired saving range / µs', f"{r['saving_range_us'][0]:.12g} to {r['saving_range_us'][1]:.12g}"),
         ('Accepted pairs', c['pairs_completed']), ('Passing quantity comparisons', c['quantity_checks']),
         ('Rejected blocks', c['rejected_blocks']), ('Accepted blocks', c['accepted_blocks']),
         ('Deferred admission snapshots', c['deferred_snapshots'])]
rows += [(f'Accepted {minutes}-minute load range / 1', f'{lo:.2f}–{hi:.2f}')
         for minutes, (lo, hi) in c['accepted_load_ranges'].items()]
rows += [('First accepted timing / UTC', c['first_timing_utc']), ('Last accepted timing / UTC', c['last_timing_utc'])]
(FIGURES / 'dedicated-integrated-r4.md').write_text(
    '| Prepared-call quantity | r4 observation |\n|---|---:|\n' +
    ''.join(f'| {label} | {value} |\n' for label, value in rows))
fig, ax = plt.subplots(figsize=(7.2, 3.6), layout='constrained')
ax.plot(range(1, len(pairs) + 1), [float(x['ratio']) for x in pairs], 'o',
        color='#0072B2', label='Paired B/A block medians')
ax.axhline(math.exp(-r['C']), color='#D55E00', linestyle='--', label='exp(−C): strict saving threshold')
ax.axhline(r['R'], color='#555555', linestyle=':', label='Median R')
ax.set_xlabel('Paired observation (1)')
ax.set_ylabel('Candidate / baseline wall time (1)')
ax.legend(fontsize=8)
ax.grid(axis='y', alpha=.2)
fig.savefig(FIGURES / 'dedicated-integrated-r4-saving.png', dpi=180)
plt.close(fig)

"""Plot retained block medians and ranges, with historical competitor references."""
import csv
from datetime import datetime
from datetime import timedelta
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from summarize import LABELS

HERE = Path(__file__).resolve().parents[1]
GROUPS = {
    'ionic': ('Ionic liquid–liquid equilibrium', list(LABELS)[:2] + ['khudaida/sequential-8']),
    'mea': ('Reacting MEA equilibrium', list(LABELS)[2:4] + ['mea-nominally-unloaded/353.15-K']),
    'vle': ('Saturation and bubble points', ['C1-300', 'C2-0.5', 'methane-bubble']),
    'lle-flash': ('Neutral liquid splits and TP flash', ['C4-298.15', 'C4-323.15', 'methane-flash']),
    'density': ('Liquid density at changing pressures', ['A-propane-liquid']),
    'f1': ('Frozen MEA F1 residual and Jacobian', ['mea-f1/residual-jacobian']),
    'discovery': ('Methane/ethane phase discovery', [op for op in LABELS if op.startswith('methane-discover/')]),
    'cold-boundaries': ('Cold saturation and high-temperature boundaries', [op for op in LABELS if op.startswith(('propane-saturation/', 'methane-bubble/', 'methane-dew/'))]),
}


def main():
    r = HERE / 'results'
    campaign = json.loads((r / 'campaign.json').read_text())
    versions = [v for v in campaign['environments'] if not v['version'].endswith('-ma57')]
    ma57 = campaign['current'] + '-ma57'
    dates = {v['version']: datetime.fromisoformat(v['date']) for v in versions}
    stats = {(s['version'], s['operation']): s for s in csv.DictReader((r / 'statistics.csv').open())}
    headline = {s['operation']: s for s in csv.DictReader((r / 'headline.csv').open())}
    refs = {(s['operation'], s['library']): float(s['seconds']) for s in csv.DictReader((r / 'competitors.csv').open())}
    audit_refs = {(s['operation'], s['library']): float(s['seconds']) for s in csv.DictReader((r / 'competitors-578.csv').open())}
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False,
                         'svg.fonttype': 'none', 'savefig.bbox': 'tight'})
    for name, (title, operations) in GROUPS.items():
        scale = 1 if name == 'f1' else 1e3
        unit = 's' if name == 'f1' else 'ms'
        height = 2.9*len(operations)+1.2
        fig, axes = plt.subplots(len(operations), 2, figsize=(11.6, height),
                                 squeeze=False, sharey='row', width_ratios=[1.1, 1])
        for row_axes, op in zip(axes, operations):
            rows = [(dates[v['version']], stats[(v['version'], op)]) for v in versions if (v['version'], op) in stats]
            x = [d for d, s in rows]
            y = [float(s['median_seconds'])*scale for d, s in rows]
            lo = [m-float(s['minimum_seconds'])*scale for m, (d, s) in zip(y, rows)]
            hi = [float(s['maximum_seconds'])*scale-m for m, (d, s) in zip(y, rows)]
            # Discrete releases: points and ranges, with no interpolated performance curve.
            h = headline[op]
            for ax in row_axes:
                if (ma57, op) in stats:
                    m = stats[(ma57, op)]
                    y_ma57 = float(m['median_seconds'])*scale
                    ax.errorbar([dates[campaign['current']]], [y_ma57],
                        yerr=[[y_ma57-float(m['minimum_seconds'])*scale], [float(m['maximum_seconds'])*scale-y_ma57]],
                        fmt='s', color='#CC79A7', capsize=5, markersize=5,
                        label='Current MA57 (HSL): median and range', zorder=5)
                ax.errorbar(x, y, yerr=[lo, hi], fmt='o', color='#0072B2', capsize=4,
                            markersize=5, label='Engine: median and range of 3 blocks', zorder=4)
                ax.axhspan(float(h['aa_min_seconds'])*scale, float(h['aa_max_seconds'])*scale,
                           color='#999999', alpha=.23, label="Current MUMPS A=A range (9 blocks)")
                for library, color, dash in [('Clapeyron', '#D55E00', '--'), ('FeOs', '#009E73', ':')]:
                    reference = refs.get((op, library), audit_refs.get((op, library)))
                    if reference is not None:
                        source = '#533' if (op, library) in refs else '#578'
                        ax.axhline(reference*scale, color=color, linestyle=dash,
                                   linewidth=1.7, label=f'{library}: retained {source} median')
                ax.set_yscale('log')
                ax.yaxis.set_major_locator(ticker.LogLocator(base=10, subs=(1, 2, 5)))
                ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda value, position: f'{value:g}'))
                ax.yaxis.set_minor_formatter(ticker.NullFormatter())
                ax.grid(axis='y', which='major', alpha=.22)
            left, right = row_axes
            left.set_ylabel(f'Wall time / {unit} (log)')
            left.set_title(LABELS[op] if len(operations) > 1 else 'October 6–19 overview', loc='left', fontsize=12)
            left.set_xlim(min(dates.values())-timedelta(hours=7), max(dates.values())+timedelta(hours=7))
            left.xaxis.set_major_locator(mdates.DayLocator(interval=2))
            left.xaxis.set_major_formatter(mdates.DateFormatter('%b %d', tz=x[0].tzinfo))
            recent = [(v, d) for v, d in dates.items() if d.day >= 15]
            right.set_xlim(min(d for v, d in recent)-timedelta(hours=3), max(d for v, d in recent)+timedelta(hours=4))
            right.xaxis.set_major_locator(mdates.DayLocator())
            right.xaxis.set_major_formatter(mdates.DateFormatter('%b %d', tz=x[0].tzinfo))
            right.set_title('October 15–19 detail', loc='left', fontsize=10)
        axes[-1, 0].set_xlabel('Version date, 2026')
        axes[-1, 1].set_xlabel('Version date (America/Denver), 2026')
        fig.suptitle(title, x=.09, y=1-.1/height, ha='left', fontsize=16)
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(.07, .28/height), ncol=2,
                   fontsize=8, frameon=False)
        fig.text(.07, .08/height, 'HSL, a collection of Fortran codes for large-scale scientific computation. https://www.hsl.rl.ac.uk/', fontsize=8)
        fig.subplots_adjust(left=.09, right=.98, top=1-.75/height, bottom=1.55/height, hspace=.85, wspace=.12)
        fig.savefig(HERE / 'figures' / f'{name}.svg', metadata={'Date': None})
        inspection = HERE.parents[1] / 'build/issue-549/figure-inspection'
        inspection.mkdir(parents=True, exist_ok=True)
        fig.savefig(inspection / f'{name}.png', dpi=140, metadata={'Software': 'matplotlib'})
        plt.close(fig)
    print(f'Wrote {len(GROUPS)} operation-group charts from retained statistics, headlines and competitor medians')


if __name__ == '__main__':
    main()

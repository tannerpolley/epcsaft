"""Plot retained physical comparisons only; Haslam predictions are not recomputed."""
import csv
import json
from pathlib import Path
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
rounding = json.loads((ROOT / 'results/data/haslam-2020-table8.json').read_text())['printed_rounding_absolute_tolerance']
with (ROOT / 'results/data/haslam-2020-table8.csv').open() as stream:
    rows = list(csv.DictReader(stream))
if not any(r['status'] == 'evaluated' for r in rows):
    raise SystemExit('No native predictions available; refusing an empty comparison figure')
figure, axes = plt.subplots(3, 3, figsize=(12, 10), sharex=True)
for axis, salt in zip(axes.flat, dict.fromkeys(r['salt'] for r in rows), strict=True):
    for prop, symbol, color in [('gamma_pm_m', 'γ± (molality)', '#D55E00'), ('osmotic_coefficient', 'Φ', '#0072B2')]:
        group = [r for r in rows if r['salt'] == salt and r['observable'] == prop]
        valid = [r for r in group if r['status'] == 'evaluated']
        axis.scatter([float(r['molality_mol_kg']) for r in group], [float(r['literature_value']) for r in group],
                     marker='o', facecolors='none', edgecolors=color, label=f'Hamer–Wu {symbol}')
        axis.scatter([float(r['molality_mol_kg']) for r in valid], [float(r['model_value']) for r in valid],
                     marker='x', color=color, label=f'Current ePC-SAFT {symbol}')
    axis.set(title=salt, xscale='log', ylabel='Dimensionless coefficient')
    if not any(r['status'] == 'evaluated' for r in rows if r['salt'] == salt):
        axis.text(0.5, 0.5, 'No ePC-SAFT LiI prediction\nMissing source-backed pair', ha='center', transform=axis.transAxes)
    axis.grid(alpha=0.2)
figure.legend(*axes.flat[0].get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(0.5, 0.94), ncols=2)
figure.suptitle('ePC-SAFT comparison with Hamer–Wu; Haslam Table 8 values are not recomputed')
figure.supxlabel('Formula-unit molality / mol kg⁻¹ water')
figure.tight_layout(rect=(0, 0.03, 1, 0.87))
output = ROOT / 'results/rendered'
output.mkdir(exist_ok=True)
figure.savefig(output / 'haslam-2020-table8.png', dpi=180)
plt.close(figure)
figure, axes = plt.subplots(3, 3, figsize=(12, 9), sharex=True)
for axis, salt in zip(axes.flat, dict.fromkeys(r['salt'] for r in rows), strict=True):
    for prop, label, color in [('gamma_pm_m', 'γ± (molality)', '#D55E00'), ('osmotic_coefficient', 'Φ', '#0072B2')]:
        group = [r for r in rows if r['salt'] == salt and r['observable'] == prop and r['status'] == 'evaluated']
        axis.scatter([float(r['molality_mol_kg']) for r in group], [float(r['residual']) for r in group], color=color, label=label)
    axis.axhspan(-rounding, rounding, color='grey', alpha=0.2)
    axis.axhline(0, color='grey', linewidth=0.5)
    axis.set(title=salt, xscale='log', ylabel='ePC-SAFT − Hamer–Wu')
figure.suptitle('Coefficient residuals; grey band is the fixed printed-rounding interval')
figure.supxlabel('Formula-unit molality / mol kg⁻¹ water')
figure.legend(*axes.flat[0].get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(0.5, 0.95), ncols=2)
figure.tight_layout(rect=(0, 0.03, 1, 0.90))
figure.savefig(output / 'haslam-2020-table8-residuals.png', dpi=180)
plt.close(figure)

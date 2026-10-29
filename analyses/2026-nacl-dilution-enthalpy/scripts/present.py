"""Read retained dilution values to prepare the notebook tables and figure."""
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
identity = json.loads((HERE/'results/run-identity.json').read_text())
for name, expected in identity['input_sha256'].items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == expected, name
rows = list(csv.DictReader((HERE/'results/comparison.csv').open()))
chords = json.loads((HERE/'results/chords.json').read_text())
summary = json.loads((HERE/'results/summary.json').read_text())
states = json.loads((HERE/'results/native-states.json').read_text())
qualification = json.loads((HERE/'results/qualification.json').read_text())
heat = list(csv.DictReader((HERE/'results/heat-capacity.csv').open()))
assert len(rows) == sum(g['rows'] for g in summary['groups'][:-1])


def table(name, header, data):
    text = '| ' + ' | '.join(header) + ' |\n| ' + ' | '.join(['---']*len(header)) + ' |\n'
    text += ''.join('| ' + ' | '.join(map(str, row)) + ' |\n' for row in data)
    (HERE/f'figures/{name}.md').write_text(text)


def number(value, precision=3):
    return f'{float(value):.{precision}f}' if value is not None and value != '' else 'unavailable'


table('conditions', ['Condition', 'Value'], [
    ['Temperatures / K', ', '.join(number(g['T_K'], 2) for g in summary['groups'][:-1])],
    ['Nominal pressure / Pa', identity['P_Pa']], ['Water molar mass / kg mol⁻¹', identity['MW_kg_mol']],
    ['Salt reporting inventory / mol', 1], ['Source rows / selected / excluded',
     f"{len(rows)+len(summary['excluded_rows'])} / {len(rows)} / {len(summary['excluded_rows'])}"],
    ['Selected initial molality range / mol kg⁻¹', f"{min(float(r['mi_mol_kg']) for r in rows)}–{max(float(r['mi_mol_kg']) for r in rows)}"],
    ['Selected final molality range / mol kg⁻¹', f"{min(float(r['mf_mol_kg']) for r in rows)}–{max(float(r['mf_mol_kg']) for r in rows)}"],
    ['Heat-capacity reference molality / mol kg⁻¹', identity['heat_capacity_reference_m_mol_kg']],
    ['Manufactured ideal thermal interval / K', str(identity['manufactured_thermal_interval_K'])]])
table('procedure', ['Numerical declaration', 'Value'], [
    *[[k.replace('_', ' '), f'{v:.8g}'] for k, v in identity['criteria'].items()],
    *[[setting+' '+key, f'{v:.8g}'] for setting, options in identity['settings'].items() for key, v in options.items()],
    ['Cartesian subinterval counts per coordinate', ', '.join(map(str, identity['grids']))],
    ['Solution / water seed temperature offsets / K', f"{identity['solution_seed_delta_K']} / {identity['water_seed_delta_K']}"],
    ['Thermal finite-difference steps / K', str(identity['finite_difference_steps_K'])],
    ['γ₁₆ extensive-operand rounding factor', f"{identity['gamma16']:.12g}"]])
table('summary', ['T / K', 'Rows / chords', 'Rows beyond / overlap / unresolved',
                  'Chords beyond / overlap / unresolved', 'Mean signed / J mol⁻¹', 'MAE / J mol⁻¹',
                  'RMSE / J mol⁻¹', 'Maximum absolute error / J mol⁻¹'],
      [['All' if g['T_K'] is None else number(g['T_K'], 2), f"{g['rows']} / {g['chords']}",
        ' / '.join(str(g['rows_'+k]) for k in ['disjoint', 'overlap', 'unresolved']),
        ' / '.join(str(g['chords_'+k]) for k in ['disjoint', 'overlap', 'unresolved']),
        *[number(g.get(k)) for k in ['mean_signed_error_J_mol', 'mean_absolute_error_J_mol',
                                    'RMSE_J_mol', 'max_absolute_error_J_mol']]] for g in summary['groups']])
table('signs', ['T / K', 'Positive errors', 'Negative errors', 'Benchmark / J mol⁻¹'],
      [[number(g['T_K'], 2), g['positive_errors'], g['negative_errors'],
        number(next(c['uniform_benchmark_J_mol'] for c in chords if c['T_K'] == g['T_K']), 6)]
       for g in summary['groups'][:-1]])
table('comparison', ['Source row', 'T / K', 'Initial → final / mol kg⁻¹', 'Measured / J mol⁻¹',
                     'Predicted / J mol⁻¹', 'Signed error / J mol⁻¹', 'Nonnegative gap / J mol⁻¹',
                     'Against study-chosen benchmark'],
      [[r['table']+':'+r['source_row'], number(r['T_K'], 2), r['mi_mol_kg']+' → '+r['mf_mol_kg'],
        *[number(r[k]) for k in ['observed_J_mol', 'predicted_J_mol', 'signed_nominal_error_J_mol',
                                'nonnegative_interval_gap_J_mol']],
        {'disjoint': 'exceeds 0.003RT', 'overlap': 'overlap', 'unresolved': 'unresolved'}[r['classification']]] for r in rows])
table('intervals', ['T / K', 'Initial → final / mol kg⁻¹', 'Expanded sampled interval / J mol⁻¹',
                    'b numerical / J mol⁻¹', 'b sampling / J mol⁻¹', 'Last grid n', 'ΔCp / J mol⁻¹ K⁻¹', 'Exact ∂TΔH / J mol⁻¹ K⁻¹'],
      [[number(c['T_K'], 2), f"{c['mi_mol_kg']} → {c['mf_mol_kg']}",
        '['+', '.join(number(v, 6) for v in c['interval_J_mol'])+']' if c['sampling_resolved'] else 'unresolved',
        f"{c.get('b_num_J_mol', 0):.3e}" if c['sampling_resolved'] else 'unavailable',
        number(c.get('b_sample_J_mol'), 6), c['levels'][-1]['n'] if c['levels'] else 'unavailable',
        number(c.get('delta_Cp_J_mol_K'), 6), number(c.get('exact_dH_dT_J_mol_K'), 6)] for c in chords])
table('heat-capacity', ['T / K', 'm / mol kg⁻¹', 'Cp difference from reference / J mol⁻¹ K⁻¹',
                        'Exact temperature action / J mol⁻¹ K⁻²'],
      [[number(r['T_K'], 2), r['m_mol_kg'], number(r.get('delta_Cp_from_1m_J_mol_K'), 6),
        number(r.get('exact_d_delta_Cp_dT_J_mol_K2'), 6)] for r in heat])
table('identities', ['T / K', 'Initial → final / mol kg⁻¹', 'Total − residual / J mol⁻¹',
                     'Changed ideal H / J mol⁻¹', 'Changed ideal Cp / J mol⁻¹ K⁻¹', '∂TH − ΔCp / J mol⁻¹ K⁻¹', 'Pass'],
      [[number(r['T_K'], 2), f"{r['mi_mol_kg']} → {r['mf_mol_kg']}",
        *[f'{r[k]:.3e}' if k in r else 'unavailable' for k in ['total_minus_residual_J_mol',
         'changed_record_H_J_mol', 'changed_record_Cp_J_mol_K', 'thermal_identity_J_mol_K']], r['passed']]
       for r in qualification['identities']])
table('thermal-ladder', ['Initial → final / mol kg⁻¹', 'Exact / J mol⁻¹ K⁻¹', 'FD coarse / J mol⁻¹ K⁻¹',
                         'FD fine / J mol⁻¹ K⁻¹', 'Fine − coarse / J mol⁻¹ K⁻¹', 'Exact − Richardson / J mol⁻¹ K⁻¹', 'Pass'],
      [[f"{r['mi_mol_kg']} → {r['mf_mol_kg']}", *[f'{r[k]:.9g}' if k in r else 'unavailable'
        for k in ['exact_J_mol_K', 'FD_h_0_1_K', 'FD_h_0_05_K', 'fine_minus_coarse_J_mol_K',
                  'exact_minus_Richardson_J_mol_K']], r['passed']] for r in qualification['thermal_ladder']])
diags = [s[k] for s in states for k in ['molality_diagnostics', 'temperature_diagnostics'] if k in s]
table('numerics', ['Quantity', 'Retained value'], [
    ['Native failures in final calculation', summary['native_failures']], ['State records', summary['state_records']],
    ['Cartesian sample records across refinements', summary['samples']],
    ['Maximum state pressure closure / Pa', f"{max(abs(s['pressure_error_Pa']) for s in states if 'pressure_error_Pa' in s):.6e}"],
    ['Maximum action pressure closure / Pa', f"{max(abs(d['pressure_error_Pa']) for d in diags):.6e}"],
    ['Minimum mechanical pressure slope / Pa m³ mol⁻¹', f"{min(s['pressure_density_slope'] for s in states if 'pressure_density_slope' in s):.6f}"],
    ['Maximum molality-action association diagnostic', f"{max(s['molality_diagnostics']['association_residual'] for s in states if 'molality_diagnostics' in s):.6e}"],
    ['Maximum temperature-action association diagnostic', f"{max(s['temperature_diagnostics']['association_residual'] for s in states if 'temperature_diagnostics' in s):.6e}"],
    ['Maximum empirical numerical expansion / J mol⁻¹', f"{max(c['b_num_J_mol'] for c in chords if c['sampling_resolved']):.6e}"],
    ['Maximum empirical sampling expansion / J mol⁻¹', f"{max(c['b_sample_J_mol'] for c in chords if c['sampling_resolved']):.6e}"]])
fig, axes = plt.subplots(1, 3, figsize=(12, 4.4), sharey=True, layout='constrained')
for ax, g in zip(axes, summary['groups'][:-1], strict=True):
    group = [r for r in rows if float(r['T_K']) == g['T_K']]
    cs = [c for c in chords if c['T_K'] == g['T_K'] and c['nominal_J_mol'] is not None]
    ax.errorbar([float(r['mi_mol_kg']) for r in group], [float(r['observed_J_mol']) for r in group],
        yerr=[float(r['benchmark_J_mol'])+float(r['rounding_J_mol']) for r in group],
        fmt='o', color='#0072B2', ms=4, capsize=2, label='Measured ± benchmark + rounding')
    for c in cs:
        yerr = [[c['nominal_J_mol']-c['interval_J_mol'][0]], [c['interval_J_mol'][1]-c['nominal_J_mol']]] if c['sampling_resolved'] else None
        ax.errorbar([c['mi_mol_kg']], [c['nominal_J_mol']], yerr=yerr, fmt='s',
                    color='#D55E00', mfc='white', ms=5, capsize=3)
    ax.plot([], [], 's', color='#D55E00', mfc='white', label='Prediction ± expanded sampled range')
    ax.axhline(0, color='0.65', lw=.7)
    ax.set(title=f"{g['T_K']:.2f} K", xlabel='Initial molality / mol kg⁻¹')
axes[0].set_ylabel('Dilution enthalpy / J mol⁻¹ NaCl')
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='outside lower center', ncol=2, fontsize=9)
fig.savefig(HERE/'figures/dilution-enthalpy.png', dpi=180)
plt.close(fig)

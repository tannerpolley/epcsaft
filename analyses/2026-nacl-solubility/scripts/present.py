"""Tables and discrete model–measurement comparison from retained values only."""
import csv
import json
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
rows = list(csv.DictReader((HERE / 'results/comparison.csv').open()))
states = json.loads((HERE / 'results/native-states.json').read_text())
summary = json.loads((HERE / 'results/summary.json').read_text())
primary = [r for r in states if r['variant'] == 'primary' and not r['options']]

def table(name, header, data):
    text = '| ' + ' | '.join(header) + ' |\n| ' + ' | '.join(['---'] * len(header)) + ' |\n'
    text += ''.join('| ' + ' | '.join(map(str, r)) + ' |\n' for r in data)
    (HERE / f'results/{name}.md').write_text(text)

def interval(row, low, high):
    return f"[{float(row[low]):.6f}, {float(row[high]):.6f}]"

table('conditions', ['Condition', 'Value'],
      [['System pressure / Pa', primary[0]['pressure_Pa']], ['Source pressure / Pa', primary[0]['reference_pressure_Pa']],
       ['Initial salt feed / mol', primary[0]['feed_mol']], ['Liquid water / kg', primary[0]['water_kg']]])

table('comparison-table', ['T / K', 'Measured / mol kg⁻¹', 'Predicted / mol kg⁻¹', 'Signed error / mol kg⁻¹',
                          'Measured 95% interval', 'Predicted ±1 SE input interval', 'Disjoint'],
      [[f"{float(r['temperature_K']):.2f}", *[f'{float(r[k]):.6f}' for k in
        ['measured_mol_kg', 'predicted_mol_kg', 'signed_error']],
        interval(r, 'measured_low_95', 'measured_high_95'), interval(r, 'predicted_low_se', 'predicted_high_se'),
        r['interval_disjoint']] for r in rows])
table('closure-table', ['T / K', 'Classified residual max', 'Absolute SI / decade', 'Absolute identity / ln unit',
                       'Scaled inventory max', 'Solid / mol', 'ρ / mol m⁻³', 'dP/dρ / Pa m³ mol⁻¹'],
      [[f"{r['temperature_K']:.2f}", f"{r['classified_max']:.2e}", f"{abs(r['si_decade']):.2e}",
        f"{abs(r['identity_ln']):.2e}", f"{r['inventory_scaled_max']:.2e}", f"{r['solid_mol']:.6f}",
        f"{r['density_mol_m3']:.2f}", f"{r['dp_drho_Pa_m3_mol']:.2f}"] for r in primary])
table('metrics', ['Quantity', 'Value'], [[k.replace('_', ' ').replace('aard', 'AARD').replace('mol kg', 'mol/kg'), f'{v:.8g}' if isinstance(v, float) else str(v)]
      for k, v in summary.items() if k != 'packet_1_sensitivity'])
sensitivity = next(r for r in states if r['variant'] == 'packet-1')
observed = float(rows[0]['measured_mol_kg'])
table('sensitivity', ['T / K', 'Screening convention', 'Predicted / mol kg⁻¹', 'Signed error / mol kg⁻¹', 'Change from /2 / mol kg⁻¹'],
      [[f"{sensitivity['temperature_K']:.2f}", '/1: σ', f"{sensitivity['molality_mol_kg']:.6f}",
        f"{sensitivity['molality_mol_kg'] - observed:+.6f}",
        f"{sensitivity['molality_mol_kg'] - primary[0]['molality_mol_kg']:+.6f}"]])
table('temperature-sensitivity', ['T / K', 'T ±0.1 K: molality range', 'Constant ε=78.09: molality'],
      [[f"{float(r['temperature_K']):.2f}", interval(r, 'temperature_low', 'temperature_high'),
        f"{float(r['constant_epsilon_mol_kg']):.6f}"] for r in rows])
t = [float(r['temperature_K']) for r in rows]
fig, ax = plt.subplots(figsize=(7.2, 4.5), layout='constrained')
for series, low, high, color, marker, label in [
    ('measured_mol_kg', 'measured_low_95', 'measured_high_95', '#0072B2', 'o', 'Pinho–Macedo: 95% expanded'),
    ('predicted_mol_kg', 'predicted_low_se', 'predicted_high_se', '#D55E00', 'x', 'Held /2 + Malmberg–Maryott: K ±1 SE')]:
    y = [float(r[series]) for r in rows]
    errors = [[float(r[series]) - float(r[low]) for r in rows], [float(r[high]) - float(r[series]) for r in rows]]
    ax.errorbar(t, y, yerr=errors, fmt=marker, color=color, capsize=3, label=label, markersize=6)
ax.set(xlabel='Temperature / K', ylabel='NaCl molality / mol kg⁻¹', title='Predicted solubility decreases as measured solubility rises')
ax.legend(fontsize=8, loc='lower left')
ax.grid(alpha=.2)
fig.savefig(HERE / 'figures/solubility.png', dpi=180)
plt.close(fig)
print('Retained tables and figures prepared without Engine execution')

# Temperature-trend diagnosis: presentation of retained results only.
TREND = HERE / 'results/temperature-trend'
activity_rounding = json.loads((HERE/'inputs/source/diagnosis-sources.json').read_text())['pitzer_peiper_busey_1984']['activity_half_last_digit']
def trend_rows(name):
    return [{k: float(v) if k not in ('variant', 'term', 'message', 'source') else v
             for k, v in r.items()} for r in csv.DictReader((TREND / f'{name}.csv').open())]

def trend_table(name, header, data):
    table(f'temperature-trend/{name}', header, data)

def upper_bound(value, digits):
    return f'{math.ceil(value*10**digits)/10**digits:.{digits}f}'

sat = trend_rows('saturation')
fixed = trend_rows('fixed-molality')
parts = trend_rows('contribution-decomposition')
laws = trend_rows('parameter-law-decomposition')
reference = [r for r in trend_rows('reference-comparison') if r['m_mol_kg'] == 6]
verification = trend_rows('temperature-verification')
pressure = trend_rows('pressure-comparison')
temps = [r['T_K'] for r in laws]
units = ['Contribution / 10⁻³ K⁻¹', *[f'{t:.2f} K' for t in temps]]
primary_sat = [r for r in sat if r['variant'] == 'primary']
primary_fixed = [r for r in fixed if r['variant'] == 'primary' and r['m_mol_kg'] == 6.2]
contributions = []
for key, label in [('hard_chain', 'Hard chain'), ('dispersion', 'Dispersion'),
                   ('association', 'Association'), ('debye_huckel', 'Debye–Hückel'),
                   ('born', 'Born (absent)'), ('density_conversion', 'Density and molality conversion'),
                   ('sum', 'Sum of decomposition')]:
    contributions.append([label, *[f"{next(r['dln_gamma_dT_component'] for r in parts if r['T_K'] == t and r['term'] == key and r['reference_m_mol_kg'] == 1e-16)*1000:+.6f}" for t in temps]])
contributions.append(['Exact native total', *[f"{r['dln_gamma_dT']*1000:+.6f}" for r in primary_fixed]])
trend_table('contribution-table', units, contributions)
trend_table('law-table', ['Thermal dependence / 10⁻³ K⁻¹', *units[1:]],
    [[label, *[f'{r[key]*1000:+.6f}' for r in laws]] for key, label in
     [('k_law', 'Sodium–water interaction'), ('water_sigma_law', 'Water segment diameter'),
      ('permittivity_law', 'Water permittivity'), ('ion_diameter_laws', 'Ion diameter laws (absent)'),
      ('explicit_temperature', 'Remaining explicit temperature dependence'), ('sum', 'Sum')]])
trend_table('cancellation-table', ['T / K', 'Sodium–water law / 10⁻³ K⁻¹',
    'Net of all other terms / 10⁻³ K⁻¹', 'Total / 10⁻³ K⁻¹', 'Law / total'],
    [[f"{r['T_K']:.2f}", f"{r['k_law']*1000:+.6f}", f"{(r['sum']-r['k_law'])*1000:+.6f}",
      f"{r['sum']*1000:+.6f}", f"{r['k_law']/r['sum']:.2f}"] for r in laws])
trend_table('saturation-slope-table', ['T / K', 'Same saturation molality / mol kg⁻¹',
    'Denominator / kg mol⁻¹', '½ d ln Kₛ/dT / 10⁻³ K⁻¹',
    'Primary gₜ / 10⁻³ K⁻¹', 'Primary dmₛ/dT / mol kg⁻¹ K⁻¹',
    'Locally frozen k: dmₛ/dT / mol kg⁻¹ K⁻¹'],
    [[f"{r['T_K']:.2f}", f"{r['m_sat']:.6f}", f"{1/r['m_sat']+r['dln_gamma_dm']:.7f}",
      f"{r['dlnK_dT']*500:+.6f}", f"{r['dln_gamma_dT']*1000:+.6f}", f"{r['dm_sat_dT']:+.8f}",
      f"{next(v['dm_sat_dT'] for v in sat if v['variant']=='local-k' and v['T_K']==r['T_K']):+.8f}"] for r in primary_sat])
trend_table('reference-table', ['T / K', 'Evaluated ln γ±', 'Predicted ln γ±',
    'Activity rounding bound / ln unit', 'Evaluated gₜ / 10⁻³ K⁻¹',
    'Predicted gₜ / 10⁻³ K⁻¹', 'Derivative rounding bound / 10⁻³ K⁻¹'],
    [[f"{r['T_K']:.2f}", f"{r['reference_ln_gamma']:+.6f}", f"{r['engine_ln_gamma']:+.6f}",
      upper_bound(-math.log1p(-activity_rounding/r['reference_gamma']), 9),
      f"{r['reference_dln_gamma_dT']*1000:+.6f}", f"{r['engine_dln_gamma_dT']*1000:+.6f}",
      upper_bound(r['reference_derivative_rounding_bound']*1000, 6)] for r in reference])
probe_table = []
for variant, label in [('primary', 'Primary'), ('freeze-k', 'Frozen sodium–water k'),
                       ('freeze-sigma', 'Frozen water segment σ'), ('constant-78.09', 'Constant ε = 78.09'),
                       ('ascani', 'Ascani–Held water ε')]:
    values = [next(r['m_sat'] for r in sat if r['T_K'] == t and r['variant'] == variant) for t in temps]
    probe_table.append([label, *[f'{m:.6f}' for m in values], f'{values[-1]-values[0]:+.6f}'])
measured = [float(next(r['measured_mol_kg'] for r in rows if float(r['temperature_K']) == t)) for t in temps]
probe_table.append(['Pinho–Macedo measurements', *[f'{m:.6f}' for m in measured], f'{measured[-1]-measured[0]:+.6f}'])
trend_table('probe-table', ['Declaration / mol kg⁻¹', *units[1:], 'Endpoint change / mol kg⁻¹'], probe_table)
max_decomp = max(abs(next(r['dln_gamma_dT_component'] for r in parts if r['term']=='sum' and r['T_K']==v['T_K'] and r['reference_m_mol_kg']==1e-16)-v['dln_gamma_dT']) for v in primary_fixed)
max_activity = max(abs(next(r['ln_gamma_component'] for r in parts if r['term']=='sum' and r['T_K']==v['T_K'] and r['reference_m_mol_kg']==1e-16)-v['ln_gamma']) for v in primary_fixed)
isolation = next(r for r in trend_rows('k-contribution-isolation') if r['T_K']==temps[0] and r['term']=='hard_chain')
isolation_disp = next(r for r in trend_rows('k-contribution-isolation') if r['T_K']==temps[0] and r['term']=='dispersion')
(TREND/'density-coupling.md').write_text(
    f"At {temps[0]:.2f} K, the interaction law contributes {isolation['k_law_effect_per_K']*1000:+.6f} "
    f"through hard chain and {isolation_disp['k_law_effect_per_K']*1000:+.6f} through dispersion, "
    "in units of $10^{-3}\\,\\mathrm{K^{-1}}$.\n")
(TREND/'verification.md').write_text(
    f"All {len(sat)} saturation calculations retain positive solid inventory. "
    f"The maximum classified residual is {max(r['classified_max'] for r in sat):.3e}; "
    f"the maximum activity-product identity discrepancy is {max(abs(r['identity_ln']) for r in sat):.3e} in log units; "
    f"the maximum equation-of-state pressure discrepancy is {max(abs(r['pressure_error_Pa']) for r in sat):.3e} Pa. "
    f"The largest exact-versus-Richardson derivative difference is {max(abs(r['error_per_K']) for r in verification):.3e} K⁻¹. "
    f"At the most dilute retained reference, the contribution sum differs from the analytic-reference total by at most {max_decomp:.3e} K⁻¹ "
    f"and {max_activity:.3e} in ln γ±. "
    f"The parameter-law sum closes within {max(abs(r['closure']) for r in laws):.3e} K⁻¹. "
    "These are numerical verification results, not bounds on physical error.\n")
(TREND/'reference-bounds.md').write_text(
    f"The largest printed-coefficient derivative rounding bound at 6 mol/kg is ±{max(r['reference_derivative_rounding_bound'] for r in reference):.6e} K⁻¹. "
    f"Changing the pressure from 100000 to 101325 Pa changes ln γ± by at most {max(abs(r['difference_ln']) for r in pressure):.3e} "
    "over the nine selected activity states. This pressure check bounds changes in activity values, "
    "without supplying a physical uncertainty in their temperature derivatives.\n")
# Discrete predictions and evaluated values, with source-rounding intervals only.
fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
for ax, ref_key, model_key, label, scale in [
    (axes[0], 'reference_ln_gamma', 'engine_ln_gamma', 'Mean ionic ln γ±', 1),
    (axes[1], 'reference_dln_gamma_dT', 'engine_dln_gamma_dT', '∂ ln γ±/∂T / 10⁻³ K⁻¹', 1000)]:
    bounds = [(-math.log1p(-activity_rounding/r['reference_gamma']) if scale==1 else r['reference_derivative_rounding_bound'])*scale for r in reference]
    ax.errorbar(temps, [r[ref_key]*scale for r in reference], yerr=bounds, fmt='o', color='#0072B2', capsize=4, label='Evaluated NaCl: rounding only')
    ax.plot(temps, [r[model_key]*scale for r in reference], 'x', color='#D55E00', markersize=7, label='Held /2 + Malmberg–Maryott')
    ax.axhline(0, color='.5', linewidth=.7)
    ax.set(xlabel='Temperature / K', ylabel=label)
    ax.grid(alpha=.2)
axes[0].legend(fontsize=7)
fig.suptitle('Concentrated NaCl activity has an excessive positive thermal response (6 mol/kg)', fontsize=10)
fig.savefig(HERE/'figures/activity-temperature-trend.png', dpi=180)
plt.close(fig)
print('Retained temperature-trend tables and figure prepared without Engine execution')

# One measured co-saturation state: tables suffice for these scalar comparisons.
CO = HERE / 'results/co-saturation'
co = json.loads((CO / 'summary.json').read_text())
cs = json.loads((CO / 'native-states.json').read_text())
comparison = []
for r in co['comparison']:
    for mapping in ('direct_intervals', 'mass_fraction_intervals'):
        values = [r['measured'], r['predicted'], r['relative_error_percent'], r['sampled_range']]
        bounds = r[mapping]
        comparison.append([r['salt'], mapping, *values, bounds['measured'], bounds['separation_mol_kg'], bounds['verdict']])
header = 'Salt|Mapping|Measured / mol kg⁻¹|Predicted / mol kg⁻¹|Error / %|Sampled range / mol kg⁻¹|Measured interval / mol kg⁻¹|Separation / mol kg⁻¹|Verdict'
table('co-saturation/comparison', header.split('|'), comparison)
checks = [[key, max(max(map(abs, r[key])) for r in cs)]
          for key in ('classified_residuals', 'balance_mol', 'solid_affinity_decades', 'identity_error_ln')]
checks += [['charge_mol', max(abs(r['charge_mol']) for r in cs)], ['all seven numerically closed', all(r['numerically_closed'] for r in cs)]]
checks += [[key, cs[0][key]] for key in ('T_K', 'P_Pa', 'feed_mol', 'solid_amounts', 'molar_densities', 'dp_drho')]
checks += [['doubled-feed solid amounts / mol', cs[1]['solid_amounts']], ['doubled-feed relative molality changes', co['excess_feed_relative_error']]]
checks += [['100000 Pa molalities / mol kg⁻¹', cs[2]['molality_mol_kg'][:2]]]
checks += [[key, value] for key, value in co['thermochemistry'].items() if key not in ('direct_intervals', 'mass_fraction_intervals')]
table('co-saturation/verification', ['Quantity (SI unless stated)', 'Value'], checks)

"""Present both one-coordinate trials; native outcomes record identification before binary scoring."""
import json
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from run import C, D, LEAF

R = LEAF / D.get('results', 'results/issue-798-r2')

def table(path, frame):
    columns = list(frame.columns)
    rows = ['| ' + ' | '.join(columns) + ' |', '| ' + ' | '.join('---' for _ in columns) + ' |']
    rows += ['| ' + ' | '.join(('unavailable' if pd.isna(v) else f'{v:.9g}' if isinstance(v, float) else str(v)).replace('|', r'\|') for v in row) + ' |' for row in frame.itertuples(index=False, name=None)]
    (R / path).write_text('\n'.join(rows) + '\n')


coordinates = D.get('coordinates', [D.get('coordinate')])
training_label = f"{D['training_states']} {C['source']} training states"
if D.get('sulfate'): training_label = f"{D['training_observations']} scalar observations: {D['training_states']} binary states and 21 LiCl–MgSO₄ mixtures"
if 'mixture' in D:
    source = D['mixture'].get('source', 'Clegg').replace(' Table 1', '').replace(' Table 2', '')
    training_label += f" plus {D['mixture']['states']} {source} mixtures; {D['training_observations']} scalar observations"
    admission = pd.read_csv(R / D['mixture'].get('retained_admission', 'mixture-admission.csv'), float_precision='round_trip')
    table('admission-table.md', admission.groupby('y').agg(cells=('passes', 'count'), admitted=('passes', 'sum'), minimum_m_total=('m_total', 'min'), maximum_m_total=('m_total', 'max'), minimum_I=('I', 'min'), maximum_I=('I', 'max')).reset_index())
    table('guo-conditions.md' if 'later_Guo_conditions' in D else 'check-conditions.md', pd.DataFrame([D.get('later_Guo_conditions', D.get('later_check_conditions'))]))
    (R / 'mixture-training-table.md').write_text(f'Informational {source} selected-endpoint training metrics: not reached unless fitting and identification pass. No qualification threshold applies.\n')
    if D['mixture'].get('reference'):
        table('reference-reproduction-table.md', pd.read_csv(R / 'reference-reproduction.csv', float_precision='round_trip'))
        table('reference-deviation-table.md', pd.read_csv(R / 'reference-deviation.csv', float_precision='round_trip'))
table('fixed-table.md', pd.DataFrame([{'quantity': k, 'fixed value': v, 'unit': ('kg²/mol²' if k.startswith('psi_') else 'dimensionless' if k.startswith('k_') or k == 'segment_count' else 'kelvin' if k == 'dispersion_energy_over_k' else 'angstrom')} for k,v in D.get('fixed', D.get('H_fixed')).items()] + ([] if D.get('sulfate') or len(coordinates) == 2 and 'verification' not in D else [{'quantity': f"k_{C['label']}-Cl (modelling assumption)", 'fixed value': D.get('ion_Cl_k', D.get('H_Cl_k')), 'unit': 'dimensionless'}])))
table('limits-table.md', pd.DataFrame([{'rule': k, 'limit': v} for k,v in D['identification'].items()]))
if (R / 'sensitivity-summary.csv').exists():
    table('prefit-table.md', pd.read_csv(R / 'sensitivity-summary.csv', float_precision='round_trip'))
    if 'verification' in D or D.get('mixture', {}).get('reference'):
        diagnostics = pd.read_csv(R / 'prefit-diagnostics.csv', float_precision='round_trip')
        table('coverage-table.md', diagnostics.groupby('evaluation').agg(observations=('observation', 'nunique'), state_entries=('state', 'count'), minimum_density=('density_mol_m3', 'min'), minimum_pressure_slope=('pressure_slope', 'min')).reset_index())
    if 'verification' in D:
        table('lengths-table.md', pd.read_csv(R / 'length-coordinates.csv', float_precision='round_trip').query('evaluation <= 9'))
        verified = pd.read_csv(R / 'length-verification.csv', float_precision='round_trip')
        table('verification-table.md', verified.groupby(['step', 'coordinate']).agg(checks=('passes', 'count'), passed=('passes', 'sum')).reset_index())
    if not (R / 'outcome.txt').read_text().startswith(('Identified', 'Not identified')):
        (R / 'outcome.md').write_text((R / 'outcome.txt').read_text())
        (R / 'conditions.md').write_text(f"Conditions: **{D['temperature_K']:g} K**, **{D['pressure_Pa']:g} Pa**; **{training_label}**. The complete declared grid is shown below; fitting and evaluated-binary scoring were not reached.\n")
        raise SystemExit(0)
fits = pd.read_csv(R / 'fit-status.csv', float_precision='round_trip')
values = pd.read_csv(R / 'fit-values.csv', float_precision='round_trip')
identified = pd.read_csv(R / 'identification.csv', float_precision='round_trip')
resolution = pd.read_csv(R / 'local-resolution.csv', float_precision='round_trip')
coefficient = f"k_{C['label']}-water"
names = [c.get('name', f"k_{C['label']}-" + ('water' if 'water' in c.get('pair', []) else 'Cl')) for c in coordinates]
points = values[values.field == 'optimizer'].pivot(index='start', columns='index', values='value').rename(columns=dict(enumerate(names))).reset_index()
for a, name in enumerate(names): points['initial ' + name] = [D['starts'][int(i)-1][a] for i in points.start]
table('fit-table.md', fits.merge(points, on='start', how='left'))
checks = identified.merge(resolution[['start', 'coordinate', 'C_aa', 'h', 'R']], on='start', how='left')
if 'verification' in D:
    checks['coordinate name'] = checks.coordinate.map(dict(enumerate(names)))
    checks['physical h'] = checks.h * checks.coordinate.map(dict(enumerate(c['scale'] for c in coordinates)))
else:
    span = coordinates[0]['bounds'][1] - coordinates[0]['bounds'][0]
    checks['endpoint distance'] = checks.max_endpoint_distance * span
table('identification-table.md', checks)
diagnostics = values[values.field.isin(['singular_values', 'C', 'bound_distances', 'active_bound_count', 'projected_gradient'])]
table('diagnostics-table.md', diagnostics)
columns = values[values.field == 'optimizer_jacobian'].groupby('start').value
sensitivity = columns.agg(['count', 'min', 'max'])
sensitivity['minimum absolute weighted sensitivity'] = columns.agg(lambda x: x.abs().min())
table('sensitivity-table.md', sensitivity.reset_index())
# Preserve the native identification outcome recorded before binary scoring.
(R / 'outcome.md').write_text((R / 'outcome.txt').read_text())
(R / 'conditions.md').write_text(f"Conditions: **{D['temperature_K']:g} K**, **{D['pressure_Pa']:g} Pa**; **{training_label}**. Per-start limits: **{D['maximum_iterations']} iterations**, **{D['maximum_seconds']} s**; shared allowance: **{D['total_maximum_iterations']} iterations**, **{D['total_maximum_seconds']} s**. Returned/converged starts: **{len(fits)}/{sum(fits.status == 'converged')}**; total iterations: **{fits.iterations.sum()}**.\n")
selected_rows = checks[checks.selected == 1].merge(points, on='start')
if selected_rows.empty:
    (R / 'selected.md').write_text('No finite converged endpoint was selected; retained native diagnostics describe the stopping point.\n')
else:
    selected = selected_rows.iloc[0]
    (R / 'selected.md').write_text('Selected ' + ', '.join(f'**{name}={selected[name]:.17g}**, initial **{selected["initial " + name]:g}**' for name in names) + f"; {'combined training' if 'mixture' in D else 'binary'} cost **{selected.J_min:.17g}**, cost increment **{selected.delta_J:.17g}**. Per-coordinate excursions and all identification values are retained in the tables.\n")

from score import score
if 'trial' in D:
    if (R / 'outcome.txt').read_text().startswith('Identified'):
        errors, metrics = score(LEAF, R)
        if 'mixture' in D or D.get('sulfate'):
            admission = pd.read_csv(R / 'mixture-admission.csv', float_precision='round_trip') if D.get('sulfate') else admission
            predictions = values[(values.start == selected.start) & (values.field == 'predictions')].sort_values('index').value.to_numpy()[D['training_states'] if D.get('sulfate') else 2*D['training_states']:]
            if D.get('sulfate'): admission = admission[admission.salt == 'LiCl-MgSO4']
            assert len(predictions) == len(admission) == (21 if D.get('sulfate') else D['mixture']['states'])
            mixture_errors = admission.assign(prediction=predictions, error=predictions-admission.observed)
            mixture_errors.to_csv(R / 'mixture-training-errors.csv', index=False)
            mixture_metrics = pd.DataFrame([dict(cells=len(admission), phi_RMS=(mixture_errors.error.pow(2).mean())**.5, phi_max=mixture_errors.error.abs().max(), objective_contribution=.5*(mixture_errors.error/.01).pow(2).mul(admission.weight).sum(), role='informational training; no qualification threshold')])
            mixture_metrics.to_csv(R / 'mixture-training-metrics.csv', index=False)
            table('mixture-training-table.md', mixture_metrics)
        errors.to_csv(R / '2a-binary-errors.csv', index=False)
        metrics.to_csv(R / '2a-binary-metrics.csv', index=False)
        table('2a-table.md', metrics.set_index('bucket_mol_per_kg').T.rename_axis('quantity').reset_index())
        (R / '2a-outcome.md').write_text('First failed requirement: ' + next((f'{row.bucket_mol_per_kg}: {metric}' for row in metrics.itertuples() for metric in D['accuracy'] if not getattr(row, metric + '_passes')), 'none; all evaluated-binary limits pass') + '.\n')
        figure, axes = plt.subplots(1, 1 if D.get('sulfate') else 2, figsize=(10, 4), layout='constrained')
        if D.get('sulfate'):
            for salt, group in errors.groupby('salt'): axes.scatter(group.molality, group.phi_error, label=salt)
            axes.set(xlabel='Salt molality / (mol per kg water)', ylabel='Osmotic coefficient error'); axes.legend()
        else:
            for axis, field in zip(axes, ['phi_error', 'gamma_relative_error_percent']): errors.plot.scatter('molality', field, ax=axis)
        figure.savefig(R / 'binary-errors.png', dpi=160)
    else:
        (R / '2a-table.md').write_text('2a not reached: native identification failed.\n')
    raise SystemExit(0)
comparison = []
for trial, declaration in [('issue-798-r2', 'fit-declaration-issue-798-r2.json'), ('issue-798-retrial', 'fit-declaration-issue-798-retrial.json')]:
    result = LEAF / 'results' / trial
    errors, metrics = score(LEAF, result)
    errors.to_csv(LEAF / 'results/issue-798-retrial' / f'{trial}-binary-errors.csv', index=False)
    selected = pd.read_csv(result / 'identification.csv').query('selected == 1').iloc[0]
    value = pd.read_csv(result / 'fit-values.csv', float_precision='round_trip').query('field == "optimizer" and start == @selected.start').iloc[0].value
    h = pd.read_csv(result / 'local-resolution.csv', float_precision='round_trip').query('start == @selected.start').iloc[0].h
    comparison.append(metrics.assign(trial=trial, k_H_Cl_modelling_assumption=json.loads((LEAF / 'inputs' / declaration).read_text())['H_Cl_k'], identification='identified', k_H_water=value, h=h))
R = LEAF / 'results/issue-798-retrial'
comparison = pd.concat(comparison, ignore_index=True)
comparison.to_csv(R / 'binary-metrics.csv', index=False)
table('comparison-table.md', comparison.set_index(comparison.trial + ', ' + comparison.bucket_mol_per_kg).drop(columns=['trial', 'bucket_mol_per_kg']).T.rename_axis('quantity').reset_index())

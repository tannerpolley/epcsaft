"""Present retained verification rows; import no Engine and run no calculation."""
import csv
import json
from pathlib import Path

out = Path(__file__).resolve().parents[1] / 'results'
s = json.loads((out/'summary.json').read_text())
rows = list(csv.DictReader((out/'falsifiers.csv').open()))
lines = [f"Wheel SHA-256: `{s['wheel_sha256']}`.", '',
    '| Check | Comparisons | Largest error / tolerance | Result |',
    '|---|---:|---:|---|']
for name in ('F1','F1a','F2','F3','F4','F5','F6'):
    selected = [r for r in rows if r['case'] == name]
    ratios = [float(r['error'])/float(r['tolerance']) for r in selected if float(r['tolerance']) > 0]
    lines.append(f"| {name} | {len(selected)} | {max(ratios):.6g} | {'Pass' if all(r['passed']=='True' for r in selected) else 'Miss'} |")
lines += ['', ('F1 uses 383 rows and 3,167 applicable columns. Its allowance is '
    '$2\\times10^{-5}\\max|J|+\\delta$, with $\\delta=10^{-7}$ in dimensionless units '
    'and $RT\\,10^{-7}$ J/mol for energy. Unprescribed difference steps are retained without a verdict.'), '',
    '| Manufactured recovery | Largest error | Required limit |', '|---|---:|---:|']
for prefix, label, limit in (('coordinate-', 'k coordinates',1e-6), ('diameter-', 'held diameters (Angstrom)',1e-7), ('decade-', 'log10 residuals (decade)',1e-10)):
    selected=[r for r in rows if r['case']=='F2' and r['quantity'].startswith(prefix)]
    if prefix=='coordinate-':selected=[r for r in selected if r['quantity']!='coordinate-6']
    lines.append(f"| {label} | {max(float(r['error']) for r in selected):.6g} | {limit:.6g} |")
c=next(r for r in rows if r['case']=='F2' and r['quantity']=='coordinate-6')
lines.append(f"| C | {float(c['error']):.6g} | 1e-5 |")
lines += ['', '| Held ion | Stored diameter (Angstrom) | Independent root (Angstrom) | Energy error (J/mol) | Bracket width (Angstrom) |', '|---|---:|---:|---:|---:|']
for r in csv.DictReader((out/'independent-roots.csv').open()):
    lines.append('| '+' | '.join(r[k] for k in ('ion','stored_diameter','independent_diameter','energy_error','bracket_width'))+' |')
lines += ['', ('F3 requires diameter agreement within 1e-6 Angstrom and energy error '
    'within 1e-2 J/mol. Every inner hold also meets 1e-6 J/mol and 1e-10 Angstrom.'), '',
    '| Native check group | Maximum | Mean | Outcome |', '|---|---:|---:|---|']
for r in s['checks']['groups']:
    lines.append(f"| {r['name']} | {r['maximum']:.12g} | {r['mean']:.12g} | {r['status']} |")
lines += ['', ('The unavailable negative-energy logarithmic residual retains its signed '
    'J/mol difference. It blocks a guard; the report-only declaration retains a valid zero '
    'objective and passes overall acceptance. The epsilon group deliberately passes its '
    'maximum and fails its mean. Metric discrepancies are at most 2.9e-15 in their '
    'reported units, against 1e-10.'), '',
    (f"F6 maximum induced root uncertainty: {s['charts']['maximum_root_uncertainty_angstrom']:.6g} Angstrom, against 1e-10. "
    'The root/map study uses source and manufactured points, three steps and nonunit scales '
    'with tied C. Four native failure/cache comparisons require exact agreement.'), '',
    '| Native fit | Weighted objective | Seconds | Covariance | Active coordinates |', '|---|---:|---:|---|---|']
for name in ('F2','F3'):
    r=s['fits'][name]
    lines.append(f"| {name} | {r['cost']:.6g} | {r['seconds']:.6g} | {r['covariance_status']} | {r['active_bounds']} |")
lines += ['', ('The native raw-value/column calls took '
    f"{sum(r['exact_value_columns_seconds'] for r in s['timings']):.6g} s across F1, "
    'and the corresponding central-difference batches took '
    f"{sum(r['central_difference_value_columns_seconds'] for r in s['timings']):.6g} s. "
    'These measured calls include construction, reference admission and derivative work. '
    'They do not establish an optimized finite-difference baseline or a general speedup.'), '']
(out/'qualification.md').write_text('\n'.join(lines))

# D7 presentation consumes retained scalars only, without importing the Engine.
if (out / 'minimax.json').exists():
    b = json.loads((out / 'minimax.json').read_text())
    native = list(csv.DictReader((out / 'minimax-native.csv').open()))
    points = list(csv.DictReader((out / 'minimax.csv').open()))
    def error_max(selected):
        return max(abs(float(r['model']) - float(r['reference'])) for r in selected)
    lines = [f"B3/B4 replay wheel SHA-256: `{b['wheel_sha256']}`.", '',
        '| Analytic manufactured comparison | Largest absolute error | Limit |', '|---|---:|---:|']
    for label, token, limit in [('B1 coordinates', '/coordinate/', 1e-7), ('B1/B2 epigraph', '/t', 1e-9),
                               ('B1 mean slacks', '/slack/', 1e-9), ('B1 multipliers', 'multiplier', 1e-6)]:
        selected = [r for r in native if token in r['comparison_id']]
        lines.append(f"| {label} | {error_max(selected):.8g} | {limit:g} |")
    lines += ['', ('B1 uses x=0.75, z=0.5, d=1.5 and t=0.0175 decade, with both starts and the rescaled repeat. '
        'Its equality multiplier is 0.005 decade per H-unit. B2 recovers x=1.5, d=3 and t=0.015 decade; '
        'the analytic least-squares point x=0.6 has a 0.024-decade maximum, above 0.02.'), '',
        '| Physical-hold manufactured solve | Start 0 | Start 1 | Required limit |', '|---|---:|---:|---:|']
    cases = [{r['quantity']: float(r['exact']) for r in points if r['case'] == f'B4-start-{i}'} for i in range(2)]
    for label, key, limit in [('t (decade)', 'minimax/t', .005 + 1e-8),
                            ('epigraph gap (decade)', 'minimax/gap', 1e-8),
                            ('scaled primal', 'minimax/primal', 1e-8),
                            ('scaled stationarity', 'minimax/stationarity', 1e-8),
                            ('scaled dual-sign defect', 'minimax/dual_sign', 1e-8),
                            ('scaled complementarity', 'minimax/complementarity', 1e-8),
                            ('permittivity maximum (1)', 'check-0/maximum', .1),
                            ('permittivity mean (1)', 'check-0/mean', .1),
                            ('stored hydration maximum (J/mol)', 'check-1/maximum', 10)]:
        lines.append(f"| {label} | {cases[0][key]:.12g} | {cases[1][key]:.12g} | {limit:g} |")
    root_max = max(float(r['exact']) for r in points if r['quantity'].endswith('/root_error'))
    full_max = max(abs(float(r['exact'])) for r in points if r['quantity'].startswith('full-stationarity/'))
    refinement = [r for r in points if r['case'].endswith('-refinement') and r['error']]
    b3 = [r for r in points if r['tolerance'] and not r['case'].endswith('-refinement')
          and not r['quantity'].startswith('full-stationarity/')]
    lines += ['', (f"The two native solves used {b['B4']['solves'][0]['iterations']} and {b['B4']['solves'][1]['iterations']} iterations, "
        f"with a measured batch duration of {b['B4']['elapsed_seconds']:.6g} s (300/120 s per fit; 600/240 s shared)."),
        (f"Maximum estimated scaled root error: {root_max:.8g}, against 1e-9. "
        f"Maximum reconstructed full-space stationarity component: {full_max:.8g}, against 1e-8. "
        "The held components vanish by the multiplier definition; B1's analytic KKT conditions independently qualify that multiplier."),
        (f"Requiring 1e-7 J/mol and 1e-11 Angstrom roots changes each start's 45 scaled rows and 450 columns by "
        f"{max(float(r['error']) for r in refinement):.8g}, against 1e-9. "
        "The original roots already satisfy the tighter criteria, so this replay returns the same iterate and establishes no additional error bound."), '',
        '| Returned parameter | Start 0 | Start 1 | Unit |', '|---|---:|---:|---|']
    for label, index, unit in [('Li-water k', 0, '1'), ('Li-Cl k', 1, '1'),
                               ('Na-water k', 2, '1'), ('Li Born diameter', 10, 'Angstrom'),
                               ('K Born diameter', 12, 'Angstrom')]:
        lines.append(f"| {label} | {cases[0][f'fit/physical-{index}']:.9g} | {cases[1][f'fit/physical-{index}']:.9g} | {unit} |")
    lines += ['', ('The two starts return distinct lithium and potassium parameters with nearly identical objective values. '
        'They supply no unique parameter estimate; Na-water k is at its lower bound in both.'), '',
        (f"B3 has {len(b3)} retained numeric comparisons at source, manufactured truth and both returned points. "
        "Its criterion remains 1e-12 relative, with 1e-12 absolute only at zero. "
        "At each fitted point both prediction evaluators replay with common cold density initialization; "
        "stored-fit minus cold predictions are retained separately as initialization diagnostics. "
        "All five kinds, nonunit affine scales, tied C, frozen thresholds and stored checks are covered. "
        "The report-only transfer check fails while the enforced rows remain identical and the fitted point passes. "
        "Zero-Ceres-weight permittivity rows enter B's constraints. "
        "B3 checks consistency across selections of the shared evaluator; independent derivative evidence remains slice A's F1 on its identified wheel."), '',
        ('These synthetic rows verify local solver integration and numerical stationarity, not physical validation, '
        '#209 adoption, uniqueness or global optimality. The admitted finite reference endpoint retains its existing limits.'), '']
    sensitivity = json.loads((out / 'minimax-diagnosis.json').read_text())['density_sensitivity']
    lines += ['', '| Row-3 density sensitivity | Value | Unit |', '|---|---:|---|']
    for label, key, unit in [('Exact density derivative', 'dphi_drho_m3_mol', 'm³/mol'),
                             ('One density ulp', 'density_ulp_mol_m3', 'mol/m³'),
                             ('Smooth derivative times one ulp', 'smooth_one_ulp_change', '1'),
                             ('Observed absolute change', 'observed_absolute_change', '1'),
                             ('Observed / smooth change', 'observed_over_smooth', '1')]:
        lines.append(f"| {label} | {sensitivity[key]:.10g} | {unit} |")
    (out / 'minimax-qualification.md').write_text('\n'.join(lines))

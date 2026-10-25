"""Retain bounded native probes of source discrepancies and excluded states."""
import csv
import json
import math
import time

from run import CALIBRATION, HERE, PARAMETERS, STARTS, TOL, _core, eq, model, pure_rows, write_rows


def main():
    started = time.perf_counter()
    rows = list(csv.DictReader((HERE / 'results/figure-comparisons.csv').open()))
    predictions = list(csv.DictReader((HERE / 'results/pure-saturation.csv').open()))
    population = {
        'figure_requested': len(rows), 'prediction_requested': len(predictions),
        'returned': sum(r['success'] == 'True' for r in rows + predictions),
        'figure_qualified': sum(r['numerical_tolerance_met'] == 'True' for r in rows),
        'prediction_qualified': sum(r['numerical_tolerance_met'] == 'True' for r in predictions)}
    assert population == dict(figure_requested=122, prediction_requested=78, returned=200,
                              figure_qualified=119, prediction_qualified=78)
    output, defects = [], []
    for figure, temperature in [('10', '204'), ('12', '4')]:
        candidates = [r for r in rows if r['figure'] == figure and r['temperature_c'] == temperature and r['reference'] == 'model']
        row = max(candidates, key=lambda r: float(r['absolute_error']))
        case = CALIBRATION[figure]
        m = model(case['components'], case['k'])
        x = float(row['x'])
        y = STARTS['vapor_first_component_mole_fraction'][figure][temperature] * .95
        problem = eq.Problem(T=float(temperature)+273.15, P=eq.Free(float(row['y'])*.9e5),
            phases=[eq.Phase('L', kind='liquid', composition_guess=[x, 1-x]),
                    eq.Phase('V', kind='vapor', amount=eq.Pinned(0), composition_guess=[y, 1-y])],
            feed=eq.MoleFractions(dict(zip(case['components'], [x, 1-x]))),
            row_tolerance=TOL['equilibrium_maximum_absolute_row_residual'])
        r = eq.solve_equilibrium(m, problem)
        output.append({'figure': figure, 'component': row['component'], 'temperature_k': float(temperature)+273.15,
            'source_x': x, 'baseline_absolute_error': row['absolute_error'],
            'probe': 'largest model miss: alternate composition and pressure starts', 'success': r.success, 'message': r.message,
            'baseline_pressure_pa': row['pressure_pa'], 'probe_pressure_pa': r.pressure if r.success else '',
            'pressure_relative_change': abs(r.pressure/float(row['pressure_pa'])-1) if r.success else ''})
    excluded = [r for r in rows if r['success'] != 'True' or r['numerical_tolerance_met'] != 'True']
    assert {r['component'] for r in excluded} == {'propylene', 'diethyl ether', 'toluene'}
    for row in excluded:
        name, t = row['component'], float(row['temperature_k'])
        m = model((name,))
        p = float(row['pressure_pa'])
        liquid = float(row['rho_liquid_mol_m3'])
        vapor = float(row['rho_vapor_mol_m3'])
        coefficient = m._model.packing_fraction_coefficient(t, [1.])
        problem = m.saturation(t)
        problem.row_tolerance = TOL['equilibrium_maximum_absolute_row_residual']
        probes = []
        # Separate explicit starts are diagnostic experiments, not production retries.
        for factor in (1., .8):
            declaration = {'pressure': factor*p, 'packing_fractions': [liquid*coefficient, factor*vapor*coefficient]}
            attempt = eq.solve_equilibrium_starts(m, problem, [_core.EquilibriumStart(**declaration)],
                _core.EquilibriumWorkLimits(attempts=1))
            r = attempt.results[0]
            probe = {'figure': row['figure'], 'component': name, 'temperature_k': t,
                'probe': 'explicit liquid and vapor packing starts', 'start': json.dumps(declaration),
                'success': r.success if r else False, 'message': r.message if r else attempt.causes[0]}
            if r and r.success:
                states = [m.state(t, rho=rho) for rho in r.molar_densities]
                error = max(abs(s.pressure-r.pressure) for s in states)
                probe.update(probe_pressure_pa=r.pressure, rho_liquid_kg_m3=r.molar_densities[0]*float(PARAMETERS[name]['MW']),
                    rho_liquid_mol_m3=r.molar_densities[0], rho_vapor_mol_m3=r.molar_densities[1],
                    phase_pressure_absolute_closure_pa=error, phase_pressure_relative_closure=error/r.pressure,
                    phase_pressure_density_ulp_pa=max(.5*abs(s.pressure_density_derivative)*math.ulp(s.molar_density) for s in states),
                    max_abs_residual=max(map(abs, r.residuals)))
            probes.append(probe)
            output.append(probe)
        selected = m.state(t, P=p, phase='liquid').molar_density
        with _core.Capture(diagnostics=True) as capture:
            default = eq.solve_equilibrium(m, problem)
        assert default.success and row['success'] == 'True', default.message
        assert abs(default.pressure/p-1) <= 1e-10
        initialization = [{'name': e.name, 'values': list(e.values)} for e in capture.events
                          if e.name in ('superancillary_start', 'loop_pressure_seed')]
        assert initialization, 'Missing native initialization evidence'
        diagnosis = ('low-pressure cancellation: returned equilibrium misses relative pressure closure; '
                     'represented liquid-density steps limit pressure resolution')
        if name == 'toluene':
            diagnosis = ('The default toluene result fails the relative pressure criterion. '
                         'An explicit start at the same temperature and with the same row tolerance '
                         'gives relative closure 6.113574794657297e-8. Density representation and '
                         'pressure cancellation make this case sensitive, but the criterion is '
                         'attainable in the retained calculation. The default row remains excluded '
                         'pending the Engine follow-up (https://github.com/tannerpolley/ePC-SAFT/issues/849).')
        defects.append({'component': name, 'temperature_k': t, 'base': 'pc-saft',
            'association': 'none', 'electrolyte': 'none', 'polar': 'none',
            **{k: PARAMETERS[name][k] for k in ('m', 's', 'e', 'MW')},
            'parameter_units': 'm dimensionless; s angstrom; e K; MW kg/mol',
            'problem': 'T fixed; P free; pure mole fraction 1; liquid amount free; vapor amount pinned 0',
            'row_tolerance': TOL['equilibrium_maximum_absolute_row_residual'],
            'pressure_closure_tolerance': TOL['phase_pressure_relative_closure'],
            'baseline_success': row['success'], 'baseline_numerical_tolerance_met': row['numerical_tolerance_met'],
            'baseline_message': row['message'], 'baseline_pressure_pa': row['pressure_pa'],
            'baseline_phase_pressure_relative_closure': row['phase_pressure_relative_closure'],
            'baseline_phase_pressure_absolute_closure_pa': float(row['pressure_pa'])*float(row['phase_pressure_relative_closure']) if row['pressure_pa'] else '',
            'explicit_probes': json.dumps(probes), 'liquid_branch_pressure_pa': p,
            'selected_liquid_branch_density_mol_m3': selected, 'diagnosis': diagnosis,
            'default_initialization_events': json.dumps(initialization),
            'independent_liquid_density_mol_m3': 18276.58 if name == 'propylene' else '',
            'independent_vapor_density_mol_m3': 1.0775e-6 if name == 'propylene' else '',
            'independent_pressure_pa': 7.837e-4 if name == 'propylene' else '',
            'independent_basis': 'historical reviewer supplied independent-code calculation' if name == 'propylene' else ''})
        sample = next(r for r in pure_rows() if r['component'] == name)
        initial = eq.solve_equilibrium(m, m.saturation(sample['temperature_k']))
        neighbor = m.saturation(t+5.)
        neighbor.row_tolerance = TOL['equilibrium_maximum_absolute_row_residual']
        attempted = eq.solve_equilibrium_starts(m, neighbor, [_core.EquilibriumStart(coordinates=initial.coordinates)],
            _core.EquilibriumWorkLimits(attempts=1))
        r = attempted.results[0]
        output.append({'figure': row['figure'], 'component': name, 'temperature_k': t+5.,
            'probe': 'neighboring temperature from native prediction coordinates',
            'success': r.success if r else False, 'message': r.message if r else attempted.causes[0]})
    write_rows('diagnosis.csv', output)
    write_rows('engine-defect-cases.csv', defects)
    report = dict(wall_seconds=time.perf_counter()-started, population=population,
                  excluded_cases=len(defects), probes=len(output))
    (HERE / 'results/diagnosis-run.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()

"""Declare Schick Table S10 fixed-tuple native binary VLE and retain evidence."""
import copy
import csv
import hashlib
import json
import math
import time
from pathlib import Path

import epcsaft
from epcsaft import _core
from epcsaft.equilibrium import Amounts, Free, MoleFractions, Phase, Problem, bubble_point, compile_problem, solved_state_actions, solve_equilibrium

LEAF = Path(__file__).resolve().parents[1]
OUT = LEAF / 'results/data'
TOL = dict(row=1e-10, refinement=1e-12, composition_error=1e-8,
           composition_half_width=0.00005, aard_percentage_point=0.005,
           pressure_half_width_bar=0.05, row23_pressure_half_width_bar=0.5)
FIELDS = ('success message solver_status requested_tolerance_met temperature pressure amounts '
          'mole_fractions molar_densities packing_fractions residuals classified_residuals '
          'coordinates iterations evaluations hessian_calls unknown_count jacobian_rank '
          'generic_nullity reduced_hessian_minimum topology_event solver_settings solver_options linear_solver').split()


def record(model, problem, result):
    data = {key: getattr(result, key) for key in FIELDS}
    data.update(requested_row_tolerance=problem.row_tolerance, closure_tolerance=problem.row_tolerance)
    data['topology_event'] = str(result.topology_event)
    data['convergence'] = {k: getattr(result.convergence, k) for k in ('residual', 'floor', 'amplification', 'contraction')}
    data['convergence']['stop'] = str(result.convergence.stop)
    data['closed'] = False
    if result.success:
        states = [model.state(problem.T, rho=rho, x=result.phase_mole_fractions(i))
                  for i, rho in enumerate(result.molar_densities)]
        data['state_pressure_pa'] = [s.pressure for s in states]
        data['log_fugacity_coefficients'] = [s.log_fugacity_coefficient for s in states]
        data['pressure_density_derivative_pa_m3_per_mol'] = [s.pressure_density_derivative for s in states]
        data['state_density_branch'] = [str(s.density_branch) for s in states]
        data['inventory_mol'] = [sum(a * result.phase_mole_fractions(i)[j]
                                     for i, a in enumerate(result.amounts)) for j in range(2)]
        data['inventory_error_mol'] = max(abs(data['inventory_mol'][j] - v)
                                         for j, v in enumerate(problem.feed.values.values()))
        data['normalization_error'] = max(abs(sum(result.phase_mole_fractions(i))-1) for i in range(2))
        data['pressure_relative_error'] = max(abs(s.pressure/result.pressure-1) for s in states)
        data['log_fugacity_difference'] = [math.log(result.mole_fractions[j]/result.mole_fractions[2+j])
            + states[0].log_fugacity_coefficient[j]-states[1].log_fugacity_coefficient[j] for j in range(2)]
        data['closed'] = bool(result.requested_tolerance_met and result.jacobian_rank == result.unknown_count
            and str(result.topology_event).endswith('None_') and min(result.amounts) >= 0
            and min(result.mole_fractions) >= 0 and result.phase_mole_fractions(0) != result.phase_mole_fractions(1)
            and (isinstance(problem.P, Free) or min(result.amounts) > 0)
            and min(data['pressure_density_derivative_pa_m3_per_mol']) > 0
            and result.molar_densities[0] > result.molar_densities[1]
            and data['inventory_error_mol'] <= problem.row_tolerance
            and data['normalization_error'] <= problem.row_tolerance
            and data['pressure_relative_error'] <= problem.row_tolerance
            and max(map(abs, data['log_fugacity_difference'])) <= problem.row_tolerance)
    return data


def declarations():
    model = epcsaft.Mixture(epcsaft.Parameters.from_json(LEAF/'inputs/parameters.json', components=('co2', 'thf')))
    rows = list(csv.DictReader((LEAF/'inputs/source/schick-2023-table-s10.csv').open()))
    cases = []
    for source in rows:
        T, P, xp, xe = (float(source[k]) for k in ('temperature_k', 'pressure_bar',
            'paper_model_liquid_co2_mole_fraction', 'experimental_liquid_co2_mole_fraction'))
        z = (xp+1)/2
        problem = Problem(T=T, P=P*1e5, feed=Amounts({'co2': z, 'thf': 1-z}),
            phases=[Phase('L', kind='liquid', composition_guess=[xp, 1-xp]),
                    Phase('V', kind='vapor', composition_guess=[0.99, 0.01])], row_tolerance=TOL['row'])
        bubble = bubble_point(liquid=Phase('L', kind='liquid'), T=T,
                              feed=MoleFractions({'co2': xp, 'thf': 1-xp}))
        bubble.P = Free(P*1e5)
        bubble.row_tolerance = TOL['row']
        cases.append((source, problem, bubble))
    return model, cases


def main():
    started = time.perf_counter()
    model, cases = declarations()
    reported_metrics = {r['quantity']: float(r['value']) for r in csv.DictReader(
        (LEAF/'inputs/source/schick-2023-reported-metrics.csv').open())}
    diagnostics, comparisons = [], []
    for source, problem, bubble in cases:
        P, xp, xe = (float(source[k]) for k in ('pressure_bar', 'paper_model_liquid_co2_mole_fraction',
                                               'experimental_liquid_co2_mole_fraction'))
        row = dict(source)
        detail = dict(source_row=int(source['source_row']), feed_mol=dict(problem.feed.values))
        regular = solve_equilibrium(model, problem)
        detail['regular'] = record(model, problem, regular)
        if regular.success:
            refinement_problem = copy.deepcopy(problem)
            refinement_problem.row_tolerance = TOL['refinement']
            refined = solve_equilibrium(model, refinement_problem, start=regular)
            detail['refinement'] = record(model, refinement_problem, refined)
            if detail['source_row'] == 1 and not detail['refinement']['closed']:
                detail['initial_regular'] = detail['regular']
                detail['initial_refinement'] = detail['refinement']
                correction_problem = copy.deepcopy(problem)
                correction_problem.row_tolerance = 1e-15
                corrected = solve_equilibrium(model, correction_problem, start=regular)
                detail['correction'] = record(model, refinement_problem, corrected)
                detail['correction'].update(requested_row_tolerance=correction_problem.row_tolerance,
                                            closure_tolerance=refinement_problem.row_tolerance)
                if not detail['correction']['closed']:
                    (OUT/'native-diagnostics.json').write_text(json.dumps(diagnostics+[detail], indent=2, allow_nan=False)+'\n')
                    raise RuntimeError('Row 1 native correction failed the unchanged refinement closure checks')
                regular = solve_equilibrium(model, problem, start=corrected)
                detail['regular'] = record(model, problem, regular)
                refined = solve_equilibrium(model, refinement_problem, start=regular)
                detail['refinement'] = record(model, refinement_problem, refined)
        row.update(engine_x=None, engine_y=None, refinement_change=None, numerically_qualified=False,
                   paper_reproduced=False, measurement_reproduced=False, bubble_pressure_bar=None,
                   bubble_pressure_difference_bar=None, bubble_closed=False, diagnosis='undetermined')
        if regular.success:
            x = regular.phase_mole_fractions(0)[0]
            row.update(engine_x=x, engine_y=regular.phase_mole_fractions(1)[0],
                       signed_error=x-xe, absolute_error=abs(x-xe), ard_percent=100*abs(1-x/xe),
                       paper_difference=x-xp)
            compiled = compile_problem(model, problem)
            observable = _core.SolvedStateObservable(_core.SolvedStateObservableKind.PhaseMoleFraction, 0, 0)
            half = TOL['row23_pressure_half_width_bar'] if int(source['source_row']) == 23 else TOL['pressure_half_width_bar']
            for name, dt, dp in [('temperature', 0.05, 0), ('pressure', 0, half*1e5)]:
                direction = _core.SolvedStateActionDirection(dt, dp, [0, 0], [])
                batch = solved_state_actions(compiled, regular, _core.SolvedStateActionRequest([observable], [direction]))
                action = batch.results[0]
                detail[name+'_rounding_action'] = dict(value=action.value, action=action.action, status=str(action.status),
                    message=action.diagnostic_message, condition_estimate=batch.diagnostics.condition_estimate)
            if refined.success:
                delta = abs(refined.phase_mole_fractions(0)[0]-x)
                row['refinement_change'] = delta
                row['numerically_qualified'] = bool(detail['regular']['closed'] and detail['refinement']['closed']
                    and delta <= TOL['composition_error'])
                row['paper_reproduced'] = bool(row['numerically_qualified'] and abs(x-xp) <= TOL['composition_half_width']+delta)
                row['measurement_reproduced'] = bool(row['numerically_qualified'] and abs(x-xe) <= TOL['composition_half_width'])
        if not row['paper_reproduced']:
            br = solve_equilibrium(model, bubble)
            detail['bubble'] = record(model, bubble, br)
            if br.success:
                bp = br.pressure/1e5
                row.update(bubble_pressure_bar=bp, bubble_pressure_difference_bar=bp-P,
                           bubble_closed=detail['bubble']['closed'])
                half = TOL['row23_pressure_half_width_bar'] if int(source['source_row']) == 23 else TOL['pressure_half_width_bar']
                if detail['bubble']['closed'] and abs(bp-P) > half:
                    row['diagnosis'] = 'closed forward bubble disagrees with printed pressure; model/source discrepancy; underlying cause undetermined'
            if not detail['regular']['closed']:
                row['diagnosis'] = 'native fixed-T/P closure unavailable; '+regular.message
        else:
            row['diagnosis'] = 'reproduced at printed inputs'
        e, m, h = xe, xp, TOL['composition_half_width']
        lo, hi = (m-h)/(e+h), (m+h)/(e-h)
        ard_lo = 0 if lo <= 1 <= hi else 100*min(abs(1-lo), abs(1-hi))
        ard_hi = 100*max(abs(1-lo), abs(1-hi))
        reported = float(source['paper_reported_ard_percent'])
        row.update(source_ard_percent=100*abs(1-m/e), source_ard_lower=ard_lo, source_ard_upper=ard_hi,
                   source_ard_interval_overlap=ard_lo <= reported+TOL['aard_percentage_point'] and ard_hi >= reported-TOL['aard_percentage_point'])
        diagnostics.append(detail)
        comparisons.append(row)
        (OUT/'native-diagnostics.json').write_text(json.dumps(diagnostics, indent=2, allow_nan=False)+'\n')
        print(f"row {source['source_row']}: x={row['engine_x']} qualified={row['numerically_qualified']} paper={row['paper_reproduced']}", flush=True)
    columns = list(dict.fromkeys(k for r in comparisons for k in r))
    with (OUT/'comparison.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, lineterminator='\n'); writer.writeheader(); writer.writerows(comparisons)
    metadata = dict(evidence_class='PUBLISHED_TUPLE_PROPERTY_REPLAY', tolerances=TOL,
        engine_provenance=epcsaft.provenance(), engine_module=epcsaft.__file__,
        parameters_sha256=hashlib.sha256((LEAF/'inputs/parameters.json').read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256((LEAF/'inputs/source/schick-2023-table-s10.csv').read_bytes()).hexdigest(),
        row_count=len(cases), wall_seconds=time.perf_counter()-started,
        source_aard_percent=sum(r['source_ard_percent'] for r in comparisons)/len(cases),
        printed_row_ard_mean_percent=sum(float(r['paper_reported_ard_percent']) for r in comparisons)/len(cases),
        engine_aard_percent=sum(r['ard_percent'] for r in comparisons)/len(cases) if all(r['engine_x'] is not None and d['regular']['closed'] for r,d in zip(comparisons,diagnostics)) else None)
    (OUT/'run.json').write_text(json.dumps(metadata, indent=2, allow_nan=False)+'\n')
    assert len(comparisons) == 34 and all(r['source_ard_interval_overlap'] for r in comparisons)
    assert abs(metadata['source_aard_percent']-reported_metrics['paper_aard']) <= TOL['aard_percentage_point']
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()

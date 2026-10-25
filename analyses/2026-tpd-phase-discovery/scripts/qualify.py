"""Installed-wheel #223 phase-count and tangent-plane falsifiers; native search only."""
import argparse
import copy
import itertools
import json
import math
import sys
from pathlib import Path
import time

import baseline as b
import costs
from epcsaft import trace

eq = b.eq
COSTS = []
FAILURES = []
UNRESOLVED = []


def require(condition, message):
    if not condition:
        FAILURES.append(message)
        raise RuntimeError(message)


def quantities(case, values):
    for name, value in values.items():
        b.record(case, name, value)


def report_rows(case, report):
    b.fields(case, report, ('completed', 'message', 'candidate_coordinates', 'standard_over_rt', 'plane_over_rt'))
    quantities(case+'/work', report.work)
    for i, phase in enumerate(report.present_phases):
        quantities(f'{case}/present-{i}', phase)
    for i, trial in enumerate(report.trials):
        name = f'{case}/trial-{i}'
        b.fields(name, trial, ('seed', 'evidence', 'message', 'raw_coordinates', 'raw_residuals', 'composition',
                             'split_coordinates', 'basis', 'row_transformation', 'donor', 'resolved_descent'))
        b.fields(name+'/phase', trial.phase, ('name', 'kind', 'support'))
        b.state(name+'/optimizer', trial.optimizer)
        quantities(name+'/quantities', trial.quantities)
        for j, split in enumerate(trial.split_checks):
            quantities(f'{name}/finite-split-{j}', split)


def trace_work(captured):
    spans = captured.spans
    root_evaluations = 0
    for span in spans:
        if span.name != 'evaluate':
            continue
        parent = span.parent
        while parent is not None:
            if spans[parent].name == 'solve_density':
                root_evaluations += 1
                break
            parent = spans[parent].parent
    return dict(driver_invocations=sum(e.name == 'hessian_calls' for e in captured.events),
                driver_spans=sum(s.name in ('solve_equilibrium', 'tangent_plane_trial') for s in spans),
                trial_invocations=sum(s.name == 'tangent_plane_trial' for s in spans),
                root_calls=sum(s.name == 'solve_density' for s in spans), root_evaluations=root_evaluations)


def timed(case, role, operation):
    first = None
    for sample in range(3):
        result = error = None
        start = time.perf_counter()
        with trace() as captured:
            try:
                result = operation()
            except Exception as caught:
                error = caught
        row = dict(case=case, role=role, sample=sample, seconds=time.perf_counter()-start,
                   trial_spans=sum(s.name == 'tangent_plane_trial' for s in captured.spans),
                   driver_finalizations=sum(e.name == 'hessian_calls' for e in captured.events),
                   hessian_calls=sum(sum(e.values) for e in captured.events if e.name == 'hessian_calls'),
                   iterations=sum(s.convergence.iterations for s in captured.spans if s.convergence and s.name in ('solve_equilibrium', 'tangent_plane_trial', 'discover_phases')),
                   evaluations=sum(s.convergence.evaluations for s in captured.spans if s.convergence and s.name in ('solve_equilibrium', 'tangent_plane_trial', 'discover_phases')),
                   trace_work=trace_work(captured))
        native_work = None if result is None else result.work if hasattr(result, 'work') else result.discovery.work if result.discovery is not None else None
        if native_work is not None:
            row['native_work'] = dict(native_work)
            for name in ('iterations', 'evaluations', 'hessian_calls'):
                if name in native_work:
                    row[name] = native_work[name]
        COSTS.append(row)
        if error is not None:
            b.record(f'{case}/{role}/sample-{sample}', 'exception_type', type(error).__name__)
            b.record(f'{case}/{role}/sample-{sample}', 'message', str(error))
            for i, span in enumerate(captured.spans):
                b.record(f'{case}/{role}/sample-{sample}/span-{i}', 'name', span.name)
                b.record(f'{case}/{role}/sample-{sample}/span-{i}', 'parent', span.parent)
            raise error
        if role == 'declared':
            require(row['trial_spans'] == 0, f'{case}: declared call implicitly ran a tangent-plane trial')
            require(result.discovery is None, f'{case}: declared call returned a discovery calculation')
        if sample == 0:
            first = result
            if hasattr(result, 'trials'):
                report_rows(case+'/'+role, result)
            else:
                b.state(case+'/'+role, result)
                b.record(case+'/'+role, 'message', result.message)
                if result.discovery is not None:
                    discovery_rows(case+'/'+role, result.discovery)
    return first


def discovery_rows(case, discovery):
    quantities(case+'/work', discovery.work)
    b.record(case, 'declaration_roles', discovery.declaration_roles)
    b.record(case, 'declaration_attempt_coordinates', discovery.declaration_attempt_coordinates)
    b.state(case+'/candidate', discovery.candidate)
    b.record(case+'/candidate', 'message', discovery.candidate.message)
    b.record(case, 'declaration_attempt_count', len(discovery.declaration_attempt_coordinates))
    for i, coordinates in enumerate(discovery.declaration_attempt_coordinates):
        b.record(f'{case}/attempt-{i}', 'coordinate_count', len(coordinates))
    for i, declaration in enumerate(discovery.declarations):
        b.state(f'{case}/declaration-{i}', declaration)
        b.record(f'{case}/declaration-{i}', 'message', declaration.message)
    for i, report in enumerate(discovery.rounds):
        report_rows(f'{case}/round-{i}', report)


def check_report(case, report, expect_descent=None, charges=None):
    require(report.completed, f'{case}: incomplete trial enumeration: {report.message}')
    witnesses = sum(t.resolved_descent for t in report.trials)
    if expect_descent is not None:
        require(bool(witnesses) == expect_descent, f'{case}: expected descent={expect_descent}, found {witnesses}')
    for phase in report.present_phases:
        require(abs(phase['present_identity_defect_over_rt']) <= 1e-10, f'{case}: present-phase residual identity exceeds 1e-10')
    for index, trial in enumerate(report.trials):
        q = trial.quantities
        if 'charge_residual' in q:
            require(abs(q['charge_residual']) <= 1e-12, f'{case}: reconstructed trial is not electroneutral')
            require(abs(q['normalization_residual']) <= 1e-12, f'{case}: reconstructed trial is not normalized')
        if charges and 'Vapor' in str(trial.phase.kind):
            require(all(charges[i] == 0 for i in trial.phase.support), f'{case}: vapor trial includes an ion')
        condensed_refusal = trial.message.endswith('charged trial outside the condensed domain of density-independent permittivity')
        if condensed_refusal:
            require(trial.evidence == 'unresolved' and 'tpd_over_rt' not in q and not trial.resolved_descent,
                    f'{case}/{trial.seed}: condensed-domain refusal has admitted evidence, TPD or resolved descent')
        elif trial.optimizer.success:
            require(abs(q['trial_identity_defect_over_rt']) <= 1e-10, f'{case}: stationary trial dual identity exceeds 1e-10')
        if trial.resolved_descent:
            require(q['tpd_over_rt'] < 0 and q['donor_slope_over_rt'] < -q['descent_error_estimate_over_rt'], f'{case}: unresolved slope admitted')
            require(q['split_delta_g_over_rt'] < -q['split_error_estimate_over_rt'], f'{case}: finite split lacks resolved descent')
        if trial.evidence in ('unresolved', 'indeterminate', 'negative_plane_without_positive_donor'):
            UNRESOLVED.append(dict(case=case, trial=index, seed=trial.seed, evidence=trial.evidence, message=trial.message))


def compare_phases(case, found, reference, composition_tolerance, energy_tolerance, density_tolerance=None):
    require(found.phase_count == reference.phase_count, f'{case}: wrong phase count {found.phase_count}')
    n = len(reference.mole_fractions)//reference.phase_count
    # At most three phases in these pinned cases: compare all label permutations, never search thermodynamics.
    orders = list(itertools.permutations(range(found.phase_count)))
    permutation = min(orders, key=lambda order: max(abs(found.mole_fractions[order[p]*n+i]-reference.mole_fractions[p*n+i]) for p in range(reference.phase_count) for i in range(n)))
    x_error = max(abs(found.mole_fractions[permutation[p]*n+i]-reference.mole_fractions[p*n+i]) for p in range(reference.phase_count) for i in range(n))
    amount_error = max(abs(found.amounts[permutation[p]]-reference.amounts[p]) for p in range(reference.phase_count))
    if case.startswith('case-7'):
        x_error = max(abs(found.mole_fractions[permutation[p]*n+i]/reference.mole_fractions[p*n+i]-1) if reference.mole_fractions[p*n+i] else abs(found.mole_fractions[permutation[p]*n+i]) for p in range(reference.phase_count) for i in range(n))
        amount_error = max(abs(found.amounts[permutation[p]]/reference.amounts[p]-1) for p in range(reference.phase_count) if reference.amounts[p])
    g_error = abs(found.objective-reference.objective)/max(1., abs(reference.objective))
    b.record(case+'/comparison', 'composition_max_abs', x_error)
    b.record(case+'/comparison', 'g_relative', g_error)
    b.record(case+'/comparison', 'amount_error', amount_error)
    require(x_error <= composition_tolerance and g_error <= energy_tolerance, f'{case}: phase comparison exceeds its issue tolerance ({x_error}, {g_error})')
    if case.startswith('case-7'):
        require(amount_error <= 1e-9, f'{case}: relative amount error exceeds 1e-9 ({amount_error})')
    if density_tolerance is not None:
        rho_error = max(abs(found.molar_densities[permutation[p]]/reference.molar_densities[p]-1) for p in range(reference.phase_count))
        b.record(case+'/comparison', 'density_max_relative', rho_error)
        require(rho_error <= density_tolerance, f'{case}: density comparison exceeds its issue tolerance')


def run_case(case, model, problem, expected_phases, expect_candidate_descent, *, test_only=False, lower_energy=False, charges=None):
    b.INPUTS[case] = dict(components=model.component_ids, **b.problem_input(problem))
    reference = timed(case, 'declared', lambda: eq.solve_equilibrium(model, problem))
    require(reference.success and reference.requested_tolerance_met, f'{case}: declared reference failed: {reference.message}')
    compiled = eq.compile_problem(model, problem)
    reference_report = timed(case, 'explicit-test', lambda: eq.stability_test(compiled, reference, eq.StabilityOptions()))
    check_report(case+'/declared-reference', reference_report, None if test_only else lower_energy, charges)
    if test_only:
        return
    test_reference = len(problem.phases) == 1 or case.startswith('case-4') or expect_candidate_descent is None
    candidate_problem = problem if test_reference else copy.copy(problem)
    if not test_reference:
        candidate_problem.phases = (eq.Phase('candidate', kind='liquid'),)
    candidate = reference if test_reference else eq.solve_equilibrium(model, candidate_problem)
    b.state(case+'/explicit-test-candidate', candidate)
    require(candidate.success, f'{case}: one-liquid candidate failed: {candidate.message}')
    candidate_compiled = eq.compile_problem(model, candidate_problem)
    candidate_report = reference_report if test_reference else timed(case, 'candidate-test', lambda: eq.stability_test(candidate_compiled, candidate, eq.StabilityOptions()))
    if not test_reference and expect_candidate_descent is not None:
        check_report(case, candidate_report, expect_candidate_descent, charges)
    request = copy.copy(problem)
    request.phases = eq.Discover()
    request.topology_declared = False
    found = timed(case, 'discover', lambda: eq.solve_equilibrium(model, request))
    require(found.success and found.requested_tolerance_met, f'{case}: discovery failed: {found.message}')
    if case.startswith('case-3'):
        phase_counts = [len(report.present_phases) for report in found.discovery.rounds]
        successful_additions = [result.phase_count for role, result in zip(found.discovery.declaration_roles, found.discovery.declarations) if role == 'addition' and result.success]
        b.record(case+'/addition-evidence', 'round_phase_counts', phase_counts)
        b.record(case+'/addition-evidence', 'successful_addition_phase_counts', successful_additions)
        require(phase_counts == [1, 2, 3] and 2 in successful_additions and 3 in successful_additions,
                f'{case}: native record does not show two successive additions ({phase_counts}, {successful_additions})')
    require(found.phase_count == expected_phases or lower_energy, f'{case}: expected {expected_phases} phases, found {found.phase_count}')
    if lower_energy:
        require(found.objective < reference.objective, f'{case}: discovery does not lower the declared three-liquid energy')
    else:
        compare_phases(case, found, reference, 1e-9, 1e-9 if case.startswith(('case-6', 'case-7')) else 1e-10)
    # The returned native declaration is the fixed topology consumed by later callers.
    final_compiled = eq.compile_problem(model, found.discovery.declaration)
    final = eq.stability_test(final_compiled, found, eq.StabilityOptions())
    report_rows(case+'/final-explicit-test', final)
    check_report(case+'/final', final, False, charges)
    redeclared = eq.solve_equilibrium(model, found.discovery.declaration, found)
    b.state(case+'/redeclared', redeclared)
    require(redeclared.success, f'{case}: returned native declaration cannot be solved')


def neutral():
    references = {case: (model, problem) for case, model, problem in b.declared_cases() if case.startswith(('case-1/', 'case-3/'))}
    for case, (model, problem) in references.items():
        run_case(case, model, problem, len(problem.phases), True if case.startswith('case-1/') else None)
    for case, species, feed in [('case-2/water-ethanol', ('water', 'ethanol'), (.5, .5)), ('case-2/dilute-isobutanol', ('water', 'isobutanol'), (.995, .005))]:
        model = b.Mixture(b.Parameters.from_json(b.KH_PACKET, components=species))
        problem = eq.Problem(T=293.15, P=1e5, row_tolerance=1e-10, feed=eq.Amounts(dict(zip(species, feed))), phases=[eq.Phase('L', kind='liquid')])
        run_case(case, model, problem, 1, False)
    packet = b.ROOT/'analyses/phase-stability-cutover/inputs/phase-stability-neutral-ternary.json'
    model = b.Mixture(b.Parameters.from_json(packet))
    b.INPUTS['case-4/parameters'] = dict(path=str(packet), sha256=b.digest(packet))
    for pressure in (1e7, 1e5):
        problem = eq.Problem(T=300., P=pressure, row_tolerance=1e-10, feed=eq.Amounts(dict.fromkeys(model.component_ids, 1./3)),
            phases=[eq.Phase(f'L{i}', kind='liquid', composition_guess=x) for i, x in enumerate(((.90, .06, .04), (.11, .83, .06), (.14, .12, .74)))])
        run_case(f'case-4/three-liquid-{pressure:g}-Pa', model, problem, 3, pressure == 1e5, lower_energy=pressure == 1e5)


def gibbs_scaling(reference_path):
    require(reference_path is not None, 'gibbs-scaling requires the retained one-mol neutral reference CSV')
    reference_path = reference_path.resolve()
    reference_metadata = json.loads(reference_path.with_suffix('.json').read_text())
    require(reference_metadata['failure'] is None, 'one-mol neutral qualification did not pass')
    unit = {r['quantity']: r for r in b.csv.DictReader(reference_path.open()) if r['case'] == 'case-1/two-liquid/declared'}
    for row in unit.values():
        if row['kind'] == 'float':
            require(b.struct.pack('>d', float(row['value'])).hex() == row['binary64'], 'one-mol reference lost binary64 identity')
    value = lambda key: float(unit[key]['value'])
    require(int(value('phase_count')) == 2 and abs(sum(value(f'amounts/{p}') for p in range(2))-1.) <= 1e-12, 'reference is not the one-mol two-liquid case')
    case, scale = 'case-1/100-mol', 100.
    model, problem = next((model, problem) for name, model, problem in b.declared_cases() if name == 'case-1/two-liquid')
    reference_input = reference_metadata['inputs']['case-1/two-liquid']
    require(reference_input['T'] == problem.T and reference_input['P'] == problem.P and reference_input['feed'] == problem.feed.values,
            'one-mol reference physical conditions differ')
    problem.feed = eq.Amounts({name: scale*amount for name, amount in problem.feed.values.items()})
    b.INPUTS[case] = dict(components=model.component_ids, **b.problem_input(problem),
        parameters=dict(path=str(b.KH_PACKET.relative_to(b.ROOT)), sha256=b.digest(b.KH_PACKET)),
        one_mol_reference=dict(csv=str(reference_path), sha256=b.digest(reference_path), metadata_sha256=b.digest(reference_path.with_suffix('.json'))),
        purpose='Extensive addition-error weighting and feed scaling; numerical verification, not an application prediction')
    request = copy.copy(problem)
    request.phases, request.topology_declared = eq.Discover(), False
    b.INPUTS[case]['phases'] = 'Discover(default trial families and budgets)'
    b.INPUTS[case]['scope_limit'] = 'Discovery only; the 100-mol cold declaration and returned-declaration reuse keep the split since #243 (engine/tests/equilibrium/test_feed_scale.py).'
    unit_order = sorted(range(2), key=lambda p: value(f'mole_fractions/{2*p}'))
    found = timed(case, 'discover', lambda: eq.solve_equilibrium(model, request))
    require(found.success and found.requested_tolerance_met and found.phase_count == 2, f'{case}: discovery failed: {found.message}')
    order = sorted(range(2), key=lambda p: found.mole_fractions[2*p])
    x_error = max(abs(found.mole_fractions[2*p+i]-value(f'mole_fractions/{2*q+i}')) for p, q in zip(order, unit_order) for i in range(2))
    amount_error = max(abs(found.amounts[p]-scale*value(f'amounts/{q}')) for p, q in zip(order, unit_order))
    density_error = max(abs(found.molar_densities[p]/value(f'molar_densities/{q}')-1) for p, q in zip(order, unit_order))
    g_error = abs(found.objective-scale*value('objective'))/max(1., abs(scale*value('objective')))
    inventory_error = max(abs(sum(found.amounts[p]*found.mole_fractions[2*p+i] for p in range(2))-problem.feed.values[name]) for i, name in enumerate(model.component_ids))
    quantities(case+'/discover/scaling', dict(feed_total_mol=scale, composition_max_abs=x_error,
        amount_scaling_max_abs_mol=amount_error, component_conservation_max_abs_mol=inventory_error,
        density_max_relative=density_error, g_relative=g_error))
    require(x_error <= 1e-9 and g_error <= 1e-10 and max(amount_error, inventory_error) <= 1e-12*scale,
            f'{case}: scaling exceeds the owner criteria ({x_error}, {g_error}, {amount_error}, {inventory_error})')
    for i, report in enumerate(found.discovery.rounds):
        check_report(f'{case}/round-{i}', report, False if i == len(found.discovery.rounds)-1 else None)
    first = found.discovery.rounds[0]
    require(len(first.present_phases) == 1, 'scaling margin requires the first one-phase addition')
    initial_index, initial = min(((i, result) for i, (role, result) in enumerate(zip(found.discovery.declaration_roles, found.discovery.declarations))
                                 if role == 'initial' and result.success), key=lambda item: item[1].objective)
    require(initial.phase_count == 1 and abs(initial.amounts[0]-scale) <= 1e-12*scale,
            'retained initial candidate is not the specified 100-mol single phase')
    q = first.present_phases[0]
    epsilon = q['association_energy_estimate_over_rt']+q['euler_defect_over_rt']+q['summation_estimate_over_rt']
    key = 'addition_gibbs_estimate_over_rt_mol'
    recorded = [(i, t.quantities[key]) for i, t in enumerate(first.trials) if key in t.quantities]
    require(recorded and math.isfinite(epsilon) and epsilon > 0, 'positive first-addition energy estimate was not retained')
    index, actual = recorded[0]
    expected = initial.amounts[0]*epsilon
    weighted_error, unweighted_error = abs(actual/expected-1), abs(epsilon/expected-1)
    quantities(case+'/addition-margin', dict(per_mol_energy_estimate_over_rt=epsilon,
        expected_extensive_estimate_over_rt_mol=expected, weighted_relative_error=weighted_error,
        old_unweighted_relative_error=unweighted_error))
    b.record(case+'/addition-margin', 'native_operand', f'{case}/discover/round-0/trial-{index}/quantities/{key}')
    b.record(case+'/addition-margin', 'weight_operand', f'{case}/discover/declaration-{initial_index}/amounts/0')
    require(weighted_error <= 1e-13 and unweighted_error > 1e-13, 'addition margin does not falsify missing amount weighting')


def condensed_case(case, model, problem, *, neutral_only=False):
    b.INPUTS.setdefault(case, {}).update(components=model.component_ids, **b.problem_input(problem))
    reference = timed(case, 'declared', lambda: eq.solve_equilibrium(model, problem))
    require(reference.success and reference.requested_tolerance_met, f'{case}: declared reference failed: {reference.message}')
    compiled = eq.compile_problem(model, problem)
    report = timed(case, 'explicit-test', lambda: eq.stability_test(compiled, reference, eq.StabilityOptions()))
    check_report(case+'/explicit-test', report)
    if case.startswith(('case-6/', 'case-7/')):
        charged = [t for t in report.trials if len(t.phase.support) == len(model.component_ids)]
        expected = {'positive-neutral-base', 'near-pure-0'}
        if len(model.component_ids) > 3:
            expected.update(('near-pure-1', 'near-pure-2'))
        require({t.seed for t in charged} == expected, f'{case}: charged trial seeds missing or changed')
        for trial in charged:
            unavailable = case.startswith('case-7/') and trial.seed == 'near-pure-0'
            if unavailable:
                require(trial.evidence == 'unresolved' and 'tpd_over_rt' not in trial.quantities
                        and trial.message.endswith('charged trial outside the condensed domain of density-independent permittivity'),
                        f'{case}/{trial.seed}: expected the condensed-domain refusal without TPD')
            else:
                require(trial.evidence != 'unresolved' and math.isfinite(trial.quantities.get('tpd_over_rt', math.nan)),
                        f'{case}/{trial.seed}: expected an admitted charged trial with finite TPD')
    request = copy.copy(problem)
    request.phases, request.topology_declared = eq.Discover(), False
    if isinstance(problem.P, (float, int)):
        found = timed(case, 'discover', lambda: eq.solve_equilibrium(model, request))
        require(found.discovery is not None, f'{case}: discovery retained no trial reports')
        for index, report in enumerate(found.discovery.rounds):
            check_report(f'{case}/discover/round-{index}', report)
        if not found.success:
            UNRESOLVED.append(dict(case=case, operation='discover', message=found.message))
    if neutral_only:
        species = tuple(name for name in model.component_ids if name not in ('sodium-cation', 'chloride-anion'))
        pairs = [eq.Phase('neutral-liquid', kind='liquid', support=species), eq.Phase('neutral-vapor', kind='vapor', support=species)]
        request.phases = eq.Discover(trials=pairs)
        # Reuse the existing Phase -> native-options transport, with no numerical calls.
        options = request._native(model).discovery
        report = timed(case, 'neutral-only-test', lambda: eq.stability_test(compiled, reference, options))
        check_report(case+'/neutral-only', report)
        neutral_indices = {model.component_ids.index(name) for name in species}
        require(report.trials and all(set(t.phase.support) <= neutral_indices for t in report.trials), f'{case}: replacement family included an ion or ran no trials')
        found = timed(case, 'neutral-only-discover', lambda: eq.solve_equilibrium(model, request))
        require(found.success and found.discovery is not None, f'{case}: neutral-only discovery failed: {found.message}')
        require(all(set(t.phase.support) <= neutral_indices for r in found.discovery.rounds for t in r.trials), f'{case}: discovery replaced the neutral-only family')
        quantities(case+'/neutral-only-comparison', dict(delta_g_over_rt_mol=found.objective-reference.objective,
            discovered_phase_count=found.phase_count, declared_phase_count=reference.phase_count))


def ionic():
    b.INPUTS['ionic/parameters'] = {str(p.relative_to(b.ROOT)): b.digest(p) for p in (b.KH_PACKET, b.MEA_PACKET)}
    for case, model, problem in b.declared_cases():
        if case.startswith(('case-6/', 'case-7/')):
            condensed_case(case, model, problem, neutral_only=case.startswith('case-6/'))
    species = ('water', 'sodium-cation', 'chloride-anion')
    model = b.Mixture(b.Parameters.from_json(b.KH_PACKET, components=species))
    problem = eq.Problem(T=298.15, P=1e5, row_tolerance=1e-10, feed=eq.Amounts(dict(zip(species, (.98, .01, .01)))), phases=[eq.Phase('L', kind='liquid')])
    condensed_case('case-6/aqueous-salt-one-liquid', model, problem)


def reacting_checks(case, model, result, gibbs, nu, *, expected=None):
    require(result.success and result.requested_tolerance_met, f'{case}: solve failed: {result.message}')
    largest = max(abs(value) for value in result.classified_residuals)
    quantities(case+'/checks', dict(classified_residual_max_abs=largest))
    require(largest <= 1e-10, f'{case}: strict rows exceed 1e-10 ({largest})')
    totals = [sum(result.amounts[p]*result.mole_fractions[p*3+i] for p in range(result.phase_count)) for i in range(3)]
    b.record(case+'/checks', 'component_amounts_mol', totals)
    if expected is not None:
        amount_error = max(abs(a-bb) for a, bb in zip(totals, expected))
        quantities(case+'/checks', dict(amount_max_abs_error_mol=amount_error))
        require(amount_error <= 1e-8, f'{case}: amount error exceeds 1e-8 mol ({amount_error})')
    for ion in (1, 2):
        # B=(W+nu_W Na, W+nu_W Cl); the feed has one mol of each ion.
        defect = totals[0]+nu[0]*totals[ion]-nu[0]
        quantities(case+'/checks', {f'moiety-{ion}_error_mol': defect})
        require(abs(defect) <= 1e-12, f'{case}: moiety error exceeds 1e-12 mol ({defect})')
    rt = b._core.GAS_CONSTANT_J_PER_MOL_K*result.temperature
    for p in range(result.phase_count):
        x = result.mole_fractions[p*3:(p+1)*3]
        norm = sum(x)-1
        charge = x[1]-x[2]
        quantities(f'{case}/phase-{p}/checks', dict(normalization_residual=norm, charge_residual=charge))
        require(abs(norm) <= 1e-12 and abs(charge) <= 1e-12, f'{case}: phase normalization/charge exceeds 1e-12')
        missing = [model.component_ids[i] for i, coefficient in enumerate(nu) if coefficient != 0 and x[i] <= 0]
        if missing or result.amounts[p] <= 0:
            # The shared classified chemical rows above still apply; log-affinity does not.
            b.record(f'{case}/phase-{p}/independent-affinity', 'availability', 'structural_zero' if missing else 'zero_amount_phase')
            b.record(f'{case}/phase-{p}/independent-affinity', 'unavailable_participants', missing)
            b.record(f'{case}/phase-{p}/independent-affinity', 'affinity_over_rt', None)
            continue
        physical = model.state(result.temperature, rho=result.molar_densities[p], x=x)
        mu = [math.log(physical.molar_density*x[i])+physical.residual_chemical_potential_over_rt[i]+gibbs[i]/rt for i in range(3)]
        affinity = sum(a*bb for a, bb in zip(nu, mu))
        b.record(f'{case}/phase-{p}/independent-affinity', 'availability', 'available')
        quantities(f'{case}/phase-{p}/independent-affinity', dict(affinity_over_rt=affinity, pressure_relative_residual=abs(physical.pressure/result.pressure-1)))
        require(abs(affinity) <= 1e-10, f'{case}: independent affinity exceeds 1e-10 ({affinity})')


def reacting(reference_path=None):
    species = ('water', 'sodium-cation', 'chloride-anion')
    path = b.ROOT/'analyses/2026-dielectric-alternatives/inputs/selections/density-temperature.json'
    model = b.Mixture(b.Parameters.from_json(path, components=species))
    x = (.98, .01, .01)
    rt = b._core.GAS_CONSTANT_J_PER_MOL_K*298.15
    frozen_k = {}
    provenance = dict(path=str(path.relative_to(b.ROOT)), sha256=b.digest(path), construction='equilibrium_checks.cpp manufactured charged reference; one Liquid x* state, frozen g')
    if reference_path is None:
        central = model.state(298.15, P=1e5, x=x, phase='liquid')
        gibbs = [-rt*(math.log(central.molar_density*x[i])+central.residual_chemical_potential_over_rt[i]) for i in range(3)]
        quantities('reacting/reference', dict(T_K=298.15, P_Pa=1e5, rho_mol_m3=central.molar_density, x=x, g_j_per_mol=gibbs, R_j_per_mol_K=b._core.GAS_CONSTANT_J_PER_MOL_K))
    else:
        prior = json.loads(reference_path.with_suffix('.json').read_text())
        require(prior['inputs']['reacting/reference']['sha256'] == b.digest(path), 'frozen reacting parameter selection differs')
        provenance['reference_wheel_sha256'] = prior['identity']['wheel_sha256']
        values = {}
        for row in b.csv.DictReader(reference_path.open()):
            value = float(row['value'])
            require(b.struct.pack('>d', value).hex() == row['binary64'], 'frozen reacting binary64 operand differs')
            values[row['case'], row['quantity']] = value
        require(values['reacting/reference', 'T_K'] == 298.15 and values['reacting/reference', 'P_Pa'] == 1e5 and
                values['reacting/reference', 'R_j_per_mol_K'] == b._core.GAS_CONSTANT_J_PER_MOL_K and
                tuple(values['reacting/reference', f'x/{i}'] for i in range(3)) == x, 'frozen reacting reference conditions differ')
        gibbs = [values['reacting/reference', f'g_j_per_mol/{i}'] for i in range(3)]
        if ('reacting/fixed-TP/inputs', 'equilibrium_constant') in values:
            frozen_k[(2., -1., -1.)] = values['reacting/fixed-TP/inputs', 'equilibrium_constant']
        provenance.update(operands_csv=str(reference_path.resolve()), operands_sha256=b.digest(reference_path),
                          method='Continue with the exact already-frozen g/K operands after the constructor mapping correction; no second reference-state evaluation')
    b.INPUTS['reacting/reference'] = provenance
    for pinned, nu, expected in [(False, (2., -1., -1.), (1.96, .02, .02)), (True, (1., -1., -1.), (98/99, 1/99, 1/99))]:
        case = 'reacting/pinned-free-pressure' if pinned else 'reacting/fixed-TP'
        k = frozen_k[nu] if nu in frozen_k else math.exp(-sum(a*bb for a, bb in zip(nu, gibbs))/rt)
        quantities(case+'/inputs', dict(stoichiometry=nu, equilibrium_constant=k, expected_amounts_mol=expected))
        problem = eq.Problem(T=298.15, P=eq.Free(1e5) if pinned else 1e5, row_tolerance=1e-10,
            feed=eq.Amounts(dict(zip(species, (0., 1., 1.)))), standard_gibbs_j_per_mol=gibbs,
            phases=[eq.Phase('L', kind='liquid', composition_guess=x, amount=eq.Pinned(50/99) if pinned else eq.Free())],
            reactions=[eq.Reaction(dict(zip(species, nu)), name='manufactured-water-salt', equilibrium_constant=k)])
        b.INPUTS[case] = dict(components=species, **b.problem_input(problem), reaction_operands_csv_case=case+'/inputs', standard_operands_csv_case='reacting/reference')
        reference = timed(case, 'declared', lambda: eq.solve_equilibrium(model, problem))
        reacting_checks(case+'/declared', model, reference, gibbs, nu, expected=expected)
        compiled = eq.compile_problem(model, problem)
        report = timed(case, 'explicit-test', lambda: eq.stability_test(compiled, reference, eq.StabilityOptions()))
        check_report(case+'/declared-test', report, charges=(0, 1, -1))
        require(any(any(i in (1, 2) for i in trial.phase.support) for trial in report.trials), f'{case}: no charged trial was performed')
        standards_error = max(abs(a-bb/rt) for a, bb in zip(report.standard_over_rt, gibbs))
        quantities(case+'/snapshot', dict(standard_max_abs_error_over_rt=standards_error, solved_T_K=reference.temperature, solved_P_Pa=reference.pressure))
        require(len(report.standard_over_rt) == 3 and standards_error <= 1e-10, f'{case}: snapshot standards differ from frozen g')
        if pinned:
            continue
        request = copy.copy(problem)
        request.phases, request.topology_declared = eq.Discover(), False
        found = timed(case, 'discover', lambda: eq.solve_equilibrium(model, request))
        reacting_checks(case+'/discovered', model, found, gibbs, nu)
        delta = found.objective-reference.objective
        quantities(case+'/discovery-comparison', dict(delta_g_over_rt_mol=delta, observed_phase_count=found.phase_count))
        require(delta <= 1e-9*max(1., abs(reference.objective)), f'{case}: discovery increased G/RT ({delta})')
        returned = found.discovery.declaration
        final_compiled = eq.compile_problem(model, returned)
        final = timed(case, 'returned-test', lambda: eq.stability_test(final_compiled, found, eq.StabilityOptions()))
        check_report(case+'/returned-test', final, charges=(0, 1, -1))
        redeclared = timed(case, 'returned-declared', lambda: eq.solve_equilibrium(model, returned, found))
        reacting_checks(case+'/redeclared', model, redeclared, gibbs, nu)
        g_error = abs(redeclared.objective-found.objective)/max(1., abs(found.objective))
        quantities(case+'/returned-comparison', dict(g_relative=g_error))
        require(g_error <= 1e-9, f'{case}: returned declaration G differs by {g_error}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('neutral', 'ionic', 'reacting', 'permittivity', 'gibbs-scaling'), required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reacting-reference', type=Path, help='Reuse exact frozen g/K CSV operands with unchanged physical selection and conditions')
    parser.add_argument('--neutral-reference', type=Path, help='Retained one-mol case-1 neutral CSV for the exact 100-mol scaling check')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if (args.output/f'{args.stage}.csv').exists() or (args.output/f'{args.stage}.json').exists():
        raise RuntimeError('output already contains retained rows')
    failure = None
    distribution = b.metadata.distribution('epcsaft')
    direct = json.loads(distribution.read_text('direct_url.json'))
    wheel = Path(direct['url'].removeprefix('file://'))
    identity = dict(wheel=str(wheel), wheel_sha256=b.digest(wheel), script_sha256=b.digest(__file__), baseline_script_sha256=b.digest(b.__file__),
                    source_mtimes={str(p.relative_to(b.ROOT)): p.stat().st_mtime_ns for p in (b.ROOT/'engine/native/equilibrium/gibbs_minimizer.cpp', b.ROOT/'engine/native/equilibrium/request_compile.cpp')},
                    thread_environment={name: b.os.environ.get(name) for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS')})
    if args.stage == 'permittivity':
        identity['permittivity_script_sha256'] = b.digest(Path(__file__).with_name('permittivity.py'))
    try:
        if args.stage == 'gibbs-scaling':
            gibbs_scaling(args.neutral_reference)
        elif args.stage == 'reacting':
            reacting(args.reacting_reference)
        elif args.stage == 'permittivity':
            import permittivity
            permittivity.qualify(sys.modules[__name__])
        else:
            globals()[args.stage]()
    except Exception as error:
        failure = f'{type(error).__name__}: {error}'
        raise
    finally:
        b.csv_write(args.output/f'{args.stage}.csv', b.ROWS)
        costs.write(args.output/f'{args.stage}.json', dict(identity=identity, inputs=b.INPUTS, costs=COSTS, failure=failure, failures=FAILURES, unresolved=UNRESOLVED,
            limits=['numerical verification of the supplied model and finite trial set', 'estimated uncertainty is not a rigorous bound',
                'Condensed-domain refusal applies per trial; unsuccessful discovery remains unresolved.',
                'no physical validation or global stability certificate']))


if __name__ == '__main__':
    main()

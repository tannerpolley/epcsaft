"""#223 acceptance 10: explicit declared calculations, with binary64 replay operands."""
import argparse
import ast
import csv
import hashlib
import gzip
import importlib.metadata as metadata
import importlib.util
import json
import math
import os
from pathlib import Path
import struct
import subprocess
import sys
import time

import epcsaft
from epcsaft import Mixture, Parameters, _core, equilibrium as eq, regression as r

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parents[1]
KH = ROOT / 'analyses/2026-khudaida-idaes-settler'
KH_PACKET = KH / 'inputs/data/khudaida-2026-figure-2-electrolyte-lle/7/parameters/parameters.json'
MEA_PACKET = ROOT / 'data/packets/mea-co2-h2o-nine-species-estimation/3/parameters/parameters.json'
GROSS_PACKET = ROOT / 'data/packets/gross-2001-methane-ethane/3/parameters/parameters.json'
FIT_PACKET = ROOT / 'data/packets/figiel-2025-reference-electrolytes/4/parameters/parameters.json'
STATE_FIELDS = ('success', 'solver_status', 'temperature', 'pressure', 'amounts', 'mole_fractions',
                'molar_densities', 'packing_fractions', 'multipliers', 'constraint_multipliers',
                'residuals', 'classified_residuals', 'coordinates', 'phase_count',
                'iterations', 'evaluations', 'hessian_calls', 'objective', 'reduced_hessian_minimum',
                'unknown_count', 'jacobian_rank', 'generic_nullity', 'requested_tolerance_met',
                'topology_event', 'topology_phase', 'topology_pair')
ROWS = []
INPUTS = {}
COSTS = []


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def record(case, quantity, value):
    if isinstance(value, (list, tuple)):
        for i, item in enumerate(value):
            record(case, f'{quantity}/{i}', item)
    else:
        kind = 'none' if value is None else 'bool' if isinstance(value, bool) else 'int' if isinstance(value, int) else 'float' if isinstance(value, float) else 'text'
        ROWS.append(dict(case=case, quantity=quantity, kind=kind,
                         value=repr(value) if kind != 'text' else str(value),
                         binary64=struct.pack('>d', value).hex() if kind == 'float' else ''))


def fields(case, obj, names):
    for name in names:
        value = getattr(obj, name)
        # Records keep the one-entry T and P lists of earlier wheels (#609).
        record(case, name, ([] if value is None else [value]) if name in ('temperature', 'pressure') else value)


def state(case, result):
    fields(case, result, STATE_FIELDS)
    fields(case + '/convergence', result.convergence,
           ('residual', 'floor', 'amplification', 'contraction', 'iterations', 'evaluations', 'stop'))
    record(case + '/convergence', 'criterion', result.convergence.criterion and result.convergence.criterion.value)
    for i, row in enumerate(result.rows):
        fields(f'{case}/row/{i}', row, ('eqid', 'kind', 'phase', 'index'))


def problem_input(problem):
    return dict(T=problem.T if isinstance(problem.T, float) else vars_intensive(problem.T),
                P=problem.P if isinstance(problem.P, float) else vars_intensive(problem.P),
                feed=dict(problem.feed.values), row_tolerance=problem.row_tolerance,
                phases=[dict(name=p.name, kind=p.kind, support=p.support, amount=vars_intensive(p.amount),
                             composition_guess=p.composition_guess, packing_guess=p.packing_guess) for p in problem.phases])


def vars_intensive(value):
    return dict(kind=type(value).__name__, value=value.value if isinstance(value, eq.Pinned) else value.guess)


def solve(case, model, problem):
    INPUTS[case] = dict(components=model.component_ids, **problem_input(problem))
    start = time.perf_counter()
    result = eq.solve_equilibrium(model, problem)
    COSTS.append(dict(case=case, seconds=time.perf_counter()-start, iterations=result.iterations,
                      evaluations=result.evaluations, hessian_calls=result.hessian_calls))
    state(case, result)
    if not result.success or result.requested_tolerance_met is not True:
        raise RuntimeError(f'{case}: {result.message}')
    compiled = eq.compile_problem(model, problem)
    observables = [eq.SolvedStateObservable(eq.SolvedStateObservableKind.PhaseComponentAmount, p, i)
                   for p in range(result.phase_count) for i in range(len(model.component_ids))]
    direction = eq.SolvedStateActionDirection(1., 0., [0.] * len(model.component_ids), [])
    for order in (1, 2):
        batch = eq.solved_state_actions(compiled, result, eq.SolvedStateActionRequest(observables, [direction]*order))
        for i, item in enumerate(batch.results):
            fields(f'{case}/T-action-{order}/{i}', item, ('value', 'action', 'status', 'order', 'unit'))
        fields(f'{case}/T-action-{order}/diagnostics', batch.diagnostics,
               ('equation_count', 'state_dimension', 'maximum_energy_order', 'state_factorizations', 'linear_solves',
                'rank', 'condition_estimate', 'maximum_linear_residual', 'central_residual', 'central_floor'))
    return result


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def declared_cases():
    model = Mixture(Parameters.from_json(KH_PACKET, components=('water', 'isobutanol')))
    yield ('case-1/two-liquid', model, eq.Problem(T=293.15, P=1e5, row_tolerance=1e-10,
          feed=eq.Amounts({'water': .6, 'isobutanol': .4}),
          phases=[eq.Phase('Lw', kind='liquid', composition_guess=(.999, .001)),
                  eq.Phase('Lo', kind='liquid', composition_guess=(.3, .7))]))
    model = Mixture(Parameters.from_json(KH_PACKET, components=('water', 'ethanol', 'isobutanol')))
    yield ('case-3/VLLE', model, eq.Problem(T=364.2, P=1e5, row_tolerance=1e-10,
          feed=eq.Amounts({'water': .8, 'ethanol': .02, 'isobutanol': .18}),
          phases=[eq.Phase('Lw', kind='liquid', composition_guess=(.98, .01, .01)),
                  eq.Phase('Lo', kind='liquid', composition_guess=(.5, .03, .47)),
                  eq.Phase('V', kind='vapor', composition_guess=(.713829, .039282, .246889))]))
    fixed = load('run_fixed_settler', KH / 'scripts/run_fixed_settler.py')
    model = Mixture(Parameters.from_json(KH_PACKET, components=fixed.COMPONENTS))
    yield ('case-6/figure-2-point-7', model, fixed.engine_problem(fixed.feed_rows()[7]))
    tree = ast.parse((ROOT / 'engine/tests/equilibrium/test_trace_reliability.py').read_text())
    constants = {node.targets[0].id: ast.literal_eval(node.value) for node in tree.body
                 if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name)
                 and node.targets[0].id in ('_MEA', '_MEA_REACTIONS')}
    species = constants['_MEA']
    ln_k = (-605.348, -587.8, -948.2312, -25.1014, -298.6415)
    reactions = [eq.Reaction(nu, correlation=eq.ReactionLogPolynomial(a, 0., 0., 0., 318.15, 250., 450., eq.EOS_STANDARD_STATE_ID))
                 for nu, a in zip(constants['_MEA_REACTIONS'], ln_k)]
    INPUTS['case-7/reactions'] = dict(species=species, stoichiometry=constants['_MEA_REACTIONS'], ln_k=ln_k,
                                   T_ref=318.15, T_min=250., T_max=450., standard_state=eq.EOS_STANDARD_STATE_ID)
    mapping = json.loads(MEA_PACKET.read_text())
    # Packet 3 has no CO2/MEA cross edge; declared to keep that model (#569), not a physical-absence claim.
    mapping['topology']['no_cross_association'] = [{
        'component_id_a': 'carbon-dioxide', 'component_id_b': 'monoethanolamine',
        'provenance': {'source_id': mapping['document_id'], 'locator': 'packet 3: no CO2/MEA edge'}}]
    model = Mixture(Parameters.from_mapping(mapping, components=species))
    feed = eq.Amounts({'carbon-dioxide': .25557798094241313*9.7292821, 'monoethanolamine': 9.7292821, 'water': 76.9751925})
    phases = [eq.Phase('L', kind='liquid'), eq.Phase('V', kind='vapor', support=species[:3], amount=eq.Pinned(0.))]
    yield ('case-7/bubble', model, eq.Problem(T=318.15, P=eq.Free(7398.16), phases=phases, feed=feed, reactions=reactions, row_tolerance=1e-10))
    phases = [eq.Phase('L', kind='liquid'), eq.Phase('V', kind='vapor', support=species[:3])]
    yield ('case-7/flash-5918.5286-Pa', model, eq.Problem(T=318.15, P=5918.5286, phases=phases, feed=feed, reactions=reactions, row_tolerance=1e-10))
    yield ('case-7/liquid-8877.7929-Pa', model, eq.Problem(T=318.15, P=8877.7929, phases=[eq.Phase('L', kind='liquid')], feed=feed, reactions=reactions, row_tolerance=1e-10))


def declared():
    for case in declared_cases():
        solve(*case)


def trace_rows(captured):
    global ROWS
    values, ROWS = ROWS, []
    try:
        for i, span in enumerate(captured.spans):
            record(f'span-{i}', 'name', span.name)
            record(f'span-{i}', 'parent', span.parent)
            if span.convergence is not None:
                fields(f'span-{i}/convergence', span.convergence, ('iterations', 'evaluations', 'residual', 'stop'))
        for i, event in enumerate(captured.events):
            fields(f'event-{i}', event, ('name', 'span', 'iteration', 'values'))
        return ROWS
    finally:
        ROWS = values


def csv_write(path, rows):
    with (gzip.open(path, 'wt', newline='') if path.suffix == '.gz' else path.open('w', newline='')) as stream:
        writer = csv.DictWriter(stream, fieldnames=('case', 'quantity', 'kind', 'value', 'binary64'))
        writer.writeheader()
        writer.writerows(rows)


def continuation():
    model = Mixture(Parameters.from_json(GROSS_PACKET, components=('methane', 'ethane')))
    def bubble(T):
        return eq.Problem(T=T, P=eq.Free(1e5), phases=[eq.Phase('L', kind='liquid'), eq.Phase('V', kind='vapor', amount=eq.Pinned(0.))], feed=eq.Amounts({'methane': .5, 'ethane': .5}))
    for step, steps in ((.5, 2), (.25, 4)):
        case = f'continuation/bubble-180-to-181-K/step-{step}'
        INPUTS[case] = dict(start=problem_input(bubble(180.)), target=problem_input(bubble(181.)),
                           initial_step=step, minimum_step=step, max_steps=steps, max_retries_per_step=1)
        start = time.perf_counter()
        result = eq.continue_equilibrium(model, bubble(180.), bubble(181.), eq.Continuation(initial_step=step, minimum_step=step, max_steps=steps, max_retries_per_step=1))
        COSTS.append(dict(case=case, seconds=time.perf_counter()-start))
        fields(case, result, ('success', 'stop', 'path_values', 'reached_path', 'attempted_path', 'event_phase', 'event_pair', 'attempts', 'accepted_steps', 'retries', 'cause_stop'))
        for i, item in enumerate(result.states):
            state(f'{case}/state-{i}', item)
        if result.last_result is not None:
            state(case+'/last-result', result.last_result)


def fitting():
    p = Parameters.from_json(FIT_PACKET, components=('water', 'sodium-cation', 'chloride-anion'))
    cs = [r.coordinate(p, [('k_ij', pair, 1)], origin=0, scale=.2, bounds=(-1, 1)) for pair in (('sodium-cation', 'water'), ('sodium-cation', 'chloride-anion'))]
    cs += [r.coordinate(p, [('ion_suppression_coefficient', 'sodium-cation', 1), ('ion_suppression_coefficient', 'chloride-anion', 1)], origin=0, scale=2, bounds=(3, 9))]
    energy = r.electrolyte(p, 'solvation_gibbs_energy', 298.15, 1e5, {'water': 1}, ions=('sodium-cation', 'chloride-anion'), stoichiometry=(1, 1), molality=0, species='sodium-cation', observed=-381930, form='difference', scale=1)
    holds = [r.hydration_hold(energy, bracket=(1, 8), energy_tolerance=1e-6, diameter_tolerance=1e-10, maximum_iterations=50)]
    physical_inputs = [-.25, .9, 6.01]
    physical_cs = []
    for coordinate in cs:
        unit = r.Coordinate()
        unit.targets = coordinate.targets
        unit.origin, unit.scale = 0, 1
        physical_cs.append(unit)
    declared = r.evaluate_coordinates(p, physical_cs, physical_inputs, [])
    native = r.CoordinateSet()
    native.coordinates = cs
    truth = native.start(declared.parameters)
    templates = [r.electrolyte(p, kind, 298.15, 1e5, {'water': 1}, ions=('sodium-cation', 'chloride-anion'), stoichiometry=(1, 1), molality=m, observed=1, form='log_ratio') for m in (.1, 1, 4.5) for kind in ('mean_ionic_activity_coefficient', 'osmotic_coefficient')]
    data = _core.observations_at(p._native, cs, truth, holds, templates, [0]*6)
    INPUTS['fitting/219-held-NaCl'] = dict(source='engine/tests/regression/test_fit.py::test_native_held_electrolyte_recovery_checks_and_columns', physical_inputs=physical_inputs, truth=truth, observations=[dict(row=i, kind=ob.kind.name, molality=ob.molality, observed=ob.observed) for i, ob in enumerate(data)], weights=[1/(.02*math.log(10))**2]*6, controls=dict(function_tolerance=1e-15, gradient_tolerance=1e-14, parameter_tolerance=1e-14))
    observation = _core.evaluate_observations(p._native, cs, truth, holds, data)
    fields('fitting/219-observations', observation, ('predictions', 'gradient', 'optimizer_gradient', 'statuses'))
    controls = r.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0,
                             function_tolerance=1e-15, gradient_tolerance=1e-14, parameter_tolerance=1e-14)
    start = time.perf_counter()
    result = r.fit(p, cs, data, holds=holds, controls=controls, weights=INPUTS['fitting/219-held-NaCl']['weights'])
    COSTS.append(dict(case='fitting/219-held-NaCl', seconds=time.perf_counter()-start))
    fields('fitting/219-held-NaCl', result, ('optimizer', 'status', 'usable', 'physical', 'predictions', 'residuals', 'weighted_residuals', 'physical_jacobian', 'optimizer_jacobian', 'statuses', 'iterations', 'residual_evaluations', 'jacobian_evaluations', 'trial_failures', 'iteration_costs', 'initial_cost', 'final_cost'))
    fields('fitting/219-held-NaCl/covariance', result.covariance, ('rank', 'singular_values', 'optimizer', 'physical', 'variance_factor', 'projected_gradient_norm', 'status'))
    if result.status != 'converged':
        raise RuntimeError(result.message)


def asl():
    sys.path.insert(0, str(KH / 'scripts'))
    live = load('khudaida_live', KH / 'scripts/run_live_settler.py')
    reference = live.Reference()
    INPUTS['asl/anchor'] = dict(F0=live.F0.tolist(), continuation=dict(initial_step=.25, minimum_step=.25, max_steps=4, max_retries_per_step=1), record_sha256=hashlib.sha256(reference.record.encode()).hexdigest(), library_sha256=digest(reference.library))
    state('asl/anchor/declared', reference.anchor)
    for label, F in [('anchor', live.F0), ('start-0.9', live.F0*(1., 1., .9, .9)), ('start-1.1', live.F0*(1., 1., 1.1, 1.1))]:
        INPUTS['asl/'+label] = dict(feed=F.tolist())
        start = time.perf_counter()
        public, diagnostics = reference.public(F)
        native = reference.asl(F)
        COSTS.append(dict(case='asl/'+label, seconds=time.perf_counter()-start))
        for surface, values in (('public', public), ('ASL', native)):
            for name, value in zip(('amounts', 'gradient', 'hessian'), values):
                record('asl/'+label+'/'+surface, name, value.tolist())
        for name, value in diagnostics.items():
            record('asl/'+label+'/diagnostics', name, value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--section', choices=('declared', 'continuation', 'fitting', 'asl'), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    target = args.output / f'{args.section}.csv'
    if target.exists() or (args.output/f'{args.section}.json').exists():
        raise RuntimeError(f'refusing to overwrite {target}')
    sys.path.insert(0, str(ROOT / 'tools'))
    from engine_identity import loaded_identity
    prepared_environment = loaded_identity()
    identity = dict(**prepared_environment, engine_import=epcsaft.__file__,
                    source_head=prepared_environment['source_revision'],
                    thread_environment={name: os.environ.get(name) for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS')},
                    dependencies={name: prepared_environment['version'] if name == 'epcsaft' else metadata.version(name)
                                  for name in ('epcsaft', 'numpy', 'scipy', 'pyomo', 'idaes-pse')},
                    physical_inputs={str(p.relative_to(ROOT)): digest(p) for p in (KH_PACKET, MEA_PACKET, GROSS_PACKET, FIT_PACKET)},
                    script_sha256=digest(__file__), source_mtimes={str(p.relative_to(ROOT)): p.stat().st_mtime_ns for p in (ROOT/'engine/native/equilibrium/gibbs_minimizer.cpp', ROOT/'engine/native/equilibrium/request_compile.cpp')})
    failure = None
    captured = None
    try:
        with epcsaft.trace() as captured:
            globals()[args.section]()
    except Exception as error:
        failure = f'{type(error).__name__}: {error}'
        raise
    finally:
        csv_write(target, ROWS)
        trace_path = args.output / f'{args.section}-trace.csv.gz'
        if captured is not None:
            csv_write(trace_path, trace_rows(captured))
        trace_counts = None if captured is None else dict(tangent_plane_trials=sum(s.name == 'tangent_plane_trial' for s in captured.spans),
            driver_finalizations=sum(e.name == 'hessian_calls' for e in captured.events),
            hessian_calls=sum(sum(e.values) for e in captured.events if e.name == 'hessian_calls'))
        (args.output / f'{args.section}.json').write_text(json.dumps(dict(identity=identity, inputs=INPUTS, costs=COSTS, failure=failure, trace_counts=trace_counts,
            limits=['declared topology numerical replay; no physical validation']), indent=2)+'\n')


def declarations():
    return list(declared_cases())


if __name__ == '__main__':
    main()

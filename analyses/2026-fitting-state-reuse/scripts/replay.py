"""Replay retained MEA physical inputs through the installed Engine wheel."""
import csv
import hashlib
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import epcsaft
from epcsaft import equilibrium as eq, regression as reg

ROOT = Path(__file__).resolve().parents[3]
STUDY = Path(__file__).resolve().parents[1]
INPUT = STUDY / 'inputs/mea-initial.json'


def declaration(value):
    if isinstance(value, list):
        return [declaration(item) for item in value]
    if not isinstance(value, dict) or 'type' not in value:
        return value
    kind, fields = value['type'], value['fields']
    fields = {key: declaration(item) for key, item in fields.items()}
    if kind in ('Amounts', 'MoleFractions', 'MassFractions'):
        return getattr(eq, kind)(fields['values'])
    if kind == 'Reaction':
        if fields['correlation'] is not None:
            fields['correlation'] = eq.ReactionLogPolynomial(**fields['correlation'])
        if fields['reference'] is not None:
            ref = fields['reference']
            fields['reference'] = eq.ReactionReference(
                getattr(eq.ReferenceSourceBasis, ref['source_basis']), ref['reference_pressure_pa'])
    return getattr(eq, kind)(**fields)


def declared(data):
    """The retained parameter mapping with its CO2/MEA omission declared (#569): no cross edge, not a
    physical-absence claim."""
    data['parameters']['topology']['no_cross_association'] = [{
        'component_id_a': 'carbon-dioxide', 'component_id_b': 'monoethanolamine',
        'provenance': {'source_id': 'mea-initial', 'locator': 'inputs/mea-initial.json: no CO2/MEA edge'}}]
    return data['parameters']


def main():
    mode, output_name = sys.argv[1:]
    if mode not in ('initial', 'short', 'affine'):
        raise ValueError('mode must be initial, short or affine')
    data = json.loads(INPUT.read_text())
    parameters = epcsaft.Parameters.from_mapping(declared(data))
    coordinates = []
    for record in data['coordinates']:
        coordinate = reg.Coordinate()
        targets = []
        for row in record['targets']:
            target = reg.CoordinateTarget()
            target.family = getattr(epcsaft._core.ParameterFamily, row['family'])
            target.indices, target.coefficient = row['indices'], row['coefficient']
            targets.append(target)
        coordinate.targets = targets
        for key in ('origin', 'scale', 'lower', 'upper'): setattr(coordinate, key, record[key])
        coordinates.append(coordinate)
    if mode == 'affine':
        coordinates[0].origin = 1e6  # resolve the finite affine round trip at the initial physical point
    states = data['states']
    if mode in ('short', 'affine'):
        states = [state for state in states if len(state['targets']) > 1][:3] + states[:3]
    observations, weights, identities = [], [], []
    kinds = {'PartialPressure': 'partial_pressure', 'PhaseMoleFraction': 'mole_fraction'}
    for state in states:
        problem = declaration(state['problem'])
        for row in state['targets']:
            form = 'log_ratio' if row['form'] == 'LogRatio' else 'difference'
            observations.append(reg.observation(parameters, kinds[row['kind']], problem,
                observed=row['observed'], form=form, scale=None if form == 'log_ratio' else row['scale'],
                mass_basis=row['mass_basis'], phase=problem.phases[row['phase']].name,
                species=[parameters.component_ids[i] for i in row['species']]))
            weights.append(row['weight'])
            identities.append(row['identity'])
    if '--declare-only' in sys.argv:
        return
    wheel = Path(os.environ['EPCSAFT_REPLAY_WHEEL']) if 'EPCSAFT_REPLAY_WHEEL' in os.environ else next((ROOT / 'build/environment-wheel').glob('*.whl'))
    subprocess.run([sys.executable, str(ROOT / 'tools/verify-installed-wheel.py'),
                    '--wheel', str(wheel)], check=True, cwd=ROOT)
    controls = reg.FitControls(maximum_iterations=3 if mode == 'short' else 1,
                               maximum_elapsed_time_seconds=120.0 if mode == 'short' else 1e-6)
    load_before = os.getloadavg()
    usage_before = resource.getrusage(resource.RUSAGE_SELF)
    wall, cpu = time.perf_counter(), time.process_time()
    with epcsaft.trace() as trace:
        result = reg.fit(parameters, coordinates, observations, weights=weights, controls=controls)
    elapsed, cpu_elapsed = time.perf_counter()-wall, time.process_time()-cpu
    usage_after = resource.getrusage(resource.RUSAGE_SELF)
    load_after = os.getloadavg()
    if 'EPCSAFT_REPLAY_PROFILE' in os.environ:
        from performance import profile_rows
        profile_rows(trace, Path(os.environ['EPCSAFT_REPLAY_PROFILE']))
    scalar = []
    if mode == 'short':
        coordinate_set = reg.CoordinateSet()
        coordinate_set.coordinates = coordinates
        model = epcsaft.Mixture(reg.fitted_parameters(parameters, result),
                               active_parameters=coordinate_set.active())
        for item in observations:
            prediction, gradient = reg.predict(model, item), reg.differentiate(model, item)
            if prediction.status != reg.ObservationStatus.Available or gradient.status != reg.ObservationStatus.Available:
                raise RuntimeError(f'scalar reference unavailable: {prediction.message}; {gradient.message}')
            scalar.append([prediction.value] + list(gradient.value))
    output = Path(output_name)
    output.parent.mkdir(parents=True, exist_ok=True)
    width = len(coordinates)
    with output.with_suffix('.csv').open('w', newline='') as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(['identity', 'prediction', 'residual', 'weighted_residual', 'status'] +
                        [f'physical_jacobian_{i}' for i in range(width)] +
                        [f'optimizer_jacobian_{i}' for i in range(width)] +
                        (['scalar_prediction'] + [f'scalar_physical_jacobian_{i}' for i in range(width)] if scalar else []))
        for i, name in enumerate(identities):
            writer.writerow([name, result.predictions[i], result.residuals[i],
                result.weighted_residuals[i], str(result.statuses[i])] +
                result.physical_jacobian[i*width:(i+1)*width] +
                result.optimizer_jacobian[i*width:(i+1)*width] + (scalar[i] if scalar else []))
    reuse = [event.values[0] for event in trace.events if event.name == 'reuse']
    summary = dict(wheel_sha256=hashlib.sha256(wheel.read_bytes()).hexdigest(),
        input_sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(), mode=mode,
        states=len(states), targets=len(observations), wall_s=elapsed, cpu_s=cpu_elapsed,
        user_s=usage_after.ru_utime-usage_before.ru_utime, system_s=usage_after.ru_stime-usage_before.ru_stime,
        load_before=load_before, load_after=load_after,
        initial_cost=result.initial_cost, final_cost=result.final_cost,
        physical=result.physical, status=result.status, message=result.message,
        residual_evaluations=result.residual_evaluations,
        jacobian_evaluations=result.jacobian_evaluations, trial_failures=result.trial_failures,
        iterations=result.iterations, residual_to_jacobian_promotions=reuse.count(3.0))
    output.write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary), flush=True)


def declarations():
    return main()


if __name__ == '__main__':
    main()

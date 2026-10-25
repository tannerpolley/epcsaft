"""Verify the #557 coordinate pair using frozen Engine-local declarations."""
import copy
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys

import epcsaft
from epcsaft import _core, regression as reg

ROOT = Path(__file__).resolve().parents[3]
LEAF = Path(__file__).resolve().parents[1]


def decode(value):
    if isinstance(value, list):
        return [decode(v) for v in value]
    if not isinstance(value, dict):
        return value
    if 'enum' in value:
        return getattr(getattr(_core, value['enum']), value['name'])
    cls = getattr(_core, value['type'])
    fields = {key: decode(v) for key, v in value['fields'].items()}
    try:
        obj = cls()
    except TypeError:
        return cls(**fields)
    for key, field in fields.items():
        setattr(obj, key, field)
    return obj


def declarations():
    path = LEAF / 'inputs/kim-svendsen-interval.json'
    identity = json.loads((path.parent / 'kim-svendsen-identity.json').read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest() == identity['declaration_sha256']
    kim = json.loads(path.read_text())
    raw = gzip.decompress((ROOT / 'tools/fixtures/mea-f1.json.gz').read_bytes())
    assert hashlib.sha256(raw).hexdigest() == identity['f1_uncompressed_sha256']
    packet = json.loads(raw)
    coordinates = decode(packet['coordinates']) + [
        reg.reaction_coordinate('R5', 'b', origin=-2677.91*math.log(10), scale=80*math.log(10),
                                bounds=(-3000*math.log(10), -2400*math.log(10))),
        reg.reaction_enthalpy_shift('R5', pivot_temperature=313.15, origin=0, scale=2500, bounds=(-10000, 10000)),
        reg.reaction_enthalpy_shift('R2', pivot_temperature=313.15, origin=0, scale=2500, bounds=(-10000, 10000))]
    rows = [packet['observations'][0], packet['observations'][136], kim['observation']]
    for row in rows:
        for name in ('problem', 'final_problem'):
            for reaction in row['fields'][name]['fields']['reactions']:
                if reaction['fields']['name'] == 'R5':
                    reaction['fields']['correlation']['fields']['b'] = -2677.91*math.log(10)
    # The F1-lineage input has no CO2/MEA cross edge; declared to keep that model (#569), not a physical-absence claim.
    kim['parameters']['topology']['no_cross_association'] = [{
        'component_id_a': 'carbon-dioxide', 'component_id_b': 'monoethanolamine',
        'provenance': {'source_id': 'kim-svendsen-interval', 'locator': 'inputs/kim-svendsen-interval.json: no CO2/MEA edge'}}]
    return epcsaft.Parameters.from_mapping(kim['parameters']), coordinates, rows, kim, identity


def main():
    parameters, coordinates, source, kim, identity = declarations()
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'build/issue-557'
    out.mkdir(parents=True, exist_ok=True)
    records, replays, centers = [], [], []
    for temperature in (293.15, 313.15, 353.15, 393.15):
        inputs = copy.deepcopy(source)
        for row in inputs:
            row['fields']['temperature'] = temperature
            for name in ('problem', 'final_problem'):
                row['fields'][name]['fields']['temperature']['fields']['value'] = temperature
        rows = decode(inputs)
        thermal = decode(kim['thermal_extrapolation']['thermochemistry'] if temperature == 293.15 else kim['thermochemistry'])
        chart = _core.CoordinateSet()
        chart.coordinates = coordinates
        start = chart.start(parameters._native, rows, thermal)
        def evaluate(point, declared=rows, coords=coordinates):
            with epcsaft.trace(diagnostics=True) as capture:
                value = _core.evaluate_observations(parameters._native, coords, point, [], declared, thermal)
            assert all(s == _core.ObservationStatus.Available for s in value.statuses)
            residual = max(s.convergence.residual for s in capture.trace.spans
                           if s.name == 'solve_equilibrium' and s.convergence is not None)
            return value, residual
        for plain in (2677.91, 2990.0):
            for shift in (0.0, 250.0):
                point = list(start)
                point[7:] = [(2677.91-plain)/80, shift/2500, -shift/2500]
                exact, central_residual = evaluate(point)
                laws = exact.chart.reaction_laws
                centers.append(dict(T_K=temperature, plain_a_k_K=plain, pre_applied_J_mol=0.0,
                    incremental_J_mol=[shift, -shift], total_J_mol=[shift, -shift], physical=exact.chart.physical,
                    laws={name: {key: getattr(law, key) for key in kim['observation']['fields']['problem']['fields']['reactions'][4]['fields']['correlation']['fields']}
                          for name, law in laws.items()}))
                replay = decode(inputs)
                for row in replay:
                    for problem in (row.problem, row.final_problem):
                        for reaction in problem.reactions:
                            if reaction.name in laws:
                                reaction.correlation = laws[reaction.name]
                model = _core.HelmholtzModel(exact.chart.parameters, exact.chart.active, thermal)
                repeated = [_core.predict(model, row).value for row in replay]
                assert repeated == exact.predictions
                replays.append(dict(T_K=temperature, plain_a_k_K=plain, shift_J_mol=shift,
                                    maximum_error=max(abs(a-b) for a, b in zip(repeated, exact.predictions))))
                for column, label, step, du in ((7, 'R5_a_k', 2.0, -1/80), (8, 'R5_h', 25.0, 1/2500), (9, 'R2_h', 25.0, 1/2500)):
                    rungs = []
                    for rung in (step, step/2):
                        plus, minus = list(point), list(point)
                        plus[column] += rung*du
                        minus[column] -= rung*du
                        vp, rp = evaluate(plus)
                        vm, rm = evaluate(minus)
                        rungs.append((rung, vp, vm, rp, rm))
                    for row, kind in enumerate(('pressure', 'species', 'kim-addition')):
                        slope = 1/exact.predictions[row] if rows[row].form == _core.ResidualForm.LogRatio else 1/rows[row].scale
                        for quantity, factor, scale in (('prediction_physical', 1.0, du), ('prediction_optimizer', 1.0, 1.0),
                                                       ('residual_optimizer', slope, 1.0), ('weighted_residual_optimizer', slope*math.sqrt((1, 4, .25)[row]), 1.0)):
                            analytic = (exact.residual_optimizer_gradient[row*10+column]*factor/slope
                                        if 'residual' in quantity else exact.optimizer_gradient[row*10+column]*scale)
                            differences = [(vp.predictions[row]-vm.predictions[row])/(2*h*du)*factor*scale
                                           for h, vp, vm, rp, rm in rungs]
                            if 'residual' in quantity and rows[row].form == _core.ResidualForm.LogRatio:
                                differences = [(math.log(vp.predictions[row])-math.log(vm.predictions[row]))/(2*h*du)*factor/slope
                                               for h, vp, vm, rp, rm in rungs]
                            richardson = (4*differences[1]-differences[0])/3
                            error = abs(richardson-analytic)
                            passed = error <= max(1e-8, 1e-5*max(abs(analytic), abs(richardson)))
                            records.append(dict(T_K=temperature, plain_a_k_K=plain, shift_J_mol=shift, coordinate=label, kind=kind,
                                quantity=quantity, prediction=exact.predictions[row], exact=analytic, step=step, half_step=step/2,
                                D_h=differences[0], D_half=differences[1], raw_error_h=abs(differences[0]-analytic),
                                raw_error_half=abs(differences[1]-analytic), richardson=richardson, absolute_error=error,
                                relative_error=error/max(abs(analytic), abs(richardson), 1e-300), passed=passed,
                                signal_h=abs(rungs[0][1].predictions[row]-rungs[0][2].predictions[row]),
                                central_residual=central_residual, residual_plus=rungs[0][3], residual_minus=rungs[0][4],
                                residual_half_plus=rungs[1][3], residual_half_minus=rungs[1][4]))
                print(temperature, plain, shift, 'completed', flush=True)
    for name, values in (('columns', records), ('replay', replays)):
        with (out / f'{name}.csv').open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(values[0]))
            writer.writeheader()
            writer.writerows(values)
    (out / 'native-laws.json').write_text(json.dumps(centers, indent=2)+'\n')
    (out / 'input-identities.json').write_text(json.dumps(identity, indent=2)+'\n')
    assert all(row['passed'] for row in records), 'Richardson comparison failed; retain and diagnose'


if __name__ == '__main__':
    main()

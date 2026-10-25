"""Native #658 calibration with the complete frozen joint-check inputs."""
import hashlib
import json
import math
from pathlib import Path
import sys

import epcsaft
from epcsaft import Parameters, regression as r

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'analyses/2026-lithium-chloride-refit/scripts'))
import fit as li
import evaluate as checks
from trace_activity import Activities, COMPONENTS as TRACE_COMPONENTS

OUT = HERE / 'results/cation-mixing'
D = json.loads((HERE / 'inputs/cation-mixing-declaration.json').read_text())
POINT = json.loads((HERE / 'inputs/cation-mixing-reference.json').read_text())
MIXING = {'solvent': 'water', 'cations': ['lithium-cation', 'magnesium-cation'], 'anion': 'chloride-anion',
          'coefficients': [{'identity': 'mixing/Li-Mg-Cl/theta', 'family': 'cation_pair_mixing',
                            'value': {'magnitude': -0.2205, 'unit': 'kilogram / mole'}},
                           {'identity': 'mixing/Li-Mg-Cl/psi', 'family': 'cation_triple_mixing',
                            'value': {'magnitude': 0.024, 'unit': 'kilogram ** 2 / mole ** 2'}}]}


def document(physical, mixing=True):
    doc = json.loads((ROOT / f"data/packets/{li.D['base_packet']}/parameters/parameters.json").read_text())
    values = dict(zip(li.LABELS, physical, strict=True))
    for owner, coefficients in [(c['component_id'], c.get('coefficients', [])) for c in doc['components']] + [
            ('+'.join(sorted((q['component_id_a'], q['component_id_b']))), q['coefficients']) for q in doc['pairs']]:
        for coefficient in coefficients:
            if (label := f"{coefficient['family']}/{owner}") in values:
                coefficient['value']['magnitude'] = values[label]
    if mixing:
        doc['cation_mixing'] = MIXING
    return doc


def parameters(physical, mixing=True, components=li.COMPONENTS):
    return Parameters.from_mapping(document(physical, mixing), components=components)


def observations(p):
    training = li.training(p)
    for row in training:
        row['observation'].scale = D['observation_scales'][row['source']]
    rows = [q for q in checks.mixture_rows('Yao') if q['bucket'] in ('I<=3', '3<I<=8.3')]
    obs = [q['observation'] for q in training]
    for q in rows:
        obs.append(r.mixture_osmotic_coefficient(p, li.T, li.P, li.composition(p, {
            li.ION: q['a'], 'magnesium-cation': q['b'], 'chloride-anion': q['a'] + 2*q['b']}),
            observed=q['phi'], scale=D['observation_scales'][q['bucket']]))
    assert len(obs) == 96
    return obs


def score(p, physical):
    water = checks.Water(p)
    nominal, rounds, buckets = checks.held_out(water, 'Yao', 'native-mixing')
    checks.fit.write(OUT / 'yao-nominal.csv', nominal)
    rounds = [dict(original_passed=q.pop('passed'), original_mixing_passed=q.pop('mixing_passed'), **q,
        amended_passed=q['RMS'] <= .0105 and q['maximum_absolute_error'] <= .027 if q['bucket'] != 'report-only' else None) for q in rounds]
    checks.fit.write(OUT / 'yao-rounding.csv', rounds)
    yao = {b['bucket']: {k: b[k] for k in ('RMS', 'maximum_absolute_error', 'worst_RMS', 'worst_maximum')}
           for b in buckets if b['bucket'] in ('I<=3', '3<I<=8.3')}
    for bucket, values in yao.items():
        scored = [v for v in rounds if v['bucket'] == bucket]
        values['passing_realizations'] = sum(v['RMS'] <= .0105 and v['maximum_absolute_error'] <= .027 for v in scored)
        values['realizations'] = len(scored)
    hw, eg = [], []
    for q in li.training(p):
        pred = water.predict(li=q['m_LiCl'], na=q['m_NaCl'])
        if q['source'] == 'Hamer-Wu':
            hw.append(pred['phi'] - q['observation'].observed)
        else:
            eg.append((q['m_LiCl']+q['m_NaCl'], pred['aw']-math.exp(-p.molar_masses[0]*2*(q['m_LiCl']+q['m_NaCl'])*q['observed'])))
    trace = Activities(parameters(physical, components=TRACE_COMPONENTS))
    transfers, density = checks.transfers(p)
    frozen, _ = checks.transfers(li.parameters())
    pv = {v['state']: checks.combinations(v['log10_activity']) for v in json.loads((li.RESULTS/'pitzer.json').read_text())['states']}
    lithium = []
    for key, value in transfers.items():
        state, name = key.split('/')
        lithium.append({'row': key, 'check': name, 'engine_transfer': value,
                        'difference_decade': value - (pv[state][name]-pv['fig7'][name]),
                        'change_from_frozen_/1': value-frozen[key]})
    return {'Yao': yao, 'Hamer_Wu': checks.metrics(hw),
            'El_Guendouzi_aw': {b: checks.metrics([e for I, e in eg if (I <= 3) == low])
                               for b, low in [('I<=3', True), ('I>3', False)]},
            'F3': [dict(trace.row(m), binary_route_error_decade=trace.binary_check(m)) for m in (1, 3, 6)],
            'Lithium': lithium, 'density': density}


def declarations():
    p = parameters(D['start_physical'])
    return p, li.coordinates(p), observations(p)


if __name__ == '__main__':
    p = parameters(D['start_physical'])
    obs = observations(p)
    result = r.fit(p, li.coordinates(p), obs, controls=r.FitControls(**D['controls']), starts=[D['start_optimizer']])[0]
    sigma, dh, born, kw, kcl = result.optimizer
    physical = [sigma, .88*sigma, dh, born, kw, kcl]
    record = {k: getattr(result, k) for k in ('status', 'message', 'usable', 'iterations', 'residual_evaluations',
              'jacobian_evaluations', 'trial_failures', 'initial_cost', 'final_cost')}
    record.update(physical=physical, optimizer=list(result.optimizer), available_rows=sum(str(s)=="ObservationStatus.Available" for s in result.statuses),
                  covariance={k: getattr(result.covariance, k) for k in ('status', 'rank', 'active_bounds', 'singular_values', 'projected_gradient_norm')},
                  failures=[f.message for f in result.failures],
                  extension=epcsaft._core.__file__, extension_sha256=hashlib.sha256(Path(epcsaft._core.__file__).read_bytes()).hexdigest())
    li.write(OUT/'training.csv', [{'row':i+1, 'observed':ob.observed, 'scale':ob.scale, 'prediction':v,
        'residual':u, 'weighted_residual':w} for i,(ob,v,u,w) in enumerate(zip(obs,result.predictions,
        result.residuals,result.weighted_residuals,strict=True))])
    (OUT/'native-fit.json').write_text(json.dumps(record, indent=2, default=list)+'\n')
    fitted = r.fitted_parameters(p, result)
    assert isinstance(fitted._document, epcsaft._core.ParameterDocument), fitted._document
    scores = score(fitted, physical)
    (OUT/'native-scores.json').write_text(json.dumps(scores, indent=2)+'\n')
    print(json.dumps(scores, indent=2), flush=True)

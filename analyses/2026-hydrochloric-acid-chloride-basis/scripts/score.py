"""Score only retained Hamer–Wu predictions, following design-stage2 section 2a.i."""
import math
import pandas as pd
from run import C, D


def score(leaf, result):
    assert (result / 'outcome.txt').read_text().startswith('Identified'), '2a requires identification'
    decisions = pd.read_csv(result / 'identification.csv')
    selected = decisions[decisions.selected == 1]
    assert len(selected) == 1
    values = pd.read_csv(result / 'fit-values.csv', float_precision='round_trip')
    values = values[(values.stage == selected.iloc[0].stage) & (values.start == selected.iloc[0].start)]
    predictions = values[values.field == 'predictions'].sort_values('index')
    statuses = values[values.field == 'observation_status'].sort_values('index')
    assert predictions['index'].tolist() == statuses['index'].tolist() == list(range(D.get('training_observations', 2 * D['training_states'])))
    assert all(math.isfinite(v) for v in predictions.value) and all(statuses.value == 1)  # ObservationStatus::Available, observations.hpp
    if D.get('sulfate'):
        targets = pd.read_csv(result / 'mixture-admission.csv', float_precision='round_trip')
        assert targets.observation.tolist() == predictions['index'].tolist()
        targets['phi'] = targets.observed
        targets['prediction'] = predictions.value.to_numpy()
        targets = targets[targets.salt != 'LiCl-MgSO4'].copy()
    else:
        predictions = predictions.iloc[:2 * D['training_states']]
        targets = pd.read_csv(leaf / 'inputs' / C['training'], float_precision='round_trip').loc[lambda x: x.molality <= C['cutoff']].copy()
    assert len(targets) == D['training_states']
    if not D.get('sulfate'): assert sum(targets.molality <= C['split']) == C['low_states']
    targets['phi_error'] = (targets.prediction if D.get('sulfate') else predictions.value.to_numpy()[::2]) - targets.phi
    if not D.get('sulfate'): targets['gamma_relative_error_percent'] = 100 * (predictions.value.to_numpy()[1::2] / targets.gamma - 1)
    rows = []
    limits = D.get('accuracy', dict(phi_RMS=.01, phi_max=.03, gamma_MARD_percent=2, gamma_max_percent=5))
    buckets = [(name, targets.bucket == name) for name in ['Na-low','Na-high','Mg-low','Mg-high']] if D.get('sulfate') else [(f"<={C['split']}", targets.molality <= C['split']), (f"{C['split']}-{C['cutoff']}", targets.molality > C['split'])]
    for bucket, mask in buckets:
        group = targets[mask]
        if D.get('sulfate'): assert len(group) == dict(zip(['Na-low','Na-high','Mg-low','Mg-high'], [5,4,9,5]))[bucket]
        rms = math.sqrt(sum(group.phi_error ** 2) / len(group))
        maximum = group.phi_error.abs().max()
        if not D.get('sulfate'): mard, relative_max = group.gamma_relative_error_percent.abs().mean(), group.gamma_relative_error_percent.abs().max()
        metrics = dict(phi_RMS=rms, phi_max=maximum)
        if not D.get('sulfate'): metrics.update(gamma_MARD_percent=mard, gamma_max_percent=relative_max)
        rows.append(dict(bucket_mol_per_kg=bucket, states=len(group), phi_RMS=rms, phi_RMS_limit=limits['phi_RMS'], phi_max=maximum, phi_max_limit=limits['phi_max'], passes=all(v <= limits[k] for k,v in metrics.items())))
        if not D.get('sulfate'):
            rows[-1] = dict(bucket_mol_per_kg=bucket, states=len(group), phi_RMS=rms, phi_RMS_limit=limits['phi_RMS'], phi_max=maximum, phi_max_limit=limits['phi_max'], gamma_MARD_percent=mard, gamma_MARD_limit_percent=limits['gamma_MARD_percent'], gamma_max_percent=relative_max, gamma_max_limit_percent=limits['gamma_max_percent'], passes=rows[-1]['passes'])
        if 'trial' in D: rows[-1].update({k + '_passes': v <= limits[k] for k,v in metrics.items()})
    return targets, pd.DataFrame(rows)

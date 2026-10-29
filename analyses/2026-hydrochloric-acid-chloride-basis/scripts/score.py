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


def profile_scores(result):
    """Score every retained attempt; comparison rows are never selection inputs."""
    rows = pd.read_csv(result / 'training.csv', float_precision='round_trip')
    attempts = pd.read_csv(result / 'attempts.csv', float_precision='round_trip')
    metrics = []
    for (point, start), group in rows.groupby(['point','start'], sort=False):
        for quantity, source, limits in [(q,source,(D['profile']['mean_limits'][i],D['profile']['max_limits'][i])) for i,(q,source) in enumerate([('phi','binary'),('gamma','binary'),('aw','Messnaoui')])]:
            for high in [False,True]:
                frame = group[(group.source == source) & (group.observable == quantity) & ((group.molality > 3) == high)]
                error = frame.prediction-frame.observed if quantity != 'gamma' else 100*(frame.prediction/frame.observed-1)
                mean = error.abs().mean() if quantity == 'gamma' else math.sqrt(error.pow(2).mean())
                metrics.append(dict(point=point,start=start,quantity=quantity,bucket='high' if high else 'low',states=len(frame),mean=mean,worst=error.abs().max(),mean_limit=limits[0],worst_limit=limits[1],verified=frame.verified.eq(1).all()))
        frame=group[group.source == 'Jiang-I5']
        for role,subset in [('mixtures',frame[frame.background>0]),('endpoint',frame[frame.background==0]),('line',frame)]:
            error=100*(subset.prediction/subset.observed-1)
            metrics.append(dict(point=point,start=start,quantity='Jiang-I5',bucket=role,states=len(subset),mean=error.abs().mean(),worst=error.abs().max(),signed_mean=error.mean(),verified=subset.verified.eq(1).all()))
    pd.DataFrame(metrics).to_csv(result/'metrics.csv',index=False)
    ordinary=[]
    for point, frame in attempts[attempts.point>0].groupby('point',sort=False):
        finite=frame[(frame.converged==1) & frame.cost.notna()]
        if finite.empty:
            continue
        best=finite.sort_values(['cost','start'],kind='stable').iloc[0].copy()
        comparable=frame.cost <= best.cost+max(1e-10,1e-4*abs(best.cost))
        best['complete']=bool(len(frame)==3 and frame.converged.eq(1).all() and frame.verified.eq(1).all() and frame.u.between(-50+1e-6,1-1e-6,inclusive='neither').all() and ((frame.loc[comparable,'u']-best.u).abs()<=1e-6).all())
        ordinary.append(best)
    profile=pd.DataFrame(ordinary)
    profile.to_csv(result/'profile.csv',index=False)
    return rows,attempts,profile


def select_profile(attempts, profile):
    successful=attempts[(attempts.point>0) & attempts.converged.eq(1) & attempts.verified.eq(1) & (attempts.B<=1) & (attempts.V<=1)]
    if not successful.empty:
        minimum=successful.cost.min()
        comparable=successful[successful.cost<=minimum+max(1e-10,1e-4*abs(minimum))].copy()
        return comparable.assign(abs_c=comparable.c.abs()).sort_values(['abs_c','worst_I5','point','start'],kind='stable').iloc[0]
    complete=profile[profile.complete].copy()
    eligible=complete[complete.B<=1]
    candidates=eligible if not eligible.empty else complete
    return candidates.assign(abs_c=candidates.c.abs()).sort_values((['V'] if not eligible.empty else ['B','V'])+['cost','abs_c','point'],kind='stable').iloc[0]

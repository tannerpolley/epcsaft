"""Rescore the selected endpoint, then evaluate required adjusted binary water rows."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path
import pandas as pd
from run import D, LEAF, ROOT
from score import score

if D.get('profile'):
    import os
    import time
    from score import profile_scores, select_profile
    R=Path(sys.argv[2]).resolve()
    rows,attempts,profile=profile_scores(R)
    frozen=pd.read_csv(R/'frozen.csv',float_precision='round_trip').iloc[0]
    original=frozen.to_dict()
    selection=select_profile(attempts,profile)
    assert (selection.point,selection.start)==(frozen.point,frozen.start)
    if not ((profile.point==frozen.point) & (profile.start==frozen.start)).any():
        profile=pd.concat([profile,selection.to_frame().T],ignore_index=True)
    executable=ROOT/'build/native-drivers/activity/driver'
    environment={**os.environ,'LD_LIBRARY_PATH':str(ROOT/'build/engine-meson/shared')}
    summaries=[]
    for endpoint in profile.itertuples(index=False):
        folder=R/f'point-{int(endpoint.point)}'; folder.mkdir()
        for source,mixtures,endpoints in [('jiang1996-hcl-li',27,3),('roy1980-hcl-mg',47,11),('khoo1977-hcl-mg',25,5),('khoo1977-hcl-ca',25,5)]:
            destination=folder/source; destination.mkdir()
            began=time.monotonic()
            result=subprocess.run([str(executable),str(LEAF),str(destination),f'source-{source}:100000',str(endpoint.u),str(endpoint.c)],env=environment)
            frame=pd.read_csv(destination/'source.csv',float_precision='round_trip')
            if source=='jiang1996-hcl-li': frame=frame[frame.ionic_strength!=5].copy()
            assert len(frame)==mixtures+endpoints
            frame['relative_error_percent']=100*(frame.gamma_prediction/frame.gamma_observed-1)
            frame.to_csv(destination/'errors.csv',index=False)
            for endpoint_role,group in frame.groupby(frame.y_background==0):
                for label,g in [('whole',group),*((str(i),g) for i,g in group.groupby('ionic_strength'))]:
                    errors=g.relative_error_percent
                    complete=g.passes.eq(1).all() and errors.notna().all()
                    summaries.append(dict(point=endpoint.point,start=endpoint.start,source=source,role='endpoint' if endpoint_role else 'mixture',group=label,states=len(g),complete=complete,mean=errors.abs().mean() if complete else None,worst=errors.abs().max() if complete else None,signed_mean=errors.mean() if complete else None,wall_seconds=time.monotonic()-began,returncode=result.returncode))
    spans=[]
    for folder in R.glob('point-*'):
        for source in folder.iterdir():
            frame=pd.read_csv(source/'errors.csv',float_precision='round_trip')
            for strength,group in frame.groupby('ionic_strength'):
                ordered=group.sort_values('y_background'); lo,hi=ordered.iloc[0],ordered.iloc[-1]
                spans.append(dict(point=folder.name,source=source.name,I=strength,y_min=lo.y_background,y_max=hi.y_background,D=hi.log_gamma-lo.log_gamma))
    pd.DataFrame(spans).to_csv(R/'transfer-spans.csv',index=False)
    pd.DataFrame(summaries).to_csv(R/'transfer.csv',index=False)
    assert pd.read_csv(R/'frozen.csv',float_precision='round_trip').iloc[0].to_dict()==original
    unavailable=pd.DataFrame(summaries).assign(complete=False,mean=float('nan'),worst=float('nan'),signed_mean=float('nan'))
    repeated=select_profile(attempts,profile)
    assert (repeated.point,repeated.start)==(frozen.point,frozen.start)
    (R/'selection-invariance.json').write_text(json.dumps(dict(before=original,after_point=float(repeated.point),after_start=float(repeated.start),unavailable_summary_memberships=int(unavailable.states.sum()),equal=True))+'\n')
    selected=R/f'point-{int(frozen.point)}'
    replay=R/'replay'; replay.mkdir()
    subprocess.run([str(executable),str(LEAF),str(replay),'profile-replay',str(frozen.u),str(frozen.c)],env=environment,check=True)
    expected=rows[(rows.point==frozen.point) & (rows.start==frozen.start)].reset_index(drop=True)
    current=pd.read_csv(replay/'training.csv',float_precision='round_trip')
    assert len(current)==len(expected)==84 and current.verified.eq(1).all()
    differences=(current.prediction-expected.prediction).abs()
    tolerance=(expected.prediction.abs()*1e-10).clip(lower=1e-12)
    phi=expected.observable=='phi'
    tolerance.loc[phi]=tolerance.loc[phi].combine(1e-12/(2*.0180153*expected.loc[phi,'molality']),max)
    pd.DataFrame(dict(difference=differences,allowance=tolerance,passes=differences<=tolerance)).to_csv(replay/'training-comparison.csv',index=False)
    assert (differences<=tolerance).all() and current.cause.fillna('').equals(expected.cause.fillna(''))
    comparisons=[]
    for source in ['jiang1996-hcl-li','roy1980-hcl-mg','khoo1977-hcl-mg','khoo1977-hcl-ca']:
        destination=replay/source; destination.mkdir()
        began=time.monotonic()
        result=subprocess.run([str(executable),str(LEAF),str(destination),f'source-{source}:100000',str(frozen.u),str(frozen.c)],env=environment)
        a=pd.read_csv(selected/source/'source.csv',float_precision='round_trip'); b=pd.read_csv(destination/'source.csv',float_precision='round_trip')
        assert a.shape==b.shape and a.isna().equals(b.isna())
        for column in a:
            if pd.api.types.is_numeric_dtype(a[column]):
                delta=(a[column]-b[column]).abs(); tolerance=(a[column].abs()*1e-10).clip(lower=1e-12)
                assert ((delta<=tolerance) | (a[column].isna() & b[column].isna())).all(),column
            else: assert a[column].fillna('').equals(b[column].fillna('')),column
        comparisons.append(dict(source=source,rows=len(a),equal=True,wall_seconds=time.monotonic()-began,returncode=result.returncode))
    pd.DataFrame(comparisons).to_csv(replay/'summary.csv',index=False)
    shared={}
    for path in replay.glob('*/*.csv'):
        retained=selected/path.relative_to(replay)
        if retained.exists() and path.read_bytes()==retained.read_bytes():
            shared[str(path.relative_to(replay))]=dict(retained=str(retained.relative_to(R)),sha256=hashlib.sha256(path.read_bytes()).hexdigest()); path.unlink()
    (replay/'equivalence.json').write_text(json.dumps(shared,indent=2)+'\n')
    eligible=profile[profile.complete & (profile.B<=1)]
    metric='V' if not eligible.empty else 'B'
    edge=profile[profile.c.isin([-10,-9.75])]
    edge_incomplete=len(edge)!=2 or not edge.complete.all()
    improving=not edge_incomplete and float(edge[edge.c==-10][metric].iloc[0])<(1-1e-6)*float(edge[edge.c==-9.75][metric].iloc[0])
    boundary=frozen.c in [-10,0] or improving
    disposition=1 if frozen.success else 3 if len(profile)!=attempts[attempts.point>0].point.nunique() or not profile.complete.all() or edge_incomplete or boundary else 2
    (R/'disposition.json').write_text(json.dumps(dict(disposition=disposition,edge='incomplete' if edge_incomplete else 'boundary-limited' if boundary else 'interior',transfer_complete=all(row['complete'] for row in summaries),optional_pressure='not run under optional pressure allowance'))+'\n')
    raise SystemExit(0)

R = LEAF / 'results/binary-water-qualification'
R.mkdir(exist_ok=False)
D['accuracy'] = dict(phi_RMS=.011, phi_max=.03, gamma_MARD_percent=2, gamma_max_percent=5)
errors, metrics = score(LEAF, LEAF / 'results/issue-798-r2')
assert errors.to_csv(index=False) == (LEAF / 'results/issue-798-retrial/issue-798-r2-binary-errors.csv').read_text()
metrics.to_csv(R / 'binary-metrics.csv', index=False)
assert metrics.passes.all(), '2a.i failed; no later evaluation permitted'
subprocess.run([sys.executable, str(LEAF / 'scripts/run.py'), str(LEAF / 'inputs/fit-declaration-issue-798-r2.json'), str(R / 'native'), '--binary-water'], cwd=ROOT, check=True)
water = pd.read_csv(R / 'native/water.csv', float_precision='round_trip')
assert len(water) == 14 and water.molality.between(.2, 6).all()
rows = []
for name, group in water.groupby(water.molality > 3):
    rms, maximum = (group.error.pow(2).mean())**.5, group.error.abs().max()
    available = group.aw_prediction.notna().all() and group.cause.isna().all()
    rows.append(dict(bucket='3-6' if name else '<=3', states=len(group), aw_RMS=rms, RMS_limit=.003, aw_max=maximum, max_limit=.006, passes=bool(available and rms <= .003 and maximum <= .006)))
pd.DataFrame(rows).to_csv(R / 'water-metrics.csv', index=False)
(R / 'input-identities.json').write_text(json.dumps({name: hashlib.sha256((LEAF / 'inputs' / name).read_bytes()).hexdigest() for name in ['messnaoui2020-hcl.csv','fit-declaration-issue-798-r2.json']}, indent=2) + '\n')
if not all(row['passes'] for row in rows):
    (R / 'outcome.txt').write_text('2a.i passes with osmotic RMS limit 0.011. First failure: 2a.ii adjusted Messnaoui water activity. Report-only binary comparisons, 2b numerical checks and reserved qualification, and 2c preservation/materialization not reached. No packet written.\n')
    sys.exit(1)
executable = ROOT / 'build/native-drivers/activity/driver'
environment = {**__import__('os').environ, 'LD_LIBRARY_PATH': str(ROOT / 'build/engine-meson/shared')}
def evaluate(destination, mode):
    destination.mkdir()
    subprocess.run([str(executable), str(LEAF), str(destination), mode], env=environment, check=True)
evaluate(R / 'report-only', 'report-binary')
evaluate(R / 'numerical', 'numerical')
evaluate(R / 'cold', 'numerical')
warm = pd.read_csv(R / 'numerical/numerical.csv', float_precision='round_trip')
cold = pd.read_csv(R / 'cold/numerical.csv', float_precision='round_trip')
assert len(warm) == 14 and warm.passes.eq(1).all() and warm.equals(cold)
training = pd.read_csv(R / 'cold/training.csv', float_precision='round_trip')
retained = pd.read_csv(LEAF / 'results/issue-798-r2/fit-values.csv', float_precision='round_trip').query('stage == 1 and start == 3 and field == "predictions"').sort_values('index').value.to_numpy()
current = training[['phi','gamma']].to_numpy().reshape(-1)
water_mass = __import__('run').declarations().molar_masses[0]
allowance = [max(1e-10*abs(v), 1e-12/(water_mass*2*m) if i%2 == 0 else 1e-12) for i,(v,m) in enumerate(zip(retained, training.molality.repeat(2)))]
difference = abs(current-retained)
pd.DataFrame(dict(observation=range(60), retained=retained, cold=current, absolute_difference=difference, allowance=allowance, passes=difference<=allowance)).to_csv(R / 'cold/training-comparison.csv', index=False)
assert all(difference <= allowance)
shared = {}
for directory, mode, filename in [('native', 'binary-water', 'water.csv'), ('report-only', 'report-binary', 'report.csv')]:
    destination = R / 'cold' / directory
    evaluate(destination, mode)
    original = R / directory / filename
    replay = destination / filename
    a, b = [pd.read_csv(path, float_precision='round_trip') for path in [original, replay]]
    assert a.shape == b.shape and a.isna().equals(b.isna())
    comparisons = []
    for column in a:
        if pd.api.types.is_numeric_dtype(a[column]):
            allowance = (a[column].abs()*1e-10).clip(lower=1e-12)
            if column == 'phi':
                allowance = allowance.combine(1e-12/(water_mass*2*a.molality), max)
            difference = (b[column]-a[column]).abs()
            passes = (difference <= allowance) | (a[column].isna() & b[column].isna())
            assert passes.all(), column
            comparisons.extend(dict(row=i, quantity=column, retained=a[column][i], cold=b[column][i], absolute_difference=difference[i], allowance=allowance[i], passes=passes[i]) for i in a.index)
        else:
            assert a[column].fillna('').equals(b[column].fillna('')), column
    pd.DataFrame(comparisons).to_csv(destination / 'comparison.csv', index=False)
    shared[directory] = dict(states=len(a), availability_and_causes_equal=True, retained='../' + directory + '/' + filename, sha256=hashlib.sha256(original.read_bytes()).hexdigest(), replay_sha256=hashlib.sha256(replay.read_bytes()).hexdigest())
    if replay.read_bytes() == original.read_bytes():
        replay.unlink()

for path in (R / 'cold').glob('*.csv'):
    original = R / 'numerical' / path.name
    if original.exists() and path.read_bytes() == original.read_bytes():
        shared[path.name] = dict(retained='../numerical/' + path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        path.unlink()
(R / 'cold/equivalence.json').write_text(json.dumps(shared, indent=2) + '\n')
for source, mixtures, endpoints in [('roy1980-hcl-mg',47,11), ('khoo1977-hcl-mg',25,5), ('jiang1996-hcl-li',36,4), ('khoo1977-hcl-ca',25,5)]:
    folder = R / source
    folder.mkdir()
    observed = pd.read_csv(LEAF / 'inputs' / (source + '.csv'), float_precision='round_trip')
    (folder / 'input.sha256').write_text(hashlib.sha256((LEAF / 'inputs' / (source + '.csv')).read_bytes()).hexdigest() + '\n')
    metrics, frames = [], []
    for pressure in [100000,101325]:
        destination = folder / str(pressure)
        try:
            evaluate(destination, f'source-{source}:{pressure}')
            evaluate(destination / 'cold', f'source-{source}:{pressure}')
        except subprocess.CalledProcessError:
            (R / 'outcome.txt').write_text(f'First numerical failure: {source} at {pressure} Pa; later sources, preservation and packet write not reached.\n')
            raise
        frame = pd.read_csv(destination / 'source.csv', float_precision='round_trip')
        cold = pd.read_csv(destination / 'cold/source.csv', float_precision='round_trip')
        assert frame.equals(cold) and frame.passes.eq(1).all()
        assert len(frame) == mixtures + endpoints and sum(frame.y_background > 0) == mixtures and sum(frame.y_background == 0) == endpoints
        frame['relative_error_percent'] = 100*(frame.gamma_prediction/frame.gamma_observed-1)
        frame['locator'] = observed.locator
        frame.to_csv(destination / 'errors.csv', index=False)
        frames.append(frame)
        for endpoint, group in frame.groupby(frame.y_background == 0):
            for name, rows in [('whole',group), *(((str(i),g) for i,g in group.groupby('ionic_strength')) if not endpoint else [])]:
                metrics.append(dict(pressure_Pa=pressure, role='endpoint' if endpoint else 'mixture', group=name, states=len(rows), MARD_percent=rows.relative_error_percent.abs().mean(), max_percent=rows.relative_error_percent.abs().max(), signed_mean_percent=rows.relative_error_percent.mean()))
    pd.DataFrame(metrics).to_csv(folder / 'pressure-metrics.csv', index=False)
    (frames[1].select_dtypes('number')-frames[0].select_dtypes('number')).to_csv(folder / 'pressure-difference.csv', index=False)

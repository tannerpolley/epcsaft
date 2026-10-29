"""Retained block medians, paired ratios and release-owned physical-value checks."""
import csv
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import statistics as st

HERE = Path(__file__).resolve().parents[1]; OUT = HERE/'results'
spec = importlib.util.spec_from_file_location('release_values', HERE.parents[1]/'analyses/2026-release-computational-cost/scripts/verify.py')
v = importlib.util.module_from_spec(spec); spec.loader.exec_module(v)


def physical(sample, arm):
    return v.values('explicit-discovery', 'Engine' if arm.startswith('Engine') else arm, sample)


def check(reference, actual, arm):
    if not arm.startswith('Engine'):
        return all(r['passed'] for r in v.compare('explicit-discovery',arm,reference,actual))
    a,b = v.flat(reference),v.flat(actual)
    return a.keys()==b.keys() and all(math.isclose(a[k],b[k],rel_tol=1e-10,abs_tol=1e-12) for k in a)


def write(name, rows):
    assert rows, name
    with (OUT/name).open('w') as f:
        w = csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def main():
    selected, admissions = {}, []
    for path in sorted((OUT/'raw').glob('measurement-*.jsonl.gz'),key=lambda p:int(p.name.split('-')[1].split('.')[0])):
        records = list(map(json.loads,gzip.open(path,'rt')))
        identities = next(r['identities'] for r in records if r['kind']=='identities')
        matched = all(identities[arm]['image_sha256']==hashlib.sha256((OUT/f'image-{image}.json').read_bytes()).hexdigest() for arm,image in [('Engine-A','A'),('Engine-control','A'),('Engine-B','B'),('Engine-B-control','B')])
        blocks = [r for r in records if r['kind']=='block']
        for op in {r['operation'] for r in blocks}:
            rows = [r for r in blocks if r['operation']==op]
            quiet = all(s['load'][0]<4 and not s['compilers'] for r in rows for s in (r['before'],r['after']))
            admissions.append(dict(path=str(path.relative_to(HERE)),operation=op,quiet=quiet,blocks=len(rows),matched_images=matched))
            if not matched: continue
            if quiet: selected[op] = rows
            else: selected.pop(op,None)
    assert set(selected)=={'ordinary','prepared'}
    pre = json.loads((OUT/'preflight.json').read_text()); blocks, checks, times = [], [], {}
    for op,rows in selected.items():
        assert len(rows)==18
        reference = physical(pre[op]['Engine-A']['samples'][0],'Engine-A')
        for r in rows:
            reply,arm,rep = r['reply'],r['arm'],r['repetition']
            assert 'error' not in reply and len(reply['samples'])==25, reply
            for i,sample in enumerate(reply['samples']):
                ok = check(reference,physical(sample,arm),arm)
                checks.append(dict(operation=op,arm=arm,repetition=rep,sample=i,passed=ok)); assert ok
            wall = [s['seconds']*1e6 for s in reply['samples']]
            cpu = [s['cpu_seconds']*1e6 for s in reply['samples'] if 'cpu_seconds' in s]
            assert (op,arm,rep) not in times
            times[op,arm,rep] = st.median(wall)
            blocks.append(dict(operation=op,arm=arm,repetition=rep,wall_median_us=st.median(wall),
                wall_min_us=min(wall),wall_max_us=max(wall),cpu_median_us=st.median(cpu) if cpu else None,
                prepare_us=reply['prepare_seconds']*1e6 if reply.get('prepare_seconds') is not None else None,
                load_before=r['before']['load'][0],load_after=r['after']['load'][0]))
    pairs = []
    for op in selected:
        for rep in range(3):
            for label,num,den in [('A-control/A','Engine-control','Engine-A'),('B-control/B','Engine-B-control','Engine-B')]+[
                ('B/'+den,b,den) for b in ('Engine-B','Engine-B-control') for den in ('Engine-A','FeOs','Clapeyron')]:
                pairs.append(dict(operation=op,comparison=label.replace('Engine-A','A'),repetition=rep,arm=num,ratio=times[op,num,rep]/times[op,den,rep]))
    ratios,statistics = [],[]
    for op,comparison in sorted({(p['operation'],p['comparison']) for p in pairs}):
        values = [p['ratio'] for p in pairs if (p['operation'],p['comparison'])==(op,comparison)]
        control = min(p['ratio'] for p in pairs if p['operation']==op and p['comparison']=='A-control/A')
        ratios.append(dict(operation=op,comparison=comparison,pairs=len(values),median=st.median(values),minimum=min(values),
            maximum=max(values),meets_target=st.median(values)<=.8 and max(values)<control if comparison=='B/A' else None))
    for op in selected:
        for name,arms in [('A',('Engine-A','Engine-control')),('B',('Engine-B','Engine-B-control')),('FeOs',('FeOs',)),('Clapeyron',('Clapeyron',))]:
            group = [r for r in blocks if r['operation']==op and r['arm'] in arms]
            wall = [r['wall_median_us'] for r in group]; cpu = [r['cpu_median_us'] for r in group if r['cpu_median_us'] is not None]
            statistics.append(dict(operation=op,arm=name,blocks=len(group),wall_median_us=st.median(wall),wall_min_us=min(wall),
                wall_max_us=max(wall),cpu_median_us=st.median(cpu) if cpu else None))
    for name,rows in [('selection.csv',admissions),('blocks.csv',blocks),('cost-value-checks.csv',checks),('paired-ratios.csv',pairs),('ratios.csv',ratios),('statistics.csv',statistics)]: write(name,rows)
    print(json.dumps(dict(ratios=ratios,statistics=statistics),indent=2))


if __name__=='__main__': main()

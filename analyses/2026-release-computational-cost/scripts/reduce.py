"""Reduce fresh, quiet, value-matched blocks; retain every exclusion and sample."""
import csv
import gzip
import json
from pathlib import Path
import statistics

from verify import compare, flat, values

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / 'results'
MAXIMUM_LOAD = json.loads((HERE / 'inputs/inputs.json').read_text())['quiet_host']['maximum_load1_exclusive']


def write(name, rows):
    assert rows, name
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with (OUT / name).open('w') as stream:
        writer = csv.DictWriter(stream, keys, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def quantile(v, p):
    v = sorted(v)
    at = (len(v) - 1) * p
    lo = int(at)
    return v[lo] + (v[min(lo + 1, len(v) - 1)] - v[lo]) * (at - lo)


def main():
    reference = json.loads((OUT / 'preflight-values.json').read_text())
    blocks, samples, host, checks = [], [], [], []
    def started(path):
        with gzip.open(path, 'rt') as stream:
            return json.loads(next(stream))['utc']
    for path in sorted((OUT / 'raw').glob('measurement-*.jsonl.gz'), key=started):
        for line in gzip.open(path, 'rt'):
            r = json.loads(line)
            if r['kind'] != 'block':
                continue
            operation, arm = r['operation'], r['arm']
            tool = 'Engine' if arm.startswith('Engine') else arm
            key = f'{path.stem}/{operation}/{r["repetition"]}/{arm}'
            quiet = all(s['load'][0] < MAXIMUM_LOAD and not s['compilers'] for s in (r['before'], r['after']))
            for side in ('before', 'after'):
                s = r[side]
                host.append(dict(block=key, operation=operation, repetition=r['repetition'], arm=arm, side=side, utc=s['utc'],
                                 load1=s['load'][0], load5=s['load'][1], load15=s['load'][2],
                                 compile_processes=json.dumps(s['compile_processes']), outside_fits=json.dumps(s['outside_fits']),
                                 disposition='quiet' if quiet else 'disturbed; excluded and retry required'))
            valid = True
            times = []
            try:
                for index, sample in enumerate(r['reply']['samples']):
                    result = values(operation, tool, sample)
                    a, b = flat(reference[operation][tool]), flat(result)
                    assert a.keys() == b.keys()
                    worst = max(abs(a[k] - b[k]) / max(1., abs(a[k]), abs(b[k])) for k in a)
                    passed = worst <= 1e-10
                    checks.append(dict(block=key, sample=index, maximum_scaled_difference=worst, tolerance=1e-10, passed=passed))
                    valid &= passed
                    seconds = sample['seconds'] / sample.get('per_call', 1)
                    assert seconds > 0
                    times.append(seconds)
                    samples.append(dict(block=key, operation=operation, repetition=r['repetition'], arm=arm, sample=index,
                                        seconds=seconds, cpu_seconds=sample.get('cpu_seconds'), quiet=quiet, value_matched=passed))
                assert times
            except Exception as error:
                valid = False
                checks.append(dict(block=key, passed=False, cause=r['reply'].get('error') or str(error)))
            blocks.append(dict(block=key, operation=operation, repetition=r['repetition'], arm=arm, quiet=quiet,
                               value_matched=valid, samples=len(times), median_seconds=statistics.median(times) if times else None))
    assert blocks, 'no measurement blocks'
    write('timed-value-checks.csv', checks)
    summary = []
    groups = sorted({(r['operation'], r['arm']) for r in blocks})
    for operation, arm in groups:
        accepted = [r for r in blocks if (r['operation'], r['arm']) == (operation, arm) and r['quiet'] and r['value_matched']]
        # A later successful repetition supersedes a disturbed one; never select by speed.
        selected = {r['repetition']: r for r in accepted}
        medians = [r['median_seconds'] for r in selected.values()]
        if len(medians) != 3:
            print('RETRY REQUIRED', operation, arm, 'quiet repetitions:', sorted(selected))
            continue
        selected_keys = {r['block'] for r in selected.values()}
        times = [s['seconds'] for s in samples if s['block'] in selected_keys]
        summary.append(dict(operation=operation, arm=arm, blocks=len(medians), samples=len(times),
                            median_seconds=statistics.median(medians), minimum_block_median_seconds=min(medians),
                            maximum_block_median_seconds=max(medians), sample_p10_seconds=quantile(times, .1),
                            sample_p90_seconds=quantile(times, .9), selected_blocks=json.dumps(sorted(selected_keys))))
    write('statistics.csv', summary)
    selected_keys = {key for r in summary for key in json.loads(r['selected_blocks'])}
    for r in blocks:
        r['selected'] = r['block'] in selected_keys
        r['disposition'] = 'selected quiet block' if r['selected'] else 'disturbed; excluded' if not r['quiet'] else 'superseded by complete retry'
    for r in samples:
        r['selected'] = r['block'] in selected_keys
    for r in host:
        if r['disposition'] == 'quiet':
            r['disposition'] = 'selected quiet block' if r['block'] in selected_keys else 'superseded by complete retry'
    write('host.csv', host)
    write('samples.csv', samples)
    write('blocks.csv', blocks)
    by_key = {(r['operation'], r['arm']): r for r in summary}
    ratios = []
    for operation, arm in groups:
        if arm not in ('FeOs', 'Clapeyron') or (operation, arm) not in by_key:
            continue
        if not all((operation, a) in by_key for a in ('Engine-A', 'Engine-control')):
            continue
        paired, control = [], []
        for rep in range(3):
            selected = {a: next(r for r in reversed(blocks) if r['operation'] == operation and r['arm'] == a and r['repetition'] == rep and r['quiet'] and r['value_matched'])
                        for a in ('Engine-A', 'Engine-control', arm)}
            for a in ('Engine-A', 'Engine-control'):
                paired.append(selected[a]['median_seconds'] / selected[arm]['median_seconds'])
            control.append(selected['Engine-A']['median_seconds'] / selected['Engine-control']['median_seconds'])
        ratios.append(dict(operation=operation, comparator=arm, package_over_comparator=statistics.median(paired),
                           minimum_paired_ratio=min(paired), maximum_paired_ratio=max(paired), paired_ratios=json.dumps(paired),
                           minimum_control_ratio=min(control), maximum_control_ratio=max(control)))
    for case in ('C2-0.5', 'C4-298.15'):
        exact, fd = f'{case}/actions-exact', f'{case}/actions-fd'
        if not all((op, a) in by_key for op in (exact, fd) for a in ('Engine-A', 'Engine-control')):
            continue
        paired = []
        for rep in range(3):
            for arm in ('Engine-A', 'Engine-control'):
                row = {op: next(r for r in reversed(blocks) if r['operation'] == op and r['arm'] == arm and r['repetition'] == rep and r['quiet'] and r['value_matched']) for op in (exact, fd)}
                paired.append(row[exact]['median_seconds'] / row[fd]['median_seconds'])
        ratios.append(dict(operation=f'{case}/actions-on-solved-state', comparator='package central differences',
                           package_over_comparator=statistics.median(paired), minimum_paired_ratio=min(paired), maximum_paired_ratio=max(paired), paired_ratios=json.dumps(paired)))
    write('ratios.csv', ratios)
    print('quiet value-matched blocks:', sum(r['quiet'] and r['value_matched'] for r in blocks), 'of', len(blocks))
    for r in ratios:
        print(r['operation'], r['comparator'], r['package_over_comparator'], (r['minimum_paired_ratio'], r['maximum_paired_ratio']))


if __name__ == '__main__':
    main()

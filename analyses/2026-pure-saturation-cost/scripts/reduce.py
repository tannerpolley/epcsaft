"""Read every retained cost block; reject disturbed operations and compare all returned values."""
import csv
import gzip
import json
import math
import os
from pathlib import Path
import statistics as st
from epcsaft.equilibrium import (ComparisonEligibility, ComparisonMode, NumericalEstimate, Quantity,
                                  compare_equilibrium_quantity)

HERE = Path(os.environ.get('EPCSAFT_COST_OUTPUT', str(Path(__file__).resolve().parents[1])))


def write(name, rows):
    with (HERE / 'results' / name).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def main():
    root = HERE / 'results'
    records, selected = [], {}
    estimate, eligibility = NumericalEstimate(), ComparisonEligibility()
    for path in sorted((root / 'raw').glob('measurement-*.jsonl.gz'), key=lambda p: int(p.name.split('-')[1].split('.')[0])):
        rows = [json.loads(line) for line in gzip.open(path, 'rt')]
        blocks = [r for r in rows if r['kind'] == 'block']
        for operation in {r['operation'] for r in blocks}:
            part = [r for r in blocks if r['operation'] == operation]
            quiet = all(s['load'][0] < 4 and not s['compilers'] for r in part for s in (r['before'], r['after']))
            records.append(dict(path=str(path.relative_to(HERE)), operation=operation, quiet=quiet, blocks=len(part)))
            if quiet:
                selected[operation] = part
            else:
                selected.pop(operation, None)
    assert set(selected) == {'ordinary', 'prepared'}
    pre = json.loads((root / 'preflight.json').read_text())
    ref = pre['operations']['ordinary']['Engine-A']['samples'][0]['value']
    reference = [ref['P'], *[p['rho'] for p in ref['phases']]]
    block_rows, checks, times = [], [], {}
    for operation, rows in selected.items():
        assert len(rows) == 18
        for row in rows:
            reply = row['reply']; assert 'error' not in reply, reply
            samples = reply['samples']; assert len(samples) == 25
            arm, rep = row['arm'], row['repetition']
            walls = [s['seconds'] * 1e6 for s in samples]
            cpus = [s['cpu_seconds'] * 1e6 for s in samples if 'cpu_seconds' in s]
            median = st.median(walls)
            times[operation, arm, rep] = median
            block_rows.append(dict(operation=operation, arm=arm, repetition=rep, wall_median_us=median,
                wall_min_us=min(walls), wall_max_us=max(walls), cpu_median_us=st.median(cpus) if cpus else None,
                prepare_us=reply.get('prepare_seconds') * 1e6 if reply.get('prepare_seconds') is not None else None,
                load_before=row['before']['load'][0], load_after=row['after']['load'][0]))
            for sample_index, sample in enumerate(samples):
                value = sample['value']; actual = [value['P'], *(value['densities'] if 'densities' in value else [p['rho'] for p in value['phases']])]
                for quantity, av, bv in zip(('pressure_Pa', 'liquid_density_mol_m3', 'vapor_density_mol_m3'), reference, actual):
                    error = abs(av-bv)/max(abs(av), abs(bv))
                    limit = 1e-10 if arm.startswith('Engine') else 1e-7
                    comparison = compare_equilibrium_quantity(
                        Quantity.Linear if quantity == 'pressure_Pa' else Quantity.Positive,
                        av, bv, limit * max(abs(av), abs(bv)), estimate, estimate, eligibility, ComparisonMode.Strict)
                    checks.append(dict(operation=operation, arm=arm, repetition=rep, sample=sample_index,
                        quantity=quantity, A=av, value=bv, relative_error=error, limit=limit, judgment=comparison.judgment, cause=comparison.cause, passed=comparison.judgment != 'failed'))
    for row in selected['prepared']:
        for sample in row['reply']['samples']:
            value = sample['value']
            if row['arm'].startswith('Engine'):
                expected = pre['operations']['prepared']['Engine-A']['samples'][0]['value']
                assert value == expected, (row['arm'], value, expected)
            if 'compositions' in value:
                assert value['compositions'] == [[1.0], [1.0]]
            for phase, expected in zip(value.get('phases', []), ref['phases']):
                assert phase['phase'] == expected['phase'] and phase['x'] == expected['x']
                pressure = phase.get('P', phase.get('pressure'))
                if pressure is not None:
                    target = expected.get('P', expected.get('pressure', ref['P']))
                    limit = 1e-10 if row['arm'].startswith('Engine') else 1e-7
                    assert abs(pressure-target) <= limit * max(abs(pressure), abs(target))
    assert all(r['passed'] for r in checks)
    ratio_rows = []
    for operation in selected:
        for comparator in ('FeOs', 'Clapeyron'):
            for rep in range(3):
                for arm in ('Engine-B', 'Engine-B-control'):
                    ratio_rows.append(dict(operation=operation, comparison='B/' + comparator, repetition=rep, arm=arm,
                        ratio=times[operation, arm, rep]/times[operation, comparator, rep]))
        for rep in range(3):
            ratio_rows.append(dict(operation=operation, comparison='B-control/B', repetition=rep, arm='Engine-B-control',
                ratio=times[operation, 'Engine-B-control', rep]/times[operation, 'Engine-B', rep]))
            ratio_rows.append(dict(operation=operation, comparison='A-control/A', repetition=rep, arm='Engine-control',
                ratio=times[operation, 'Engine-control', rep]/times[operation, 'Engine-A', rep]))
            for arm, baseline in (('Engine-B', 'Engine-A'), ('Engine-B-control', 'Engine-control')):
                ratio_rows.append(dict(operation=operation, comparison='B/A', repetition=rep, arm=arm,
                    ratio=times[operation, arm, rep]/times[operation, baseline, rep]))
    summary = []
    for operation, comparison in sorted({(r['operation'], r['comparison']) for r in ratio_rows}):
        values = [r['ratio'] for r in ratio_rows if (r['operation'], r['comparison']) == (operation, comparison)]
        summary.append(dict(operation=operation, comparison=comparison, pairs=len(values), median=st.median(values),
            minimum=min(values), maximum=max(values), target=None))
    statistics = []
    for operation in selected:
        for name, arms in [('A', ('Engine-A','Engine-control')), ('B', ('Engine-B','Engine-B-control')), ('FeOs', ('FeOs',)), ('Clapeyron', ('Clapeyron',))]:
            group = [r for r in block_rows if r['operation'] == operation and r['arm'] in arms]
            if not group: continue
            walls = [r['wall_median_us'] for r in group]; cpus = [r['cpu_median_us'] for r in group if r['cpu_median_us'] is not None]
            statistics.append(dict(operation=operation, arm=name, blocks=len(group), wall_median_us=st.median(walls),
                wall_min_us=min(walls), wall_max_us=max(walls), cpu_median_us=st.median(cpus) if cpus else None))
    verdicts = []
    for operation in selected:
        ratios = [r['ratio'] for r in ratio_rows if r['operation'] == operation and r['comparison'] == 'B/A']
        controls = [abs(math.log(r['ratio'])) for r in ratio_rows if r['operation'] == operation
                    and r['comparison'] in ('A-control/A', 'B-control/B')]
        savings = [times[operation, baseline, rep] - times[operation, arm, rep]
                   for rep in range(3) for arm, baseline in (('Engine-B', 'Engine-A'), ('Engine-B-control', 'Engine-control'))]
        R, C = st.median(ratios), max(controls)
        b = next(r for r in statistics if r['operation'] == operation and r['arm'] == 'B')
        feos = [r['ratio'] for r in ratio_rows if r['operation'] == operation and r['comparison'] == 'B/FeOs']
        verdicts.append(dict(operation=operation, R=R, C=C, negative_log_R=-math.log(R), saving_resolved=-math.log(R)>C,
            paired_savings_us=savings, saving_median_us=st.median(savings), saving_range_us=[min(savings),max(savings)],
            below_historical_feos=b['wall_median_us']<24.5, below_fresh_feos=st.median(feos)<1))
    (root/'saving.json').write_text(json.dumps(verdicts,indent=2)+'\n')
    for name, rows in [('selection.csv',records), ('blocks.csv',block_rows), ('checks.csv',checks), ('paired-ratios.csv',ratio_rows), ('ratios.csv',summary), ('statistics.csv',statistics)]:
        write(name,rows)
    print(json.dumps(dict(ratios=summary,absolute_times=statistics,verdicts=verdicts),indent=2))


if __name__ == '__main__':
    main()

"""Reduce retained timed returns; no Engine imports or model execution."""
from collections import defaultdict
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics
import runpy
import zipfile

HERE = Path(__file__).resolve().parents[1]
R = HERE / 'results'
LABELS = {
    'case-6/figure-2-point-7/cold': 'Ionic LLE, cold',
    'case-6/figure-2-point-7/accepted': 'Ionic LLE, accepted coordinates',
    'case-7/flash-5918.5286-Pa/cold': 'MEA flash, cold',
    'case-7/flash-5918.5286-Pa/accepted': 'MEA flash, accepted coordinates',
    'khudaida/sequential-8': 'Khudaida, eight-point sequence',
    'mea-nominally-unloaded/353.15-K': 'Nominally unloaded MEA',
    'C1-300': 'Propane saturation, 300 K',
    'C2-0.5': 'Ethanol/butane bubble, x = 0.5',
    'methane-bubble': 'Methane/ethane bubble, 180 K',
    'C4-298.15': 'Water/pentanol LLE, 298.15 K',
    'C4-323.15': 'Water/pentanol LLE, 323.15 K',
    'methane-flash': 'Methane/ethane flash, 180 K',
    'A-propane-liquid': 'Propane liquid density, per call',
    'mea-f1/residual-jacobian': 'MEA F1 residual + Jacobian, serial',
}
EXTENSION = json.loads((HERE / 'inputs/extension.json').read_text())
for op, q in EXTENSION.items():
    LABELS[op] = (f"Methane/ethane Discover, {q['T_K']:g} K, {q['P_Pa']/1e6:g} MPa" if q['kind'] == 'discover'
                  else f"Propane saturation, {q['T_K']:g} K, cold" if q['kind'] == 'saturation'
                  else f"Methane/ethane {q['kind']}, {q['T_K']:g} K")
PHYSICAL = ('temperature', 'pressure', 'amounts', 'mole_fractions', 'molar_densities',
            'predictions', 'gradient', 'optimizer_gradient', 'residual_optimizer_gradient')


def read(name):
    return list(csv.DictReader((R / name).open()))


def write(name, rows, fields=None):
    opener = gzip.open if name.endswith('.gz') else open
    with opener(R / name, 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    campaign = json.loads((R / 'campaign.json').read_text())
    current = campaign['current']
    ma57 = current + '-ma57'
    controls = campaign['controls']
    library_checks = 0
    for label, identity in campaign['identities'].items():
        version = controls.get(label, label)
        env = next(v for v in campaign['environments'] if v['version'] == version)
        with zipfile.ZipFile(env['wheel']) as wheel:
            for path, loaded_hash in identity['loaded_libraries'].items():
                member = next(n for n in wheel.namelist() if Path(n).name == Path(path).name)
                with wheel.open(member) as stream:
                    assert hashlib.file_digest(stream, 'sha256').hexdigest() == loaded_hash, (label, path)
                library_checks += 1
    versions = [v['version'] for v in campaign['environments']]
    operations = list(LABELS)
    timings = read('timings.csv')
    by_block, by_version = defaultdict(list), defaultdict(list)
    for row in timings:
        value = float(row['wall_seconds'])
        assert math.isfinite(value) and value > 0
        by_block[(row['version'], row['operation'], int(row['block']))].append(value)
    for (version, operation, block), values in by_block.items():
        by_version[(version, operation)].append(statistics.median(values))
    statistics_rows = []
    for (version, operation), medians in by_version.items():
        statistics_rows.append(dict(version=version, operation=operation, blocks=len(medians),
            median_seconds=statistics.median(medians), minimum_seconds=min(medians), maximum_seconds=max(medians)))
    write('statistics.csv', statistics_rows)
    competitors = {(r['operation'], r['library']): float(r['seconds']) for r in read('competitors.csv')}
    headline = []
    for operation in operations:
        def median(version):
            return statistics.median(by_version[(version, operation)]) if (version, operation) in by_version else None
        earliest, morning, today = map(median, ('1218dc78', '5ad73fb9', current))
        aa = by_version[('control', operation)] + by_version[(current, operation)]
        aa_ma57 = by_version[('control-ma57', operation)] + by_version[(ma57, operation)]
        headline.append(dict(operation=operation, label=LABELS[operation], earliest_seconds=earliest,
            morning_seconds=morning, last_seconds=median('cf26a17b'), oct16_seconds=median('21e761c3'),
            oct17_seconds=median('cc63a135'), today_seconds=today, ma57_seconds=median(ma57),
            earliest_over_today=earliest/today if earliest and today else None,
            morning_over_today=morning/today if morning and today else None,
            today_over_clapeyron=today/competitors[(operation, 'Clapeyron')] if today and (operation, 'Clapeyron') in competitors else None,
            today_over_feos=today/competitors[(operation, 'FeOs')] if today and (operation, 'FeOs') in competitors else None,
            aa_min_seconds=min(aa) if aa else None, aa_max_seconds=max(aa) if aa else None,
            aa_range_over_median=(max(aa)-min(aa))/statistics.median(aa) if aa else None,
            ma57_aa_min_seconds=min(aa_ma57) if aa_ma57 else None,
            ma57_aa_max_seconds=max(aa_ma57) if aa_ma57 else None))
    write('headline.csv', headline)

    references = {(r['version'], r['operation'], r['point'], r['quantity'], int(r['index'])): r
                  for r in csv.DictReader(gzip.open(R / 'retained-values.csv.gz', 'rt'))}
    # The harness owns F1 quantity scales. Reading its frozen inputs here performs no calculation.
    f1 = json.load(gzip.open(HERE.parents[1] / 'tools/fixtures/mea-f1.json.gz', 'rt'))
    key_harness = runpy.run_path(str(HERE.parents[1] / 'tools/time-key-solves'), run_name='history_comparison')
    def differences(actual, expected, operation, quantity):
        if operation.startswith(('case-', 'khudaida/', 'mea-nominally', 'mea-f1/')):
            return key_harness['scaled_differences'](actual, expected,
                case=key_harness['F1'] if operation.startswith('mea-f1/') else None, quantity=quantity)
        return [abs(a-b)/max(abs(a), abs(b), 1.) for a, b in zip(actual, expected, strict=True)]
    first, aggregate, fields = {}, {}, []
    diagnostics, density, invariants = [], [], []
    work = defaultdict(set)
    with gzip.open(R / 'returns.jsonl.gz', 'rt') as handle:
        for line in handle:
            record = json.loads(line)
            if 'unavailable' in record:
                continue
            version = controls.get(record['version'], record['version'])
            operation = record['operation']
            key = (version, operation)
            a = aggregate.setdefault(key, dict(version=version, operation=operation, checked_fields=0,
                differing_fields=0, maximum_scaled_difference=0., repeat_maximum_scaled_difference=0.,
                unreferenced_fields=0, samples=0, maximum_raw_residual=0., maximum_classified_residual=0.))
            for sample, value in enumerate(record['values']):
                a['samples'] += 1
                if 'missing_work_diagnostics' in value:
                    diagnostics.append(dict(block=record['block'], version=record['version'],
                        missing=','.join(value['missing_work_diagnostics']),
                        **{q: value['rows'][0].get(q, '') for q in ('reference_computations', 'reference_columns', 'density_solves', 'equilibrium_iterations')}))
                for row in value['rows']:
                    counters = key_harness['F1_COUNTS'] if operation.startswith('mea-f1/') else key_harness['COUNTS']
                    if all(q in row for q in counters):
                        work[(version, operation, row['point'])].add(tuple(row[q] for q in counters))
                    for q in ('raw_residual', 'classified_residual'):
                        a['maximum_' + q] = max(a['maximum_' + q], row.get(q, 0.))
                    if 'mole_fractions' in row:
                        n = len(row['mole_fractions']) // len(row['amounts'])
                        normalization = max(abs(sum(row['mole_fractions'][i:i+n])-1) for i in range(0, len(row['mole_fractions']), n))
                        assert normalization <= 1e-10 and min(row['amounts']) >= 0
                        invariants.append(dict(block=record['block'], sample=sample, version=record['version'],
                            operation=operation, point=row['point'], composition_sum_error=normalization,
                            minimum_amount=min(row['amounts'])))
                    if 'density_oracle_relative_error' in row:
                        error = row['density_oracle_relative_error']
                        assert error <= 1e-10
                        density.append(dict(block=record['block'], sample=sample, version=record['version'],
                            point=row['point'], relative_error=error))
                    for quantity in PHYSICAL:
                        actual_values = row.get(quantity, [])
                        keys = [(version, operation, row['point'], quantity, i) for i in range(len(actual_values))]
                        reference_keys = [(current if version == ma57 else version, *k[1:]) for k in keys]
                        reference_values = [float(references[k]['value']) if k in references else actual
                                            for k, actual in zip(reference_keys, actual_values)]
                        repeats = [first.setdefault(k, actual) for k, actual in zip(keys, actual_values)]
                        deltas = differences(actual_values, reference_values, operation, quantity)
                        repeat_deltas = differences(actual_values, repeats, operation, quantity)
                        for i, actual in enumerate(actual_values):
                            assert math.isfinite(actual)
                            a['repeat_maximum_scaled_difference'] = max(a['repeat_maximum_scaled_difference'],
                                repeat_deltas[i])
                            ref = references.get(reference_keys[i])
                            if ref is None:
                                a['unreferenced_fields'] += 1
                                continue
                            expected = float(ref['value'])
                            delta = deltas[i]
                            a['checked_fields'] += 1
                            a['differing_fields'] += delta > float(ref['tolerance'])
                            a['maximum_scaled_difference'] = max(a['maximum_scaled_difference'], delta)
                            fields.append(dict(block=record['block'], sample=sample, version=record['version'],
                                operation=operation, point=row['point'], quantity=quantity, index=i,
                                actual=actual, reference=expected, scaled_difference=delta,
                                tolerance=ref['tolerance'], within_tolerance=delta <= float(ref['tolerance'])))
    write('return-checks.csv.gz', fields)
    write('return-summary.csv', list(aggregate.values()))
    write('f1-work-diagnostics.csv', diagnostics)
    write('density-checks.csv.gz', density)
    write('composition-checks.csv.gz', invariants)
    changes = []
    first_vectors = defaultdict(dict)
    for (v, op, point, quantity, i), actual in first.items():
        first_vectors[(v, op, point, quantity)][i] = actual
    for version in versions:
        for operation in operations:
            ds = []
            for (v, op, point, quantity), actual in first_vectors.items():
                if v != version or op != operation:
                    continue
                current_values = first_vectors.get((current, op, point, quantity))
                if current_values is not None:
                    ds.extend(differences([actual[i] for i in sorted(actual)],
                        [current_values[i] for i in sorted(current_values)], operation, quantity))
            if ds:
                changes.append(dict(version=version, operation=operation, compared_fields=len(ds),
                                    maximum_scaled_difference_to_today=max(ds)))
    write('state-changes.csv', changes)
    solver_difference = max(r['maximum_scaled_difference_to_today'] for r in changes if r['version'] == ma57)
    availability = {(r['version'], r['operation'], int(r['repetition'])) for r in read('unavailable.csv')}
    for v in versions:
        for op in operations:
            found = by_version.get((v, op), [])
            assert len(found) == 3 or all((v, op, rep) in availability for rep in range(3)), (v, op, len(found))
    assert all(len(by_version[(control, op)]) == 6 for control in controls for op in operations)
    work_checks = []
    for left, right in (('cc63a135', current), (current, ma57)):
        for (v, op, point), counts in work.items():
            if v != left or not op.startswith(('case-', 'khudaida/', 'mea-nominally', 'mea-f1/')):
                continue
            other = work.get((right, op, point), set())
            work_checks.append(dict(left=left, right=right, operation=op, point=point,
                left_counts=json.dumps(sorted(counts)), right_counts=json.dumps(sorted(other)),
                equal=len(counts) == 1 and counts == other))
    assert len(work_checks) == 28, 'missing key-solve work observations'
    write('work-count-checks.csv', work_checks)
    comparisons = []
    for h in headline:
        for label, field in (('1218dc78 to current MUMPS', 'earliest_seconds'),
                             ('5ad73fb9 to current MUMPS', 'morning_seconds'),
                             ('cf26a17b to current MUMPS', 'last_seconds'),
                             ('cc63a135 to current MUMPS', 'oct17_seconds'),
                             ('current MUMPS to MA57', 'today_seconds')):
            a = h[field]
            b = h['ma57_seconds'] if label.endswith('MA57') else h['today_seconds']
            width = h['aa_max_seconds'] - h['aa_min_seconds']
            if label.endswith('MA57'):
                width = max(width, h['ma57_aa_max_seconds'] - h['ma57_aa_min_seconds'])
            comparisons.append(dict(operation=h['operation'], comparison=label, earlier_seconds=a, later_seconds=b,
                difference_seconds=a-b if a and b else None, aa_range_seconds=width,
                beyond_aa_range=abs(a-b)>width if a and b else None))
    write('timing-comparisons.csv', comparisons)
    summary = dict(measured_samples=len(timings), measured_blocks=len(by_block),
        loaded_library_checks=library_checks,
        repeat_tolerance=1e-10, repeats_within_tolerance=all(a['repeat_maximum_scaled_difference'] <= 1e-10 for a in aggregate.values()),
        maximum_solver_arm_difference=solver_difference, solver_arms_within_tolerance=solver_difference <= 1e-10,
        work_count_comparisons=len(work_checks), work_counts_equal=all(r['equal'] for r in work_checks),
        work_count_differences=sum(not r['equal'] for r in work_checks), f1_observations=len(f1['observations']),
        f1_coordinates=len(f1['coordinates']),
        unavailable_blocks=len(availability), prior_value_checks=len(fields),
        differing_fields=sum(a['differing_fields'] for a in aggregate.values()),
        maximum_repeat_difference=max(a['repeat_maximum_scaled_difference'] for a in aggregate.values()),
        maximum_density_relative_error=max(r['relative_error'] for r in density),
        maximum_composition_sum_error=max(r['composition_sum_error'] for r in invariants),
        load_min=min(float(r['load_before']) for r in timings),
        load_max=max(max(float(r['load_before']), float(r['load_after'])) for r in timings),
        gates=sorted({r['gate'] for r in timings}),
        declaration_sha256={v['version']: hashlib.sha256(
            (Path(v['source'])/'tools/time-key-solves').read_bytes()).hexdigest() for v in campaign['environments']})
    (R / 'verification.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()

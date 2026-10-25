"""Measure native work on frozen Engine fitting and Khudaida inputs."""
import csv
import hashlib
import importlib.util
import json
import os
import resource
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
STUDY = Path(__file__).resolve().parents[1]


def write_csv(path, header, rows):
    with path.open('w', newline='') as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def profile_rows(trace, output):
    spans = trace.spans
    exclusive = [s.wall_us for s in spans]
    totals = defaultdict(lambda: [0, 0.0, 0.0])
    for s in spans:
        if s.parent is not None:
            exclusive[s.parent] -= s.wall_us
    for s, own in zip(spans, exclusive, strict=True):
        row = totals[s.name]
        row[0] += 1
        row[1] += s.wall_us / 1e6
        row[2] += own / 1e6
    write_csv(output, ['span', 'calls', 'inclusive_s', 'exclusive_s'],
              [[name] + values for name, values in sorted(totals.items())])


def equilibrium_run(output):
    import epcsaft
    from epcsaft import equilibrium as eq
    source = ROOT / 'analyses/2026-khudaida-idaes-settler/scripts/run_fixed_settler.py'
    spec = importlib.util.spec_from_file_location('fixed_settler', source)
    case = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = case
    spec.loader.exec_module(case)
    model = epcsaft.Mixture(epcsaft.Parameters.from_json(case.PARAMETERS, components=case.COMPONENTS))
    rows, previous = [], None
    load_before = os.getloadavg()
    usage_before = resource.getrusage(resource.RUSAGE_SELF)
    wall, cpu = time.perf_counter(), time.process_time()
    with epcsaft.trace() as trace:
        for point, feed in sorted(case.feed_rows().items()):
            result = eq.solve_equilibrium(model, case.engine_problem(feed), previous)
            if not result.success or not result.requested_tolerance_met:
                raise RuntimeError(f'Khudaida point {point}: {result.message}')
            for field in ('amounts', 'mole_fractions', 'molar_densities', 'packing_fractions',
                          'residuals', 'coordinates', 'multipliers'):
                for index, value in enumerate(getattr(result, field)):
                    rows.append([point, field, index, value])
            compiled = eq.compile_problem(model, case.engine_problem(feed))
            observables = [eq.SolvedStateObservable(eq.SolvedStateObservableKind.PhaseMoleFraction, phase, i)
                           for phase in range(2) for i in range(5)]
            direction = eq.SolvedStateActionDirection(1.0, 0.0, [0.0] * 5, [])
            for order in (1, 2):
                batch = eq.solved_state_actions(compiled, result,
                    eq.SolvedStateActionRequest(observables, [direction] * order))
                for index, item in enumerate(batch.results):
                    if item.status != eq.SolvedStateActionStatus.Available:
                        raise RuntimeError(f'Khudaida A{order}: {item.status}')
                    rows.append([point, f'mole_fraction_temperature_A{order}', index, item.action])
            previous = result
    usage_after = resource.getrusage(resource.RUSAGE_SELF)
    summary = dict(wall_s=time.perf_counter()-wall, cpu_s=time.process_time()-cpu, points=8,
        user_s=usage_after.ru_utime-usage_before.ru_utime, system_s=usage_after.ru_stime-usage_before.ru_stime,
        load_before=load_before, load_after=os.getloadavg())
    write_csv(output.with_suffix('.csv'), ['point', 'quantity', 'index', 'value'], rows)
    profile_rows(trace, Path(os.environ.get('EPCSAFT_REPLAY_PROFILE', str(output.with_suffix('.profile.csv')))))
    output.write_text(json.dumps(summary)+'\n')
    print(json.dumps(summary), flush=True)


def repeat(workload, label, wheel):
    subprocess.run(['uv', 'pip', 'install', '--python', sys.executable, '--no-deps', str(wheel)], check=True)
    subprocess.run([sys.executable, str(ROOT / 'tools/verify-installed-wheel.py'),
                    '--wheel', str(wheel)], check=True, cwd=ROOT)
    output = STUDY / 'results'
    timing_rows, span_rows, summaries = [], [], []
    wheel_hash = hashlib.sha256(wheel.read_bytes()).hexdigest()
    with tempfile.TemporaryDirectory(prefix='epcsaft-129-') as temporary:
        reference = None
        for run in range(1, 4):
            result = Path(temporary) / f'{workload}-{run}.json'
            profile = result.with_suffix('.profile.csv')
            env = dict(os.environ, EPCSAFT_REPLAY_WHEEL=str(wheel), EPCSAFT_REPLAY_PROFILE=str(profile))
            command = ([sys.executable, str(STUDY / 'scripts/replay.py'), 'initial', str(result)]
                       if workload == 'fit' else [sys.executable, __file__, 'equilibrium-once', str(result)])
            subprocess.run(command, check=True, env=env, cwd=ROOT)
            values = result.with_suffix('.csv').read_bytes()
            if reference is None:
                reference = values
                values_hash = hashlib.sha256(values).hexdigest()
                values_file = next((p for p in output.glob('*.csv')
                    if p.read_bytes().splitlines() == values.splitlines()),
                    output / f'{label}-{workload}-values.csv')
                if not values_file.exists():
                    values_file.write_bytes(values)
                values_hash = hashlib.sha256(values_file.read_bytes()).hexdigest()
            elif reference != values:
                raise RuntimeError('Numerical rows changed between timing repeats')
            summary = json.loads(result.read_text())
            summaries.append(summary)
            timing_rows.append([label, workload, run, wheel_hash, summary['wall_s'], summary['cpu_s']])
            with profile.open() as handle:
                span_rows.extend([[label, workload, run] + row for row in list(csv.reader(handle))[1:]])
    (output / f'{label}-{workload}-summary.json').write_text(json.dumps(dict(
        wheel_sha256=wheel_hash, values_file=str(values_file.relative_to(STUDY)),
        values_sha256=values_hash, repeats_identical=True,
        calculation={k: v for k, v in summaries[0].items() if k not in ('wall_s', 'cpu_s', 'message', 'wheel_sha256')},
        termination_messages=[s.get('message', 'eight accepted points and available actions') for s in summaries]), indent=2)+'\n')
    write_csv(output / f'{label}-{workload}-timing.csv',
              ['wheel', 'workload', 'repeat', 'wheel_sha256', 'wall_s', 'cpu_s'], timing_rows)
    write_csv(output / f'{label}-{workload}-profile.csv',
              ['wheel', 'workload', 'repeat', 'span', 'calls', 'inclusive_s', 'exclusive_s'], span_rows)


def paired(workload, pair, baseline, candidate):
    output = STUDY / 'results'
    timings, profiles, records = [], [], []
    with tempfile.TemporaryDirectory(prefix='epcsaft-129-paired-') as temporary:
        for label, wheel in [('baseline', baseline), ('candidate', candidate)]:
            subprocess.run(['uv', 'pip', 'install', '--python', sys.executable, '--no-deps', str(wheel)], check=True)
            subprocess.run([sys.executable, str(ROOT / 'tools/verify-installed-wheel.py'),
                            '--wheel', str(wheel)], check=True, cwd=ROOT)
            result = Path(temporary) / f'{label}.json'
            profile = result.with_suffix('.profile.csv')
            env = dict(os.environ, EPCSAFT_REPLAY_WHEEL=str(wheel), EPCSAFT_REPLAY_PROFILE=str(profile))
            command = ([sys.executable, str(STUDY / 'scripts/replay.py'), 'initial', str(result)]
                       if workload == 'fit' else [sys.executable, __file__, 'equilibrium-once', str(result)])
            subprocess.run(command, check=True, env=env, cwd=ROOT)
            values = result.with_suffix('.csv').read_bytes()
            reference = output / ('corrected-initial.csv' if workload == 'fit' else 'before-equilibrium-values.csv')
            if values.splitlines() != reference.read_bytes().splitlines():
                raise RuntimeError(f'{label} numerical values disagree with immutable baseline')
            summary = json.loads(result.read_text())
            cpu = summary['user_s'] + summary['system_s']
            timings.append([pair, label, hashlib.sha256(wheel.read_bytes()).hexdigest(),
                summary['wall_s'], cpu, summary['user_s'], summary['system_s'],
                *summary['load_before'], *summary['load_after']])
            with profile.open() as handle:
                profiles.extend([[pair, label] + row for row in list(csv.reader(handle))[1:]])
            records.append(dict(pair=pair, wheel=label, numerical_rows_identical=True,
                values_file=str(reference.relative_to(STUDY)),
                termination=summary.get('message', 'eight accepted states and available actions')))
    write_csv(output / f'paired-{workload}-{pair}-timing.csv',
        ['pair', 'wheel', 'wheel_sha256', 'wall_s', 'cpu_s', 'user_s', 'system_s',
         'load_before_1m', 'load_before_5m', 'load_before_15m',
         'load_after_1m', 'load_after_5m', 'load_after_15m'], timings)
    write_csv(output / f'paired-{workload}-{pair}-profile.csv',
        ['pair', 'wheel', 'span', 'calls', 'inclusive_s', 'exclusive_s'], profiles)
    (output / f'paired-{workload}-{pair}-summary.json').write_text(json.dumps(records, indent=2)+'\n')


if __name__ == '__main__':
    if sys.argv[1] == 'paired':
        paired(sys.argv[2], int(sys.argv[3]), Path(sys.argv[4]).resolve(), Path(sys.argv[5]).resolve())
    elif sys.argv[1] == 'equilibrium-once':
        equilibrium_run(Path(sys.argv[2]))
    else:
        repeat(sys.argv[1], sys.argv[2], Path(sys.argv[3]).resolve())

"""One serial interleaved campaign on retained wheels, with current A=A controls."""
import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import selectors
import subprocess
import time

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
BUILD = ROOT / 'build/issue-549'
OPERATIONS = ('case-6/figure-2-point-7/cold', 'case-6/figure-2-point-7/accepted',
              'case-7/flash-5918.5286-Pa/cold', 'case-7/flash-5918.5286-Pa/accepted',
              'khudaida/sequential-8', 'mea-nominally-unloaded/353.15-K',
              'C1-300', 'C2-0.5', 'methane-bubble', 'C4-298.15', 'C4-323.15',
              'methane-flash', 'A-propane-liquid', 'mea-f1/residual-jacobian')
THREADS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
           'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')
CURRENT = '1cd3ab15'
HSL = Path.home() / '.cache/epcsaft-hsl/coinhsl-lib/libhsl.so'


def utc():
    return datetime.now(timezone.utc).isoformat()


def snapshot():
    compilers = subprocess.check_output(['ps', '-eo', 'comm='], text=True).splitlines()
    return dict(utc=utc(), load=list(os.getloadavg()),
                compilers=sorted({c.strip() for c in compilers if c.strip() in
                                  ('cc1plus', 'cc1', 'clang', 'clang++', 'gcc', 'g++', 'ninja')}))


def gate():
    begun = time.monotonic()
    while True:
        s = snapshot()
        # Prefer quiet conditions whenever they are already available; otherwise use the authorized relaxed gate.
        if s['load'][0] < 12 and not s['compilers']:
            return s | dict(mode='quiet' if s['load'][0] < 2 else 'relaxed', waited_seconds=time.monotonic()-begun)
        if time.monotonic()-begun > 3600:
            raise RuntimeError('load/compiler gate unavailable for one hour')
        print('waiting for load < 12 and no compilers:', s, flush=True)
        time.sleep(30)


def receive(process, timeout=600):
    with selectors.DefaultSelector() as selector:
        selector.register(process.stdout, selectors.EVENT_READ)
        if not selector.select(timeout):
            raise TimeoutError(f'worker {process.pid} exceeded {timeout} seconds')
    line = process.stdout.readline()
    if not line:
        raise RuntimeError(f'worker {process.pid} exited {process.poll()}')
    return json.loads(line)


def loaded_libraries(pid):
    paths = {line.split()[-1] for line in Path(f'/proc/{pid}/maps').read_text().splitlines()
             if 'libepcsaft_' in line or '/epcsaft/_core.' in line}
    return {path: hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in sorted(paths)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cpu', type=int, default=31)
    args = parser.parse_args()
    assert args.cpu in os.sched_getaffinity(0)
    environments = json.loads((BUILD / 'environments.json').read_text())
    current = next(v for v in environments if v['version'] == CURRENT)
    assert HSL.is_file(), HSL
    environments.append(current | dict(version=CURRENT + '-ma57'))
    controls = {'control': CURRENT, 'control-ma57': CURRENT + '-ma57'}
    operations = (*OPERATIONS, *json.loads((HERE / 'inputs/extension.json').read_text()))
    environment = os.environ.copy()
    environment.update({k: '1' for k in THREADS})
    workers, logs, identities = {}, {}, {}
    raw_path = HERE / 'results/returns.jsonl.gz'
    if raw_path.exists():
        raise FileExistsError(f'refusing to overwrite completed or partial campaign: {raw_path}')
    started = utc()
    timing_rows, failures, blocks = [], [], []
    try:
        for v in [*environments, current | dict(version='control'), current | dict(version='control-ma57')]:
            label = v['version']
            worker_environment = environment.copy()
            worker_environment.pop('EPCSAFT_HSL_LIBRARY', None)
            if label.endswith('-ma57'):
                worker_environment['EPCSAFT_HSL_LIBRARY'] = str(HSL)
            logs[label] = (BUILD / f'{label}.stderr.log').open('w')
            p = subprocess.Popen(['taskset', '-c', str(args.cpu), v['python'], '-u',
                                  str(HERE / 'scripts/worker.py'), v['source']],
                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=logs[label],
                                 text=True, env=worker_environment)
            workers[label] = p
            identities[label] = receive(p)
            if label in (CURRENT, 'control'):
                assert identities[label]['linear_solver'] == 'mumps', identities[label]
            elif label.endswith('-ma57'):
                assert identities[label]['linear_solver'] == 'ma57', identities[label]
            identities[label]['loaded_libraries'] = loaded_libraries(p.pid)
            assert identities[label]['affinity'] == [args.cpu] and identities[label]['nice'] == 10
            extension = Path(identities[label]['extension'])
            identities[label]['extension_sha256'] = hashlib.sha256(extension.read_bytes()).hexdigest()
        with gzip.open(raw_path, 'wt') as raw:
            for rep in range(3):
                for op_index, operation in enumerate(operations):
                    labels = [v['version'] for v in environments]
                    shift = (rep * 3 + op_index) % len(labels)
                    labels = labels[shift:] + labels[:shift]
                    if rep % 2:
                        labels.reverse()
                    labels = [*controls, *labels, *reversed(controls)]
                    for ordinal, label in enumerate(labels):
                        before = gate()
                        block = len(blocks)
                        samples = 5 if operation in OPERATIONS[6:13] else 1
                        p = workers[label]
                        p.stdin.write(json.dumps(dict(operation=operation, samples=samples)) + '\n')
                        p.stdin.flush()
                        reply = receive(p)
                        after = snapshot()
                        record = dict(block=block, repetition=rep, ordinal=ordinal, version=label,
                                      operation=operation, before=before, after=after, **reply)
                        raw.write(json.dumps(record, allow_nan=False) + '\n')
                        raw.flush()
                        blocks.append({k: record[k] for k in ('block', 'repetition', 'ordinal', 'version', 'operation', 'before', 'after')})
                        if 'unavailable' in reply:
                            failures.append(dict(block=block, repetition=rep, version=label,
                                                 operation=operation, reason=reply['unavailable']))
                            print(rep, label, operation, reply['unavailable'][:130], flush=True)
                        else:
                            for sample, value in enumerate(reply['values']):
                                assert value['seconds'] > 0 and value['rows']
                                timing_rows.append(dict(block=block, repetition=rep, ordinal=ordinal, version=label,
                                    operation=operation, sample=sample, wall_seconds=value['seconds'],
                                    cpu_seconds=value['cpu_seconds'], load_before=before['load'][0],
                                    load_after=after['load'][0], gate=before['mode']))
                            print(rep, label, operation, f"{reply['values'][0]['seconds']:.6g} s", flush=True)
    finally:
        for label, process in workers.items():
            if process.stdin and not process.stdin.closed:
                process.stdin.close()
            try:
                process.wait(timeout=20)
            except subprocess.TimeoutExpired:
                process.terminate()
                process.wait(timeout=20)
            logs[label].close()
        for name, rows in [('timings.csv', timing_rows), ('unavailable.csv', failures)]:
            if rows:
                with (HERE / 'results' / name).open('w', newline='') as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator='\n')
                    writer.writeheader()
                    writer.writerows(rows)
        (HERE / 'results/campaign.json').write_text(json.dumps(dict(started=started, ended=utc(),
            host=platform.uname()._asdict(), cpu_count=os.cpu_count(), cpu=args.cpu,
            cpu_model=next(line.split(':', 1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')),
            thread_environment={k: environment[k] for k in THREADS}, environments=environments,
            controls=controls, current=CURRENT, hsl_library=str(HSL),
            unavailable_wheels=json.loads((BUILD / 'unavailable-wheels.json').read_text()),
            identities=identities, scheduled_blocks=len(blocks), repetitions=3, full_fit='omitted: serial F1 only',
            time_boundary='native solve calls only; F1 uses the key-solve capture boundary; density per call over 32 changing pressures'), indent=2) + '\n')


if __name__ == '__main__':
    main()

"""Reuse #578's persistent workers, interleaving and timers for current main.

Only declarations, locations and reserved-core per-block admission differ.
"""
import argparse
import ast
import csv
from datetime import datetime, timedelta, timezone
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import time

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
BUILD = ROOT / 'build/issue-929'
SOURCE = ROOT / 'analyses/2026-audit-timing-campaign/scripts/campaign.py'
spec = importlib.util.spec_from_file_location('existing_campaign', SOURCE)
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)
c.HERE, c.ROOT, c.BUILD = HERE, ROOT, BUILD
c.INPUTS = json.loads((HERE / 'inputs/inputs.json').read_text())
c.IMAGE = json.loads((HERE / 'results/image.json').read_text()) if (HERE / 'results/image.json').exists() else {}
c.COMPILERS.update({'c++', 'clang-18', 'clang++-18', 'clang-tidy-18'})
QUIET = c.INPUTS['quiet_host']

original_snapshot = c.snapshot
def snapshot():
    s = original_snapshot()
    rows = subprocess.check_output(['ps', '-eo', 'pid=,ppid=,ni=,comm=,args='], text=True).splitlines()
    s['compile_processes'] = [row.strip() for row in rows if row.split()[3] in c.COMPILERS]
    s['cpu'] = 31
    s['affinity'] = sorted(os.sched_getaffinity(0))
    s['nice'] = os.getpriority(os.PRIO_PROCESS, 0)
    s['core_reserved'] = core_reserved()
    return s


def admit(allowance):
    """Retain deferred admissions; never use the old campaign's loaded-host fallback."""
    begun = time.monotonic()
    record = HERE / 'results/deferred-host.jsonl'
    started = ADMISSION_START
    deadline = started + timedelta(hours=QUIET['maximum_wait_hours'])
    while True:
        s = snapshot()
        if datetime.now(timezone.utc) >= deadline:
            raise TimeoutError('Eight-hour admission deadline reached')
        if admissible(s):
            return s | dict(gate='quiet', waited=time.monotonic() - begun)
        with (HERE / 'results/deferred-host.jsonl').open('a') as stream:
            stream.write(json.dumps(s | dict(disposition='deferred before timing')) + '\n')
        print('DEFERRED', s['utc'], 'load1', s['load'][0], 'compiler PIDs', [row.split()[0] for row in s['compile_processes']], flush=True)
        if datetime.now(timezone.utc) >= deadline:
            raise TimeoutError('No complete quiet set within the owner-authorized waiting interval; deferred handoff required')
        time.sleep(QUIET['poll_seconds'])


class Clapeyron(c.Clapeyron):
    def __init__(self):
        self.label, self.count = 'Clapeyron', 0
        (BUILD / 'requests').mkdir(exist_ok=True)
        ready = next(row for row in (BUILD / 'clapeyron.log').read_text().splitlines() if row.startswith('READY '))
        self.identity = json.loads(ready[6:])
        self.identity.update(affinity=sorted(os.sched_getaffinity(self.identity['pid'])), nice=os.getpriority(os.PRIO_PROCESS, self.identity['pid']))
        assert self.identity['julia'] == '1.13.1' and self.identity['threads'] == self.identity['blas_threads'] == 1
        assert self.identity['affinity'] == [31] and self.identity['nice'] == 10


c.snapshot, c.admit, c.Clapeyron = snapshot, admit, Clapeyron
# Keep the existing interleaving/timing function; require quiet admission before EVERY arm.
tree = ast.parse(SOURCE.read_text())
fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'run')
for n in ast.walk(fn):
    if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'before' for t in n.targets):
        n.value = ast.parse('admit(0)', mode='eval').body
exec(compile(ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[])), str(SOURCE), 'exec'), c.__dict__)


def blocks():
    r = dict(samples=25, warmups=5)
    yield 'property-bundle', {t: dict(op='t1-bundle', state='O-P01/propane-dilute', **r) for t in ('Engine', 'FeOs', 'Clapeyron')}
    yield 'pressure-state', dict(Engine=dict(op='t2-sequence', case='A-propane-liquid', mode='cold', samples=5, warmups=2),
                               **{t: dict(op='t2-sequence', case='A-propane-liquid', samples=5, warmups=2) for t in ('FeOs', 'Clapeyron')})
    for case in ('C1-300', 'C2-0.5', 'C4-298.15'):
        yield case, dict(Engine=dict(op='release-equilibrium', case=case, **r),
                         **{t: dict(op='t3', case=case, **r) for t in ('FeOs', 'Clapeyron')})
    yield 'C1-300-prepared', dict(Engine=dict(op='prepared-equilibrium', case='C1-300', **r),
        **{t: dict(op='t3', case='C1-300', **r) for t in ('FeOs', 'Clapeyron')})
    yield 'electrolyte-eight-feeds', dict(Engine=dict(op='release-electrolyte', samples=5, warmups=2))
    yield 'explicit-discovery', dict(Engine=dict(op='t4', system='methane_ethane', mode='discover', point=0, samples=5, warmups=2),
                                    FeOs=dict(op='t4', system='methane_ethane', T=190., P=1e6, samples=5, warmups=2),
                                    Clapeyron=dict(op='t4', system='methane_ethane', T=190., P=1e6, method='default', samples=5, warmups=2))
    for case in ('C2-0.5', 'C4-298.15'):
        for mode in ('exact', 'fd'):
            yield f'{case}/actions-{mode}', dict(Engine=dict(op='release-action', case=case, mode=mode, samples=10, warmups=2))


def host_identity():
    paths = [SOURCE, ROOT / 'analyses/2026-audit-timing-campaign/scripts/engine_worker.py',
             ROOT / 'analyses/2026-audit-timing-campaign/scripts/feos_worker.py',
             ROOT / 'analyses/2026-audit-timing-campaign/scripts/clapeyron_server.jl',
             HERE / 'scripts/definitions.jl']
    return dict(system=platform.uname()._asdict(), lscpu=subprocess.check_output(['lscpu', '-J'], text=True),
                numerical_threads={key: '1' for key in c.THREADS}, nice=10, cpu=31,
                methods_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                host=snapshot())


def acquire():
    global ADMISSION_START
    ADMISSION_START = datetime.now(timezone.utc)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['preflight', 'measure'])
    parser.add_argument('--name', default='measurement-1')
    args = parser.parse_args()
    plan = list(blocks())
    if args.mode == 'preflight':
        # Numerical answers and compilation warmups are deliberately outside measurement.
        workers = {}
        try:
            for name in ('Engine-A', 'FeOs'):
                workers[name] = c.Worker(name, [c.IMAGE['python'] if name == 'Engine-A' else str(c.FEOS), '-u',
                    str(HERE / 'scripts' / ('engine_worker.py' if name == 'Engine-A' else 'feos_worker.py'))],
                    31, BUILD / f'preflight-{name}.log')
            workers['Clapeyron'] = Clapeyron()
            assert workers['FeOs'].identity['version'] == '0.10.1'
            admit(0)
            out = dict(host=host_identity(), identities={name: w.identity for name, w in workers.items()}, operations={})
            for operation, requests in plan:
                out['operations'][operation] = {}
                for name, worker in workers.items():
                    key = 'Engine' if name == 'Engine-A' else name
                    if key in requests:
                        reply = worker.send(requests[key] | dict(samples=1, warmups=1))
                        out['operations'][operation][key] = reply
                        print(operation, key, reply.get('error', 'returned'), flush=True)
                (HERE / 'results/preflight.json').write_text(json.dumps(out, indent=2) + '\n')
        finally:
            for worker in workers.values():
                worker.close()
        return
    # checks.csv is produced from every preflight quantity; no ratio may bypass a failed match.
    if not (HERE / 'results/checks.csv').exists():
        raise RuntimeError('run preflight and its value checks before measurement')
    checks = list(csv.DictReader((HERE / 'results/checks.csv').open()))
    matched = {(r['operation'], r['comparator']) for r in checks}
    failed = {(r['operation'], r['comparator']) for r in checks if r['passed'] != 'True'}
    for operation, requests in plan:
        for comparator in ('FeOs', 'Clapeyron'):
            if comparator in requests and ((operation, comparator) not in matched or (operation, comparator) in failed):
                del requests[comparator]
        if '/actions-' in operation:
            fd_operation = operation.replace('/actions-exact', '/actions-fd')
            if (fd_operation, 'package central differences') not in matched or (fd_operation, 'package central differences') in failed:
                requests.clear()
    plan = [(operation, requests) for operation, requests in plan if requests]
    if any(op == 'C1-300-prepared' and len(req) != 3 for op, req in plan):
        plan = [(op, req) for op, req in plan if op != 'C1-300-prepared']
    plan = [(op, req) for op, req in plan if (op, 'Engine') not in failed and not (
        op == 'electrolyte-eight-feeds' and any(o == op for o, _ in failed))]
    # Keep the same worker processes across complete operation retries.
    original_worker = c.Worker
    pool = {}
    def persistent(label, command, cpu, log):
        if label not in pool:
            w = original_worker(label, command, cpu, log)
            pool[label] = (w, w.close)
            w.close = lambda: None
        return pool[label][0]
    c.Worker = persistent
    c.close_pool = lambda: [close() for _, close in pool.values()]
    all_plan = plan
    retry = 0
    while plan:
        out = HERE / 'results/raw' / f'{args.name}{"-retry-" + str(retry) if retry else ""}.jsonl.gz'
        if out.exists():
            raise FileExistsError(out)
        c.run('main-cost', plan, out, 31, ['Engine-A', 'Engine-control', 'FeOs', 'Clapeyron'], repetitions=3)
        disturbed = set()
        for line in gzip.open(out, 'rt'):
            r = json.loads(line)
            if r['kind'] == 'block' and any(not admissible(s) for s in (r['before'], r['after'])):
                disturbed.add(r['operation'])
        disturbed |= {op.replace('/actions-fd', '/actions-exact') for op in disturbed if '/actions-fd' in op}
        disturbed |= {op.replace('/actions-exact', '/actions-fd') for op in disturbed if '/actions-exact' in op}
        plan = [(op, req) for op, req in all_plan if op in disturbed]
        retry += 1
        if plan:
            print('RETRY disturbed operations', sorted(disturbed), flush=True)


def core_reserved():
    base = Path('/sys/fs/cgroup/user.slice/user-1000.slice/user@1000.service')
    def cpus(text):
        result = set()
        for part in filter(None, text.strip().split(',')):
            lo, _, hi = part.partition('-')
            result.update(range(int(lo), int(hi or lo) + 1))
        return result
    return all({15, 31}.isdisjoint(cpus((base / s / 'cpuset.cpus.effective').read_text()))
               for s in ('app.slice', 'background.slice', 'session.slice', 'agents.slice'))


def admissible(s):
    return s['core_reserved'] and not s['compilers'] and s['affinity'] == [31] and s['nice'] == 10


def main():
    from provenance import check
    check()
    import selectors
    BUILD.mkdir(parents=True, exist_ok=True)
    fifo = BUILD / 'clapeyron.fifo'
    if not fifo.exists():
        os.mkfifo(fifo)
    env = os.environ | {k: '1' for k in c.THREADS} | {'JULIA_DEPOT_PATH': str(Path.home() / '.cache/epcsaft-comparators/222/depot')}
    julia = Path.home() / '.cache/epcsaft-comparators/222/julia-1.13.1/bin/julia'
    command = ['taskset', '-c', '31', str(julia), '--startup-file=no', '--project=' + str(ROOT / 'analyses/2026-competitor-cross-check/inputs'), str(HERE / 'scripts/clapeyron_server.jl'), str(fifo)]
    log_path = BUILD / 'clapeyron.log'
    if log_path.exists():
        log_path.rename(BUILD / f'clapeyron-{time.time_ns()}.log')
    with log_path.open('w') as log:
        server = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log, text=True, env=env, cwd=ROOT)
        try:
            with selectors.DefaultSelector() as sel:
                sel.register(server.stdout, selectors.EVENT_READ)
                while True:
                    if not sel.select(600):
                        raise TimeoutError('Clapeyron startup exceeded 600 seconds')
                    line = server.stdout.readline()
                    if not line:
                        raise RuntimeError('Clapeyron startup failed; see clapeyron.log')
                    log.write(line); log.flush()
                    if line.startswith('READY '):
                        identity = json.loads(line[6:])
                        assert '0.6.29' in json.dumps(identity), identity
                        break
            acquire()
            check()
        finally:
            if hasattr(c, 'close_pool'):
                c.close_pool()
            server.terminate()
            try: server.wait(timeout=30)
            except subprocess.TimeoutExpired:
                server.kill(); server.wait(timeout=30)


if __name__ == '__main__':
    main()

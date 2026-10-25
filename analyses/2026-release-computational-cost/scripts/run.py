"""Reuse #578's persistent workers, interleaving and timers for the named release.

Only declarations, locations and strict per-block quiet-host admission differ.
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
BUILD = ROOT / 'build/issue-681'
SOURCE = ROOT / 'analyses/2026-audit-timing-campaign/scripts/campaign.py'
spec = importlib.util.spec_from_file_location('existing_campaign', SOURCE)
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)
c.HERE, c.ROOT, c.BUILD = HERE, ROOT, BUILD
c.INPUTS = json.loads((HERE / 'inputs/inputs.json').read_text())
c.IMAGE = json.loads((HERE / 'results/image.json').read_text())
c.COMPILERS.add('c++')
QUIET = c.INPUTS['quiet_host']

original_snapshot = c.snapshot
def snapshot():
    s = original_snapshot()
    rows = subprocess.check_output(['ps', '-eo', 'pid=,ppid=,ni=,comm=,args='], text=True).splitlines()
    s['compile_processes'] = [row.strip() for row in rows if row.split()[3] in c.COMPILERS]
    s['cpu'] = 31
    return s


def admit(allowance):
    """Retain deferred admissions; never use the old campaign's loaded-host fallback."""
    begun = time.monotonic()
    record = HERE / 'results/deferred-host.jsonl'
    started = datetime.fromisoformat(json.loads(record.open().readline())['utc']) if record.exists() else datetime.now(timezone.utc)
    deadline = started + timedelta(hours=QUIET['maximum_wait_hours'])
    while True:
        s = snapshot()
        if s['load'][0] < QUIET['maximum_load1_exclusive'] and not s['compilers']:
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
             Path('/home/tnnrpolley21/Workspaces/Engineering/ePC-SAFT/build/issue-533/scripts/definitions.jl')]
    return dict(system=platform.uname()._asdict(), lscpu=subprocess.check_output(['lscpu', '-J'], text=True),
                numerical_threads={key: '1' for key in c.THREADS}, nice=10, cpu=31,
                methods_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                host=snapshot())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['preflight', 'measure'])
    parser.add_argument('--name', default='measurement-1')
    args = parser.parse_args()
    plan = list(blocks())
    if args.mode == 'preflight':
        # Numerical answers and compilation warmups are deliberately outside measurement.
        workers = {name: c.Worker(name, [c.IMAGE['python'] if name == 'Engine-A' else str(c.FEOS), '-u',
                    str(HERE / 'scripts' / ('engine_worker.py' if name == 'Engine-A' else 'feos_worker.py'))],
                    31, BUILD / f'preflight-{name}.log') for name in ('Engine-A', 'FeOs')}
        workers['Clapeyron'] = Clapeyron()
        try:
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
    retry = 0
    while plan:
        out = HERE / 'results/raw' / f'{args.name}{"-retry-" + str(retry) if retry else ""}.jsonl.gz'
        if out.exists():
            raise FileExistsError(out)
        c.run('release-cost', plan, out, 31, ['Engine-A', 'Engine-control', 'FeOs', 'Clapeyron'], repetitions=3)
        disturbed = set()
        for line in gzip.open(out, 'rt'):
            r = json.loads(line)
            if r['kind'] == 'block' and any(s['load'][0] >= QUIET['maximum_load1_exclusive'] or s['compilers'] for s in (r['before'], r['after'])):
                disturbed.add(r['operation'])
        plan = [(operation, requests) for operation, requests in plan if operation in disturbed]
        retry += 1
        if plan:
            print('RETRY disturbed operations', sorted(disturbed), flush=True)


if __name__ == '__main__':
    main()

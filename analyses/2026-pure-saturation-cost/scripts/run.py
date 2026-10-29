"""C1-300 only: release workers, timers, interleaving and quiet admission."""
import ast
import gzip
import importlib.util
import json
import os
from pathlib import Path
import selectors
import subprocess
import time

STUDY = Path(__file__).resolve().parents[1]
HERE = Path(os.environ.get('EPCSAFT_COST_OUTPUT', str(STUDY)))
ROOT = STUDY.parents[1]
BUILD = HERE / 'cost-workers'
BUILD.mkdir(parents=True, exist_ok=True)
if HERE != STUDY and not (HERE / 'scripts').exists():
    (HERE / 'scripts').symlink_to(STUDY / 'scripts', target_is_directory=True)
source = ROOT / 'analyses/2026-release-computational-cost/scripts/run.py'
spec = importlib.util.spec_from_file_location('release_campaign', source)
r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r)
c = r.c
r.HERE, r.BUILD = HERE, BUILD
c.HERE, c.BUILD = HERE, BUILD
r.QUIET = r.QUIET | dict(poll_seconds=30)
images = {arm: json.loads((HERE / f'results/image-{arm}.json').read_text()) for arm in ('A', 'B')}
c.IMAGE = images['A']
original_worker = c.Worker


class Worker(original_worker):
    def __init__(self, label, command, cpu, log):
        if label.startswith('Engine'):
            arm = 'B' if label.startswith('Engine-B') else 'A'
            image = BUILD / f'guard-{arm}/results/image.json'
            os.environ['EPCSAFT_COST_ARM'] = arm
            os.environ['EPCSAFT_COST_IMAGE'] = str(image)
            os.environ['EPCSAFT_COST_SOURCE'] = str(Path(images[arm]['python']).parents[2])
            command = [images[arm]['python'], '-u', str(STUDY / 'scripts/engine_worker.py')]
        super().__init__(label, command, cpu, log)


c.Worker = Worker
# Preserve the shared loop; add B's two independent workers to its existing interleaving.
tree = ast.parse(c.SOURCE.read_text()) if hasattr(c, 'SOURCE') else ast.parse(r.SOURCE.read_text())
fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'run')
for node in ast.walk(fn):
    if isinstance(node, ast.If) and ast.unparse(node.test) == "label in ('Engine-A', 'Engine-control')":
        node.test = ast.parse("label.startswith('Engine')", mode='eval').body
    if isinstance(node, ast.Assign):
        if any(isinstance(t, ast.Name) and t.id == 'before' for t in node.targets):
            node.value = ast.parse('admit(0)', mode='eval').body
        if any(isinstance(t, ast.Name) and t.id == 'arms' for t in node.targets):
            node.value = ast.parse("['Engine-A', 'Engine-B', 'FeOs', 'Clapeyron', 'Engine-B-control', 'Engine-control']", mode='eval').body
exec(compile(ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[])), str(r.SOURCE), 'exec'), c.__dict__)


original_snapshot = r.snapshot

def snapshot():
    from datetime import datetime
    from zoneinfo import ZoneInfo
    s = original_snapshot()
    local = datetime.now(ZoneInfo('America/Denver'))
    s['protected_window'] = local.date().isoformat() == '2026-10-28' and 165 <= local.hour * 60 + local.minute < 330
    s['compilers'] = s['compilers'] or s['protected_window']
    return s

r.snapshot = snapshot

def quiet_blocks(path):
    disturbed = set()
    for line in gzip.open(path, 'rt'):
        row = json.loads(line)
        if row['kind'] == 'block' and any(s['load'][0] >= r.QUIET['maximum_load1_exclusive'] or s['compilers'] for s in (row['before'], row['after'])):
            disturbed.add(row['operation'])
    return disturbed


def values(reply):
    if 'error' in reply:
        raise RuntimeError(reply['error'])
    return [[sample['value']['P'], *(sample['value']['densities'] if 'densities' in sample['value'] else [p['rho'] for p in sample['value']['phases']])] for sample in reply['samples']]


def main():
    fifo = BUILD / 'clapeyron.fifo'
    if not fifo.exists():
        os.mkfifo(fifo)
    env = os.environ | {k: '1' for k in c.THREADS} | {'JULIA_DEPOT_PATH': str(Path.home() / '.cache/epcsaft-comparators/222/depot')}
    julia = Path.home() / '.cache/epcsaft-comparators/222/julia-1.13.1/bin/julia'
    command = ['taskset', '-c', '31', str(julia), '--startup-file=no', '--project=' + str(ROOT / 'analyses/2026-competitor-cross-check/inputs'), str(ROOT / 'analyses/2026-release-computational-cost/scripts/clapeyron_server.jl'), str(fifo)]
    log = (BUILD / 'clapeyron.log').open('w')
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
                    break
        # The FIFO worker writes progress to this same task-owned log after startup.
        c.Clapeyron = r.Clapeyron
        c.snapshot, c.admit = r.snapshot, r.admit
        plan = [(name, dict(Engine=dict(op=op, case='C1-300', samples=25, warmups=5),
                    **{tool: dict(op='t3', case='C1-300', samples=25, warmups=5) for tool in ('FeOs', 'Clapeyron')}))
                for name, op in [('ordinary', 'release-equilibrium'), ('prepared', 'prepared-equilibrium')]]
        arms = ['Engine-A', 'Engine-control', 'Engine-B', 'Engine-B-control', 'FeOs', 'Clapeyron']
        workers = {name: Worker(name, [str(c.FEOS), '-u', str(STUDY / 'scripts/feos_worker.py')], 31, BUILD / f'preflight-{name}.log')
                   for name in arms if name != 'Clapeyron'}
        workers['Clapeyron'] = c.Clapeyron()
        preflight = dict(identities={name: w.identity for name, w in workers.items()}, operations={})
        try:
            r.admit(0)
            reference = None
            for operation, requests in plan:
                preflight['operations'][operation] = {}
                for name, worker in workers.items():
                    reply = worker.send(requests['Engine' if name.startswith('Engine') else name] | dict(samples=1, warmups=1))
                    preflight['operations'][operation][name] = reply
                    (HERE / 'results/preflight-progress.json').write_text(json.dumps(preflight, indent=2) + '\n')
                    observed = values(reply)[0]
                    print('Preflight', operation, name, observed, flush=True)
                    if reference is None:
                        reference = observed
                    limit = 1e-10 if name.startswith('Engine') else 1e-7
                    assert all(abs(a-b) <= limit * max(abs(a), abs(b)) for a, b in zip(reference, observed)), (name, observed, reference)
            (HERE / 'results/preflight.json').write_text(json.dumps(preflight, indent=2) + '\n')
        finally:
            for worker in workers.values():
                worker.close()
        retry = 0
        while plan:
            out = HERE / f'results/raw/measurement-{retry}.jsonl.gz'
            if out.exists():
                raise FileExistsError(out)
            c.run('pure-saturation-cost', plan, out, 31, arms, repetitions=3)
            disturbed = quiet_blocks(out)
            plan = [(operation, requests) for operation, requests in plan if operation in disturbed]
            retry += 1
        (HERE / 'results/host.json').write_text(json.dumps(r.host_identity(), indent=2) + '\n')
    finally:
        server.terminate()
        server.wait(timeout=30)
        log.close()


if __name__ == '__main__':
    (HERE / 'results/raw').mkdir(exist_ok=True)
    main()

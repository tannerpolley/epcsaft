"""Release interleaving and strict quiet admission for ordinary/prepared A/B discovery."""
import ast
import gzip
import importlib.util
import json
import os
from pathlib import Path
import selectors
import subprocess

HERE = Path(__file__).resolve().parents[1]; ROOT = HERE.parents[1]
MAIN = Path(subprocess.check_output(['git','rev-parse','--path-format=absolute','--git-common-dir'], cwd=ROOT, text=True).strip()).parent
BUILD = MAIN/'build/issue-855/build-r2/cost-workers'; BUILD.mkdir(parents=True, exist_ok=True)
spec = importlib.util.spec_from_file_location('release_campaign', ROOT/'analyses/2026-release-computational-cost/scripts/run.py')
r = importlib.util.module_from_spec(spec); spec.loader.exec_module(r)
c = r.c; r.HERE = c.HERE = HERE; r.BUILD = c.BUILD = BUILD
r.QUIET = r.QUIET | dict(poll_seconds=30)
# A resumed campaign gets its own waiting window; retained deferrals stay intact.
r.CAMPAIGN_STARTED = r.datetime.now(r.timezone.utc)
admission = next(n for n in ast.parse(Path(r.__file__).read_text()).body if isinstance(n, ast.FunctionDef) and n.name == 'admit')
next(n for n in admission.body if isinstance(n, ast.Assign) and n.targets[0].id == 'started').value = ast.Name(id='CAMPAIGN_STARTED', ctx=ast.Load())
exec(compile(ast.fix_missing_locations(ast.Module(body=[admission], type_ignores=[])), str(r.__file__), 'exec'), r.__dict__)
c.admit = r.admit
images = {arm:json.loads((HERE/f'results/image-{arm}.json').read_text()) for arm in ('A','B')}
c.IMAGE = images['A']; original = c.Worker

class Worker(original):
    def __init__(self, label, command, cpu, log):
        if label.startswith('Engine'):
            arm = 'B' if label.startswith('Engine-B') else 'A'
            os.environ['EPCSAFT_COST_IMAGE'] = str(MAIN/f'build/issue-855/build-r{2 if arm == "A" else 1}/guard-{arm}/results/image.json')
            os.environ['EPCSAFT_COST_SOURCE'] = str(Path(images[arm]['python']).parents[2])
            command = [images[arm]['python'], '-u', str(HERE/'scripts/engine_worker.py')]
        elif label == 'FeOs':
            command = [str(c.FEOS), '-u', str(ROOT/'analyses/2026-release-computational-cost/scripts/feos_worker.py')]
        super().__init__(label, command, cpu, log)

c.Worker = Worker
fn = next(n for n in ast.parse(r.SOURCE.read_text()).body if isinstance(n, ast.FunctionDef) and n.name == 'run')
for node in ast.walk(fn):
    if isinstance(node, ast.If) and ast.unparse(node.test) == "label in ('Engine-A', 'Engine-control')":
        node.test = ast.parse("label.startswith('Engine')", mode='eval').body
    if isinstance(node, ast.Assign):
        if any(isinstance(t, ast.Name) and t.id == 'before' for t in node.targets):
            node.value = ast.parse('admit(0)', mode='eval').body
        if any(isinstance(t, ast.Name) and t.id == 'arms' for t in node.targets):
            node.value = ast.parse("['Engine-A','Engine-B','FeOs','Clapeyron','Engine-B-control','Engine-control']", mode='eval').body
exec(compile(ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[])), str(r.SOURCE), 'exec'), c.__dict__)

def main():
    from reduce import physical, check
    (BUILD/'controller.pid').write_text(str(os.getpid()))
    c.admit(0)
    fifo = BUILD/'clapeyron.fifo'
    if not fifo.exists(): os.mkfifo(fifo)
    env = os.environ | {k:'1' for k in c.THREADS} | dict(JULIA_DEPOT_PATH=str(Path.home()/'.cache/epcsaft-comparators/222/depot'))
    command = ['taskset','-c','31',str(Path.home()/'.cache/epcsaft-comparators/222/julia-1.13.1/bin/julia'),
        '--startup-file=no','--project='+str(ROOT/'analyses/2026-competitor-cross-check/inputs'),
        str(ROOT/'analyses/2026-release-computational-cost/scripts/clapeyron_server.jl'),str(fifo)]
    with (BUILD/'clapeyron.log').open('w') as log:
        server = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log, text=True, env=env, cwd=ROOT)
        try:
            with selectors.DefaultSelector() as selector:
                selector.register(server.stdout, selectors.EVENT_READ)
                while True:
                    if not selector.select(600): raise TimeoutError('Clapeyron startup')
                    line = server.stdout.readline()
                    if not line: raise RuntimeError('Clapeyron exited; see its log')
                    log.write(line); log.flush()
                    if line.startswith('READY '): break
            c.Clapeyron = r.Clapeyron
            common = dict(system='methane_ethane', samples=25, warmups=5)
            plan = [(mode, dict(Engine=dict(op=op, mode='discover', point=0, **common),
                FeOs=dict(op='t4',T=190.,P=1e6,**common), Clapeyron=dict(op='t4',T=190.,P=1e6,method='default',**common)))
                for mode,op in [('ordinary','t4'),('prepared','prepared-discovery')]]
            arms = ['Engine-A','Engine-control','Engine-B','Engine-B-control','FeOs','Clapeyron']
            workers, pre = {}, {}
            try:
                for a in arms:
                    workers[a] = c.Clapeyron() if a == 'Clapeyron' else Worker(a, [], 31, BUILD/f'preflight-{a}.log')
                for operation, requests in plan:
                    pre[operation] = {a:w.send(requests['Engine' if a.startswith('Engine') else a] | dict(samples=1,warmups=1)) for a,w in workers.items()}
                    reference = physical(pre[operation]['Engine-A']['samples'][0], 'Engine-A')
                    for a, reply in pre[operation].items():
                        assert 'error' not in reply, reply
                        assert check(reference, physical(reply['samples'][0],a),a), (operation,a,reply)
                (HERE/'results/preflight.json').write_text(json.dumps(pre,indent=1)+'\n')
            finally:
                for w in workers.values(): w.close()
            attempt = 1 + max((int(p.name.split('-')[1].split('.')[0]) for p in (HERE/'results/raw').glob('measurement-*.jsonl.gz')), default=-1)
            while plan:
                path = HERE/f'results/raw/measurement-{attempt}.jsonl.gz'
                assert not path.exists(), path
                c.run('explicit-discovery-cost',plan,path,31,arms,repetitions=3)
                records = [json.loads(s) for s in gzip.open(path,'rt')]
                disturbed = {v['operation'] for v in records if v['kind']=='block' and any(s['load'][0]>=4 or s['compilers'] for s in (v['before'],v['after']))}
                plan = [(op,req) for op,req in plan if op in disturbed]; attempt += 1
        finally:
            server.terminate(); server.wait(timeout=30)

if __name__ == '__main__':
    (HERE/'results/raw').mkdir(exist_ok=True)
    main()

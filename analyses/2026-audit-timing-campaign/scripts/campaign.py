"""#578 campaign controller: one section per invocation, run inside a named user unit with nice 10.

Engine A and Engine control load the same frozen image (results/image.json) through the guard.
FeOs runs in its own interpreter; Clapeyron is the one persistent Julia service reached through a
FIFO. Every block records host load before and after; nothing is dropped on failure.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import selectors
import subprocess
import time

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
BUILD = ROOT / 'build/issue-578'
INPUTS = json.loads((HERE / 'inputs/inputs.json').read_text())
IMAGE = json.loads((HERE / 'results/image.json').read_text())
FEOS = Path.home() / '.cache/epcsaft-comparators/222/feos/bin/python'
THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
           'VECLIB_MAXIMUM_THREADS', 'RAYON_NUM_THREADS', 'JULIA_NUM_THREADS')
COMPILERS = {'cc', 'cc1', 'cc1plus', 'gcc', 'g++', 'clang', 'clang++', 'rustc', 'cargo', 'ninja', 'samu', 'make', 'cmake',
             'meson', 'ld', 'ld.gold', 'ld.lld', 'collect2', 'mold'}


def utc():
    return datetime.now(timezone.utc).isoformat()


def outside_fits():
    """PIDs of native fit workers or fit drivers not started through this campaign's guarded launcher."""
    found = []
    for proc in Path('/proc').iterdir():
        try:
            cmd = (proc / 'cmdline').read_bytes().replace(b'\0', b' ')
        except OSError:
            continue
        if (b'_serve_fit' in cmd or b'replay_fit.py' in cmd) and b'fit_launch.py' not in cmd:
            found.append(int(proc.name))
    return found


def busy_outside(interval=1.0, threshold=0.5):
    """Processes outside this campaign's epcsaft-578 units that used more than `threshold` CPU over `interval` s."""
    def ticks():
        out = {}
        for proc in Path('/proc').iterdir():
            try:
                if 'epcsaft-578' in (proc / 'cgroup').read_text():
                    continue
                fields = (proc / 'stat').read_text().rsplit(')', 1)[1].split()
                out[int(proc.name)] = (int(fields[11]) + int(fields[12]), (proc / 'comm').read_text().strip())
            except (OSError, ValueError, IndexError):
                continue
        return out
    hz = os.sysconf('SC_CLK_TCK')
    first = ticks()
    time.sleep(interval)
    second = ticks()
    busy = [(pid, name, (t - first[pid][0]) / hz / interval) for pid, (t, name) in second.items() if pid in first]
    return sorted(([pid, name, round(cpu, 2)] for pid, name, cpu in busy if cpu > threshold), key=lambda r: -r[2])


def snapshot():
    comms = subprocess.check_output(['ps', '-eo', 'comm='], text=True).split()
    mem = dict(line.split(':') for line in Path('/proc/meminfo').read_text().splitlines())
    loadavg = Path('/proc/loadavg').read_text().split()
    return dict(utc=utc(), load=[float(v) for v in loadavg[:3]], runnable_total=loadavg[3], outside_fits=outside_fits(),
                compilers=sorted({c for c in comms if c in COMPILERS}),
                mem_available_kb=int(mem['MemAvailable'].split()[0]), swap_free_kb=int(mem['SwapFree'].split()[0]))


def admit(allowance):
    """Strict gate: load1 < 2 and no compiler. Wait at most `allowance` seconds, then record the block as loaded."""
    begun = time.monotonic()
    while True:
        s = snapshot()
        if not s['compilers'] and not s['outside_fits'] and s['load'][0] < 2.0:
            return s | dict(gate='quiet', waited=time.monotonic() - begun)
        if time.monotonic() - begun >= allowance and not s['compilers'] and not s['outside_fits']:
            return s | dict(gate='loaded', waited=time.monotonic() - begun)
        if time.monotonic() - begun >= 1800:
            raise RuntimeError(f'compilers {s["compilers"]} or outside fits {s["outside_fits"]} present for 30 min')
        time.sleep(min(60, max(1, allowance)))


class Worker:
    def __init__(self, label, command, cpu, log):
        env = os.environ | {k: '1' for k in THREADS}
        if label.startswith('Engine'):
            env['CAMPAIGN_FIT_LAUNCHER'] = str(BUILD / 'fit-python')
        prefix = ['taskset', '-c', str(cpu)] if cpu is not None else []
        self.label, self.log = label, open(log, 'a')
        self.process = subprocess.Popen(prefix + command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                        stderr=self.log, text=True, env=env, cwd=ROOT)
        self.identity = self.receive(600)['identity']

    def receive(self, timeout):
        with selectors.DefaultSelector() as sel:
            sel.register(self.process.stdout, selectors.EVENT_READ)
            if not sel.select(timeout):
                raise TimeoutError(f'{self.label} exceeded {timeout} s')
        line = self.process.stdout.readline()
        if not line:
            raise RuntimeError(f'{self.label} exited {self.process.poll()}')
        return json.loads(line)

    def send(self, req, timeout=1800):
        self.process.stdin.write(json.dumps(req) + '\n')
        self.process.stdin.flush()
        return self.receive(timeout)

    def close(self):
        self.process.stdin.close()
        try:
            self.process.wait(timeout=30)
        except subprocess.TimeoutExpired:
            self.process.terminate()
            self.process.wait(timeout=30)
        self.log.close()


class Clapeyron:
    """Client of the persistent epcsaft-578-clapeyron service."""

    def __init__(self):
        self.label, self.count = 'Clapeyron', 0
        (BUILD / 'requests').mkdir(exist_ok=True)
        show = subprocess.check_output(['systemctl', '--user', 'show', 'epcsaft-578-clapeyron', '-p', 'MainPID', '-p',
                                        'ActiveState', '-p', 'MemoryMax', '-p', 'CPUAffinity', '-p', 'Nice'], text=True)
        ready = next(l for l in (BUILD / 'clapeyron.log').read_text().splitlines() if l.startswith('READY '))
        self.identity = dict(line.split('=', 1) for line in show.split()) | json.loads(ready[6:]) | dict(unit='epcsaft-578-clapeyron.service')

    def send(self, req, timeout=1800):
        self.count += 1
        stem = BUILD / f'requests/{os.getpid()}-{self.count}'
        out = Path(f'{stem}.out.json')
        Path(f'{stem}.json').write_text(json.dumps(req | dict(out=str(out))))
        with open(BUILD / 'clapeyron.fifo', 'w') as fifo:
            fifo.write(f'{stem}.json\n')
        begun = time.monotonic()
        while not out.exists():
            if time.monotonic() - begun > timeout:
                raise TimeoutError('Clapeyron request timed out')
            time.sleep(0.002)
        reply = json.loads(out.read_text())
        out.unlink()
        Path(f'{stem}.json').unlink()
        return reply

    def close(self):
        pass


# ---------------------------------------------------------------- sections: lists of (operation id, {tool: request})

def t7():
    for case in INPUTS['t7']:
        heavy = case.startswith('mea-')
        yield case, dict(Engine=dict(op='key', case=case, warmups=0 if heavy else 1, samples=1 if heavy else 3))


def t3():
    for case in INPUTS['t3']:
        r = dict(warmups=5, samples=25)
        yield case, dict(Engine=dict(op='key', case=case) | r, FeOs=dict(op='t3', case=case) | r, Clapeyron=dict(op='t3', case=case) | r)


def t1():
    r = dict(warmups=5, samples=25)
    for s in INPUTS['t1']:
        tools = dict(Engine=dict(op='t1-bundle', state=s['id']) | r)
        tools |= {t: dict(op='t1-bundle', state=s['id']) | r for t in s['comparators'] if not (t == 'Clapeyron' and s['model'] == 'water-nacl')}
        yield f"{s['id']}/bundle", tools
        if s['neutral_property']:
            for order in ('Value', 'Gradient', 'Hessian'):
                yield f"{s['id']}/helmholtz-{order}", dict(Engine=dict(op='t1-neutral', state=s['id'], order=order) | r,
                                                           Clapeyron=dict(op='t1-derivs', state=s['id'], order=order) | r)
        for q in s['requests']:
            if not q['inadmissible_chart']:
                yield f"{s['id']}/action-{'+'.join(q['names'])}", dict(Engine=dict(op='t1-action', state=s['id'], names=q['names']) | r)


def t2(roots):
    for c in INPUTS['t2']['ordinary']:
        r = dict(warmups=5, samples=5)
        yield f"{c['id']}/pressure-cold", dict(Engine=dict(op='t2-sequence', case=c['id'], mode='cold') | r,
                                               FeOs=dict(op='t2-sequence', case=c['id']), Clapeyron=dict(op='t2-sequence', case=c['id']))
        yield f"{c['id']}/pressure-anchored", dict(Engine=dict(op='t2-sequence', case=c['id'], mode='anchored') | r)
        yield f"{c['id']}/fixed-density", dict(Engine=dict(op='t2-sequence', case=c['id'], mode='fixed-density', roots=roots[c['id']]) | r)
    for c in INPUTS['t2']['ethane_loop']:
        r = dict(warmups=5, samples=25)
        yield f"{c['id']}/pressure-cold", dict(Engine=dict(op='t2-loop', case=c['id'], mode='cold') | r,
                                               FeOs=dict(op='t2-loop', case=c['id']) | r, Clapeyron=dict(op='t2-loop', case=c['id']) | r)
        yield f"{c['id']}/pressure-anchored", dict(Engine=dict(op='t2-loop', case=c['id'], mode='anchored') | r)
        yield f"{c['id']}/fixed-density", dict(Engine=dict(op='t2-loop', case=c['id'], mode='fixed-density', roots=roots[c['id']]) | r)
        for combo in ([0], [1], [0, 0], [0, 1], [1, 1]):
            yield f"{c['id']}/pressure-action-{combo}", dict(Engine=dict(op='t2-actions', case=c['id'], combo=combo) | r)


def t4():
    r = dict(warmups=2, samples=10)
    for system in ('methane_ethane', 'methanol_hexane'):
        d = INPUTS['t4'][system]
        points = d['points'] if system == 'methane_ethane' else [[T, d['P']] for T in d['temperatures']]
        for i, (T, P) in enumerate(points):
            yield f'{system}/{T}-{P}/discover', dict(Engine=dict(op='t4', system=system, mode='discover', point=i) | r,
                                                     FeOs=dict(op='t4', system=system, T=T, P=P) | r,
                                                     Clapeyron=dict(op='t4', system=system, T=T, P=P, method='default') | r)
            yield f'{system}/{T}-{P}/michelsen', dict(Clapeyron=dict(op='t4', system=system, T=T, P=P, method='michelsen') | r)
        for mode in ('discover', 'declared-cold', 'declared-sequence'):
            yield f'{system}/sequence/{mode}', dict(Engine=dict(op='t4', system=system, mode=mode, warmups=1, samples=5))


def t6():
    for mode in ('reused', 'fresh'):
        yield f'nacl/{mode}', dict(Engine=dict(op='t6', mode=mode, warmups=1, samples=5))
    yield 'nacl/clapeyron-curve', dict(Clapeyron=dict(op='t6-curve', warmups=1, samples=5))


def t8():
    import itertools
    r = dict(warmups=5, samples=25)
    for case, n in (('isomerisation', 3), ('bubble-C2', 4), ('ionic-lle', 4)):
        for mode in ('first-batch', 'first-separate'):
            for count in (1, 2, 'all'):
                yield f'{case}/{mode}-{count}', dict(Engine=dict(op='t8', case=case, mode=mode, n=count) | r)
        for k in (2, 3):
            for combo in itertools.combinations_with_replacement(range(n), k):
                yield f'{case}/order-{k}-{list(combo)}', dict(Engine=dict(op='t8', case=case, mode='higher', combo=list(combo)) | r)


def t9_serial():
    for family, heavy in (('propane', False), ('electrolyte', False), ('F1', True)):
        r = dict(warmups=0, samples=1) if heavy else dict(warmups=1, samples=3)
        for shift in (0.0, 1e-4):
            yield f'{family}/evaluate-shift-{shift}', dict(Engine=dict(op='t9-evaluate', family=family, shift=shift) | r)


def t10():
    for function in ('property', 'equilibrium'):
        for fgh in (0, 1, 2):
            for controls in ('same', 'alternating'):
                r = dict(warmups=2, samples=10)
                yield f'{function}/fgh{fgh}/{controls}', dict(
                    Engine=dict(op='t10-node', function=function, fgh=fgh, controls=controls) | r)
                yield f'{function}/fgh{fgh}/{controls}/native', dict(
                    Engine=dict(op='t10-native', function=function, fgh=fgh, controls=controls) | r)


def t11(states):
    for name in INPUTS['t11']['sequences']:
        yield f'sequence/{name}', dict(Engine=dict(op='t11', sequence=name, samples=3))
    if states:
        yield 'replay/t4-two-phase', dict(Engine=dict(op='t11-replay', states=states, passes=50, samples=3))


def run_cohorts(out, cpu_count, allowance):
    """T10 concurrency: one fresh Engine process per cohort and arm, so record slots never carry over."""
    raw = gzip.open(out, 'at')
    try:
        cohorts = [(f, n, recs) for f in ('property', 'equilibrium') for n, recs in ((1, 'same'), (16, 'same'), (16, 'distinct'),
                                                                                     (cpu_count, 'same'))]
        for rep in range(2):
            for function, threads, records in cohorts:
                for arm in (('Engine-A', 'Engine-control') if rep == 0 else ('Engine-control', 'Engine-A')):
                    host = admit(allowance)
                    w = Worker(arm, [IMAGE['python'], '-u', str(HERE / 'scripts/engine_worker.py')], None, BUILD / 'T10-cohorts.log')
                    try:
                        req = dict(op='t10-concurrent', function=function, threads=threads, records=records, calls=20, fgh=2)
                        reply = w.send(req)
                        raw.write(json.dumps(dict(kind='block', section='T10-cohorts', operation=f'{function}/{threads}-{records}',
                                                  repetition=rep, arm=arm, request=req, before=host, after=snapshot(),
                                                  identity=w.identity, reply=reply)) + '\n')
                        raw.flush()
                        print('T10 cohort', rep, arm, function, threads, records, reply.get('error') or
                              round(reply['samples'][0]['nodes_per_second'], 1), flush=True)
                    finally:
                        w.close()
    finally:
        raw.close()


def run_fits(out, allowance):
    """T9 fits: unpinned native workers, strict admission (bounded wait) before every fit, end observations retained."""
    raw = gzip.open(out, 'at')
    ncpu = len(os.sched_getaffinity(0))
    w = Worker('Engine-A', [IMAGE['python'], '-u', str(HERE / 'scripts/engine_worker.py')], None, BUILD / 'T9-fits-A.log')
    c = Worker('Engine-control', [IMAGE['python'], '-u', str(HERE / 'scripts/engine_worker.py')], None, BUILD / 'T9-fits-control.log')
    arms = {'Engine-A': w, 'Engine-control': c}
    try:
        raw.write(json.dumps(dict(kind='identities', section='T9-fits', ncpu=ncpu, identities={k: v.identity for k, v in arms.items()})) + '\n')
        plan = [(f, workers, dict(iterations=1, seconds=60)) for f in ('propane', 'electrolyte', 'F1') for workers in (1, max(1, ncpu // 2), ncpu)]
        plan += [('F1', 0, dict(iterations=40, seconds=120))] * 3
        for index, (family, workers, budget) in enumerate(plan):
            order = ('Engine-A', 'Engine-control') if index % 2 == 0 else ('Engine-control', 'Engine-A')
            for arm in order:
                before = admit(allowance) | dict(busy_outside=busy_outside())
                req = dict(op='t9-fit', family=family, workers=workers, **budget)
                t = time.perf_counter()
                try:
                    reply = arms[arm].send(req, timeout=900)
                except Exception as error:
                    reply = dict(error=f'{type(error).__name__}: {error}')
                after = snapshot() | dict(busy_outside=busy_outside())
                raw.write(json.dumps(dict(kind='block', section='T9-fits', operation=f'{family}/workers-{workers}/iterations-{budget["iterations"]}',
                                          repetition=index, arm=arm, request=req, before=before, after=after,
                                          block_wall=time.perf_counter() - t, reply=reply)) + '\n')
                raw.flush()
                v = reply.get('samples', [{}])[0]
                print('T9 fit', family, workers, budget['iterations'], arm, before['load'][0], '->', after['load'][0],
                      reply.get('error') or (round(v['seconds'], 3), v['value']['status'], v['value']['iterations']), flush=True)
        for label, worker in arms.items():
            raw.write(json.dumps(dict(kind='recheck', section='T9-fits', arm=label, reply=worker.send(dict(op='recheck')))) + '\n')
    finally:
        for worker in arms.values():
            worker.close()
        raw.close()


def run(section, blocks, out, cpu, tools, repetitions=3, allowance=0):
    workers = {}
    raw = gzip.open(out, 'at')
    try:
        for label in tools:
            if label in ('Engine-A', 'Engine-control'):
                workers[label] = Worker(label, [IMAGE['python'], '-u', str(HERE / 'scripts/engine_worker.py')], cpu, BUILD / f'{section}-{label}.log')
            elif label == 'FeOs':
                workers[label] = Worker(label, [str(FEOS), '-u', str(HERE / 'scripts/feos_worker.py')], cpu, BUILD / f'{section}-FeOs.log')
            elif label == 'Clapeyron':
                workers[label] = Clapeyron()
        scripts = {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted((HERE / 'scripts').iterdir()) if f.is_file()}
        raw.write(json.dumps(dict(kind='identities', section=section, utc=utc(), cpu=cpu, scripts=scripts,
                                  identities={k: w.identity for k, w in workers.items()})) + '\n')
        start = admit(allowance)
        raw.write(json.dumps(dict(kind='section-start', section=section, host=start)) + '\n')
        print(section, 'start gate', start['gate'], start['load'], flush=True)
        blocks = list(blocks)
        for rep in range(repetitions):
            for index, (operation, requests) in enumerate(blocks):
                arms = ['Engine-A'] + [t for t in ('FeOs', 'Clapeyron') if t in requests and t in workers] + ['Engine-control']
                if rep % 2:
                    arms.reverse()
                for ordinal, arm in enumerate(arms):
                    tool = 'Engine' if arm.startswith('Engine') else arm
                    if arm not in workers or tool not in requests:
                        continue
                    req = requests[tool]
                    before = snapshot()
                    t = time.perf_counter()
                    try:
                        reply = workers[arm].send(req)
                    except Exception as error:
                        reply = dict(error=f'{type(error).__name__}: {error}')
                    wall = time.perf_counter() - t
                    after = snapshot()
                    raw.write(json.dumps(dict(kind='block', section=section, operation=operation, repetition=rep, ordinal=ordinal,
                                              arm=arm, request=req, before=before, after=after, block_wall=wall, reply=reply)) + '\n')
                    if 'error' in reply:
                        print(section, rep, operation, arm, 'ERROR', reply['error'][:200], flush=True)
            raw.flush()
            print(section, 'repetition', rep, 'done', snapshot()['load'], flush=True)
        for label, w in workers.items():
            if label.startswith('Engine'):
                raw.write(json.dumps(dict(kind='recheck', section=section, arm=label, reply=w.send(dict(op='recheck')))) + '\n')
    finally:
        for w in workers.values():
            w.close()
        raw.close()


def prepare_t2(cpu):
    """Untimed: accepted cold roots for the fixed-density diagnostic and native root diagnostics, retained once."""
    path = HERE / 'results/raw/T2-prepare.json'
    if path.exists():
        return json.loads(path.read_text())['roots']
    w = Worker('Engine-A', [IMAGE['python'], '-u', str(HERE / 'scripts/engine_worker.py')], cpu, BUILD / 'T2-prepare.log')
    try:
        roots, replies = {}, {}
        for c in INPUTS['t2']['ordinary']:
            replies[c['id']] = w.send(dict(op='t2-sequence', case=c['id'], mode='cold'))
            roots[c['id']] = [r['rho'] for r in replies[c['id']]['samples'][0]['value']]
        for c in INPUTS['t2']['ethane_loop']:
            replies[c['id']] = w.send(dict(op='t2-loop', case=c['id'], mode='cold'))
            roots[c['id']] = [r['rho'] for r in replies[c['id']]['samples'][0]['value']]
        replies['diagnostics'] = w.send(dict(op='t2-diagnostics'))
        path.write_text(json.dumps(dict(identity=w.identity, roots=roots, replies=replies)) + '\n')
    finally:
        w.close()
    return roots


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('section')
    parser.add_argument('--cpu', type=int, default=31)
    parser.add_argument('--allowance', type=float, default=300, help='quiet-host wait at section start, seconds')
    parser.add_argument('--name', help='output name when a section is rerun; earlier records are retained')
    args = parser.parse_args()
    out = HERE / f'results/raw/{args.name or args.section}.jsonl.gz'
    out.parent.mkdir(exist_ok=True)
    if out.exists():
        raise FileExistsError(f'{out} exists; retain it and name a new section')
    engines = ['Engine-A', 'Engine-control']
    if args.section == 'T7':
        run('T7', t7(), out, args.cpu, engines)
    elif args.section == 'T3':
        run('T3', t3(), out, args.cpu, engines + ['FeOs', 'Clapeyron'], allowance=args.allowance)
    elif args.section == 'T1-verify':
        run('T1-verify', ((s['id'], dict(Engine=dict(op='t1-verify', state=s['id']))) for s in INPUTS['t1']), out, args.cpu,
            ['Engine-A'], repetitions=1)
    elif args.section == 'T1':
        run('T1', t1(), out, args.cpu, engines + ['FeOs', 'Clapeyron'], allowance=args.allowance)
    elif args.section == 'T2':
        roots = prepare_t2(args.cpu)
        run('T2', t2(roots), out, args.cpu, engines + ['FeOs', 'Clapeyron'], allowance=args.allowance)
    elif args.section == 'T4':
        run('T4', t4(), out, args.cpu, engines + ['FeOs', 'Clapeyron'], allowance=args.allowance)
    elif args.section == 'T6':
        run('T6', t6(), out, args.cpu, engines + ['Clapeyron'], allowance=args.allowance)
    elif args.section == 'T8-verify':
        run('T8-verify', ((c, dict(Engine=dict(op='t8-verify', case=c))) for c in ('isomerisation', 'bubble-C2', 'ionic-lle')),
            out, args.cpu, ['Engine-A'], repetitions=1)
    elif args.section == 'T8':
        run('T8', t8(), out, args.cpu, engines, allowance=args.allowance)
    elif args.section == 'T9-serial':
        run('T9-serial', t9_serial(), out, args.cpu, engines, allowance=args.allowance)
    elif args.section == 'T9-fits':
        run_fits(out, args.allowance)
    elif args.section == 'T10':
        run('T10', t10(), out, args.cpu, engines, allowance=args.allowance)
    elif args.section == 'T10-cohorts':
        run_cohorts(out, len(os.sched_getaffinity(0)), args.allowance)
    elif args.section == 'T11':
        states = []
        t4raw = HERE / 'results/raw/T4.jsonl.gz'
        if t4raw.exists():
            for line in gzip.open(t4raw, 'rt'):
                r = json.loads(line)
                if r.get('arm') == 'Engine-A' and r.get('operation') == 'methane_ethane/sequence/discover' and r['repetition'] == 0:
                    for row in r['reply']['samples'][0]['value']:
                        if row.get('phases') == 2:
                            n = 2
                            states += [(row['T'], row['molar_densities'][j], row['mole_fractions'][n * j:n * j + n]) for j in range(2)]
        run('T11', t11(states), out, args.cpu, engines, allowance=args.allowance)
        run('T11-check', [('fresh-model-check', dict(Engine=dict(op='t11-check')))], HERE / 'results/raw/T11-check.jsonl.gz',
            args.cpu, ['Engine-A'], repetitions=1)
    elif args.section == 'T14-replay':
        run('T14-replay', [('benzene-toluene-ph-round-trip', dict(Engine=dict(op='t14-ph')))], out, args.cpu, ['Engine-A'], repetitions=1)
    else:
        raise SystemExit(f'unknown section {args.section}')


if __name__ == '__main__':
    main()

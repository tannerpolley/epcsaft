"""Persistent FeOs 0.10.1 comparator worker for #578 (same request protocol as engine_worker.py).

Parameter conversions are the audit probes' and the #222 cross-check's; values are compared
with the Engine before any time is ranked.
"""
import importlib.metadata
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time
import traceback

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
sys.path = [p for p in sys.path if Path(p or '.').resolve() != Path(__file__).parent]  # never shadow feos
import numpy as np  # noqa: E402
import feos  # noqa: E402
import si_units as si  # noqa: E402

feos.set_num_threads(1)
INPUTS = json.loads((HERE / 'inputs/inputs.json').read_text())
T1 = {r['id']: r for r in INPUTS['t1']}
R = feos.Contributions.Residual
RHO = si.MOL / si.METER**3
CROSS = json.loads((ROOT / 'analyses/2026-competitor-cross-check/inputs/cases.json').read_text())
DENSITY = json.loads((ROOT / 'analyses/2026-density-root-benchmark/inputs/cases.json').read_text())


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


cc = load(ROOT / 'analyses/2026-competitor-cross-check/scripts/feos.py', 'feos_conversion578')
sys.path.insert(0, str(ROOT / 'analyses/2026-density-root-benchmark/scripts'))
dc = load(ROOT / 'analyses/2026-density-root-benchmark/scripts/competitors.py', 'density_competitors578')


def record(c):
    """Opus feos_p01.py:rec, unchanged."""
    kw = dict(m=c['m'], sigma=c['sigma'], epsilon_k=c['eps'])
    if 'assoc' in c:
        kw.update(association_sites=[dict(epsilon_k_ab=c['assoc'][0], kappa_ab=c['assoc'][1], na=1.0, nb=1.0)])
    if 'mu' in c:
        kw.update(mu=c['mu'])
    return feos.PureRecord(feos.Identifier(name=c['cid']), c['M'] * 1000.0, **kw)


def eos(key, kij=0.0):
    if key == 'water-nacl':  # Opus feos_p07.py, unchanged inputs
        I = feos.Identifier
        water = feos.PureRecord(I(name='water'), 18.0153, m=1.2047, sigma=2.7927, epsilon_k=353.95, z=0.0,
                                permittivity_record={'ExperimentalData': {'data': [[298.15, 78.09]]}},
                                association_sites=[dict(epsilon_k_ab=2425.7, kappa_ab=0.0451, na=1.0, nb=1.0)])
        na = feos.PureRecord(I(name='sodium-cation'), 22.98977, m=1.0, sigma=2.8232, epsilon_k=230.0, z=1.0,
                             permittivity_record={'ExperimentalData': {'data': [[298.15, 8.0]]}})
        cl = feos.PureRecord(I(name='chloride-anion'), 35.453, m=1.0, sigma=2.756, epsilon_k=170.0, z=-1.0,
                             permittivity_record={'ExperimentalData': {'data': [[298.15, 8.0]]}})
        binaries = [feos.BinaryRecord(I(name='water'), I(name='sodium-cation'), k_ij=0.00045485),
                    feos.BinaryRecord(I(name='water'), I(name='chloride-anion'), k_ij=-0.25),
                    feos.BinaryRecord(I(name='sodium-cation'), I(name='chloride-anion'), k_ij=0.317)]
        return feos.EquationOfState.epcsaft(feos.Parameters.from_records([water, na, cl], binaries), epcsaft_variant='advanced')
    comps = [INPUTS['components'][n] for n in INPUTS['models'][key]]
    if len(comps) == 1:
        params = feos.Parameters.new_pure(record(comps[0]))
    else:
        params = feos.Parameters.new_binary([record(c) for c in comps], k_ij=kij) if kij else feos.Parameters.new_binary([record(c) for c in comps])
    return feos.EquationOfState.pcsaft(params)


def bundle(st, T):
    return dict(P=st.pressure() / si.PASCAL, h_res=st.molar_enthalpy(R) / (si.JOULE / si.MOL),
                cv_res=st.molar_isochoric_heat_capacity(R) / (si.JOULE / si.MOL / si.KELVIN),
                lnphi=list(st.ln_phi()), a_res=st.molar_helmholtz_energy(R) / (si.RGAS * T * si.KELVIN))


def op_t1_bundle(req):
    r = T1[req['state']]
    e = eos(r['model'])
    T, rho, x = r['T'], r['rho'], np.array(r['x'])
    return lambda: dict(value=bundle(feos.State(e, T * si.KELVIN, density=rho * RHO, composition=x), T))


def op_t2_sequence(req):
    case = next(c for c in DENSITY['cases'] if c['id'] == req['case'])

    def call():
        out = dc.feos_case(DENSITY, case, cc, feos.State, np, si.KELVIN, si.PASCAL, RHO)
        unique = [s for s in out['samples'] if s['mode'] == dc.PRIMARY_MODE]
        seconds = sum(s['elapsed_us'] for s in unique if s['elapsed_us'] is not None) / 1e6
        return dict(seconds=seconds, per_call=len(unique), value=unique)
    return call


def op_t2_loop(req):
    loop = next(c for c in INPUTS['t2']['ethane_loop'] if c['id'] == req['case'])
    e = eos('ethane')

    def call():
        rows, ns = [], 0
        for P in loop['pressures']:
            t = time.perf_counter_ns()
            try:
                st = feos.State(e, loop['T'] * si.KELVIN, pressure=P * si.PASCAL, density_initialization=loop['branch'])
                ns += time.perf_counter_ns() - t
                rows.append(dict(rho=st.density / RHO, P=st.pressure() / si.PASCAL, dPdrho=st.dp_drho() / (si.PASCAL / RHO)))
            except Exception as error:
                ns += time.perf_counter_ns() - t
                rows.append(dict(error=f'{type(error).__name__}: {error}'))
        return dict(seconds=ns / 1e9, per_call=len(rows), value=rows)
    return call


def op_t3(req):
    name = req['case']
    if name == 'methane-bubble':
        q = dict(id=name, case='methane-ethane', operation='bubble', T=180., x=[.5, .5])
        inputs = CROSS | {'models': CROSS['models'] | {'methane-ethane': DENSITY['comparator_inputs']['models']['methane-ethane']}}
    else:
        q, inputs = next(q for q in CROSS['requests'] if q['id'] == name), CROSS
    prepared = cc.prepare(inputs, q)

    def call():
        out, _ = cc.evaluate(q, prepared)
        out.pop('reference_state', None)
        return dict(value=out)
    return call


def op_t4(req):
    """FeOs tp_flash without an initial state (its own stability-analysis start); a stable feed is a returned error."""
    system = req['system']
    e = eos('methane-ethane') if system == 'methane_ethane' else eos('methanol-hexane', kij=INPUTS['t4']['methanol_hexane']['kij'])
    feed = np.array([.5, .5]) * si.MOL

    def call():
        try:
            f = feos.PhaseEquilibrium.tp_flash(e, req['T'] * si.KELVIN, req['P'] * si.PASCAL, feed)
        except Exception as error:
            stable = 'No phase split' in str(error)  # FeOs found the feed stable; anything else is a failed flash
            return dict(value=dict(phases=1 if stable else None, error=f'{type(error).__name__}: {error}'))
        return dict(value=dict(phases=2, compositions=[f.liquid.molefracs.tolist(), f.vapor.molefracs.tolist()],
                               densities=[f.liquid.density / RHO, f.vapor.density / RHO], vapor_fraction=f.vapor_phase_fraction))
    return call


OPS = {name[3:].replace('_', '-'): f for name, f in globals().items() if name.startswith('op_')}


def clean(x):
    if isinstance(x, float):
        return x if math.isfinite(x) else None
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.generic):
        return clean(x.item())
    return x


def main():
    replies = os.fdopen(os.dup(1), 'w')
    os.dup2(2, 1)
    print(json.dumps(dict(identity=dict(pid=os.getpid(), library='FeOs', version=importlib.metadata.version('feos'),
                                        python=sys.executable, affinity=sorted(os.sched_getaffinity(0)),
                                        nice=os.getpriority(os.PRIO_PROCESS, 0)))), file=replies, flush=True)
    calls = {}
    for line in sys.stdin:
        req = json.loads(line)
        try:
            key = json.dumps({k: v for k, v in req.items() if k not in ('samples', 'warmups', 'block')}, sort_keys=True)
            prepare = None
            if key not in calls:
                t = time.perf_counter_ns()
                calls[key] = OPS[req['op']](req)
                prepare = (time.perf_counter_ns() - t) / 1e9
            f = calls[key]
            for _ in range(req.get('warmups', 0)):
                f()
            samples = []
            for _ in range(req.get('samples', 1)):
                t, c = time.perf_counter_ns(), time.process_time_ns()
                try:
                    out = f()
                except Exception as error:  # a failed call is retained, never dropped
                    out = dict(error=f'{type(error).__name__}: {error}')
                out.setdefault('seconds', (time.perf_counter_ns() - t) / 1e9)
                out.setdefault('cpu_seconds', (time.process_time_ns() - c) / 1e9)
                samples.append(clean(out))
            reply = dict(samples=samples, prepare_seconds=prepare)
        except Exception as error:
            traceback.print_exc()
            reply = dict(error=f'{type(error).__name__}: {error}')
        print(json.dumps(reply, allow_nan=False), file=replies, flush=True)


if __name__ == '__main__':
    main()

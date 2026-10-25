"""Freeze every #578 input before timing: states, directions, verification steps and sequences.

Reads the audit probe declarations (main checkout build/audit-2026-10-15, hashed) and repository
inputs; writes inputs/inputs.json. Pure arithmetic: no Engine import, no thermodynamic call.
"""
import ast
import hashlib
import itertools
import json
import math
from pathlib import Path
import random

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
AUDIT = Path('/home/tnnrpolley21/Workspaces/Engineering/ePC-SAFT/build/audit-2026-10-15')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def epc_components():
    """The Opus probe component constants (Gross & Sadowski 2001/2002, Gross & Vrabec 2006)."""
    tree = ast.parse((AUDIT / 'opus/probes/epc.py').read_text())
    out = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Call) and getattr(node.value.func, 'id', '') == 'dict':
            out[node.targets[0].id] = {k.arg: ast.literal_eval(k.value) for k in node.value.keywords}
    return out


def steps(T, rho, x, directions, order, thermal=None):
    """Design R1 rule: deterministic admissible steps, one per direction in a k-th order stencil."""
    out = []
    for dT, drho, dx in directions:
        b = [v for v in (0.02 / abs(dT) if dT else None, 1e-4 * rho / abs(drho) if drho else None,
                         0.002 / max(map(abs, dx)) if any(dx) else None) if v is not None]
        h = min(b)
        if dT:
            margin = min([T] + ([min(T - lo, hi - T) for lo, hi in thermal] if thermal else []))
            h = min(h, margin / (4 * order * abs(dT)))
        if drho:
            h = min(h, rho / (4 * order * abs(drho)))
        for xi, di in zip(x, dx):
            if di:
                h = min(h, xi / (4 * order * abs(di)))
        out.append(h)
    return out


def stencil_points(T, rho, x, directions, h):
    points = []
    for signs in itertools.product((-1.0, 1.0), repeat=len(directions)):
        pT = T + sum(s * hj * d[0] for s, hj, d in zip(signs, h, directions))
        prho = rho + sum(s * hj * d[1] for s, hj, d in zip(signs, h, directions))
        px = [xi + sum(s * hj * d[2][i] for s, hj, d in zip(signs, h, directions)) for i, xi in enumerate(x)]
        points.append(dict(signs=signs, T=pT, rho=prho, x=px))
    return points


def admissible(T, rho, x, point, charges):
    ok = (all(math.isfinite(v) for v in (point['T'], point['rho'], *point['x'])) and point['T'] > 0 and point['rho'] > 0
          and all(p > 0 for p, c in zip(point['x'], x) if c > 0) and all(p == 0 for p, c in zip(point['x'], x) if c == 0)
          and abs(sum(point['x']) - 1) < 1e-12)
    if charges:
        ok = ok and abs(sum(z * v for z, v in zip(charges, point['x']))) <= 1e-12
    return ok


def t1(comps):
    """18 audit states plus the nine-species MEA liquid and the associating-binary action fixture."""
    density = json.loads((ROOT / 'analyses/2026-density-root-benchmark/inputs/cases.json').read_text())
    mea = next(c for c in density['cases'] if c['id'] == 'B-MEA-liquid')
    rows = [
        ('O-P01/propane-dilute', 'propane', 300., 500., [1.]), ('O-P01/propane-liquid', 'propane', 300., 11000., [1.]),
        ('O-P01/methanol-liquid', 'methanol', 300., 24000., [1.]), ('O-P01/methanol-vapor', 'methanol', 450., 300., [1.]),
        ('O-P01/acetone-polar', 'acetone-polar', 300., 13000., [1.]),
        ('O-P01/methane-ethane', 'methane-ethane', 230., 8000., [.4, .6]),
        ('O-P01/methanol-hexane', 'methanol-hexane', 320., 7000., [.3, .7]),
        ('O-P07/nacl-0.01', 'water-nacl', 298.15, 55000., [.98, .01, .01]),
        ('O-P07/nacl-0.05', 'water-nacl', 298.15, 52000., [.9, .05, .05]),
        ('O-P07/nacl-0.10', 'water-nacl', 298.15, 50000., [.8, .1, .1]),
        ('A-Q01/300-5000', 'gross-methane-ethane', 300., 5000., [.4, .6]),
        ('A-Q01/180-200', 'gross-methane-ethane', 180., 200., [.5, .5]),
        ('A-Q01/180-18000', 'gross-methane-ethane', 180., 18000., [.5, .5]),
        ('A-Q01/300-dilute', 'gross-methane-ethane', 300., 1e-10, [.4, .6]),
        ('Sol-1/ethane-100', 'ethane', 300., 100., [1.]), ('Sol-1/ethane-10000', 'ethane', 300., 10000., [1.]),
        ('Sol-1/ethane-loop', 'ethane', 308.8, 5930.399614704038, [1.]),
        ('O-S1/MEA-liquid', 'mea', mea['T_K'], mea['states'][0]['oracle_rho_mol_m3'], mea['x']),
        ('A-Q05/associating-binary', 'associating-binary', 313.15, 42441.09427799721, [.25, .75]),
    ]
    charges = {'water-nacl': [0, 1, -1]}
    out = []
    for sid, model, T, rho, x in rows:
        n = len(x)
        directions = [('T', (1., 0., [0.] * n)), ('rho', (0., 1., [0.] * n))]
        directions += [(f'x{i}-x{n - 1}', (0., 0., [1. if j == i else -1. if j == n - 1 else 0. for j in range(n)]))
                       for i in range(n - 1)]
        top = 3 if model == 'associating-binary' else 2
        requests = []
        for order in range(1, top + 1):
            for combo in itertools.combinations_with_replacement(range(len(directions)), order):
                ds = [directions[j][1] for j in combo]
                z = charges.get(model)
                inadmissible = bool(z) and any(abs(sum(a * b for a, b in zip(z, d[2]))) > 0 for d in ds)
                levels = []
                if not inadmissible:
                    h = steps(T, rho, x, ds, order)
                    for scale in (1., .5):
                        hs = [v * scale for v in h]
                        pts = stencil_points(T, rho, x, ds, hs)
                        levels.append(dict(h=hs, admissible=all(admissible(T, rho, x, p, z) for p in pts), points=len(pts)))
                requests.append(dict(order=order, names=[directions[j][0] for j in combo],
                                     inadmissible_chart='charge-violating composition direction' if inadmissible else None,
                                     levels=levels))
        comparators = [] if model in ('mea', 'associating-binary') else ['FeOs', 'Clapeyron']
        out.append(dict(id=sid, model=model, T=T, rho=rho, x=x, directions=dict(directions), requests=requests,
                        comparators=comparators, neutral_property=n == 1 and model in ('propane', 'ethane')))
    return out


def t2():
    density = json.loads((ROOT / 'analyses/2026-density-root-benchmark/inputs/cases.json').read_text())
    ordinary = [dict(id=c['id'], model=c['model'], T=c['T_K'], x=c['x'], branch=c['branch'],
                     states=[dict(index=s['index'], P=s['P_Pa'], oracle_rho=s['oracle_rho_mol_m3'],
                                  oracle_roots=s['oracle_stable_root_count'], floor=s['main_pressure_floor'])
                             for s in c['states']])
                for c in density['cases'] if c['id'] in ('A-propane-liquid', 'A-propane-vapor', 'A-C2-liquid', 'A-methane-ethane-vapor')]
    loop = [dict(id=f'A-Q02/{T}-{branch}', T=T, branch=branch, pressures=[P * (1 + j * 1e-8) for j in (0, -2, -1, 1, 2)])
            for T, P in ((308.8, 5147500.), (308.0, 5057915.461281232)) for branch in ('vapor', 'liquid')]
    return dict(ordinary=ordinary, ethane_loop=loop, actions=dict(directions=[[1., 0.], [0., 1.]], orders=[1, 2]))


def main():
    out = HERE / 'inputs/inputs.json'
    if out.exists():
        raise FileExistsError(f'{out} is frozen')
    comps = epc_components()
    molalities = [10 ** (math.log10(.001) + j * (math.log10(6) - math.log10(.001)) / 29) for j in range(30)]
    rng = random.Random(578)
    multiset = [275., 325., 375.] * 100
    rng.shuffle(multiset)
    inputs = dict(
        sources={str(p): sha(p) for p in (AUDIT / 'opus/probes/epc.py', AUDIT / 'opus/probes/p07/engine.py',
                                          AUDIT / 'opus/probes/p03/engine.py', AUDIT / 'opus/probes/p04/engine.py',
                                          AUDIT / 'opus/probes/p10/engine.py', AUDIT / 'astra/probes/engine_probes.py',
                                          AUDIT / 'astra/probes/extended_engine_probes.py', AUDIT / 'sol/probes/engine_cases.py',
                                          ROOT / 'analyses/2026-density-root-benchmark/inputs/cases.json',
                                          ROOT / 'analyses/2026-competitor-cross-check/inputs/cases.json')},
        components=comps,
        models={'propane': ['PROPANE'], 'methanol': ['METHANOL'], 'acetone-polar': ['ACETONE_POLAR'],
                'methane-ethane': ['METHANE', 'ETHANE'], 'methanol-hexane': ['METHANOL', 'NHEXANE'],
                'gross-methane-ethane': ['METHANE', 'ETHANE'], 'ethane': ['ETHANE']},
        tolerances=dict(action_abs=2e-8, action_rel=2e-5, helmholtz_scaled=1e-6, property_rel=1e-8,
                        floors=dict(P=1., h=1., cv=1., ar=1., lnphi=1.), density_rel=1e-10, same_image=1e-10,
                        equilibrium_rel=1e-7, composition_abs=1e-8),
        t1=t1(comps), t2=t2(),
        t3=['C1-300', 'C2-0.5', 'methane-bubble'],
        t4=dict(methane_ethane=dict(x=[.5, .5], points=[[200., 3e6], [220., 4e6], [240., 5e6], [190., 1e6], [250., 6.5e6], [200., .5e6]]),
                methanol_hexane=dict(x=[.5, .5], kij=.03, P=1e5, temperatures=[250., 270., 290., 300., 310.])),
        t6=dict(T=298.15, P=1e5, request='C3-0.88sigma', molalities=molalities),
        t7=['case-6/figure-2-point-7/cold', 'case-6/figure-2-point-7/accepted', 'case-7/flash-5918.5286-Pa/cold',
            'case-7/flash-5918.5286-Pa/accepted', 'khudaida/sequential-8', 'mea-nominally-unloaded/353.15-K',
            'mea-f1/residual-jacobian', 'mea-cp/hilliard-353.15-K-0.358'],
        t11=dict(temperatures=[275., 325., 375.], rho=100., x=[.4, .6], intervals=[[250., 300.], [300., 350.], [350., 400.]],
                 sequences={'two-repeated': [275., 325.] * 150, 'three-repeated': [275., 325., 375.] * 100,
                            'three-permuted-578': multiset, 'three-partitioned': sorted(multiset)}),
    )
    out.write_text(json.dumps(inputs, indent=1) + '\n')
    t1rows = inputs['t1']
    print(len(t1rows), 'T1 states;', sum(len(r['requests']) for r in t1rows), 'action requests;',
          sum(1 for r in t1rows for q in r['requests'] if q['inadmissible_chart']), 'inadmissible;',
          sum(1 for r in t1rows for q in r['requests'] for lv in q['levels'] if not lv['admissible']), 'nonadmissible step levels')


if __name__ == '__main__':
    main()

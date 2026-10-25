"""Persistent guarded Engine worker for #578: one JSON request per stdin line, one JSON reply per line.

All thermodynamic work is native Engine calls on frozen declarations (inputs/inputs.json and the
repository declaration owners). Python times calls and records returns; finite differences only
verify exact derivatives, outside timing.
"""
import ast
import importlib.util
from importlib.machinery import SourceFileLoader
import itertools
import json
import math
import os
from pathlib import Path
import sys
import time
import traceback

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
AUDIT = Path('/home/tnnrpolley21/Workspaces/Engineering/ePC-SAFT/build/audit-2026-10-15')
INPUTS = json.loads((HERE / 'inputs/inputs.json').read_text())
T1 = {r['id']: r for r in INPUTS['t1']}
CACHE = {}


def load(path, name, skip=()):
    """Load a repository/audit module, dropping top-level statements named in skip (AST node types or names)."""
    path = Path(path)
    if not skip:
        loader = SourceFileLoader(name, str(path))
        spec = importlib.util.spec_from_loader(name, loader)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        loader.exec_module(module)
        return module
    tree = ast.parse(path.read_text())
    tree.body = [n for n in tree.body if not skip(n)]
    module = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, loader=None))
    module.__file__ = str(path)
    sys.modules[name] = module
    exec(compile(tree, str(path), 'exec'), module.__dict__)
    return module


def clean(x):
    if isinstance(x, float):
        return x if math.isfinite(x) else None
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if hasattr(x, 'name') and type(x).__module__.startswith('epcsaft'):
        return x.name
    return x


def cached(key, build):
    if key not in CACHE:
        CACHE[key] = build()
    return CACHE[key]


def counters(model):
    c = model.counters()
    return dict(forward=c.forward, reverse=c.reverse, hessian_calls=c.hessian_calls, recordings=c.recordings)


# ---------------------------------------------------------------- declarations

def epc():
    return cached('epc', lambda: load(AUDIT / 'opus/probes/epc.py', 'epc578'))


def p07_document():
    # P07's top level declares the document, then evaluates states; keep only the declarations.
    def skip(n):
        evaluates = isinstance(n, ast.For) and {getattr(t, 'id', None) for t in ast.walk(n.target)} & {'rho', 'r'}
        return evaluates or isinstance(n, ast.Expr) or isinstance(n, ast.Assign) and getattr(n.targets[0], 'id', '') in ('mix', 'rows')
    sys.path.insert(0, str(AUDIT / 'opus/probes'))
    return load(AUDIT / 'opus/probes/p07/engine.py', 'p07_578', skip).doc


def density_declarations():
    # The density benchmark declarations run a rebuild probe at import; the guard forbids it, so drop it.
    return cached('density', lambda: load(ROOT / 'analyses/2026-density-root-benchmark/scripts/declarations.py', 'density578',
                                          lambda n: isinstance(n, ast.For)))


def model(key):
    import epcsaft

    def build():
        if key in ('gross-methane-ethane', 'ethane'):
            components = ('methane', 'ethane') if key == 'gross-methane-ethane' else ('ethane',)
            p = epcsaft.Parameters.from_json(ROOT / 'data/packets/gross-2001-methane-ethane/3/parameters/parameters.json',
                                             components=components)
        elif key in INPUTS['models']:
            comps = [dict(INPUTS['components'][n]) for n in INPUTS['models'][key]]
            for c in comps:
                if 'assoc' in c:
                    c['assoc'] = tuple(c['assoc'])
            p = epc().parameters([epc().component(**c) for c in comps], polar=any('mu' in c for c in comps))
        elif key == 'water-nacl':
            p = epcsaft.Parameters.from_mapping(p07_document())
        elif key == 'mea':
            return density_declarations().models()['MEA']
        elif key == 'associating-binary':
            sys.path.insert(0, str(ROOT / 'engine/tests'))
            from native_inputs import associating_binary
            p = associating_binary()
        elif key in ('C1', 'C2', 'methane-ethane'):
            return density_declarations().models()[key]
        else:
            raise KeyError(key)
        return (epcsaft.Pure if len(p.component_ids) == 1 else epcsaft.Mixture)(p)
    return cached(('model', key), build)


# ---------------------------------------------------------------- T7 / T3: key solves and represented VLE

def history_worker():
    return cached('history', lambda: load(ROOT / 'analyses/2026-engine-performance-history/scripts/worker.py', 'history578'))


def key_tool():
    return cached('key', lambda: load(ROOT / 'tools/time-key-solves', 'key578'))


def op_key(req):
    name = req['case']
    if name in key_tool().FROZEN:
        def frozen():
            result = key_tool().frozen_result(name, 'evaluation')
            row = result['rows'][0]
            return dict(seconds=result['seconds'], cpu_seconds=result['cpu_seconds'],
                        work={k: row[k] for k in key_tool().F1_COUNTS}, value={k: row[k] for k in key_tool().F1_VALUES})
        return frozen
    call = history_worker().make(name, ROOT)

    def solve():
        result = call()
        rows = result['rows']
        work = {k: sum(r[k] for r in rows) for k in ('iterations', 'evaluations', 'hessian_calls')}
        return dict(seconds=result['seconds'], cpu_seconds=result['cpu_seconds'], work=work, value=rows)
    return solve


# ---------------------------------------------------------------- T1: fixed-density properties and actions

def t1_state(sid):
    r = T1[sid]
    return model(r['model']), r['T'], r['rho'], r['x']


def bundle(s):
    return dict(P=s.pressure, h_res=s.residual_enthalpy, cv_res=s.residual_isochoric_heat_capacity,
                lnphi=s.log_fugacity_coefficient, a_res=s.residual_helmholtz)


def op_t1_bundle(req):
    m, T, rho, x = t1_state(req['state'])
    return lambda: dict(value=bundle(m.state(T, rho=rho, x=x)))


def directions(sid, names):
    from epcsaft import DensityActionDirection
    d = T1[sid]['directions']
    return [DensityActionDirection(d[n][0], d[n][1], d[n][2]) for n in names]


def action_record(batch):
    r = batch.results[0]
    d = batch.diagnostics
    return dict(value=r.value, action=r.action, status=r.status.name, energy_order=r.energy_order,
                maximum_energy_order=batch.maximum_energy_order,
                diagnostics={k: getattr(d, k) for k in ('site_rank', 'site_factorizations', 'site_solves', 'association_residual',
                                                         'thermal_evaluations', 'recordings', 'directional_paths', 'forward_sweeps',
                                                         'pressure_factorizations', 'pressure_solves')})


def op_t1_action(req):
    from epcsaft import PropertyObservable
    m, T, rho, x = t1_state(req['state'])
    ds = directions(req['state'], req['names'])
    obs = [PropertyObservable.ResidualEnthalpy]
    return lambda: dict(value=action_record(m.density_actions(T, rho, x=x, observables=obs, directions=ds)))


def op_t1_neutral(req):
    from epcsaft import _core
    m, T, rho, _ = t1_state(req['state'])
    order = getattr(_core.DerivativeOrder, req['order'])
    prop = _core.NeutralProperty.ResidualHelmholtz

    def call():
        r = m._model.neutral_property(T, rho, prop, order)
        return dict(value=dict(value=r.value, gradient=list(r.gradient), hessian=list(r.hessian)))
    return call


def stencil(f, T, rho, x, dirs, h):
    total = 0.0
    for signs in itertools.product((-1.0, 1.0), repeat=len(dirs)):
        pT = T + sum(s * hj * d[0] for s, hj, d in zip(signs, h, dirs))
        prho = rho + sum(s * hj * d[1] for s, hj, d in zip(signs, h, dirs))
        px = [xi + sum(s * hj * d[2][i] for s, hj, d in zip(signs, h, dirs)) for i, xi in enumerate(x)]
        total += math.prod(signs) * f(pT, prho, px)
    return total / (2 ** len(dirs) * math.prod(h))


def op_t1_verify(req):
    """Exact actions against frozen-step stencils (design R1); untimed."""
    from epcsaft import PropertyObservable, _core
    sid = req['state']
    m, T, rho, x = t1_state(sid)
    tol = INPUTS['tolerances']
    raw = T1[sid]['directions']

    def call():
        rows = []
        for q in T1[sid]['requests']:
            row = dict(order=q['order'], names=q['names'])
            try:
                batch = m.density_actions(T, rho, x=x, observables=[PropertyObservable.ResidualEnthalpy],
                                          directions=directions(sid, q['names']))
                rec = action_record(batch)
                row.update(action=rec['action'], status=rec['status'])
            except Exception as error:
                row.update(action=None, status=f'refused: {type(error).__name__}: {error}')
            if q['inadmissible_chart'] or 'charge balance' in str(row['status']):
                row.update(verdict='inadmissible chart', cause=q['inadmissible_chart'] or row['status'])
            elif row['action'] is None:
                row.update(verdict='native unavailable', cause=row['status'])
            else:
                dirs = [raw[n] for n in q['names']]
                levels = []
                for lv in q['levels']:
                    try:
                        fd = stencil(lambda a, b, c: m.state(a, rho=b, x=c).residual_enthalpy, T, rho, x, dirs, lv['h'])
                    except Exception as error:
                        levels.append(dict(h=lv['h'], fd=None, cause=f'{type(error).__name__}: {error}'))
                        continue
                    err = abs(row['action'] - fd)
                    allowed = max(tol['action_abs'], tol['action_rel'] * max(abs(row['action']), abs(fd)))
                    levels.append(dict(h=lv['h'], fd=fd, error=err, allowed=allowed, passed=err <= allowed))
                row['levels'] = levels
                row['verdict'] = 'verified' if levels and all(lv.get('passed') for lv in levels) else 'unresolved verification'
            rows.append(row)
        if T1[sid]['neutral_property']:
            exact = m._model.neutral_property(T, rho, _core.NeutralProperty.ResidualHelmholtz, _core.DerivativeOrder.Hessian)
            value = lambda a, b: m._model.neutral_property(a, b, _core.NeutralProperty.ResidualHelmholtz, _core.DerivativeOrder.Value).value
            grad = lambda a, b: m._model.neutral_property(a, b, _core.NeutralProperty.ResidualHelmholtz, _core.DerivativeOrder.Gradient).gradient
            for scale in (1., .5):
                hT, hr = min(.02, T / 4) * scale, min(1e-4 * rho, rho / 4) * scale
                fd_g = [(value(T + hT, rho) - value(T - hT, rho)) / (2 * hT), (value(T, rho + hr) - value(T, rho - hr)) / (2 * hr)]
                gT = [(a - b) / (2 * hT) for a, b in zip(grad(T + hT, rho), grad(T - hT, rho))]
                gr = [(a - b) / (2 * hr) for a, b in zip(grad(T, rho + hr), grad(T, rho - hr))]
                fd_h = [gT[0], gT[1], gr[0], gr[1]]
                for kind, ex, fd in (('gradient', list(exact.gradient), fd_g), ('hessian', list(exact.hessian), fd_h)):
                    errs = [abs(a - b) / max(1, abs(a), abs(b)) for a, b in zip(ex, fd)]
                    rows.append(dict(order=kind, names=['T', 'rho'], step_scale=scale, exact=ex, fd=fd, scaled_error=max(errs),
                                     verdict='verified' if max(errs) <= tol['helmholtz_scaled'] else 'unresolved verification'))
        return dict(value=rows)
    return call


# ---------------------------------------------------------------- T2: density closure

def t2_case(case_id):
    return next(c for c in INPUTS['t2']['ordinary'] if c['id'] == case_id)


def density_record(s):
    c = s.density_convergence
    return dict(rho=s.molar_density, P=s.pressure, dPdrho=s.pressure_density_derivative, branch=s.density_branch.name,
                roots=s.stable_root_count,
                convergence=None if c is None else dict(residual=c.residual, iterations=c.iterations, evaluations=c.evaluations,
                                                        stop=c.stop.name))


def op_t2_sequence(req):
    """32 changing pressures on one branch: cold, anchored on the previous accepted root, or fixed density."""
    c = t2_case(req['case'])
    m = model(c['model'])
    mode = req['mode']
    T, x, branch = c['T'], c['x'], c['branch']
    roots = req.get('roots')

    def call():
        rows, ns = [], 0
        previous = None
        for i, st in enumerate(c['states']):
            t = time.perf_counter_ns()
            if mode == 'fixed-density':
                s = m.state(T, rho=roots[i], x=x)
            else:
                s = m.state(T, P=st['P'], x=x, phase=branch, anchor=previous if mode == 'anchored' else None)
            ns += time.perf_counter_ns() - t
            rows.append(density_record(s))
            previous = s
        return dict(seconds=ns / 1e9, per_call=len(rows), value=rows)
    return call


def op_t2_loop(req):
    loop = next(c for c in INPUTS['t2']['ethane_loop'] if c['id'] == req['case'])
    m = model('ethane')
    roots = req.get('roots')

    def call():
        rows, ns, previous = [], 0, None
        for i, P in enumerate(loop['pressures']):
            t = time.perf_counter_ns()
            if req['mode'] == 'fixed-density':
                s = m.state(loop['T'], rho=roots[i])
            else:
                s = m.state(loop['T'], P=P, phase=loop['branch'], anchor=previous if req['mode'] == 'anchored' else None)
            ns += time.perf_counter_ns() - t
            rows.append(density_record(s))
            previous = s
        return dict(seconds=ns / 1e9, per_call=len(rows), value=rows)
    return call


def op_t2_actions(req):
    from epcsaft import PressureActionDirection, PropertyObservable
    loop = next(c for c in INPUTS['t2']['ethane_loop'] if c['id'] == req['case'])
    m = model('ethane')
    ds = [PressureActionDirection(*INPUTS['t2']['actions']['directions'][j]) for j in req['combo']]

    def call():
        b = m.pressure_actions(loop['T'], loop['pressures'][0], phase=loop['branch'],
                               observables=[PropertyObservable.ResidualEnthalpy], directions=ds)
        rec = action_record(b)
        rec.update(rho=b.molar_density, density_action=list(b.density_action), density_mixed_action=b.density_mixed_action,
                   density_status=b.density_status.name)
        return dict(value=rec)
    return call


def op_t2_diagnostics(req):
    """Native root diagnostics for the narrow-loop states, once, outside timing."""
    from epcsaft import _core
    m = model('ethane')

    def call():
        rows = []
        for loop in INPUTS['t2']['ethane_loop']:
            if loop['branch'] != 'vapor':
                continue
            T, P = loop['T'], loop['pressures'][0]
            hi = m.molar_density_upper_bound(T)
            nodes = [hi * i / 4096 for i in range(1, 4096)]
            d = m.density_diagnostics(T, P, nodes=nodes, density_tolerance=1e-10, log_tolerance=1e-12, max_evaluations=3 * 4096,
                                      max_seconds=60, adaptive=False, finite_upper=hi, evaluation_order=_core.DerivativeOrder.PressureSlope)
            rows.append(dict(T=T, P=P, roots=[dict(rho=r.molar_density, residual=r.residual, slope=r.pressure_density_derivative,
                                                   sign=r.mechanical_sign) for r in d.pressure_roots],
                             spinodals=[s.molar_density for s in d.spinodals], unresolved=len(d.unresolved),
                             brackets=[(b.lower, b.upper) for b in d.declared_brackets], evaluations=d.evaluations, stop=d.stop))
        return dict(value=rows)
    return call


# ---------------------------------------------------------------- T4: unknown phase count

def t4_model(system):
    if system == 'methane_ethane':
        return model('methane-ethane')
    kij = INPUTS['t4']['methanol_hexane']['kij']

    def build():
        import epcsaft
        comps = [dict(INPUTS['components'][n]) for n in ('METHANOL', 'NHEXANE')]
        comps[0]['assoc'] = tuple(comps[0]['assoc'])
        return epcsaft.Mixture(epc().parameters([epc().component(**c) for c in comps], kij={('methanol', 'n-hexane'): kij}))
    return cached(('model', 'methanol-hexane-kij'), build)


def t4_points(system):
    d = INPUTS['t4'][system]
    return d['points'] if system == 'methane_ethane' else [[T, d['P']] for T in d['temperatures']]


def equilibrium_record(r):
    out = dict(success=r.success, criterion_met=r.requested_tolerance_met, phases=r.phase_count, amounts=list(r.amounts),
               mole_fractions=list(r.mole_fractions), molar_densities=list(r.molar_densities), pressure=[] if r.pressure is None else [r.pressure],
               temperature=[] if r.temperature is None else [r.temperature], iterations=r.iterations, evaluations=r.evaluations, hessian_calls=r.hessian_calls,
               max_residual=max(map(abs, r.residuals)) if r.residuals else None, message=r.message[:300])
    if r.discovery is not None:
        d = r.discovery
        out['discovery'] = dict(work=dict(d.work), rounds=[dict(r.work) for r in d.rounds], declarations=len(d.declarations),
                                roles=list(d.declaration_roles), attempts=len(d.declaration_attempt_coordinates))
    return out


def op_t4(req):
    """Discover at each point, the accepted discovered topology solved cold, or that topology along the sequence."""
    from epcsaft.equilibrium import Discover, MoleFractions, Problem, solve_equilibrium
    system, mode = req['system'], req['mode']
    m = t4_model(system)
    ids = m.component_ids
    points = t4_points(system)
    if 'point' in req:
        points = [points[req['point']]]
    z = MoleFractions(dict(zip(ids, (.5, .5))))
    declared = None
    if mode != 'discover':
        declared = [solve_equilibrium(m, Problem(phases=Discover(), T=T, P=P, feed=z)).discovery.declaration for T, P in points]

    def call():
        rows, ns, previous = [], 0, None
        for i, (T, P) in enumerate(points):
            t = time.perf_counter_ns()
            try:
                if mode == 'discover':
                    r = solve_equilibrium(m, Problem(phases=Discover(), T=T, P=P, feed=z))
                else:
                    start = previous if mode == 'declared-sequence' and previous is not None and previous.phase_count == len(declared[i].phases) else None
                    r = solve_equilibrium(m, declared[i], start)
                ns += time.perf_counter_ns() - t
                rows.append(dict(T=T, P=P, **equilibrium_record(r)))
                previous = r if r.success else None
            except Exception as error:
                ns += time.perf_counter_ns() - t
                rows.append(dict(T=T, P=P, error=f'{type(error).__name__}: {error}'))
                previous = None
        return dict(seconds=ns / 1e9, per_call=len(rows), value=rows)
    return call


# ---------------------------------------------------------------- T6: one reference, 30-point activity curve

def t6_model():
    import epcsaft
    cc = cached('cc', lambda: load(ROOT / 'analyses/2026-competitor-cross-check/scripts/engine.py', 'cc578'))
    inputs = json.loads((ROOT / 'analyses/2026-competitor-cross-check/inputs/cases.json').read_text())
    q = next(r for r in inputs['requests'] if r['id'] == INPUTS['t6']['request'] + '-1')
    spec = inputs['models'][q['case']]
    ids = [c['id'] for c in spec['components']]
    return epcsaft.Mixture(epcsaft.Parameters.from_mapping(cc.document(spec, q), components=ids)), ids


def op_t6(req):
    """Each call builds a fresh model (untimed) so no arm inherits a warm reference from an earlier block."""
    T, P, molalities = INPUTS['t6']['T'], INPUTS['t6']['P'], INPUTS['t6']['molalities']

    def call():
        m, ids = t6_model()
        reference = lambda: m.reference(T, P, solvent_fractions=[1, 0, 0], ions=(ids[1], ids[2]), stoichiometry=(1, 1))
        rows, ref_ns, curve_ns = [], 0, 0
        if req['mode'] == 'reused':
            t = time.perf_counter_ns()
            ref = reference()
            ref_ns = time.perf_counter_ns() - t
            for mj in molalities:
                t = time.perf_counter_ns()
                a = ref.activity(mj)
                curve_ns += time.perf_counter_ns() - t
                rows.append(dict(m=mj, gamma=a.mean_ionic_molality_scale, phi=a.osmotic_coefficient))
        else:
            for mj in molalities:
                t = time.perf_counter_ns()
                ref = reference()
                a = ref.activity(mj)
                curve_ns += time.perf_counter_ns() - t
                rows.append(dict(m=mj, gamma=a.mean_ionic_molality_scale, phi=a.osmotic_coefficient))
        c = counters(m)
        return dict(seconds=(ref_ns + curve_ns) / 1e9, reference_seconds=ref_ns / 1e9, curve_seconds=curve_ns / 1e9,
                    counters=c, value=dict(points=rows, reference_error=ref.convergence_error,
                                           reference_molality=ref.molality))
    return call


# ---------------------------------------------------------------- T8: solved-state actions on one central solve

def t8_case(name):
    from epcsaft import equilibrium as q

    def build():
        K = q.SolvedStateObservableKind
        if name == 'isomerisation':
            sys.path[:0] = [str(ROOT / 'engine'), str(ROOT / 'engine/tests')]
            test = load(ROOT / 'engine/tests/equilibrium/test_native_transport.py', 'transport578')
            fn = next(n for n in ast.parse(Path(test.__file__).read_text()).body
                      if isinstance(n, ast.FunctionDef) and n.name == 'test_thermal_actions_of_the_logistic_isomerisation_are_analytic')
            fn.body = [n for n in fn.body if not isinstance(n, (ast.For, ast.Expr))] + ast.parse('return model, problem, analytic').body
            fn.name, fn.args.args = 'declaration578', []
            exec(compile(ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[])), test.__file__, 'exec'), test.__dict__)
            m, problem, analytic = test.declaration578()
            obs = [q.SolvedStateObservable(K.PhaseMoleFraction, 0, i) for i in range(2)]
            dirs = [q.SolvedStateActionDirection(1., 0., [0., 0.], []), q.SolvedStateActionDirection(0., 1e5, [0., 0.], []),
                    q.SolvedStateActionDirection(0., 0., [.01, -.01], [])]
            return dict(model=m, problem=problem, observables=obs, directions=dirs, analytic=analytic)
        if name == 'bubble-C2':
            hh = load(ROOT / 'analyses/2026-head-to-head-comparisons/scripts/h2h_engine.py', 'h2h578')
            tree = ast.parse(Path(hh.__file__).read_text())
            fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'sensitivity_case')
            fn.body = fn.body[:next(i for i, n in enumerate(fn.body) if isinstance(n, ast.FunctionDef) and n.name == 'exact')]
            fn.body += ast.parse('return model, prob, central, observables, [d for _, d in directions]').body
            exec(compile(ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[])), hh.__file__, 'exec'), hh.__dict__)
            inputs = json.loads((ROOT / 'analyses/2026-competitor-cross-check/inputs/cases.json').read_text())
            m, problem, central, obs, dirs = hh.sensitivity_case(inputs, 'C2-0.5')
            return dict(model=m, problem=problem, observables=obs, directions=dirs)
        if name == 'ionic-lle':
            import epcsaft
            names = ('water', 'isobutanol', 'sodium-cation', 'chloride-anion')
            m = epcsaft.Mixture(epcsaft.Parameters.from_json(
                ROOT / 'data/packets/khudaida-2026-figure-2-electrolyte-lle/7/parameters/parameters.json', components=names))
            problem = q.Problem(T=293.15, P=1e5, row_tolerance=1e-10, feed=q.Amounts(dict(zip(names, (.59, .39, .01, .01)))),
                                phases=[q.Phase('Lw', kind='liquid', composition_guess=(.97, .01, .01, .01)),
                                        q.Phase('Lo', kind='liquid', composition_guess=(.3, .698, .001, .001))])
            obs = [q.SolvedStateObservable(kind, j, i) for kind in (K.PhaseComponentAmount, K.PhaseMoleFraction)
                   for j in range(2) for i in range(4)]
            dirs = [q.SolvedStateActionDirection(1., 0., [0.] * 4, []), q.SolvedStateActionDirection(0., 1e5, [0.] * 4, []),
                    q.SolvedStateActionDirection(0., 0., [-.002, 0., .001, .001], []),
                    q.SolvedStateActionDirection(0., 0., [.01, -.01, 0., 0.], [])]
            return dict(model=m, problem=problem, observables=obs, directions=dirs)
        raise KeyError(name)
    case = cached(('t8', name), build)
    if 'central' not in case:
        case['central'] = q.solve_equilibrium(case['model'], case['problem'])
        case['compiled'] = q.compile_problem(case['model'], case['problem'])
    return case


def solved_record(batch):
    d = batch.diagnostics
    return dict(actions=[r.action for r in batch.results], status=[r.status.name for r in batch.results],
                order=[r.order for r in batch.results],
                diagnostics={k: getattr(d, k) for k in ('equation_count', 'state_dimension', 'maximum_energy_order', 'state_factorizations',
                                                         'linear_solves', 'rank', 'condition_estimate', 'maximum_linear_residual',
                                                         'central_residual', 'central_floor')})


def op_t8(req):
    from epcsaft import equilibrium as q
    c = t8_case(req['case'])
    dirs = c['directions']
    if req['mode'] == 'first-batch':
        chosen = dirs if req['n'] == 'all' else dirs[:req['n']]
        return lambda: dict(value=[solved_record(b) for b in q.solved_state_first_actions(c['compiled'], c['central'], c['observables'], chosen)])
    if req['mode'] == 'first-separate':
        chosen = dirs if req['n'] == 'all' else dirs[:req['n']]
        return lambda: dict(value=[solved_record(q.solved_state_actions(c['compiled'], c['central'], q.SolvedStateActionRequest(c['observables'], [d])))
                                   for d in chosen])
    request = q.SolvedStateActionRequest(c['observables'], [dirs[j] for j in req['combo']])
    return lambda: dict(value=solved_record(q.solved_state_actions(c['compiled'], c['central'], request)))


def op_t8_verify(req):
    """Outside timing: batch equals separate first actions; isomerisation A1-A3 against its analytic law; every order-k
    action with an input (T, P or feed) direction against centered differences of the exact order-(k-1) action on
    states re-solved at y +/- h d (design T8 targets 1e-6/1e-5/1e-4 relative, 1e-10 absolute)."""
    import dataclasses
    import itertools as it
    from epcsaft import equilibrium as q
    c = t8_case(req['case'])
    dirs, obs = c['directions'], c['observables']
    targets = {1: 1e-6, 2: 1e-5, 3: 1e-4}

    def moved(d, h):
        p = c['problem']
        changes = {}
        if d.temperature_k:
            changes['T'] = p.T + h * d.temperature_k
        if d.pressure_pa:
            changes['P'] = p.P + h * d.pressure_pa
        if any(d.feed_amounts_mol):
            values = p.feed.values
            changes['feed'] = type(p.feed)({k: v + h * dv for (k, v), dv in zip(values.items(), d.feed_amounts_mol)})
        return dataclasses.replace(p, **changes)

    def exact(problem, combo):
        central = q.solve_equilibrium(c['model'], problem)
        if not central.success:
            raise RuntimeError(central.message)
        rest = [dirs[j] for j in combo] or [dirs[0]]
        b = q.solved_state_actions(q.compile_problem(c['model'], problem), central, q.SolvedStateActionRequest(obs, rest))
        return [r.value if not combo else r.action for r in b.results]

    def call():
        rows = []
        batch = q.solved_state_first_actions(c['compiled'], c['central'], obs, dirs)
        for j, d in enumerate(dirs):
            sep = q.solved_state_actions(c['compiled'], c['central'], q.SolvedStateActionRequest(obs, [d]))
            a, b = [r.action for r in batch[j].results], [r.action for r in sep.results]
            scale = max([abs(v) for v in b if v is not None] + [1e-300])
            rows.append(dict(check='batch-equals-separate', combo=[j], scaled_difference=max(
                abs(x - y) / scale for x, y in zip(a, b) if x is not None and y is not None)))
        if 'analytic' in c:
            for k in (1, 2, 3):
                b = q.solved_state_actions(c['compiled'], c['central'], q.SolvedStateActionRequest(obs[1:], [dirs[0]] * k))
                value, ref = b.results[0].action, c['analytic'][k - 1]
                rows.append(dict(check='analytic', combo=[0] * k, action=value, reference=ref, error=abs(value - ref),
                                 allowed=2e-10 * abs(ref) + 1e-14, passed=abs(value - ref) <= 2e-10 * abs(ref) + 1e-14))
        inputs = [j for j, d in enumerate(dirs) if not any(d.active_parameters)]
        for k in (1, 2, 3):
            for combo in it.combinations_with_replacement(range(len(dirs)), k):
                lead = next((j for j in combo if j in inputs), None)
                b = q.solved_state_actions(c['compiled'], c['central'], q.SolvedStateActionRequest(obs, [dirs[j] for j in combo]))
                actions = [r.action for r in b.results]
                row = dict(check='re-solve', combo=list(combo), actions=actions, status=[r.status.name for r in b.results])
                if lead is None:
                    row.update(verdict='no independent replay (parameter directions only)')
                    rows.append(row)
                    continue
                rest = list(combo)
                rest.remove(lead)
                levels = []
                for h in (1e-2, 5e-3):
                    try:
                        hi, lo = exact(moved(dirs[lead], h), rest), exact(moved(dirs[lead], -h), rest)
                    except Exception as error:
                        levels.append(dict(h=h, cause=f'{type(error).__name__}: {error}'))
                        continue
                    fd = [(x - y) / (2 * h) for x, y in zip(hi, lo)]
                    errs = [abs(a - f) - max(1e-10, targets[k] * max(abs(a), abs(f))) for a, f in zip(actions, fd) if a is not None]
                    levels.append(dict(h=h, fd=fd, worst_excess=max(errs), passed=max(errs) <= 0))
                row.update(lead=lead, levels=levels, verdict='verified' if any(l.get('passed') for l in levels) else 'unresolved verification')
                rows.append(row)
        return dict(value=rows)
    return call


# ---------------------------------------------------------------- T9: fitting throughput

def t9_declaration(family):
    from epcsaft import _core, regression as rg

    def build():
        if family == 'F1':
            tool = key_tool()
            fx = tool.FROZEN[tool.F1]
            data = fx['data']
            # Reuse the key-solve decoder by evaluating its nested decode() on this packet.
            fn = next(n for n in ast.parse(Path(tool.__file__).read_text()).body if isinstance(n, ast.FunctionDef) and n.name == 'frozen_result')
            dec = next(n for n in fn.body if isinstance(n, ast.FunctionDef) and n.name == 'decode')
            ns = dict(_core=_core)
            exec(compile(ast.fix_missing_locations(ast.Module(body=[dec], type_ignores=[])), tool.__file__, 'exec'), ns)
            import epcsaft
            params = epcsaft.Parameters.from_mapping(data['parameters'])
            coords, observations, thermal = (ns['decode'](data[k]) for k in ('coordinates', 'observations', 'thermochemistry'))
            return dict(parameters=params._native, coordinates=coords, observations=observations, thermal=thermal,
                        weights=data['weights'], holds=[])
        if family == 'propane':
            import csv
            from epcsaft.equilibrium import Amounts, Free, Phase, Pinned, Problem
            path = ROOT / 'data/packets/nist-srd69-pure-saturation-properties/1/observations/propane-reference-saturation.csv'
            rows = [r for r in csv.DictReader(path.open()) if 150.0 <= float(r['T_K']) <= 360.0]
            start = dict(INPUTS['components']['PROPANE'], m=2.2, sigma=3.5, eps=200.0)
            p0 = epc().parameters([epc().component(**start)])
            obs = []
            for r in rows:
                sat = Problem(phases=[Phase('L', kind='liquid'), Phase('V', kind='vapor', amount=Pinned(0.0))], T=float(r['T_K']),
                              P=Free(float(r['p_sat_Pa'])), feed=Amounts({'propane': 1.0}))
                obs.append(rg.observation(p0, 'pressure', sat, observed=float(r['p_sat_Pa']), form='log_ratio'))
                obs.append(rg.observation(p0, 'phase_density', sat, observed=float(r['rho_sat_liq_kg_m3']), form='log_ratio',
                                          mass_basis=True, phase='L'))
            coords = [rg.coordinate(p0, [('segment_count', 'propane', 1.0)], origin=0.0, scale=1.0, bounds=(1.0, 5.0)),
                      rg.coordinate(p0, [('segment_diameter', 'propane', 1.0)], origin=0.0, scale=1.0, bounds=(2.5, 5.0)),
                      rg.coordinate(p0, [('dispersion_energy_over_k', 'propane', 1.0)], origin=0.0, scale=100.0, bounds=(100.0, 400.0))]
            return dict(parameters=p0._native, coordinates=coords, observations=obs, thermal=None, weights=[], holds=[],
                        rows=len(rows))
        if family == 'electrolyte':
            sys.path[:0] = [str(ROOT / 'engine'), str(ROOT / 'engine/tests')]
            test = load(ROOT / 'engine/tests/regression/test_fit.py', 'test_fit578')
            fn = next(n for n in ast.parse(Path(test.__file__).read_text()).body
                      if isinstance(n, ast.FunctionDef) and n.name == 'test_native_held_electrolyte_recovery_checks_and_columns')
            cut = next(i for i, n in enumerate(fn.body) if isinstance(n, ast.Assign) and getattr(n.targets[0], 'id', '') == 'controls')
            fn.body = fn.body[:cut] + ast.parse('return p, cs, observations, holds, [1 / (0.02 * math.log(10)) ** 2] * 6').body
            fn.name = 'declaration578'
            exec(compile(ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[])), test.__file__, 'exec'), test.__dict__)
            p, cs, observations, holds, weights = test.declaration578()
            return dict(parameters=p._native, coordinates=cs, observations=observations, thermal=None, weights=weights, holds=holds)
        raise KeyError(family)
    d = cached(('t9', family), build)
    if 'start' not in d:
        chart = _core.CoordinateSet()
        chart.coordinates = d['coordinates']
        d['start'] = list(chart.start(d['parameters'], d['observations'], d['thermal']) if d['thermal'] is not None
                          else chart.start(d['parameters'], d['observations']))
    return d


def op_t9_evaluate(req):
    from epcsaft import _core
    d = t9_declaration(req['family'])
    x = list(d['start'])
    x[0] += req.get('shift', 0.0)

    def call():
        t, c = time.perf_counter_ns(), time.process_time_ns()
        r = _core.evaluate_observations(d['parameters'], d['coordinates'], x, d['holds'], d['observations'], d['thermal'])
        wall, cpu = time.perf_counter_ns() - t, time.process_time_ns() - c
        return dict(seconds=wall / 1e9, cpu_seconds=cpu / 1e9,
                    value=dict(statuses=[s.name for s in r.statuses], predictions=list(r.predictions), gradient=list(r.gradient),
                               optimizer_gradient=list(r.optimizer_gradient)))
    return call


def op_t9_fit(req):
    """One native fit with the declared worker count; wall, parent and child CPU, and peak memory."""
    import resource
    from epcsaft import _core
    d = t9_declaration(req['family'])

    def call():
        controls = _core.FitControls(maximum_iterations=req['iterations'], maximum_elapsed_time_seconds=req['seconds'], workers=req['workers'])
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        t, c = time.perf_counter_ns(), time.process_time_ns()
        r = _core.fit(d['parameters'], d['thermal'], d['coordinates'], d['observations'], d['weights'], controls, d['holds'])
        wall, cpu = time.perf_counter_ns() - t, time.process_time_ns() - c
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        child = after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime
        return dict(seconds=wall / 1e9, cpu_seconds=cpu / 1e9, child_cpu_seconds=child,
                    parent_maxrss_kb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, child_maxrss_kb=after.ru_maxrss,
                    value=dict(status=r.status, usable=r.usable, message=r.message[:300], iterations=r.iterations,
                               residual_evaluations=r.residual_evaluations, jacobian_evaluations=r.jacobian_evaluations,
                               trial_failures=r.trial_failures if isinstance(r.trial_failures, int) else len(r.trial_failures), initial_cost=r.initial_cost, final_cost=r.final_cost,
                               optimizer=list(r.optimizer), physical=list(r.physical), predictions=list(r.predictions),
                               statuses=sorted({s.name for s in r.statuses}), covariance_rank=r.covariance.rank))
    return call

# ---------------------------------------------------------------- T10: ASL callbacks and record concurrency

PROPERTY, EQUILIBRIUM = 'epcsaft_pressure_property', 'epcsaft_equilibrium_phase_component_amount'
T10 = dict(property=dict(controls=[300., 1e6, .4, .6], neighbor=[300.01, 1e6 + 1., .4, .6], rows=[(r, 1) for r in range(4)]),
           equilibrium=dict(controls=[.7, .3], neighbor=[.7, .285], rows=[(p, c) for p in range(2) for c in range(2)]))


def nasa9():
    sys.path[:0] = [str(ROOT / 'engine'), str(ROOT / 'engine/tests')]
    from native_inputs import p3_forms
    return p3_forms()['nasa9']


def t10_record(function, k):
    """Property: methane/ethane with the T10 NASA9 record on [250, 400 + k] K. Equilibrium: check_records.py record B_k."""
    import epcsaft
    from epcsaft import IdealCorrelation, IdealInterval, ThermochemistryRecord, equilibrium as eq

    def build():
        if function == 'property':
            rec = ThermochemistryRecord(1e5, [[IdealInterval(250., 400. + k, True, True, IdealCorrelation(nasa9(), 1e5))] for _ in range(2)])
            p = epcsaft.Parameters.from_json(ROOT / 'data/packets/gross-2001-methane-ethane/3/parameters/parameters.json',
                                             components=('methane', 'ethane'))
            return epcsaft.eos.property_asl_record(epcsaft.Mixture(p, thermochemistry=rec))
        both = ('water', 'isobutanol')
        m = epcsaft.Mixture(epcsaft.Parameters.from_json(
            ROOT / 'data/packets/khudaida-2026-figure-2-electrolyte-lle/7/parameters/parameters.json', components=list(both)))
        problem = eq.Problem(phases=[eq.Phase('L1', kind='liquid', support=both, composition_guess=(0.46, 0.54)),
                                     eq.Phase('L2', kind='liquid', support=both, composition_guess=(0.99, 0.01))],
                             T=298.15, P=1.0e5, feed=eq.Amounts({'water': .7, 'isobutanol': .3}))
        options = eq.Continuation(initial_step=0.25, minimum_step=0.25, max_steps=4 + k, max_retries_per_step=1)
        return eq.asl_record(m, problem, options, [0.0, 0.0], [[1.0, 0.0], [0.0, 1.0]], [.7, .3], 1, 0.1)
    return cached(('record', function, k), build)


def asl_function(name):
    def build():
        from pyomo.environ import ExternalFunction
        image = json.loads((HERE / 'results/image.json').read_text())
        ext = ExternalFunction(library=str(Path(image['build_dir']) / 'libepcsaft_asl.so'), function=name)
        ext.load_library()
        return ext._known_functions[name][0]
    return cached(('asl', name), build)


def asl_node(function, record, controls, fgh):
    """One complete process node: every row of the record at these controls, as Pyomo would call them."""
    import ctypes
    from pyomo.core.base.external import _ARGLIST
    fn = asl_function(PROPERTY if function == 'property' else EQUILIBRIUM)
    n = len(controls)
    rows, per_row = [], []
    for a, b in T10[function]['rows']:
        packed = _ARGLIST([*controls, record, str(a), str(b)], fgh=fgh)
        t = time.perf_counter_ns()
        value = fn(ctypes.byref(packed))
        per_row.append((time.perf_counter_ns() - t) / 1e9)
        rows.append(dict(value=value, gradient=[packed.derivs[i] for i in range(n)] if fgh >= 1 else None,
                         hessian=[packed.hes[i] for i in range(n * (n + 1) // 2)] if fgh == 2 else None,
                         errmsg=packed.Errmsg.decode() if packed.Errmsg else ''))
    return rows, per_row


def op_t10_node(req):
    record = t10_record(req['function'], req.get('record', 0))
    spec = T10[req['function']]
    sequence = [spec['controls']] if req['controls'] == 'same' else [spec['controls'], spec['neighbor']]

    def call():
        out, per_row, ns = [], [], 0
        for controls in sequence:
            t = time.perf_counter_ns()
            rows, timings = asl_node(req['function'], record, controls, req['fgh'])
            ns += time.perf_counter_ns() - t
            out.append(rows)
            per_row.append(timings)
        return dict(seconds=ns / 1e9, per_call=len(sequence), per_row_seconds=per_row, value=out)
    return call


def op_t10_native(req):
    """The same nodes from native calls: pressure_actions mapped to (T, P, w) by the adapter's normalization chain,
    or the declared equilibrium solved and differentiated along the identity feed map."""
    import epcsaft
    from epcsaft import PressureActionDirection, PropertyObservable, ThermochemistryRecord, IdealInterval, IdealCorrelation
    from epcsaft import equilibrium as eq
    fgh = req['fgh']
    if req['function'] == 'property':
        rec = ThermochemistryRecord(1e5, [[IdealInterval(250., 400., True, True, IdealCorrelation(nasa9(), 1e5))] for _ in range(2)])
        m = epcsaft.Mixture(epcsaft.Parameters.from_json(ROOT / 'data/packets/gross-2001-methane-ethane/3/parameters/parameters.json',
                                                         components=('methane', 'ethane')), thermochemistry=rec)
        obs = [PropertyObservable.TotalEnthalpy, PropertyObservable.LogFugacityCoefficient]

        def node(q):
            T, P, w = q[0], q[1], q[2:]
            s = sum(w)
            x = [v / s for v in w]
            basis = [PressureActionDirection(1., 0., [0., 0.]), PressureActionDirection(0., 1., [0., 0.])]
            basis += [PressureActionDirection(0., 0., [(1. if i == j else 0.) - x[i] for i in range(2)]) for j in range(2)]
            n = len(q)

            def rows_of(b, k):
                return [b.density_action[0] if k == 1 else b.density_mixed_action] + [r.action for r in b.results]
            first = m.pressure_actions(T, P, x=x, phase='vapor', observables=obs, directions=[basis[0]])
            values = [first.molar_density] + [r.value for r in first.results]
            grads = [[0.] * n for _ in range(4)]
            for k in range(n if fgh >= 1 else 0):
                b = m.pressure_actions(T, P, x=x, phase='vapor', observables=obs, directions=[basis[k]])
                for row, a in enumerate(rows_of(b, 1)):
                    grads[row][k] = a / s ** (k >= 2)
            hes = [[] for _ in range(4)]
            for j in range(n if fgh == 2 else 0):
                for i in range(j + 1):
                    b = m.pressure_actions(T, P, x=x, phase='vapor', observables=obs, directions=[basis[i], basis[j]])
                    for row, a in enumerate(rows_of(b, 2)):
                        hes[row].append(a / s ** ((i >= 2) + (j >= 2)) - ((grads[row][i] + grads[row][j]) / s if i >= 2 else 0.))
            return [dict(value=values[r], gradient=grads[r] if fgh >= 1 else None, hessian=hes[r] if fgh == 2 else None) for r in range(4)]
    else:
        both = ('water', 'isobutanol')
        m = epcsaft.Mixture(epcsaft.Parameters.from_json(
            ROOT / 'data/packets/khudaida-2026-figure-2-electrolyte-lle/7/parameters/parameters.json', components=list(both)))
        K = eq.SolvedStateObservableKind
        obs = [eq.SolvedStateObservable(K.PhaseComponentAmount, p, c) for p, c in T10['equilibrium']['rows']]
        dirs = [eq.SolvedStateActionDirection(0., 0., [1., 0.], []), eq.SolvedStateActionDirection(0., 0., [0., 1.], [])]

        def node(q):
            problem = eq.Problem(phases=[eq.Phase('L1', kind='liquid', support=both, composition_guess=(0.46, 0.54)),
                                         eq.Phase('L2', kind='liquid', support=both, composition_guess=(0.99, 0.01))],
                                 T=298.15, P=1.0e5, feed=eq.Amounts({'water': q[0], 'isobutanol': q[1]}))
            central = eq.solve_equilibrium(m, problem)
            compiled = eq.compile_problem(m, problem)
            values = [central.amounts[p] * central.mole_fractions[2 * p + c] for p, c in T10['equilibrium']['rows']]
            grads = [[None, None] for _ in obs]
            if fgh >= 1:
                for k, b in enumerate(eq.solved_state_first_actions(compiled, central, obs, dirs)):
                    for row, r in enumerate(b.results):
                        grads[row][k] = r.action
            hes = [[] for _ in obs]
            for j in range(2 if fgh == 2 else 0):
                for i in range(j + 1):
                    b = eq.solved_state_actions(compiled, central, eq.SolvedStateActionRequest(obs, [dirs[i], dirs[j]]))
                    for row, r in enumerate(b.results):
                        hes[row].append(r.action)
            return [dict(value=values[r], gradient=grads[r] if fgh >= 1 else None, hessian=hes[r] if fgh == 2 else None,
                         success=central.success) for r in range(len(obs))]
    spec = T10[req['function']]
    sequence = [spec['controls']] if req['controls'] == 'same' else [spec['controls'], spec['neighbor']]

    def call():
        out, ns = [], 0
        for q in sequence:
            t = time.perf_counter_ns()
            out.append(node(q))
            ns += time.perf_counter_ns() - t
        return dict(seconds=ns / 1e9, per_call=len(sequence), value=out)
    return call


def op_t10_concurrent(req):
    """N caller threads, each with its own argument storage, alternating neighbor controls so every node is computed."""
    from concurrent.futures import ThreadPoolExecutor
    function, threads, calls = req['function'], req['threads'], req['calls']
    spec = T10[function]
    records = [t10_record(function, 0 if req['records'] == 'same' else k) for k in range(threads)]

    def caller(k):
        out = []
        for i in range(calls):
            rows, _ = asl_node(function, records[k], spec['controls'] if i % 2 == 0 else spec['neighbor'], req['fgh'])
            out.append(rows)
        return out

    def call():
        for k in range(threads):  # prepare every record before timing (one admitted slot each)
            asl_node(function, records[k], spec['controls'], 0)
        t = time.perf_counter_ns()
        with ThreadPoolExecutor(threads) as pool:
            results = list(pool.map(caller, range(threads)))
        wall = (time.perf_counter_ns() - t) / 1e9
        nodes = threads * calls
        return dict(seconds=wall, nodes=nodes, nodes_per_second=nodes / wall,
                    value=dict(first=[r[0] for r in results], last=[r[-1] for r in results],
                               errors=sum(1 for r in results for node in r for row in node if row['errmsg'])))
    return call


# ---------------------------------------------------------------- T11: three thermal signatures on one model

def t11_model():
    def build():
        import epcsaft
        from epcsaft import IdealCorrelation, IdealInterval, ThermochemistryRecord
        sys.path[:0] = [str(ROOT / 'engine'), str(ROOT / 'engine/tests')]
        from native_inputs import ETHANE, parameters
        lows = INPUTS['t11']['intervals']
        species = [[IdealInterval(lo, hi, True, i == len(lows) - 1, IdealCorrelation(nasa9(), 1e5)) for i, (lo, hi) in enumerate(lows)]
                   for _ in range(2)]
        return epcsaft.Mixture(parameters((ETHANE, ETHANE), [('k_ij', (0, 1), 0.0)]), thermochemistry=ThermochemistryRecord(1e5, species))
    return build()


def op_t11(req):
    """One reusable model per block; total enthalpy and its first T action at each call of the declared sequence."""
    from epcsaft import DensityActionDirection, PropertyObservable
    d = INPUTS['t11']
    sequence = d['sequences'][req['sequence']]
    obs, direction = [PropertyObservable.TotalEnthalpy], [DensityActionDirection(1., 0., [0., 0.])]

    def call():
        m = t11_model()
        for T in sequence[:len(set(sequence))]:  # one complete warm-up pass over the distinct states
            m.density_actions(T, d['rho'], x=d['x'], observables=obs, directions=direction)
        before = counters(m)
        ns, recordings, values = 0, [], {}
        for T in sequence:
            r0 = m.counters().recordings
            t = time.perf_counter_ns()
            b = m.density_actions(T, d['rho'], x=d['x'], observables=obs, directions=direction)
            ns += time.perf_counter_ns() - t
            recordings.append(m.counters().recordings - r0)
            values.setdefault(T, (b.results[0].value, b.results[0].action, b.results[0].status.name))
        after = counters(m)
        return dict(seconds=ns / 1e9, per_call=len(sequence), counters={k: after[k] - before[k] for k in after},
                    value=dict(recording_calls=sum(1 for r in recordings if r), values=values))
    return call


def op_t11_check(req):
    """Fresh-model values at each state, and the declared NASA9 ideal enthalpy h = 1000 + 30 (T - 300) J/mol."""
    d = INPUTS['t11']

    def call():
        rows = []
        for T in d['temperatures']:
            s = t11_model().state(T, rho=d['rho'], x=d['x'])
            rows.append(dict(T=T, total_enthalpy=s.total_enthalpy, ideal_enthalpy=s.ideal_enthalpy,
                             analytic_ideal=1000. + 30. * (T - 300.), status=s.thermochemistry_status.name))
        return dict(value=rows)
    return call


def op_t11_replay(req):
    """Retained T4 two-phase methane/ethane phase states, alternated on one reusable model."""
    states = req['states']

    def call():
        m = model('methane-ethane')
        before = counters(m)
        t = time.perf_counter_ns()
        for _ in range(req['passes']):
            for T, rho, x in states:
                m.state(T, rho=rho, x=x)
        ns = time.perf_counter_ns() - t
        after = counters(m)
        return dict(seconds=ns / 1e9, per_call=req['passes'] * len(states), counters={k: after[k] - before[k] for k in after})
    return call


# ---------------------------------------------------------------- T14: replay of the audit's benzene/toluene PH round trip

def op_t14_ph(req):
    """Fable p10b top-level statements from its first benzene declaration on (loops and try blocks dropped): the
    hydrogen/toluene cold and guessed flashes and part (3), unchanged declarations, packets from this checkout; untimed."""
    path = AUDIT / 'fable/probes/p10b_ph_ps_basis.py'
    tree = ast.parse(path.read_text())
    start = next(n.lineno for n in tree.body if isinstance(n, ast.Assign) and 'benzene' in ast.unparse(n.value))
    stop = next(n.lineno for n in tree.body if isinstance(n, ast.Expr) and 'write_csv' in ast.unparse(n))
    keep = [n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef)) and 'common' not in ast.unparse(n)
            or isinstance(n, ast.Assign) and (n.lineno < start and getattr(n.targets[0], 'id', '') in ('NASA9', 'R_SOURCE', 'rows')
                                               or start <= n.lineno < stop)
            or isinstance(n, ast.If) and start <= n.lineno < stop or isinstance(n, ast.Expr) and start <= n.lineno < stop]

    def call():
        ns = dict(MAIN=ROOT, packet=lambda name, version: ROOT / 'data/packets' / name / str(version) / 'parameters/parameters.json')
        exec(compile(ast.fix_missing_locations(ast.Module(body=keep, type_ignores=[])), str(path), 'exec'), ns)
        return dict(value=ns['rows'])
    return call


OPS = {name[3:].replace('_', '-'): f for name, f in globals().items() if name.startswith('op_')}


def declarations():
    """Construct every T1/T2/T7 model without solving (repository declaration check)."""
    return [model(r['model']) for r in INPUTS['t1']] + [model(c['model']) for c in INPUTS['t2']['ordinary']] + [model('ethane')]


def main():
    sys.path.insert(0, str(HERE / 'scripts'))
    import guard
    identity = guard.install(HERE / 'results/image.json', launcher=os.environ.get('CAMPAIGN_FIT_LAUNCHER'))
    replies = os.fdopen(os.dup(1), 'w')
    os.dup2(2, 1)  # native diagnostics stay off the reply pipe
    print(json.dumps(dict(identity=identity, threads={k: v for k, v in os.environ.items() if k.endswith('_THREADS')})),
          file=replies, flush=True)
    calls = {}
    for line in sys.stdin:
        req = json.loads(line)
        try:
            if req['op'] == 'recheck':
                reply = dict(maps=guard.recheck(HERE / 'results/image.json'))
                print(json.dumps(reply), file=replies, flush=True)
                continue
            key = json.dumps({k: v for k, v in req.items() if k not in ('samples', 'warmups', 'block')}, sort_keys=True)
            if key not in calls:
                t = time.perf_counter_ns()
                calls[key] = OPS[req['op']](req)
                prepare_ns = time.perf_counter_ns() - t
            else:
                prepare_ns = None
            f = calls[key]
            m = model(T1[req['state']]['model']) if req.get('state') in T1 else None
            for _ in range(req.get('warmups', 0)):
                f()
            samples = []
            for _ in range(req.get('samples', 1)):
                before = counters(m) if m else None
                t, c = time.perf_counter_ns(), time.process_time_ns()
                out = f()
                wall, cpu = time.perf_counter_ns() - t, time.process_time_ns() - c
                out.setdefault('seconds', wall / 1e9)
                out.setdefault('cpu_seconds', cpu / 1e9)
                if m:
                    after = counters(m)
                    out['counters'] = {k: after[k] - before[k] for k in after}
                samples.append(clean(out))
            reply = dict(samples=samples, prepare_seconds=None if prepare_ns is None else prepare_ns / 1e9)
        except Exception as error:
            traceback.print_exc()
            reply = dict(error=f'{type(error).__name__}: {error}')
        print(json.dumps(reply, allow_nan=False), file=replies, flush=True)


if __name__ == '__main__':
    main()

"""Release declarations for the existing #578 persistent worker and timers."""
import ast
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
RELEASE = ROOT / 'build/issue-681/release'
sys.path.insert(0, str(ROOT / 'analyses/2026-audit-timing-campaign/scripts'))
import guard

# Load definitions without importing the Engine; the existing worker installs its guard first.
import importlib.util
spec = importlib.util.spec_from_file_location('release_worker', ROOT / 'analyses/2026-audit-timing-campaign/scripts/engine_worker.py')
w = importlib.util.module_from_spec(spec)
spec.loader.exec_module(w)
w.HERE, w.ROOT = HERE, RELEASE
w.INPUTS = json.loads((HERE / 'inputs/inputs.json').read_text())
w.T1 = {r['id']: r for r in w.INPUTS['t1']}


def cross():
    return w.cached('release-cross', lambda: w.load(RELEASE / 'analyses/2026-competitor-cross-check/scripts/engine.py', 'release_cross'))


def model(key):
    def build():
        import epcsaft
        if key in ('propane', 'C1'):
            return epcsaft.Pure(epcsaft.Parameters.from_json(RELEASE / 'data/packets/gross-2001-propane/3/parameters/parameters.json'))
        if key == 'methane-ethane':
            return epcsaft.Mixture(epcsaft.Parameters.from_json(RELEASE / 'data/packets/gross-2001-methane-ethane/3/parameters/parameters.json', components=('methane', 'ethane')))
        raise KeyError(key)
    return w.cached(('model', key), build)


def equilibrium(req):
    cc = cross()
    inputs = json.loads((RELEASE / 'analyses/2026-competitor-cross-check/inputs/cases.json').read_text())
    request = next(r for r in inputs['requests'] if r['id'] == req['case'])
    prepared = cc.prepare(inputs, request)
    return lambda: dict(value=cc.call(request, prepared))


def electrolyte(req):
    """Reuse the head-to-head's constant-permittivity document, feeds and declared starts."""
    from epcsaft import Mixture, Parameters, equilibrium as q
    hh = w.load(RELEASE / 'analyses/2026-head-to-head-comparisons/scripts/h2h_engine.py', 'release_h2h')
    m = Mixture(Parameters.from_mapping(hh.khudaida_document('constant-permittivity'), components=hh.KH_IDS))
    problems = [(point, q.Problem(T=293.15, P=1e5, feed=q.Amounts(dict(zip(hh.KH_IDS, feed))), row_tolerance=1e-10, phases=[
        q.Phase('phase-0', kind='liquid', support=hh.KH_IDS, composition_guess=(.80, .05, .05, .05, .05)),
        q.Phase('phase-1', kind='liquid', support=hh.KH_IDS, composition_guess=(.20, .10, .6998, .0001, .0001))]))
        for point, feed in hh.khudaida_feeds().items()]
    assert len(problems) == 8

    def call():
        rows = []
        for point, p in problems:
            r = q.solve_equilibrium(m, p)
            if not r.success or not r.requested_tolerance_met:
                raise RuntimeError(f'feed {point}: {r.message}; residuals={r.residuals}')
            rows.append(dict(point=point, **w.equilibrium_record(r)))
        return dict(value=rows)
    return call


def sensitivity(ident):
    """Keep the existing head-to-head declarations and central-difference closures, omit its timers."""
    from epcsaft import equilibrium as q
    def build():
        hh = w.load(RELEASE / 'analyses/2026-head-to-head-comparisons/scripts/h2h_engine.py', 'release_sensitivity')
        tree = ast.parse(Path(hh.__file__).read_text())
        fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'sensitivity_case')
        fn.body = [n for n in fn.body if not isinstance(n, (ast.For, ast.Return)) and not (
            isinstance(n, ast.Assign) and any(isinstance(v, ast.Call) and isinstance(v.func, ast.Name) and v.func.id == 'timed' for v in ast.walk(n)))]
        fn.body += ast.parse('return dict(model=model, problem=prob, central=central, observables=observables, directions=directions, observe=observe, fd=finite_differences, steps=steps)').body
        exec(compile(ast.fix_missing_locations(ast.Module(body=[fn], type_ignores=[])), hh.__file__, 'exec'), hh.__dict__)
        inputs = json.loads((RELEASE / 'analyses/2026-competitor-cross-check/inputs/cases.json').read_text())
        c = hh.sensitivity_case(inputs, ident)
        if not c['central'].success or not c['central'].requested_tolerance_met:
            raise RuntimeError(c['central'].message)
        c['compiled'] = q.compile_problem(c['model'], c['problem'])
        return c
    return w.cached(('sensitivity', ident), build)


def action(req):
    from epcsaft import equilibrium as q
    c = sensitivity(req['case'])
    if req['mode'] == 'fd':
        return lambda: dict(value=c['fd']())
    def call():
        batches = q.solved_state_first_actions(c['compiled'], c['central'], c['observables'], [d for _, d in c['directions']])
        return dict(value={name: [r.action for r in b.results] for (name, _), b in zip(c['directions'], batches)},
                    diagnostics=[w.solved_record(b) for b in batches], central=c['observe'](c['central']), steps=c['steps'])
    return call


w.model = model
w.OPS.update({'release-equilibrium': equilibrium, 'release-electrolyte': electrolyte, 'release-action': action})


def declarations():
    """Construct the repository-owned declarations without the historical image or a solve."""
    from epcsaft import Mixture, Parameters, equilibrium as q
    cc = w.load(ROOT / 'analyses/2026-competitor-cross-check/scripts/engine.py', 'cost_declaration_cross')
    inputs = json.loads((ROOT / 'analyses/2026-competitor-cross-check/inputs/cases.json').read_text())
    declared = [cc.prepare(inputs, r, declare_only=True)[:2] for r in inputs['requests'] if r['id'] in ('C1-300', 'C2-0.5', 'C4-298.15')]
    m = Mixture(Parameters.from_json(ROOT / 'data/packets/gross-2001-methane-ethane/3/parameters/parameters.json', components=('methane', 'ethane')))
    declared.append((m, q.Problem(T=190., P=1e6, feed=q.MoleFractions({'methane': .5, 'ethane': .5}), phases=q.Discover())))
    hh = w.load(ROOT / 'analyses/2026-head-to-head-comparisons/scripts/h2h_engine.py', 'cost_declaration_h2h')
    m = Mixture(Parameters.from_mapping(hh.khudaida_document('constant-permittivity'), components=hh.KH_IDS))
    for feed in hh.khudaida_feeds().values():
        declared.append((m, q.Problem(T=293.15, P=1e5, feed=q.Amounts(dict(zip(hh.KH_IDS, feed))), row_tolerance=1e-10, phases=[
            q.Phase('phase-0', kind='liquid', support=hh.KH_IDS, composition_guess=(.80, .05, .05, .05, .05)),
            q.Phase('phase-1', kind='liquid', support=hh.KH_IDS, composition_guess=(.20, .10, .6998, .0001, .0001))])))
    return declared


if __name__ == '__main__':
    w.main()

"""Time commit-matched native declarations; retain every measured return."""
import ast
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time
import traceback
from types import ModuleType
from importlib.machinery import SourceFileLoader


def load(path, name):
    loader = SourceFileLoader(name, str(path))
    spec = importlib.util.spec_from_loader(name, loader)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    loader.exec_module(module)
    return module


def fixed_inputs(path):
    # Keep the original Engine declarations; omit the unrelated IDAES process model.
    tree = ast.parse(path.read_text())
    nodes = [n for n in tree.body if isinstance(n, (ast.Assign, ast.AnnAssign))
             or isinstance(n, ast.FunctionDef) and n.name in ('feed_rows', 'engine_problem')
             or isinstance(n, ast.Import) and all(a.name in ('csv', 'math', 'json', 'epcsaft') for a in n.names)
             or isinstance(n, ast.ImportFrom) and n.module in ('__future__', 'pathlib', 'epcsaft')]
    module = ModuleType('run_fixed_settler')
    module.__file__ = str(path)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), module.__dict__)
    sys.modules[module.__name__] = module
    return module


def state(result):
    if not result.success or result.requested_tolerance_met is not True:
        raise RuntimeError(f'{result.message}; residuals={result.residuals}')
    out = {k: getattr(result, k) for k in ('temperature', 'pressure', 'amounts', 'mole_fractions',
           'molar_densities', 'iterations', 'evaluations', 'hessian_calls')}
    # Retained wheels before #609 return one-entry lists; later ones return numbers.
    out.update({k: v if isinstance(v := out[k], list) else [v] for k in ('temperature', 'pressure')})
    out.update(raw_residual=max(map(abs, result.residuals)),
               classified_residual=max(map(abs, result.classified_residuals)),
               linear_solver=getattr(result, 'linear_solver', 'not exposed by this wheel'))
    if getattr(result, 'discovery', None) is not None:
        out['discovery_work'] = dict(result.discovery.work)
    for key in ('temperature', 'pressure', 'amounts', 'mole_fractions', 'molar_densities'):
        assert out[key] and all(map(math.isfinite, out[key])), key
    return out


def make(operation, source):
    import epcsaft
    from epcsaft import equilibrium as eq
    if operation == 'mea-f1/residual-jacobian':
        key = load(source / 'tools/time-key-solves', 'key_f1')
        if not hasattr(key, 'f1_result') and not hasattr(key, 'frozen_result'):
            # Use the same frozen F1 packet, but do not claim that old wheels expose later diagnostics.
            key = load(source.parent / 'cf26a17b/tools/time-key-solves', 'key_f1_historical')
        # Missing later work counters do not make a complete older numerical evaluation unavailable.
        tree = ast.parse(Path(key.__file__).read_text())
        name = 'f1_result' if hasattr(key, 'f1_result') else 'frozen_result'
        function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
        function.body = [n for n in function.body if not (isinstance(n, ast.If)
                         and ast.unparse(n.test) == 'not all(counts.values())')]
        exec(compile(ast.Module(body=[function], type_ignores=[]), str(key.__file__), 'exec'), key.__dict__)
        def f1():
            result = key.f1_result('evaluation') if name == 'f1_result' else key.frozen_result(key.F1, 'evaluation')
            result['missing_work_diagnostics'] = [name for name in key.F1_COUNTS
                                                   if not result['rows'][0].get(name)]
            return result
        return f1
    if operation.startswith(('case-', 'khudaida/', 'mea-nominally')):
        key = load(source / 'tools/time-key-solves', 'key')
        b = load(source / 'analyses/2026-tpd-phase-discovery/scripts/baseline.py', 'baseline')
        original_load = b.load
        b.load = lambda name, path: fixed_inputs(path) if name == 'run_fixed_settler' else original_load(name, path)
        if operation.startswith('mea-nominally'):
            calls = [(operation, *key.unloaded_declaration())]
        else:
            declarations = {name: (m, p) for name, m, p in b.declared_cases() if name in key.ORDINARY}
            if operation == 'khudaida/sequential-8':
                fixed = sys.modules['run_fixed_settler']
                m = epcsaft.Mixture(epcsaft.Parameters.from_json(b.KH_PACKET, components=fixed.COMPONENTS))
                calls = [(f'khudaida/point-{i}', m, fixed.engine_problem(fixed.feed_rows()[i])) for i in range(1, 9)]
            else:
                base, mode = operation.rsplit('/', 1)
                calls = [(operation, *declarations[base])]
        accepted = None
        if operation.endswith('/accepted'):
            _, m, p = calls[0]
            r = eq.solve_equilibrium(m, p)
            state(r)
            accepted = eq.EquilibriumStart()
            accepted.coordinates = r.coordinates

        def call():
            previous, rows, ns, cpu_ns = None, [], 0, 0
            for point, m, p in calls:
                t, c = time.perf_counter_ns(), time.process_time_ns()
                if accepted is None:
                    result = eq.solve_equilibrium(m, p, previous)
                else:
                    attempts = eq.solve_equilibrium_starts(m, p, [accepted], eq.EquilibriumWorkLimits())
                    result = attempts.results[0]
                ns += time.perf_counter_ns() - t
                cpu_ns += time.process_time_ns() - c
                if result is None:
                    raise RuntimeError(f'no accepted state: {attempts.causes}')
                rows.append(dict(point=point, **state(result)))
                previous = result
            return dict(seconds=ns / 1e9, cpu_seconds=cpu_ns / 1e9, rows=rows)
        return call

    if operation in extension_operations():
        model, problem = extension_declaration(operation, source)
    else:
        model, problem = None, None
    cc = load(source / 'analyses/2026-competitor-cross-check/scripts/engine.py', 'cross_check')
    inputs = json.loads((source / 'analyses/2026-competitor-cross-check/inputs/cases.json').read_text())
    if model is not None:
        pass
    elif operation in ('C1-300', 'C2-0.5', 'C4-298.15', 'C4-323.15'):
        q = next(q for q in inputs['requests'] if q['id'] == operation)
        prepared = cc.prepare(inputs, q)
        model, problem = prepared[:2]
    elif operation.startswith('methane-'):
        path = source / 'analyses/2026-binary-methane-ethane-tp-flash/scripts/performance.py'
        flash = load(path, 'methane_flash')
        model = epcsaft.Mixture(epcsaft.Parameters.from_json(flash.PACKET, components=flash.COMPONENTS))
        if operation == 'methane-flash':
            problem = flash.declare_inputs()
        else:
            # The bubble declaration in #533 uses these exact default guesses.
            problem = eq.bubble_point(T=180., liquid=eq.Phase('L', kind='liquid', composition_guess=(.5, .5)),
                                     feed=eq.MoleFractions({'methane': .5, 'ethane': .5}))
    elif operation == 'A-propane-liquid':
        q = next(q for q in inputs['requests'] if q['id'] == 'C1-300')
        model = cc.prepare(inputs, q)[0]
        declaration = json.loads((Path(__file__).resolve().parents[1] / 'results/density-input.json').read_text())

        def density():
            rows, ns, cpus = [], 0, 0
            for st in declaration['states']:
                t, c = time.perf_counter_ns(), time.process_time_ns()
                value = model.state(st['T_K'], P=st['P_Pa'], phase='liquid')
                ns += time.perf_counter_ns() - t
                cpus += time.process_time_ns() - c
                rho = value.molar_density
                assert math.isfinite(rho) and rho > 0
                rows.append(dict(point=str(st['index']), molar_densities=[rho], pressure=[value.pressure],
                                 density_oracle_relative_error=abs(rho / st['oracle_rho_mol_m3'] - 1)))
            # Cost per density call over the same changing-pressure sequence as #533.
            return dict(seconds=ns / len(rows) / 1e9, cpu_seconds=cpus / len(rows) / 1e9, rows=rows)
        return density
    else:
        raise ValueError(operation)

    def neutral():
        t, c = time.perf_counter_ns(), time.process_time_ns()
        result = eq.solve_equilibrium(model, problem)
        ns, cpus = time.perf_counter_ns() - t, time.process_time_ns() - c
        return dict(seconds=ns / 1e9, cpu_seconds=cpus / 1e9, rows=[dict(point=operation, **state(result))])
    return neutral


def extension_operations():
    return json.loads((Path(__file__).resolve().parents[1] / 'inputs/extension.json').read_text())


def extension_declaration(operation, source):
    import epcsaft
    from epcsaft import equilibrium as eq
    request = extension_operations()[operation]
    if request['kind'] == 'saturation':
        cc = load(source / 'analyses/2026-competitor-cross-check/scripts/engine.py', 'cold_propane')
        inputs = json.loads((source / 'analyses/2026-competitor-cross-check/inputs/cases.json').read_text())
        q = next(q for q in inputs['requests'] if q['id'] == 'C1-300') | dict(T=request['T_K'])
        return cc.prepare(inputs, q)[:2]
    packet = source / 'data/packets/gross-2001-methane-ethane/3/parameters/parameters.json'
    model = epcsaft.Mixture(epcsaft.Parameters.from_json(packet, components=['methane', 'ethane']))
    feed = eq.MoleFractions({'methane': .5, 'ethane': .5})
    if request['kind'] == 'discover':
        problem = eq.Problem(phases=eq.Discover(), T=request['T_K'], P=request['P_Pa'], feed=feed)
    elif request['kind'] == 'bubble':
        problem = eq.bubble_point(T=request['T_K'],
            liquid=eq.Phase('L', kind='liquid', composition_guess=(.5, .5)), feed=feed)
    else:
        problem = eq.dew_point(T=request['T_K'],
            vapor=eq.Phase('V', kind='vapor', composition_guess=(.5, .5)), feed=feed)
    return model, problem


def main():
    source = Path(sys.argv[1])
    replies = os.fdopen(os.dup(1), 'w')
    os.dup2(2, 1)
    import epcsaft
    from epcsaft import _core
    solver = _core.linear_solver() if hasattr(_core, 'linear_solver') else 'not exposed by this wheel'
    print(json.dumps(dict(pid=os.getpid(), affinity=sorted(os.sched_getaffinity(0)),
                          nice=os.getpriority(os.PRIO_PROCESS, 0), epcsaft=epcsaft.__file__,
                          extension=_core.__file__, linear_solver=solver)), file=replies, flush=True)
    calls = {}
    for line in sys.stdin:
        req = json.loads(line)
        try:
            name = req['operation']
            if name not in calls:
                calls[name] = make(name, source)
                calls[name]()  # one complete, untimed warmup per operation and process
            values = [calls[name]() for _ in range(req['samples'])]
            for value in values:
                value['linear_solver'] = solver
                for row in value['rows']:
                    row.setdefault('linear_solver', solver)
            response = dict(values=values)
        except Exception as error:
            traceback.print_exc()
            response = dict(unavailable=f'{type(error).__name__}: {error}')
        print(json.dumps(response, allow_nan=False), file=replies, flush=True)


def declarations():
    source = Path(__file__).resolve().parents[3]
    cc = load(source / 'analyses/2026-competitor-cross-check/scripts/engine.py', 'history_declarations_cc')
    inputs = json.loads((source / 'analyses/2026-competitor-cross-check/inputs/cases.json').read_text())
    declared = [cc.prepare(inputs, q)[:2] for q in inputs['requests']
                if q['id'] in ('C1-300', 'C2-0.5', 'C4-298.15', 'C4-323.15')]
    b = load(source / 'analyses/2026-tpd-phase-discovery/scripts/baseline.py', 'history_declarations_b')
    original_load = b.load
    b.load = lambda name, path: fixed_inputs(path) if name == 'run_fixed_settler' else original_load(name, path)
    declared.extend((m, p) for name, m, p in b.declared_cases())
    key = load(source / 'tools/time-key-solves', 'history_declarations_key')
    declared.append(key.unloaded_declaration())
    flash = load(source / 'analyses/2026-binary-methane-ethane-tp-flash/scripts/performance.py', 'history_declarations_flash')
    from epcsaft import Mixture, Parameters
    declared.append((Mixture(Parameters.from_json(flash.PACKET, components=flash.COMPONENTS)), flash.declare_inputs()))
    declared.extend(extension_declaration(op, source) for op in extension_operations())
    return declared


if __name__ == '__main__':
    main()

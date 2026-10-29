"""Guarded release T4 worker, with native request mapping outside prepared timing."""
import importlib.util
import os
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location('release_engine', ROOT/'analyses/2026-release-computational-cost/scripts/engine_worker.py')
e = importlib.util.module_from_spec(spec); spec.loader.exec_module(e)
e.RELEASE = Path(os.environ.get('EPCSAFT_COST_SOURCE', str(ROOT)))
e.w.ROOT = e.RELEASE
if 'EPCSAFT_COST_IMAGE' in os.environ: e.w.HERE = Path(os.environ['EPCSAFT_COST_IMAGE']).parents[1]


def declarations():
    from epcsaft import equilibrium as q
    model = e.w.t4_model('methane_ethane')
    T, P = e.w.t4_points('methane_ethane')[0]
    problem = q.Problem(T=T, P=P, feed=q.MoleFractions(dict.fromkeys(model.component_ids, .5)), phases=q.Discover())
    return [(model, problem)]


def prepared(req):
    from epcsaft import equilibrium as q
    model, problem = declarations()[0]
    T, P = problem.T, problem.P
    native = problem._native(model)
    def call():
        start = time.perf_counter_ns()
        r = q.solve_equilibrium(model, native)
        seconds = (time.perf_counter_ns()-start)/1e9
        return dict(seconds=seconds, per_call=1, value=[dict(T=T, P=P, **e.w.equilibrium_record(r))])
    return call


e.w.OPS['prepared-discovery'] = prepared
if __name__ == '__main__':
    os.environ.pop('CAMPAIGN_FIT_LAUNCHER', None)
    e.w.main()

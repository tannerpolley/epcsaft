"""FeOs 0.10.1 side of the 2026-10-08 head-to-head: C2 parameter gradients (native AD) and advanced-variant Born terms."""
import importlib.util, json, statistics, sys, time
from pathlib import Path
CC = Path(__file__).resolve().parents[2] / '2026-competitor-cross-check'
spec = importlib.util.spec_from_file_location('feos_case', CC / 'scripts/feos.py'); fc = importlib.util.module_from_spec(spec); spec.loader.exec_module(fc)
import numpy as np
from feos import BinaryRecord, Contributions, EquationOfState, Identifier, Parameters, PureRecord, State
from si_units import KELVIN, METER, MOL, RGAS
RHO = MOL / METER**3


def timed(f, warm=3, samples=10):
    t = time.perf_counter_ns(); value = f(); first = time.perf_counter_ns() - t
    for _ in range(warm): f()
    ns = []
    for _ in range(samples):
        t = time.perf_counter_ns(); f(); ns.append(time.perf_counter_ns() - t)
    return value, dict(first_ms=first / 1e6, median_ms=statistics.median(ns) / 1e6, samples=samples)


def born_eos(comps, form):
    records = []
    for c in comps:
        kw = dict(m=c['m'], sigma=c['sigma'], epsilon_k=c['epsilon'], z=c['z'])
        if c['id'] == 'water': kw['association_sites'] = [dict(kappa_ab=0.04509, epsilon_k_ab=2425.7, na=1., nb=1.)]
        kw['permittivity_record'] = ({'ExperimentalData': {'data': [[298., c['relative_permittivity']], [299., c['relative_permittivity']]]}} if form == 'mole-fraction-mixing' else
                                     {'PerturbationTheory': {k: c[k] for k in ('dipole_scaling', 'polarizability_scaling', 'correlation_integral_parameter')}})
        records.append(PureRecord(Identifier(name=c['id']), 18.0153 if c['id'] == 'water' else 23.0, **kw))
    ids = [Identifier(name=c['id']) for c in comps]
    binaries = [BinaryRecord(ids[0], ids[1], k_ij=-0.3), BinaryRecord(ids[0], ids[2], k_ij=-0.3), BinaryRecord(ids[1], ids[2], k_ij=0.8)]
    return EquationOfState.epcsaft(Parameters.from_records(records, binaries), epcsaft_variant='advanced')


def main():
    engine = json.loads(Path(sys.argv[1]).read_text()); inputs = json.loads((CC / 'inputs/cases.json').read_text()); out = {}
    req = next(x for x in inputs['requests'] if x['id'] == 'C2-gradient')
    prepared = fc.prepare(inputs, req)
    t = time.time(); value, cost = timed(lambda: fc.gradient(prepared[0], prepared[2], req)); out['C2_gradient'] = dict(value=value, cost=cost, wall_s=time.time() - t)
    rows = []; t = time.time()
    for form in ('mole-fraction-mixing', 'perturbation'):
        eos = born_eos(engine['born_comparators']['components'], form)
        for st in engine['born_states']:
            s = State(eos, temperature=st['T'] * KELVIN, density=st['rho'] * RHO, composition=np.array(st['x']))
            contrib = {n: v / (RGAS * s.temperature) for n, v in s.residual_molar_helmholtz_energy_contributions()}
            rows.append(dict(form=form, m=st['m'], a_born=float(contrib.get('Born')), a_ionic=float(contrib.get('Ionic')),
                             a_residual=float(s.molar_helmholtz_energy(Contributions.Residual) / (RGAS * s.temperature))))
    out['born'] = rows; out['born_wall_s'] = time.time() - t
    Path(sys.argv[2]).write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()

"""#243 matrix: the six declared baseline.py cases at their original feed, then cold and from the reference
state at 1, 100 and 10000 mol, compared with --reference original-states.json under the accepted quantity
tolerances (42 rows). The reference start declares the reference phase amounts times the feed ratio, its phase
compositions, packing fractions and free T/P; the Engine builds the coordinates. Native Engine solves only.
usage: physical_comparison.py OUTPUT_DIR [--reference FILE]"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import time

from epcsaft import _core, equilibrium as eq

import baseline

cases = list(baseline.declared_cases())
parser = argparse.ArgumentParser()
parser.add_argument('output', type=Path)
parser.add_argument('--reference', type=Path)
args = parser.parse_args()
folder = args.output
folder.mkdir(parents=True, exist_ok=False)
rows = []
fields = ('success', 'message', 'solver_status', 'phase_count', 'iterations', 'evaluations',
          'hessian_calls', 'objective', 'coordinates', 'amounts', 'mole_fractions', 'molar_densities',
          'packing_fractions', 'multipliers', 'constraint_multipliers', 'residuals',
          'classified_residuals', 'temperature', 'pressure', 'requested_tolerance_met')


def run(name, lane, model, problem, start=None, previous=None):
    begun = time.perf_counter()
    with _core.Capture() as capture:
        if start is None:
            result = eq.solve_equilibrium(model, problem, previous)
        else:
            attempts = eq.solve_equilibrium_starts(model, problem, [start], eq.EquilibriumWorkLimits())
            result = attempts.results[0]
            if result is None:
                raise RuntimeError(str(attempts.causes))
    row = {key: getattr(result, key) for key in fields} | {key: [] if (v := getattr(result, key)) is None else [v] for key in ('temperature', 'pressure')}
    N = (sum(problem.feed.values.values()) if isinstance(problem.phases, eq.Discover)
         else eq.compile_problem(model, problem).feed_total)
    row.update(case=name, lane=lane, N=N, topology=str(result.topology_event), seconds=time.perf_counter() - begun,
               inventory=max((abs(value) for value, item in zip(result.residuals, result.rows)
                              if item.kind == _core.RowKind.Moiety), default=0.))
    rows.append(row)
    trace = [dict(name=e.name, span=e.span, iteration=e.iteration, values=e.values) for e in capture.events]
    trace_name = f'{len(rows):02d}-{name.replace("/", "-")}-{lane}-N{N:g}.json'
    (folder / trace_name).write_text(json.dumps(trace, indent=2) + '\n')
    print(name, lane, N, result.solver_status, result.success,
          result.iterations, result.evaluations, result.hessian_calls, result.objective, flush=True)
    return result, row


references = json.loads(args.reference.read_text()) if args.reference else {}
frozen = {}
for name, model, problem in cases:
    result, row = run(name, 'reference' if not args.reference else 'original', model, problem)
    frozen[name] = row
    if not result.success:
        raise RuntimeError(f'{name}: original feed failed: {result.message}')
    if args.reference:
        source = references[name]
        ids, x = list(model.component_ids), source['mole_fractions']
        for target in (1., 100., 10000.):
            trial = eq.Problem(T=problem.T, P=problem.P, feed=eq.Amounts({
                key: value * (target / source['N']) for key, value in problem.feed.values.items()}),
                phases=problem.phases, reactions=problem.reactions, row_tolerance=problem.row_tolerance)
            for lane in ('cold', 'accepted'):
                start = None if lane == 'cold' else eq.EquilibriumStart(
                    amounts=[amount * (target / source['N']) for amount in source['amounts']],
                    compositions=[[x[p * len(ids) + ids.index(s)] for s in phase.support or ids]
                                  for p, phase in enumerate(problem.phases)],
                    packing_fractions=source['packing_fractions'],
                    temperature=source['temperature'][0] if isinstance(problem.T, eq.Free) else None,
                    pressure=source['pressure'][0] if isinstance(problem.P, eq.Free) else None)
                run(name, lane, model, trial, start)

# The specific diagnosis uses an accepted state at each feed size: discovery at
# 100 mol, the cold split at 1 and 10000 mol. Keep this separate from the 42 rows.
name, model, original = cases[0]
for N in (1., 100., 10000.):
    problem = eq.Problem(T=original.T, P=original.P, feed=eq.Amounts({'water': .6*N, 'isobutanol': .4*N}),
                         phases=original.phases, row_tolerance=original.row_tolerance)
    cold, _ = run('neutral-diagnosis', 'cold', model, problem)
    if N == 100.:
        discovery = eq.Problem(T=problem.T, P=problem.P, feed=problem.feed, phases=eq.Discover(), row_tolerance=1e-10)
        accepted, _ = run('neutral-diagnosis', 'discovery', model, discovery)
        problem = accepted.discovery.declaration
    else:
        accepted = cold
    if accepted.success:
        run('neutral-diagnosis', 'accepted', model, problem, previous=accepted)

(folder / 'rows.json').write_text(json.dumps(rows, indent=2) + '\n')
(folder / 'original-states.json').write_text(json.dumps(frozen, indent=2) + '\n')
(folder / 'environment.json').write_text(json.dumps(dict(extension=_core.__file__,
    nice=os.getpriority(os.PRIO_PROCESS, 0), threads={k: os.getenv(k) for k in
    ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')}), indent=2) + '\n')
if args.reference:
    comparisons = []
    for row in rows:
        if row['case'] == 'neutral-diagnosis':
            continue
        ref = references[row['case']]
        def difference(key, relative=False, normalized=False):
            a, b = row[key], ref[key]
            if len(a) != len(b) or not a:
                return None
            if normalized:
                return max(abs(x/row['N'] - y/ref['N']) for x, y in zip(a, b))
            if relative:
                return max(abs(x/y - 1.) for x, y in zip(a, b))
            return max(abs(x-y) for x, y in zip(a, b))
        # Amounts and the Gibbs objective are extensive: compared per mol of declared feed.
        f, f_ref = row['objective'] / row['N'], ref['objective'] / ref['N']
        errors = dict(dx=difference('mole_fractions'), dq=difference('amounts', normalized=True),
                      drho=difference('molar_densities', relative=True), dT=difference('temperature', relative=True),
                      dP=difference('pressure', relative=True), df=abs(f - f_ref) / max(1., abs(f_ref)),
                      inventory=row.get('inventory'), classified=max(map(abs, row['classified_residuals']), default=math.inf))
        passes = row['success'] and row['requested_tolerance_met'] and row['topology'] == ref['topology'] and all(
            errors[k] is not None and math.isfinite(errors[k]) and errors[k] <= limit for k, limit in
            dict(dx=1e-9, dq=1e-12, drho=1e-9, dT=1e-9, dP=1e-9, df=1e-10, inventory=1e-12, classified=1e-10).items())
        comparisons.append(dict(case=row['case'], lane=row['lane'], N=row['N'], passes=passes,
                                I=row['iterations'], E=row['evaluations'], H=row['hessian_calls'],
                                objective=row['objective'], solver_status=row['solver_status'], topology=row['topology'],
                                requested_tolerance_met=row['requested_tolerance_met'], message=row['message'], **errors))
    with (folder / 'comparisons.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=comparisons[0])
        writer.writeheader()
        writer.writerows(comparisons)
    summary = dict(total=len(comparisons), passes=sum(r['passes'] for r in comparisons),
                   maxima={k: max(r[k] for r in comparisons if r[k] is not None) for k in errors},
                   failures=[r for r in comparisons if not r['passes']])
    (folder / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))

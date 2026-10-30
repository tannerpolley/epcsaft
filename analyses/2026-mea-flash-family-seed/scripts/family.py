"""The 73-row reacting MEA flash family through public solve_equilibrium on the running Engine.

Usage: python family.py <label> <out.jsonl>. Python declares inputs and reads results only; points run in parallel
processes, one native solve each. Declaration: case-7 flash of analyses/2026-tpd-phase-discovery/scripts/baseline.py
with loading, T and P varied (feed MEA 9.7292821, water 76.9751925, CO2 = loading x MEA mol; row tolerance 1e-10).
Same declaration and record fields as the #350/#356 driver (build/issue-356-build/family.py in the #356 archive).
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'tools'))
from engine_identity import loaded_identity

sys.path.insert(0, str(ROOT / 'analyses/2026-tpd-phase-discovery/scripts'))
import baseline as b  # noqa: E402
import epcsaft  # noqa: E402
from epcsaft import equilibrium as eq  # noqa: E402

MEA, WATER = 9.7292821, 76.9751925
POINTS = [tuple(float(v) for v in line.split())
          for line in (Path(__file__).parent / 'family-points.txt').read_text().splitlines() if line.strip()]


def declare(index):
    loading, T, P = POINTS[index]
    model, template = next((m, p) for name, m, p in b.declared_cases() if name == 'case-7/flash-5918.5286-Pa')
    return model, eq.Problem(T=T, P=P, phases=template.phases, reactions=template.reactions, row_tolerance=1e-10,
                             feed=eq.Amounts({'carbon-dioxide': loading * MEA, 'monoethanolamine': MEA, 'water': WATER}))


def solve(index):
    loading, T, P = POINTS[index]
    model, problem = declare(index)
    row = dict(index=index, loading=loading, T=T, P=P)
    start = time.perf_counter()
    try:
        r = eq.solve_equilibrium(model, problem)
        row.update(success=r.success, met=r.requested_tolerance_met, status=r.solver_status, message=r.message,
                   iterations=r.iterations, evaluations=r.evaluations, hessian_calls=r.hessian_calls,
                   topology=str(r.topology_event), amounts=list(r.amounts), mole_fractions=list(r.mole_fractions),
                   molar_densities=list(r.molar_densities), residual=r.convergence.residual)
    except Exception as error:  # a raised solve is a retained outcome, not a skipped row
        row.update(success=False, message=f'raised {type(error).__name__}: {error}')
    row['seconds'] = time.perf_counter() - start
    return row


def main():
    label, out = sys.argv[1], Path(sys.argv[2])
    identity = loaded_identity()
    provenance = {key: identity[key] for key in ('wheel', 'wheel_sha256', 'source_revision', 'source_archive_sha256', 'extension_sha256', 'native_library_sha256')}
    with ProcessPoolExecutor(os.cpu_count()) as pool:
        rows = list(pool.map(solve, range(len(POINTS))))
    out.write_text(''.join(json.dumps(dict(label=label, engine_provenance=provenance, wheel_sha256=identity['wheel_sha256'], **row)) + '\n'
                           for row in rows))


def declarations():
    return [declare(i) for i in range(len(POINTS))]


if __name__ == '__main__':
    main()

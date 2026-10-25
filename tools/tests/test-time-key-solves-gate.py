"""Planted work, state and missing-state changes against the time-key-solves gate (#442)."""
import importlib.machinery
import importlib.util
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

loader = importlib.machinery.SourceFileLoader('time_key_solves', str(Path(__file__).resolve().parents[1] / 'time-key-solves'))
tool = importlib.util.module_from_spec(importlib.util.spec_from_loader(loader.name, loader))
loader.exec_module(tool)
CASE = 'case-6/figure-2-point-7/cold'
X = [0.9, 0.1, 1e-20]  # mole fractions with a trace entry


def rows(counts_b=(10, 15, 11), fractions_b=X, case=CASE, empty_b=(), skip_b=()):
    """Both wheels solve every case once; wheel B differs only in `case` as planted."""
    work, states = [], []
    for name in tool.CASES:
        for wheel in 'AB':
            if wheel == 'B' and name in skip_b:
                continue
            planted = wheel == 'B' and name == case
            c = counts_b if planted else (10, 15, 11)
            work.append(dict(wheel=wheel, case=name, sample='warmup', point=name,
                             **dict(zip(tool.count_keys(name), (*c, 12) if name in tool.FROZEN else c, strict=True))))
            for quantity in tool.quantities(name):
                values = fractions_b if planted and quantity == 'mole_fractions' else X[:1] if name == tool.CP else X
                if wheel == 'B' and (name, quantity) in empty_b:
                    values = []
                states.extend(dict(wheel=wheel, case=name, sample='warmup', point=name, quantity=quantity,
                                   index=i, value=v) for i, v in enumerate(values))
    return work, states


above = [0.9 * (1 + 2e-10), 0.1, 1e-20]
for label, (work, states), expect_change, fails in (
        ('identical', rows(), False, False),
        ('one extra iteration', rows(counts_b=(11, 15, 11)), False, True),
        ('extra iteration, --expect-change', rows(counts_b=(11, 15, 11)), True, False),
        ('extra iteration outside case-6/case-7', rows(counts_b=(11, 15, 11), case='khudaida/sequential-8'), False, False),
        ('state 2e-10 above', rows(fractions_b=above), False, True),
        ('state 5e-11 below', rows(fractions_b=[0.9 * (1 + 5e-11), 0.1, 1e-20]), False, False),
        ('trace entry doubled, below its floor', rows(fractions_b=[0.9, 0.1, 2e-20]), False, False),
        ('state above, --expect-change', rows(fractions_b=above), True, True),
        ('no rows at all', ([], []), False, True),
        ('no state rows', (rows()[0], []), False, True),
        ('wheel B quantity empty', rows(empty_b=[(CASE, 'pressure')]), False, True),
        ('wheel B case not solved', rows(skip_b=[CASE]), False, True)):
    assert bool(tool.gate(work, states, expect_change)) is fails, label
print('time-key-solves gate plants: ok')

# The F1 performance claim is no increased work, rather than exact equality.
for key in tool.F1_COUNTS:
    work, states = rows()
    planted = next(row for row in work if row['wheel'] == 'B' and row['case'] == tool.F1)
    planted[key] += 1
    assert tool.gate(work, states), key
    assert not tool.gate(work, states, expect_change=True), key
    planted[key] -= 2
    assert not tool.gate(work, states), f'{key} decrease'

frozen = tool.FROZEN[tool.F1]['identity']['frozen']
ceiling = 1.25 * frozen['serial_cpu_seconds']
evaluation = dict(cpu_seconds=ceiling, seconds=ceiling, reference_taylor_evaluations=664, work=frozen['work'])
fit = dict(seconds=120, critical_cpu_seconds=10, status='converged')
assert not tool.ceilings(tool.F1, {'A': evaluation, 'B': evaluation}, [('A', fit), ('B', fit)])
assert tool.ceilings(tool.F1, {'B': dict(evaluation, cpu_seconds=ceiling * 1.001)})
invalid_control = subprocess.run([sys.executable, '-c',
    f"import runpy; tool = runpy.run_path({loader.path!r}); "
    f"tool['ceilings']({tool.F1!r}, {{'A': {dict(evaluation, cpu_seconds=ceiling * 1.007)!r}}})"],
    capture_output=True, text=True)
assert invalid_control.returncode == 1 and 'invalid host measurement: rerun on a quieter host' in invalid_control.stderr
assert tool.ceilings(tool.F1, {}, [('B', dict(fit, seconds=120.001))])
# Only an arm doing the frozen work measures the host; a candidate B is held to the ceiling whatever its work.
other = dict(evaluation, cpu_seconds=10 * ceiling, work=dict(frozen['work'], reference_columns=1))
assert not tool.ceilings(tool.F1, {'A': other, 'B': evaluation})
assert tool.ceilings(tool.F1, {'A': evaluation, 'B': other})
assert not tool.ceilings(tool.CP, {'A': other, 'B': other}, control=True)
assert tool.gate(*rows(skip_b=[tool.F1]), expect_change=True)
print('F1 A=A, planted work increases, fixed ceilings and missing evidence: ok')

# A failed attempt followed by a successful attempt can share one trace span.
trace = SimpleNamespace(
    spans=[SimpleNamespace(name='solve_equilibrium', convergence=SimpleNamespace(iterations=3))],
    events=[SimpleNamespace(name='ipopt_iteration', span=0, iteration=i) for i in (0, 1, 2, 2, 4, 0, 1, 3)])
assert tool.equilibrium_iterations(trace) == 7

work, states = rows()
work.extend(dict(row, sample='timer', reference_computations=row['reference_computations'] + (row['wheel'] == 'B'))
            for row in list(work) if row['case'] == tool.F1)
states.extend(dict(row, sample='timer') for row in list(states) if row['case'] == tool.F1)
assert tool.gate(work, states, expect_change=True)

controls = [({'A': evaluation}, {'A': fit, 'B': fit}),
            ({'A': evaluation, 'B': evaluation}, {'A': fit, 'B': dict(fit, critical_cpu_seconds=float('nan'))}),
            ({'A': evaluation, 'B': evaluation}, {'A': fit, 'B': dict(fit, status='failed')})]
assert all(tool.ceilings(tool.F1, evaluations, fits.items(), control=True) for evaluations, fits in controls)

work, states = rows()
next(row for row in states if row['case'] == tool.F1 and row['wheel'] == 'B')['value'] = float('nan')
assert tool.gate(work, states, expect_change=True)
print('F1 attempt totals, timer invariance, invalid controls and nonfinite retained states: ok')

"""Planted work, state and missing-state changes against the time-key-solves gate (#442)."""
import importlib.machinery
import importlib.util
from pathlib import Path

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
                             iterations=c[0], evaluations=c[1], hessian_calls=c[2]))
            for quantity in tool.PHYSICAL:
                values = fractions_b if planted and quantity == 'mole_fractions' else X
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

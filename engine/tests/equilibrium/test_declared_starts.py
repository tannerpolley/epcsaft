"""Ordered declared starts: analytic exchange amounts and a shared work quota (#238)."""
from typing import Any

from epcsaft import Mixture
from epcsaft import equilibrium as eq
from ..native_inputs import ETHANE, parameters


def test_declared_starts_close_analytic_exchange_and_bound_work(evidence: Any) -> None:
    p = parameters((ETHANE,) * 6, pairs=[("k_ij", (i, j), 0.) for i in range(6) for j in range(i+1, 6)])
    model = Mixture(p)
    problem = eq.Problem(T=200., P=5e6,
        phases=[eq.Phase("aqueous", kind="liquid", support=("s0", "s2", "s4")),
                eq.Phase("organic", kind="liquid", support=("s1", "s3", "s5"))],
        feed=eq.Amounts({"s0": 1., "s1": 1., "s2": .3, "s4": .3}), row_tolerance=1e-10,
        reactions=[eq.Reaction({"s2": -1., "s3": 1.}, equilibrium_constant=2., name="exchange_a"),
                   eq.Reaction({"s4": -1., "s5": 1.}, equilibrium_constant=.5, name="exchange_b")])
    exact = ([10/13, 1/13, 2/13], [10/13, 2/13, 1/13])
    packing = []
    for support, x in zip(((0, 2, 4), (1, 3, 5)), exact, strict=True):
        full = [0.] * 6
        for i, value in zip(support, x, strict=True):
            full[i] = value
        packing.append(model._model.packing_fraction_coefficient(200., full)
                       * model.state(200., P=5e6, x=full, phase="liquid").molar_density)

    def seed(x: Any, eta: Any) -> eq.EquilibriumStart:
        out = eq.EquilibriumStart()
        out.compositions = x
        out.packing_fractions = eta
        out.amounts = [1.3, 1.3]
        return out

    invalid = seed(exact, [-1., packing[1]])
    positive = seed(exact, packing)
    later = seed(([.8, .08, .12], [.8, .12, .08]), packing)
    solved = eq.solve_equilibrium_starts(model, problem, [invalid, positive, later], eq.EquilibriumWorkLimits())
    assert solved.selected == 1, solved.causes
    assert solved.results[0] is None and solved.results[2] is None
    result = solved.results[1]
    assert result is not None
    for phase, support in enumerate(((0, 2, 4), (1, 3, 5))):
        error = max(abs(result.amounts[phase] - 1.3),
                    *(abs(result.mole_fractions[6*phase+i] - x) for i, x in zip(support, exact[phase], strict=True)))
        assert error <= 1e-9
        evidence.compare(f"declared-starts/exchange/{phase}", model=error, reference=0., unit="mol or mole fraction",
                         criterion={"rule": "symmetric_close", "atol": 1e-9})
    raw = max(abs(r) for r in result.residuals)
    assert raw <= 1e-10
    for a, b, target in ((2, 3, 2.), (4, 5, .5)):
        ratio = result.amounts[1]*result.mole_fractions[6+b] / (result.amounts[0]*result.mole_fractions[a])
        assert abs(ratio-target) <= 1e-8
    only = eq.solve_equilibrium_starts(model, problem, [invalid], eq.EquilibriumWorkLimits())
    assert only.selected is None and only.iterations == 0
    for name in ("iterations", "seconds"):
        limits = eq.EquilibriumWorkLimits()
        setattr(limits, name, 0)
        stopped = eq.solve_equilibrium_starts(model, problem, [positive], limits)
        assert stopped.selected is None and stopped.budget_exhausted and stopped.results[0] is None
    limits = eq.EquilibriumWorkLimits()
    limits.iterations = 1
    stopped = eq.solve_equilibrium_starts(model, problem, [invalid, later, positive], limits)
    assert stopped.selected is None and stopped.budget_exhausted and stopped.iterations <= 1
    assert stopped.results[1] is not None and stopped.results[2] is None
    unfinished = max(abs(r) for r in stopped.results[1].residuals)
    assert unfinished > 1e-10
    evidence.observe("declared-starts/work", raw_rows=raw, limited_raw_rows=unfinished,
                     limited_iterations=stopped.iterations, selected=solved.selected)

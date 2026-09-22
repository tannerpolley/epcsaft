"""#87/#90 trace-state reliability: the analytic trace solution below the former composition box
and the intensive accuracy of trace rows after the amount-weighted optimizer stop."""
import math
from typing import Any

import pytest

from epcsaft import Mixture
from epcsaft.equilibrium import (
    Amounts, EOS_STANDARD_STATE_ID, Phase, Problem, Reaction, ReactionLogPolynomial,
    SolvedStateActionDirection, SolvedStateActionRequest, SolvedStateActionStatus,
    SolvedStateObservable, SolvedStateObservableKind as Kind,
    compile_problem, solve_equilibrium, solved_state_actions,
)
from tests.native_inputs import ETHANE, parameters, same


@pytest.mark.parametrize("total", [1.0, 10.0])
@pytest.mark.parametrize("equilibrium_constant", [1e-4, 1e-8, 1e-12])
def test_trace_isomerisation_meets_intensive_screens(evidence: Any, equilibrium_constant: float, total: float) -> None:
    """Identical-ethane A <=> B at 350 K and 200000 Pa has x_B = K/(1+K) exactly, with the logistic
    A1-A3 references of analyses/2026-trace-action-diagnosis. At K = 1e-8 the solution lies below the
    former 1e-8 log-composition box of the driver, which rejected both feeds (classified residuals
    2.98e-10 and 1.09e-10 against the 1e-10 request despite Ipopt SUCCESS). Screens are the retained
    diagnosis screens: composition 1e-6 relative, actions 5e-6, 5e-5 and 5e-4 relative, no absolute floor."""
    model = Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), 0.0)]))
    temperature = 350.0
    l1 = -1200 / temperature**2 + 0.4 / temperature + 0.002
    l2 = 2400 / temperature**3 - 0.4 / temperature**2
    l3 = -7200 / temperature**4 + 0.8 / temperature**3
    p = equilibrium_constant / (1 + equilibrium_constant)
    references = [p, p * (1 - p) * l1, p * (1 - p) * (l2 + (1 - 2 * p) * l1**2),
                  p * (1 - p) * (l3 + 3 * (1 - 2 * p) * l1 * l2 + (1 - 6 * p + 6 * p * p) * l1**3)]
    a = math.log(equilibrium_constant) - 1200 / temperature - 0.4 * math.log(temperature / 300) - 0.002 * temperature
    law = ReactionLogPolynomial(a, 1200, 0.4, 0.002, 300, 250, 450, EOS_STANDARD_STATE_ID)
    problem = Problem(phases=[Phase("V", kind="vapor")], T=temperature, P=2e5,
                      feed=Amounts({"s0": total / 2, "s1": total / 2}), row_tolerance=1e-10,
                      reactions=[Reaction({"s0": -1.0, "s1": 1.0}, correlation=law)])
    central = solve_equilibrium(model, problem)
    assert central.success and central.requested_tolerance_met and central.convergence.residual <= 1e-10, central.message
    label = f"trace/{equilibrium_constant:g}/{total:g}"
    inputs = {"K": equilibrium_constant, "feed_mol": total, "T_K": temperature, "P_Pa": 2e5,
              "classified_residual": central.convergence.residual,
              "unscaled_residual_max": max(abs(value) for value in central.residuals),
              "solver_status": central.solver_status, "iterations": central.iterations}
    same(evidence, f"{label}/x_B", central.mole_fractions[1], references[0], atol=0.0, rtol=1e-6, unit="1",
         inputs=inputs, basis="Identity", method="logistic-analytic",
         claim="trace mole fraction equals K/(1+K) without a composition box")
    if equilibrium_constant <= 1e-8:
        assert central.mole_fractions[1] < 1e-8  # the analytic solution sits below the former box
    compiled = compile_problem(model, problem)
    observable = SolvedStateObservable(Kind.PhaseMoleFraction, 0, 1)
    direction = SolvedStateActionDirection(1.0, 0.0, [0.0, 0.0], [])
    for order, tolerance in enumerate((5e-6, 5e-5, 5e-4), 1):
        item = solved_state_actions(compiled, central, SolvedStateActionRequest([observable], [direction] * order)).results[0]
        assert item.status == SolvedStateActionStatus.Available, item.status
        same(evidence, f"{label}/A{order}", item.action, references[order], atol=0.0, rtol=tolerance, unit=f"K^-{order}",
             inputs=inputs, basis="Identity", method="logistic-analytic",
             claim=f"order-{order} temperature action of the trace mole fraction matches the logistic chain")
    same(evidence, f"{label}/material_dual", central.constraint_multipliers[0], -total * central.multipliers[0],
         atol=0.0, rtol=1e-8, unit="mol", inputs=inputs, basis="Identity", method="normalized-constraint",
         claim="the extensive objective multiplier is -N times the normalized balance multiplier")

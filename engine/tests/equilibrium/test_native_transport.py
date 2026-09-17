"""B1-B3 (#56-#59): public reaction targets, solved-state actions and result diagnostics."""
import math
from typing import Any

import pytest

from epcsaft import Mixture, ParameterError
from epcsaft.equilibrium import (
    EOS_STANDARD_STATE_ID,
    Amounts,
    Phase,
    Problem,
    Reaction,
    ReactionLogPolynomial,
    solve_equilibrium,
)
from tests.native_inputs import parameters, same

NU = (-1.0, -1.0, 1.0)
B, T0 = 1500.0, 298.15  # manufactured EOS-basis ln K = b (1/T - 1/T0), K(T0) = 1


def reactive_model() -> Mixture:
    """equilibrium_support::synthetic_reactive_ternary(): A + B = C with two self-associating species."""
    pairs = [("k_ij", (0, 1), -0.045), ("k_ij", (0, 2), 0.045), ("k_ij", (1, 2), -0.025)]
    for index in ((0, 1), (2, 3)):
        pairs += [("association_energy_over_k", index, 2500.0), ("association_volume", index, 0.05)]
    species = ((30.0, 2.40, 3.20, 200.0), (31.0, 1.08, 3.00, 400.0), (59.0, 2.80, 3.80, 280.0))
    return Mixture(parameters(species, pairs, ((0, 0, 1, 1), (1, 1, 1, 1))))


def reactive_problem(temperature: float, **target: Any) -> Problem:
    return Problem(
        phases=[Phase("L1", kind="liquid", composition_guess=(0.05, 0.90, 0.05)),
                Phase("L2", kind="liquid", composition_guess=(0.05, 0.05, 0.90))],
        T=temperature, P=1.0e5, feed=Amounts({"s0": 0.05, "s1": 0.475, "s2": 0.475}),
        reactions=[Reaction({"s0": -1.0, "s1": -1.0, "s2": 1.0}, name="A+B=C", **target)],
    )


def correlation(identity: str = EOS_STANDARD_STATE_ID) -> ReactionLogPolynomial:
    return ReactionLogPolynomial(-B / T0, B, 0.0, 0.0, T0, 250.0, 450.0, identity)


def quotients(model: Mixture, result: Any) -> list[float]:
    """sum_i nu_i ln a_i per phase, a_i = (x_i rho / (1 mol/m^3)) exp(mu_res_i / RT), from public States."""
    values = []
    for phase in range(2):
        x = result.mole_fractions[3 * phase:3 * phase + 3]
        state = model.state(result.temperature[0], rho=result.molar_densities[phase], x=x)
        values.append(sum(nu * (math.log(x[i] * state.molar_density) + state.residual_chemical_potential_over_rt[i])
                          for i, nu in enumerate(NU)))
    return values


def test_correlated_reaction_target_governs_the_solved_quotient(evidence: Any) -> None:
    model = reactive_model()
    central = solve_equilibrium(model, reactive_problem(T0, correlation=correlation()))
    assert central.success
    for phase, value in enumerate(quotients(model, central)):
        same(evidence, f"e5-10/manifold/value/phase{phase}", value)
        assert abs(value) <= 1.0e-9  # the retained central reaction-affinity bound at K(T0) = 1
    # One kelvin away the correlated solve follows ln K(T); a scalar K(T0) = 1 loses that response.
    warm = T0 + 1.0
    expected = B * (1.0 / warm - 1.0 / T0)
    followed = solve_equilibrium(model, reactive_problem(warm, correlation=correlation()))
    frozen = solve_equilibrium(model, reactive_problem(warm, equilibrium_constant=1.0))
    assert followed.success and frozen.success
    for phase in range(2):
        evidence.compare(f"correlated-quotient/phase{phase}", model=quotients(model, followed)[phase], reference=expected,
                         criterion={"type": "scalar", "rule": "symmetric_close", "atol": 1e-9, "rtol": 0.0})
        assert abs(quotients(model, followed)[phase] - expected) <= 1.0e-9
        assert abs(quotients(model, frozen)[phase]) <= 1.0e-9 < abs(expected)
    with pytest.raises(ParameterError, match="reference-basis-mismatch"):
        solve_equilibrium(model, reactive_problem(T0, correlation=correlation("pressure-basis-P0")))
    with pytest.raises(ParameterError, match="must be one of equilibrium_constant or affinity_j_per_mol or correlation"):
        reactive_problem(T0, correlation=correlation(), equilibrium_constant=1.0)

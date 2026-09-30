"""#187: two independent exchange constants recovered from synthetic component-amount ratios."""
import math
from typing import Any

import numpy as np
import pytest
from scipy.optimize import least_squares  # type: ignore[import-untyped]

from epcsaft import Mixture, _core, regression
from epcsaft.equilibrium import Amounts, Phase, Problem, Reaction, solve_equilibrium
from ..native_inputs import ETHANE, parameters


@pytest.mark.slow
def test_two_exchange_constants_match_known_values_and_finite_difference_fit(evidence: Any) -> None:
    # Six identical neutral ethane pseudo-species. Two fixed-support liquids have
    # separate conserved solvents; A and B exchange independently between them.
    # Equal residual potentials cancel in each exchange law: K = x_org / x_aq.
    p = parameters((ETHANE,) * 6, pairs=[("k_ij", (i, j), 0.0) for i in range(6) for j in range(i + 1, 6)])
    model = Mixture(p)
    known = np.array([0.4, -0.2])

    def problem(ln_k: Any) -> Problem:
        return Problem(T=200.0, P=5.0e6,
            phases=[Phase("aqueous", kind="liquid", support=("s0", "s2", "s4"), composition_guess=(0.8, 0.12, 0.08)),
                    Phase("organic", kind="liquid", support=("s1", "s3", "s5"), composition_guess=(0.8, 0.12, 0.08))],
            feed=Amounts({"s0": 1.0, "s1": 1.0, "s2": 0.3, "s4": 0.2}),
            reactions=[Reaction({"s2": -1.0, "s3": 1.0}, name="exchange_a", equilibrium_constant=math.exp(ln_k[0])),
                       Reaction({"s4": -1.0, "s5": 1.0}, name="exchange_b", equilibrium_constant=math.exp(ln_k[1]))])

    def observations(ln_k: Any) -> list[Any]:
        return [regression.observation(p, "amount_ratio", problem(ln_k), phase="organic", species=(organic,),
            denominator_phase="aqueous", denominator_species=(aqueous,),
            observed=1.0, form="log_ratio", mass_basis=False)
            for aqueous, organic in (("s2", "s3"), ("s4", "s5"))]

    state = solve_equilibrium(model, problem(known))
    assert state.success and state.topology_event == _core.TopologyEvent.None_, state.message
    for i, (aqueous, organic, total) in enumerate(((2, 3, 0.3), (4, 5, 0.2))):
        law = math.log(state.mole_fractions[6 + organic] / state.mole_fractions[aqueous])
        assert abs(law - known[i]) <= 1e-9
        balance = state.amounts[0] * state.mole_fractions[aqueous] + state.amounts[1] * state.mole_fractions[6 + organic]
        assert abs(balance - total) <= 1e-9
        evidence.observe(f"reaction-c/law/{i}", model=law, reference=float(known[i]),
                         inputs={"T_K": 200.0, "P_Pa": 5e6, "balance_mol": balance})
    data = []
    for row in observations(known):
        prediction = regression.predict(model, row)
        assert prediction.value is not None and prediction.status == regression.ObservationStatus.Available
        data.append(prediction.value)
    start = np.array([-0.1, 0.1])
    rows = observations(start)
    for row, observed in zip(rows, data):
        row.observed = observed
    coordinates = [regression.reaction_coordinate(name, "a", origin=0.0, scale=1.0, bounds=(-2.0, 2.0))
                   for name in ("exchange_a", "exchange_b")]
    fitted = regression.fit(p, coordinates, rows)
    assert fitted.status == "converged" and fitted.usable and fitted.trial_failures == 0, fitted.message
    assert fitted.covariance.status == "available" and fitted.covariance.rank == 2
    assert fitted.final_cost <= 1e-14 and max(abs(r) for r in fitted.residuals) <= 1e-8

    def objective(ln_k: Any) -> Any:
        values = []
        for row, observed in zip(observations(ln_k), data):
            row.observed = observed
            prediction = regression.predict(model, row)
            assert prediction.residual is not None and prediction.status == regression.ObservationStatus.Available
            values.append(prediction.residual)
        return np.array(values)

    external = least_squares(objective, start, jac="3-point", diff_step=1e-3, bounds=(-2.0, 2.0),
                             xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=40)
    assert external.success and max(abs(r) for r in external.fun) <= 1e-8
    for i, value in enumerate(fitted.physical):
        for label, reference in (("known", known[i]), ("scipy", external.x[i])):
            error = abs(value - reference)
            evidence.compare(f"reaction-c/{label}/{i}", model=value, reference=float(reference), unit="ln K",
                inputs={"T_K": 200.0, "P_Pa": 5e6, "engine_cost": fitted.final_cost,
                        "scipy_cost": external.cost, "engine_residuals": fitted.residuals, "ratios": data},
                criterion={"rule": "symmetric_close", "atol": 1e-6, "rtol": 0.0})
            assert error <= 1e-6


def test_direct_log_constants_recover_coefficients_with_absolute_weights(evidence: Any) -> None:
    from epcsaft.equilibrium import ReactionLogPolynomial, EOS_STANDARD_STATE_ID

    p = parameters((ETHANE,) * 2, pairs=[("k_ij", (0, 1), 0.0)])
    rows = []
    temperatures = (200.0, 250.0, 300.0)
    weights = (4.0, 9.0, 16.0)
    for T in temperatures:
        problem = Problem(T=T, P=5e6, phases=[Phase("liquid", kind="liquid")],
            feed=Amounts({"s0": 0.5, "s1": 0.5}),
            reactions=[Reaction({"s0": -1.0, "s1": 1.0}, name="isomerisation",
                correlation=ReactionLogPolynomial(-1.25, 90.0, 0.0, 0.0, 250.0, 150.0, 350.0, EOS_STANDARD_STATE_ID))])
        row = regression.observation(p, "pressure", problem, observed=-1.5 + 150.0 / T,
                                     form="difference", scale=1.0, mass_basis=False)
        row.kind = regression.ObservationKind.LogEquilibriumConstant
        row.reaction = "isomerisation"
        rows.append(row)
    coordinates = [regression.reaction_coordinate("isomerisation", "a", origin=0.0, scale=1.0),
                   regression.reaction_coordinate("isomerisation", "b", origin=0.0, scale=250.0)]
    fitted = regression.fit(p, coordinates, rows, weights=weights)
    assert fitted.usable and fitted.status == "converged", fitted.message
    assert fitted.covariance.rank == 2 and fitted.covariance.status == "available"
    for k, known in enumerate((-1.5, 150.0)):
        assert abs(fitted.physical[k] - known) <= 1e-6
        evidence.compare(f"ln-k/python/recovery/{k}", model=fitted.physical[k], reference=known,
                         unit="native", criterion={"rule": "symmetric_close", "atol": 1e-6, "rtol": 0.0})
    for i, (T, weight) in enumerate(zip(temperatures, weights)):
        for k, basis in enumerate((1.0, 250.0 / T)):
            expected = math.sqrt(weight) * basis
            actual = fitted.optimizer_jacobian[i * 2 + k]
            assert abs(actual - expected) <= 1e-12
            evidence.compare(f"ln-k/python/weighted-column/{i}/{k}", model=actual, reference=expected,
                             unit="1", criterion={"rule": "symmetric_close", "atol": 1e-12, "rtol": 0.0})

"""#187: two independent exchange constants recovered from synthetic component-amount ratios."""
import math
from typing import Any

from epcsaft import Mixture, _core, regression
from epcsaft.equilibrium import Amounts, Phase, Problem, Reaction, solve_equilibrium
from tests.native_inputs import ETHANE, budget, parameters


def test_two_exchange_constants_recover_analytic_amount_ratios(evidence: Any) -> None:
    # Six identical neutral ethane pseudo-species. Two fixed-support liquids have
    # separate conserved solvents; A and B exchange independently between them.
    # Equal residual potentials cancel in each exchange law: K = x_org / x_aq.
    p = parameters((ETHANE,) * 6, pairs=[("k_ij", (i, j), 0.0) for i in range(6) for j in range(i + 1, 6)])
    model = Mixture(p)
    known = (math.log(2.0), -math.log(2.0))

    def problem(ln_k: Any) -> Problem:
        return Problem(T=200.0, P=5.0e6,
            phases=[Phase("aqueous", kind="liquid", support=("s0", "s2", "s4"), composition_guess=(0.8, 0.12, 0.08)),
                    Phase("organic", kind="liquid", support=("s1", "s3", "s5"), composition_guess=(0.8, 0.12, 0.08))],
            feed=Amounts({"s0": 1.0, "s1": 1.0, "s2": 0.3, "s4": 0.3}),
            reactions=[Reaction({"s2": -1.0, "s3": 1.0}, name="exchange_a", equilibrium_constant=math.exp(ln_k[0])),
                       Reaction({"s4": -1.0, "s5": 1.0}, name="exchange_b", equilibrium_constant=math.exp(ln_k[1]))])

    def observations(ln_k: Any) -> list[Any]:
        return [regression.observation(p, "amount_ratio", problem(ln_k), phase="organic", species=(organic,),
            denominator_phase="aqueous", denominator_species=(aqueous,),
            observed=ratio, form="log_ratio", mass_basis=False)
            for aqueous, organic, ratio in (("s2", "s3", 2.0), ("s4", "s5", 0.5))]

    # #234 analytic targets are inputs, independent of Engine predictions.
    # The 1e-6 coefficient limit is retained from #187; closure at 1e-9 and
    # log-ratio residuals at 1e-8 independently constrain the solved quantities.
    rows = observations((-0.1, 0.1))
    coordinates = [regression.reaction_coordinate(name, "a", origin=0.0, scale=1.0, bounds=(-2.0, 2.0))
                   for name in ("exchange_a", "exchange_b")]
    fitted = regression.fit(p, coordinates, rows, controls=budget())
    assert fitted.status == "converged" and fitted.usable and fitted.trial_failures == 0, fitted.message
    assert fitted.covariance.status == "available" and fitted.covariance.rank == 2
    assert fitted.final_cost <= 1e-14 and max(abs(r) for r in fitted.residuals) <= 1e-8

    for i, value in enumerate(fitted.physical):
        evidence.compare(f"reaction-c/analytic/{i}", model=value, reference=known[i], unit="ln K",
            inputs={"T_K": 200.0, "P_Pa": 5e6, "engine_cost": fitted.final_cost,
                    "engine_residuals": fitted.residuals, "ratios": [2.0, 0.5],
                    "rank": fitted.covariance.rank},
            criterion={"rule": "symmetric_close", "atol": 1e-6, "rtol": 0.0})
        assert abs(value - known[i]) <= 1e-6

    # Check the fitted equilibrium against independent amounts and fractions.
    state = solve_equilibrium(model, problem(fitted.physical))
    assert state.success and state.topology_event == _core.TopologyEvent.None_, state.message
    for phase, fractions in enumerate(((10/13, 1/13, 2/13), (10/13, 2/13, 1/13))):
        assert abs(state.amounts[phase] - 1.3) <= 1e-9
        evidence.compare(f"reaction-c/phase-amount/{phase}", model=state.amounts[phase],
            reference=1.3, unit="mol", criterion={"rule": "symmetric_close", "atol": 1e-9, "rtol": 0.0})
        for species, fraction in zip(((0, 2, 4), (1, 3, 5))[phase], fractions):
            actual = state.mole_fractions[6 * phase + species]
            assert abs(actual - fraction) <= 1e-9
            assert abs(state.amounts[phase] * actual - 1.3 * fraction) <= 1e-9
            evidence.compare(f"reaction-c/fraction/{phase}/{species}", model=actual,
                reference=fraction, unit="1", criterion={"rule": "symmetric_close", "atol": 1e-9, "rtol": 0.0})
    for i, (aqueous, organic, ratio) in enumerate(((2, 3, 2.0), (4, 5, 0.5))):
        aqueous_amount = state.amounts[0] * state.mole_fractions[aqueous]
        organic_amount = state.amounts[1] * state.mole_fractions[6 + organic]
        law = math.log(state.mole_fractions[6 + organic] / state.mole_fractions[aqueous])
        assert abs(law - known[i]) <= 1e-9
        assert abs(aqueous_amount + organic_amount - 0.3) <= 1e-9
        assert abs(math.log((organic_amount / aqueous_amount) / ratio)) <= 1e-8
        evidence.observe(f"reaction-c/law/{i}", model=law, reference=known[i],
            inputs={"T_K": 200.0, "P_Pa": 5e6, "balance_mol": aqueous_amount + organic_amount,
                    "amount_ratio": organic_amount / aqueous_amount})


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
    fitted = regression.fit(p, coordinates, rows, weights=weights, controls=budget())
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

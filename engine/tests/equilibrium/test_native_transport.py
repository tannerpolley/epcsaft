"""B1-B3 (#56-#59): public reaction targets, solved-state actions and result diagnostics."""
import gc
import math
from typing import Any, cast

import pytest

from epcsaft import ActiveParameter, Mixture, ParameterError, ParameterFamily, PropertyObservable, Pure, Stop, _core
from epcsaft.equilibrium import (
    EOS_STANDARD_STATE_ID,
    Amounts,
    Free,
    Phase,
    Pinned,
    Problem,
    Reaction,
    ReactionLogPolynomial,
    SolvedStateActionDirection,
    SolvedStateActionRequest,
    SolvedStateActionStatus,
    SolvedStateObservable,
    SolvedStateObservableKind as Kind,
    compile_problem,
    solve_equilibrium,
    solved_state_actions,
)
from tests.native_inputs import ETHANE, ROWS, aqueous_electrolyte, p3_forms, parameters, record, same, surface_model

CRITERION = (5.0e-6, 5.0e-5, 5.0e-4)  # the native E5 order-k action criteria

NU = (-1.0, -1.0, 1.0)
B, T0 = 1500.0, 298.15  # manufactured EOS-basis ln K = b (1/T - 1/T0), K(T0) = 1


def reactive_model() -> Mixture:
    """equilibrium_support::synthetic_reactive_ternary(): A + B = C with two self-associating species."""
    pairs = [("k_ij", (0, 1), -0.045), ("k_ij", (0, 2), 0.045), ("k_ij", (1, 2), -0.025)]
    for index in ((0, 1), (2, 3)):
        pairs += [("association_energy_over_k", index, 2500.0), ("association_volume", index, 0.05)]
    species = ((30.0, 2.40, 3.20, 200.0), (31.0, 1.08, 3.00, 400.0), (59.0, 2.80, 3.80, 280.0))
    # The native ternary has no s0/s1 cross edge; declared explicitly since #569.
    none = [("no_cross_association", (0, 1), 0.0, "dimensionless")]
    return Mixture(parameters(species, pairs, ((0, 0, 1, 1), (1, 1, 1, 1)), extra=none))


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
        state = model.state(result.temperature, rho=result.molar_densities[phase], x=x)
        terms = []
        for i, nu in enumerate(NU):
            if nu == 0:
                continue
            potential = state.residual_chemical_potential_over_rt[i]
            assert potential is not None
            terms.append(nu * (math.log(x[i] * state.molar_density) + potential))
        values.append(sum(terms))
    return values


def test_correlated_reaction_target_governs_the_solved_quotient(evidence: Any) -> None:
    model = reactive_model()
    central = solve_equilibrium(model, reactive_problem(T0, correlation=correlation()))
    assert central.success
    for phase, value in enumerate(quotients(model, central)):
        name = f"e5-10/manifold/value/phase{phase}"  # formed here from public States: the affinity bound applies
        same(evidence, name, value, ROWS[name], 1.0e-9, 0.0)
        assert abs(value) <= 1.0e-9  # the retained central reaction-affinity bound at K(T0) = 1
    # One kelvin away the correlated solve follows ln K(T); a scalar K(T0) = 1 loses that response.
    warm = T0 + 1.0
    expected = B * (1.0 / warm - 1.0 / T0)
    followed = solve_equilibrium(model, reactive_problem(warm, correlation=correlation()))
    frozen = solve_equilibrium(model, reactive_problem(warm, equilibrium_constant=1.0))
    assert followed.success and frozen.success
    for phase in range(2):
        same(evidence, f"correlated-quotient/phase{phase}", quotients(model, followed)[phase], expected, 1e-9, 0.0)
        assert abs(quotients(model, frozen)[phase]) <= 1.0e-9 < abs(expected)
    with pytest.raises(ParameterError, match="reference-basis-mismatch"):
        solve_equilibrium(model, reactive_problem(T0, correlation=correlation("pressure-basis-P0")))
    with pytest.raises(ParameterError, match="admits at most one equilibrium_constant, affinity_j_per_mol or correlation"):
        reactive_problem(T0, correlation=correlation(), equilibrium_constant=1.0)


def actions(model: Any, problem: Problem, observables: list[Any], directions: list[tuple[float, float, list[float]]]) -> Any:
    """One compile, one central solve, one exact action batch; directions are (dT / K, dP / Pa, dn_feed / mol)."""
    central = solve_equilibrium(model, problem)
    assert central.success
    request = SolvedStateActionRequest(observables, [SolvedStateActionDirection(*item, []) for item in directions])
    return solved_state_actions(compile_problem(model, problem), central, request)


def test_thermal_actions_of_the_logistic_isomerisation_are_analytic(evidence: Any) -> None:
    """E5-3: identical species with ln K(T) give x_B = K / (1 + K); bound 2e-10 |reference| + 1e-14."""
    T = 350.0
    a = math.log(2.0) - (1200.0 / T + 0.4 * math.log(T) + 0.002 * T)  # K(350 K) = 2
    law = ReactionLogPolynomial(a + 0.4 * math.log(300.0), 1200.0, 0.4, 0.002, 300.0, 250.0, 450.0, EOS_STANDARD_STATE_ID)
    model = Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), 0.0)]))
    problem = Problem(phases=[Phase("V", kind="vapor")], T=T, P=2.0e5, feed=Amounts({"s0": 0.5, "s1": 0.5}),
                      reactions=[Reaction({"s0": -1.0, "s1": 1.0}, correlation=law)])
    p = 2.0 / 3.0
    L1, L2, L3 = (-1200.0 / T**2 + 0.4 / T + 0.002, 2400.0 / T**3 - 0.4 / T**2, -7200.0 / T**4 + 0.8 / T**3)
    analytic = (p * (1 - p) * L1, p * (1 - p) * (L2 + (1 - 2 * p) * L1**2),
                p * (1 - p) * (L3 + 3 * (1 - 2 * p) * L1 * L2 + (1 - 6 * p + 6 * p * p) * L1**3))
    for order in (1, 2, 3):
        result = actions(model, problem, [SolvedStateObservable(Kind.PhaseMoleFraction, 0, 1)], [(1.0, 0.0, [0.0, 0.0])] * order).results[0]
        assert result.order == order
        same(evidence, f"logistic/A{order}", result.action, analytic[order - 1], 2.0e-10 * abs(analytic[order - 1]) + 1.0e-14, 0.0)
        same(evidence, f"e5-3/isomerisation/A{order}", result.action)


def test_log_activity_actions_follow_the_correlated_reaction(evidence: Any) -> None:
    """E5-10 on A + B = C: sum_i nu_i D^k ln a_i equals d^k ln K / dT^k within CRITERION[k-1] sum |nu_i D^k ln a_i|."""
    model = reactive_model()
    problem = reactive_problem(T0, correlation=correlation())
    logs = [SolvedStateObservable(Kind.PhaseLogActivity, phase, i) for phase in range(2) for i in range(3)]
    dT, dP = (1.0, 0.0, [0.0] * 3), (0.0, 1.0e5, [0.0] * 3)
    cases = (("temperature-temperature-temperature", [dT, dT, dT], -6.0 * B / T0**4), ("pressure-temperature", [dP, dT], 0.0))
    for label, directions, reference in cases:
        batch = actions(model, problem, logs, directions)
        assert batch.diagnostics.rank == batch.diagnostics.state_dimension
        for phase in range(2):
            terms = [nu * batch.results[3 * phase + i].action for i, nu in enumerate(NU)]
            bound = CRITERION[len(directions) - 1] * sum(abs(term) for term in terms)
            name = f"e5-10/manifold/{label}/phase{phase}"
            if reference == 0.0:  # the solve resolves this exact zero: twice the recording route's distance from it (#605)
                same(evidence, name, sum(terms), ROWS[name], 2 * max(abs(ROWS[f"e5-10/manifold/{label}/phase{p}"]) for p in range(2)))
            else:
                same(evidence, name, sum(terms))
            same(evidence, f"manifold/{label}/phase{phase}", sum(terms), reference, bound, 0.0)
            assert abs(sum(terms) + terms[2] - reference) > bound  # nu = (-1, -1, 2) fails the identity's own criterion


def test_ionic_feed_and_caloric_actions_match_native(evidence: Any) -> None:
    """E5-5: water, Na+ and Cl- along charge-safe T, P and feed directions, orders one to three."""
    model = Mixture(aqueous_electrolyte(), thermochemistry=record(*[p3_forms()["polynomial"]] * 3))
    problem = Problem(phases=[Phase("L", kind="liquid")], T=298.15, P=1.0e5, feed=Amounts({"s0": 0.98, "s1": 0.01, "s2": 0.01}))
    directions = [(0.1, 1000.0, [-0.002, 0.001, 0.001]), (-0.05, -500.0, [0.001, -0.0005, -0.0005]), (0.02, 250.0, [-0.0004, 0.0002, 0.0002])]
    names = ("ResidualEnthalpy", "ResidualIsobaricHeatCapacity", "TotalEnthalpy", "TotalIsobaricHeatCapacity")
    observables = [SolvedStateObservable(Kind.PhaseProperty, 0, 0, getattr(PropertyObservable, name)) for name in names]
    observables += [SolvedStateObservable(kind, 0, i) for kind in (Kind.PhaseLogActivity, Kind.PhaseComponentAmount) for i in range(3)]
    ids = [f"property/{name}" for name in names] + [f"log_activity[{i}]" for i in range(3)]
    for order in (1, 2, 3):
        batch = actions(model, problem, observables, directions[:order])
        for name, result in zip(ids, batch.results):
            same(evidence, f"e5-5/electrolyte/phase[0].{name}/A{order}", result.action)
        for i, result in enumerate(batch.results[7:]):  # one phase holds the whole feed: n_i = n_feed,i exactly
            assert result.unit == "mol" and result.status == SolvedStateActionStatus.Available
            same(evidence, f"ionic-amount/{i}/A{order}", result.action, directions[0][2][i] if order == 1 else 0.0)
    # Without the thermal record the total row is typed and empty while its residual sibling keeps the same action.
    partial = actions(Mixture(aqueous_electrolyte()), problem, observables[:3], directions[:1]).results
    assert partial[2].status == SolvedStateActionStatus.MissingThermochemistry
    assert partial[2].value is None and partial[2].action is None and partial[2].diagnostic_message
    same(evidence, "partial/residual-sibling", partial[0].action, ROWS["e5-5/electrolyte/phase[0].property/ResidualEnthalpy/A1"])
    with pytest.raises(ParameterError, match="solved-state action direction"):  # feed width is the component count
        actions(model, problem, observables[:1], [(0.0, 0.0, [0.0, 0.0])])


def surface_flash(pressure: float | Free, vapor: Any = None, **options: Any) -> Problem:
    """The C++ example's associating binary at 353.15 K: its two-phase flash, or its bubble point with a pinned vapor."""
    return Problem(phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=vapor)], T=353.15, P=pressure,
                   feed=Amounts({"s0": 0.25, "s1": 0.75}), **options)


def test_two_phase_actions_match_the_cpp_example(evidence: Any) -> None:
    """Declaration-order phases; thermal and radial-feed actions of amount, ln a, x, total h and density."""
    observables = [SolvedStateObservable(Kind.PhaseComponentAmount, 0, 1), SolvedStateObservable(Kind.PhaseLogActivity, 0, 1),
                   SolvedStateObservable(Kind.PhaseMoleFraction, 0, 1),
                   SolvedStateObservable(Kind.PhaseProperty, 0, 0, PropertyObservable.TotalEnthalpy),
                   SolvedStateObservable(Kind.PhaseMolarDensity, 0)]
    thermal = actions(surface_model(), surface_flash(1.8e4), observables, [(1.0, 0.0, [0.0, 0.0])])
    radial = actions(surface_model(), surface_flash(1.8e4), observables, [(0.0, 0.0, [0.25, 0.75])])
    for index, unit in enumerate(("mol", "1", "1", "J/mol", "mol/m^3")):
        assert thermal.results[index].unit == unit
        # equivalence standard, CODING_STANDARDS.md, PR #311; retain the absolute floor near zero.
        name = f"solved/thermal/{index}"
        same(evidence, name, thermal.results[index].action, ROWS[name], max(1e-12, 1e-10 * abs(ROWS[name])), 0.0)
        same(evidence, f"solved/radial/{index}", radial.results[index].action)
    same(evidence, "solved/condition", thermal.diagnostics.condition_estimate)
    amount, activity = radial.results[0], radial.results[1]  # Euler: D n_i[n_feed] = n_i and D ln a_i[n_feed] = 0
    assert abs(amount.action - amount.value) <= CRITERION[0] * abs(amount.value) and abs(activity.action) <= 1.0e-12
    # Lifetime: the compiled system alone keeps its native model alive after the Python owner is gone.
    compiled = compile_problem(surface_model(), surface_flash(1.8e4))
    gc.collect()
    request = SolvedStateActionRequest(observables[4:], [SolvedStateActionDirection(1.0, 0.0, [0.0, 0.0], [])])
    central = solve_equilibrium(surface_model(), surface_flash(1.8e4))
    same(evidence, "solved/thermal/4", solved_state_actions(compiled, central, request).results[0].action)


def test_result_diagnostics_match_the_cpp_example(evidence: Any) -> None:
    model = surface_model()
    split = solve_equilibrium(model, surface_flash(1.8e4))
    for name, values in (("amounts", split.amounts), ("mole_fractions", split.mole_fractions),
                         ("multipliers", split.multipliers), ("constraint_multipliers", split.constraint_multipliers)):
        for index, value in enumerate(values):
            same(evidence, f"flash/{name}/{index}", value)
    # Two bases at the returned state: reduced-energy balance coordinates and reconstructed
    # extensive-objective duals (here material-balance and normalization rows, in rows order).
    implemented = [row.kind for row in split.rows if row.kind not in (_core.RowKind.Chemical, _core.RowKind.Pressure)]
    assert implemented == [_core.RowKind.Moiety] * 2 + [_core.RowKind.Normalization] * 2
    assert len(split.multipliers) == 2 and len(split.constraint_multipliers) == len(implemented)
    assert all(abs(a - b) > 1.0 for a, b in zip(split.multipliers, split.constraint_multipliers))
    same(evidence, "bubble/pressure", solve_equilibrium(model, surface_flash(Free(3.0e4), Pinned(0.0))).pressure)
    # An unattainable row tolerance: no success, the achieved finite residuals and the owner's stop stay inspectable.
    floor = solve_equilibrium(model, surface_flash(Free(3.0e4), Pinned(0.0), row_tolerance=1.0e-18))
    same(evidence, "floor/residual", floor.convergence.residual)
    assert not floor.success and floor.requested_tolerance_met is False and floor.convergence.stop == Stop.converged_to_floor
    assert floor.solver_status == "STOP_AT_TINY_STEP"
    assert floor.residuals and all(math.isfinite(value) for value in floor.residuals) and floor.temperature is None and floor.pressure is None


def test_active_parameter_and_thermal_actions_match_native(evidence: Any) -> None:
    """E5-1: a vapor-liquid ternary along T, P, feed and (k_01, m_0) directions; E5-7: pure-vapor calorics."""
    species = ((16.043, 1.0, 3.7039, 150.03), ETHANE, (44.097, 2.0020, 3.6184, 208.11))
    active = [ActiveParameter(ParameterFamily.Kij, (1, 3)), ActiveParameter(ParameterFamily.SegmentCount, (0,))]
    model = Mixture(parameters(species, [("k_ij", pair, 0.0) for pair in ((0, 1), (0, 2), (1, 2))]), active_parameters=active)
    problem = Problem(phases=[Phase("V", kind="vapor"), Phase("L", kind="liquid")], T=180.0, P=1.0e5,
                      feed=Amounts({f"s{i}": 1.0 / 3.0 for i in range(3)}))
    moves = [(0.2, 2000.0, [-0.01, 0.006, 0.004], [0.001, 0.01]), (-0.1, -1000.0, [0.003, -0.004, 0.001], [-0.0005, 0.004]),
             (0.05, 500.0, [0.002, 0.001, -0.003], [0.0002, -0.002])]
    kinds = ((Kind.PhaseComponentAmount, "amount"), (Kind.PhaseLogActivity, "log_activity"))
    observables = [SolvedStateObservable(kind, phase, i) for kind, _ in kinds for phase in range(2) for i in range(3)]
    names = [f"phase[{phase}].{name}[{i}]" for _, name in kinds for phase in range(2) for i in range(3)]
    central, compiled = solve_equilibrium(model, problem), compile_problem(model, problem)
    for order in (1, 2, 3):
        request = SolvedStateActionRequest(observables, [SolvedStateActionDirection(*item) for item in moves[:order]])
        for name, result in zip(names, solved_state_actions(compiled, central, request).results):
            same(evidence, f"e5-1/neutral-ternary/{name}/A{order}", result.action)
    # Euler homogeneity (E5-10): n_feed radial with a non-radial v; amounts are degree one, ln a degree zero.
    radial, skew = SolvedStateActionDirection(0.0, 0.0, [1.0 / 3.0] * 3, [0.0, 0.0]), SolvedStateActionDirection(0.0, 0.0, [-0.01, 0.006, 0.004], [0.0, 0.0])
    first = solved_state_actions(compiled, central, SolvedStateActionRequest(observables, [skew])).results
    mixed = solved_state_actions(compiled, central, SolvedStateActionRequest(observables, [radial, skew])).results
    for name, one, two in zip(names, first, mixed):
        slope = cast(float, one.action)
        if "amount" in name:  # D^2 n_i[n_feed, v] = 0
            same(evidence, f"euler/mixed/{name}", two.action, 0.0, 1.0e-12 * max(1.0, abs(slope)), 0.0)
        else:  # D^2 ln a_i[n_feed, v] = -D ln a_i[v], with a nonzero first-order response
            assert abs(slope) > 1.0e-3
            same(evidence, f"euler/mixed/{name}", two.action, -slope, CRITERION[1] * abs(slope), 0.0)
    fixed = SolvedStateActionRequest(observables[:1], [SolvedStateActionDirection(*moves[0][:3], [0.0, 0.0])])
    lost = solved_state_actions(compiled, central, fixed).results[0].action  # without the parameter increments
    assert lost is not None and abs(lost - ROWS["e5-1/neutral-ternary/phase[0].amount[0]/A1"]) > 1.0e-4
    # A same-layout result of another problem keeps its mapped value and is typed, not differentiated.
    other = solved_state_actions(compile_problem(model, Problem(phases=problem.phases, T=181.0, P=1.0e5, feed=problem.feed)),
                                 central, fixed).results[0]
    assert other.status == SolvedStateActionStatus.DefiningEquationsNotSatisfied and other.action is None
    same(evidence, "typed/defining-equations/value", other.value, central.amounts[0] * central.mole_fractions[0])

    caloric_names = ("IdealEnthalpy", "TotalEnthalpy", "IdealIsobaricHeatCapacity", "TotalIsobaricHeatCapacity")
    calorics = [SolvedStateObservable(Kind.PhaseProperty, 0, 0, getattr(PropertyObservable, name)) for name in caloric_names]
    vapor = Pure(parameters((ETHANE,)), thermochemistry=record(p3_forms()["shomate"]))
    for order in (1, 2, 3):
        steps = [(20.0, 0.0, [0.0]), (-10.0, 0.0, [0.0]), (5.0, 0.0, [0.0])][:order]
        batch = actions(vapor, Problem(phases=[Phase("V", kind="vapor")], T=320.0, P=1.0e5, feed=Amounts({"s0": 1.0})), calorics, steps)
        for name, result in zip(caloric_names, batch.results):
            # equivalence standard, CODING_STANDARDS.md, PR #311; retain the absolute floor near zero.
            row = f"e5-7/thermal-320/phase[0].property/{name}/A{order}"
            same(evidence, row, result.action, ROWS[row], max(1e-12, 1e-10 * abs(ROWS[row])), 0.0)
    for temperature, status in ((250.0, SolvedStateActionStatus.ThermalBoundary), (249.0, SolvedStateActionStatus.OutsideInterval)):
        edge = actions(vapor, Problem(phases=[Phase("V", kind="vapor")], T=temperature, P=1.0e5, feed=Amounts({"s0": 1.0})),
                       calorics, [(20.0, 0.0, [0.0])]).results
        assert all(item.status == status and item.action is None for item in edge)
        assert (edge[0].value is not None) == (temperature == 250.0)  # the owned endpoint keeps its central value

"""Record-anchored calorics on a referenced reactive state (EqID: reference_species_calorics).

Identical ethane 2A ⇌ B: A carries an ideal-gas record, B has none and is completed from the
EOS-basis reaction enthalpy and heat capacity. The oracle takes C' and C'' from the pure-fluid
identities RT² dΛ/dT = -(h_res + RT) and d²Λ/dT² = -(cp_res + R)/(RT²) + 2(h_res + RT)/(RT³),
not from the Engine's reference ladder.
"""
import math
from typing import Any

import pytest

from epcsaft import (
    ActiveParameter, GAS_CONSTANT_J_PER_MOL_K as R, ParameterFamily, IdealCorrelation, IdealInterval, IdealPolynomial, Mixture,
    PropertyObservable as Property, Pure, ThermochemistryRecord,
)
from epcsaft.equilibrium import (
    Amounts, NeutralReference, Phase, Problem, Reaction, ReactionLogPolynomial, ReactionReference,
    ReferenceSourceBasis as Basis, SolvedStateActionDirection as Direction, SolvedStateActionRequest,
    SolvedStateActionStatus as Status, SolvedStateObservable, SolvedStateObservableKind as Kind,
    compile_problem, solve_equilibrium, solved_state_actions, solved_state_first_actions,
)
from tests.native_inputs import ETHANE, aqueous_electrolyte, parameters, same

T0, P0 = 260.0, 3.0e6
LAW = ReactionLogPolynomial(math.log(0.3), 200.0, 3.0, -0.004, 300.0, 250.0, 350.0, "analytic-source")
CP_A = (52.0, 0.12)  # ideal-gas cp(T) = 52 + 0.12 (T - 300) J/(mol K), anchored h = s = 0 at 300 K
FEED = {"s0": 2.0, "s1": 0.2}


def _record(*species: list[IdealInterval]) -> ThermochemistryRecord:
    return ThermochemistryRecord(1.0e5, list(species))


A_RECORD = [IdealInterval(200.0, 400.0, True, True, IdealCorrelation(IdealPolynomial(CP_A, 300.0, 0.0, 0.0), 1.0e5))]


def _problem(temperature: float, reference_pressure: float | None, feed: dict[str, float] = FEED) -> Problem:
    return Problem(
        T=temperature, P=P0, feed=Amounts(feed), phases=[Phase("L", kind="liquid")],
        neutral_reference=NeutralReference(solvent="s0", solute_molality_ratios={"s1": 1.0},
                                           standard_molality_mol_per_kg=1.0),
        reactions=[Reaction({"s0": -2.0, "s1": 1.0}, correlation=LAW, reference=ReactionReference(
            Basis.RawMoleFractionInfiniteDilution, reference_pressure))],
    )


def _model(*record: list[IdealInterval]) -> Mixture:
    return Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), 0.0)]), thermochemistry=_record(*record))


H, CP = (SolvedStateObservable(Kind.PhaseProperty, 0, 0, p) for p in (Property.TotalEnthalpy,
                                                                     Property.TotalIsobaricHeatCapacity))
AMOUNT = SolvedStateObservable(Kind.PhaseAmount, 0)
LN_A = SolvedStateObservable(Kind.PhaseLogActivity, 0, 0)


def _values(model: Mixture, problem: Problem, observables: list[SolvedStateObservable],
            direction: Direction) -> list[Any]:
    result = solve_equilibrium(model, problem)
    assert result.success, result.message
    return list(solved_state_actions(compile_problem(model, problem), result,
                                     SolvedStateActionRequest(observables, [direction])).results)


@pytest.mark.parametrize("reference_pressure", [None, 5.0e6])
def test_reference_calorics_recover_pure_fluid_chain(evidence: Any, reference_pressure: float | None) -> None:
    """Frozen h and cp: omitting C' (source Δh), C'', the source ΔCp l2 or a density/pressure RT
    shift each moves h or cp by >1 J/mol or >0.1 J/(mol K)."""
    t, rt = T0, R * T0
    pure = Pure(parameters((ETHANE,)))
    system = pure.state(t, P=P0, phase="liquid")
    ref = pure.state(t, P=reference_pressure or P0, phase="liquid")
    assert ref.residual_isobaric_heat_capacity is not None and system.residual_isobaric_heat_capacity is not None
    lam1 = -(ref.residual_enthalpy + rt) / (rt * t)
    lam2 = -(ref.residual_isobaric_heat_capacity + R) / (rt * t) + 2 * (ref.residual_enthalpy + rt) / (rt * t * t)
    l1 = -LAW.b / t**2 + LAW.c / t + LAW.d
    l2 = 2 * LAW.b / t**3 - LAW.c / t**2
    c1, c2 = -lam1, -lam2  # C = Δν Λ with Δν = -1
    dt = t - 300.0
    h_a, cp_a = CP_A[0] * dt + CP_A[1] * dt * dt / 2, CP_A[0] + CP_A[1] * dt
    h_b = 2 * (h_a - rt) + rt * t * (l1 + c1) + rt
    cp_b = 2 * (cp_a - R) + 2 * rt * (l1 + c1) + rt * t * (l2 + c2) + R
    model = _model(A_RECORD, [])
    result = solve_equilibrium(model, _problem(t, reference_pressure))
    assert result.success, result.message
    x_a, x_b = result.mole_fractions
    items = _values(model, _problem(t, reference_pressure), [H, CP], Direction(0.0, 0.0, [0.0, 0.0], []))
    name = f"reference-calorics/{'system' if reference_pressure is None else 'fixed'}"
    inputs = {"T_K": t, "P_Pa": P0, "P_ref_Pa": reference_pressure}
    same(evidence, name + "/h", items[0].value, x_a * h_a + x_b * h_b + system.residual_enthalpy,
         rtol=1e-9, atol=1e-6, inputs=inputs, basis="Identity", method="pure-fluid-residual-enthalpy",
         claim="total h carries the EOS-basis reaction enthalpy RT^2 (l1 + C')")
    same(evidence, name + "/cp", items[1].value, x_a * cp_a + x_b * cp_b + system.residual_isobaric_heat_capacity,
         rtol=1e-9, atol=1e-8, inputs=inputs, basis="Identity", method="pure-fluid-residual-heat-capacity",
         claim="frozen cp carries 2RT(l1 + C') + RT^2 (l2 + C'')")


def test_reference_enthalpy_actions_match_resolves_and_gibbs_helmholtz(evidence: Any) -> None:
    """Feed and T actions of nH equal centered re-solves, and the apparent partial enthalpy
    satisfies h_A - d(nH)/dn_A = RT² d ln f_A/dT at fixed P and feed (the absorption heat)."""
    model = _model(A_RECORD, [])
    problem = _problem(T0, None)
    result = solve_equilibrium(model, problem)
    assert result.success, result.message
    compiled = compile_problem(model, problem)
    feed, thermal = solved_state_first_actions(compiled, result, [H, AMOUNT, LN_A], [
        Direction(0.0, 0.0, [1.0, 0.0], []), Direction(1.0, 0.0, [0.0, 0.0], [])])
    assert all(item.status == Status.Available for item in feed.results + thermal.results)
    h, n, dh, dn, th, tn, t_ln_a = (math.nan if v is None else v for v in (  # all Available, asserted above
        feed.results[0].value, feed.results[1].value, feed.results[0].action, feed.results[1].action,
        thermal.results[0].action, thermal.results[1].action, thermal.results[2].action))
    partial, heat = n * dh + h * dn, n * th + h * tn
    dt = T0 - 300.0
    h_a = CP_A[0] * dt + CP_A[1] * dt * dt / 2
    inputs = {"T_K": T0, "P_Pa": P0}
    same(evidence, "reference-calorics/gibbs-helmholtz", h_a - partial,
         R * T0 * T0 * t_ln_a + R * T0, rtol=1e-8, atol=0.0, inputs=inputs,
         basis="Identity", method="gibbs-helmholtz-log-activity",
         claim="caloric apparent partial enthalpy equals the equilibrium's T-derivative of ln f")

    def extensive(temperature: float, extra: float) -> float:
        values = _values(model, _problem(temperature, None, {"s0": 2.0 + extra, "s1": 0.2}), [H, AMOUNT],
                         Direction(0.0, 0.0, [0.0, 0.0], []))
        return float(values[0].value * values[1].value)

    for label, action, step, shift in (("feed", partial, 0.2, (0.0, 1.0)), ("T", heat, 0.05, (1.0, 0.0))):
        coarse, fine = ((extensive(T0 + s * shift[0], s * shift[1]) - extensive(T0 - s * shift[0], -s * shift[1]))
                        / (2 * s) for s in (2 * step, step))
        # O(step^2) toward the action above re-solve rounding (a few ulps of nH over the step; along the feed
        # nH is linear here to rounding, so its truncation is below this floor at every step on both starts)
        assert abs(fine - action) <= max(abs(coarse - action) / 2, 16 * math.ulp(n * h) / step)
        same(evidence, f"reference-calorics/d(nH)/{label}/resolve", action, fine, rtol=1e-6, atol=0.0,
             inputs={**inputs, "step": step}, basis="Identity", method="centered-resolve-difference",
             claim="extensive enthalpy action equals the re-solved slope")


def test_reference_calorics_refuse_missing_orders_and_undetermined_species() -> None:
    model = _model(A_RECORD, [])
    problem = _problem(T0, None)
    result = solve_equilibrium(model, problem)
    compiled = compile_problem(model, problem)
    residual = SolvedStateObservable(Kind.PhaseProperty, 0, 0, Property.ResidualEnthalpy)
    thermal, pressure = solved_state_first_actions(compiled, result, [H, CP, residual], [
        Direction(1.0, 0.0, [0.0, 0.0], []), Direction(0.0, 1.0e5, [0.0, 0.0], [])])
    assert thermal.results[0].status == Status.Available
    assert thermal.results[1].status == Status.ReferenceUnavailable and thermal.results[1].value is not None
    assert pressure.results[0].status == Status.ReferenceUnavailable and pressure.results[0].value is not None
    assert pressure.results[2].status == Status.Available and pressure.results[2].action is not None
    twice = Direction(0.0, 0.0, [1.0, 0.0], [])
    second = solved_state_actions(compiled, result, SolvedStateActionRequest([H], [twice, twice])).results[0]
    assert second.status == Status.ReferenceUnavailable and second.value is not None
    for record in ((A_RECORD, A_RECORD), ([], [])):  # overdetermined-inconsistent; element gauge undetermined
        other = _model(*record)
        solved = solve_equilibrium(other, problem)
        item = solved_state_actions(compile_problem(other, problem), solved, SolvedStateActionRequest(
            [H, residual], [Direction(0.0, 0.0, [0.0, 0.0], [])])).results
        assert item[0].status == Status.ReferenceUnavailable and item[0].value is None
        assert item[1].status == Status.Available and item[1].value is not None


def test_reference_parameter_actions_refuse_missing_orders() -> None:
    """A k_ij direction reaches only order-one non-caloric actions: h would need d2C/dTdtheta."""
    model = Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), 0.0)]), thermochemistry=_record(A_RECORD, []),
                    active_parameters=[ActiveParameter(ParameterFamily.Kij, (1, 2))])
    problem = _problem(T0, None)
    result = solve_equilibrium(model, problem)
    compiled = compile_problem(model, problem)
    fraction = SolvedStateObservable(Kind.PhaseMoleFraction, 0, 1)
    k = Direction(0.0, 0.0, [0.0, 0.0], [1.0])
    caloric, composition = solved_state_actions(compiled, result, SolvedStateActionRequest([H, fraction], [k])).results
    assert caloric.status == Status.ReferenceUnavailable and caloric.value is not None and caloric.action is None
    assert composition.status == Status.Available and composition.action is not None
    second = solved_state_actions(compiled, result, SolvedStateActionRequest([fraction], [k, k])).results[0]
    assert second.status == Status.ReferenceUnavailable and second.value is not None and second.action is None


def _ionic(temperature: float, reference_pressure: float | None) -> Problem:
    law = ReactionLogPolynomial(math.log(1e-8), 0.0, 0.0, 0.0, 298.15, 250.0, 350.0, "analytic-source")
    return Problem(
        T=temperature, P=1.0e5, feed=Amounts({"s0": 1.0, "s1": 1e-6, "s2": 1e-6}), phases=[Phase("L", kind="liquid")],
        neutral_reference=NeutralReference(solvent="s0", solute_molality_ratios={"s1": 1.0, "s2": 1.0},
                                           standard_molality_mol_per_kg=1.0),
        reactions=[Reaction({"s0": -1.0, "s1": 1.0, "s2": 1.0}, correlation=law, reference=ReactionReference(
            Basis.RawMoleFractionInfiniteDilution, reference_pressure))],
    )


@pytest.mark.parametrize("reference_pressure", [None, 2.0e5])
def test_ionic_reference_calorics_reproduce_native_rows_and_resolves(
    evidence: Any, reference_pressure: float | None,
) -> None:
    """Aqueous NaCl, H2O ⇌ Na+ + Cl-: water carries a record, the ions are completed up to the
    charge gauge. The installed wheel reproduces the native total h and cp rows, and the
    equilibrium heat capacity d(nH)/dT|_P, which carries the ionic second slope, equals
    centered re-solves."""
    water = [IdealInterval(200.0, 500.0, True, True, IdealCorrelation(
        IdealPolynomial([33.59], 298.15, -241826.0, 188.83), 1.0e5))]
    model = Mixture(aqueous_electrolyte(), thermochemistry=_record(water, [], []))
    zero = Direction(0.0, 0.0, [0.0, 0.0, 0.0], [])
    items = _values(model, _ionic(298.15, reference_pressure), [H, CP], zero)
    name = f"reference-calorics/ionic/{'system' if reference_pressure is None else 'fixed'}"
    same(evidence, name + "/h", items[0].value)
    same(evidence, name + "/cp", items[1].value)
    thermal = _values(model, _ionic(298.15, reference_pressure), [H, AMOUNT], Direction(1.0, 0.0, [0.0, 0.0, 0.0], []))
    assert all(item.status == Status.Available for item in thermal)
    h, n, dh, dn = (math.nan if v is None else v for v in (
        thermal[0].value, thermal[1].value, thermal[0].action, thermal[1].action))

    def extensive(temperature: float) -> float:
        values = _values(model, _ionic(temperature, reference_pressure), [H, AMOUNT], zero)
        return float(values[0].value * values[1].value)

    action = n * dh + h * dn
    coarse, fine = ((extensive(298.15 + s) - extensive(298.15 - s)) / (2 * s) for s in (0.1, 0.05))
    assert abs(fine - action) <= abs(coarse - action) / 2  # O(step^2) toward the action
    same(evidence, name + "/d(nH)/T/resolve", action, fine, rtol=1e-6, atol=0.0,
         inputs={"T_K": 298.15, "P_Pa": 1.0e5, "P_ref_Pa": reference_pressure, "step": 0.05},
         basis="Identity", method="centered-resolve-difference",
         claim="electrolyte equilibrium heat capacity equals the re-solved enthalpy slope")

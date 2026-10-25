"""A nonzero reaction mole change falsifies density and molality reference signs."""
from functools import cache
import math
from typing import Any

import pytest

from epcsaft import ActiveParameter, GAS_CONSTANT_J_PER_MOL_K as R, Mixture, ParameterFamily, Pure, _core
from epcsaft.equilibrium import (
    Amounts, Continuation, ContinuationStop, Free, NeutralReference, Phase, Pinned, Problem, Reaction, ReactionLogPolynomial,
    ReactionReference, ReferenceSourceBasis as Basis, SolvedStateActionDirection,
    SolvedStateActionRequest, SolvedStateActionStatus as Status, SolvedStateObservable,
    SolvedStateObservableKind as Kind,
    compile_problem, continue_equilibrium, problem_from_mapping, solve_equilibrium, solved_state_actions,
    solved_state_first_actions,
)
from tests.native_inputs import ETHANE, aqueous_electrolyte, parameters, same


@pytest.mark.parametrize("reference_pressure", [None, 2.0e5])
def test_species_reference_switch_preserves_density_target(evidence: Any, reference_pressure: float | None) -> None:
    """EqID: neutral_reaction_reference_transfer; transform K with independent native fugacities."""
    from epcsaft.equilibrium import SpeciesReference
    temperature, pressure = 260.0, 1.0e5
    for shift in (0.0, 0.01):
        solute = (ETHANE[0], ETHANE[1], ETHANE[2], ETHANE[3] + shift)
        model = Mixture(parameters((ETHANE, solute), [("k_ij", (0, 1), 0.0)]))
        source_p = reference_pressure or pressure
        liquid = model.state(temperature, P=source_p, x=(1.0, 0.0), phase="liquid")
        vapor = model.state(temperature, P=source_p, x=(0.0, 1.0), phase="vapor")
        assert liquid.log_fugacity_coefficient is not None and vapor.log_fugacity_coefficient is not None
        liquid_coefficient, vapor_coefficient = liquid.log_fugacity_coefficient[1], vapor.log_fugacity_coefficient[1]
        assert liquid_coefficient is not None and vapor_coefficient is not None
        delta = vapor_coefficient - liquid_coefficient
        results = []
        for mixed in (False, True):
            reference = NeutralReference(solvent="s0", solute_molality_ratios={"s1": 1.0},
                                         standard_molality_mol_per_kg=1.0,
                                         species_reference={"s1": SpeciesReference.PureGas} if mixed else None)
            law = ReactionLogPolynomial(math.log(0.3) - (delta if mixed else 0.0), 0.0, 0.0, 0.0,
                                        300.0, 250.0, 350.0, "manufactured-reference-switch")
            problem = Problem(T=temperature, P=pressure, feed=Amounts({"s0": 2.0, "s1": 0.2}),
                              phases=[Phase("L", kind="liquid")], neutral_reference=reference,
                              reactions=[Reaction({"s0": -2.0, "s1": 1.0}, correlation=law,
                                                  reference=ReactionReference(Basis.RawMoleFractionInfiniteDilution,
                                                                              reference_pressure))])
            result = solve_equilibrium(model, problem)
            assert result.success, result.message
            assert max(abs(r) for r in result.residuals) < 1e-10
            results.append(result)
        left, right = results
        c_left, c_right = (r.neutral_reference_diagnostics[0].terminal_contraction for r in results)
        assert c_left is not None and c_right is not None
        name = f"species-reference/switch/{reference_pressure}/{shift}"
        same(evidence, name + "/density-target", c_right - delta, c_left, rtol=0.0, atol=1e-10,
             method="native-fugacity-convention-switch", claim="ln K + C is invariant after a reference switch")
        same(evidence, name + "/x_B", right.mole_fractions[1], left.mole_fractions[1], rtol=0.0, atol=1e-10,
             method="native-fugacity-convention-switch", claim="transformed references preserve solved composition")


def test_pure_gas_reference_selects_vapor_at_ethane_saturation(evidence: Any) -> None:
    """At saturation the fugacities coincide; the pressure slope distinguishes branches."""
    from epcsaft.equilibrium import SpeciesReference
    temperature = 180.0
    pure = Pure(parameters((ETHANE,)))
    saturation = solve_equilibrium(pure, pure.saturation(temperature))
    assert saturation.success, saturation.message
    pressure = saturation.pressure
    model = Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), 0.0)]))
    problem = Problem(T=temperature, P=Free(pressure), feed=Amounts({"s0": 2.0, "s1": 0.2}),
                      phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))],
                      neutral_reference=NeutralReference(solvent="s0", solute_molality_ratios={"s1": 1.0},
                                                         standard_molality_mol_per_kg=1.0,
                                                         species_reference={"s1": SpeciesReference.PureGas}),
                      reactions=[Reaction({"s0": -2.0, "s1": 1.0}, correlation=ReactionLogPolynomial(
                          math.log(0.3), 0.0, 0.0, 0.0, 180.0, 170.0, 200.0, "manufactured-saturation"),
                          reference=ReactionReference(Basis.RawMoleFractionInfiniteDilution))])
    result = solve_equilibrium(model, problem)
    assert result.success, result.message
    slope = result.neutral_reference_diagnostics[0].terminal_contraction_pressure_derivative
    assert slope is not None
    vapor = pure.state(temperature, P=result.pressure, phase="vapor")
    liquid = pure.state(temperature, P=result.pressure, phase="liquid")
    # Remove the solvent's -2 v_liquid/(RT) contribution from the reaction contraction.
    gas_slope = slope + 2.0 / (liquid.molar_density * R * temperature)
    expected = 1.0 / (vapor.molar_density * R * temperature)
    same(evidence, "species-reference/ethane-vapor-pressure-slope", gas_slope, expected,
         rtol=1e-9, atol=0.0, method="native-saturation-branch-volume",
         claim="the pure-gas source-pressure slope is v_vapor/(RT)")
    assert abs(gas_slope - 1.0 / (liquid.molar_density * R * temperature)) > expected / 2.0


def test_pure_gas_reference_refuses_compressed_liquid_ethane() -> None:
    """At 180 K and 3 MPa pure ethane has only its liquid root; a declared vapor must refuse."""
    from epcsaft.equilibrium import SpeciesReference
    model = Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), 0.0)]))
    problem = Problem(T=180.0, P=3.0e6, feed=Amounts({"s0": 2.0, "s1": 0.2}), phases=[Phase("L", kind="liquid")],
                      neutral_reference=NeutralReference(solvent="s0", solute_molality_ratios={"s1": 1.0},
                                                         standard_molality_mol_per_kg=1.0,
                                                         species_reference={"s1": SpeciesReference.PureGas}),
                      reactions=[Reaction({"s0": -2.0, "s1": 1.0}, correlation=ReactionLogPolynomial(
                          math.log(0.3), 0.0, 0.0, 0.0, 180.0, 170.0, 200.0, "manufactured-compressed-liquid"),
                          reference=ReactionReference(Basis.RawMoleFractionInfiniteDilution))])
    with pytest.raises(_core.ParameterError, match="pure reference species 1: vapor branch not established by the finite search"):
        solve_equilibrium(model, problem)


@pytest.mark.parametrize("temperature,pressure", [(260.0, 3e6), (280.0, 5e6)])
@pytest.mark.parametrize("common", [False, True])
def test_neutral_reference_recovers_mole_fraction_law(
    evidence: Any, temperature: float, pressure: float, common: bool,
) -> None:
    """Identical EOS species: 2A ⇌ B gives x_B = K_x(1-x_B)^2 at any regular T/P."""
    model = Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), 0.0)]))
    constant = 0.3
    expected = 4 * constant / (1 + math.sqrt(1 + 4 * constant)) ** 2
    ln_k = math.log(constant) - (math.log(ETHANE[0] / 1000) if common else 0.0)
    basis = Basis.CommonMolalityInfiniteDilution if common else Basis.RawMoleFractionInfiniteDilution
    law = ReactionLogPolynomial(ln_k, 0.0, 0.0, 0.0, 300.0, 250.0, 350.0, "analytic-source")
    problem = Problem(
        T=temperature, P=pressure, feed=Amounts({"s0": 2.0, "s1": 0.2}),
        phases=[Phase("L", kind="liquid")],
        neutral_reference=NeutralReference(solvent="s0", solute_molality_ratios={"s1": 1.0},
                                           standard_molality_mol_per_kg=1.0),
        reactions=[Reaction({"s0": -2.0, "s1": 1.0}, correlation=law,
                            reference=ReactionReference(basis))],
    )
    if common:
        problem = problem_from_mapping({
            "T": temperature, "P": pressure, "feed": {"basis": "amounts", "values": {"s0": 2.0, "s1": 0.2}},
            "phases": [{"name": "L", "kind": "liquid"}],
            "neutral_reference": {"solvent": "s0", "solute_molality_ratios": {"s1": 1.0},
                                  "standard_molality_mol_per_kg": 1.0},
            "reactions": [{"stoichiometry": {"s0": -2.0, "s1": 1.0}, "correlation": law,
                           "reference": {"source_basis": "CommonMolalityInfiniteDilution"}}],
        })
    result = solve_equilibrium(model, problem)
    assert result.success, result.message
    inputs = {"T_K": temperature, "P_Pa": pressure, "K_x": constant, "common_molality": common}
    name = f"neutral-reference/{temperature:g}/{common}"
    same(evidence, name + "/x_B", result.mole_fractions[1], expected, rtol=1e-8, atol=0.0,
         inputs=inputs, basis="Identity", method="2A-to-B-analytic",
         claim="density-reference and source molality conversion recover the mole-fraction law")
    # Independent pressure/fugacity identity, not the reference owner's density expression.
    state = model.state(temperature, P=pressure, x=(0.5, 0.5), phase="liquid")
    assert state.log_fugacity_coefficient is not None
    solvent_fugacity = state.log_fugacity_coefficient[0]
    assert solvent_fugacity is not None
    contraction = -solvent_fugacity - math.log(pressure / (R * temperature))
    same(evidence, name + "/contraction", result.neutral_reference_diagnostics[0].terminal_contraction,
         contraction, rtol=1e-9, atol=1e-10, inputs=inputs, basis="Identity",
         method="fugacity-pressure-identity", claim="Lambda includes the standard-density factor")
    direction = SolvedStateActionDirection(0.0, 0.0, [2.0, 0.2], [])
    item = solved_state_actions(
        compile_problem(model, problem), result,
        SolvedStateActionRequest([SolvedStateObservable(Kind.PhaseMoleFraction, 0, 1)], [direction]),
    ).results[0]
    same(evidence, name + "/feed_homogeneity", item.action, 0.0, rtol=0.0, atol=1e-10,
         inputs=inputs, basis="Identity", method="extensive-feed-scaling",
         claim="referenced equilibrium composition remains invariant under total-feed scaling")


@pytest.mark.parametrize("constant", [0.3, 1e-10])
def test_reactive_bubble_recovers_pure_saturation(evidence: Any, constant: float) -> None:
    """Identical EOS species separate source chemistry from the phase-equilibrium oracle."""
    temperature = 180.0
    pure = Pure(parameters((ETHANE,)))
    pure_result = solve_equilibrium(pure, pure.saturation(temperature))
    assert pure_result.success, pure_result.message
    assert max(abs(value) for value in pure_result.residuals) < 1e-10
    saturation = pure_result.pressure
    liquid = pure.state(temperature, P=saturation, phase="liquid")
    vapor = pure.state(temperature, P=saturation, phase="vapor")
    assert liquid.log_fugacity_coefficient is not None
    assert vapor.log_fugacity_coefficient is not None
    liquid_fugacity = liquid.log_fugacity_coefficient[0]
    vapor_fugacity = vapor.log_fugacity_coefficient[0]
    assert liquid_fugacity is not None and vapor_fugacity is not None
    assert abs(liquid_fugacity - vapor_fugacity) < 1e-8
    model = Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), 0.0)]))
    expected = 4 * constant / (1 + math.sqrt(1 + 4 * constant)) ** 2
    name = f"reference-bubble/{constant:g}"
    problem = Problem(
        T=temperature, P=Free(1e5), feed=Amounts({"s0": 2.0, "s1": 0.2}),
        phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))],
        neutral_reference=NeutralReference(solvent="s0", solute_molality_ratios={"s1": 1.0},
                                           standard_molality_mol_per_kg=1.0),
        reactions=[Reaction({"s0": -2.0, "s1": 1.0}, correlation=ReactionLogPolynomial(
            math.log(constant), 0.0, 0.0, 0.0, 300.0, 170.0, 350.0, "analytic-source"),
            reference=ReactionReference(Basis.RawMoleFractionInfiniteDilution))],
    )
    result = solve_equilibrium(model, problem)
    assert result.success, result.message
    assert max(abs(value) for value in result.residuals) < 1e-10
    same(evidence, name + "/pressure", result.pressure, saturation, rtol=1e-7,
         method="native-pure-saturation", claim="identical-species reactive bubble has pure-EOS pressure")
    same(evidence, name + "/liquid-composition", result.mole_fractions[1], expected,
         rtol=1e-8, atol=0.0, method="2A-to-B-analytic", claim="source mole-fraction law at solved pressure")
    reference = result.neutral_reference_diagnostics[0]
    liquid = pure.state(temperature, P=saturation, phase="liquid")
    same(evidence, name + "/reference-pressure", reference.reference_pressure_pa, saturation,
         rtol=1e-7, method="native-pure-saturation", claim="diagnostics use final pressure, not the seed")
    same(evidence, name + "/pressure-first", reference.terminal_contraction_pressure_derivative,
         -1.0 / (liquid.molar_density * R * temperature), rtol=1e-8, atol=0.0,
         method="identical-species-chemical-potential", claim="nonideal reference C_P=-V/(RT)")
    same(evidence, name + "/pressure-second",
         reference.terminal_contraction_pressure_second_derivative,
         1.0 / (liquid.molar_density**2 * R * temperature * liquid.pressure_density_derivative),
         rtol=1e-6, atol=0.0, method="identical-species-chemical-potential",
         claim="reference curvature follows the density response")


def test_free_pressure_mixes_system_and_fixed_pressure_references(evidence: Any) -> None:
    """Identical species: 2A ⇌ B at system P and 2A ⇌ C at fixed P0 share one bubble.

    With C(P) = -ln φ(P) - ln(P/RT), row r satisfies x_r/x_A^2 = K_r exp(C(P_r) - C(P)),
    so only the fixed-P0 row carries the pure-liquid factor between P0 and P_sat.
    """
    temperature, reference_pressure = 180.0, 1.0e6
    pure = Pure(parameters((ETHANE,)))
    model = Mixture(parameters((ETHANE, ETHANE, ETHANE),
                               [("k_ij", (0, 1), 0.0), ("k_ij", (0, 2), 0.0), ("k_ij", (1, 2), 0.0)]))
    constants = (0.3, 0.2)
    law = [ReactionLogPolynomial(math.log(k), 0.0, 0.0, 0.0, 300.0, 170.0, 350.0, "analytic-source")
           for k in constants]
    problem = Problem(
        T=temperature, P=Free(1e5), feed=Amounts({"s0": 2.0, "s1": 0.2, "s2": 0.1}),
        phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))],
        neutral_reference=NeutralReference(solvent="s0", solute_molality_ratios={"s1": 1.0, "s2": 1.0},
                                           standard_molality_mol_per_kg=1.0),
        reactions=[
            Reaction({"s0": -2.0, "s1": 1.0}, correlation=law[0],
                     reference=ReactionReference(Basis.RawMoleFractionInfiniteDilution)),
            Reaction({"s0": -2.0, "s2": 1.0}, correlation=law[1],
                     reference=ReactionReference(Basis.RawMoleFractionInfiniteDilution,
                                                 reference_pressure)),
        ],
    )
    result = solve_equilibrium(model, problem)
    assert result.success and result.pressure is not None, result.message
    saturation = result.pressure

    def contraction(pressure: float) -> float:
        coefficient = pure.state(temperature, P=pressure, phase="liquid").log_fugacity_coefficient
        assert coefficient is not None
        solvent_fugacity = coefficient[0]
        assert solvent_fugacity is not None
        return -solvent_fugacity - math.log(pressure / (R * temperature))

    effective = (constants[0], constants[1] * math.exp(contraction(reference_pressure) - contraction(saturation)))
    total = sum(effective)
    solvent = (-1 + math.sqrt(1 + 4 * total)) / (2 * total)
    inputs = {"T_K": temperature, "P0_Pa": reference_pressure, "K": constants}
    for row, species in enumerate((1, 2)):
        same(evidence, f"reference-mixed-pressure/x{species}", result.mole_fractions[species],
             effective[row] * solvent**2, rtol=1e-8, atol=0.0, inputs=inputs, basis="Identity",
             method="identical-species-two-reaction-analytic",
             claim="free-pressure solve keeps system-pressure and fixed-P0 reference rows distinct")
    diagnostics = result.neutral_reference_diagnostics
    same(evidence, "reference-mixed-pressure/system-pressure", diagnostics[0].reference_pressure_pa,
         saturation, rtol=1e-12, method="solved-pressure", claim="system reference follows solved pressure")
    same(evidence, "reference-mixed-pressure/fixed-pressure", diagnostics[1].terminal_contraction,
         contraction(reference_pressure), rtol=1e-9, atol=1e-10, method="pure-fugacity-identity",
         claim="fixed-P0 reference is independent of the solved pressure")


@pytest.mark.parametrize("reference_pressure", [None, 5.0e6])
def test_temperature_reference_chain_recovers_identical_species_law(
    evidence: Any, reference_pressure: float | None,
) -> None:
    """Identical species: 2A ⇌ B obeys ln K_eff = ln K_src + C(T,P_ref) - C(T,P) with
    RT² dC/dT|_P = h_res + RT and dC/dP = -V/(RT), both from the pure-fluid state owner.

    A system-pressure row cancels the mixture term exactly, so dx_B/dT is source-driven;
    a fixed row keeps the residual-enthalpy difference.  Omitting psi_nT or the density
    slope changes d ln K_eff/dT by ~2e-2 or ~5e-2 per K.
    """
    temperature, pressure = 260.0, 3.0e6
    law = ReactionLogPolynomial(math.log(0.3), 200.0, 3.0, -0.004, 300.0, 250.0, 350.0, "analytic-source")
    ln_k = (law.a + law.b / temperature + law.c * math.log(temperature / law.reference_temperature)
            + law.d * temperature)
    source_slope = -law.b / temperature**2 + law.c / temperature + law.d
    pure = Pure(parameters((ETHANE,)))
    system = pure.state(temperature, P=pressure, phase="liquid")
    fixed = pure.state(temperature, P=reference_pressure or pressure, phase="liquid")
    rt = R * temperature

    def contraction(state: Any, p: float) -> float:
        assert state.log_fugacity_coefficient is not None
        return -float(state.log_fugacity_coefficient[0]) - math.log(p / rt)

    ln_k_eff = ln_k + contraction(fixed, reference_pressure or pressure) - contraction(system, pressure)
    slope_t = source_slope + (fixed.residual_enthalpy - system.residual_enthalpy) / (rt * temperature)
    slope_p = 0.0 if reference_pressure is None else 1.0 / (system.molar_density * rt)
    k = math.exp(ln_k_eff)
    x = 4 * k / (1 + math.sqrt(1 + 4 * k)) ** 2
    dx_dlnk = k * (1 - x) ** 2 / (1 + 2 * k * (1 - x))
    model = Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), 0.0)]))
    problem = Problem(
        T=temperature, P=pressure, feed=Amounts({"s0": 2.0, "s1": 0.2}),
        phases=[Phase("L", kind="liquid")],
        neutral_reference=NeutralReference(solvent="s0", solute_molality_ratios={"s1": 1.0},
                                           standard_molality_mol_per_kg=1.0),
        reactions=[Reaction({"s0": -2.0, "s1": 1.0}, correlation=law, reference=ReactionReference(
            Basis.RawMoleFractionInfiniteDilution, reference_pressure))],
    )
    result = solve_equilibrium(model, problem)
    assert result.success, result.message
    compiled = compile_problem(model, problem)
    observable = [SolvedStateObservable(Kind.PhaseMoleFraction, 0, 1)]
    name = f"reference-temperature/{'system' if reference_pressure is None else 'fixed'}"
    inputs = {"T_K": temperature, "P_Pa": pressure, "P_ref_Pa": reference_pressure}
    same(evidence, name + "/x_B", result.mole_fractions[1], x, rtol=1e-8, atol=0.0, inputs=inputs,
         basis="Identity", method="2A-to-B-analytic", claim="reference value at the declared pressure")
    for label, d_t, d_p in (("T", 1.0, 0.0), ("T+P", 1.0, 2.0e5)):
        direction = SolvedStateActionDirection(d_t, d_p, [0.0, 0.0], [])
        item = solved_state_actions(compiled, result, SolvedStateActionRequest(observable, [direction])).results[0]
        same(evidence, f"{name}/dx_B/{label}", item.action, dx_dlnk * (slope_t * d_t + slope_p * d_p),
             rtol=1e-8, atol=0.0, inputs=inputs, basis="Identity",
             method="identical-species-residual-enthalpy",
             claim="first-order T/P action carries dC/dT|_P and C_P of the reference")
    # Driven feed A3 has native numerical witnesses; moving-T A3 is outside #225.


def _ionic_problem(temperature: float, pressure: float, reference_pressure: float | None) -> Problem:
    """Manufactured formula reaction H2O ⇌ Na+ + Cl- (charge neutral, not atom balanced) on
    the native equilibrium_checks aqueous case."""
    law = ReactionLogPolynomial(math.log(1e-8), 0.0, 0.0, 0.0, 298.15, 250.0, 350.0, "analytic-source")
    return Problem(
        T=temperature, P=pressure, feed=Amounts({"s0": 1.0, "s1": 1e-6, "s2": 1e-6}),
        phases=[Phase("L", kind="liquid")],
        neutral_reference=NeutralReference(solvent="s0", solute_molality_ratios={"s1": 1.0, "s2": 1.0},
                                           standard_molality_mol_per_kg=1.0),
        reactions=[Reaction({"s0": -1.0, "s1": 1.0, "s2": 1.0}, correlation=law, reference=ReactionReference(
            Basis.RawMoleFractionInfiniteDilution, reference_pressure))],
    )


@pytest.mark.parametrize("reference_pressure", [None, 2.0e5])
def test_ionic_temperature_action_matches_resolved_difference(
    evidence: Any, reference_pressure: float | None,
) -> None:
    """The native temperature-slope case: first-order T and T/P actions of x_Na equal centered
    re-solves; the step pair shows the O(step^2) truncation below the 1e-6 gate."""
    temperature, pressure = 298.15, 1.0e5
    model = Mixture(aqueous_electrolyte())

    def sodium(d_t: float, d_p: float) -> float:
        result = solve_equilibrium(model, _ionic_problem(temperature + d_t, pressure + d_p, reference_pressure))
        assert result.success, result.message
        return float(result.mole_fractions[1])

    problem = _ionic_problem(temperature, pressure, reference_pressure)
    result = solve_equilibrium(model, problem)
    assert result.success, result.message
    compiled = compile_problem(model, problem)
    observable = [SolvedStateObservable(Kind.PhaseMoleFraction, 0, 1)]
    name = f"reference-temperature/ionic/{'system' if reference_pressure is None else 'fixed'}"
    for label, d_t, d_p in (("T", 1.0, 0.0), ("T+P", 1.0, 1.0e5)):
        direction = SolvedStateActionDirection(d_t, d_p, [0.0, 0.0, 0.0], [])
        item = solved_state_actions(compiled, result, SolvedStateActionRequest(observable, [direction])).results[0]
        step = 0.05
        coarse, fine = ((sodium(h * d_t, h * d_p) - sodium(-h * d_t, -h * d_p)) / (2 * h) for h in (2 * step, step))
        assert item.action is not None
        assert abs(fine - item.action) <= abs(coarse - item.action) / 2  # O(step^2) toward the action
        same(evidence, f"{name}/dx_Na/{label}", item.action)  # the native row of the same solve
        same(evidence, f"{name}/dx_Na/{label}/resolve", item.action, fine, rtol=1e-6, atol=0.0,
             inputs={"T_K": temperature, "P_Pa": pressure, "P_ref_Pa": reference_pressure, "step": step},
             basis="Identity", method="centered-resolve-difference",
             claim="first-order referenced T/P action equals the re-solved composition slope")


@pytest.mark.parametrize("reference_pressure", [None, 2.0e5])
def test_ionic_mixed_temperature_actions_match_differenced_first_actions(
    evidence: Any, reference_pressure: float | None,
) -> None:
    """D2 n_Na[e_T, e] for e = water feed, P and T equals the first action D n_Na[e] differenced
    along T over re-solves (Richardson, O(step^2) shown). The system-pressure row fails without
    s_TP (EqID: standard_state_transfer_temperature_pressure_derivative) and T-T without s_TT."""
    temperature, pressure = 298.15, 1.0e5
    model = Mixture(aqueous_electrolyte())
    observable = [SolvedStateObservable(Kind.PhaseComponentAmount, 0, 1)]
    thermal = SolvedStateActionDirection(1.0, 0.0, [0.0, 0.0, 0.0], [])

    @cache
    def solved(d_t: float) -> tuple[Any, Any]:
        problem = _ionic_problem(temperature + d_t, pressure, reference_pressure)
        result = solve_equilibrium(model, problem)
        assert result.success, result.message
        return compile_problem(model, problem), result

    def action(d_t: float, *directions: SolvedStateActionDirection) -> float:
        compiled, result = solved(d_t)
        assert result.success, result.message
        item = solved_state_actions(compiled, result,
                                    SolvedStateActionRequest(observable, list(directions))).results[0]
        assert item.status == Status.Available and item.action is not None, item.diagnostic_message
        return float(item.action)

    name = f"reference-temperature/ionic/{'system' if reference_pressure is None else 'fixed'}"
    for label, other in (("feed", SolvedStateActionDirection(0.0, 0.0, [1.0, 0.0, 0.0], [])),
                         ("P", SolvedStateActionDirection(0.0, 1.0e5, [0.0, 0.0, 0.0], [])), ("T", thermal)):
        exact = action(0.0, thermal, other)
        coarse, fine = ((action(h, other) - action(-h, other)) / (2 * h) for h in (0.2, 0.1))
        # O(step^2) toward the action; the 1e-6 term is the asserted resolution. The near-zero P action
        # (-3.03e-13) scatters by up to 3.2e-19 for row tolerances 1e-10 to 1e-13, its floating-point floor (#314)
        assert abs(fine - exact) <= abs(coarse - exact) / 2 + 1e-6 * abs(exact)
        same(evidence, f"{name}/d2n_Na/T,{label}/resolve", exact, (4 * fine - coarse) / 3, rtol=1e-6, atol=1e-18,
             inputs={"T_K": temperature, "P_Pa": pressure, "P_ref_Pa": reference_pressure, "steps": [0.2, 0.1]},
             basis="Identity", method="richardson-differenced-first-action",
             claim="referenced mixed temperature second action equals the T derivative of the first action")


@pytest.mark.parametrize("reference_pressure", [None, 2.0e5])
def test_referenced_continuation_follows_the_cold_solves(reference_pressure: float | None) -> None:
    """The continuation tangent carries the reference T/P chain (the two reference pressures move
    d ln x_Na/dq by ~2%): it equals centered cold re-solves, and the T/P path lands on them."""
    model = Mixture(aqueous_electrolyte())
    start, end = (298.15, 1.0e5), (318.15, 3.0e5)

    def at(q: float) -> Any:
        result = solve_equilibrium(model, _ionic_problem(
            (1 - q) * start[0] + q * end[0], (1 - q) * start[1] + q * end[1], reference_pressure))
        assert result.success, result.message
        return result

    central, step = at(0.0), 1e-3
    lever = [end[0] - start[0], end[1] - start[1], 0.0, 0.0, 0.0]
    tangent = _core.sensitivity(compile_problem(model, _ionic_problem(*start, reference_pressure)),
                                central.coordinates, lever)
    resolved = [(u - v) / (2 * step) for u, v in zip(at(step).coordinates, at(-step).coordinates)]
    assert max(abs(t - r) for t, r in zip(tangent, resolved)) <= 1e-6  # re-solve noise ~4e-8
    result = continue_equilibrium(
        model, _ionic_problem(*start, reference_pressure), _ionic_problem(*end, reference_pressure),
        Continuation(initial_step=0.25, minimum_step=0.25, max_steps=4, max_retries_per_step=1),
    )
    assert result.success and result.stop == ContinuationStop.reached_target, result.message
    assert result.path_values[-1] == 1.0
    for q, state in zip(result.path_values, result.states):
        assert max(abs(a / b - 1.0) for a, b in zip(state.mole_fractions, at(q).mole_fractions)) <= 1e-8


@pytest.mark.parametrize("pure_solvent", [False, True])
def test_steep_permittivity_temperature_action_matches_resolved_difference(evidence: Any, pure_solvent: bool) -> None:
    """A steep declared permittivity, d ln(eps T)/d ln T ~ -75, makes the T slope's sqrt(m) ladder
    tail 73 times the constant-eps one, past 5e-5 at the 1e-14 mol/kg rung. The sqrt(m) Richardson
    limit (#353) admits it, and the T action of x_Na equals centred re-solves."""
    from epcsaft.equilibrium import SpeciesReference
    model = Mixture(aqueous_electrolyte(
        ("permittivity_term_coefficient", (0, 0, 1), -20.0, "1 / kelvin"),
        ("permittivity_reference_temperature", (0,), 298.15, "kelvin")))
    temperature = 298.15

    def problem(d_t: float) -> Problem:
        item = _ionic_problem(temperature + d_t, 1.0e5, None)
        if pure_solvent:
            assert item.neutral_reference is not None
            item.neutral_reference.species_reference = {"s0": SpeciesReference.PureGas}
        return item

    def sodium(d_t: float) -> float:
        result = solve_equilibrium(model, problem(d_t))
        assert result.success, result.message
        return float(result.mole_fractions[1])

    result = solve_equilibrium(model, problem(0.0))
    assert result.success, result.message
    thermal = SolvedStateActionDirection(1.0, 0.0, [0.0, 0.0, 0.0], [])
    item = solved_state_actions(compile_problem(model, problem(0.0)), result, SolvedStateActionRequest(
        [SolvedStateObservable(Kind.PhaseMoleFraction, 0, 1)], [thermal])).results[0]
    assert item.status == Status.Available and item.action is not None, item.diagnostic_message
    # eps moves 1.3% per 0.05 K: the centred difference keeps 7e-5 at 0.025 K, its Richardson pair 2e-8.
    step = 0.025
    coarse, fine = ((sodium(h) - sodium(-h)) / (2 * h) for h in (2 * step, step))
    assert abs(fine - item.action) <= abs(coarse - item.action) / 2  # O(step^2) toward the action
    same(evidence, f"reference-temperature/ionic/steep-permittivity/{'pure-gas' if pure_solvent else 'default'}/resolve",
         item.action, (4 * fine - coarse) / 3, rtol=1e-6, atol=0.0,
         inputs={"T_K": temperature, "P_Pa": 1.0e5, "dlneps_dT_per_K": -20.0 / 78.09, "steps": [2 * step, step]},
         basis="Identity", method="richardson-centered-resolve-difference",
         claim="steep-permittivity T action from the sqrt(m) Richardson reference equals the re-solved slope")


def _kij_problem(mode: str, shift: float = 0.0) -> Problem:
    """Near-identical ethane 2A ⇌ B at 180 K: a system-P row at 3 MPa, a fixed 5 MPa row, or a
    free-pressure bubble (system-P row). `shift` moves the source ln K for the dx/d ln K probe."""
    law = ReactionLogPolynomial(math.log(0.3) + shift, 0.0, 0.0, 0.0, 300.0, 170.0, 350.0, "analytic-source")
    bubble = mode == "bubble"
    problem = Problem(
        T=180.0, P=Free(1e5) if bubble else 3.0e6, feed=Amounts({"s0": 2.0, "s1": 0.2}),
        phases=[Phase("L", kind="liquid")] + ([Phase("V", kind="vapor", amount=Pinned(0.0))] if bubble else []),
        neutral_reference=NeutralReference(solvent="s0", solute_molality_ratios={"s1": 1.0},
                                           standard_molality_mol_per_kg=1.0),
        reactions=[Reaction({"s0": -2.0, "s1": 1.0}, correlation=law, reference=ReactionReference(
            Basis.RawMoleFractionInfiniteDilution, 5.0e6 if mode == "fixed" else None))],
    )
    return problem


@pytest.mark.parametrize("family", ["k_ij", "m_solvent"])
@pytest.mark.parametrize("mode", ["system", "fixed", "bubble"])
def test_parameter_reference_chain_matches_resolves(evidence: Any, mode: str, family: str) -> None:
    """EqID standard_state_transfer_parameter_derivative: the theta action of x_B carries
    ds/dtheta = -dC/dtheta of the reference, whose infinite-dilution B depends on theta. Centered
    re-solves rebuild the model, the reference ladder and the solve (Richardson over 2h, h removes
    the O(h²) truncation, ~1e-4 at h = 1e-3). k_AB moves only the dilute solute, so its omitted
    reference term (dx_B/d ln K_src)(dC/dk) is about twice the action. The solvent segment count
    also moves the terminal volume: dropping C_V V_theta = -C_P p_theta would shift the action by
    (dx_B/d ln K_src) C_P p_theta, measured about 180 times the action. Both parameters stay active, so
    one set of reference jets serves both directions (#152); the batch equals separate requests."""
    center = 0.02 if family == "k_ij" else ETHANE[1]

    def model(theta: float) -> Mixture:
        kij, segments = (theta, ETHANE[1]) if family == "k_ij" else (0.02, theta)
        return Mixture(parameters(((ETHANE[0], segments, *ETHANE[2:]), ETHANE), [("k_ij", (0, 1), kij)]),
                       active_parameters=[ActiveParameter(ParameterFamily.Kij, (1, 2)), ActiveParameter(ParameterFamily.SegmentCount, (0,))])

    def solve(theta: float, shift: float = 0.0) -> Any:
        result = solve_equilibrium(model(theta), _kij_problem(mode, shift))
        assert result.success, result.message
        return result

    problem = _kij_problem(mode)
    result = solve(center)
    compiled, observable = compile_problem(model(center), problem), [SolvedStateObservable(Kind.PhaseMoleFraction, 0, 1)]
    directions = [SolvedStateActionDirection(0.0, 0.0, [0.0, 0.0], unit) for unit in ([1.0, 0.0], [0.0, 1.0])]
    batch = [item.results[0] for item in solved_state_first_actions(compiled, result, observable, directions)]
    alone = [solved_state_actions(compiled, result, SolvedStateActionRequest(observable, [d])).results[0] for d in directions]
    assert [item.action for item in batch] == [item.action for item in alone]
    item = batch[0 if family == "k_ij" else 1]
    assert item.status == Status.Available and item.action is not None, item.diagnostic_message
    step = 5.0e-4
    far, near = ([solve(center + h), solve(center - h)] for h in (2 * step, step))
    coarse, fine = ((pair[0].mole_fractions[1] - pair[1].mole_fractions[1]) / (2 * h)
                    for pair, h in ((far, 2 * step), (near, step)))
    assert abs(fine - item.action) <= abs(coarse - item.action) / 2  # O(step^2) toward the action
    inputs = {"T_K": 180.0, "mode": mode, family: center, "step": step}
    same(evidence, f"reference-parameter/{mode}/dx_B/d{family}/resolve", item.action, (4 * fine - coarse) / 3,
         rtol=1e-6, atol=0.0, inputs=inputs, basis="Identity", method="richardson-centered-resolve-difference",
         claim=f"first-order referenced {family} action equals the re-solved composition slope")
    if family == "k_ij":
        contraction = [solved.neutral_reference_diagnostics[0].terminal_contraction for solved in near]
        dc_dk = (contraction[0] - contraction[1]) / (2 * step)
        dx_dlnk = (solve(center, 1e-4).mole_fractions[1] - solve(center, -1e-4).mole_fractions[1]) / 2e-4
        assert abs(dx_dlnk * dc_dk) > 0.5 * abs(item.action)  # measured ~2x the action: omission flips its sign

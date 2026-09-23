"""A nonzero reaction mole change falsifies density and molality reference signs."""
import math
from typing import Any

import pytest

from epcsaft import GAS_CONSTANT_J_PER_MOL_K as R, Mixture, Pure
from epcsaft.equilibrium import (
    Amounts, Free, NeutralReference, Phase, Pinned, Problem, Reaction, ReactionLogPolynomial,
    ReactionReference, ReferenceSourceBasis as Basis, SolvedStateActionDirection,
    SolvedStateActionRequest, SolvedStateObservable, SolvedStateObservableKind as Kind,
    compile_problem, problem_from_mapping, solve_equilibrium, solved_state_actions,
    solved_state_first_actions,
)
from tests.native_inputs import ETHANE, aqueous_electrolyte, parameters, same


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
    contraction = -state.log_fugacity_coefficient[0] - math.log(pressure / (R * temperature))
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
    lower, upper = 30000.0, 150000.0
    for _ in range(40):
        pressure = (lower + upper) / 2
        liquid = pure.state(temperature, P=pressure, phase="liquid")
        vapor = pure.state(temperature, P=pressure, phase="vapor")
        assert liquid.log_fugacity_coefficient is not None
        assert vapor.log_fugacity_coefficient is not None
        if liquid.log_fugacity_coefficient[0] > vapor.log_fugacity_coefficient[0]:
            lower = pressure
        else:
            upper = pressure
    saturation = (lower + upper) / 2
    assert liquid.log_fugacity_coefficient is not None and vapor.log_fugacity_coefficient is not None
    assert abs(liquid.log_fugacity_coefficient[0] - vapor.log_fugacity_coefficient[0]) < 1e-8
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
    same(evidence, name + "/pressure", result.pressure[0], saturation, rtol=1e-7,
         method="pure-fugacity-bisection", claim="identical-species reactive bubble has pure-EOS pressure")
    same(evidence, name + "/liquid-composition", result.mole_fractions[1], expected,
         rtol=1e-8, atol=0.0, method="2A-to-B-analytic", claim="source mole-fraction law at solved pressure")
    reference = result.neutral_reference_diagnostics[0]
    liquid = pure.state(temperature, P=saturation, phase="liquid")
    same(evidence, name + "/reference-pressure", reference.reference_pressure_pa, saturation,
         rtol=1e-7, method="pure-fugacity-bisection", claim="diagnostics use final pressure, not the seed")
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
    assert result.success, result.message
    saturation = result.pressure[0]

    def contraction(pressure: float) -> float:
        coefficient = pure.state(temperature, P=pressure, phase="liquid").log_fugacity_coefficient
        assert coefficient is not None
        return -coefficient[0] - math.log(pressure / (R * temperature))

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
        return -state.log_fugacity_coefficient[0] - math.log(p / rt)

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
    twice = SolvedStateActionDirection(1.0, 0.0, [0.0, 0.0], [])
    second = solved_state_actions(compiled, result, SolvedStateActionRequest(observable, [twice, twice])).results[0]
    assert second.status.name == "ReferenceUnavailable" and second.action is None


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
        return result.mole_fractions[1]

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


def test_refused_temperature_slope_keeps_value_and_feed_direction() -> None:
    """A steep declared permittivity, d ln eps/d ln T ~ -76, multiplies the slope's terminal
    change past 5e-5 while the value criterion holds: only the T direction is refused."""
    model = Mixture(aqueous_electrolyte(
        ("permittivity_term_coefficient", (0, 0, 1), -20.0, "1 / kelvin"),
        ("permittivity_reference_temperature", (0,), 298.15, "kelvin")))
    problem = _ionic_problem(298.15, 1.0e5, None)
    result = solve_equilibrium(model, problem)
    assert result.success, result.message
    observable = [SolvedStateObservable(Kind.PhaseMoleFraction, 0, 1)]
    thermal, feed = solved_state_first_actions(compile_problem(model, problem), result, observable, [
        SolvedStateActionDirection(1.0, 0.0, [0.0, 0.0, 0.0], []),
        SolvedStateActionDirection(0.0, 0.0, [1.0, 0.0, 0.0], [])])
    assert thermal.results[0].status.name == "ReferenceUnavailable" and thermal.results[0].action is None
    assert feed.results[0].status.name == "Available" and feed.results[0].action is not None
    assert thermal.results[0].value == feed.results[0].value == result.mole_fractions[1]

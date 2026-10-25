"""#87/#90 trace-state reliability: the analytic trace solution below the former composition box
and the intensive accuracy of trace rows after the amount-weighted optimizer stop; #154/#157 cold
reacting-liquid starts that recover a manufactured speciation with trace water ions; #159 stability
decisions that trace curvature below rounding cannot flip."""
import json
import math
from pathlib import Path
from typing import Any

import pytest


from epcsaft import (
    GAS_CONSTANT_J_PER_MOL_K, IdealCorrelation, IdealInterval, IdealPolynomial, IdealShomate, Mixture, Parameters,
    PropertyObservable, ThermochemistryRecord,
)
from epcsaft.equilibrium import (
    Amounts, EOS_STANDARD_STATE_ID, Free, NeutralReference, Phase, Problem, Reaction, ReactionLogPolynomial,
    ReactionReference, ReferenceSourceBasis, SolvedStateActionDirection, SolvedStateActionRequest,
    SolvedStateActionStatus, SolvedStateObservable, SolvedStateObservableKind as Kind,
    compile_problem, solve_equilibrium, solved_state_actions, solved_state_first_actions,
)
from tests.native_inputs import ETHANE, parameters, same
from tests.packets import parameters_json


def test_trace_ethanol_retains_two_liquids(evidence: Any) -> None:
    # #650/#652 failing declaration: trace inventory must not drain an initially enriched phase away.
    model = Mixture(Parameters.from_json(parameters_json("khudaida-2026-figure-2-electrolyte-lle", 7),
                                        components=("water", "ethanol", "isobutanol")))
    problem = Problem(T=293.15, P=1e5, row_tolerance=1e-10,
                      feed=Amounts({"water": .6, "ethanol": 1e-19, "isobutanol": .4}), phases=[
                          Phase("Lw", kind="liquid", composition_guess=(.98, .001, .019)),
                          Phase("Lo", kind="liquid", composition_guess=(.3, .001, .699))])
    result = solve_equilibrium(model, problem)
    assert result.success and result.requested_tolerance_met, result.message
    assert result.amounts[0] > 0. and result.amounts[1] > 0.
    same(evidence, "trace-ethanol/two-liquids/classified-residual", result.convergence.residual, 0.,
         atol=1e-10, rtol=0., unit="1", basis="Identity", method="shared-equilibrium-rows",
         claim="the original trace-ethanol declaration retains two liquids and meets its classified rows")


@pytest.mark.parametrize("total", [1.0, 10.0])
@pytest.mark.parametrize("equilibrium_constant", [1e-4, 1e-8, 1e-12, 1e-16])
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
    # First-order optimality of the extensive objective with the returned constraint multipliers, from the
    # compiled objective gradient and row Jacobian rather than the driver's own reconstruction: a wrong
    # feed scaling of the material multiplier or a nonzero normalization multiplier leaves a residual of
    # order |grad F| (0.4 to 142 here) instead of rounding.
    n = len(central.coordinates)
    gradient = compiled.objective_gradient(list(central.coordinates))
    jacobian = compiled.jacobian(list(central.coordinates))
    constraint_rows = [k for k, item in enumerate(compiled.rows) if str(item.kind).split(".")[-1] not in ("Chemical", "Pressure")]
    assert len(constraint_rows) == len(central.constraint_multipliers)
    kkt = max(abs(gradient[i] + sum(mu * jacobian[row * n + i] for mu, row in zip(central.constraint_multipliers, constraint_rows)))
              for i in range(compiled.primal_size))
    # equivalence standard, CODING_STANDARDS.md, PR #311: scale by the extensive objective gradient.
    objective_scale = max(abs(value) for value in gradient)
    same(evidence, f"{label}/kkt_residual", kkt, 0.0, atol=max(1e-12, 1e-10 * objective_scale), rtol=0.0, unit="1",
         inputs={**inputs, "objective_gradient_max": objective_scale}, basis="Identity",
         method="objective-gradient-plus-jacobian",
         claim="the returned constraint multipliers make the extensive objective stationary on the constraint rows")
    # Stationarity and curvature at the returned point, with the standard-state gauge made explicit so
    # the identity is independent of the driver's internal gauge: g0 = (0, -RT ln K) with the consistent
    # scalar K reproduces the same state. The balance multipliers act through the compiled
    # balance basis, read from the row Jacobian's multiplier columns (-B_ji on a chemical row), so
    # ln(x_i rho) + mu_res_i/RT + g0_i/RT = sum_j lambda_j B_ji for both species from public state quantities.
    explicit = Problem(phases=[Phase("V", kind="vapor")], T=temperature, P=2e5,
                       feed=Amounts({"s0": total / 2, "s1": total / 2}), row_tolerance=1e-10,
                       reactions=[Reaction({"s0": -1.0, "s1": 1.0}, equilibrium_constant=equilibrium_constant)],
                       standard_gibbs_j_per_mol=(0.0, -GAS_CONSTANT_J_PER_MOL_K * temperature * math.log(equilibrium_constant)))
    gauged = solve_equilibrium(model, explicit)
    assert gauged.success and math.isclose(gauged.mole_fractions[1], central.mole_fractions[1], rel_tol=1e-12)
    compiled_explicit = compile_problem(model, explicit)
    n = len(gauged.coordinates)
    jacobian = compiled_explicit.jacobian(list(gauged.coordinates))
    state = model.state(temperature, rho=gauged.molar_densities[0], x=tuple(gauged.mole_fractions))
    for species, standard in ((0, 0.0), (1, -math.log(equilibrium_constant))):
        row = next(k for k, item in enumerate(compiled_explicit.rows) if "Chemical" in str(item.kind) and item.index == species)
        basis = [-jacobian[row * n + compiled_explicit.primal_size + j] for j in range(compiled_explicit.rank)]
        residual_potential = state.residual_chemical_potential_over_rt[species]
        assert residual_potential is not None
        potential = math.log(gauged.mole_fractions[species] * state.molar_density) + residual_potential + standard
        same(evidence, f"{label}/stationarity/{species}", sum(m * b for m, b in zip(gauged.multipliers, basis)), potential,
             atol=1e-12, rtol=1e-12, unit="1", inputs={**inputs, "balance_basis": basis}, basis="Identity",
             method="public-state-chemical-potential",
             claim="the returned balance multipliers reproduce the species chemical potential over RT at the returned point")
    # Curvature certificate at the returned point: positive, reproduced by a fresh solve started there to 1e-9
    # (the pre-correction certificate with optimizer duals differed by 5.9e-8 at K = 1e-16), and on the
    # physical scale of the trace amount N x_B, the log-chart curvature of the trace coordinate. The
    # reproduction row checks consistency only; the scale row falsifies the feed scaling of the multipliers
    # at N = 10 (mu/N or mu N would give 4.8e-7 or -3.7e-6 against 1e-7), while at N = 1 the scalings coincide.
    warm = solve_equilibrium(model, explicit, start=gauged)
    assert gauged.reduced_hessian_minimum > 0.0 and warm.success
    same(evidence, f"{label}/curvature", warm.reduced_hessian_minimum, gauged.reduced_hessian_minimum, atol=0.0, rtol=1e-9,
         unit="1", inputs=inputs, basis="Relationship", method="warm-resolve",
         claim="the reduced-Hessian minimum certified at the returned point is reproduced from that point")
    same(evidence, f"{label}/curvature_scale", gauged.reduced_hessian_minimum, total * gauged.mole_fractions[1], atol=0.0,
         rtol=2.0 * equilibrium_constant + 1e-9, unit="1", inputs=inputs, basis="Relationship", method="trace-amount-scale",
         claim="the minimum reduced curvature of the isomerisation is the trace amount to first order in K")


def test_large_feed_trace_recovers_the_analytic_state(evidence: Any) -> None:
    """At 1e6 mol the former tiny-step state missed x_B by 0.96% and curvature by 1.9%."""
    total, temperature, k = 1e6, 350.0, 1e-8
    model = Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), 0.0)]))
    a = math.log(k) - 1200 / temperature - 0.4 * math.log(temperature / 300) - 0.002 * temperature
    law = ReactionLogPolynomial(a, 1200, 0.4, 0.002, 300, 250, 450, EOS_STANDARD_STATE_ID)
    problem = Problem(phases=[Phase("V", kind="vapor")], T=temperature, P=2e5,
                      feed=Amounts({"s0": total / 2, "s1": total / 2}),
                      reactions=[Reaction({"s0": -1.0, "s1": 1.0}, correlation=law)])
    result = solve_equilibrium(model, problem)
    assert result.success and result.requested_tolerance_met, result.message
    inputs = {"feed_mol": total, "K": k, "T_K": temperature, "P_Pa": 2e5}
    same(evidence, "large-feed-trace/x_B", result.mole_fractions[1], k / (1 + k),
         atol=0.0, rtol=1e-6, unit="1", inputs=inputs, basis="Identity",
         method="logistic-analytic", claim="trace equilibrium is independent of total feed")
    compiled = compile_problem(model, problem)
    point = list(result.coordinates[:compiled.primal_size])
    gradient, jacobian = compiled.objective_gradient(point), compiled.jacobian(point)
    rows = [i for i, row in enumerate(compiled.rows)
            if str(row.kind).split(".")[-1] not in ("Chemical", "Pressure")]
    kkt = max(abs(gradient[j] + sum(mu * jacobian[i * compiled.primal_size + j]
              for mu, i in zip(result.constraint_multipliers, rows)))
              for j in range(compiled.primal_size)) / total
    same(evidence, "large-feed-trace/kkt_per_feed", kkt, 0.0, atol=1e-10, rtol=0.0,
         unit="1", inputs=inputs, basis="Identity", method="objective-gradient-plus-jacobian",
         claim="per-feed stationarity meets the existing intensive criterion")
    assert result.reduced_hessian_minimum > 0.0
    same(evidence, "large-feed-trace/curvature", result.reduced_hessian_minimum,
         total * k / (1 + k), atol=0.0, rtol=2 * k + 1e-9, unit="1", inputs=inputs,
         basis="Relationship", method="trace-amount-scale",
         claim="large-feed reduced curvature has the analytic trace-amount scale")


def test_trace_isomerisation_below_the_composition_cutoff(evidence: Any) -> None:
    """#266: the A <=> B declaration above at K = 1e-24 puts x_B below EQUILIBRIUM_BOX_FLOOR (1e-20). The
    polish rejected every positive trial below that cutoff, so main stalled at raw rows 9.21; Newton reaches
    K/(1+K) to 2.4e-14 relative."""
    model = Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), 0.0)]))
    temperature, k = 350.0, 1e-24
    a = math.log(k) - 1200 / temperature - 0.4 * math.log(temperature / 300) - 0.002 * temperature
    law = ReactionLogPolynomial(a, 1200, 0.4, 0.002, 300, 250, 450, EOS_STANDARD_STATE_ID)
    result = solve_equilibrium(model, Problem(phases=[Phase("V", kind="vapor")], T=temperature, P=2e5,
                                              feed=Amounts({"s0": 0.5, "s1": 0.5}), row_tolerance=1e-10,
                                              reactions=[Reaction({"s0": -1.0, "s1": 1.0}, correlation=law)]))
    assert result.mole_fractions, result.message
    same(evidence, "trace/1e-24/x_B", result.mole_fractions[1], k / (1 + k), atol=0.0, rtol=1e-9, unit="1",
         inputs={"K": k, "feed_mol": 1.0, "T_K": temperature, "P_Pa": 2e5, "message": result.message},
         basis="Identity", method="logistic-analytic", claim="the polish reaches a positive root below the cutoff")


@pytest.mark.parametrize("inventory", [1e-9, 1e-12, 1e-13])
def test_trace_inventory_below_the_composition_box(evidence: Any, inventory: float) -> None:
    """#172: identical ethanes, 1 mol inert s0 and s1 <=> s2 (K(350 K) = 3) fed as `inventory` mol s1, so
    x_s2 = inventory K / ((1 + inventory)(1 + K)) and dx_s2/dT = x_s2 (1 - p) l1 exactly. The uniform 1e-8 box
    demanded 2e-8 mol of s1 + s2, far above the inventory, and the subproblem was infeasible (wheel b1367b32).
    The s1 + s2 moiety row resolves the inventory to rounding; the orthonormal balance basis mixed it into the
    s0 row and returned the 1e-12 mol state 2.8e-4 relative off."""
    model = Mixture(parameters((ETHANE,) * 3, [("k_ij", pair, 0.0) for pair in ((0, 1), (0, 2), (1, 2))]))
    temperature, k = 350.0, 3.0
    a = math.log(k) - 1200 / temperature - 0.4 * math.log(temperature / 300) - 0.002 * temperature
    law = ReactionLogPolynomial(a, 1200, 0.4, 0.002, 300, 250, 450, EOS_STANDARD_STATE_ID)
    problem = Problem(phases=[Phase("V", kind="vapor")], T=temperature, P=2e5, row_tolerance=1e-10,
                      feed=Amounts({"s0": 1.0, "s1": inventory, "s2": 0.0}),
                      reactions=[Reaction({"s1": -1.0, "s2": 1.0}, correlation=law)])
    result = solve_equilibrium(model, problem)
    assert result.success and result.requested_tolerance_met, result.message
    assert result.reduced_hessian_minimum > 0.0
    p = k / (1 + k)
    x_s2 = inventory / (1 + inventory) * p
    inputs = {"inventory_mol": inventory, "K": k, "T_K": temperature, "P_Pa": 2e5,
              "classified_residual": result.convergence.residual, "iterations": result.iterations}
    same(evidence, f"trace-inventory/{inventory:g}/x_s2", result.mole_fractions[2], x_s2, atol=0.0, rtol=1e-12,
         unit="1", inputs=inputs, basis="Identity", method="logistic-analytic",
         claim="a species whose whole inventory lies below the box times the phase amount reaches its analytic state")
    item = solved_state_actions(compile_problem(model, problem), result, SolvedStateActionRequest(
        [SolvedStateObservable(Kind.PhaseMoleFraction, 0, 2)], [SolvedStateActionDirection(1.0, 0.0, [0.0] * 3, [])])).results[0]
    assert item.status == SolvedStateActionStatus.Available, item.status
    l1 = -1200 / temperature**2 + 0.4 / temperature + 0.002
    same(evidence, f"trace-inventory/{inventory:g}/A1", item.action, x_s2 * (1 - p) * l1, atol=0.0, rtol=5e-6,
         unit="K^-1", inputs=inputs, basis="Identity", method="logistic-analytic",
         claim="the temperature action of the trace-inventory species matches the logistic chain")


def test_moiety_basis_of_an_oligomerization_network(evidence: Any) -> None:
    """#172: identical ethanes, 3 s0 <=> s2, s0 + 2 s2 <=> s3, s2 + 2 s3 <=> s4 (K = 2) beside inert s1. The one
    mass moiety sum M_i n_i (M = 1, 3, 7, 17) and s1 must be conserved. Echelon pivots chosen on the SVD null
    space admitted a dependent species by rounding, gave 1e14 basis entries and a rank refusal."""
    model = Mixture(parameters((ETHANE,) * 5, [("k_ij", (i, j), 0.0) for i in range(5) for j in range(i + 1, 5)]))
    result = solve_equilibrium(model, Problem(
        phases=[Phase("V", kind="vapor")], T=300.0, P=1e5, feed=Amounts({"s0": 1.0, "s1": 1.0, "s2": 0.0, "s3": 0.0, "s4": 0.0}),
        reactions=[Reaction({"s0": -3.0, "s2": 1.0}, equilibrium_constant=2.0),
                   Reaction({"s0": -1.0, "s2": -2.0, "s3": 1.0}, equilibrium_constant=2.0),
                   Reaction({"s2": -1.0, "s3": -2.0, "s4": 1.0}, equilibrium_constant=2.0)]))
    assert result.success and result.requested_tolerance_met, result.message
    n = [result.amounts[0] * x for x in result.mole_fractions]
    same(evidence, "oligomerization/mass-moiety", sum(m * v for m, v in zip((1, 0, 3, 7, 17), n)), 1.0, atol=0.0, rtol=1e-12,
         unit="mol", inputs={"K": 2.0, "T_K": 300.0, "P_Pa": 1e5}, basis="Identity", method="moiety-balance",
         claim="the echelon balance basis conserves the network's mass moiety")


@pytest.mark.parametrize("total", [1.0, 10.0])
def test_ionic_constraint_multipliers_are_extensive(evidence: Any, total: float) -> None:
    """Khudaida Figure 2 point 1: nonzero electrochemical offset, two stable liquid phases.
    The independent KKT and public chemical potentials falsify charge-dual sign and amount
    scaling. The 10 mol case distinguishes phase amount from phase fraction. Existing packet
    physics is consumed unchanged; this fixed state does not qualify the full LLE curve."""
    names = ("water", "ethanol", "isobutanol", "sodium-cation", "chloride-anion")
    model = Mixture(Parameters.from_json(
        parameters_json("khudaida-2026-figure-2-electrolyte-lle", 7), components=names))
    feed = [0.7295582904360919, 0.013336215598851255, 0.20495308450174893,
            0.02607620473165392, 0.02607620473165392]
    starts = [(0.6416650425482523, [0.9166627632724768, 0.0010121592710021154,
               0.001360946512171668, 0.04048206547217478, 0.04048206547217478]),
              (0.35833495745174776, [0.39451311265799627, 0.03540471872289053,
               0.5695224773809328, 0.0002798456190902245, 0.0002798456190902245])]
    problem = Problem(T=293.15, P=1e5, feed=Amounts(dict(zip(names, [total * n for n in feed]))),
                      phases=[Phase(f"L{p}", kind="liquid", amount=Free(), composition_guess=x)
                              for p, (n, x) in enumerate(starts)], row_tolerance=1e-10)
    compiled = compile_problem(model, problem)
    result = solve_equilibrium(model, problem)
    assert result.success and result.requested_tolerance_met, result.message
    assert result.reduced_hessian_minimum > 0.0
    point = list(result.coordinates[:compiled.primal_size])
    gradient, jacobian = compiled.objective_gradient(point), compiled.jacobian(point)
    rows = [(i, row) for i, row in enumerate(compiled.rows)
            if str(row.kind).split(".")[-1] not in ("Chemical", "Pressure")]
    electro = next(k for k, (_, row) in enumerate(rows) if str(row.kind).split(".")[-1] == "Electro")
    phase = rows[electro][1].phase
    inputs = {"feed_mol": total, "T_K": 293.15, "P_Pa": 1e5,
              "classified_residual": result.convergence.residual,
              "electrochemical_offset": result.multipliers[compiled.rank]}
    label = f"ionic-multiplier/{total:g}"

    def stationarity(multipliers: list[float]) -> float:
        return max(abs(gradient[j] + sum(mu * jacobian[i * compiled.primal_size + j]
                   for mu, (i, _) in zip(multipliers, rows))) for j in range(compiled.primal_size))

    same(evidence, f"{label}/kkt", stationarity(list(result.constraint_multipliers)), 0.0,
         atol=1e-12, rtol=0.0, unit="1", inputs=inputs, basis="Identity",
         method="objective-gradient-plus-jacobian", claim="ionic constraint duals make the extensive objective stationary")
    wrong = list(result.constraint_multipliers)
    wrong[electro] = -(result.amounts[phase] / total) * result.multipliers[compiled.rank]
    if total == 10.0:
        assert stationarity(wrong) > 1e-6  # phase fractions cannot replace extensive phase amounts
    potentials = []
    for p in range(2):
        x = tuple(result.mole_fractions[5 * p:5 * (p + 1)])
        state = model.state(293.15, rho=result.molar_densities[p], x=x)
        phase_potentials = []
        for i in range(5):
            residual_potential = state.residual_chemical_potential_over_rt[i]
            assert residual_potential is not None
            phase_potentials.append(math.log(x[i] * state.molar_density) + residual_potential)
        potentials.append(phase_potentials)
        same(evidence, f"{label}/charge/{p}", x[3] - x[4], 0.0, atol=1e-12, rtol=0.0,
             unit="1", inputs=inputs, basis="Identity", method="public-composition",
             claim="each declared liquid phase is electroneutral")
    delta = [potentials[phase][i] - potentials[1 - phase][i] for i in range(5)]
    same(evidence, f"{label}/salt-transfer", delta[3] + delta[4], 0.0, atol=1e-12, rtol=0.0,
         unit="1", inputs=inputs, basis="Identity", method="public-state-chemical-potential",
         claim="neutral salt chemical potential is equal across the two liquids")
    same(evidence, f"{label}/electrochemical-offset", (delta[3] - delta[4]) / 2,
         result.multipliers[compiled.rank], atol=1e-12, rtol=0.0, unit="1", inputs=inputs,
         basis="Identity", method="public-state-chemical-potential",
         claim="the row electrochemical offset matches the charged species chemical potentials")


_MEA = ("carbon-dioxide", "monoethanolamine", "water", "protonated-monoethanolamine", "carbamate-anion",
        "bicarbonate-anion", "carbonate-anion", "hydronium-cation", "hydroxide-anion")
_MEA_REACTIONS = (  # R1-R5: water, CO2 and bicarbonate dissociation, carbamate hydrolysis, MEAH+ dissociation
    {"water": -2.0, "hydronium-cation": 1.0, "hydroxide-anion": 1.0},
    {"carbon-dioxide": -1.0, "water": -2.0, "bicarbonate-anion": 1.0, "hydronium-cation": 1.0},
    {"water": -1.0, "bicarbonate-anion": -1.0, "carbonate-anion": 1.0, "hydronium-cation": 1.0},
    {"monoethanolamine": 1.0, "water": -1.0, "carbamate-anion": -1.0, "bicarbonate-anion": 1.0},
    {"protonated-monoethanolamine": -1.0, "water": -1.0, "monoethanolamine": 1.0, "hydronium-cation": 1.0})
_MEA_FEED = [1e-05, 1e-05, 8e-05, 1e-05, 0.00011, 1e-05]  # ion amounts [mol] after CO2, MEA, water


def _mea_mapping(path: Path) -> dict[str, Any]:
    """An MEA document with its historical CO2/MEA omission declared explicitly (#569): no cross edge, not a
    claim that the physical interaction is absent."""
    mapping: dict[str, Any] = json.loads(path.read_text())
    mapping["topology"].setdefault("no_cross_association", []).append({
        "component_id_a": "carbon-dioxide", "component_id_b": "monoethanolamine",
        "provenance": {"source_id": mapping["document_id"], "locator": f"{path.relative_to(Path(__file__).resolve().parents[3])}: no CO2/MEA edge"}})
    return mapping


@pytest.mark.parametrize("label, temperature, carbon_dioxide, guess, speciation", [
    ("bottinger-046", 313.15, 0.34290000000000004,
     [0.00012044382032884565, 0.03820530851007986, 0.88492326364512, 0.03838499164509782, 0.03561505266280674,
      0.0027319402457567526, 1.899928767640581e-05, 1.0976016899728266e-11, 1.7215753829197593e-10],
     [1.1167290436645361e-06, 0.037675712115251, 0.8817275121233354, 0.04210421915583975, 0.032438812298105475,
      0.0024383205817141316, 0.0036127794353711146, 7.803077515596226e-11, 1.5274833087670273e-06]),
    ("matin-016", 293.15, 0.4549,
     [0.00033720310764111846, 0.01727031113072209, 0.8809380252760517, 0.05074929564786032, 0.04416142155735913,
      0.006499612403774798, 4.4130824078475495e-05, 6.971446826142033e-12, 4.5540883406069727e-11],
     [2.0195659951573153e-06, 0.015800147213209173, 0.8689944151250493, 0.06414532316923494, 0.03227317187141063,
      0.00569706237942454, 0.013087228408673619, 8.297500114780182e-11, 6.321840276008174e-07]),
])
def test_reacting_liquid_cold_start_recovers_manufactured_speciation(
    evidence: Any, label: str, temperature: float, carbon_dioxide: float, guess: list[float], speciation: list[float]
) -> None:
    """#154/#157: a cold one-liquid MEA solve from the application's unspeciated guess reaches the speciation.
    The speciation is exact by construction: each density-basis ln K_r = sum_i nu_ri (ln(x_i rho) + mu_i^res/RT)
    at x (packet 3 EOS with the applications' Born runtime defaults), and x carries the feed's balances. The
    684b213a start, which seeded the guess directly, exhausted the 200-iteration budget on both states."""
    mapping = _mea_mapping(parameters_json("mea-co2-h2o-nine-species-estimation", 4))
    for family in mapping["model_families"]:
        if family["kind"] == "electrolyte":
            family.update(c_shell=1.0, c_dielectric=1.0)
    model = Mixture(Parameters.from_mapping(mapping, components=_MEA))
    state = model.state(temperature, P=101325.0, x=speciation, phase="liquid")
    ln_k = []
    for nu in _MEA_REACTIONS:
        terms = []
        for species, x, mu in zip(_MEA, speciation, state.residual_chemical_potential_over_rt):
            coefficient = nu.get(species, 0.0)
            if coefficient == 0:
                continue
            assert mu is not None
            terms.append(coefficient * (math.log(x * state.molar_density) + mu))
        ln_k.append(sum(terms))
    reactions = [Reaction(nu, correlation=ReactionLogPolynomial(a, 0.0, 0.0, 0.0, temperature, 250.0, 450.0,
                                                                EOS_STANDARD_STATE_ID)) for nu, a in zip(_MEA_REACTIONS, ln_k)]
    problem = Problem(phases=[Phase("L", kind="liquid", composition_guess=guess)], T=temperature, P=101325.0,
                      feed=Amounts(dict(zip(_MEA, [carbon_dioxide, 0.99998, 7.910947365207347, *_MEA_FEED]))),
                      reactions=reactions)
    result = solve_equilibrium(model, problem)
    assert result.success, result.message
    for species, value, reference in zip(_MEA, result.mole_fractions, speciation):
        same(evidence, f"cold-start/{label}/{species}", value, reference, atol=0.0, rtol=1e-8, unit="1",
             inputs={"T_K": temperature, "P_Pa": 101325.0, "iterations": result.iterations}, basis="Identity",
             method="manufactured-speciation", claim="a cold reacting-liquid solve recovers the manufactured speciation")


@pytest.mark.parametrize("pressure, loading", [(105159.47802968061, 0.25557798094241313),
                                               (109766.9773907787, 0.3259285401314349)])
def test_trace_curvature_below_rounding_does_not_decide_stability(pressure: float, loading: float) -> None:
    """#159: packet-3 MEA liquid at 318.15 K with the absorber's fixed density-basis ln K and CO2 near 1e-15.
    The raw reduced-Hessian minimum (-2.9e-13 and -1e-16 to -3e-16) lies under its rounding n eps ||S|| = 8e-11, and
    wheels 684b213a, 4105e2a4 and main 99e333cd all rejected both states on its sign. The equilibrated inertia
    (minimum 1.4e-3 against a band below 1e-13, analyses/2026-stability-noise-band) certifies them on every solve."""
    model = Mixture(Parameters.from_mapping(_mea_mapping(parameters_json("mea-co2-h2o-nine-species-estimation", 3)), components=_MEA))
    ln_k = (-605.348, -587.8, -948.2312, -25.1014, -298.6415)
    reactions = [Reaction(nu, correlation=ReactionLogPolynomial(a, 0.0, 0.0, 0.0, 318.15, 250.0, 450.0,
                                                                EOS_STANDARD_STATE_ID)) for nu, a in zip(_MEA_REACTIONS, ln_k)]
    problem = Problem(phases=[Phase("L", kind="liquid")], T=318.15, P=pressure, reactions=reactions,
                      feed=Amounts({"carbon-dioxide": loading * 9.7292821, "monoethanolamine": 9.7292821,
                                    "water": 76.9751925}))
    first, second = solve_equilibrium(model, problem), solve_equilibrium(model, problem)
    assert first.success and second.success, first.message
    assert first.mole_fractions == second.mole_fractions and first.mole_fractions[0] < 1e-12


def test_mea_liquid_reaches_carbon_dioxide_root_below_the_cutoff() -> None:
    """#266: packet-3 MEA liquid at loading 1e-6, 318.15 K and 8877.7929 Pa has neutral CO2 at x = 1.1e-21, below
    EQUILIBRIUM_BOX_FLOOR. Every full polish step contracted the raw rows but was rejected by that cutoff; the
    half steps reduced Newton to bisection in ln x and main stalled at raw rows 27.68. #322: the action there was
    refused by that cutoff; the carbon amounts move with the CO2 feed, sum_C dn_i/d ln n_CO2 = n_CO2."""
    model = Mixture(Parameters.from_mapping(_mea_mapping(parameters_json("mea-co2-h2o-nine-species-estimation", 3)), components=_MEA))
    reactions = [Reaction(nu, correlation=ReactionLogPolynomial(a, 0.0, 0.0, 0.0, 318.15, 250.0, 450.0,
                 EOS_STANDARD_STATE_ID)) for nu, a in zip(_MEA_REACTIONS, (-605.348, -587.8, -948.2312, -25.1014, -298.6415))]
    carbon_dioxide = 1e-6 * 9.7292821
    problem = Problem(phases=[Phase("L", kind="liquid")], T=318.15, P=8877.7929, row_tolerance=1e-10, reactions=reactions,
                      feed=Amounts({"carbon-dioxide": carbon_dioxide, "monoethanolamine": 9.7292821, "water": 76.9751925}))
    result = solve_equilibrium(model, problem)
    assert result.success and result.requested_tolerance_met, result.message
    assert max(abs(value) for value in result.residuals) <= 1e-10 and 0.0 < result.mole_fractions[0] < 1e-20
    items = solved_state_actions(compile_problem(model, problem), result, SolvedStateActionRequest(
        [SolvedStateObservable(Kind.PhaseComponentAmount, 0, i) for i in (0, 4, 5, 6)],
        [SolvedStateActionDirection(0.0, 0.0, [carbon_dioxide] + [0.0] * 8, [])])).results
    carbon = [item.action for item in items if item.status == SolvedStateActionStatus.Available and item.action is not None]
    assert len(carbon) == 4 and math.isclose(sum(carbon), carbon_dioxide, rel_tol=1e-10, abs_tol=0.0), items[0].status


@pytest.mark.parametrize("k_ij, temperature", [(0.2, 200.0), (0.5, 150.0)])
def test_symmetric_isomer_curvature_maximum_is_rejected(k_ij: float, temperature: float) -> None:
    """#159: identical-ethane A <=> B (K = 1) started at the symmetric stationary point x = 1/2, where a large
    k_ij makes the reaction direction a Gibbs maximum. Its equilibrated minimum (-0.052 and -0.65) lies
    ten orders outside the 5e-15 rounding band: rejected as not positive definite, never indeterminate.
    Rounding may carry the stronger maximum (0.5, 150 K) to one of its two symmetric minima near a pure
    isomer (x_B = 6.2e-6 on the #172 moiety basis); a success must be that certified minimum, never x = 1/2."""
    model = Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), k_ij)]))
    result = solve_equilibrium(model, Problem(phases=[Phase("L", kind="liquid")], T=temperature, P=5e6,
                                              feed=Amounts({"s0": 0.5, "s1": 0.5}),
                                              reactions=[Reaction({"s0": -1.0, "s1": 1.0}, equilibrium_constant=1.0)]))
    if result.success:
        assert k_ij == 0.5 and min(result.mole_fractions) < 1e-4 and result.reduced_hessian_minimum > 0.0
        return
    assert "is not positive definite" in result.message
    assert str(result.convergence.stop).split(".")[-1] == "converged" and result.reduced_hessian_minimum < -0.1


def test_zero_feed_ions_allow_provisional_flash_start() -> None:
    """#265: packet-3 neutral feed at loading 1e-3 reaches an accepted flash that conserves C, H, N, O and charge."""
    model = Mixture(Parameters.from_mapping(_mea_mapping(parameters_json("mea-co2-h2o-nine-species-estimation", 3)), components=_MEA))
    reactions = [Reaction(nu, correlation=ReactionLogPolynomial(a, 0., 0., 0., 318.15, 250., 450.,
                 EOS_STANDARD_STATE_ID)) for nu, a in zip(_MEA_REACTIONS, (-605.348, -587.8, -948.2312, -25.1014, -298.6415))]
    feed = {"carbon-dioxide": 1e-3 * 9.7292821, "monoethanolamine": 9.7292821, "water": 76.9751925}
    problem = Problem(T=318.15, P=5918.5286, row_tolerance=1e-10, reactions=reactions, feed=Amounts(feed),
                      phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", support=_MEA[:3])])
    result = solve_equilibrium(model, problem)
    assert result.success, result.message
    count = len(_MEA)
    returned = [math.fsum(result.amounts[p] * result.mole_fractions[p * count + i] for p in range(result.phase_count))
                for i in range(count)]
    fed = [feed.get(species, 0.0) for species in _MEA]
    # Species formulas in _MEA order: C, H, N, O and charge number.
    for row in ((1, 2, 0, 2, 3, 1, 1, 0, 0), (0, 7, 2, 8, 6, 1, 0, 3, 1), (0, 1, 0, 1, 1, 0, 0, 0, 0),
                (2, 1, 1, 1, 3, 3, 3, 1, 1), (0, 0, 0, 1, -1, -1, -2, 1, -1)):
        balance = math.fsum(a * (n - f) for a, n, f in zip(row, returned, fed))
        assert abs(balance) <= 1e-12 * sum(feed.values()), (row, balance)


def test_nominally_unloaded_mea_resolves_trace_ion_stationarity(evidence: Any) -> None:
    """#261: 7 mol MEA/kg water at CO2 loading 1e-5, 353.15 K and 101325 Pa (engine/data/reference/
    mea-nominally-unloaded, Amine-Thermodynamics #124 record). CO2 (4.9e-13) and H3O+ (1.6e-12) lie far below
    the optimizer box, so the polish starts at raw rows near 10. The former three-step cap stopped a still
    quadratic contraction with every ionic stationarity row near -1.8e-9 z_i^2 under a 7.6e-12 classified
    residual; continuing while raw rows exceed the criterion reaches rounding at the fourth step. The equilibrium heat-capacity action at the polished state
    matches the centered enthalpy difference of re-solved neighbours."""
    root = Path(__file__).resolve().parents[2] / "data" / "reference"
    request = json.loads((root / "mea-nominally-unloaded" / "request.json").read_text())
    system = request["reaction_system"]
    mapping = _mea_mapping(root / "mea-nominally-unloaded" / "parameters.json")
    for family in mapping["model_families"]:
        if family["kind"] == "electrolyte":
            family.update(c_shell=1.0, c_dielectric=1.0)
    calorics = json.loads((root / "mea-zero-inventory" / "ideal-calorics.json").read_text())
    p0, components = calorics["reference_pressure_pa"], calorics["components"]
    intervals: list[list[IdealInterval]] = [[] for _ in _MEA]
    for name in ("carbon-dioxide", "water"):
        lower, upper = components[name]["range_k"]
        intervals[_MEA.index(name)] = [IdealInterval(lower, upper, True, True, IdealCorrelation(
            IdealShomate(tuple(components[name]["coefficients"]), components[name]["formation_enthalpy_j_per_mol"]), p0))]
    mea = components["monoethanolamine"]  # Cp(T) re-centred on the h = s = 0 anchor
    centre = mea["enthalpy_entropy_anchor_k"]
    shifted = [sum(a * math.comb(j, k) * centre ** (j - k) for j, a in enumerate(mea["coefficients"]) if j >= k)
               for k in range(4)]
    lower, upper = mea["range_k"]
    intervals[1] = [IdealInterval(lower, upper, True, True, IdealCorrelation(IdealPolynomial(shifted, centre, 0.0, 0.0), p0))]
    model = Mixture(Parameters.from_mapping(mapping), thermochemistry=ThermochemistryRecord(p0, intervals))
    reactions = []
    for record in system["engine_reactions"]:
        law, reference = record["engine_correlation"], record["engine_reference"]
        reactions.append(Reaction(dict(zip(_MEA, record["stoichiometry"])), name=record["reaction_id"],
            correlation=ReactionLogPolynomial(*(law[key] for key in (
                "a", "b", "c", "d", "reference_temperature", "temperature_min", "temperature_max", "standard_state_id"))),
            reference=ReactionReference(getattr(ReferenceSourceBasis, reference["source_basis"]),
                                        reference_pressure_pa=reference["reference_pressure_pa"])))
    # The application's neutral reference: neutral solutes 1, source ion pairs summed, scaled by 1/8.
    ratios = dict(zip(_MEA[:2] + _MEA[3:], (0.125, 0.125, 0.625, 0.25, 0.125, 0.125, 0.125, 0.125)))
    guess = request["continuation"]["state"]["phases"][0]["mole_fractions"]

    def solve(temperature: float) -> tuple[Problem, Any]:
        problem = Problem(phases=[Phase("L", kind="liquid", support=_MEA, composition_guess=guess)], T=temperature,
                          P=101325.0, feed=Amounts(dict(zip(_MEA, system["feed_amounts_mol"]))), reactions=reactions,
                          neutral_reference=NeutralReference(solvent="water", solute_molality_ratios=ratios,
                                                             standard_molality_mol_per_kg=1.0), topology_declared=True)
        result = solve_equilibrium(model, problem)
        assert result.success, result.message
        return problem, result

    def actions(problem: Problem, result: Any, d_t: float) -> list[tuple[Any, Any]]:
        observables = [SolvedStateObservable(Kind.PhaseProperty, 0, 0, PropertyObservable.TotalEnthalpy),
                       SolvedStateObservable(Kind.PhaseAmount, 0, 0)]
        (batch,) = solved_state_first_actions(compile_problem(model, problem), result, observables,
                                              [SolvedStateActionDirection(d_t, 0.0, [0.0] * len(_MEA), [])])
        assert all(item.status == SolvedStateActionStatus.Available for item in batch.results)
        return [(item.value, item.action) for item in batch.results]

    problem, result = solve(353.15)
    raw = max(abs(value) for value in result.residuals)
    inputs = {"T_K": 353.15, "P_Pa": 101325.0, "classified_residual": result.convergence.residual,
              "unscaled_residual_max": raw, "residuals": list(result.residuals)}
    same(evidence, "mea-unloaded/unscaled-rows", raw, 0.0, atol=1e-9, rtol=0.0, unit="1", inputs=inputs,
         basis="Identity", method="newton-polish", claim="trace ionic stationarity meets the unscaled 1e-9 gate")
    mass = math.fsum(n * m for n, m in zip(system["feed_amounts_mol"], system["molar_masses_kg_per_mol"]))
    (h, dh), (n, dn) = actions(problem, result, 1.0)
    heat_capacity = (n * dh + h * dn) / mass / 1000.0
    enthalpy = []
    for temperature in (353.10, 353.20):  # the enthalpy at re-solved neighbours verifies the exact action
        neighbour, solved = solve(temperature)
        (h, _), (n, _) = actions(neighbour, solved, 0.0)
        enthalpy.append(h * n)
    same(evidence, "mea-unloaded/equilibrium-cp", heat_capacity, (enthalpy[1] - enthalpy[0]) / 0.1 / mass / 1000.0,
         atol=0.0, rtol=1e-6, unit="kJ/(kg K)", inputs=inputs, basis="Identity", method="centered-enthalpy-difference",
         claim="the equilibrium heat-capacity action equals the re-solved enthalpy slope")
    # #563: at 298.15 K the trace-ion phase of this neutral feed carries sum z x = -3.7e-16 (3.6e-12 of
    # sum |z| x); the absolute charge rule admits it back into the public state.
    _, cold = solve(298.15)
    for p in range(cold.phase_count):
        model.state(298.15, rho=cold.molar_densities[p], x=cold.phase_mole_fractions(p))

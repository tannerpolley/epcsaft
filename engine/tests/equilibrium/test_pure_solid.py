"""Halite in water at 298.15 K and 101325 Pa: the falsifier F1-F6 of the pure stoichiometric solid (#573, ADR 0013).

Held 2014 strategy-2 water/Na+/Cl-; Clarke and Glew 1985 (doi:10.1063/1.555730) Table 34 dissolution constants and
Table 17A solubility 6.1476 mol/kg; halite molar volume 27.1 cm3/mol from PHREEQC pitzer.dat (phreeqpython 1.6.2).
"""
from functools import cache
import dataclasses
import json
import math
from typing import Any

import pytest

from epcsaft import ActiveParameter, GAS_CONSTANT_J_PER_MOL_K as R, Mixture, ParameterFamily, Parameters, _core
from epcsaft import regression as reg
from epcsaft.equilibrium import (
    Free, Molalities, NeutralReference, Phase, Pinned, Problem, ReactionLogPolynomial, ReactionReference,
    ReferenceSourceBasis, Solid, SolvedStateActionDirection, SolvedStateActionRequest, SolvedStateObservable,
    SolvedStateObservableKind as Kind, compile_problem, solid_activities, solve_equilibrium, solved_state_actions,
    solved_state_first_actions,
)
from tests.native_inputs import same
from tests.packets import parameters_json

T, P = 298.15, 101325.0
SPECIES = ("water", "sodium-cation", "chloride-anion")
# Clarke and Glew 1985 Table 34, 25 C: dissolution Gibbs energy, enthalpy (J/mol) and heat capacity (J/(mol K)).
DG, DH, DCP = -9038.7, 3858.4, -134.03
LN_K = -DG / (R * T)
C_LAW = DCP / R
B_LAW = -(DH - T * DCP) / R
LAW = ReactionLogPolynomial(LN_K - B_LAW / T, B_LAW, C_LAW, 0.0, T, 288.15, 308.15, "molality-infinite-dilution")
V_SOLID = 27.1e-6
MEASURED = 6.1476  # mol/kg, Table 17A


def held_parameters(shift: tuple[str, float] = ("", 0.0)) -> Parameters:
    """The Held 2014 set with one document coefficient shifted by an absolute amount."""
    document = json.loads(parameters_json("held-2014-strategy2", 1).read_text())
    for record in (*document["components"], *document["pairs"]):
        for item in record["coefficients"]:
            item["value"]["magnitude"] += shift[1] if item["identity"] == shift[0] else 0.0
    return Parameters.from_mapping(document, components=list(SPECIES))


def held(shift: tuple[str, float] = ("", 0.0), active: tuple[ActiveParameter, ...] = ()) -> Mixture:
    return Mixture(held_parameters(shift), active_parameters=list(active))


@cache
def model() -> Mixture:
    return held()


def law(a: float = LAW.a) -> ReactionLogPolynomial:
    return ReactionLogPolynomial(a, LAW.b, LAW.c, LAW.d, LAW.reference_temperature, 288.15, 308.15, "molality-infinite-dilution")


def problem(molality: float, amount: Free | Pinned = Free(), temperature: float = T, pressure: float = P,
            correlation: ReactionLogPolynomial = LAW) -> Problem:
    halite = Solid("halite", formula={"sodium-cation": 1.0, "chloride-anion": 1.0}, correlation=correlation,
                   reference=ReactionReference(ReferenceSourceBasis.CommonMolalityInfiniteDilution, P),
                   molar_volume_m3_per_mol=V_SOLID, amount=amount)
    return Problem(T=temperature, P=pressure, phases=[Phase("L", kind="liquid")], solids=[halite],
                   feed=Molalities({"sodium-cation": molality, "chloride-anion": molality}, solvent=("water",),
                                   solvent_mass_kg=1.0),
                   neutral_reference=NeutralReference(solvent="water", standard_molality_mol_per_kg=1.0,
                                                      solute_molality_ratios={"sodium-cation": 1.0, "chloride-anion": 1.0}))


def solve(molality: float, amount: Free | Pinned = Free(), temperature: float = T, pressure: float = P) -> _core.EquilibriumResult:
    result = solve_equilibrium(model(), problem(molality, amount, temperature, pressure))
    assert result.success, result.message
    assert max(abs(r) for r in result.classified_residuals) <= 1e-10
    return result


def liquid_molality(result: _core.EquilibriumResult) -> float:
    return result.amounts[0] * result.mole_fractions[1]  # mol Na+ per the feed's 1 kg of water


def reference_log_ip(molality: float) -> float:
    """2 ln(m gamma_pm) - ln K from the independent Reference.activity owner; both use one infinite-dilution limit."""
    reference = model().reference(T, P, solvent_fractions=[1.0, 0.0, 0.0], ions=SPECIES[1:], stoichiometry=(1, 1))
    return 2.0 * math.log(molality * reference.activity(molality).mean_ionic_molality_scale) - LN_K


def test_halite_solubility_and_complementarity(evidence: Any) -> None:
    """F1 physical check of the Held 2014 set, F2 activity identity and F3 complementarity."""
    saturated = solve(6.5)
    m_s = liquid_molality(saturated)
    same(evidence, "halite/F1/solubility", m_s, MEASURED, rtol=0.030, atol=0.0)
    same(evidence, "halite/F2/ln-ion-product", reference_log_ip(m_s), 0.0, rtol=0.0, atol=1e-9)
    assert saturated.solid_affinity_decades[0] == pytest.approx(0.0, abs=1e-10)
    assert solid_activities(saturated) == [1.0]
    excess = solve(7.5)
    same(evidence, "halite/F3/feed-invariant-molality", liquid_molality(excess), m_s, rtol=1e-8, atol=0.0)
    for feed, result in ((6.5, saturated), (7.5, excess)):
        same(evidence, f"halite/F3/excess/{feed}", result.solid_amounts[0] + liquid_molality(result), feed,
             rtol=0.0, atol=1e-10 * feed)
    dilute = solve(5.0)
    assert dilute.solid_amounts == [0.0]
    assert solid_activities(dilute) == [None]
    assert dilute.solid_affinity_decades[0] < 0.0


@pytest.mark.parametrize(("molality", "phreeqc"), [(5.0, -0.29986498), (6.0, -0.03325842)])
def test_pinned_halite_saturation_index(evidence: Any, molality: float, phreeqc: float) -> None:
    """F4a exact Pinned(0) identity; F4b Held ePC-SAFT against PHREEQC pitzer.dat (model against model)."""
    result = solve(molality, Pinned(0.0))
    index = result.solid_affinity_decades[0]
    same(evidence, f"halite/F4a/{molality}", index, reference_log_ip(molality) / math.log(10.0), rtol=0.0, atol=1e-9)
    same(evidence, f"halite/F4b/{molality}", index, phreeqc, rtol=0.0, atol=0.04)


@pytest.mark.parametrize("molality", [6.5, 7.5])
def test_halite_polymorph(evidence: Any, molality: float) -> None:
    """R4: a free same-formula solid with ln K + 0.01 stays absent at SI = -0.01/ln 10; halite is unchanged."""
    halite = problem(molality).solids[0]
    shifted = ReactionLogPolynomial(LAW.a + 0.01, LAW.b, LAW.c, LAW.d, LAW.reference_temperature, 288.15, 308.15,
                                    "molality-infinite-dilution")
    copy = dataclasses.replace(halite, name="halite-b", correlation=shifted)
    result = solve_equilibrium(model(), dataclasses.replace(problem(molality), solids=[halite, copy]))
    assert result.success, result.message
    single = solve(molality)
    same(evidence, f"halite/R4/{molality}/m_s", liquid_molality(result), liquid_molality(single), rtol=1e-8, atol=0.0)
    same(evidence, f"halite/R4/{molality}/n_s", result.solid_amounts[0], single.solid_amounts[0], rtol=1e-8, atol=0.0)
    assert result.solid_amounts[1] == 0.0
    same(evidence, f"halite/R4/{molality}/SI", result.solid_affinity_decades[1], -0.01 / math.log(10.0),
         rtol=0.0, atol=1e-9)


def actions(result: _core.EquilibriumResult, molality: float, observables: list[_core.SolvedStateObservable],
            directions: list[tuple[float, float, list[float]]], amount: Free | Pinned = Free()) -> list[float]:
    compiled = compile_problem(model(), problem(molality, amount))
    request = SolvedStateActionRequest(observables, [SolvedStateActionDirection(t, p, feed, []) for t, p, feed in directions])
    batch = solved_state_actions(compiled, result, request)
    assert all(item.action is not None for item in batch.results), [item.status for item in batch.results]
    return [item.action for item in batch.results if item.action is not None]


def test_free_solid_amount_label() -> None:
    """F10: the free halite amount has its own column, after the liquid's coordinates and before the multipliers."""
    compiled = compile_problem(model(), problem(6.5))
    labels = compiled.coordinate_labels()
    assert len(labels) == len(compiled.generic_point()) and labels[labels.index("ln eta[L]") + 1] == "n_solid[halite] / mol"


def test_present_branch_actions(evidence: Any) -> None:
    """F5: exact dm_s/dT and water-feed actions of order 1 and 2 against central differences of native solves."""
    center = solve(6.5)
    sodium, solid = SolvedStateObservable(Kind.PhaseComponentAmount, 0, 1), SolvedStateObservable(Kind.SolidAmount, 0)
    dm_dt, = actions(center, 6.5, [sodium], [(1.0, 0.0, [0.0, 0.0, 0.0])])
    h = 0.1
    fd = (liquid_molality(solve(6.5, temperature=T + h)) - liquid_molality(solve(6.5, temperature=T - h))) / (2 * h)
    same(evidence, "halite/F5/dm_dT", dm_dt, fd, rtol=1e-6, atol=0.0)
    water = [1.0, 0.0, 0.0]
    first = actions(center, 6.5, [sodium, solid], [(0.0, 0.0, water)])
    second = actions(center, 6.5, [sodium, solid], [(0.0, 0.0, water), (0.0, 0.0, water)])
    # Adding water w kg dissolves m_s w mol: the liquid molality and the solid's chemical potential stay fixed.
    mw = 1.0 / center.amounts[0] / center.mole_fractions[0]  # kg/mol from the feed's 1 kg of water
    same(evidence, "halite/F5/dn_Na_dn_w", first[0], liquid_molality(center) * mw, rtol=1e-6, atol=0.0)
    same(evidence, "halite/F5/dn_s_dn_w", first[1], -first[0], rtol=1e-12, atol=1e-14)
    for name, value in zip(("Na", "s"), second):
        same(evidence, f"halite/F5/second/{name}", value, 0.0, rtol=0.0, atol=1e-10)


# Document identity, active parameter (k_ij by flat pair indices) and its Held 2014 value.
PARAMETERS = (("component/sodium-cation/dispersion_energy_over_k", ActiveParameter(ParameterFamily.DispersionEnergyOverK, (1,)), 230.0),
              ("component/chloride-anion/segment_diameter", ActiveParameter(ParameterFamily.SigmaConstant, (2,)), 2.756),
              ("pair/chloride-anion/water/k_ij", ActiveParameter(ParameterFamily.Kij, (2, 6)), -0.25))


def test_present_branch_parameter_actions(evidence: Any) -> None:
    """F5: exact dm_s/dtheta, including the solid's dC_s/dtheta, against central differences of rebuilt-model solves."""
    mixture = held(active=tuple(item for _, item, _ in PARAMETERS))
    center = solve_equilibrium(mixture, problem(6.5))
    assert center.success, center.message
    directions = [SolvedStateActionDirection(0.0, 0.0, [0.0] * 3, [float(j == k) for j in range(3)]) for k in range(3)]
    batches = solved_state_first_actions(compile_problem(mixture, problem(6.5)), center,
                                         [SolvedStateObservable(Kind.PhaseComponentAmount, 0, 1)], directions)
    for (identity, _, value), batch in zip(PARAMETERS, batches):
        exact = batch.results[0].action
        assert exact is not None, batch.results[0].status
        for step, rtol in ((1e-3, 1e-5), (1e-4, 1e-7), (1e-5, 1e-7)):
            h = step * max(1.0, abs(value))
            fd = (liquid_molality(solve_equilibrium(held((identity, h)), problem(6.5)))
                  - liquid_molality(solve_equilibrium(held((identity, -h)), problem(6.5)))) / (2 * h)
            same(evidence, f"halite/F5/dm_dtheta/{identity}/{step}", exact, fd, rtol=rtol, atol=0.0)


def test_halite_pressure_reference(evidence: Any) -> None:
    """F6: dg_s/dP = V_s exactly at fixed T; the transfer at P0 adds no ion partial molar volume."""
    center = solve(6.0, Pinned(0.0))
    observables = [SolvedStateObservable(Kind.SolidSaturationIndex, 0), SolvedStateObservable(Kind.PhaseLogActivity, 0, 1),
                   SolvedStateObservable(Kind.PhaseLogActivity, 0, 2)]
    d_index, d_sodium, d_chloride = actions(center, 6.0, observables, [(0.0, 1.0, [0.0, 0.0, 0.0])], Pinned(0.0))
    # ln 10 SI = sum_i nu_i ln a_i - ln K^rho(T, P0) - V_s (P - P0)/RT on a nonreacting liquid.
    same(evidence, "halite/F6/dg_dP", R * T * (d_sodium + d_chloride - math.log(10.0) * d_index), V_SOLID,
         rtol=1e-9, atol=0.0)
    h = 5.0e4
    fd = (solve(6.0, Pinned(0.0), pressure=P + h).solid_affinity_decades[0]
          - solve(6.0, Pinned(0.0), pressure=P - h).solid_affinity_decades[0]) / (2 * h)
    same(evidence, "halite/F6/dSI_dP", d_index, fd, rtol=1e-6, atol=0.0)


def test_solid_refusals() -> None:
    """A solid without its source pressure, or with a free intensive, is refused with a named cause."""
    unset = problem(6.5)
    assert unset.solids
    unset.solids[0].reference = ReactionReference(ReferenceSourceBasis.CommonMolalityInfiniteDilution)
    with pytest.raises(ValueError, match="reference_pressure_pa"):
        solve_equilibrium(model(), unset)
    free = problem(6.5)
    free.P = Free(P)
    with pytest.raises(ValueError, match="solids require fixed T and P"):
        solve_equilibrium(model(), free)


# Solubility as a fit observation: n_Na/n_w = m_s M_w in the liquid, the existing amount ratio with M_w in the value.
MW = 18.01528e-3  # kg/mol, water


def solubility(molality: float, temperature: float = T, correlation: ReactionLogPolynomial = LAW,
               observed: float = 0.0) -> _core.Observation:
    return reg.observation(held_parameters(), "amount_ratio", problem(molality, temperature=temperature, correlation=correlation),
                           observed=observed, form="difference", scale=1e-3, phase="L", species=["sodium-cation"],
                           denominator_phase="L", denominator_species=["water"])


def ratio(result: _core.EquilibriumResult) -> float:
    assert result.success, result.message
    return result.mole_fractions[1] / result.mole_fractions[0]


def test_solubility_observation(evidence: Any) -> None:
    """The saturated amount ratio equals the direct native solve at a 10 mol/kg feed."""
    prediction = reg.predict(model(), solubility(10.0))
    assert prediction.status == _core.ObservationStatus.Available, prediction.message
    assert prediction.value is not None
    same(evidence, "halite/solubility/observation", prediction.value, ratio(solve(10.0)), rtol=1e-12, atol=0.0)
    same(evidence, "halite/solubility/molality", prediction.value / MW, 6.29890460910619, rtol=1e-12, atol=0.0)


def test_solubility_jacobian(evidence: Any) -> None:
    """Exact EOS and solubility-product columns of the observation against central differences of native solves."""
    p = held_parameters()
    targets = (("dispersion_energy_over_k", "sodium-cation"), ("segment_diameter", "chloride-anion"), ("k_ij", ("chloride-anion", "water")))
    coordinates = [reg.coordinate(p, [(family, target, 1.0)], origin=value, scale=1.0)
                   for (family, target), (_, _, value) in zip(targets, PARAMETERS)]
    coordinates.append(reg.reaction_coordinate("halite", "a", origin=LAW.a, scale=1.0))
    evaluation = reg.differentiate(p, coordinates, [solubility(10.0)], optimizer=[0.0] * 4)
    assert evaluation.statuses == [_core.ObservationStatus.Available]

    def check(name: str, exact: float, move: Any, steps: tuple[tuple[float, float], ...]) -> None:
        for h, rtol in steps:
            fd = (ratio(move(h)) - ratio(move(-h))) / (2 * h)
            same(evidence, f"halite/solubility/jacobian/{name}/{h:.3g}", exact / MW, fd / MW, rtol=rtol, atol=0.0)

    for k, (identity, _, value) in enumerate(PARAMETERS):
        check(identity, evaluation.gradient[k], lambda h, identity=identity: solve_equilibrium(held((identity, h)), problem(10.0)),
              tuple((step * max(1.0, abs(value)), rtol) for step, rtol in ((1e-3, 1e-5), (1e-4, 1e-7), (1e-5, 1e-7))))
    check("a", evaluation.gradient[3], lambda h: solve_equilibrium(model(), problem(10.0, correlation=law(LAW.a + h))),
          ((1e-2, 1e-5), (1e-3, 1e-7), (1e-4, 1e-7)))


TEMPERATURES = (288.15, 298.15, 308.15)


def solubility_fit(feed: float, workers: int) -> Any:
    observed = [reg.predict(model(), solubility(10.0, t)).value for t in TEMPERATURES]
    rows = [solubility(feed, t, law(LAW.a + 0.05), y or 0.0) for t, y in zip(TEMPERATURES, observed)]
    return reg.fit(held_parameters(), [reg.reaction_coordinate("halite", "a", origin=LAW.a, scale=0.1)], rows,
                   controls=reg.FitControls(maximum_iterations=50, maximum_elapsed_time_seconds=120.0, workers=workers))


@pytest.mark.parametrize("workers", [1, 2])
def test_solubility_product_fit(evidence: Any, workers: int) -> None:
    """A manufactured three-temperature solubility fit recovers the law's a from a + 0.05, in process and in workers."""
    result = solubility_fit(10.0, workers)
    assert result.usable, result.message
    same(evidence, f"halite/solubility/fit/{workers}", result.physical[0], LAW.a, rtol=0.0, atol=1e-8)


def test_undersaturated_solubility_refused() -> None:
    """At 6.35 mol/kg and a + 0.05 halite is absent: the row is a named topology event, never the feed molality."""
    prediction = reg.predict(model(), solubility(6.35, correlation=law(LAW.a + 0.05)))
    assert prediction.status == _core.ObservationStatus.TopologyEvent and "solid halite is absent" in prediction.message
    result = solubility_fit(6.35, 1)
    assert not result.usable
    assert {(f.status, "halite" in f.message) for f in result.failures} == {(_core.ObservationStatus.TopologyEvent, True)}

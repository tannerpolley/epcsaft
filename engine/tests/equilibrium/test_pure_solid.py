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

from epcsaft import GAS_CONSTANT_J_PER_MOL_K as R, Mixture, Parameters, _core
from epcsaft.equilibrium import (
    Free, Molalities, NeutralReference, Phase, Pinned, Problem, ReactionLogPolynomial, ReactionReference,
    ReferenceSourceBasis, Solid, SolvedStateActionDirection, SolvedStateActionRequest, SolvedStateObservable,
    SolvedStateObservableKind as Kind, compile_problem, solve_equilibrium, solved_state_actions,
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


@cache
def model() -> Mixture:
    document = json.loads(parameters_json("held-2014-strategy2", 1).read_text())
    return Mixture(Parameters.from_mapping(document, components=list(SPECIES)))


def problem(molality: float, amount: Free | Pinned = Free(), temperature: float = T, pressure: float = P) -> Problem:
    halite = Solid("halite", formula={"sodium-cation": 1.0, "chloride-anion": 1.0}, correlation=LAW,
                   reference=ReactionReference(ReferenceSourceBasis.CommonMolalityInfiniteDilution, P),
                   molar_volume_m3_per_mol=V_SOLID, amount=amount)
    return Problem(T=temperature, P=pressure, phases=[Phase("L", kind="liquid")], solids=[halite],
                   feed=Molalities({"sodium-cation": molality, "chloride-anion": molality}, solvent=("water",),
                                   solvent_mass_kg=1.0),
                   neutral_reference=NeutralReference(solvent="water", standard_molality_mol_per_kg=1.0,
                                                      solute_molality_ratios={"sodium-cation": 1.0, "chloride-anion": 1.0}))


def solve(molality: float, amount: Free | Pinned = Free(), **state: float) -> _core.EquilibriumResult:
    result = solve_equilibrium(model(), problem(molality, amount, **state))
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
    same(evidence, "halite/F1/solubility", m_s, MEASURED, rtol=0.030, atol=0.0,
         method="native solid-liquid solve", claim="Held 2014 set within the 0.02-decade target at saturation")
    same(evidence, "halite/F2/ln-ion-product", reference_log_ip(m_s), 0.0, rtol=0.0, atol=1e-9,
         method="Reference.activity at the solved molality", claim="the saturation row imposes 2 ln(m g) = ln K")
    assert saturated.solid_affinity_decades[0] == pytest.approx(0.0, abs=1e-10)
    excess = solve(7.5)
    same(evidence, "halite/F3/feed-invariant-molality", liquid_molality(excess), m_s, rtol=1e-8, atol=0.0,
         method="native solves at 6.5 and 7.5 mol/kg", claim="the saturated liquid is feed invariant")
    for feed, result in ((6.5, saturated), (7.5, excess)):
        same(evidence, f"halite/F3/excess/{feed}", result.solid_amounts[0] + liquid_molality(result), feed,
             rtol=0.0, atol=1e-10 * feed, method="element balance", claim="n_s is the feed excess")
    dilute = solve(5.0)
    assert dilute.solid_amounts == [0.0]
    assert dilute.solid_affinity_decades[0] < 0.0


@pytest.mark.parametrize(("molality", "phreeqc"), [(5.0, -0.29986498), (6.0, -0.03325842)])
def test_pinned_halite_saturation_index(evidence: Any, molality: float, phreeqc: float) -> None:
    """F4a exact Pinned(0) identity; F4b Held ePC-SAFT against PHREEQC pitzer.dat (model against model)."""
    result = solve(molality, Pinned(0.0))
    index = result.solid_affinity_decades[0]
    same(evidence, f"halite/F4a/{molality}", index, reference_log_ip(molality) / math.log(10.0), rtol=0.0, atol=1e-9,
         method="Reference.activity", claim="SI = [2 ln(m g) - ln K]/ln 10")
    same(evidence, f"halite/F4b/{molality}", index, phreeqc, rtol=0.0, atol=0.04,
         method="PHREEQC 3 pitzer.dat via phreeqpython 1.6.2", claim="2 x the 0.02-decade activity target")


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
    same(evidence, f"halite/R4/{molality}/m_s", liquid_molality(result), liquid_molality(single), rtol=1e-8, atol=0.0,
         method="native solve with two free solids", claim="the lower-g_s solid takes the amount")
    same(evidence, f"halite/R4/{molality}/n_s", result.solid_amounts[0], single.solid_amounts[0], rtol=1e-8, atol=0.0,
         method="native solve with two free solids", claim="the lower-g_s solid takes the amount")
    assert result.solid_amounts[1] == 0.0
    same(evidence, f"halite/R4/{molality}/SI", result.solid_affinity_decades[1], -0.01 / math.log(10.0),
         rtol=0.0, atol=1e-9, method="exact ln K shift", claim="the absent polymorph keeps a negative affinity")


def actions(result: _core.EquilibriumResult, molality: float, observables: list[_core.SolvedStateObservable],
            directions: list[tuple[float, float, list[float]]], amount: Free | Pinned = Free()) -> list[float]:
    compiled = compile_problem(model(), problem(molality, amount))
    request = SolvedStateActionRequest(observables, [SolvedStateActionDirection(t, p, feed, []) for t, p, feed in directions])
    batch = solved_state_actions(compiled, result, request)
    assert all(item.action is not None for item in batch.results), [item.status for item in batch.results]
    return [item.action for item in batch.results if item.action is not None]


def test_present_branch_actions(evidence: Any) -> None:
    """F5: exact dm_s/dT and water-feed actions of order 1 and 2 against central differences of native solves."""
    center = solve(6.5)
    sodium, solid = SolvedStateObservable(Kind.PhaseComponentAmount, 0, 1), SolvedStateObservable(Kind.SolidAmount, 0)
    dm_dt, = actions(center, 6.5, [sodium], [(1.0, 0.0, [0.0, 0.0, 0.0])])
    h = 0.1
    fd = (liquid_molality(solve(6.5, temperature=T + h)) - liquid_molality(solve(6.5, temperature=T - h))) / (2 * h)
    same(evidence, "halite/F5/dm_dT", dm_dt, fd, rtol=1e-6, atol=0.0, method="central difference, 0.1 K",
         claim="exact implicit actions on the present branch")
    water = [1.0, 0.0, 0.0]
    first = actions(center, 6.5, [sodium, solid], [(0.0, 0.0, water)])
    second = actions(center, 6.5, [sodium, solid], [(0.0, 0.0, water), (0.0, 0.0, water)])
    # Adding water w kg dissolves m_s w mol: the liquid molality and the solid's chemical potential stay fixed.
    mw = 1.0 / center.amounts[0] / center.mole_fractions[0]  # kg/mol from the feed's 1 kg of water
    same(evidence, "halite/F5/dn_Na_dn_w", first[0], liquid_molality(center) * mw, rtol=1e-6, atol=0.0,
         method="saturated-liquid balance", claim="the liquid dissolves m_s M_w per mol of water")
    same(evidence, "halite/F5/dn_s_dn_w", first[1], -first[0], rtol=1e-12, atol=1e-14,
         method="sodium balance", claim="the solid supplies what the liquid takes")
    for name, value in zip(("Na", "s"), second):
        same(evidence, f"halite/F5/second/{name}", value, 0.0, rtol=0.0, atol=1e-10,
             method="saturated amounts are linear in the water feed", claim="second feed action")


def test_halite_pressure_reference(evidence: Any) -> None:
    """F6: dg_s/dP = V_s exactly at fixed T; the transfer at P0 adds no ion partial molar volume."""
    center = solve(6.0, Pinned(0.0))
    observables = [SolvedStateObservable(Kind.SolidSaturationIndex, 0), SolvedStateObservable(Kind.PhaseLogActivity, 0, 1),
                   SolvedStateObservable(Kind.PhaseLogActivity, 0, 2)]
    d_index, d_sodium, d_chloride = actions(center, 6.0, observables, [(0.0, 1.0, [0.0, 0.0, 0.0])], Pinned(0.0))
    # ln 10 SI = sum_i nu_i ln a_i - ln K^rho(T, P0) - V_s (P - P0)/RT on a nonreacting liquid.
    same(evidence, "halite/F6/dg_dP", R * T * (d_sodium + d_chloride - math.log(10.0) * d_index), V_SOLID,
         rtol=1e-9, atol=0.0, method="exact pressure actions", claim="dg_s/dP = V_s")
    h = 5.0e4
    fd = (solve(6.0, Pinned(0.0), pressure=P + h).solid_affinity_decades[0]
          - solve(6.0, Pinned(0.0), pressure=P - h).solid_affinity_decades[0]) / (2 * h)
    same(evidence, "halite/F6/dSI_dP", d_index, fd, rtol=1e-6, atol=0.0, method="central difference, 0.5 bar",
         claim="the reference transfer stays at P0 when the system pressure moves")


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

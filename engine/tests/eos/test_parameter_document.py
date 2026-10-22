"""#608: parameter documents built from several packets equal the hand-edited documents they replace.
Input assembly only; no physical-validation claim."""
import json
from typing import Any

import pytest

from epcsaft import Mixture, ParameterError, Parameters
from epcsaft import equilibrium as eq
from tests.packets import PACKETS


def native(parameters: Parameters) -> tuple[Any, ...]:
    """Everything build_parameters reads: components, sites, choices and every entry in order."""
    d: Any = parameters._document
    return (tuple(d.component_ids), tuple(d.site_component), tuple(d.site_multiplicity), d.base_choice,
            d.association_choice, d.electrolyte_choice, d.ionic_dispersion_choice, d.permittivity_choice,
            d.polar_choice, tuple((e.identity, e.field, tuple(e.index), e.magnitude, e.unit) for e in d.entries))


def packet(name: str) -> dict[str, Any]:
    result: dict[str, Any] = json.loads((PACKETS / name / "parameters/parameters.json").read_text())
    return result


def test_methanol_hexane_from_two_packets_equals_the_hand_edited_document() -> None:
    alc, alk = packet("feos-gross-sadowski-2002/1"), packet("gross-sadowski-2001/3")
    alc["components"] += [c for c in alk["components"] if c["component_id"] == "hexane"]
    alc["pairs"] = [{"component_id_a": "methanol", "component_id_b": "hexane", "coefficients": [
        {"identity": "pair/methanol/hexane/k_ij", "family": "k_ij", "value": {"magnitude": 0.03, "unit": "dimensionless"}}]}]
    hand = Parameters.from_mapping(alc, components=["methanol", "hexane"])
    built = Parameters.from_packets(PACKETS, [("feos-gross-sadowski-2002/1", ["methanol"]), ("gross-sadowski-2001/3", ["hexane"])],
                                    k_ij={("methanol", "hexane"): 0.03})
    assert native(built) == native(hand)
    feed = eq.MoleFractions({"methanol": 0.5, "hexane": 0.5})
    splits = [eq.solve_equilibrium(Mixture(p), eq.Problem(phases=eq.Discover(), T=250.0, P=1.0e5, feed=feed)) for p in (hand, built)]
    assert splits[0].phase_count == 2
    assert splits[1].mole_fractions == splits[0].mole_fractions  # bitwise


def test_set_k_ij_replaces_the_stored_pair() -> None:
    doc = packet("gross-2001-methane-ethane/3")
    doc["pairs"][0]["coefficients"][0]["value"]["magnitude"] = 0.02
    hand = Parameters.from_mapping(doc, components=["methane", "ethane"])
    built = Parameters.from_packets(PACKETS, [("gross-2001-methane-ethane/3", ["methane", "ethane"])],
                                    k_ij={("methane", "ethane"): 0.02})
    assert native(built) == native(hand)
    assert [(e[0], e[3]) for e in native(built)[-1] if e[1] == "k_ij"] == [("pair/ethane/methane/k_ij", 0.02)]


def test_polar_and_nonpolar_packets_are_refused_by_their_model_family() -> None:
    with pytest.raises(ParameterError, match="polar model .* conflicts"):
        Parameters.from_packets(PACKETS, [("gross-vrabec-2006/1", ["acetone"]), ("gross-sadowski-2001/3", ["hexane"])])


def test_cation_mixing_with_zero_dodecane_preserves_aqueous_state() -> None:
    doc = packet("figiel-2025-chloride-mixture-refit/2")
    names = [c["component_id"] for c in doc["components"]]
    aqueous = Parameters.from_mapping(doc)
    doc["components"].append(next(c for c in packet("gross-sadowski-2001/3")["components"] if c["component_id"] == "dodecane"))
    doc["components"][-1]["coefficients"].append({"identity": "component/dodecane/relative_permittivity",
        "family": "relative_permittivity", "value": {"magnitude": 2.0, "unit": "dimensionless"}})
    doc["components"][-1]["coefficients"].append({"identity": "component/dodecane/solvation_factor",
        "family": "solvation_factor", "value": {"magnitude": 1.0, "unit": "dimensionless"}})
    for name in names:
        doc["pairs"].append({"component_id_a": name, "component_id_b": "dodecane", "coefficients": [
            {"identity": f"pair/{name}/dodecane/k_ij", "family": "k_ij", "value": {"magnitude": 0.0, "unit": "dimensionless"}}]})
        doc["topology"]["no_cross_association"].append({"component_id_a": name, "component_id_b": "dodecane",
            "provenance": {"locator": "#781 explicit inactive cross association"}})
    extended = Parameters.from_mapping(doc)
    amounts = [1 / aqueous.molar_masses[0], 0.2, 0.1, 0.4, 0.1, 1.35, 0.05]
    x = [n / sum(amounts) for n in amounts]
    old = Mixture(aqueous).state(298.15, rho=30000, x=x)
    new = Mixture(extended).state(298.15, rho=30000, x=[*x, 0.0])
    assert old.residual_helmholtz == new.residual_helmholtz
    assert old.cation_mixing == new.cation_mixing != 0

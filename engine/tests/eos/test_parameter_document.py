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

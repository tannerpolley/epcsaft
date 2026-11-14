"""#608: parameter documents built from several packets equal the hand-edited documents they replace.
Input assembly only; no physical-validation claim."""
import json
from typing import Any

import pytest

from epcsaft import Mixture, ParameterError, Parameters, _core
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


def test_entries_read_before_reassignment_stay_valid() -> None:
    """#1100: entries are copies; a list read earlier once dangled after the setter replaced the vector and aborted the process."""
    d: Any = Parameters.from_mapping(packet("gross-sadowski-2001/3"), components=["methane", "ethane"])._document
    earlier = d.entries
    expected = [(e.identity, e.magnitude) for e in earlier]
    d.entries = d.entries
    d.entries = earlier
    assert [(e.identity, e.magnitude) for e in earlier] == expected == [(e.identity, e.magnitude) for e in d.entries]


def test_packet_composition_bounds() -> None:
    """Verify species molality, 2:1 ionic strength and all-ion fraction; no EOS validation claim."""
    from epcsaft import _core

    def mixture(packet_id: str, names: list[str]) -> Mixture:
        return Mixture(Parameters.from_packets(PACKETS, [(packet_id, names)]))

    refit = "figiel-2025-chloride-mixture-refit/3"
    hcl = mixture(refit, ["water", "hydronium-cation", "chloride-anion"])
    # One kilogram of solvent and stoichiometric ion amounts, normalized independently.
    water_moles = 1 / hcl.parameters.molar_masses[0]
    for molality in (8.0, 5.0):
        total = water_moles + 2 * molality
        records = hcl.outside_domains([water_moles / total, molality / total, molality / total])
        if molality == 5:
            assert records == []
        else:
            assert len(records) == 1
            assert records[0] == {"packet": refit, "domain_id": "hydrogen-chloride-298-15-k",
                "quantity": _core.DomainQuantity.SoluteMolality, "solute": "hydronium-cation",
                "limit": _core.DomainLimit.Maximum, "bound": 7.0,
                "value": pytest.approx(8.0, rel=1e-12), "unit": "mole / kilogram"}

    mgcl2 = mixture(refit, ["water", "magnesium-cation", "chloride-anion"])
    water_moles = 1 / mgcl2.parameters.molar_masses[0]
    for molality in (3.4, 3.2):
        total = water_moles + 3 * molality
        records = mgcl2.outside_domains([water_moles / total, molality / total, 2 * molality / total])
        if molality == 3.2:
            assert records == []  # I = 9.6 must not trigger the absent lithium fit's 8.3 maximum.
        else:
            assert len(records) == 1
            assert records[0] == {"packet": refit, "domain_id": "issue-575-chloride-brine-298-15-k",
                "quantity": _core.DomainQuantity.IonicStrength, "limit": _core.DomainLimit.Maximum,
                "bound": 9.9, "value": pytest.approx(10.2, rel=1e-12), "unit": "mole / kilogram"}

    nacl = mixture("figiel-yu-held-2025/4", ["water", "sodium-cation", "chloride-anion"])
    records = nacl.outside_domains([0.6, 0.2, 0.2])
    assert len(records) == 1
    assert records[0]["quantity"] == _core.DomainQuantity.IonMoleFraction
    assert records[0]["limit"] == _core.DomainLimit.Maximum
    assert records[0]["bound"] == 0.38
    assert records[0]["value"] == pytest.approx(0.4, rel=1e-12)
    assert records[0]["unit"] == "dimensionless"

    held = "held-2014-strategy2/2"
    nacl = mixture(held, ["water", "sodium-cation", "chloride-anion"])
    water_moles = 1 / nacl.parameters.molar_masses[0]
    total = water_moles + 2 * 0.0005
    records = nacl.outside_domains([water_moles / total, 0.0005 / total, 0.0005 / total])
    assert len(records) == 1
    assert records[0] == {"packet": held, "domain_id": "held-2014-chlorides-298-15-k-1-bar",
        "quantity": _core.DomainQuantity.SoluteMolality, "solute": "sodium-cation",
        "limit": _core.DomainLimit.Minimum, "bound": 0.001,
        "value": pytest.approx(0.0005, rel=1e-12), "unit": "mole / kilogram"}

    # A's sodium bound must not apply when B supplies sodium, although sodium is in the final model.
    other = "held-2014-strategy2/1"
    # Cross-packet pairs are explicit inputs; retain B's published coefficients unchanged.
    cross_pairs = {(p["component_id_a"], p["component_id_b"]): c["value"]["magnitude"]
                   for p in packet(other)["pairs"] for c in p["coefficients"]
                   if "sodium-cation" in (p["component_id_a"], p["component_id_b"]) and c["family"] == "k_ij"}
    mixed = Mixture(Parameters.from_packets(PACKETS,
        [(held, ["water", "chloride-anion"]), (other, ["sodium-cation"])], k_ij=cross_pairs))
    water_moles = 1 / mixed.parameters.molar_masses[0]
    total = water_moles + 2 * 5.0
    records = mixed.outside_domains([water_moles / total, 5.0 / total, 5.0 / total])
    assert len(records) == 1
    assert records[0] == {"packet": other, "domain_id": "held-2014-chlorides-298-15-k-1-bar",
        "quantity": _core.DomainQuantity.SoluteMolality, "solute": "sodium-cation",
        "limit": _core.DomainLimit.Maximum, "bound": 4.5,
        "value": pytest.approx(5.0, rel=1e-12), "unit": "mole / kilogram"}


def test_cation_pair_coordinate_round_trip_preserves_the_other_record() -> None:
    from copy import deepcopy
    from epcsaft import _core
    from epcsaft import regression as reg
    doc = packet("figiel-2025-chloride-mixture-refit/3")
    inherited = deepcopy(doc["cation_mixing"])
    added = deepcopy(inherited)
    added["cations"] = ["sodium-cation", "hydronium-cation"]
    for c in added["coefficients"]:
        c["identity"] = "H-Na/" + c["family"]
        c["value"]["magnitude"] = 0.0
    doc["cation_mixing"] = [inherited, added]  # reverse of canonical pair order
    start = Parameters.from_mapping(doc)
    assert start._document is not None
    coordinate = reg.coordinate(start, [("cation_pair_mixing", ("sodium-cation", "hydronium-cation"), 1)],
                                origin=0, scale=1, bounds=(-1, 1))
    chart = reg.evaluate_coordinates(start, [coordinate], [-0.11], [])
    written = _core.fitted_document(start._document, chart.parameters)
    # Native read/write/read preserves every untargeted coefficient, including the other cation pair.
    target = next(e for e in written.entries if e.identity == "H-Na/cation_pair_mixing")
    assert target.magnitude == -0.11
    before = [(e.identity, e.magnitude) for e in start._document.entries if not e.identity.startswith("H-Na/")]
    after = [(e.identity, e.magnitude) for e in written.entries if not e.identity.startswith("H-Na/")]
    assert before == after
    assert native_document_entries(_core.fitted_document(written, _core.build_parameters(written))) == native_document_entries(written)
    duplicated = _core.fitted_document(start._document, chart.parameters)
    entry = next(e for e in duplicated.entries if e.identity == "H-Na/cation_pair_mixing")
    entry.identity += "/duplicate"
    duplicated.entries = [*duplicated.entries, entry]
    with pytest.raises(ParameterError, match="duplicate H-Na/"):
        _core.build_parameters(duplicated)
    with pytest.raises(ParameterError, match="duplicate cation_mixing pair"):
        Parameters.from_mapping({**doc, "cation_mixing": [added, added]})
    with pytest.raises(ParameterError, match="unknown cation pair"):
        reg.active_parameter(start, "cation_pair_mixing", ("hydronium-cation", "calcium-cation"))


def native_document_entries(document: _core.ParameterDocument) -> list[tuple[object, ...]]:
    return [(e.identity, e.field, tuple(e.index), e.magnitude, e.unit) for e in document.entries]

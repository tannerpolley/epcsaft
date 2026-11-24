"""#608: parameter documents built from several packets equal the hand-edited documents they replace.
Input assembly only; no physical-validation claim."""
import json
import math
import warnings
from pathlib import Path
from typing import Any

import pytest

from epcsaft import DomainWarning, Mixture, ParameterError, Parameters, Pure, _core
from epcsaft import equilibrium as eq
from epcsaft.parameters import parameter_document
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
    def mixture(packet_id: str, names: list[str]) -> Mixture:
        return Mixture(Parameters.from_packets(PACKETS, [(packet_id, names)]))

    def fields(r: _core.OutsideDomain) -> dict[str, Any]:
        return dict(packet=r.packet, domain_id=r.domain_id, variable=r.variable, side=r.side,
                    bound=r.bound, value=r.value, unit=r.unit, **({"solute": r.solute} if r.solute else {}))

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
            assert fields(records[0]) == {"packet": refit, "domain_id": "hydrogen-chloride-298-15-k",
                "solute": "hydronium-cation",
                "variable": "solute_molality", "side": "above",
                "bound": 7.0,
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
            assert fields(records[0]) == {"packet": refit, "domain_id": "issue-575-chloride-brine-298-15-k",

                "variable": "ionic_strength", "side": "above",
                "bound": 9.9, "value": pytest.approx(10.2, rel=1e-12), "unit": "mole / kilogram"}

    nacl = mixture("figiel-yu-held-2025/4", ["water", "sodium-cation", "chloride-anion"])
    records = nacl.outside_domains([0.6, 0.2, 0.2])
    assert len(records) == 1
    assert records[0].variable == "ion_mole_fraction"
    assert records[0].side == "above"
    assert records[0].bound == 0.38
    assert records[0].value == pytest.approx(0.4, rel=1e-12)
    assert records[0].unit == "dimensionless"

    held = "held-2014-strategy2/2"
    nacl = mixture(held, ["water", "sodium-cation", "chloride-anion"])
    water_moles = 1 / nacl.parameters.molar_masses[0]
    total = water_moles + 2 * 0.0005
    records = nacl.outside_domains([water_moles / total, 0.0005 / total, 0.0005 / total])
    assert len(records) == 1
    assert fields(records[0]) == {"packet": held, "domain_id": "held-2014-chlorides-298-15-k-1-bar",
        "solute": "sodium-cation",
        "bound": 0.001, "variable": "solute_molality", "side": "below",
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
    assert fields(records[0]) == {"packet": other, "domain_id": "held-2014-chlorides-298-15-k-1-bar",
        "solute": "sodium-cation",
        "bound": 4.5, "variable": "solute_molality", "side": "above",
        "value": pytest.approx(5.0, rel=1e-12), "unit": "mole / kilogram"}


@pytest.fixture
def held_document_wide_model() -> Mixture:
    ids = ["water", "sodium-cation", "chloride-anion"]
    document = parameter_document(PACKETS, [("held-2014-strategy2/2", ids)])
    document["document_id"] = "document-conditions"
    document["domains"] = [{"domain_id": "document-conditions", "kind": "reported-conditions",
                           "temperature_min": 298.15, "temperature_max": 298.15, "temperature_unit": "kelvin",
                           "pressure_min": 1, "pressure_max": 1, "pressure_unit": "bar"}]
    return Mixture(Parameters.from_mapping(document))


def test_state_domain_warning(held_document_wide_model: Mixture) -> None:
    model = held_document_wide_model
    water = 1 / model.parameters.molar_masses[0]
    inside = [water / (water + 2), 1 / (water + 2), 1 / (water + 2)]
    outside = [water / (water + 10), 5 / (water + 10), 5 / (water + 10)]
    with warnings.catch_warnings():
        warnings.simplefilter("error", DomainWarning)
        assert model.state(298.15, P=100000, x=inside, phase="liquid").domain_violations is None
        with pytest.raises(DomainWarning, match="temperature=300 above 298.15 kelvin"):
            model.state(300, P=100000, x=inside, phase="liquid")

    with pytest.warns(DomainWarning) as caught:
        state = model.state(300, P=90000, x=outside, phase="liquid")
    assert len(caught) == 1
    assert state.domain_violations is not None
    by_variable = {r.variable: r for r in state.domain_violations}
    assert by_variable.keys() == {"temperature", "pressure", "solute_molality"}
    for key, value, bound, side, domain_id in (
        ("temperature", 300, 298.15, "above", "document-conditions"),
        ("pressure", 90000, 100000, "below", "document-conditions"),
        ("solute_molality", 5, 4.5, "above", "held-2014-chlorides-298-15-k-1-bar"),
    ):
        record = by_variable[key]
        assert (record.bound, record.side, record.packet, record.domain_id) == (
            bound, side, "held-2014-strategy2/2" if key == "solute_molality" else "document-conditions", domain_id)
        assert record.value == pytest.approx(value)

    with pytest.warns(DomainWarning) as caught:
        batch = model.states([298.15, 300, -1], P=100000, x=inside, phase="liquid")
    assert len(caught) == 1
    assert batch.domain_violations is not None
    assert [(r.row, r.variable) for r in batch.domain_violations] == [(1, "temperature")]

    with pytest.warns(DomainWarning, match="pressure=90000 below 100000 pascal"):
        actions = model.pressure_actions(298.15, 90000, x=inside, phase="liquid",
                                         observables=[_core.PropertyObservable.ResidualEnthalpy], directions=[])
    assert actions.domain_violations is not None
    assert [r.variable for r in actions.domain_violations] == ["pressure"]


def test_domain_warning_equilibrium_phases() -> None:
    document = packet("gross-sadowski-2001/3")
    document["domains"].append({"domain_id": "pressure-test", "kind": "fit-range", "pressure_max": 1, "pressure_unit": "bar"})
    model = Mixture(Parameters.from_mapping(document, components=["methane", "ethane"]))
    feed = eq.Amounts({"methane": 0.5, "ethane": 0.5})
    with pytest.warns(DomainWarning) as caught:
        result = eq.tp_flash(model, T=190, P=1e6, feed=feed)
    assert len(caught) == 1
    assert result.success and result.phase_count == 2
    assert result.domain_violations is not None
    assert [(r.phase, r.variable, r.value, r.bound) for r in result.domain_violations] == [
        (0, "pressure", 1e6, 1e5), (1, "pressure", 1e6, 1e5)]


def test_unknown_domain_reports_nothing() -> None:
    document = packet("gross-2001-methane-ethane/3")
    document["domains"] = [d for d in document["domains"] if d["kind"] == "unknown"]
    model = Mixture(Parameters.from_mapping(document))
    with warnings.catch_warnings():
        warnings.simplefilter("error", DomainWarning)
        assert model.state(330, P=1e5, x=[0.5, 0.5]).domain_violations is None
    with pytest.raises(ParameterError, match="phase mole fractions"):
        model.outside_domains([1])

def test_materialized_domain_component_selection(tmp_path: Path) -> None:
    source = packet("gross-sadowski-2001/3")
    source["domains"].append({"domain_id": "document-temperature", "kind": "fit-range",
                              "temperature_min": 80, "temperature_unit": "kelvin"})
    path = tmp_path / "test/1/parameters/parameters.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(source))
    doc: dict[str, Any] = parameter_document(tmp_path, [("test/1", ["methane", "propane"])])
    model = Pure(Parameters.from_mapping(doc, components=["propane"]))
    assert model.outside_domains([1], T=84) == []
    assert [(r.packet, r.domain_id, r.bound) for r in model.outside_domains([1], T=79)] == [
        ("test/1", "document-temperature", 80),
    ]

    doc["materialization"][0].pop("document_wide_domains")
    legacy = Pure(Parameters.from_mapping(doc, components=["propane"]))
    assert legacy.outside_domains([1], T=79) == []


def test_materialized_domain_ids_are_packet_local(tmp_path: Path) -> None:
    names = ["methane", "ethane"]
    for name in names:
        doc = packet("gross-2001-methane-ethane/3")
        doc["domains"].append({"domain_id": "local-range", "kind": "fit-range",
                               "temperature_max": 290, "temperature_unit": "kelvin",
                               "pressure_min": 1.1, "pressure_unit": "bar"})
        if name == "ethane":
            component = next(c for c in doc["components"] if c["component_id"] == name)
            component["coefficients"][0]["provenance"] = {"domain_id": "local-range"}
        path = tmp_path / name / "1/parameters/parameters.json"
        path.parent.mkdir(parents=True)
        path.write_text(json.dumps(doc))

    saved = parameter_document(tmp_path, [(f"{name}/1", [name]) for name in names])
    for components, x in ((names, [0.5, 0.5]), (["methane"], [1])):
        model = Mixture(Parameters.from_mapping(saved, components=components))
        with pytest.warns(DomainWarning) as caught:
            state = model.state(300, P=100000, x=x, phase="vapor")
        assert len(caught) == 1 and state.domain_violations is not None
        assert [(r.packet, r.domain_id, r.variable, r.bound) for r in state.domain_violations] == [
            ("methane/1", "local-range", "temperature", 290),
            ("methane/1", "local-range", "pressure", 110000),
        ]


@pytest.mark.parametrize("bound,unit", [(1.1, "bar"), (0.11, "megapascal")])
def test_pressure_domain_at_single_point(bound: float, unit: str) -> None:
    doc = packet("figiel-2025-reference-electrolytes/3")
    doc["domains"].append({"domain_id": "document-pressure", "kind": "reported-conditions",
                           "pressure_min": bound, "pressure_max": bound, "pressure_unit": unit})
    model = Mixture(Parameters.from_mapping(doc, components=["water", "sodium-cation", "chloride-anion"]))
    x = [0.98, 0.01, 0.01]
    with warnings.catch_warnings():
        warnings.simplefilter("error", DomainWarning)
        assert model.state(295, P=110000, x=x, phase="liquid").domain_violations is None
    for target, side in ((0, "below"), (math.inf, "above")):
        adjacent = math.nextafter(110000, target)
        (record,) = model.outside_domains(x, T=295, P=adjacent)
        assert (record.variable, record.value, record.bound, record.side) == ("pressure", adjacent, 110000, side)


@pytest.mark.parametrize("override", [None, 0.0])
def test_pair_coefficient_domain_is_deferred(override: float | None) -> None:
    a, b = "carbon-dioxide", "ethane"
    packets = [("gross-sadowski-2001/3", [a, b])]
    model = Mixture(Parameters.from_packets(PACKETS, packets, k_ij=None if override is None else {(a, b): override}))
    with warnings.catch_warnings():
        warnings.simplefilter("error", DomainWarning)
        assert model.state(300, P=1e5, x=[0.5, 0.5]).domain_violations is None


@pytest.mark.parametrize("document_local", [False, True])
def test_omitted_cation_mixing_domain(document_local: bool) -> None:
    doc = packet("figiel-2025-chloride-mixture-refit/3")
    domains = doc["domains"] if document_local else doc["materialization"][0]["domains"]
    domains.append({"domain_id": "omitted-cation-mixing", "kind": "fit-range",
                    "temperature_max": 295, "temperature_unit": "kelvin"})
    doc["cation_mixing"]["coefficients"][0]["provenance"] = {"domain_id": "omitted-cation-mixing"}
    model = Mixture(Parameters.from_mapping(doc, components=["water", "sodium-cation", "chloride-anion"]))
    with warnings.catch_warnings():
        warnings.simplefilter("error", DomainWarning)
        assert model.state(298.15, P=1e5, x=[1, 0, 0], phase="liquid").domain_violations is None


@pytest.mark.parametrize("record", ["born", "ethanol-constant", "ethanol-value"])
def test_cited_temperature_pressure_domain_is_deferred(record: str) -> None:
    doc = packet("figiel-2025-reference-electrolytes/3")
    doc["domains"].append({"domain_id": "cited-conditions", "kind": "fit-range",
                           "temperature_max": 295, "temperature_unit": "kelvin",
                           "pressure_max": 0.9, "pressure_unit": "bar"})
    citation = {"domain_id": "cited-conditions"}
    if record == "born":
        born = next(r for r in doc["model_families"] if r["choice"] == "born")
        born["provenance"] = citation
    else:
        constant = next(r["constant"] for r in doc["correlations"]
                        if r["component_id"] == "ethanol" and r["family"] == "relative_permittivity")
        if record == "ethanol-value":
            constant["value"]["evidence"] = [{"original": {"provenance": [citation]}}]
        else:
            constant["provenance"] = citation

    ids = ["water", "sodium-cation", "chloride-anion"]
    for components, x in ((ids, [0.98, 0.01, 0.01]), ([*ids, "ethanol"], [0.97, 0.01, 0.01, 0.01])):
        model = Mixture(Parameters.from_mapping(doc, components=components))
        assert model.outside_domains(x, T=298.15, P=1e5) == []
        with warnings.catch_warnings():
            warnings.simplefilter("error", DomainWarning)
            assert model.state(298.15, P=1e5, x=x, phase="liquid").domain_violations is None


def test_document_domain_with_inherited_materialization() -> None:
    source = PACKETS / "figiel-2025-chloride-mixture-refit/3/parameters/parameters.json"
    model = Mixture(Parameters.from_json(source))
    water = [float(name == "water") for name in model.component_ids]
    records = model.outside_domains(water, T=300, P=90000)
    assert {(r.domain_id, r.variable, r.bound, r.side) for r in records} == {
        ("issue-209-provisional-binary-aqueous", "temperature", 298.15, "above"),
        ("issue-209-provisional-binary-aqueous", "pressure", 100000, "below"),
    }


def test_fixed_density_pressure_domain() -> None:
    doc = packet("gross-2001-propane/3")
    model = Pure(Parameters.from_mapping(doc))
    temperature, density = 300, 1000
    pressure = model.state(temperature, rho=density).pressure
    doc["domains"].append({"domain_id": "pressure-test", "kind": "fit-range",
                           "pressure_min": pressure, "pressure_max": pressure, "pressure_unit": "pascal"})
    model = Pure(Parameters.from_mapping(doc))
    with warnings.catch_warnings():
        warnings.simplefilter("error", DomainWarning)
        assert model.log_fugacity_terms(temperature, density).domain_violations is None
    with pytest.warns(DomainWarning) as caught:
        action = model.density_actions(temperature, 500, observables=[_core.PropertyObservable.ResidualEnthalpy],
                                       directions=[])
    assert len(caught) == 1 and action.domain_violations is not None
    assert [(r.variable, r.bound, r.side) for r in action.domain_violations] == [("pressure", pressure, "below")]


def test_compiled_saturation_pressure_domain() -> None:
    doc = packet("gross-2001-propane/3")
    doc["domains"].append({"domain_id": "pressure-test", "kind": "fit-range", "pressure_max": 1, "pressure_unit": "bar"})
    model = Pure(Parameters.from_mapping(doc))
    compiled = eq.compile_problem(model, model.saturation(300))
    with pytest.warns(DomainWarning) as caught:
        result = compiled.solve_saturation()
    assert len(caught) == 1 and result.domain_violations is not None
    assert [(r.phase, r.variable, r.value) for r in result.domain_violations] == [
        (0, "pressure", result.pressure), (1, "pressure", result.pressure)]
    with warnings.catch_warnings():
        warnings.simplefilter("error", DomainWarning)
        with pytest.raises(DomainWarning, match="pressure"):
            compiled.solve()

    request = eq.SolvedStateActionRequest([eq.SolvedStateObservable(eq.SolvedStateObservableKind.Pressure)],
                                         [eq.SolvedStateActionDirection(1, 0, [0], [])])
    with pytest.warns(DomainWarning) as caught:
        actions = eq.solved_state_actions(compiled, result, request)
    assert len(caught) == 1 and actions.domain_violations is not None
    assert [r.phase for r in actions.domain_violations] == [0, 1]


def test_reference_activity_domains(held_document_wide_model: Mixture) -> None:
    ids = ["water", "sodium-cation", "chloride-anion"]
    model = held_document_wide_model
    with warnings.catch_warnings():
        warnings.simplefilter("error", DomainWarning)
        assert model.reference(298.15, 100000, solvent_fractions=[1, 0, 0]).domain_violations is None
        with pytest.raises(DomainWarning, match="temperature"):
            model.reference(300, 90000, solvent_fractions=[1, 0, 0])
    with pytest.warns(DomainWarning) as caught:
        solvent_reference = model.reference(300, 90000, solvent_fractions=[1, 0, 0])
    assert len(caught) == 1 and solvent_reference.domain_violations is not None
    assert {r.variable for r in solvent_reference.domain_violations} == {"temperature", "pressure"}
    with warnings.catch_warnings():
        warnings.simplefilter("error", DomainWarning)
        (refused,) = solvent_reference.henry_constants(solutes=["sodium-cation"])
    assert refused.status == _core.NeutralReferenceStatus.NonNeutralDeclaration and refused.domain_violations is None

    with pytest.warns(DomainWarning) as caught:
        reference = model.reference(300, 90000, solvent_fractions=[1, 0, 0], ions=(ids[1], ids[2]), stoichiometry=(1, 1))
    assert len(caught) == 1 and reference.domain_violations is not None
    assert {r.variable for r in reference.domain_violations} == {"temperature", "pressure"}
    water = 1 / model.parameters.molar_masses[0]
    with pytest.warns(DomainWarning) as caught:
        activity = reference.activity(5)
    assert len(caught) == 1 and activity.domain_violations is not None
    assert {r.variable for r in activity.domain_violations} == {"temperature", "pressure", "solute_molality"}
    with pytest.warns(DomainWarning) as caught:
        solvent = reference.solvent_activity(x=[water / (water + 10), 5 / (water + 10), 5 / (water + 10)])
    assert len(caught) == 1 and solvent.domain_violations is not None
    assert {r.variable for r in solvent.domain_violations} == {"temperature", "pressure", "solute_molality"}


def test_henry_constant_pressure_domain() -> None:
    doc = packet("gross-2001-methane-ethane/3")
    doc["domains"].append({"domain_id": "pressure-test", "kind": "fit-range", "pressure_max": 1, "pressure_unit": "bar"})
    model = Mixture(Parameters.from_mapping(doc, components=["methane", "ethane"]))
    with pytest.warns(DomainWarning) as caught:
        reference = model.reference(200, 1e6, solvent_fractions=[0, 1])
    assert len(caught) == 1
    with pytest.warns(DomainWarning) as caught:
        (henry,) = reference.henry_constants(solutes=["methane"])
    assert len(caught) == 1 and henry.henry_constant_px is not None and henry.domain_violations is not None
    assert [(r.variable, r.value, r.bound) for r in henry.domain_violations] == [("pressure", 1e6, 1e5)]


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

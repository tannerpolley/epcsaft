import copy
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import pytest

from epcsaft import (DensityBranch, IdealCorrelation, IdealInterval, IdealPolynomial, Mixture, ParameterError, Parameters,
                     Pure, ThermochemistryRecord)
from tests.equilibrium.test_trace_reliability import _mea_mapping
from tests.native_inputs import packet_metadata, same
from tests.packets import parameters_json

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
_MEA_ORDER = (
    "water",
    "monoethanolamine",
    "carbon-dioxide",
    "protonated-monoethanolamine",
    "carbamate-anion",
    "bicarbonate-anion",
    "carbonate-anion",
    "hydronium-cation",
    "hydroxide-anion",
)


def close(left: float, right: float, relative: float) -> None:
    assert math.isclose(left, right, rel_tol=relative, abs_tol=relative)


def test_relabelling_and_state_identities(evidence: Any) -> None:
    packet = parameters_json("mea-co2-h2o-nine-species-estimation", 4)
    declared = _mea_mapping(packet)
    order = _MEA_ORDER
    reverse = tuple(reversed(order))
    x = (0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0)
    temperature = 313.15
    density = 3.0e4
    pressure_target = 1.0e5
    phase = "liquid"
    first = Mixture(Parameters.from_mapping(declared, components=order)).state(
        temperature, rho=density, x=x
    )
    second = Mixture(Parameters.from_mapping(declared, components=reverse)).state(
        temperature, rho=density, x=tuple(reversed(x))
    )
    first_fugacity = first.log_fugacity_coefficient
    second_fugacity = second.log_fugacity_coefficient
    fugacity_available = first_fugacity is not None and second_fugacity is not None
    mixture = Mixture(Parameters.from_mapping(declared, components=order))
    state = mixture.state(temperature, P=pressure_target, x=x, phase=phase)
    anchored = mixture.state(
        temperature, P=pressure_target, x=x, phase=phase, anchor=state
    )
    parameter_metadata = {**packet_metadata(packet, order), "relabelled_selected_components": list(reverse)}
    inputs = {
        "T_K": temperature,
        "rho_mol_m3": density,
        "pressure_target_Pa": pressure_target,
        "phase": phase,
        "composition": list(x),
        "relabelled_composition": list(reversed(x)),
        "components": list(order),
        "relabelled_components": list(reverse),
        "fugacity_available": fugacity_available,
        "_parameters": parameter_metadata,
    }

    def identity(name: str, left: float | None, right: float | None, atol: float, rtol: float, **more: Any) -> None:
        same(evidence, f"state_identity/{name}", left, right, atol, rtol, inputs={**inputs, **more}, basis="Relationship")

    # Relabelling invariance of the scalar state and reversal of the component vectors.
    for key in ("pressure", "residual_helmholtz", "compressibility_factor", "bulk_relative_permittivity"):
        identity(key, getattr(first, key), getattr(second, key), 1e-12, 1e-12)
    assert (first_fugacity is None) == (second_fugacity is None)
    euler_terms = []
    for index, component in enumerate(order):
        mirror = len(order) - index - 1
        potential = first.residual_chemical_potential_over_rt[index]
        assert potential is not None
        if x[index] != 0:
            euler_terms.append(x[index] * potential)
        identity(f"residual_chemical_potential/{index}", potential,
                 second.residual_chemical_potential_over_rt[mirror], 1e-12, 1e-12, component=component)
        if first_fugacity is not None and second_fugacity is not None:
            identity(f"log_fugacity_coefficient/{index}", first_fugacity[index], second_fugacity[mirror],
                     1e-12, 1e-12, component=component)
    identity("euler_identity", sum(euler_terms), first.residual_helmholtz + first.compressibility_factor - 1.0,
             1e-12, 1e-12, identity="sum(x_i * mu_i^res/RT) = a^res + Z - 1")
    identity("codata_pressure", first.pressure,
             first.compressibility_factor * first.molar_density * 8.31446261815324 * temperature, 0.0, 1e-12)
    identity("pressure_anchor", state.pressure, pressure_target, 1e-9, 1e-9)
    identity("anchored_density", anchored.molar_density, state.molar_density, 1e-10, 1e-10)


def test_pure_is_one_component_mixture(evidence: Any) -> None:
    packet = parameters_json("gross-2001-methane-ethane", 3)
    temperature = 300.0
    density = 1.0e4
    composition = (0.0, 1.0)
    pure_composition = (1.0,)
    pure_components = ("ethane",)
    mixture_components = ("methane", "ethane")
    pure = Pure(Parameters.from_json(packet, components=pure_components)).state(
        temperature, rho=density
    )
    mixture = Mixture(Parameters.from_json(packet, components=mixture_components)).state(
        temperature, rho=density, x=composition
    )
    parameter_metadata = packet_metadata(packet, mixture_components)
    evidence.compare(
        "pure_mixture_identity", quantity="pressure", model=float(pure.pressure),
        reference=float(mixture.pressure), unit="Pa",
        inputs={"T_K": temperature, "rho_mol_m3": density, "composition": list(composition),
                "pure_components": list(pure_components), "pure_composition": list(pure_composition),
                "mixture_components": list(mixture_components),
                "_parameters": parameter_metadata},
        criterion={"type": "scalar", "rule": "symmetric_close", "atol": 0.0, "rtol": 1e-12,
                   "scale": "max(abs(model), abs(reference))",
                   "combining_rule": "max(atol, rtol*scale)"},
        basis="Relationship", method="pure/mixture identity",
        claim="Pure pressure equals one-component mixture pressure.",
        relationship_latex=r"P_{pure}(T,\rho) = P_{mixture}(T,\rho,x_{ethane}=1)", quantity_latex="P",
        inputs_latex=r"T,\rho,\mathbf{x}",
    )
    close(pure.pressure, mixture.pressure, 1e-12)


def test_neutral_selection_is_the_ion_free_electrolyte_limit() -> None:
    # DH and Born vanish exactly without ions, so a water-only selection of an
    # electrolyte packet must equal the electrolyte selection at zero ion amount.
    packet = parameters_json("figiel-yu-held-2025", 3)
    water = Pure(Parameters.from_json(packet, components=("water",)))
    salt = Mixture(Parameters.from_json(packet, components=("water", "sodium-cation", "chloride-anion")))
    for density in (20.0, 55000.0):
        pure = water.state(298.15, rho=density)
        limit = salt.state(298.15, rho=density, x=(1.0, 0.0, 0.0))
        close(pure.pressure, limit.pressure, 1e-12)
        close(pure.residual_enthalpy, limit.residual_enthalpy, 1e-12)
        pure_potential, limit_potential = pure.residual_chemical_potential_over_rt[0], limit.residual_chemical_potential_over_rt[0]
        assert pure_potential is not None and limit_potential is not None
        close(pure_potential, limit_potential, 1e-12)
    liquid = water.state(298.15, P=1e5, phase="liquid")
    close(liquid.molar_density, salt.state(298.15, P=1e5, x=(1.0, 0.0, 0.0), phase="liquid").molar_density, 1e-12)


def test_physical_state_refuses_a_charged_feed() -> None:
    # full-eos.md section 12: an ordinary physical ionic state must be electroneutral (#563).
    packet = parameters_json("figiel-yu-held-2025", 3)
    salt = Mixture(Parameters.from_json(packet, components=("water", "sodium-cation", "chloride-anion")))
    with pytest.raises(ParameterError, match="charge sum 0.004"):
        salt.state(298.15, P=1e5, x=(0.98, 0.012, 0.008), phase="liquid")
    with pytest.raises(ParameterError, match="charge sum 0.004"):
        salt.state(298.15, rho=55000.0, x=(0.98, 0.012, 0.008))


@pytest.mark.parametrize(
    ("temperature", "pressure", "phase", "reference_density"),
    [(308.8, 5147500.0, "vapor", 5930.399614704038),
     (308.8, 5147500.0, "liquid", 6828.295025027434),
     (308.0, 5057915.461281232, "vapor", 5124.770478882947),
     (308.0, 5057915.461281232, "liquid", 7160.554342626920)],
)
def test_pressure_branch_selects_sub_cell_loop_root(
    evidence: Any, temperature: float, pressure: float, phase: str, reference_density: float,
) -> None:
    """#89: pure ethane at 308.8 K / 5147500 Pa has two stable roots 0.14 apart in log density,
    inside one 40-cell scan interval. References are the retained independent fixed-density
    bisection roots of analyses/2026-density-closure-diagnosis/results/root-controls.csv.
    At 308 K the liquid root coincides with a scan node; the independent fixed-density
    bisection roots still require two distinct branches, regardless of that coincidence."""
    packet = parameters_json("gross-2001-methane-ethane", 3)
    model = Pure(Parameters.from_json(packet, components=("ethane",)))
    state = model.state(temperature, P=pressure, phase=phase)
    evidence.compare(
        f"loop_branch_{phase}", quantity="molar_density", model=float(state.molar_density),
        reference=reference_density, unit="mol/m3",
        inputs={"T_K": temperature, "P_Pa": pressure, "phase": phase,
                "stable_root_count": state.stable_root_count,
                "density_branch": str(state.density_branch),
                "pressure_density_derivative": float(state.pressure_density_derivative),
                "_parameters": packet_metadata(packet, ("ethane",))},
        criterion={"type": "scalar", "rule": "symmetric_close", "atol": 0.0, "rtol": 1e-8,
                   "scale": "max(abs(model), abs(reference))",
                   "combining_rule": "max(atol, rtol*scale)"},
        basis="Relationship", method="independent fixed-density bisection",
        claim="A declared branch returns its extreme stable root when the loop is narrower than a scan cell.",
        relationship_latex=r"\rho_{branch}(T,P^\star)=\rho_{ref}", quantity_latex=r"\rho",
        inputs_latex=r"T,P^\star,\text{branch}",
    )
    close(state.molar_density, reference_density, 1e-8)
    assert state.density_branch == getattr(DensityBranch, phase.capitalize()) and state.stable_root_count == 2
    assert state.pressure_density_derivative > 0.0
    with pytest.raises(RuntimeError, match="ambiguous"):
        model.state(temperature, P=pressure)


def test_mapping_fails_closed_by_name() -> None:
    with pytest.raises(ParameterError, match="model/electrolyte/ionic_dispersion"):
        Parameters.from_json(
            parameters_json("figiel-2025-reference-electrolytes", 2),
            components=("water", "lithium-cation", "chloride-anion"),
        )
    with pytest.raises(ParameterError, match="pair/ethanol/methanol/k_ij"):
        Parameters.from_json(parameters_json("figiel-2025-reference-electrolytes", 3), components=("methanol", "ethanol"))
    ionic_mapping = copy.deepcopy(json.loads(parameters_json("figiel-2025-reference-electrolytes", 3).read_text()))
    ionic_mapping["pairs"].append({
        "component_id_a": "lithium-cation",
        "component_id_b": "sodium-cation",
        "coefficients": [{
            "identity": "pair/lithium-cation/sodium-cation/k_ij",
            "family": "k_ij",
            "value": {"magnitude": 0.0, "unit": "dimensionless"},
        }],
    })
    with pytest.raises(ParameterError, match="lithium-cation/sodium-cation.*exclude-same-sign-ion-pairs"):
        Parameters.from_mapping(
            ionic_mapping,
            components=("water", "lithium-cation", "sodium-cation", "chloride-anion"),
        )
    mapping = copy.deepcopy(json.loads(parameters_json("gross-2001-methane-ethane", 3).read_text()))
    mapping["pair"] = mapping.pop("pairs")  # an unread top-level field is refused, not dropped with its k_ij
    with pytest.raises(ParameterError, match="does not read: pair$"):
        Parameters.from_mapping(mapping)
    mapping = copy.deepcopy(json.loads(parameters_json("gross-2001-methane-ethane", 3).read_text()))
    mapping["components"][0]["fixed"]["molar_mass"]["value"]["unit"] = "gram"
    with pytest.raises(ParameterError, match="(?=.*component/ethane/molar_mass)(?=.*gram)"):
        Parameters.from_mapping(mapping, components=("ethane",))
    parameters = Parameters.from_json(parameters_json("gross-2001-methane-ethane", 3), components=("methane", "ethane"))
    with pytest.raises(ParameterError):
        Mixture(parameters).state(300.0, P=1e5, x=(0.5, 0.5), phase="gas")


def test_named_composition_and_molar_masses_follow_component_ids() -> None:
    # Without components=, document order (ethane, methane) sets component_ids; x=(0.3, 0.7) would mean 30 % ethane.
    parameters = Parameters.from_json(parameters_json("gross-2001-methane-ethane", 3))
    model = Mixture(parameters)
    assert parameters.component_ids == ("ethane", "methane")
    assert parameters.molar_masses == (0.03007, 0.016043)  # kg/mol, the packet's declared values
    named = model.state(200.0, P=2e6, x={"methane": 0.3, "ethane": 0.7}, phase="liquid")
    assert named.molar_density == model.state(200.0, P=2e6, x=(0.7, 0.3), phase="liquid").molar_density
    with pytest.raises(ParameterError, match="propane; expected one of ethane, methane"):
        model.state(200.0, P=2e6, x={"propane": 1.0}, phase="liquid")
    with pytest.raises(ParameterError, match="unknown species None"):  # None is a name, not a not-found marker
        model.state(200.0, P=2e6, x={"methane": 1.0, None: 1.0}, phase="liquid")  # type: ignore[dict-item]
    with pytest.raises(AttributeError):  # read from the native record, not a settable copy
        parameters.molar_masses = (1.0, 2.0)  # type: ignore[misc]



def test_declared_model_preserves_association_equations(evidence: Any) -> None:
    # Fixed-density equation control, outside the source's reported 1-bar states.
    path = parameters_json("figiel-yu-held-2025", 3)
    document = json.loads(path.read_text())
    components = ("water", "sodium-cation", "chloride-anion")
    state = Mixture(Parameters.from_mapping(document, components=components)).state(
        298.15, rho=56000.0, x=(0.98, 0.01, 0.01)
    )
    # Symmetric water 2B reduction: B X² + X - 1 = 0.
    bond_density = 536.45118457612296
    fraction = 2.0 / (1.0 + math.sqrt(1.0 + 4.0 * bond_density))
    association = 2.0 * 0.98 * (math.log(fraction) - fraction / 2.0 + 0.5)
    assert len(state.site_fractions) == 2
    assert state.site_fractions == pytest.approx((fraction, fraction), rel=1e-11, abs=1e-12)
    assert state.association == pytest.approx(association, rel=1e-11, abs=1e-12)
    evidence.compare(
        "declared_model/association", quantity="reduced association Helmholtz energy",
        model=state.association, reference=association, unit="1",
        inputs={"T_K": 298.15, "rho_mol_m3": 56000.0, "x": [0.98, 0.01, 0.01],
                "components": components, "bond_density": bond_density,
                "parameter_sha256": hashlib.sha256(path.read_bytes()).hexdigest()},
        criterion={"rule": "symmetric_close", "atol": 1e-12, "rtol": 1e-11},
        basis="Relationship", method="independent symmetric 2B scalar root",
        claim="The declared general-site choice preserves the association equations.",
    )

    relabelled = copy.deepcopy(document)
    relabelled["document_id"] = "display-only"
    for family in relabelled["model_families"]:
        family["family_id"] = "display/" + family["kind"]
    same = Mixture(Parameters.from_mapping(relabelled, components=components)).state(
        298.15, rho=56000.0, x=(0.98, 0.01, 0.01)
    )
    for field in ("pressure", "association", "residual_helmholtz",
                  "pressure_density_derivative", "site_fractions",
                  "residual_chemical_potential_over_rt"):
        assert getattr(same, field) == pytest.approx(getattr(state, field), rel=1e-12, abs=1e-12)

    for kind, choice in (("association", "none"), ("base", "simplified-pc-saft")):
        invalid = copy.deepcopy(document)
        next(item for item in invalid["model_families"] if item["kind"] == kind)["choice"] = choice
        with pytest.raises(ParameterError, match=kind):
            Parameters.from_mapping(invalid, components=components)
    for extra in ({"kind": "polar", "choice": "dipole"},
                  {"kind": "association", "choice": "general-site"}):
        invalid = copy.deepcopy(document)
        invalid["model_families"].append(extra)
        with pytest.raises(ParameterError, match=extra["kind"]):
            Parameters.from_mapping(invalid, components=components)
    malformed_first = copy.deepcopy(document)
    malformed_first["model_families"].insert(0, {"kind": "association"})
    with pytest.raises(ParameterError, match="duplicate model choice association"):
        Parameters.from_mapping(malformed_first, components=components)

    ambiguous = copy.deepcopy(document)
    water = copy.deepcopy(next(item for item in document["components"] if item["component_id"] == "water"))
    next(item for item in water["coefficients"] if item["family"] == "dispersion_energy_over_k")["value"]["magnitude"] *= 0.9
    ambiguous["components"].append(water)
    with pytest.raises(ParameterError, match="duplicate component_id water"):
        Parameters.from_mapping(ambiguous, components=components)

    unbonded = copy.deepcopy(document)
    unbonded["topology"] = {"sites": [], "edges": []}
    next(item for item in unbonded["model_families"] if item["kind"] == "association")["choice"] = "none"
    assert Mixture(Parameters.from_mapping(unbonded, components=components)).state(
        298.15, rho=56000.0, x=(0.98, 0.01, 0.01)
    ).association == 0.0

    missing_charge = copy.deepcopy(document)
    del next(item for item in missing_charge["components"] if item["component_id"] == "water")["fixed"]["charge_number"]
    with pytest.raises(ParameterError, match="missing component/water/charge_number"):
        Parameters.from_mapping(missing_charge, components=components)
    malformed = copy.deepcopy(document)
    malformed["components"][0]["component_id"] = 7
    with pytest.raises(ParameterError, match=r"components\[0\]\.component_id"):
        Parameters.from_mapping(malformed, components=components)
    with pytest.raises(ParameterError, match=r"component selection\[0\]"):
        Parameters.from_mapping(document, components=("",))
    for section, member, field, message in (
        ("correlations", 0, "component_id", r"correlation\[0\]\.component_id"),
        ("pairs", 0, "component_id_a", r"pair\[0\]\.component_id_a"),
    ):
        malformed = copy.deepcopy(document)
        malformed[section][member][field] = 7
        with pytest.raises(ParameterError, match=message):
            Parameters.from_mapping(malformed, components=components)
    malformed = copy.deepcopy(document)
    malformed["topology"]["sites"][0]["site_id"] = 7
    with pytest.raises(ParameterError, match=r"topology site\[0\]\.site_id"):
        Parameters.from_mapping(malformed, components=components)
    malformed = copy.deepcopy(document)
    malformed["topology"]["edges"][0]["endpoint_a"]["component_id"] = 7
    with pytest.raises(ParameterError, match=r"topology edge\[0\]\.endpoint_a\.component_id"):
        Parameters.from_mapping(malformed, components=components)


def test_wolbach_sandler_rule_matches_source_cross_volumes() -> None:
    """Gross2002 Eq. 3 with water sigma(T) from Khudaida Table 5.

    At 293.15 K the reference is Khudaida bundle 2, whose cross volumes independent
    historical code evaluated once from the same equation. At 353.15 K the
    reference volumes follow the water sigma correlation, so a size factor frozen
    at one temperature fails. Residual entropy falsifies the temperature chain.
    """
    rule: dict[str, Any] = json.loads(parameters_json("khudaida-2026-figure-2-electrolyte-lle", 7).read_text())
    bundle = REPOSITORY_ROOT / "analyses/reproduction/2026-khudaida/inputs/parameters/bundle/khudaida-2026-figure-2-electrolyte-lle/2/association.csv"
    stored = {(r["component_id_a"], r["site_id_a"], r["component_id_b"], r["site_id_b"], r["family"]): float(r["value"])
              for r in csv.DictReader(bundle.open())}
    pure = {c: (stored[(c, "a", c, "b", "association_energy_over_k")], stored[(c, "a", c, "b", "association_volume")])
            for c in ("water", "ethanol", "isobutanol")}
    sigma = {"ethanol": 3.1771, "isobutanol": 3.5561}

    def explicit(temperature: float | None) -> dict[str, Any]:
        mapping = copy.deepcopy(rule)
        t = 293.15 if temperature is None else temperature
        sigma["water"] = 2.7927 + 10.11 * math.exp(-0.01775 * t) - 1.417 * math.exp(-0.01146 * t)
        for edge in mapping["topology"]["edges"]:
            if edge["source"]["kind"] != "combining-rule":
                continue
            (i, s_i), (j, s_j) = ((edge[e]["component_id"], edge[e]["site_id"]) for e in ("endpoint_a", "endpoint_b"))
            key = next(k for k in ((i, s_i, j, s_j), (j, s_j, i, s_i)) if (*k, "association_volume") in stored)
            factor = (math.sqrt(sigma[i] * sigma[j]) / (0.5 * (sigma[i] + sigma[j]))) ** 3
            volume = stored[(*key, "association_volume")] if temperature is None else math.sqrt(pure[i][1] * pure[j][1]) * factor
            tag = f"association/{i}/{s_i}/{j}/{s_j}"
            edge["energy_over_k"] = {"identity": f"{tag}/energy_over_k", "value": {"magnitude": 0.5 * (pure[i][0] + pure[j][0]), "unit": "kelvin"}}
            edge["volume"] = {"identity": f"{tag}/volume", "value": {"magnitude": volume, "unit": "dimensionless"}}
            edge["source"] = {"kind": "explicit", "provenance": edge["source"]["provenance"]}
        return mapping

    composition = (0.01, 0.25, 0.25, 0.01, 0.48)
    for temperature, reference in ((293.15, explicit(None)), (353.15, explicit(353.15))):
        model = Mixture(Parameters.from_mapping(rule))
        state = model.state(temperature, rho=1.5e4, x=composition)
        frozen = Mixture(Parameters.from_mapping(reference)).state(temperature, rho=1.5e4, x=composition)
        for field in ("pressure", "association", "residual_helmholtz"):
            close(getattr(state, field), getattr(frozen, field), 1e-12)
        for left, right in zip(state.residual_chemical_potential_over_rt, frozen.residual_chemical_potential_over_rt):
            assert left is not None and right is not None
            close(left, right, 1e-12)
        step = 1e-3
        helmholtz = [8.314462618 * t * model.state(t, rho=1.5e4, x=composition).residual_helmholtz
                     for t in (temperature - step, temperature + step)]
        close(state.residual_entropy, -(helmholtz[1] - helmholtz[0]) / (2 * step), 1e-7)


def test_association_topology_is_order_independent_and_typed() -> None:
    mapping = json.loads(
        parameters_json("khudaida-2026-figure-2-electrolyte-lle", 7).read_text()
    )
    reversed_mapping = copy.deepcopy(mapping)
    reversed_mapping["topology"]["edges"].reverse()
    first_parameters = Parameters.from_mapping(mapping)
    second_parameters = Parameters.from_mapping(reversed_mapping)
    assert first_parameters.component_ids == second_parameters.component_ids
    assert first_parameters.site_ids == second_parameters.site_ids
    composition = (0.01, 0.25, 0.25, 0.01, 0.48)
    first = Mixture(first_parameters).state(293.15, rho=1.0e4, x=composition)
    second = Mixture(second_parameters).state(293.15, rho=1.0e4, x=composition)
    for field in (
        "pressure", "residual_helmholtz", "compressibility_factor", "association",
        "debye_huckel", "born", "bulk_relative_permittivity",
    ):
        assert getattr(first, field) == pytest.approx(getattr(second, field), rel=1e-12, abs=1e-12)
    assert first.residual_chemical_potential_over_rt == pytest.approx(
        second.residual_chemical_potential_over_rt, rel=1e-12, abs=1e-12
    )
    assert first.site_fractions == pytest.approx(second.site_fractions, rel=1e-12, abs=1e-12)

    missing = copy.deepcopy(mapping)
    missing["topology"]["edges"][1]["source"]["inputs"]["energy_over_k"][1] = "missing"
    with pytest.raises(ParameterError, match="missing association input missing"):
        Parameters.from_mapping(missing)
    duplicate = copy.deepcopy(mapping)
    duplicate["topology"]["edges"][5]["energy_over_k"]["identity"] = (
        duplicate["topology"]["edges"][0]["energy_over_k"]["identity"]
    )
    with pytest.raises(ParameterError, match="duplicate association identity"):
        Parameters.from_mapping(duplicate)
    non_string = copy.deepcopy(mapping)
    non_string["topology"]["edges"][0]["energy_over_k"]["identity"] = 123
    with pytest.raises(ParameterError, match="association edge identities must be strings"):
        Parameters.from_mapping(non_string)
    cross_family = copy.deepcopy(mapping)
    cross_family["topology"]["edges"][1]["source"]["inputs"]["energy_over_k"][0] = (
        cross_family["topology"]["edges"][0]["volume"]["identity"]
    )
    with pytest.raises(ParameterError, match="energy_over_k input has wrong family"):
        Parameters.from_mapping(cross_family)

    subset = Parameters.from_mapping(mapping, components=("ethanol", "water"))
    assert subset.component_ids == ("ethanol", "water")
    invalid_subset = copy.deepcopy(mapping)
    invalid_subset["topology"]["edges"][2]["endpoint_b"]["site_id"] = "missing"
    with pytest.raises(ParameterError, match="missing topology site water/missing"):
        Parameters.from_mapping(invalid_subset, components=("ethanol", "water"))


def test_density_anchor_outside_domain_keeps_native_convergence_failure() -> None:
    from epcsaft import ConvergenceError, EosError
    from tests.native_inputs import ETHANE, parameters

    model = Pure(parameters((ETHANE,)))
    anchor = model.state(200.0, rho=1.0e-15)
    with pytest.raises(ConvergenceError, match="left_domain: residual 0") as failure:
        model.state(200.0, P=1.0e5, phase="liquid", anchor=anchor)
    assert isinstance(failure.value, EosError)


# #566 derived State properties. Clapeyron BasicIdeal as one constant polynomial: Cp = 5R/2, anchors at 300 K, 1 bar.
_R = 8.31446261815324
_BASIC_IDEAL = IdealPolynomial([2.5 * _R], 300.0, 2.5 * _R * 300.0, _R * (2.5 + 2.5 * math.log(300.0) + math.log(_R / 1e5)))
_DERIVED = ("isothermal_compressibility", "isobaric_expansivity", "heat_capacity_ratio", "isentropic_compressibility",
            "speed_of_sound", "joule_thomson_coefficient")
_PARTIAL = ("partial_molar_volume", "partial_molar_enthalpy", "partial_molar_entropy")
# Clapeyron 0.6.29 VT_* and _partial_property at 300 K, the same parameters, masses, R and BasicIdeal; columns:
# P, Cp, Cv, kappa_T, alpha_P, gamma, kappa_S, w, mu_JT, then v_i, h_i, s_i per component. Retained comparator:
# #566 reference-values.json (SHA-256 25ca53f55c3646826390b4ad5a26c4377e89d901ee938994fac17ea5737c4ded).
_CLAPEYRON: dict[tuple[str, float, tuple[float, ...]], tuple[tuple[float, ...], tuple[tuple[float, ...], ...]]] = {
    ("propane", 10.0, (1.0,)): ((
        24851.063870309004, 20.92484505730069, 12.499257735892634, 4.03896269160768e-05, 0.0033680178565465794,
        1.674087013760289, 2.4126360567934118e-05, 306.5875359596401, 4.972728321514283e-05), (
        (0.1,), (6210.066348330103,), (72.75285746093314,))),
    ("propane", 12000.0, (1.0,)): ((
        17281151.61786351, 57.077344774420204, 25.826492475130674, 3.465815537595254e-09, 0.0020814388241019435,
        2.210030836722493, 1.5682204429034611e-09, 1097.757021369006, -5.48332492593646e-07), (
        (8.333333333333326e-05,), (-9596.036203472815,), (-14.042556200231374,))),
    ("methane-ethane", 200.0, (0.4, 0.6)): ((
        487070.6864253928, 21.778686538246202, 12.639631605845253, 2.103406943133779e-06, 0.0035798651369842753,
        1.723047571115487, 1.2207480387625344e-06, 409.2143553160764, 1.697979833756275e-05), (
        (0.005090799699017654, 0.004939466867321562), (6177.358948615365, 5979.632802908508),
        (55.54237630672732, 51.75357688601219))),
    ("methane-ethane", 8000.0, (0.4, 0.6)): ((
        9520487.906309793, 95.96140244639435, 16.194940470527886, 1.34988005329227e-07, 0.016945021368564746,
        5.925393959985727, 2.2781270956969778e-08, 473.63622130238076, 5.319204266593722e-06), (
        (0.00025345238566291625, 3.9365076224722084e-05), (7989.34478068556, -3789.1227552397054),
        (36.764362090305546, 1.4195930772730305))),
}


def _basic_model(name: str) -> Mixture:
    components = ("propane",) if name == "propane" else ("methane", "ethane")
    parameters = Parameters.from_json(parameters_json(f"gross-2001-{name}", 3), components=components)
    record = ThermochemistryRecord(1e5, [[IdealInterval(100.0, 1000.0, True, True, IdealCorrelation(_BASIC_IDEAL, 1e5))]
                                         for _ in components])
    return (Pure if len(components) == 1 else Mixture)(parameters, thermochemistry=record)


def test_derived_properties_match_clapeyron(evidence: Any) -> None:
    for (name, rho, x), (row, partials) in _CLAPEYRON.items():
        model = _basic_model(name)
        state = model.state(300.0, rho=rho, x=x, derived_properties=True)
        inputs = {"T_K": 300.0, "rho_mol_m3": rho, "x": list(x), "packet": f"gross-2001-{name}/3"}
        scalars = ("pressure", "total_isobaric_heat_capacity", "total_isochoric_heat_capacity", *_DERIVED)
        for key, reference in zip(scalars, row):
            same(evidence, f"derived/{name}/{rho}/{key}", getattr(state, key), reference, 0.0, 1e-10, inputs=inputs)
        for key, components in zip(_PARTIAL, partials):
            for i, value in enumerate(getattr(state, key)):
                same(evidence, f"derived/{name}/{rho}/{key}/{i}", value, components[i], 0.0, 1e-10, inputs=inputs)
        # Weighted partial molar sums recover the molar totals: v, h and s on the State's own standard pressure.
        for key, total in zip(_PARTIAL, (1.0 / state.molar_density, state.total_enthalpy, state.total_entropy)):
            same(evidence, f"derived/{name}/{rho}/sum/{key}", sum(a * b for a, b in zip(x, getattr(state, key))),
                 total, 0.0, 1e-12, inputs=inputs)
        # The pressure-closed State assembles the same values at its converged density.
        closed = model.state(300.0, P=row[0], x=x, phase="liquid" if rho > 1000.0 else "vapor", derived_properties=True)
        at_root = model.state(300.0, rho=closed.molar_density, x=x, derived_properties=True)
        for key in (*_DERIVED, *_PARTIAL):
            assert getattr(closed, key) == getattr(at_root, key), key


def test_sound_speed_reaches_the_ideal_gas_limit(evidence: Any) -> None:
    for name, x in (("propane", (1.0,)), ("methane-ethane", (0.4, 0.6))):
        model = _basic_model(name)
        mass = sum(a * b for a, b in zip(x, model.parameters.molar_masses))
        ideal = math.sqrt(5.0 / 3.0 * _R * 300.0 / mass)  # gamma_id = Cp_id/(Cp_id - R) for Cp_id = 5R/2
        speeds = [model.state(300.0, rho=rho, x=x, derived_properties=True).speed_of_sound for rho in (1e-2, 1e-4, 1e-6, 1e-8)]
        errors = [abs(speed - ideal) / ideal for speed in speeds if speed is not None]
        assert len(errors) == 4 and errors == sorted(errors, reverse=True), errors  # dilution error falls with density
        same(evidence, f"sound_limit/{name}", ideal * (1.0 + errors[-1]), ideal, 0.0, 1e-10,
             inputs={"T_K": 300.0, "rho_mol_m3": 1e-8, "x": list(x), "ladder_relative_errors": errors})


def test_derived_property_availability() -> None:
    parameters = Parameters.from_json(parameters_json("gross-2001-methane-ethane", 3), components=("methane", "ethane"))
    interval = [IdealInterval(100.0, 1000.0, True, True, IdealCorrelation(_BASIC_IDEAL, 1e5))]
    late = [IdealInterval(350.0, 1000.0, True, True, IdealCorrelation(IdealPolynomial([2.5 * _R], 500.0, 0.0, 0.0), 1e5))]
    assert Mixture(parameters).state(300.0, rho=8000.0, x=(0.4, 0.6)).partial_molar_volume is None  # not requested
    for record in (None, ThermochemistryRecord(1e5, [interval, []]), ThermochemistryRecord(1e5, [interval, late])):
        state = Mixture(parameters, thermochemistry=record).state(300.0, rho=8000.0, x=(0.4, 0.6), derived_properties=True)
        assert state.isothermal_compressibility is not None and state.isobaric_expansivity is not None
        assert all(value is not None for value in state.partial_molar_volume or [None])
        assert state.heat_capacity_ratio is state.speed_of_sound is state.joule_thomson_coefficient is None
        assert state.partial_molar_enthalpy == state.partial_molar_entropy == [None, None]
    # A zero-fraction neutral species keeps finite insertion volume and enthalpy; its entropy diverges.
    zero = _basic_model("methane-ethane").state(300.0, rho=200.0, x=(0.0, 1.0), derived_properties=True)
    assert zero.partial_molar_volume is not None and zero.partial_molar_enthalpy is not None
    assert zero.partial_molar_entropy is not None
    assert None not in (zero.partial_molar_volume[0], zero.partial_molar_enthalpy[0], zero.partial_molar_entropy[1])
    assert zero.partial_molar_entropy[0] is None
    # A nonpositive (dP/drho)_T keeps every old field and makes every new response value unavailable.
    unstable = _basic_model("propane").state(300.0, rho=5000.0, derived_properties=True)
    assert unstable.pressure_density_derivative < 0.0 and unstable.total_enthalpy is not None
    assert all(getattr(unstable, key) is None for key in _DERIVED)
    assert unstable.partial_molar_volume == unstable.partial_molar_enthalpy == [None]
    # Associated neutral liquid with absent charged slots: neutral volumes finite, IonicInsertion entries None.
    mea = Mixture(Parameters.from_mapping(_mea_mapping(parameters_json("mea-co2-h2o-nine-species-estimation", 4)), components=_MEA_ORDER))
    liquid = mea.state(313.15, P=1e5, x=(0.85, 0.10, 0.05, 0, 0, 0, 0, 0, 0), phase="liquid", derived_properties=True)
    assert liquid.partial_molar_volume is not None
    assert all(value is not None for value in liquid.partial_molar_volume[:3])
    assert liquid.partial_molar_volume[3:] == [None] * 6


def test_contributions_sum_to_residual_helmholtz() -> None:
    # Acetone at 400 K, 10000 mol/m3 is the failing dipolar split (a negative-pressure algebraic check only).
    acetone = Pure(Parameters.from_json(parameters_json("gross-vrabec-2006", 1), components=("acetone",))).state(400.0, rho=1e4)
    mea = Mixture(Parameters.from_mapping(_mea_mapping(parameters_json("mea-co2-h2o-nine-species-estimation", 4)), components=_MEA_ORDER))
    nonpolar = Mixture(Parameters.from_json(parameters_json("gross-2001-methane-ethane", 3))).state(300.0, rho=8000.0, x=(0.4, 0.6))
    ionic = mea.state(313.15, rho=3.0e4, x=(0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0))
    associated = mea.state(313.15, P=1e5, x=(0.85, 0.10, 0.05, 0, 0, 0, 0, 0, 0), phase="liquid")
    assert acetone.dipolar < 0.0 and nonpolar.dipolar == 0.0
    for state in (acetone, nonpolar, ionic, associated):
        terms = (state.hard_chain, state.dispersion, state.association, state.debye_huckel, state.born, state.dipolar)
        close(sum(terms), state.residual_helmholtz, 1e-12)

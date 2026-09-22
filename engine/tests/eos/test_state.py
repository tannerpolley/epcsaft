import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import pytest

from epcsaft import DensityBranch, Mixture, ParameterError, Parameters, Pure
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


def _parameter_metadata(path: Path, packet_id: str, version: int, components: tuple[str, ...]) -> dict[str, Any]:
    supplied_mapping = json.loads(path.read_text())
    return {
        "packet_id": packet_id,
        "packet_version": version,
        "parameter_file": {
            "path": str(path.relative_to(REPOSITORY_ROOT)),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "size_bytes": path.stat().st_size,
        },
        "effective_parameters": {
            "supplied_mapping": supplied_mapping,
            "selected_components": list(components),
        },
    }


def close(left: float, right: float, relative: float) -> None:
    assert math.isclose(left, right, rel_tol=relative, abs_tol=relative)


SYMMETRIC_STATE_CRITERION = {
    "type": "scalar",
    "rule": "symmetric_close",
    "atol": 1e-12,
    "rtol": 1e-12,
    "scale": "max(abs(model), abs(reference))",
    "combining_rule": "max(atol, rtol*scale)",
}


def test_relabelling_and_state_identities(evidence: Any) -> None:
    packet = parameters_json("mea-co2-h2o-nine-species-estimation", 3)
    order = _MEA_ORDER
    reverse = tuple(reversed(order))
    x = (0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0)
    temperature = 313.15
    density = 3.0e4
    pressure_target = 1.0e5
    phase = "liquid"
    first = Mixture(Parameters.from_json(packet, components=order)).state(
        temperature, rho=density, x=x
    )
    second = Mixture(Parameters.from_json(packet, components=reverse)).state(
        temperature, rho=density, x=tuple(reversed(x))
    )
    first_fugacity = first.log_fugacity_coefficient
    second_fugacity = second.log_fugacity_coefficient
    fugacity_available = first_fugacity is not None and second_fugacity is not None
    mixture = Mixture(Parameters.from_json(packet, components=order))
    state = mixture.state(temperature, P=pressure_target, x=x, phase=phase)
    anchored = mixture.state(
        temperature, P=pressure_target, x=x, phase=phase, anchor=state
    )
    parameter_metadata = _parameter_metadata(
        packet,
        "mea-co2-h2o-nine-species-estimation",
        3,
        order,
    )
    parameter_metadata["effective_parameters"] = {
        **parameter_metadata["effective_parameters"],
        "relabelled_selected_components": list(reverse),
    }
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
    criterion = dict(SYMMETRIC_STATE_CRITERION)
    for key, quantity, left, right, latex in (
        ("pressure", "pressure", first.pressure, second.pressure, "P"),
        ("residual_helmholtz", "residual Helmholtz energy", first.residual_helmholtz, second.residual_helmholtz, "a^{res}"),
        ("compressibility_factor", "compressibility factor", first.compressibility_factor, second.compressibility_factor, "Z"),
        ("relative_permittivity", "bulk relative permittivity", first.bulk_relative_permittivity, second.bulk_relative_permittivity, "\\varepsilon_r"),
    ):
        evidence.compare(
            f"state_identity/{key}", quantity=quantity, model=float(left), reference=float(right),
            unit="1" if key not in {"pressure"} else "Pa", inputs=inputs,
            criterion=criterion, basis="Relationship", method="component relabelling",
            claim="The scalar state quantity is invariant under component relabelling.",
            relationship_latex=fr"{latex}(x_1,\ldots,x_n) = {latex}(x_n,\ldots,x_1)",
            quantity_latex=latex,
            inputs_latex=r"T,\rho,\mathbf{x},\mathrm{components}",
        )
    for index, component in enumerate(order):
        reverse_index = len(order) - index - 1
        evidence.compare(
            f"state_identity/residual_chemical_potential/{index}",
            quantity="residual chemical potential over RT",
            model=float(first.residual_chemical_potential_over_rt[index]),
            reference=float(second.residual_chemical_potential_over_rt[reverse_index]),
            unit="1",
            inputs={**inputs, "component": component, "relabelled_component": reverse[reverse_index]},
            coordinates={"component_index": index, "relabelled_index": reverse_index},
            criterion=criterion,
            basis="Relationship",
            method="component relabelling",
            claim="The residual chemical potential vector reverses with the component labels.",
            comparison_group_id="state_identity/residual_chemical_potential",
            point_id=component,
            relationship_latex=r"\mu_{i,\mathrm{first}}^{res}/RT = \mu_{n-1-i,\mathrm{relabelled}}^{res}/RT",
            quantity_latex=r"\mu_i^{res}/RT",
            inputs_latex=r"T,\rho,\mathbf{x},\mathrm{components}",
        )
        if first_fugacity is not None and second_fugacity is not None:
            evidence.compare(
                f"state_identity/log_fugacity_coefficient/{index}",
                quantity="log fugacity coefficient",
                model=float(first_fugacity[index]),
                reference=float(second_fugacity[reverse_index]),
                unit="1",
                inputs={**inputs, "component": component, "relabelled_component": reverse[reverse_index]},
                coordinates={"component_index": index, "relabelled_index": reverse_index},
                criterion=criterion,
                basis="Relationship",
                method="component relabelling",
                claim="The log fugacity coefficient vector reverses with the component labels.",
                comparison_group_id="state_identity/log_fugacity_coefficient",
                point_id=component,
                relationship_latex=r"\ln\phi_{i,\mathrm{first}} = \ln\phi_{n-1-i,\mathrm{relabelled}}",
                quantity_latex=r"\ln\phi_i",
                inputs_latex=r"T,\rho,\mathbf{x},\mathrm{components}",
            )
    euler_left = sum(a * b for a, b in zip(x, first.residual_chemical_potential_over_rt))
    euler_right = first.residual_helmholtz + first.compressibility_factor - 1.0
    evidence.compare(
        "state_identity/euler_identity",
        quantity="Euler residual identity",
        model=float(euler_left),
        reference=float(euler_right),
        unit="1",
        inputs={**inputs, "identity": "sum(x_i * mu_i^res/RT) = a^res + Z - 1"},
        criterion=criterion,
        basis="Relationship",
        method="residual Euler identity",
        claim="The composition-weighted residual chemical potential identity equals a^res + Z - 1.",
        relationship_latex=r"\sum_i x_i\mu_i^{res}/RT = a^{res} + Z - 1",
        quantity_latex=r"\sum_i x_i\mu_i^{res}/RT",
        inputs_latex=r"T,\rho,\mathbf{x}",
    )
    gas_constant = 8.31446261815324
    codata_pressure = first.compressibility_factor * first.molar_density * gas_constant * temperature
    evidence.compare(
        "state_identity/codata_pressure",
        quantity="pressure from CODATA gas constant",
        model=float(first.pressure),
        reference=float(codata_pressure),
        unit="Pa",
        inputs={**inputs, "gas_constant_J_mol_K": gas_constant, "gas_constant_basis": "CODATA 2018 exact value"},
        criterion={
            "type": "scalar",
            "rule": "symmetric_close",
            "atol": 0.0,
            "rtol": 1e-12,
            "scale": "max(abs(model), abs(reference))",
            "combining_rule": "max(atol, rtol*scale)",
        },
        basis="Relationship",
        method="pressure equation of state identity",
        claim="Pressure equals Z times molar density, the CODATA gas constant, and temperature.",
        relationship_latex=r"P = Z\rho RT",
        quantity_latex="P",
        inputs_latex=r"T,\rho,Z,R",
    )
    evidence.compare(
        "state_identity/pressure_anchor",
        quantity="pressure anchor",
        model=float(state.pressure),
        reference=pressure_target,
        unit="Pa",
        inputs={**inputs, "resolved_density_mol_m3": float(state.molar_density)},
        criterion={
            "type": "scalar",
            "rule": "symmetric_close",
            "atol": 1e-9,
            "rtol": 1e-9,
            "scale": "max(abs(model), abs(reference))",
            "combining_rule": "max(atol, rtol*scale)",
        },
        basis="Relationship",
        method="fixed-pressure state solve",
        claim="The pressure-specified state returns the requested pressure anchor.",
        relationship_latex=r"P(T,P_{target},x) = P_{target}",
        quantity_latex="P",
        inputs_latex=r"T,P_{target},\mathbf{x},\mathrm{phase}",
    )
    evidence.compare(
        "state_identity/anchored_density",
        quantity="anchored molar density",
        model=float(anchored.molar_density),
        reference=float(state.molar_density),
        unit="mol/m^3",
        inputs={
            **inputs,
            "anchor": "state",
            "resolved_density_mol_m3": float(anchored.molar_density),
            "anchor_density_mol_m3": float(state.molar_density),
        },
        criterion={
            "type": "scalar",
            "rule": "symmetric_close",
            "atol": 1e-10,
            "rtol": 1e-10,
            "scale": "max(abs(model), abs(reference))",
            "combining_rule": "max(atol, rtol*scale)",
        },
        basis="Relationship",
        method="anchored fixed-pressure state solve",
        claim="Reusing the converged state as an anchor preserves its molar density.",
        relationship_latex=r"\rho(T,P,x;\mathrm{anchor}=s) = \rho(s)",
        quantity_latex=r"\rho",
        inputs_latex=r"T,P,\mathbf{x},\mathrm{phase},\mathrm{anchor}",
    )
    for field in ("pressure", "residual_helmholtz", "compressibility_factor", "bulk_relative_permittivity"):
        close(getattr(first, field), getattr(second, field), 1e-12)
    assert (first.log_fugacity_coefficient is None) == (second.log_fugacity_coefficient is None)
    for index in range(len(order)):
        reverse_index = len(order) - index - 1
        close(
            first.residual_chemical_potential_over_rt[index],
            second.residual_chemical_potential_over_rt[reverse_index],
            1e-12,
        )
        if first_fugacity is not None and second_fugacity is not None:
            close(
                first_fugacity[index],
                second_fugacity[reverse_index],
                1e-12,
            )
    assert math.isclose(
        sum(a * b for a, b in zip(x, first.residual_chemical_potential_over_rt)),
        first.residual_helmholtz + first.compressibility_factor - 1.0,
        rel_tol=1e-12,
        abs_tol=1e-12,
    )
    assert math.isclose(
        first.pressure,
        first.compressibility_factor * first.molar_density * 8.31446261815324 * temperature,
        rel_tol=1e-12,
    )
    close(state.pressure, pressure_target, 1e-9)
    close(anchored.molar_density, state.molar_density, 1e-10)


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
    parameter_metadata = _parameter_metadata(
        packet,
        "gross-2001-methane-ethane",
        2,
        mixture_components,
    )
    parameter_metadata["effective_parameters"] = {
        **parameter_metadata["effective_parameters"],
        "pure_components": list(pure_components),
        "mixture_components": list(mixture_components),
    }
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
                "_parameters": _parameter_metadata(packet, "gross-2001-methane-ethane", 3, ("ethane",))},
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
    mapping["components"][0]["fixed"]["molar_mass"]["value"]["unit"] = "gram"
    with pytest.raises(ParameterError, match="(?=.*component/ethane/molar_mass)(?=.*gram)"):
        Parameters.from_mapping(mapping, components=("ethane",))
    parameters = Parameters.from_json(parameters_json("gross-2001-methane-ethane", 3), components=("methane", "ethane"))
    with pytest.raises(ParameterError):
        Mixture(parameters).state(300.0, P=1e5, x=(0.5, 0.5), phase="gas")



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


def test_association_topology_is_order_independent_and_typed() -> None:
    mapping = json.loads(
        parameters_json("khudaida-2026-figure-2-electrolyte-lle", 6).read_text()
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

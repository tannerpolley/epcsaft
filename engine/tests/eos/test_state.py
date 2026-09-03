import copy
import json
import math

import pytest

from epcsaft import Mixture, ParameterError, Parameters, Pure

from tests.packets import parameters_json


def close(left: float, right: float, relative: float) -> None:
    assert math.isclose(left, right, rel_tol=relative, abs_tol=relative)


def test_relabelling_and_state_identities() -> None:
    packet = parameters_json("mea-co2-h2o-nine-species-estimation", 3)
    order = ("water", "monoethanolamine", "carbon-dioxide", "protonated-monoethanolamine", "carbamate-anion", "bicarbonate-anion", "carbonate-anion", "hydronium-cation", "hydroxide-anion")
    reverse = tuple(reversed(order))
    x = (0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0)
    first = Mixture(Parameters.from_json(packet, components=order)).state(313.15, rho=3.0e4, x=x)
    second = Mixture(Parameters.from_json(packet, components=reverse)).state(313.15, rho=3.0e4, x=tuple(reversed(x)))
    for field in ("pressure", "residual_helmholtz", "compressibility_factor", "bulk_relative_permittivity"):
        close(getattr(first, field), getattr(second, field), 1e-12)
    assert (first.log_fugacity_coefficient is None) == (second.log_fugacity_coefficient is None)
    for index in range(len(order)):
        reverse_index = len(order) - index - 1
        close(first.residual_chemical_potential_over_rt[index], second.residual_chemical_potential_over_rt[reverse_index], 1e-12)
        if first.log_fugacity_coefficient is not None and second.log_fugacity_coefficient is not None:
            close(first.log_fugacity_coefficient[index], second.log_fugacity_coefficient[reverse_index], 1e-12)
    assert math.isclose(sum(a * b for a, b in zip(x, first.residual_chemical_potential_over_rt)), first.residual_helmholtz + first.compressibility_factor - 1.0, rel_tol=1e-12, abs_tol=1e-12)
    # CODATA 2018 exact gas constant.
    assert math.isclose(first.pressure, first.compressibility_factor * first.molar_density * 8.31446261815324 * 313.15, rel_tol=1e-12)
    mixture = Mixture(Parameters.from_json(packet, components=order))
    state = mixture.state(313.15, P=1e5, x=x, phase="liquid")
    close(state.pressure, 1e5, 1e-9)
    anchored = mixture.state(313.15, P=1e5, x=x, phase="liquid", anchor=state)
    close(anchored.molar_density, state.molar_density, 1e-10)


def test_pure_is_one_component_mixture() -> None:
    packet = parameters_json("gross-2001-methane-ethane", 2)
    pure = Pure(Parameters.from_json(packet, components=("ethane",))).state(300.0, rho=1.0e4)
    mixture = Mixture(Parameters.from_json(packet, components=("methane", "ethane"))).state(300.0, rho=1.0e4, x=(0.0, 1.0))
    close(pure.pressure, mixture.pressure, 1e-12)


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
    mapping = copy.deepcopy(json.loads(parameters_json("gross-2001-methane-ethane", 2).read_text()))
    mapping["components"][0]["fixed"]["molar_mass"]["value"]["unit"] = "gram"
    with pytest.raises(ParameterError, match="(?=.*component/ethane/molar_mass)(?=.*gram)"):
        Parameters.from_mapping(mapping, components=("ethane",))
    parameters = Parameters.from_json(parameters_json("gross-2001-methane-ethane", 2), components=("methane", "ethane"))
    with pytest.raises(ParameterError):
        Mixture(parameters).state(300.0, P=1e5, x=(0.5, 0.5), phase="gas")

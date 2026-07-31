from __future__ import annotations

from copy import deepcopy

import pytest

from epcsaft import Mixture, ParameterError, Parameters, unit_registry as u

from tests._raw_cases import RAW_CASES


def _neutral_binary(*, descriptive: bool = False) -> dict[str, object]:
    count = "segment_count" if descriptive else "m"
    diameter = "segment_diameter" if descriptive else "s"
    energy = "dispersion_energy_over_k" if descriptive else "e"
    interaction = "dispersion_interaction" if descriptive else "k_ij"
    return {
        "schema": "epcsaft.parameters",
        "schema_version": 1,
        "components": ["component-a", "component-b"],
        "parameters": {
            count: [1.5, 2.25],
            diameter: [3.4, 3.75],
            energy: [180.0, 245.0],
            interaction: [[0.0, 0.03], [0.03, 0.0]],
        },
        "options": {"relative_permittivity_formulation": "none"},
        "validity": {"kind": "unknown"},
    }


def test_from_dictionary__normalizes_aliases_and_copies_inputs() -> None:
    abbreviated = _neutral_binary()
    descriptive = _neutral_binary(descriptive=True)

    left = Parameters.from_dictionary(abbreviated)
    right = Parameters.from_dictionary(descriptive)
    abbreviated["components"][0] = "mutated"  # type: ignore[index]
    abbreviated["parameters"]["m"][0] = 999.0  # type: ignore[index]

    assert left.component_ids == ("component-a", "component-b")
    assert left.fingerprint == right.fingerprint
    state = Mixture(left).state(
        T=300.0 * u.kelvin,
        rho=1000.0 * u.mole / u.meter**3,
        x=(0.4, 0.6),
    )
    assert state.pressure.to("pascal").magnitude > 0.0


def test_from_dictionary__round_trips_through_an_external_bundle(tmp_path) -> None:
    expected = Parameters.from_dictionary(_neutral_binary())
    path = tmp_path / "bundle"

    expected.to_bundle(path)
    observed = Parameters.from_bundle(path)

    assert observed.component_ids == expected.component_ids
    assert observed.fingerprint == expected.fingerprint
    assert observed.records == expected.records


@pytest.mark.parametrize(
    ("abbreviated", "descriptive", "value"),
    (
        ("mw", "molar_mass", (0.02, 0.03)),
        ("m", "segment_count", (1.5, 2.25)),
        ("s", "segment_diameter", (3.4, 3.75)),
        ("e", "dispersion_energy_over_k", (180.0, 245.0)),
        ("z", "charge_number", (0, 0)),
        ("epsilon_r", "relative_permittivity", (2.0, 3.0)),
        ("d_born", "born_diameter", (3.0, 4.0)),
        ("f_solv", "solvation_factor", (1.0, 1.0)),
        ("v_dielectric", "dielectric_volume", (1.0e-5, 2.0e-5)),
        ("t_dielectric", "dielectric_temperature", (200.0, 250.0)),
        ("a_ion", "ion_specific_suppression_coefficient", (1.0, 2.0)),
        ("mu2", "dipole_scaling", (1.0, 2.0)),
        ("alpha", "polarizability_scaling", (1.0, 2.0)),
        ("a2", "correlation_integral_parameter", (1.0, 2.0)),
        ("k_ij", "dispersion_interaction", ((0.0, 0.02), (0.02, 0.0))),
        ("l_ij", "diameter_interaction", ((0.0, 0.01), (0.01, 0.0))),
    ),
)
def test_from_dictionary__all_parameter_aliases_are_equivalent(
    abbreviated: str, descriptive: str, value: object
) -> None:
    left = _neutral_binary()
    right = _neutral_binary()
    if abbreviated in left["parameters"]:
        del left["parameters"][abbreviated]
        del right["parameters"][abbreviated]
    left["parameters"][abbreviated] = value
    right["parameters"][descriptive] = value

    assert (
        Parameters.from_dictionary(left).fingerprint
        == Parameters.from_dictionary(right).fingerprint
    )


@pytest.mark.parametrize(
    ("abbreviated", "descriptive", "value"),
    (
        ("electrolyte", "electrolyte_formulation", "fully-dissociated-debye-huckel"),
        ("permittivity_model", "relative_permittivity_formulation", "none"),
        ("epsilon_r_ion", "ionic_region_relative_permittivity", 8.0),
        ("a_dh", "ion_fraction_suppression_coefficient", 7.01),
    ),
)
def test_from_dictionary__all_option_aliases_are_equivalent(
    abbreviated: str, descriptive: str, value: object
) -> None:
    left = _neutral_binary()
    right = _neutral_binary()
    left["options"] = {abbreviated: value, "permittivity_model": "none"}
    right["options"] = {descriptive: value, "permittivity_model": "none"}
    if abbreviated == "permittivity_model":
        left["options"] = {abbreviated: value}
        right["options"] = {descriptive: value}

    assert (
        Parameters.from_dictionary(left).fingerprint
        == Parameters.from_dictionary(right).fingerprint
    )


def test_from_dictionary__interprets_fixed_units() -> None:
    parameters = Parameters.from_dictionary(_neutral_binary())
    records = {
        (record.component_id, record.family): record
        for record in parameters.records
        if hasattr(record, "component_id") and hasattr(record, "value")
    }

    assert float(records[("component-a", "segment_diameter")].value.to("angstrom").magnitude) == 3.4
    assert (
        float(records[("component-a", "dispersion_energy_over_k")].value.to("kelvin").magnitude)
        == 180.0
    )


@pytest.mark.parametrize(("name", "case"), RAW_CASES.items(), ids=RAW_CASES)
def test_from_dictionary__resolves_each_supported_topology(name: str, case) -> None:
    model = Mixture(Parameters.from_dictionary(case))
    fractions = (
        (0.98, 0.01, 0.01)
        if "electrolyte" in name
        else tuple(1.0 / len(model.component_ids) for _ in model.component_ids)
    )

    state = model.state(
        T=298.15 * u.kelvin,
        rho=100.0 * u.mole / u.meter**3,
        x=fractions,
    )

    assert state.pressure.to("pascal").magnitude > 0.0


@pytest.mark.parametrize(
    ("mutate", "code"),
    [
        (
            lambda value: value["parameters"].update(segment_count=[1.5, 2.25]),
            "duplicate-alias",
        ),
        (lambda value: value.update(schema_version=2), "invalid-schema"),
        (lambda value: value["parameters"].update(unknown=[1.0, 2.0]), "invalid-schema"),
        (lambda value: value["parameters"].update(m=[1.5]), "invalid-shape"),
        (lambda value: value["parameters"].update(e=[180.0, float("nan")]), "invalid-value"),
        (
            lambda value: value["parameters"].update(k_ij=[[0.0, 0.1], [0.2, 0.0]]),
            "invalid-value",
        ),
        (
            lambda value: value["parameters"].update(k_ij=[[0.1, 0.0], [0.0, 0.0]]),
            "invalid-value",
        ),
        (
            lambda value: value["options"].update(permittivity_model="none"),
            "duplicate-alias",
        ),
        (lambda value: value.update(validity={}), "invalid-domain"),
    ],
)
def test_from_dictionary__fails_closed(mutate, code: str) -> None:
    value = deepcopy(_neutral_binary())
    mutate(value)

    with pytest.raises(ParameterError) as caught:
        Parameters.from_dictionary(value)

    assert caught.value.code == code

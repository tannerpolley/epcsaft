from __future__ import annotations

from copy import deepcopy

import pytest

from epcsaft import EosError, Mixture, ParameterError, Parameters, Pure
from epcsaft import unit_registry as u
from tests.eos._cases import AQUEOUS_ELECTROLYTE, NEUTRAL_BINARY, neutral_pure


def test_state_units_and_nonfinite_values_fail_closed() -> None:
    model = Pure(Parameters.from_mapping(neutral_pure("propane")))

    with pytest.raises(EosError) as invalid_unit:
        model.state(T=300.0, rho=1000.0 * u.mole / u.meter**3)
    assert invalid_unit.value.code == "invalid-state"
    with pytest.raises(EosError) as nonfinite:
        model.state(T=float("nan") * u.kelvin, rho=1000.0 * u.mole / u.meter**3)
    assert nonfinite.value.code == "invalid-state"


@pytest.mark.parametrize(
    "composition",
    ((0.5,), (0.5, 0.6), (-0.1, 1.1), (float("nan"), 0.0)),
    ids=("count", "normalization", "negative", "nonfinite"),
)
def test_malformed_compositions_fail_closed(composition) -> None:
    model = Mixture(Parameters.from_mapping(NEUTRAL_BINARY))

    with pytest.raises(EosError) as captured:
        model.state(
            T=300.0 * u.kelvin,
            rho=1.0 * u.mole / u.meter**3,
            x=composition,
        )
    assert captured.value.code == "invalid-state"


def test_native_model_enforces_resolved_eos_domain() -> None:
    native = Mixture(Parameters.from_mapping(AQUEOUS_ELECTROLYTE))._native

    with pytest.raises(ValueError):
        native.evaluate(298.14, 1000.0, (0.98, 0.01, 0.01))
    with pytest.raises(ValueError):
        native.evaluate_pressure(298.15, 100001.0, (0.98, 0.01, 0.01), "liquid")
    with pytest.raises(ValueError):
        native.evaluate(298.15, 1000.0, (0.98, 0.02, 0.0))
    with pytest.raises(ValueError):
        native.evaluate(298.15, 1000.0, (0.5, 0.25, 0.25))
    with pytest.raises(ValueError):
        native.evaluate_phase_derivatives(
            298.14, (0.98, 0.01, 0.01), 1.0e-3, (), 2
        )
    with pytest.raises(ValueError):
        native.evaluate_fixed_pressure_composition_derivatives(
            298.15, (0.98, 0.02, 0.0), 1.0e-3
        )


def test_impossible_density_request_fails_closed() -> None:
    model = Pure(Parameters.from_mapping(neutral_pure("propane")))

    with pytest.raises(EosError):
        model.state(T=300.0 * u.kelvin, rho=1.0e9 * u.mole / u.meter**3)


def test_invalid_charge_and_permittivity_declarations_fail_closed() -> None:
    charge = deepcopy(AQUEOUS_ELECTROLYTE)
    ion = next(
        component
        for component in charge["components"]
        if component["component_id"] == "sodium-cation"
    )
    ion["fixed"]["charge_number"]["value"]["magnitude"] = 0.5
    with pytest.raises((ParameterError, TypeError, ValueError)):
        Parameters.from_mapping(charge)

    permittivity = deepcopy(AQUEOUS_ELECTROLYTE)
    correlation = next(
        record
        for record in permittivity["correlations"]
        if record["family"] == "relative_permittivity"
    )
    correlation["constant"]["value"]["magnitude"] = 0.0
    with pytest.raises((ParameterError, EosError, ValueError)):
        Mixture(Parameters.from_mapping(permittivity))


def test_candidate_model_requires_every_participating_parameter_record() -> None:
    mapping = deepcopy(NEUTRAL_BINARY)
    mapping["domains"].append(
        {
            "domain_id": "engine-candidate-domain",
            "kind": "fit-range",
            "temperature_min": {"magnitude": 290.0, "unit": "kelvin"},
            "temperature_max": {"magnitude": 400.0, "unit": "kelvin"},
        }
    )
    candidate = {
        "candidate_domain_id": "engine-candidate-domain",
        "qualification": "candidate_extrapolation",
        "source_sha256": f"sha256:{'a' * 64}",
    }
    for family in mapping["model_families"]:
        family.update(candidate)
    records = [
        *(record for component in mapping["components"] for record in component["fixed"].values()),
        *(record for component in mapping["components"] for record in component["coefficients"]),
        *(record for pair in mapping["pairs"] for record in pair["coefficients"]),
    ]
    for record in records:
        record.update(candidate)
    Parameters.from_mapping(mapping)

    for field in candidate:
        records[-1].pop(field)
    with pytest.raises(ParameterError) as failure:
        Parameters.from_mapping(mapping)
    assert failure.value.code == "missing-record"

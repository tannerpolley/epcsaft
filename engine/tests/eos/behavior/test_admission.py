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


def test_impossible_density_and_pressure_requests_fail_closed() -> None:
    model = Pure(Parameters.from_mapping(neutral_pure("propane")))

    with pytest.raises(EosError):
        model.state(T=300.0 * u.kelvin, rho=1.0e9 * u.mole / u.meter**3)
    with pytest.raises(EosError):
        model.state(T=300.0 * u.kelvin, P=1.0e12 * u.pascal, phase="vapor")


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

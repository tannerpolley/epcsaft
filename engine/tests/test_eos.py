from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest

import epcsaft
from epcsaft import ParameterBundle, ParameterError
from epcsaft.eos import EPCSAFT


def test_native_model_boundary_is_ordered_strict_and_immutable() -> None:
    bundle = ParameterBundle.from_catalog(
        "gross-2001-methane-ethane", version=1
    )
    parameters = bundle.select(("methane", "ethane"))
    model = EPCSAFT(parameters)

    assert model.component_ids == ("methane", "ethane")
    assert model.parameter_fingerprint == parameters.fingerprint
    assert EPCSAFT(bundle.select(("ethane", "methane"))).component_ids == (
        "ethane",
        "methane",
    )
    with pytest.raises(FrozenInstanceError):
        model._component_ids = ("ethane", "methane")

    fixture = ParameterBundle.from_path(
        Path(__file__).parent / "data/bundles/figiel-2025-contract/1"
    )
    with pytest.raises(ParameterError) as exc_info:
        EPCSAFT(fixture.select(("water",)))
    assert exc_info.value.code == "invalid-value"
    assert dict(exc_info.value.context) == {
        "record_id": "water-association-energy",
        "family": "association_energy_over_k",
    }

    assert epcsaft.__all__ == (
        "ParameterBundle",
        "ParameterSet",
        "ParameterError",
        "unit_registry",
    )

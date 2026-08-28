from __future__ import annotations

from pathlib import Path

import pytest

import epcsaft
from epcsaft.eos import Pure
from epcsaft.regression import fit
from tests.regression._cases import (
    PURE_TRUTH,
    pure_parameters,
    pure_saturation_fit_case,
    scaled_error,
)


def _saturation_predictions(value: float) -> tuple[float, ...]:
    model = Pure(pure_parameters(value))
    predictions = []
    for temperature in (195.0, 210.0):
        saturation = model.saturation(T=temperature * epcsaft.unit_registry.kelvin)
        predictions.extend(
            (saturation.saturation_pressure_pa, saturation.liquid.molar_density_mol_m3)
        )
    return tuple(predictions)


def test_saturation_owner_preserves_pressure_density_and_centered_public_jacobian(
    tmp_path: Path,
) -> None:
    """Falsify a fit-owned saturation solve or a pressure/density row swap."""

    mapping, _, expected_predictions = pure_saturation_fit_case(tmp_path)
    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None
    assert accepted.final_physical == pytest.approx((PURE_TRUTH,), abs=1.0e-7)
    assert accepted.predictions == pytest.approx(expected_predictions, rel=2.0e-6, abs=2.0e-8)

    step = 1.0e-3
    below = _saturation_predictions(accepted.final_physical[0] - step)
    above = _saturation_predictions(accepted.final_physical[0] + step)
    centered = tuple((high - low) / (2.0 * step) for low, high in zip(below, above, strict=True))
    owner_derivatives = tuple(
        derivative * prediction / row.multiplier
        for derivative, prediction, row in zip(
            accepted.physical_jacobian,
            accepted.predictions,
            result.specification.rows,
            strict=True,
        )
    )
    assert all(
        scaled_error(actual, expected) <= 5.0e-5
        for actual, expected in zip(owner_derivatives, centered, strict=True)
    )
    assert tuple(row.selector for row in result.evaluated_rows) == (
        "saturation_pressure",
        "saturated_liquid_molar_density",
        "saturation_pressure",
        "saturated_liquid_molar_density",
    )
    assert tuple(row.unit for row in result.evaluated_rows) == (
        "pascal",
        "mole / meter ** 3",
        "pascal",
        "mole / meter ** 3",
    )

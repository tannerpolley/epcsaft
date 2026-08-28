from __future__ import annotations

import json
from pathlib import Path

import pytest

import epcsaft
from epcsaft import equilibrium
from epcsaft.eos import Pure
from epcsaft.regression import fit
from tests.regression._cases import (
    PURE_TRUTH,
    boundary_fit_case,
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


@pytest.mark.slow
def test_boundary_owner_preserves_branch_receipt_and_centered_public_jacobian(
    tmp_path: Path,
) -> None:
    """Falsify GREPE boundary owner transfer without another recovery campaign."""

    mapping, parameters, model, problem, identity, truth = boundary_fit_case(tmp_path)
    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None
    assert accepted.final_physical == pytest.approx((0.08,), abs=3.0e-7)
    assert accepted.predictions == pytest.approx(truth.values, rel=2.0e-6, abs=2.0e-8)

    step = 1.0e-5
    centered_values = []
    for value in (accepted.final_physical[0] - step, accepted.final_physical[0] + step):
        solved = equilibrium.solve(
            model,
            problem,
            active_parameters=epcsaft.ActiveParameterSet(parameters, (identity,), values=(value,)),
        )
        assert solved.status == "evaluated" and solved.values[0] is not None
        centered_values.append(float(solved.values[0]))
    centered = (centered_values[1] - centered_values[0]) / (2.0 * step)
    row = result.specification.rows[0]
    owner_derivative = accepted.physical_jacobian[0] * row.scale / row.multiplier
    assert scaled_error(owner_derivative, centered) <= 5.0e-5

    assert len(accepted.owner_receipts) == 1
    receipt = json.loads(accepted.owner_receipts[0].canonical_bytes_utf8)
    expected_branch = json.loads(json.dumps(truth.branch_signature))
    assert receipt["branch_signature"][0] == expected_branch[0]
    assert receipt["branch_signature"][1][0][0] == expected_branch[1][0][0]
    assert receipt["branch_signature"][1][0][1] == pytest.approx(
        expected_branch[1][0][1], abs=2.0e-12
    )
    assert receipt["branch_signature"][2] == expected_branch[2]
    assert receipt["certificate_scope"] == truth.certificate_scope
    assert receipt["globality_claim"] == truth.globality_claim
    assert receipt["descriptor_fingerprint"].startswith("sha256:")
    assert receipt["problem_fingerprint"].startswith("sha256:")
    metadata = result.specification.rows[0].to_mapping()["series_metadata"]
    assert metadata["aggregate_identity"] == "boundary-output-group"
    assert metadata["covariance_identity"] == "diagonal-whitening-group"

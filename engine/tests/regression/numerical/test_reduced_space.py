from __future__ import annotations

from pathlib import Path

import pytest

from epcsaft.regression import fit
from tests.regression._cases import (
    PAIR_IDENTITY,
    fixed_state_fit_case,
    pressure,
    scaled_error,
    with_parameter,
)


def test_fixed_state_owner_recovers_coordinate_and_centered_public_jacobian(
    tmp_path: Path,
) -> None:
    """Falsify owner-value, active-order, recovery, or residual-Jacobian drift."""

    mapping, input_parameters, states = fixed_state_fit_case(tmp_path)
    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None

    truth = 0.04
    coordinate_scale = 0.05
    recovery_limit = 4.0e-6
    observation_relative_tolerance = 2.0e-6
    # affine_jacobian maps dimensionless coordinate error into the
    # dimensionless scaled-difference residual used by these observations.
    assert (
        max(abs(value) for value in accepted.affine_jacobian) * recovery_limit
        < observation_relative_tolerance
    )
    assert abs(accepted.final_physical[0] - truth) / coordinate_scale <= recovery_limit
    candidate = with_parameter(
        input_parameters, PAIR_IDENTITY, accepted.final_physical[0], "recovered-fixed-state"
    )
    expected_predictions = tuple(
        pressure(candidate, temperature, density, composition)
        for temperature, density, composition, _ in states
    )
    assert accepted.predictions == pytest.approx(expected_predictions, rel=2.0e-6, abs=2.0e-8)

    step = 1.0e-5
    below = with_parameter(
        input_parameters, PAIR_IDENTITY, accepted.final_physical[0] - step, "centered-below"
    )
    above = with_parameter(
        input_parameters, PAIR_IDENTITY, accepted.final_physical[0] + step, "centered-above"
    )
    centered = tuple(
        (
            pressure(above, temperature, density, composition)
            - pressure(below, temperature, density, composition)
        )
        / (2.0 * step)
        for temperature, density, composition, _ in states
    )
    owner_derivatives = tuple(
        derivative * row.scale / row.multiplier
        for derivative, row in zip(
            accepted.physical_jacobian, result.specification.rows, strict=True
        )
    )
    assert all(
        scaled_error(actual, expected) <= 5.0e-5
        for actual, expected in zip(owner_derivatives, centered, strict=True)
    )
    assert accepted.rank == 1
    assert accepted.condition_number == pytest.approx(1.0)
    assert accepted.active_bounds == ()
    assert accepted.owner_failures == ()
    assert result.accounting.dropped_rows == 0
    assert result.accounting.rejected_row_trials == 0
    assert result.accounting.failed_row_trials == 0

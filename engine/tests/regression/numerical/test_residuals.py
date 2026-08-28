from __future__ import annotations

import math
from copy import deepcopy
from pathlib import Path
from typing import Any, cast

import pytest

from epcsaft.regression import fit
from tests.regression._cases import (
    PAIR_IDENTITY,
    complete_residuals,
    fixed_state_fit_case,
    parameter_value,
    pressure,
    scaled_error,
    sha256,
    tied_coordinate_fit_case,
    with_parameter,
)


def test_scaled_and_log_rows_apply_declared_whitening_objective_and_gradient(
    tmp_path: Path,
) -> None:
    """Falsify residual transforms, diagonal weights, F, or J.T @ r."""

    mapping, parameters, states = fixed_state_fit_case(
        tmp_path, truth=0.04, initial=-0.02, starts=(("input", -0.02),)
    )
    first_state = states[0]
    original = cast(list[dict[str, Any]], mapping["observations"])[0]
    rows = cast(list[dict[str, float]], cast(dict[str, Any], original["source"])["rows"])
    original["source"]["rows"] = rows[:1]
    original["source"]["sha256"] = sha256(rows[:1])
    original["residuals"] = {"pressure": "log_ratio"}
    original["multiplier"] = 2.0
    scaled = deepcopy(original)
    scaled["source"]["path"] = "inline-scaled-row"
    scaled["source"]["sha256"] = sha256(rows[:1])
    scaled["residuals"] = {"pressure": "scaled_difference"}
    scaled["multiplier"] = 0.25
    mapping["observations"] = [original, scaled]
    coordinate = cast(dict[str, Any], cast(dict[str, Any], mapping["active"])["coordinates"])[
        PAIR_IDENTITY
    ]
    coordinate["bounds"] = [-0.10, 0.0]

    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None
    assert accepted.final_physical == pytest.approx((0.0,), abs=1.0e-12)
    bounded_parameters = with_parameter(parameters, PAIR_IDENTITY, 0.0, "bounded-point")
    predicted = pressure(bounded_parameters, first_state[0], first_state[1], first_state[2])
    observed = first_state[3]
    expected_residuals = (
        2.0 * math.log(predicted / observed),
        0.25 * (predicted - observed) / observed,
    )
    assert accepted.predictions == pytest.approx((predicted, predicted), rel=2.0e-12)
    assert accepted.residuals == pytest.approx(expected_residuals, rel=1.0e-12, abs=1.0e-14)

    step = 1.0e-5
    below = with_parameter(parameters, PAIR_IDENTITY, -step, "residual-below")
    above = with_parameter(parameters, PAIR_IDENTITY, step, "residual-above")
    prediction_derivative = (
        pressure(above, first_state[0], first_state[1], first_state[2])
        - pressure(below, first_state[0], first_state[1], first_state[2])
    ) / (2.0 * step)
    expected_jacobian = (
        2.0 * prediction_derivative / predicted,
        0.25 * prediction_derivative / observed,
    )
    assert all(
        scaled_error(actual, expected) <= 5.0e-5
        for actual, expected in zip(accepted.physical_jacobian, expected_jacobian, strict=True)
    )
    expected_cost = 0.5 * math.fsum(value * value for value in expected_residuals)
    assert accepted.cost == pytest.approx(expected_cost, rel=1.0e-12, abs=1.0e-14)
    gradient = math.fsum(
        residual * derivative
        for residual, derivative in zip(accepted.residuals, accepted.physical_jacobian, strict=True)
    )

    def objective(value: float) -> float:
        candidate = with_parameter(parameters, PAIR_IDENTITY, value, "objective-perturbation")
        value_prediction = pressure(candidate, first_state[0], first_state[1], first_state[2])
        residuals = (
            2.0 * math.log(value_prediction / observed),
            0.25 * (value_prediction - observed) / observed,
        )
        return 0.5 * math.fsum(residual * residual for residual in residuals)

    centered_gradient = (objective(step) - objective(-step)) / (2.0 * step)
    assert scaled_error(gradient, centered_gradient) <= 5.0e-5


def test_bounds_scale_ties_and_covariance_follow_physical_coordinates(tmp_path: Path) -> None:
    """Falsify affine scaling, tied chain rules, active bounds, rank, or covariance."""

    mapping, parameters, states, identities = tied_coordinate_fit_case(tmp_path)
    pair_a, energy, pair_b = identities
    interior = fit(mapping, base_path=tmp_path)
    accepted = interior.best_usable_start
    assert accepted is not None
    assert accepted.final_physical[0] == pytest.approx(0.05, abs=2.0e-8)
    assert accepted.final_affine == pytest.approx(
        (
            accepted.final_physical[0] / 0.05,
            (accepted.final_physical[1] - 270.0) / 300.0,
        )
    )
    assert accepted.rank == accepted.free_subspace_rank == 2
    assert math.isfinite(accepted.condition_number)
    assert accepted.active_bounds == ()
    assert accepted.covariance_status == "available"

    baseline = {
        pair_a: accepted.final_physical[0],
        pair_b: accepted.final_physical[0],
        energy: accepted.final_physical[1],
    }
    coordinate_identities = ((pair_a, pair_b), (energy,))
    steps = (1.0e-5, 1.0e-3)
    physical_columns: list[tuple[float, ...]] = []
    for tied, step in zip(coordinate_identities, steps, strict=True):
        below = complete_residuals(
            parameters,
            states,
            {**baseline, **{identity: baseline[identity] - step for identity in tied}},
        )
        above = complete_residuals(
            parameters,
            states,
            {**baseline, **{identity: baseline[identity] + step for identity in tied}},
        )
        physical_columns.append(
            tuple((high - low) / (2.0 * step) for low, high in zip(below, above, strict=True))
        )
    centered_physical = tuple(
        physical_columns[column][row] for row in range(len(states)) for column in range(2)
    )
    assert all(
        scaled_error(actual, expected) <= 5.0e-5
        for actual, expected in zip(accepted.physical_jacobian, centered_physical, strict=True)
    )

    affine_columns = (
        tuple(value * 0.05 for value in physical_columns[0]),
        tuple(value * 300.0 for value in physical_columns[1]),
    )
    a = math.fsum(value * value for value in affine_columns[0])
    b = math.fsum(left * right for left, right in zip(*affine_columns, strict=True))
    d = math.fsum(value * value for value in affine_columns[1])
    determinant = a * d - b * b
    variance = (
        math.fsum(value * value for value in accepted.residuals) / accepted.degrees_of_freedom
    )
    expected_affine_covariance = (
        variance * d / determinant,
        -variance * b / determinant,
        -variance * b / determinant,
        variance * a / determinant,
    )
    assert accepted.affine_covariance == pytest.approx(
        expected_affine_covariance, rel=5.0e-5, abs=1.0e-24
    )
    assert interior.physical_covariance == pytest.approx(
        tuple(
            expected_affine_covariance[row * 2 + column]
            * (0.05, 300.0)[row]
            * (0.05, 300.0)[column]
            for row in range(2)
            for column in range(2)
        ),
        rel=5.0e-5,
        abs=1.0e-24,
    )

    coordinate = cast(dict[str, Any], cast(dict[str, Any], mapping["active"])["coordinates"])[
        pair_a
    ]
    coordinate["bounds"] = [-0.08, 0.04]
    bounded = fit(mapping, base_path=tmp_path)
    bounded_start = bounded.best_usable_start
    assert bounded_start is not None
    assert bounded_start.final_physical[0] == pytest.approx(0.04, abs=1.0e-12)
    assert bounded_start.active_bounds == (0,)
    assert bounded_start.rank == 2
    assert math.isfinite(bounded_start.condition_number)
    assert bounded_start.covariance_status == "unavailable_active_bounds"
    fitted = {identity: value for identity, value, _ in bounded.fitted_values()}
    assert fitted[pair_a] == fitted[pair_b] == pytest.approx(0.04)
    assert parameter_value(bounded.export_candidate_parameters(), pair_a) == pytest.approx(0.04)

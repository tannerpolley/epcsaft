from __future__ import annotations

import math
from copy import deepcopy
from pathlib import Path
from typing import Any, cast

import pytest

import epcsaft
from epcsaft.eos import Mixture
from epcsaft.regression import fit
from tests.regression._cases import (
    PAIR_IDENTITY,
    fixed_state_fit_case,
    parameter_value,
    pressure,
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

    block = Mixture(bounded_parameters).observable_block(
        epcsaft.ActiveParameterSet(bounded_parameters, (PAIR_IDENTITY,)),
        identity="eos.fixed-state",
        T=first_state[0] * epcsaft.unit_registry.kelvin,
        rho=first_state[1] * epcsaft.unit_registry.mole / epcsaft.unit_registry.meter**3,
        x=first_state[2],
    )
    pressure_row = block.descriptor.row_identities.index("pressure_pa")
    prediction_derivatives = block.jacobian[pressure_row]
    assert block.values[pressure_row] == pytest.approx(predicted, rel=2.0e-12)
    assert prediction_derivatives is not None
    prediction_derivative = prediction_derivatives[0]
    expected_jacobian = (
        2.0 * prediction_derivative / predicted,
        0.25 * prediction_derivative / observed,
    )
    assert accepted.physical_jacobian == pytest.approx(expected_jacobian, rel=1.0e-12, abs=1.0e-14)
    expected_affine_jacobian = tuple(value * 0.05 for value in expected_jacobian)
    assert accepted.affine_jacobian == pytest.approx(
        expected_affine_jacobian, rel=1.0e-12, abs=1.0e-14
    )
    expected_cost = 0.5 * math.fsum(value * value for value in expected_residuals)
    assert accepted.cost == pytest.approx(expected_cost, rel=1.0e-12, abs=1.0e-14)
    physical_gradient = math.fsum(
        residual * derivative
        for residual, derivative in zip(accepted.residuals, accepted.physical_jacobian, strict=True)
    )
    expected_physical_gradient = math.fsum(
        residual * derivative
        for residual, derivative in zip(expected_residuals, expected_jacobian, strict=True)
    )
    assert physical_gradient == pytest.approx(expected_physical_gradient, rel=1.0e-12, abs=1.0e-14)
    affine_gradient = math.fsum(
        residual * derivative
        for residual, derivative in zip(accepted.residuals, accepted.affine_jacobian, strict=True)
    )
    assert affine_gradient == pytest.approx(
        expected_physical_gradient * 0.05, rel=1.0e-12, abs=1.0e-14
    )


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
    candidate = parameters.with_values(
        {
            identity: epcsaft.unit_registry.Quantity(
                value,
                next(item.unit for item in parameters.parameter_specs if item.identity == identity),
            )
            for identity, value in baseline.items()
        },
        source_id="analytic-coordinate-point",
        locator="project-authored exact chain-rule oracle",
    )
    model = Mixture(candidate)
    expected_physical: list[float] = []
    for temperature, density, composition, observed in states:
        block = model.observable_block(
            epcsaft.ActiveParameterSet(candidate, identities),
            identity="eos.fixed-state",
            T=temperature * epcsaft.unit_registry.kelvin,
            rho=density * epcsaft.unit_registry.mole / epcsaft.unit_registry.meter**3,
            x=composition,
        )
        pressure_row = block.descriptor.row_identities.index("pressure_pa")
        derivatives = block.jacobian[pressure_row]
        assert derivatives is not None
        expected_physical.extend(
            ((derivatives[0] + derivatives[2]) / observed, derivatives[1] / observed)
        )
    assert accepted.physical_jacobian == pytest.approx(expected_physical, rel=1.0e-12, abs=1.0e-14)
    expected_affine = tuple(
        value * (0.05, 300.0)[column]
        for row in range(len(states))
        for column, value in enumerate(expected_physical[row * 2 : row * 2 + 2])
    )
    assert accepted.affine_jacobian == pytest.approx(expected_affine, rel=1.0e-12, abs=1.0e-14)

    affine_columns = tuple(
        tuple(expected_affine[row * 2 + column] for row in range(len(states)))
        for column in range(2)
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

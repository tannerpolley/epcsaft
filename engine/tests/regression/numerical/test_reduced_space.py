from __future__ import annotations

import math
from copy import deepcopy
from pathlib import Path

import pytest

from epcsaft import Parameters
from epcsaft.regression import fit
from tests.eos._cases import MULTINEUTRAL_REACTIVE
from tests.regression._cases import (
    PAIR_IDENTITY,
    fixed_state_fit_case,
    full_space_fit_case,
    pressure,
    scaled_error,
    with_parameter,
    write_parameters,
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


@pytest.mark.parametrize(
    "strategy", ("reduced_space", "simultaneous_sparse_full_space")
)
def test_temperature_correlated_domain_trial_is_a_structured_rejection(
    tmp_path: Path,
    strategy: str,
) -> None:
    """Falsify divergent endpoint constraints or loss of their KKT multipliers."""

    fit_mapping, *_ = full_space_fit_case(tmp_path)
    parameter_mapping = deepcopy(MULTINEUTRAL_REACTIVE)
    parameter_mapping["domains"][0].update(
        {
            "temperature_min": {"magnitude": 293.15, "unit": "kelvin"},
            "temperature_max": {"magnitude": 393.15, "unit": "kelvin"},
        }
    )
    pair = parameter_mapping["pairs"][0]
    slope_identity = "pair/cosolvent/solvent/k_ij/reciprocal_temperature_slope"
    pair["coefficients"].append(
        {
            "identity": slope_identity,
            "family": "k_ij_reciprocal_temperature_slope",
            "value": {"magnitude": 0.0, "unit": "kelvin"},
            "reference_temperature": {"magnitude": 313.15, "unit": "kelvin"},
            "provenance": deepcopy(pair["coefficients"][0]["provenance"]),
        }
    )
    parameters = Parameters.from_mapping(parameter_mapping)
    parameter_path, parameter_hash = write_parameters(
        tmp_path, parameters, "temperature-domain-parameters.json"
    )
    identities = (
        "pair/cosolvent/solvent/k_ij",
        "component/solvent/segment_diameter",
        "component/cation/segment_diameter",
        slope_identity,
    )
    origins = (-0.1957268589414558, 2.7533328530317225, 3.017079292071164, 0.0)
    scales = (1.0, 2.15, 2.15, 1000.0)
    diagnosed_affine = (
        1.1281046895595066,
        -0.5829455130380103,
        -0.7056182753819367,
        -0.1223855692929284,
    )
    diagnosed_physical = tuple(
        origin + scale * affine
        for origin, scale, affine in zip(origins, scales, diagnosed_affine, strict=True)
    )
    fit_mapping["parameters"] = {"path": parameter_path, "sha256": parameter_hash}
    fit_mapping["active"] = {
        "identities": list(identities),
        "coordinates": {
            identity: {
                "bounds": [-1.0, 1.0]
                if index == 0
                else [1.5, 5.0]
                if index < 3
                else [-1000.0, 1000.0],
                "origin": origins[index],
                "scale": scales[index],
            }
            for index, identity in enumerate(identities)
        },
        "ties": [],
    }
    fit_mapping["solver"]["evaluation_strategy"] = strategy
    fit_mapping["solver"]["multistart"]["starts"] = [
        *(
            [
                {
                    "identity": "diagnosed-domain-trial",
                    "values": list(diagnosed_physical),
                }
            ]
            if strategy == "reduced_space"
            else []
        ),
        {"identity": "valid-reference", "values": [0.0, 2.8, 2.8, 0.0]},
    ]
    fit_mapping["solver"]["controls"]["maximum_iterations"] = 2

    result = fit(fit_mapping, base_path=tmp_path)
    expected_hot_jacobian = (
        1.0,
        0.0,
        0.0,
        1000.0 * (1.0 / 393.15 - 1.0 / 313.15),
    )
    if strategy == "reduced_space":
        rejected = next(
            start
            for start in result.starts
            if start.start_identity == "diagnosed-domain-trial"
        )
        assert rejected.initial_affine == pytest.approx(diagnosed_affine, abs=2.0e-15)
        assert not rejected.usable
        assert all(math.isfinite(value) for value in rejected.final_physical)
        assert rejected.owner_failures[0].code == "eos_domain_rejection"
        assert result.accounting.rejected_row_trials >= 1
        hot_endpoint = next(
            constraint
            for constraint in rejected.parameter_domain_constraints
            if constraint["temperature_k"] == pytest.approx(393.15)
        )
        assert hot_endpoint["initial_value"] == pytest.approx(1.0119038973581267)
        assert hot_endpoint["initial_value"] > hot_endpoint["upper"]
        assert hot_endpoint["affine_jacobian"] == pytest.approx(expected_hot_jacobian)
    evaluated = next(start for start in result.starts if start.start_identity == "valid-reference")
    if strategy == "reduced_space":
        assert evaluated.owner_failures == ()
    assert all(
        constraint["initial_value"] <= constraint["upper"]
        for constraint in evaluated.parameter_domain_constraints
    )
    assert evaluated.parameter_domain_constraints[1][
        "affine_jacobian"
    ] == pytest.approx(expected_hot_jacobian)
    if strategy == "reduced_space":
        assert all(math.isfinite(value) for value in evaluated.residuals)
        assert all(math.isfinite(value) for value in evaluated.affine_jacobian)
    kkt = evaluated.kkt_diagnostics
    assert kkt["retained"]
    assert len(kkt["parameter_constraint_multipliers"]) == 2
    if strategy == "reduced_space":
        assert kkt["primal_feasibility_inf_norm"] <= 1.0e-12
    assert all(
        math.isfinite(kkt[name])
        for name in (
            "stationarity_inf_norm",
            "complementarity_inf_norm",
            "dual_sign_violation_inf_norm",
        )
    )

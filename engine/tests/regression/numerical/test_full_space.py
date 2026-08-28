from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

import epcsaft
from epcsaft import equilibrium
from epcsaft.regression import fit
from tests.regression._cases import (
    fixed_state_fit_case,
    fixed_state_observation,
    full_space_fit_case,
    pressure,
    scaled_error,
    solved_pressure_fit_case,
    with_parameter,
)


def test_maximum_iterations_is_nonusable_and_retains_final_diagnostics(
    tmp_path: Path,
) -> None:
    """Falsify max-iteration acceptance or loss of an inspectable final state."""

    reduced_mapping, *_ = fixed_state_fit_case(tmp_path)
    reduced_mapping["solver"]["controls"]["maximum_iterations"] = 1
    reduced = fit(reduced_mapping, base_path=tmp_path)
    assert all(start.solver_status == "maximum_iterations_exceeded" for start in reduced.starts)
    assert all(not start.usable for start in reduced.starts)
    assert reduced.best_usable_start is None

    full_mapping, *_ = full_space_fit_case(tmp_path)
    full_mapping["solver"]["controls"]["maximum_iterations"] = 1
    stalled = fit(full_mapping, base_path=tmp_path).starts[0]
    assert stalled.solver_status == "maximum_iterations_exceeded"
    assert not stalled.usable
    assert all(math.isfinite(value) for value in stalled.final_physical)
    diagnostic = stalled.owner_diagnostics[0]
    assert diagnostic["solver_status"]
    assert diagnostic["local_variables"]
    certificate = diagnostic["certificate"]
    assert set(certificate) == {
        "balance_inf_norm",
        "pressure_relative_inf_norm",
        "reaction_affinity_inf_norm",
        "transfer_chemical_potential_inf_norm",
        "kkt_stationarity_inf_norm",
        "phase_distinctness",
        "minimum_amount_mol",
        "minimum_bound_margin",
        "reduced_hessian_min_eigenvalue",
        "reduced_hessian_scale",
        "reduced_hessian_threshold",
    }
    assert all(value is not None for value in certificate.values())
    assert stalled.owner_failures[0].reason == diagnostic["failure_reason"]


def test_finite_owner_full_space_recovers_truth_with_complete_certificate(
    tmp_path: Path,
) -> None:
    """Falsify finite-owner recovery, derivatives, KKT evidence, or mixed-owner order."""

    mapping, parameters, model, problem, identity, truth_value, truth = full_space_fit_case(
        tmp_path
    )
    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None

    coordinate_scale = 0.2
    recovery_limit = 1.0e-6
    observation_relative_tolerance = 2.0e-6
    assert (
        max(abs(value) for value in accepted.affine_jacobian) * recovery_limit
        < observation_relative_tolerance
    )
    assert abs(accepted.final_physical[0] - truth_value) / coordinate_scale <= recovery_limit
    assert accepted.predictions == pytest.approx(
        truth.values, rel=observation_relative_tolerance, abs=2.0e-8
    )
    assert accepted.owner_failures == ()
    assert accepted.batch_outcomes[0].status == "evaluated"
    assert result.accounting.dropped_rows == 0
    assert result.accounting.rejected_row_trials == 0
    assert result.accounting.failed_row_trials == 0

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
    evidence = dict(receipt["evidence"])
    assert evidence["balance_inf_norm"] <= 1.0e-9
    assert evidence["pressure_relative_inf_norm"] <= 1.0e-8
    assert evidence["reaction_affinity_inf_norm"] <= 1.0e-6
    assert evidence["kkt_stationarity_inf_norm"] <= 1.0e-7
    expected_branch = json.loads(json.dumps(truth.branch_signature))
    assert receipt["branch_signature"] == expected_branch
    assert receipt["certificate_scope"] == truth.certificate_scope
    assert receipt["globality_claim"] == truth.globality_claim
    assert {name: outcome["status"] for name, outcome in receipt["outcomes"].items()} == {
        "solver_convergence": "passed",
        "numerical_convergence": "passed",
        "physical_certification": "passed",
        "scientific_fit_acceptance": "not_applicable",
    }
    assert result.evaluated_rows[0].identity == "manufactured-equilibrium-row"
    assert result.evaluated_rows[0].unit == problem.outputs[0].unit
    assert result.evaluated_rows[0].basis == problem.outputs[0].basis

    assert isinstance(problem.temperature, equilibrium.Fixed)
    phase = truth.phases[0]
    temperature = float(problem.temperature.value.to("kelvin").magnitude)
    direct_pressure = pressure(
        with_parameter(parameters, identity, truth_value, "mixed-owner-truth"),
        temperature,
        phase.molar_density_mol_m3,
        phase.mole_fractions,
    )
    direct = fixed_state_observation(
        (
            {
                "temperature_k": temperature,
                "molar_density_mol_per_m3": phase.molar_density_mol_m3,
                "pressure": direct_pressure,
            },
        ),
        phase.mole_fractions,
        identity="mixed-direct-owner",
    )
    observations = mapping["observations"]
    assert isinstance(observations, list)
    mapping["observations"] = [direct, *observations]

    mixed_result = fit(mapping, base_path=tmp_path)
    mixed = mixed_result.best_usable_start
    assert mixed is not None
    assert mixed.predictions == pytest.approx(
        (direct_pressure, float(truth.values[0])), rel=2.0e-6, abs=2.0e-8
    )
    assert tuple(row.selector for row in mixed_result.specification.rows) == (
        "pressure",
        "cation-fraction",
    )
    assert tuple(outcome.status for outcome in mixed.batch_outcomes) == (
        "evaluated",
        "evaluated",
    )


@pytest.mark.slow
@pytest.mark.parametrize("ideal_vapor", [True, False], ids=("ideal", "provider"))
def test_solved_pressure_full_space_derivatives_match_independent_resolves(
    tmp_path: Path,
    ideal_vapor: bool,
) -> None:
    """Falsify solved-pressure liquid response, vapor outputs, or residual scaling."""

    mapping, parameters, model, problem, identity, truth_value, _ = solved_pressure_fit_case(
        tmp_path, ideal_vapor=ideal_vapor
    )
    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None
    assert result.evaluation_strategy == "simultaneous_sparse_full_space"
    assert accepted.final_physical == pytest.approx((truth_value,), abs=3.0e-6)
    assert accepted.owner_failures == ()
    assert result.accounting.failed_row_trials == 0
    assert result.accounting.rejected_row_trials == 0
    assert tuple((row.selector, row.residual, row.scale) for row in result.specification.rows) == (
        ("pressure", "log_ratio", 0.05),
        ("vapor-solvent-fraction", "scaled_difference", 0.01),
        ("vapor-solvent-partial-pressure", "log_ratio", 0.05),
    )

    center_value = accepted.final_physical[0]
    step = 1.0e-5
    sides = tuple(
        equilibrium.solve(
            model,
            problem,
            active_parameters=epcsaft.ActiveParameterSet(
                parameters, (identity,), values=(center_value + sign * step,)
            ),
        )
        for sign in (-1.0, 1.0)
    )
    assert all(side.status == "evaluated" for side in sides)
    centered = tuple(
        (float(high) - float(low)) / (2.0 * step)
        for low, high in zip(sides[0].values, sides[1].values, strict=True)
    )
    owner_derivatives = tuple(
        residual_derivative
        * row.scale
        * (prediction if row.residual == "log_ratio" else 1.0)
        / row.multiplier
        for residual_derivative, row, prediction in zip(
            accepted.physical_jacobian,
            result.specification.rows,
            accepted.predictions,
            strict=True,
        )
    )
    for exact, reference in zip(owner_derivatives, centered, strict=True):
        assert scaled_error(exact, reference) <= 5.0e-5

    if ideal_vapor:
        center = equilibrium.solve(
            model,
            problem,
            active_parameters=epcsaft.ActiveParameterSet(
                parameters, (identity,), values=(center_value,)
            ),
        )
        assert center.status == "evaluated"
        temperature_k = float(problem.temperature.value.to("kelvin").magnitude)
        liquid = center.phases[0]
        fugacities = tuple(
            math.exp(liquid.chemical_potential_over_rt[index])
            * 8.31446261815324
            * temperature_k
            for index in (0, 1)
        )
        pressure_pa = math.fsum(fugacities)
        reference = (pressure_pa, fugacities[0] / pressure_pa, fugacities[0])
        assert center.values == pytest.approx(reference, rel=2.0e-10, abs=2.0e-8)


@pytest.mark.slow
def test_solved_pressure_full_space_certifies_a_closed_declared_parameter_bound(
    tmp_path: Path,
) -> None:
    """Falsify affine-bound rounding or final certification at a closed bound."""

    identity = "pair/cosolvent/solvent/k_ij"
    mapping, _, _, _, _, truth_value, _ = solved_pressure_fit_case(
        tmp_path,
        ideal_vapor=False,
        parameter_identity=identity,
        truth_value=0.01,
        bounds=(-0.01, 0.01),
        origin=0.0,
        scale=0.01,
    )
    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None
    assert result.evaluation_strategy == "simultaneous_sparse_full_space"
    assert accepted.final_physical == pytest.approx((truth_value,), abs=2.0e-8)
    assert accepted.owner_failures == ()


@pytest.mark.slow
def test_solved_pressure_full_space_rejects_only_physical_domain_trials(
    tmp_path: Path,
) -> None:
    """Falsify Ipopt restoration crashes when a trial leaves the EOS domain."""

    identity = "pair/cosolvent/solvent/k_ij"
    mapping, *_ = solved_pressure_fit_case(
        tmp_path,
        ideal_vapor=False,
        parameter_identity=identity,
        truth_value=0.01,
        bounds=(-0.01, 0.01),
        origin=0.0,
        scale=0.01,
    )
    coordinate = mapping["active"]["coordinates"][identity]
    coordinate.update({"bounds": [-1.0, 1.0], "scale": 1.0})
    observation = mapping["observations"][0]
    observation["targets"] = observation["targets"][:1]
    observation["targets"][0]["observed"] *= 100.0
    mapping["solver"]["controls"]["maximum_iterations"] = 20

    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert result.evaluation_strategy == "simultaneous_sparse_full_space"
    assert accepted is not None
    assert accepted.solver_status == "solve_succeeded"
    assert accepted.owner_failures == ()
    assert accepted.final_physical[0] == pytest.approx(0.0157182582, abs=2.0e-8)

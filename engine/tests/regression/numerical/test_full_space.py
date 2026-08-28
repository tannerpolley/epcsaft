from __future__ import annotations

import json
from pathlib import Path

import pytest

import epcsaft
from epcsaft import equilibrium
from epcsaft.regression import fit
from tests.regression._cases import (
    fixed_state_observation,
    full_space_fit_case,
    pressure,
    scaled_error,
    with_parameter,
)


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

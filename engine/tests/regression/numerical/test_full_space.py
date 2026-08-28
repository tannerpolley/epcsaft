from __future__ import annotations

import json
from pathlib import Path

import pytest

import epcsaft
from epcsaft import equilibrium
from epcsaft.regression import fit
from tests.regression._cases import full_space_fit_case, scaled_error


@pytest.mark.slow
def test_finite_owner_full_space_recovers_truth_with_complete_certificate(
    tmp_path: Path,
) -> None:
    """Falsify finite-owner recovery, total derivatives, rows, or KKT evidence."""

    mapping, parameters, model, problem, identity, truth_value, truth = full_space_fit_case(
        tmp_path
    )
    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None

    # The 1e-6 scaled-coordinate gate is stricter than the 2e-6 prediction gate.
    assert abs(accepted.final_physical[0] - truth_value) / 0.2 <= 1.0e-6
    assert accepted.predictions == pytest.approx(truth.values, rel=2.0e-6, abs=2.0e-8)
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

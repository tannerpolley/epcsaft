from __future__ import annotations

from pathlib import Path

import pytest

from epcsaft.regression import FitResult, fit
from tests.regression._cases import (
    PAIR_IDENTITY,
    fixed_state_fit_case,
    parameter_value,
    pressure,
)


def test_fit_result_candidate_re_evaluates_rows_and_preserves_identity(
    tmp_path: Path,
) -> None:
    """Falsify candidate, row, provenance, diagnostic, or fingerprint drift."""

    mapping, input_parameters, states = fixed_state_fit_case(tmp_path)
    mapping["acceptance"] = {"maximum_cost": 1.0e-12}
    result = fit(mapping, base_path=tmp_path)
    assert isinstance(result, FitResult)
    accepted = result.best_usable_start
    assert accepted is not None
    assert accepted is min(
        (start for start in result.starts if start.usable),
        key=lambda start: (start.cost, start.start_identity),
    )
    assert tuple(start.start_identity for start in result.starts) == (
        "input",
        "confirmation",
    )
    assert accepted.owner_failures == ()
    assert accepted.owner_receipts == ()
    assert accepted.rank == 1
    assert accepted.physical_jacobian
    assert accepted.singular_values
    outcomes = result.outcomes
    assert outcomes is not None
    assert {name: outcome["status"] for name, outcome in outcomes.items()} == {
        "solver_convergence": "converged",
        "numerical_convergence": "converged",
        "physical_certification": "certified",
        "scientific_fit_acceptance": "accepted",
    }
    assert outcomes["scientific_fit_acceptance"]["achieved"] == accepted.cost
    assert outcomes["scientific_fit_acceptance"]["threshold"] == 1.0e-12

    candidate = result.export_candidate_parameters()
    assert parameter_value(input_parameters, PAIR_IDENTITY) == pytest.approx(-0.02)
    assert parameter_value(candidate, PAIR_IDENTITY) == pytest.approx(accepted.final_physical[0])
    assert candidate.fingerprint != input_parameters.fingerprint
    re_evaluated = tuple(
        pressure(candidate, temperature, density, composition)
        for temperature, density, composition, _ in states
    )
    assert re_evaluated == pytest.approx(accepted.predictions, rel=2.0e-6, abs=2.0e-8)

    rows = result.evaluated_rows
    assert tuple(row.identity for row in rows) == result.specification.row_identities
    assert tuple(row.predicted for row in rows) == accepted.predictions
    assert tuple(row.scaled_residual for row in rows) == accepted.residuals
    assert all(row.selector == "pressure" for row in rows)
    assert all(row.unit == "pascal" and row.basis == "absolute" for row in rows)
    assert all(row.source_identity and row.source_hash.startswith("sha256:") for row in rows)
    assert all(
        row.candidate_parameter_identity == result.candidate_parameter_identity for row in rows
    )
    assert result.candidate_parameter_identity.startswith("sha256:")
    assert result.evaluated_row_table_hash.startswith("sha256:")
    assert result.artifact_fingerprint.startswith("sha256:")
    assert result.package_fingerprint.startswith("sha256:")
    assert result.digest.startswith("sha256:")
    assert result.accounting.failed_starts == 0
    assert result.accounting.rejected_row_trials == 0
    assert result.accounting.failed_row_trials == 0

    assert accepted.cost > 0.0
    mapping["acceptance"] = {"maximum_cost": accepted.cost / 2.0}
    rejected = fit(mapping, base_path=tmp_path)
    rejected_outcomes = rejected.outcomes
    assert rejected_outcomes is not None
    assert rejected_outcomes["scientific_fit_acceptance"]["status"] == "rejected"

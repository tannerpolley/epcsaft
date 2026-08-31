from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

import epcsaft
from epcsaft import equilibrium
from epcsaft.regression import FitPreflightError, fit
from tests.regression._cases import (
    PURE_IDENTITY,
    boundary_fit_case,
    fixed_state_fit_case,
    pure_saturation_fit_case,
    scaled_error,
)


def test_invalid_saturation_start_is_a_typed_non_evaluable_preflight(
    tmp_path: Path,
) -> None:
    """Falsify conversion of a physical rejection into a plausible residual."""

    mapping, _, _ = pure_saturation_fit_case(tmp_path, initial=20.0, starts=(("invalid", 20.0),))
    mapping["active"]["coordinates"][PURE_IDENTITY]["bounds"] = [20.0, 205.0]
    with pytest.raises(FitPreflightError) as failure:
        fit(mapping, base_path=tmp_path)
    assert failure.value.code == "physical_rejection"


def test_explicit_full_space_rejects_an_ineligible_fit(tmp_path: Path) -> None:
    mapping, *_ = fixed_state_fit_case(tmp_path)
    mapping["solver"]["evaluation_strategy"] = "simultaneous_sparse_full_space"

    with pytest.raises(ValueError, match="full_space evaluation is ineligible"):
        fit(mapping, base_path=tmp_path)


@pytest.mark.slow
def test_non_evaluable_boundary_start_does_not_contaminate_best_candidate(
    tmp_path: Path,
) -> None:
    """Falsify boundary transfer or contamination after a refused start."""

    identity = "component/solvent/segment_diameter"
    mapping, parameters, model, problem, _, truth = boundary_fit_case(
        tmp_path,
        parameter_identity=identity,
        truth_value=2.9,
        bounds=(2.5, 3.5),
        origin=2.9,
        scale=0.5,
        starts=(("non-evaluable", 2.5), ("truth", 2.9)),
    )
    assert isinstance(problem.temperature, equilibrium.Fixed)
    assert truth.values[0] is not None

    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None
    assert accepted.start_identity == "truth"
    assert accepted.final_physical == pytest.approx((2.9,), abs=3.0e-7)
    assert accepted.predictions == pytest.approx(truth.values, rel=2.0e-6, abs=2.0e-8)
    assert tuple(row.selector for row in result.specification.rows) == ("cation-fraction",)
    assert tuple(outcome.status for outcome in accepted.batch_outcomes) == ("evaluated",)

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
    boundary_row = result.specification.rows[0]
    owner_derivative = accepted.physical_jacobian[0] * boundary_row.scale / boundary_row.multiplier
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
    metadata = boundary_row.to_mapping()["series_metadata"]
    assert metadata["aggregate_identity"] == "boundary-output-group"
    assert metadata["covariance_identity"] == "diagonal-whitening-group"

    refused = next(start for start in result.starts if start.start_identity == "non-evaluable")
    assert not refused.usable
    assert refused.trial_failure_code is not None
    assert refused.batch_outcomes[0].status != "evaluated"
    assert not all(math.isfinite(value) for value in refused.predictions)
    assert result.accounting.dropped_rows == 0

from __future__ import annotations

import math
from pathlib import Path

import pytest

from epcsaft.regression import FitPreflightError, fit
from tests.regression._cases import (
    PURE_IDENTITY,
    boundary_fit_case,
    pure_saturation_fit_case,
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


@pytest.mark.slow
def test_non_evaluable_boundary_start_does_not_contaminate_best_candidate(
    tmp_path: Path,
) -> None:
    """Falsify fabricated rows or best-candidate contamination after owner refusal."""

    identity = "component/solvent/segment_diameter"
    mapping, _, _, _, _, _ = boundary_fit_case(
        tmp_path,
        parameter_identity=identity,
        truth_value=2.9,
        bounds=(2.5, 3.5),
        origin=2.9,
        scale=0.5,
        starts=(("non-evaluable", 2.5), ("truth", 2.9)),
    )
    result = fit(mapping, base_path=tmp_path)
    accepted = result.best_usable_start
    assert accepted is not None
    assert accepted.start_identity == "truth"
    assert accepted.final_physical == pytest.approx((2.9,), abs=3.0e-7)
    assert all(math.isfinite(value) for value in accepted.predictions)
    assert accepted.batch_outcomes[0].status == "evaluated"

    refused = next(start for start in result.starts if start.start_identity == "non-evaluable")
    assert not refused.usable
    assert refused.trial_failure_code is not None
    assert refused.batch_outcomes[0].status != "evaluated"
    assert not all(math.isfinite(value) for value in refused.predictions)
    assert result.accounting.dropped_rows == 0

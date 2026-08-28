from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from epcsaft.regression import FitResult, fit
from tests.regression._cases import (
    PAIR_IDENTITY,
    canonical_bytes,
    fixed_state_fit_case,
    parameter_value,
    pressure,
)


def test_fit_result_candidate_re_evaluates_rows_and_preserves_identity(
    tmp_path: Path,
) -> None:
    """Falsify candidate, row, provenance, diagnostic, or fingerprint drift."""

    mapping, input_parameters, states = fixed_state_fit_case(tmp_path)
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
    expected_table_hash = (
        "sha256:" + hashlib.sha256(canonical_bytes([row.to_mapping() for row in rows])).hexdigest()
    )
    assert result.evaluated_row_table_hash == expected_table_hash
    assert result.artifact_fingerprint.startswith("sha256:")
    assert result.package_fingerprint.startswith("sha256:")
    assert result.digest.startswith("sha256:")
    serialized = json.loads(result.to_json_bytes())
    assert serialized["digest"] == result.digest
    assert serialized["candidate_parameter_identity"] == result.candidate_parameter_identity
    assert serialized["accounting"] == result.accounting.to_mapping()

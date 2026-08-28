from __future__ import annotations

from pathlib import Path

import pytest

from epcsaft import _core
from epcsaft.regression import fit
from tests.regression._cases import (
    fixed_state_fit_case,
    full_space_fit_case,
    native_fit_arguments,
)


@pytest.mark.slow
def test_native_owner_selects_both_strategies_without_public_result_change(
    tmp_path: Path,
) -> None:
    """Keep the two native identities separate without a public method label."""

    reduced_mapping, _, _ = fixed_state_fit_case(tmp_path, starts=(("truth", 0.04),))
    full_mapping, *_ = full_space_fit_case(tmp_path)
    full_mapping["solver"]["multistart"]["starts"] = [{"identity": "truth", "values": [2.9]}]
    full_mapping["solver"]["controls"]["maximum_iterations"] = 2

    reduced = _core._fit_resolved(*native_fit_arguments(reduced_mapping, tmp_path))
    full = _core._fit_resolved(*native_fit_arguments(full_mapping, tmp_path))
    assert reduced["evaluation_strategy"] == "reduced_space"
    assert full["evaluation_strategy"] == "simultaneous_sparse_full_space"

    public = fit(reduced_mapping, base_path=tmp_path)
    assert not hasattr(public, "evaluation_strategy")
    assert "evaluation_strategy" not in public.to_mapping()

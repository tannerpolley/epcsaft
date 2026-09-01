from __future__ import annotations

import pytest

from epcsaft import _core


def test_farkas_certificate_accepts_a_contradiction_and_rejects_an_invalid_ray() -> None:
    accepted = _core._phase_stability_audit_farkas_certificate(
        [[1.0], [1.0]],
        [float("-inf"), 1.0],
        [0.0, float("inf")],
        [0.0],
        [1.0],
        [1.0, -1.0],
    )
    permuted = _core._phase_stability_audit_farkas_certificate(
        [[1.0], [1.0]],
        [1.0, float("-inf")],
        [float("inf"), 0.0],
        [0.0],
        [1.0],
        [-1.0, 1.0],
    )
    rejected = _core._phase_stability_audit_farkas_certificate(
        [[2.0e-10]],
        [float("-inf")],
        [-1.0],
        [float("-inf")],
        [-1.0e11],
        [1.0],
    )

    assert accepted["accepted"] is permuted["accepted"] is True
    assert accepted["contradiction_margin"] > accepted["contradiction_threshold"]
    assert accepted["contradiction_margin"] == pytest.approx(permuted["contradiction_margin"])
    assert rejected["accepted"] is False


def test_step5_kkt_certificate_accepts_valid_bounds_and_rejects_bad_closure() -> None:
    accepted = _core._phase_stability_audit_step5_kkt(
        [0.0, 0.5],
        [0.0, 0.0],
        [1.0, 1.0],
        [1.0, 0.0],
        [],
        [],
        [1.0, 0.0],
        [0.0, 0.0],
        [],
        0.0,
        1.0,
        0.0,
        True,
        True,
        True,
    )
    rejected = _core._phase_stability_audit_step5_kkt(
        [0.0, 0.5],
        [0.0, 0.0],
        [1.0, 1.0],
        [1.0, 0.0],
        [],
        [],
        [1.0, 0.0],
        [0.0, 0.0],
        [],
        1.0e-3,
        1.0,
        0.0,
        True,
        True,
        True,
    )

    assert accepted["accepted"] is True
    assert accepted["stationarity_residual_inf"] <= 1.0e-12
    assert rejected["accepted"] is False

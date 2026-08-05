from __future__ import annotations

from pathlib import Path

import epcsaft.regression._native as native
import pytest

from epcsaft import Mixture, Parameters
from epcsaft.eos import _native_transport
from epcsaft.regression.records import FIGIEL_BORN_DIAMETER_TRACER_V1
from epcsaft.regression.workflow import _born_native_payload


def _models() -> tuple[Mixture, ...]:
    bundle = (
        Path(__file__).parents[1]
        / "fixtures"
        / "data"
        / "figiel-2025-reference-electrolytes"
        / "1"
        / "parameters"
    )
    return tuple(
        Mixture(Parameters.from_bundle(bundle, components=target.component_order))
        for target in FIGIEL_BORN_DIAMETER_TRACER_V1.targets
    )


@pytest.mark.campaign
def test_engine_model_born_derivatives_match_step_halved_value_differences() -> None:
    specification = FIGIEL_BORN_DIAMETER_TRACER_V1
    models = _models()
    capsules = tuple(_native_transport(model) for model in models)
    model_fingerprints = tuple(model.parameter_fingerprint for model in models)
    payload = _born_native_payload(specification, model_fingerprints)
    trial = (2.9, 3.2, 4.3, 3.9, 4.7)
    residuals, jacobian, rows, fingerprints, compiled_identity = native.evaluate_born(
        capsules, payload, trial
    )

    assert tuple(fingerprints) == model_fingerprints
    assert tuple(compiled_identity) == payload[0]
    assert len(residuals) == 5
    assert len(jacobian) == 25
    assert all(
        jacobian[row * 5 + column] == 0.0
        for row in range(5)
        for column in range(5)
        if row != column
    )
    for index, (diameter, row) in enumerate(zip(trial, rows, strict=True)):
        step = 1.0e-4
        values = []
        for h in (step, step / 2.0):
            plus = list(trial)
            minus = list(trial)
            plus[index] = diameter + h
            minus[index] = diameter - h
            value_plus = native.evaluate_born(capsules, payload, tuple(plus))[2][index][0]
            value_minus = native.evaluate_born(capsules, payload, tuple(minus))[2][index][0]
            values.append((value_plus - value_minus) / (2.0 * h))
        derivative = float(row[1])
        tolerance = max(1.0e-3, 20.0 * abs(values[0] - values[1]), 2.0e-8 * abs(derivative))
        assert abs(derivative - values[1]) <= tolerance
        assert jacobian[index * 5 + index] == pytest.approx(
            derivative / abs(specification.targets_j_per_mol[index]), rel=0.0, abs=0.0
        )

from __future__ import annotations

import math
from pathlib import Path

import epcsaft.regression._native as native
import pytest

from epcsaft import Mixture, Parameters, native_sdk
from epcsaft import regression
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
    capsules = tuple(native_sdk(model) for model in models)
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


@pytest.mark.campaign
def test_five_ion_born_fit_accepts_observable_recovery_and_reports_parameter_deltas() -> None:
    result = regression.fit_figiel_born_diameters(models=_models())

    assert result.solver_converged
    assert result.numerically_converged
    assert result.workflow_valid
    assert result.scientifically_valid
    assert result.predictive_status == "NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF"
    assert tuple(start.name for start in result.starts) == ("primary", "lower", "upper")
    assert all(start.termination == "CONVERGENCE" for start in result.starts)
    assert all(start.solution_usable for start in result.starts)
    assert all(start.rank == 5 for start in result.starts)
    assert all(math.isfinite(start.condition_number) for start in result.starts)
    assert all(not parameter.active_bound for parameter in result.parameters)
    assert tuple(
        parameter.final_diameter_angstrom for parameter in result.parameters
    ) == pytest.approx(
        (
            2.7888130173797934,
            3.4524616464076425,
            4.147266741279482,
            4.101505615791675,
            4.476998527506598,
        ),
        rel=2.0e-12,
        abs=2.0e-12,
    )
    assert all(
        abs(row.raw_error_j_per_mol)
        <= FIGIEL_BORN_DIAMETER_TRACER_V1.observable_round_trip_j_per_mol
        for row in result.observations
    )
    assert result.confirmation_parameter_scaled_max_deltas[0] <= 1.0e-5
    assert result.confirmation_parameter_scaled_max_deltas[1] <= 1.0e-5
    # The published Table 3 values are rounded comparison anchors, not residual
    # targets or a second acceptance gate.
    assert max(abs(parameter.published_delta_angstrom) for parameter in result.parameters) > 0.0005
    assert result.failure_reasons == ()

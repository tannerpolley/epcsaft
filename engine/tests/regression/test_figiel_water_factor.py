from __future__ import annotations

import math

import pytest
from parameter_cases import aqueous_parameters

from epcsaft import Mixture, Parameters, native_sdk
from epcsaft import regression
from epcsaft.regression import _native
from epcsaft.regression.workflow import _water_factor_native_payload


def _model() -> Mixture:
    return Mixture(
        Parameters.from_dictionary(aqueous_parameters(("water", "sodium-cation", "bromide-anion")))
    )


@pytest.mark.campaign
def test_water_factor_fit_is_rank_one_and_start_confirmed() -> None:
    result = regression.fit_figiel_water_solvation_factor(model=_model())

    assert result.specification_id == (regression.FIGIEL_WATER_SOLVATION_FACTOR_V1.specification_id)
    assert tuple(start.name for start in result.starts) == ("primary", "upper")
    assert len(result.input_row_ids) == 21
    assert result.evaluated_row_ids == result.input_row_ids
    assert result.failed_row_ids == ()
    assert result.solver_converged
    assert result.numerically_converged
    assert result.physically_valid
    assert result.workflow_valid
    assert result.predictive_status == ("NOT_ADJUDICATED_NO_APPROVED_HELD_OUT_CUTOFF")
    assert math.isfinite(result.fitted_water_solvation_factor)
    assert result.fitted_water_solvation_factor == pytest.approx(
        1.5590515389548207, rel=1.0e-12, abs=1.0e-12
    )
    assert 1.0 < result.fitted_water_solvation_factor < 2.0
    assert result.start_parameter_max_abs_delta <= 1.0e-5
    assert result.failure_reasons == ()

    for start in result.starts:
        assert start.termination == "CONVERGENCE"
        assert start.solution_usable
        assert start.final_cost <= start.initial_cost
        assert start.rank == 1
        assert start.condition_number == 1.0
        assert start.complete_jacobian_column
        assert not start.active_bound
        assert len(start.rows) == 21
        assert all(
            math.isfinite(row.modeled_gamma_pm_m)
            and row.modeled_gamma_pm_m > 0.0
            and math.isfinite(row.exact_scaled_residual_derivative)
            and row.exact_scaled_residual_derivative
            == -(row.modeled_gamma_pm_m / row.observed_gamma_pm_m) * row.provider_log_derivative
            for row in start.rows
        )


@pytest.mark.campaign
def test_water_factor_exact_jacobian_matches_centered_callback_values() -> None:
    specification = regression.FIGIEL_WATER_SOLVATION_FACTOR_V1
    model = _model()
    payload = _water_factor_native_payload(specification, model.parameter_fingerprint)
    capsule = native_sdk(model)
    trial = 1.4
    center = _native.evaluate_figiel_water_factor(capsule, payload, trial)
    differences = []
    for step in (1.0e-4, 5.0e-5):
        plus = _native.evaluate_figiel_water_factor(capsule, payload, trial + step)
        minus = _native.evaluate_figiel_water_factor(capsule, payload, trial - step)
        differences.append(
            tuple(
                (plus_row[5] - minus_row[5]) / (2.0 * step)
                for plus_row, minus_row in zip(plus, minus, strict=True)
            )
        )

    for index, row in enumerate(center):
        exact = row[7]
        tolerance = max(
            1.0e-8,
            20.0 * abs(differences[0][index] - differences[1][index]),
            2.0e-8 * abs(exact),
        )
        assert exact == pytest.approx(differences[1][index], rel=0.0, abs=tolerance)

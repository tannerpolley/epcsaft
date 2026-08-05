from __future__ import annotations

import pytest
from parameter_cases import aqueous_parameters

from epcsaft import Mixture, Parameters
from epcsaft.eos import _native_transport
from epcsaft import regression
from epcsaft.regression import _native
from epcsaft.regression.workflow import _aqueous_kij_native_payload


def _models() -> tuple[Mixture, ...]:
    return tuple(
        Mixture(Parameters.from_dictionary(aqueous_parameters(("water", cation, anion))))
        for _, cation, anion, _ in regression.FIGIEL_AQUEOUS_KIJ_V1.salt_contracts
    )


@pytest.mark.campaign
def test_aqueous_kij_residual_jacobian_maps_all_exact_provider_columns() -> None:
    specification = regression.FIGIEL_AQUEOUS_KIJ_V1
    models = _models()
    capsules = tuple(_native_transport(model) for model in models)
    payload = _aqueous_kij_native_payload(
        specification, tuple(model.parameter_fingerprint for model in models)
    )

    _, jacobian_native, rows_native = _native.evaluate_figiel_kij(
        capsules, payload, specification.published_parameters
    )
    jacobian = tuple(float(value) for value in jacobian_native)
    column_map = {salt: columns for salt, _, _, columns in specification.salt_contracts}
    for row_index, row_native in enumerate(rows_native):
        modeled = float(row_native[5])
        observed = float(row_native[3])
        local = tuple(float(value) for value in row_native[7])
        expected = [0.0] * 11
        for column, derivative in zip(column_map[str(row_native[1])], local, strict=True):
            expected[column] = -(modeled / observed) * derivative
        actual = jacobian[row_index * 11 : (row_index + 1) * 11]
        assert actual == pytest.approx(expected, rel=1.0e-14, abs=1.0e-14)


@pytest.mark.campaign
def test_aqueous_kij_exact_jacobian_matches_callback_value_direction() -> None:
    specification = regression.FIGIEL_AQUEOUS_KIJ_V1
    models = _models()
    capsules = tuple(_native_transport(model) for model in models)
    payload = _aqueous_kij_native_payload(
        specification, tuple(model.parameter_fingerprint for model in models)
    )
    parameters = specification.published_parameters
    direction = tuple((-1.0 if index % 2 else 1.0) * (index + 1.0) / 11.0 for index in range(11))
    _, jacobian_native, _ = _native.evaluate_figiel_kij(capsules, payload, parameters)
    jacobian = tuple(float(value) for value in jacobian_native)
    differences: list[tuple[float, ...]] = []
    for step in (1.0e-4, 5.0e-5):
        plus = tuple(
            value + step * component for value, component in zip(parameters, direction, strict=True)
        )
        minus = tuple(
            value - step * component for value, component in zip(parameters, direction, strict=True)
        )
        residuals_plus, _, _ = _native.evaluate_figiel_kij(capsules, payload, plus)
        residuals_minus, _, _ = _native.evaluate_figiel_kij(capsules, payload, minus)
        differences.append(
            tuple(
                (float(plus_value) - float(minus_value)) / (2.0 * step)
                for plus_value, minus_value in zip(residuals_plus, residuals_minus, strict=True)
            )
        )
    for row in range(164):
        exact = sum(jacobian[row * 11 + column] * direction[column] for column in range(11))
        coarse = differences[0][row]
        fine = differences[1][row]
        tolerance = max(
            1.0e-8,
            20.0 * abs(coarse - fine),
            2.0e-8 * abs(exact),
        )
        assert abs(exact - fine) <= tolerance

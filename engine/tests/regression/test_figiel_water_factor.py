from __future__ import annotations


import pytest
from parameter_cases import aqueous_parameters

from epcsaft import Mixture, Parameters
from epcsaft.eos import _native_transport
from epcsaft import regression
from epcsaft.regression import _native
from epcsaft.regression.workflow import _water_factor_native_payload


def _model() -> Mixture:
    return Mixture(
        Parameters.from_dictionary(aqueous_parameters(("water", "sodium-cation", "bromide-anion")))
    )


@pytest.mark.campaign
def test_water_factor_exact_jacobian_matches_centered_callback_values() -> None:
    specification = regression.FIGIEL_WATER_SOLVATION_FACTOR_V1
    model = _model()
    payload = _water_factor_native_payload(specification, model.parameter_fingerprint)
    capsule = _native_transport(model)
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

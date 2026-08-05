from __future__ import annotations

import csv
import ctypes
import math
from pathlib import Path

import pytest
from parameter_dictionaries import PURE_SATURATION_PARAMETERS

import epcsaft
from epcsaft.eos import _native_transport
from epcsaft import equilibrium
from epcsaft.equilibrium import _equilibrium

ANCHORS = Path(__file__).parents[2] / "data" / "reference" / "pure_saturation_anchors.csv"


class _NativeSdkTable(ctypes.Structure):
    _fields_ = (
        ("abi_version", ctypes.c_uint32),
        ("table_size", ctypes.c_size_t),
        ("result_size", ctypes.c_size_t),
        ("model_context", ctypes.c_void_p),
        ("evaluate_pure_phase", ctypes.c_void_p),
        ("parameterized_result_size", ctypes.c_size_t),
        ("evaluate_pure_phase_parameters", ctypes.c_void_p),
        ("component_count", ctypes.c_size_t),
        ("mixture_result_size", ctypes.c_size_t),
        ("evaluate_mixture_phase", ctypes.c_void_p),
    )


_CAPSULE_NAME_BUFFERS: list[ctypes.Array[ctypes.c_char]] = []
_EVALUATORS: list[object] = []
_SDK_TABLES: list[_NativeSdkTable] = []


def _model(component: str = "methane") -> epcsaft.Mixture:
    parameters = epcsaft.Parameters.from_dictionary(PURE_SATURATION_PARAMETERS[component])
    return epcsaft.Mixture(parameters)


@pytest.mark.parametrize(
    ("component", "temperature_k", "variables"),
    (
        ("methane", 150.0, (math.log(1_000.0), math.log(22_000.0), math.log(1.0e6))),
        ("methane", 150.0, (math.log(1.0e-3), math.log(25_000.0), math.log(1.0))),
        ("ethane", 233.15, (math.log(730.0), math.log(20_000.0), math.log(1.0e6))),
        ("ethane", 233.15, (math.log(1.0e-5), math.log(28_000.0), math.log(1.0))),
        ("propane", 300.0, (math.log(500.0), math.log(14_000.0), math.log(1.0e6))),
        ("propane", 300.0, (math.log(1.0e-5), math.log(22_000.0), math.log(1.0))),
    ),
)
def test_saturation_nlp_exact_derivatives_match_independent_directional_differences(
    component: str,
    temperature_k: float,
    variables: tuple[float, float, float],
) -> None:
    model = _model(component)
    capsule = _native_transport(model)
    multipliers = (0.7, -0.4, 1.2)
    direction = (0.3, -0.2, 0.1)
    step = 2.0e-5

    def evaluate(scale: float) -> dict[str, object]:
        point = tuple(
            value + scale * delta for value, delta in zip(variables, direction, strict=True)
        )
        return _equilibrium.evaluate_nlp(
            capsule,
            temperature_k,
            model.parameter_fingerprint,
            point,
            multipliers,
        )

    lower = evaluate(-step)
    center = evaluate(0.0)
    upper = evaluate(step)
    objective_gradient = center["objective_gradient"]
    assert center["objective"] == 0.0
    assert objective_gradient == [0.0, 0.0, 0.0]
    objective_difference = (upper["objective"] - lower["objective"]) / (2.0 * step)
    assert objective_difference == pytest.approx(
        sum(objective_gradient[index] * direction[index] for index in range(3)),
        abs=1.0e-13,
    )

    jacobian = center["jacobian"]
    for row in range(3):
        residual_difference = (upper["constraints"][row] - lower["constraints"][row]) / (2.0 * step)
        jacobian_direction = sum(
            jacobian[row * 3 + column] * direction[column] for column in range(3)
        )
        assert residual_difference == pytest.approx(jacobian_direction, rel=2.0e-7, abs=2.0e-8)

    def lagrangian_gradient(payload: dict[str, object]) -> tuple[float, float, float]:
        payload_jacobian = payload["jacobian"]
        return tuple(
            sum(multipliers[row] * payload_jacobian[row * 3 + column] for row in range(3))
            for column in range(3)
        )

    lower_gradient = lagrangian_gradient(lower)
    upper_gradient = lagrangian_gradient(upper)
    hessian_lower = center["lagrangian_hessian_lower"]
    hessian = (
        (hessian_lower[0], hessian_lower[1], hessian_lower[3]),
        (hessian_lower[1], hessian_lower[2], hessian_lower[4]),
        (hessian_lower[3], hessian_lower[4], hessian_lower[5]),
    )
    for row in range(3):
        gradient_difference = (upper_gradient[row] - lower_gradient[row]) / (2.0 * step)
        hessian_direction = sum(hessian[row][column] * direction[column] for column in range(3))
        assert gradient_difference == pytest.approx(hessian_direction, rel=2.0e-6, abs=2.0e-7)


def test_saturation_matches_retained_lab_and_nist_anchors() -> None:
    with ANCHORS.open(encoding="utf-8", newline="") as stream:
        anchors = list(csv.DictReader(stream))
    assert sorted(anchor["component"] for anchor in anchors) == [
        "ethane",
        "methane",
        "propane",
    ]

    for anchor in anchors:
        component = anchor["component"]
        model = _model(component)
        result = equilibrium.saturation(
            model,
            float(anchor["T_K"]) * epcsaft.unit_registry.kelvin,
        )
        assert result.parameter_fingerprint == anchor["parameter_fingerprint"]
        assert result.saturation_pressure_pa == pytest.approx(
            float(anchor["lab_pressure_pa"]), rel=5.0e-6
        )
        assert result.vapor.molar_density_mol_m3 == pytest.approx(
            float(anchor["lab_vapor_density_mol_m3"]), rel=5.0e-6
        )
        assert result.liquid.molar_density_mol_m3 == pytest.approx(
            float(anchor["lab_liquid_density_mol_m3"]), rel=5.0e-6
        )
        assert result.saturation_pressure_pa == pytest.approx(
            float(anchor["nist_pressure_pa"]), rel=5.0e-3
        )
        liquid_density_kg_m3 = result.liquid.molar_density_mol_m3 * float(
            anchor["molar_mass_kg_mol"]
        )
        assert liquid_density_kg_m3 == pytest.approx(
            float(anchor["nist_liquid_density_kg_m3"]), rel=1.0e-2
        )

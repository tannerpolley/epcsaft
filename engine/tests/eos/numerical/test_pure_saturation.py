from __future__ import annotations

import csv
import math
from copy import deepcopy
from pathlib import Path

import pytest

import epcsaft
from epcsaft import _core
from tests.eos._raw_cases import ASSOCIATING_PURE, set_component_values
from tests.equilibrium.parameter_dictionaries import PURE_SATURATION_PARAMETERS

ANCHORS = Path(__file__).parents[3] / "data" / "reference" / "pure_saturation_anchors.csv"


def _model(component: str = "methane") -> epcsaft.Pure:
    parameters = epcsaft.Parameters.from_mapping(PURE_SATURATION_PARAMETERS[component])
    return epcsaft.Pure(parameters)


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
    multipliers = (0.7, -0.4, 1.2)
    direction = (0.3, -0.2, 0.1)
    step = 2.0e-5

    def evaluate(scale: float) -> dict[str, object]:
        point = tuple(
            value + scale * delta for value, delta in zip(variables, direction, strict=True)
        )
        return _core._evaluate_saturation_nlp(
            model._native,
            temperature_k,
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


def test_saturation_matches_independent_and_nist_anchors() -> None:
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
        result = model.saturation(
            T=float(anchor["T_K"]) * epcsaft.unit_registry.kelvin,
        )
        assert result.saturation_pressure_pa == pytest.approx(
            float(anchor["independent_pressure_pa"]), rel=5.0e-6
        )
        assert result.vapor.molar_density_mol_m3 == pytest.approx(
            float(anchor["independent_vapor_density_mol_m3"]), rel=5.0e-6
        )
        assert result.liquid.molar_density_mol_m3 == pytest.approx(
            float(anchor["independent_liquid_density_mol_m3"]), rel=5.0e-6
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


def test_low_pressure_strongly_associating_saturation_converges() -> None:
    parameters = deepcopy(ASSOCIATING_PURE)
    set_component_values(parameters, "molar_mass", (0.119163,))
    set_component_values(parameters, "segment_count", (3.675,))
    set_component_values(parameters, "segment_diameter", (3.563,))
    set_component_values(parameters, "dispersion_energy_over_k", (228.711,))
    for site in parameters["topology"]["sites"]:
        site["multiplicity"] = 2
    edge = parameters["topology"]["edges"][0]
    edge["energy_over_k"]["value"]["magnitude"] = 2046.624
    edge["volume"]["value"]["magnitude"] = 0.123858

    result = epcsaft.Pure(epcsaft.Parameters.from_mapping(parameters)).saturation(
        T=313.15 * epcsaft.unit_registry.kelvin,
    )

    # An independent scalar pressure solve on the two public density branches
    # brackets the fugacity crossing at 4.228881 Pa.  Below about 1e-6 Pa the
    # dense-branch pressure loses relative precision through cancellation.
    assert result.saturation_pressure_pa == pytest.approx(4.228881, abs=2.0e-6)
    assert result.liquid.molar_density_mol_m3 == pytest.approx(8564.24294, abs=2.0e-5)
    assert abs(result.liquid.pressure_pa - result.saturation_pressure_pa) <= 1.0e-6
    assert abs(result.vapor.pressure_pa - result.saturation_pressure_pa) <= 1.0e-6
    assert result.diagnostics.chemical_potential_absolute_residual <= 1.0e-10
    assert result.diagnostics.solver_constraint_violation <= 1.0e-6


def test_supercritical_saturation_block_rejects_without_fabricated_rows() -> None:
    parameters = epcsaft.Parameters.from_mapping(PURE_SATURATION_PARAMETERS["methane"])
    block = epcsaft.Pure(parameters).observable_block(
        epcsaft.ActiveParameterSet(parameters, ("component/methane/segment_count",)),
        identity="eos.pure-saturation",
        T=200.0 * epcsaft.unit_registry.kelvin,
    )

    assert block.status.status == "rejected"
    assert all(value is None for value in block.values)
    assert all(row is None for row in block.jacobian)
    assert block.saturation_diagnostics is not None
    diagnostics = block.saturation_diagnostics
    assert diagnostics.physical_accepted is False
    assert diagnostics.rejection_status != "none"
    assert diagnostics.exact_derivatives is False
    assert math.isfinite(diagnostics.solver_constraint_violation)
    assert tuple(branch.identity for branch in block.saturation_branches) == ("vapor", "liquid")
    assert all(branch.pressure_relative_residual is None for branch in block.saturation_branches)
    assert all(
        branch.residual_density_derivative_m3_per_mol is None
        for branch in block.saturation_branches
    )

from __future__ import annotations

import math
from copy import deepcopy

import pytest

from epcsaft import ActiveParameterSet, Parameters, Pure, SaturationStart
from epcsaft import unit_registry as u
from tests.eos._cases import ASSOCIATING_PURE, neutral_pure, set_component_values

PURE_SATURATION_PARAMETERS = {
    "methane": neutral_pure("methane"),
    "propane": neutral_pure("propane"),
}
set_component_values(PURE_SATURATION_PARAMETERS["propane"], "molar_mass", (0.044096,))


def _set_parameter(mapping: dict[str, object], identity: str, value: float) -> None:
    pending: list[object] = [mapping]
    while pending:
        item = pending.pop()
        if isinstance(item, dict):
            if item.get("identity") == identity:
                quantity = item["value"]
                assert isinstance(quantity, dict)
                quantity["magnitude"] = value
                return
            pending.extend(item.values())
        elif isinstance(item, list):
            pending.extend(item)
    raise AssertionError(f"parameter identity not found: {identity}")


def _parameter_value(mapping: dict[str, object], identity: str) -> float:
    parameters = Parameters.from_mapping(mapping)
    value = next(spec.value for spec in parameters.parameter_specs if spec.identity == identity)
    return float(getattr(value, "magnitude", value))


def _saturation_values(mapping: dict[str, object]) -> tuple[float, float, float]:
    result = Pure(Parameters.from_mapping(mapping)).saturation(T=300.0 * u.kelvin)
    return (
        result.saturation_pressure_pa,
        result.vapor.molar_density_mol_m3,
        result.liquid.molar_density_mol_m3,
    )


def _scaled_error(exact: float, finite_difference: float) -> float:
    return abs(exact - finite_difference) / max(abs(exact), abs(finite_difference), 2.0e-9)


def test_public_saturation_returns_distinct_stable_equal_potential_branches() -> None:
    model = Pure(Parameters.from_mapping(PURE_SATURATION_PARAMETERS["propane"]))
    result = model.saturation(T=300.0 * u.kelvin)

    assert result.vapor.molar_density_mol_m3 < result.liquid.molar_density_mol_m3
    assert result.diagnostics.vapor_residual_density_derivative_m3_per_mol > 0.0
    assert result.diagnostics.liquid_residual_density_derivative_m3_per_mol > 0.0
    assert result.diagnostics.vapor_pressure_relative_residual <= 1.0e-8
    assert result.diagnostics.liquid_pressure_relative_residual <= 1.0e-8
    assert result.diagnostics.chemical_potential_absolute_residual <= 1.0e-10

    vapor = model.state(
        T=300.0 * u.kelvin,
        rho=result.vapor.molar_density_mol_m3 * u.mole / u.meter**3,
    )
    liquid = model.state(
        T=300.0 * u.kelvin,
        rho=result.liquid.molar_density_mol_m3 * u.mole / u.meter**3,
    )
    assert vapor.pressure.to("pascal").magnitude == pytest.approx(
        result.saturation_pressure_pa, rel=1.0e-8
    )
    assert liquid.pressure.to("pascal").magnitude == pytest.approx(
        result.saturation_pressure_pa, rel=1.0e-8
    )
    assert math.log(result.vapor.molar_density_mol_m3) + vapor.residual_chemical_potential_over_rt[
        0
    ] == pytest.approx(result.vapor.chemical_potential_over_rt, rel=2.0e-12, abs=2.0e-12)
    assert math.log(
        result.liquid.molar_density_mol_m3
    ) + liquid.residual_chemical_potential_over_rt[0] == pytest.approx(
        result.liquid.chemical_potential_over_rt, rel=2.0e-12, abs=2.0e-12
    )
    assert result.vapor.chemical_potential_over_rt == pytest.approx(
        result.liquid.chemical_potential_over_rt, abs=1.0e-10
    )

    uncertified = model.saturation(T=100.0 * u.kelvin)
    assert uncertified.diagnostics.candidate_available
    assert not uncertified.certified
    assert uncertified.diagnostics.sensitivity_status == "not_requested"
    assert uncertified.vapor.molar_density_mol_m3 > 0.0
    assert uncertified.liquid.molar_density_mol_m3 > 0.0

    associated = Pure(Parameters.from_mapping(ASSOCIATING_PURE))
    failed_candidate = associated.saturation(T=180.0 * u.kelvin)
    assert failed_candidate.diagnostics.candidate_available
    assert not failed_candidate.diagnostics.optimizer_converged
    assert not failed_candidate.certified
    assert (
        failed_candidate.diagnostics.feasibility_measure
        > failed_candidate.diagnostics.feasibility_tolerance
    )
    neighbor = associated.saturation(T=290.0 * u.kelvin)
    candidate = associated.saturation(T=300.0 * u.kelvin)
    assert candidate.diagnostics.candidate_available
    assert candidate.certified
    assert candidate.saturation_pressure_pa > 0.0
    assert candidate.vapor.molar_density_mol_m3 < candidate.liquid.molar_density_mol_m3
    assert candidate.diagnostics.sensitivity_status == "not_requested"
    reused = associated.saturation(
        T=300.0 * u.kelvin,
        start=SaturationStart(
            neighbor.vapor.molar_density_mol_m3 * u.mole / u.meter**3,
            neighbor.liquid.molar_density_mol_m3 * u.mole / u.meter**3,
            neighbor.saturation_pressure_pa * u.pascal,
        ),
    )
    assert reused.certified
    assert reused.saturation_pressure_pa == pytest.approx(
        candidate.saturation_pressure_pa, rel=1.0e-10
    )


def test_public_saturation_parameter_derivatives_match_resolved_boundaries() -> None:
    mapping = PURE_SATURATION_PARAMETERS["propane"]
    identities = (
        "component/propane/segment_count",
        "component/propane/segment_diameter",
        "component/propane/dispersion_energy_over_k",
    )
    steps = (1.0e-5, 1.0e-5, 1.0e-3)
    parameters = Parameters.from_mapping(mapping)
    block = Pure(parameters).observable_block(
        ActiveParameterSet(parameters, identities),
        identity="eos.pure-saturation",
        T=300.0 * u.kelvin,
    )
    branches = {branch.identity: branch for branch in block.saturation_branches}
    for row, identity in ((1, "vapor"), (2, "liquid")):
        branch = branches[identity]
        density = block.values[row]
        assert density is not None
        assert branch.density_bracket_mol_per_m3[0] < density
        assert density < branch.density_bracket_mol_per_m3[1]
        assert branch.pressure_relative_residual is not None
        assert branch.pressure_relative_residual <= 1.0e-8
        assert branch.residual_density_derivative_m3_per_mol is not None
        assert branch.residual_density_derivative_m3_per_mol > 0.0

    for column, (identity, step) in enumerate(zip(identities, steps, strict=True)):
        value = _parameter_value(mapping, identity)
        lower = deepcopy(mapping)
        upper = deepcopy(mapping)
        _set_parameter(lower, identity, value - step)
        _set_parameter(upper, identity, value + step)
        numerical = tuple(
            (high - low) / (2.0 * step)
            for low, high in zip(_saturation_values(lower), _saturation_values(upper), strict=True)
        )
        for row, expected in enumerate(numerical):
            derivative = block.jacobian[row]
            assert derivative is not None
            assert _scaled_error(derivative[column], expected) <= 5.0e-5


def test_supercritical_saturation_refuses_without_fabricated_rows() -> None:
    parameters = Parameters.from_mapping(PURE_SATURATION_PARAMETERS["methane"])
    block = Pure(parameters).observable_block(
        ActiveParameterSet(parameters, ("component/methane/segment_count",)),
        identity="eos.pure-saturation",
        T=200.0 * u.kelvin,
    )

    assert block.status.status == "rejected"
    assert all(value is None for value in block.values)
    assert all(row is None for row in block.jacobian)
    assert block.saturation_diagnostics is not None
    assert block.saturation_diagnostics.physical_accepted is False
    assert block.saturation_diagnostics.rejection_status != "none"

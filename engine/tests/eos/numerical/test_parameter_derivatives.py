from __future__ import annotations

from collections.abc import Callable
from copy import deepcopy

import pytest

from epcsaft import ActiveParameterSet, Mixture, ParameterError, Parameters, Pure
from epcsaft import unit_registry as u
from tests.eos._cases import (
    ASSOCIATING_PURE,
    MULTISOLVENT_MULTIION,
    NEUTRAL_BINARY,
    neutral_binary,
    neutral_pure,
)

Values = tuple[float, ...]


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


def _perturbed(mapping: dict[str, object], identity: str, value: float) -> dict[str, object]:
    changed = deepcopy(mapping)
    _set_parameter(changed, identity, value)
    return changed


def _central_difference(
    mapping: dict[str, object],
    identity: str,
    step: float,
    evaluate: Callable[[dict[str, object]], Values],
) -> Values:
    value = _parameter_value(mapping, identity)
    lower = evaluate(_perturbed(mapping, identity, value - step))
    upper = evaluate(_perturbed(mapping, identity, value + step))
    return tuple(
        (high - low) / (2.0 * step)
        for low, high in zip(lower, upper, strict=True)
    )


def _fixed_state_values(
    mapping: dict[str, object],
    temperature_k: float,
    density: float,
    composition: tuple[float, ...],
) -> Values:
    parameters = Parameters.from_mapping(mapping)
    model = Pure(parameters) if len(composition) == 1 else Mixture(parameters)
    state = (
        model.state(
            T=temperature_k * u.kelvin,
            rho=density * u.mole / u.meter**3,
        )
        if isinstance(model, Pure)
        else model.state(
            T=temperature_k * u.kelvin,
            rho=density * u.mole / u.meter**3,
            x=composition,
        )
    )
    return (
        float(state.ares(basis="reduced")),
        float(state.pressure.to("pascal").magnitude),
        *state.residual_chemical_potential_over_rt,
    )


def _scaled_error(exact: float, finite_difference: float) -> float:
    return abs(exact - finite_difference) / max(
        abs(exact), abs(finite_difference), 2.0e-9
    )


@pytest.mark.parametrize(
    ("mapping", "temperature", "density", "composition", "identities", "steps"),
    (
        (
            neutral_pure("propane"),
            300.0,
            4000.0,
            (1.0,),
            (
                "component/propane/segment_count",
                "component/propane/segment_diameter",
                "component/propane/dispersion_energy_over_k",
            ),
            (1.0e-5, 1.0e-5, 1.0e-3),
        ),
        (
            NEUTRAL_BINARY,
            250.0,
            4000.0,
            (0.4, 0.6),
            ("pair/ethane/methane/k_ij",),
            (1.0e-5,),
        ),
        (
            ASSOCIATING_PURE,
            300.0,
            3500.0,
            (1.0,),
            (
                "association/associating-fluid/a/associating-fluid/b/energy_over_k",
                "association/associating-fluid/a/associating-fluid/b/volume",
            ),
            (1.0e-3, 1.0e-6),
        ),
        (
            MULTISOLVENT_MULTIION,
            298.15,
            2000.0,
            (0.60, 0.33, 0.02, 0.02, 0.01, 0.02),
            (
                "component/solvent/relative_permittivity",
                "component/divalent-cation/born_diameter",
                "model/relative_permittivity/ionic_region_relative_permittivity",
            ),
            (1.0e-4, 1.0e-5, 1.0e-5),
        ),
    ),
    ids=("component", "pair", "association", "electrolyte-permittivity"),
)
def test_fixed_state_parameter_matrix_matches_centered_differences(
    mapping,
    temperature,
    density,
    composition,
    identities,
    steps,
) -> None:
    parameters = Parameters.from_mapping(mapping)
    model = Pure(parameters) if len(composition) == 1 else Mixture(parameters)
    active = ActiveParameterSet(parameters, identities)
    exact = (
        model.observable_block(
            active,
            identity="eos.fixed-state",
            T=temperature * u.kelvin,
            rho=density * u.mole / u.meter**3,
        )
        if isinstance(model, Pure)
        else model.observable_block(
            active,
            identity="eos.fixed-state",
            T=temperature * u.kelvin,
            rho=density * u.mole / u.meter**3,
            x=composition,
        )
    )

    for column, (identity, step) in enumerate(zip(identities, steps, strict=True)):
        numerical = _central_difference(
            mapping,
            identity,
            step,
            lambda candidate: _fixed_state_values(
                candidate, temperature, density, composition
            ),
        )
        for row, expected in enumerate(numerical):
            derivative = exact.jacobian[row]
            assert derivative is not None
            assert _scaled_error(derivative[column], expected) <= 5.0e-5


def test_reciprocal_temperature_pair_slope_matches_centered_full_range_resolves() -> None:
    """Falsify fit/export divergence, endpoint admission, or the exact slope column."""

    mapping = deepcopy(NEUTRAL_BINARY)
    mapping["domains"] = [
        {
            "domain_id": mapping["domains"][0]["domain_id"],
            "kind": "fit-range",
            "temperature_min": {"magnitude": 293.15, "unit": "kelvin"},
            "temperature_max": {"magnitude": 393.15, "unit": "kelvin"},
        }
    ]
    pair = mapping["pairs"][0]
    slope_identity = "pair/ethane/methane/k_ij/reciprocal_temperature_slope"
    pair["coefficients"].append(
        {
            "identity": slope_identity,
            "family": "k_ij_reciprocal_temperature_slope",
            "value": {"magnitude": 120.0, "unit": "kelvin"},
            "reference_temperature": {"magnitude": 313.15, "unit": "kelvin"},
            "provenance": deepcopy(pair["coefficients"][0]["provenance"]),
        }
    )
    parameters = Parameters.from_mapping(mapping)
    model = Mixture(parameters)
    identities = (slope_identity, "pair/ethane/methane/k_ij")

    for temperature_k in (293.15, 313.15, 393.15):
        effective_k = 120.0 * (1.0 / temperature_k - 1.0 / 313.15)
        static = deepcopy(NEUTRAL_BINARY)
        _set_parameter(static, "pair/ethane/methane/k_ij", effective_k)
        assert _fixed_state_values(
            mapping, temperature_k, 4000.0, (0.4, 0.6)
        ) == pytest.approx(
            _fixed_state_values(static, temperature_k, 4000.0, (0.4, 0.6)),
            rel=2.0e-13,
            abs=2.0e-11,
        )

    def block(temperature_k: float, slope_k: float):
        return model.observable_block(
            ActiveParameterSet(parameters, identities, values=(slope_k, 0.05)),
            identity="eos.fixed-state",
            T=temperature_k * u.kelvin,
            rho=4000.0 * u.mole / u.meter**3,
            x=(0.4, 0.6),
        )

    step = 1.0e-2
    for temperature_k in (293.15, 393.15):
        center = block(temperature_k, 120.0)
        lower = block(temperature_k, 120.0 - step)
        upper = block(temperature_k, 120.0 + step)
        numerical = tuple(
            (high - low) / (2.0 * step)
            for low, high in zip(lower.values, upper.values, strict=True)
        )
        assert abs(center.jacobian[1][0]) > 2.0e-9
        assert all(
            _scaled_error(row[0], expected) <= 5.0e-5
            for row, expected in zip(center.jacobian, numerical, strict=True)
        )

    reference = block(313.15, 120.0)
    assert all(abs(row[0]) <= 2.0e-9 for row in reference.jacobian)

    invalid = deepcopy(mapping)
    invalid["pairs"][0]["coefficients"][1]["value"]["magnitude"] = 10_000.0
    with pytest.raises(ParameterError, match="invalid-domain"):
        Mixture(Parameters.from_mapping(invalid))


def test_fixed_pressure_fugacity_derivatives_include_density_closure() -> None:
    mapping = neutral_pure("propane")
    identities = (
        "component/propane/segment_diameter",
        "component/propane/dispersion_energy_over_k",
    )
    parameters = Parameters.from_mapping(mapping)
    exact = Pure(parameters).observable_block(
        ActiveParameterSet(parameters, identities),
        identity="eos.fixed-pressure-fugacity",
        T=300.0 * u.kelvin,
        P=1.0e5 * u.pascal,
        phase="vapor",
    )

    def value(candidate: dict[str, object]) -> Values:
        state = Pure(Parameters.from_mapping(candidate)).state(
            T=300.0 * u.kelvin,
            P=1.0e5 * u.pascal,
            phase="vapor",
        )
        assert state.fugacity is not None
        return state.fugacity.ln_coefficient

    for column, (identity, step) in enumerate(
        zip(identities, (1.0e-5, 1.0e-3), strict=True)
    ):
        numerical = _central_difference(mapping, identity, step, value)
        assert _scaled_error(exact.jacobian[0][column], numerical[0]) <= 5.0e-5


def test_parameter_jacobian_permutation_resolves_by_component_identity() -> None:
    identity = "pair/ethane/methane/k_ij"
    canonical_parameters = Parameters.from_mapping(NEUTRAL_BINARY)
    canonical = Mixture(canonical_parameters).observable_block(
        ActiveParameterSet(canonical_parameters, (identity,)),
        identity="eos.fixed-state",
        T=250.0 * u.kelvin,
        rho=4000.0 * u.mole / u.meter**3,
        x=(0.4, 0.6),
    )
    reversed_parameters = Parameters.from_mapping(
        neutral_binary(("ethane", "methane"))
    )
    reversed_block = Mixture(reversed_parameters).observable_block(
        ActiveParameterSet(reversed_parameters, (identity,)),
        identity="eos.fixed-state",
        T=250.0 * u.kelvin,
        rho=4000.0 * u.mole / u.meter**3,
        x=(0.6, 0.4),
    )

    # Rows are ares, pressure, then component-labelled chemical potentials.
    for canonical_row, reversed_row in enumerate((0, 1, 3, 2)):
        assert reversed_block.values[reversed_row] == pytest.approx(
            canonical.values[canonical_row], rel=2.0e-12, abs=2.0e-12
        )
        assert reversed_block.jacobian[reversed_row][0] == pytest.approx(
            canonical.jacobian[canonical_row][0], rel=2.0e-10, abs=2.0e-10
        )

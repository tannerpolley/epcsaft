from __future__ import annotations

from collections.abc import Callable
from copy import deepcopy

import pytest

from epcsaft import ActiveParameterSet, Mixture, Parameters, Pure
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

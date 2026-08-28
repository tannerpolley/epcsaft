from __future__ import annotations

import math

import pytest

from epcsaft import Mixture, Parameters, Pure
from epcsaft import unit_registry as u
from tests.eos._cases import NEUTRAL_BINARY, neutral_binary, neutral_pure


def _first_derivative(function, value: float, step: float) -> float:
    return (
        -function(value + 2.0 * step)
        + 8.0 * function(value + step)
        - 8.0 * function(value - step)
        + function(value - 2.0 * step)
    ) / (12.0 * step)


def _second_derivative(function, value: float, step: float) -> float:
    return (
        -function(value + 2.0 * step)
        + 16.0 * function(value + step)
        - 30.0 * function(value)
        + 16.0 * function(value - step)
        - function(value - 2.0 * step)
    ) / (12.0 * step * step)


def _state(mapping, composition, *, temperature=250.0, density=1000.0):
    return Mixture(Parameters.from_mapping(mapping)).state(
        T=temperature * u.kelvin,
        rho=density * u.mole / u.meter**3,
        x=composition,
    )


def test_pressure_and_compressibility_identity() -> None:
    density = 40.0
    temperature = 298.15
    state = _state(
        NEUTRAL_BINARY,
        (0.4, 0.6),
        temperature=temperature,
        density=density,
    )

    expected = state.compressibility_factor * density * 8.31446261815324 * temperature

    assert state.pressure.to("pascal").magnitude == pytest.approx(
        expected, rel=2.0e-14, abs=2.0e-8
    )


def test_chemical_potential_and_fugacity_follow_amount_derivatives() -> None:
    model = Mixture(Parameters.from_mapping(NEUTRAL_BINARY))
    temperature = 250.0
    density = 8000.0
    fractions = (0.4, 0.6)
    volume = 1.0 / density
    state = _state(
        NEUTRAL_BINARY,
        fractions,
        temperature=temperature,
        density=density,
    )
    expected = []
    for component in range(2):

        def total_helmholtz(delta: float, component_index: int = component) -> float:
            amounts = [*fractions]
            amounts[component_index] += delta
            total = math.fsum(amounts)
            perturbed = model.state(
                T=temperature * u.kelvin,
                rho=(total / volume) * u.mole / u.meter**3,
                x=tuple(value / total for value in amounts),
            )
            return total * perturbed.ares(basis="reduced")

        expected.append(_first_derivative(total_helmholtz, 0.0, 2.0e-5))

    assert state.residual_chemical_potential_over_rt == pytest.approx(
        expected, rel=2.0e-7, abs=2.0e-8
    )
    assert state.fugacity is not None
    assert state.fugacity.ln_coefficient == pytest.approx(
        tuple(value - math.log(state.compressibility_factor) for value in expected),
        rel=2.0e-7,
        abs=2.0e-8,
    )


def test_caloric_properties_follow_helmholtz_identities() -> None:
    temperature = 250.0
    density = 1000.0
    fractions = (0.4, 0.6)
    state = _state(
        NEUTRAL_BINARY,
        fractions,
        temperature=temperature,
        density=density,
    )

    def helmholtz_at(value: float) -> float:
        return _state(
            NEUTRAL_BINARY,
            fractions,
            temperature=value,
            density=density,
        ).ares(basis="reduced")

    temperature_step = temperature * 2.0e-4
    a_t = _first_derivative(helmholtz_at, temperature, temperature_step)
    a_tt = _second_derivative(helmholtz_at, temperature, temperature_step)

    def z_at_temperature(value: float) -> float:
        return _state(
            NEUTRAL_BINARY,
            fractions,
            temperature=value,
            density=density,
        ).compressibility_factor

    def z_at_density(value: float) -> float:
        return _state(
            NEUTRAL_BINARY,
            fractions,
            temperature=temperature,
            density=value,
        ).compressibility_factor

    z = state.compressibility_factor
    z_t = _first_derivative(z_at_temperature, temperature, temperature_step)
    z_rho = _first_derivative(z_at_density, density, density * 2.0e-4)
    cv = -2.0 * temperature * a_t - temperature**2 * a_tt
    expected = {
        "ures": -temperature * a_t,
        "hres": -temperature * a_t + z - 1.0,
        "sres": -state.ares(basis="reduced") - temperature * a_t,
        "cvres": cv,
        "cpres": cv + (z + temperature * z_t) ** 2 / (z + density * z_rho) - 1.0,
    }

    for method, value in expected.items():
        assert getattr(state, method)(basis="reduced") == pytest.approx(
            value, rel=3.0e-5, abs=3.0e-5
        )


def test_pure_limit_and_component_permutation_preserve_identity() -> None:
    mixture = Mixture(Parameters.from_mapping(NEUTRAL_BINARY)).state(
        T=250.0 * u.kelvin,
        rho=1000.0 * u.mole / u.meter**3,
        x=(1.0, 0.0),
    )
    pure = Pure(Parameters.from_mapping(neutral_pure("methane"))).state(
        T=250.0 * u.kelvin,
        rho=1000.0 * u.mole / u.meter**3,
    )
    forward = _state(NEUTRAL_BINARY, (0.4, 0.6))
    reverse = _state(neutral_binary(("ethane", "methane")), (0.6, 0.4))

    for field in ("hard_chain", "dispersion", "compressibility_factor"):
        assert getattr(mixture, field) == pytest.approx(
            getattr(pure, field), rel=2.0e-13, abs=2.0e-13
        )
        assert getattr(reverse, field) == pytest.approx(
            getattr(forward, field), rel=2.0e-13, abs=2.0e-13
        )
    by_id = dict(zip(("methane", "ethane"), forward.residual_chemical_potential_over_rt))
    reversed_by_id = dict(
        zip(("ethane", "methane"), reverse.residual_chemical_potential_over_rt)
    )
    assert reversed_by_id == pytest.approx(by_id, rel=2.0e-12, abs=2.0e-12)

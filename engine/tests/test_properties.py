import math

import pytest

from epcsaft import EPCSAFT, ParameterBundle, unit_registry


GAS_CONSTANT = 8.31446261815324


def _model() -> EPCSAFT:
    return EPCSAFT(
        ParameterBundle.from_catalog(
            "gross-2001-methane-ethane", version=1
        ).select(("methane", "ethane"))
    )


def _evaluate(
    model: EPCSAFT,
    temperature: float,
    density: float,
    fractions: tuple[float, ...],
):
    return model.evaluate(
        temperature=temperature * unit_registry.kelvin,
        molar_density=density * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=fractions,
    )


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


def test_residual_chemical_potentials_are_total_mole_derivatives() -> None:
    model = _model()
    temperature = 250.0
    density = 8000.0
    fractions = (0.4, 0.6)
    volume = 1.0 / density
    result = _evaluate(model, temperature, density, fractions)

    expected = []
    for component in range(2):
        def total_residual_helmholtz(delta: float) -> float:
            moles = [*fractions]
            moles[component] += delta
            total_moles = math.fsum(moles)
            perturbed_fractions = tuple(value / total_moles for value in moles)
            perturbed_density = total_moles / volume
            return total_moles * _evaluate(
                model, temperature, perturbed_density, perturbed_fractions
            ).residual_helmholtz

        expected.append(_first_derivative(total_residual_helmholtz, 0.0, 2.0e-5))

    assert result.residual_chemical_potential_over_rt == pytest.approx(
        expected, rel=2.0e-7, abs=2.0e-8
    )
    assert result.log_fugacity_coefficient == pytest.approx(
        tuple(value - math.log(result.compressibility_factor) for value in expected),
        rel=2.0e-7,
        abs=2.0e-8,
    )


def test_residual_caloric_properties_follow_helmholtz_identities() -> None:
    model = _model()
    temperature = 250.0
    density = 8000.0
    fractions = (0.4, 0.6)
    result = _evaluate(model, temperature, density, fractions)

    def helmholtz_at(value: float) -> float:
        return _evaluate(model, value, density, fractions).residual_helmholtz

    temperature_step = temperature * 2.0e-4
    a_t = _first_derivative(helmholtz_at, temperature, temperature_step)
    a_tt = _second_derivative(helmholtz_at, temperature, temperature_step)

    def z_at_temperature(value: float) -> float:
        return _evaluate(model, value, density, fractions).compressibility_factor

    def z_at_density(value: float) -> float:
        return _evaluate(model, temperature, value, fractions).compressibility_factor

    z_t = _first_derivative(z_at_temperature, temperature, temperature_step)
    z_rho = _first_derivative(z_at_density, density, density * 2.0e-4)
    z = result.compressibility_factor
    expected = {
        "residual_entropy": GAS_CONSTANT
        * (-result.residual_helmholtz - temperature * a_t),
        "residual_internal_energy": -GAS_CONSTANT * temperature**2 * a_t,
        "residual_enthalpy": GAS_CONSTANT
        * temperature
        * (-temperature * a_t + z - 1.0),
        "residual_isochoric_heat_capacity": GAS_CONSTANT
        * (-2.0 * temperature * a_t - temperature**2 * a_tt),
        "residual_isobaric_heat_capacity": GAS_CONSTANT
        * (
            -2.0 * temperature * a_t
            - temperature**2 * a_tt
            + (z + temperature * z_t) ** 2 / (z + density * z_rho)
            - 1.0
        ),
    }

    for field, magnitude in expected.items():
        actual = getattr(result, field)
        unit = "joule / mole" if field.endswith(("energy", "enthalpy")) else (
            "joule / mole / kelvin"
        )
        assert actual.to(unit).magnitude == pytest.approx(
            magnitude, rel=3.0e-5, abs=3.0e-5
        )


def test_electrolyte_formula_unit_chemical_potential_is_electroneutral() -> None:
    model = EPCSAFT(
        ParameterBundle.from_catalog(
            "figiel-2025-water-sodium-chloride", version=1
        ).select(("water", "sodium-cation", "chloride-anion"))
    )
    temperature = 298.15
    density = 40.0
    fractions = (0.9, 0.05, 0.05)
    volume = 1.0 / density
    result = _evaluate(model, temperature, density, fractions)

    def total_residual_helmholtz(delta: float) -> float:
        moles = (fractions[0], fractions[1] + delta, fractions[2] + delta)
        total_moles = math.fsum(moles)
        return total_moles * _evaluate(
            model,
            temperature,
            total_moles / volume,
            tuple(value / total_moles for value in moles),
        ).residual_helmholtz

    expected = _first_derivative(total_residual_helmholtz, 0.0, 2.0e-6)
    actual = math.fsum(result.residual_chemical_potential_over_rt[1:])

    assert actual == pytest.approx(expected, rel=3.0e-6, abs=3.0e-6)

import math

import pytest

from epcsaft import EPCSAFT, EosError, ParameterBundle, unit_registry


ELEMENTARY_CHARGE = 1.602176634e-19
VACUUM_PERMITTIVITY = 8.8541878128e-12
BOLTZMANN = 1.380649e-23
AVOGADRO = 6.02214076e23
GAS_CONSTANT = 8.31446261815324


def _electrolyte_model() -> EPCSAFT:
    return EPCSAFT(
        ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1).select(
            ("water", "sodium-cation", "chloride-anion")
        )
    )


def _debye_huckel_oracle(temperature: float, density: float) -> float:
    fractions = (0.9, 0.05, 0.05)
    charges = (0, 1, -1)
    sigmas_angstrom = (
        2.7927
        + 10.11 * math.exp(-0.01775 * temperature)
        - 1.417 * math.exp(-0.01146 * temperature),
        2.8232,
        2.7560,
    )
    ion_fraction = fractions[1] + fractions[2]
    relative_permittivity = 78.09 / (1.0 + 7.01 * ion_fraction)
    ionic_strength_sum = sum(
        fraction * charge * charge for fraction, charge in zip(fractions, charges, strict=True)
    )
    kappa = math.sqrt(
        density
        * AVOGADRO
        * ELEMENTARY_CHARGE**2
        * ionic_strength_sum
        / (VACUUM_PERMITTIVITY * relative_permittivity * BOLTZMANN * temperature)
    )
    weighted_finite_size = 0.0
    for fraction, charge, sigma in zip(fractions, charges, sigmas_angstrom, strict=True):
        if charge == 0:
            continue
        q = kappa * 0.88 * sigma * 1.0e-10
        chi = 3.0 / q**3 * (1.5 + math.log1p(q) - 2.0 * (1.0 + q) + 0.5 * (1.0 + q) ** 2)
        weighted_finite_size += fraction * charge * charge * chi
    return (
        -kappa
        * ELEMENTARY_CHARGE**2
        * weighted_finite_size
        / (12.0 * math.pi * VACUUM_PERMITTIVITY * relative_permittivity * BOLTZMANN * temperature)
    )


def test_source_backed_debye_huckel_matches_independent_si_oracle() -> None:
    result = _electrolyte_model().evaluate(
        temperature=298.15 * unit_registry.kelvin,
        molar_density=40.0 * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=(0.9, 0.05, 0.05),
    )

    assert result.debye_huckel == pytest.approx(
        _debye_huckel_oracle(298.15, 40.0), rel=2.0e-11, abs=2.0e-13
    )
    assert result.debye_huckel < 0.0


def test_debye_huckel_pressure_uses_cppad_density_derivative() -> None:
    model = _electrolyte_model()
    temperature = 298.15
    density = 40.0
    fractions = (0.9, 0.05, 0.05)
    step = density * 2.0e-4

    def residual_at(rho: float) -> float:
        return model.evaluate(
            temperature=temperature * unit_registry.kelvin,
            molar_density=rho * unit_registry.mole / unit_registry.meter**3,
            mole_fractions=fractions,
        ).residual_helmholtz

    derivative = (
        -residual_at(density + 2.0 * step)
        + 8.0 * residual_at(density + step)
        - 8.0 * residual_at(density - step)
        + residual_at(density - 2.0 * step)
    ) / (12.0 * step)
    result = model.evaluate(
        temperature=temperature * unit_registry.kelvin,
        molar_density=density * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=fractions,
    )

    assert result.pressure.to("pascal").magnitude == pytest.approx(
        (1.0 + density * derivative) * density * GAS_CONSTANT * temperature,
        rel=5.0e-7,
    )


def test_ionic_topology_has_exact_zero_charge_reduction() -> None:
    result = _electrolyte_model().evaluate(
        temperature=298.15 * unit_registry.kelvin,
        molar_density=40.0 * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=(1.0, 0.0, 0.0),
    )

    assert result.debye_huckel == 0.0


def test_debye_huckel_approaches_zero_from_below_at_infinite_dilution() -> None:
    model = _electrolyte_model()

    values = []
    for ion_fraction_each in (1.0e-8, 1.0e-12):
        result = model.evaluate(
            temperature=298.15 * unit_registry.kelvin,
            molar_density=40.0 * unit_registry.mole / unit_registry.meter**3,
            mole_fractions=(
                1.0 - 2.0 * ion_fraction_each,
                ion_fraction_each,
                ion_fraction_each,
            ),
        )
        values.append(result.debye_huckel)

    assert values[0] < 0.0
    assert values[1] < 0.0
    assert abs(values[1]) < abs(values[0]) * 1.0e-4


def test_ionic_states_enforce_charge_balance_and_source_domain() -> None:
    model = _electrolyte_model()
    density = 40.0 * unit_registry.mole / unit_registry.meter**3
    with pytest.raises(EosError):
        model.evaluate(
            temperature=298.15 * unit_registry.kelvin,
            molar_density=density,
            mole_fractions=(0.9, 0.1, 0.0),
        )
    with pytest.raises(EosError):
        model.evaluate(
            temperature=300.0 * unit_registry.kelvin,
            molar_density=density,
            mole_fractions=(0.9, 0.05, 0.05),
        )
    with pytest.raises(EosError):
        model.evaluate(
            temperature=298.15 * unit_registry.kelvin,
            molar_density=density,
            mole_fractions=(0.6, 0.2, 0.2),
        )


def test_neutral_topology_has_zero_debye_huckel_contribution() -> None:
    result = EPCSAFT(
        ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1).select(
            ("methane", "ethane")
        )
    ).evaluate(
        temperature=250.0 * unit_registry.kelvin,
        molar_density=1000.0 * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=(0.4, 0.6),
    )

    assert result.debye_huckel == 0.0

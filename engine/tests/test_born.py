import math

import pytest

from epcsaft import EPCSAFT, ParameterBundle, unit_registry


ELEMENTARY_CHARGE = 1.602176634e-19
VACUUM_PERMITTIVITY = 8.8541878128e-12
BOLTZMANN = 1.380649e-23


def _model() -> EPCSAFT:
    return EPCSAFT(
        ParameterBundle.from_catalog("figiel-2025-water-sodium-chloride", version=1).select(
            ("water", "sodium-cation", "chloride-anion")
        )
    )


def _born_oracle(temperature: float) -> float:
    fractions = (0.9, 0.05, 0.05)
    charges = (0, 1, -1)
    born_diameters_angstrom = (0.0, 3.445, 4.100)
    solvation_factors = (1.5, 1.0, 1.0)
    mixed_solvation = sum(
        fraction * factor for fraction, factor in zip(fractions, solvation_factors, strict=True)
    ) / sum(fractions)
    bulk_permittivity = 78.09 / (1.0 + 7.01 * 0.1)
    ionic_permittivity = 8.0

    radial_sum = 0.0
    for fraction, charge, diameter_angstrom in zip(
        fractions, charges, born_diameters_angstrom, strict=True
    ):
        if charge == 0:
            continue
        born_diameter = diameter_angstrom * 1.0e-10
        shell_diameter = born_diameter * (1.0 + (mixed_solvation - 1.0) / abs(charge))
        radial_sum += (
            fraction
            * charge
            * charge
            * (
                (1.0 - 1.0 / bulk_permittivity) / shell_diameter
                + (1.0 - 1.0 / ionic_permittivity) * (1.0 / born_diameter - 1.0 / shell_diameter)
            )
        )
    return (
        -(ELEMENTARY_CHARGE**2)
        * radial_sum
        / (4.0 * math.pi * VACUUM_PERMITTIVITY * BOLTZMANN * temperature)
    )


def test_corrected_born_ssm_ds_matches_independent_two_region_oracle() -> None:
    result = _model().evaluate(
        temperature=298.15 * unit_registry.kelvin,
        molar_density=40.0 * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=(0.9, 0.05, 0.05),
    )

    assert result.born_ssm_ds == pytest.approx(_born_oracle(298.15), rel=2.0e-12, abs=2.0e-12)
    assert result.born_ssm_ds < 0.0
    assert result.residual_helmholtz == pytest.approx(
        result.hard_chain
        + result.dispersion
        + result.association
        + result.debye_huckel
        + result.born_ssm_ds
    )


def test_born_has_zero_ion_reduction_and_no_density_dependence() -> None:
    model = _model()
    state = {
        "temperature": 298.15 * unit_registry.kelvin,
        "mole_fractions": (0.9, 0.05, 0.05),
    }
    dilute = model.evaluate(
        **state,
        molar_density=20.0 * unit_registry.mole / unit_registry.meter**3,
    )
    dense = model.evaluate(
        **state,
        molar_density=80.0 * unit_registry.mole / unit_registry.meter**3,
    )
    solvent = model.evaluate(
        temperature=298.15 * unit_registry.kelvin,
        molar_density=40.0 * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=(1.0, 0.0, 0.0),
    )

    assert dilute.born_ssm_ds == pytest.approx(dense.born_ssm_ds, abs=1.0e-13)
    assert solvent.born_ssm_ds == 0.0


def test_neutral_topology_has_zero_born_contribution() -> None:
    result = EPCSAFT(
        ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1).select(
            ("methane", "ethane")
        )
    ).evaluate(
        temperature=250.0 * unit_registry.kelvin,
        molar_density=1000.0 * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=(0.4, 0.6),
    )

    assert result.born_ssm_ds == 0.0

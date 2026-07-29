from epcsaft.parameters import _ParameterBundle
from decimal import Decimal, localcontext

import pytest

from epcsaft import Mixture, unit_registry


ELEMENTARY_CHARGE = 1.602176634e-19
VACUUM_PERMITTIVITY = 8.8541878128e-12
BOLTZMANN = 1.380649e-23


def _model() -> Mixture:
    return Mixture(
        _ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1).select(
            ("water", "sodium-cation", "chloride-anion")
        )
    )


def _ascani_model() -> Mixture:
    return Mixture(
        _ParameterBundle.from_catalog("ascani-2022-case-study-2-electrolyte-lle", version=1).select(
            ("water", "1-butanol", "sodium-cation", "potassium-cation", "chloride-anion")
        )
    )


def _born_oracle(temperature: float) -> float:
    fractions = (Decimal("0.999998"), Decimal("0.000001"), Decimal("0.000001"))
    charges = (0, 1, -1)
    born_diameters_angstrom = (Decimal("0"), Decimal("3.445"), Decimal("4.100"))
    solvation_factors = (Decimal("1.5"), Decimal("1"), Decimal("1"))
    mixed_solvation = sum(
        fraction * factor for fraction, factor in zip(fractions, solvation_factors, strict=True)
    ) / sum(fractions)
    ion_fraction = fractions[1] + fractions[2]
    bulk_permittivity = Decimal("78.09") / (Decimal("1") + Decimal("7.01") * ion_fraction)
    ionic_permittivity = Decimal("8")

    radial_sum = Decimal("0")
    for fraction, charge, diameter_angstrom in zip(
        fractions, charges, born_diameters_angstrom, strict=True
    ):
        if charge == 0:
            continue
        born_diameter = diameter_angstrom * Decimal("1e-10")
        shell_diameter = born_diameter * (
            Decimal("1") + (mixed_solvation - Decimal("1")) / abs(charge)
        )
        radial_sum += (
            fraction
            * charge
            * charge
            * (
                (Decimal("1") - Decimal("1") / bulk_permittivity) / shell_diameter
                + (Decimal("1") - Decimal("1") / ionic_permittivity)
                * (Decimal("1") / born_diameter - Decimal("1") / shell_diameter)
            )
        )
    with localcontext() as context:
        context.prec = 50
        value = (
            -(Decimal(str(ELEMENTARY_CHARGE)) ** 2)
            * radial_sum
            / (
                Decimal("4")
                * Decimal("3.1415926535897932384626433832795028841971693993751")
                * Decimal(str(VACUUM_PERMITTIVITY))
                * Decimal(str(BOLTZMANN))
                * Decimal(str(temperature))
            )
        )
    return float(value)


def test_corrected_born_ssm_ds_matches_independent_two_region_oracle() -> None:
    result = _model().state(
        T=298.15 * unit_registry.kelvin,
        rho=40000.0 * unit_registry.mole / unit_registry.meter**3,
        x=(0.999998, 0.000001, 0.000001),
    )

    assert result.born == pytest.approx(_born_oracle(298.15), rel=2.0e-12, abs=2.0e-12)
    assert result.born < 0.0
    assert result._residual_helmholtz == pytest.approx(
        result.hard_chain
        + result.dispersion
        + result.association
        + result.debye_huckel
        + result.born
    )


def test_born_has_zero_ion_reduction_and_no_density_dependence() -> None:
    model = _model()
    state = {
        "T": 298.15 * unit_registry.kelvin,
        "x": (0.9, 0.05, 0.05),
    }
    dilute = model.state(
        **state,
        rho=20.0 * unit_registry.mole / unit_registry.meter**3,
    )
    dense = model.state(
        **state,
        rho=80.0 * unit_registry.mole / unit_registry.meter**3,
    )
    solvent = model.state(
        T=298.15 * unit_registry.kelvin,
        rho=40.0 * unit_registry.mole / unit_registry.meter**3,
        x=(1.0, 0.0, 0.0),
    )

    assert dilute.born == pytest.approx(dense.born, abs=1.0e-13)
    assert solvent.born == 0.0


def test_neutral_topology_has_zero_born_contribution() -> None:
    result = Mixture(
        _ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1).select(
            ("methane", "ethane")
        )
    ).state(
        T=250.0 * unit_registry.kelvin,
        rho=1000.0 * unit_registry.mole / unit_registry.meter**3,
        x=(0.4, 0.6),
    )

    assert result.born == 0.0


def test_ascani_original_born_and_combined_dielectric_match_independent_oracle() -> None:
    temperature = Decimal("298.15")
    fractions = tuple(Decimal(value) for value in ("0.74", "0.20", "0.02", "0.01", "0.03"))
    masses = (Decimal("0.01801528"), Decimal("0.0741216"))
    solvent_permittivities = (
        -Decimal("105.2") * Decimal(str(float(temperature))).ln() + Decimal("677.480"),
        -Decimal("0.1077") * temperature + Decimal("49.723"),
    )
    solvent_mass = sum(
        fraction * mass for fraction, mass in zip(fractions[:2], masses, strict=True)
    )
    salt_free_permittivity = (
        sum(
            fraction * mass * permittivity
            for fraction, mass, permittivity in zip(
                fractions[:2], masses, solvent_permittivities, strict=True
            )
        )
        / solvent_mass
    )
    bulk_permittivity = salt_free_permittivity * sum(fractions[:2]) + Decimal("8") * sum(
        fractions[2:]
    )
    charges = (0, 0, 1, 1, -1)
    sigma_angstrom = (
        Decimal("0"),
        Decimal("0"),
        Decimal("2.8232"),
        Decimal("3.3417"),
        Decimal("2.7560"),
    )
    radial_sum = sum(
        fraction * charge * charge / (sigma * Decimal("1e-10"))
        for fraction, charge, sigma in zip(fractions, charges, sigma_angstrom, strict=True)
        if charge
    )
    with localcontext() as context:
        context.prec = 50
        expected = (
            -(Decimal(str(ELEMENTARY_CHARGE)) ** 2)
            * (Decimal("1") - Decimal("1") / bulk_permittivity)
            * radial_sum
            / (
                Decimal("4")
                * Decimal("3.1415926535897932384626433832795028841971693993751")
                * Decimal(str(VACUUM_PERMITTIVITY))
                * Decimal(str(BOLTZMANN))
                * temperature
            )
        )

    result = _ascani_model().state(
        T=float(temperature) * unit_registry.kelvin,
        rho=1000.0 * unit_registry.mole / unit_registry.meter**3,
        x=tuple(float(value) for value in fractions),
    )
    assert result.born == pytest.approx(float(expected), rel=3.0e-12)
    assert result.born < 0.0

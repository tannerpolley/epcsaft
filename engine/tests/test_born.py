from decimal import Decimal, localcontext

import pytest

from epcsaft import EPCSAFT, ParameterBundle, unit_registry


ELEMENTARY_CHARGE = 1.602176634e-19
VACUUM_PERMITTIVITY = 8.8541878128e-12
BOLTZMANN = 1.380649e-23


def _model() -> EPCSAFT:
    return EPCSAFT(
        ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1).select(
            ("water", "sodium-cation", "chloride-anion")
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
    result = _model().evaluate(
        temperature=298.15 * unit_registry.kelvin,
        molar_density=40000.0 * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=(0.999998, 0.000001, 0.000001),
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

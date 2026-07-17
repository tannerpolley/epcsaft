import math

import pytest

from epcsaft import EPCSAFT, ParameterBundle, unit_registry
from epcsaft.records import (
    AssociationParameterRecord,
    ComponentRecord,
    SingleParameterRecord,
    SiteRecord,
    SourceRecord,
    ValidityDomain,
)


AVOGADRO_ANGSTROM = 6.02214076e23 / 1.0e30
GAS_CONSTANT = 8.31446261815324
PROVENANCE = {
    "source_id": "association-oracle",
    "locator": "independent pure 2B analytic fixture",
    "domain_id": "fixture-domain",
}


def _pure_2b_model() -> EPCSAFT:
    component = "water"
    singles = (
        SingleParameterRecord(
            "water-m", component, "segment_count", 1.0, **PROVENANCE
        ),
        SingleParameterRecord(
            "water-sigma",
            component,
            "segment_diameter",
            3.0 * unit_registry.angstrom,
            **PROVENANCE,
        ),
        SingleParameterRecord(
            "water-epsilon",
            component,
            "dispersion_energy_over_k",
            200.0 * unit_registry.kelvin,
            **PROVENANCE,
        ),
    )
    sites = tuple(
        SiteRecord(
            f"water-site-{site}", component, site, site, 1, **PROVENANCE
        )
        for site in ("a", "b")
    )
    associations = tuple(
        AssociationParameterRecord(
            record_id,
            component,
            "a",
            component,
            "b",
            family,
            value,
            **PROVENANCE,
        )
        for record_id, family, value in (
            (
                "water-ab-energy",
                "association_energy_over_k",
                2400.0 * unit_registry.kelvin,
            ),
            ("water-ab-volume", "association_volume", 0.03),
        )
    )
    parameters = ParameterBundle.from_records(
        bundle_id="association-oracle",
        bundle_version=1,
        purpose="package-test-fixture",
        sources=(
            SourceRecord(
                "association-oracle",
                "Independent pure 2B analytic fixture",
                "package test oracle",
            ),
        ),
        domains=(ValidityDomain("fixture-domain", "unknown"),),
        components=(ComponentRecord(component),),
        singles=singles,
        sites=sites,
        associations=associations,
    ).select((component,))
    return EPCSAFT(parameters)


def _pure_2b_association(temperature: float, molar_density: float) -> float:
    sigma = 3.0
    diameter = sigma * (1.0 - 0.12 * math.exp(-3.0 * 200.0 / temperature))
    number_density = molar_density * AVOGADRO_ANGSTROM
    eta = math.pi * number_density * diameter**3 / 6.0
    contact = (
        1.0 / (1.0 - eta)
        + 1.5 * diameter * number_density * math.pi * diameter**2 / 6.0
        / (1.0 - eta) ** 2
        + 0.5
        * diameter**2
        * (number_density * math.pi * diameter**2 / 6.0) ** 2
        / (1.0 - eta) ** 3
    )
    delta = contact * math.expm1(2400.0 / temperature) * sigma**3 * 0.03
    site_fraction = 2.0 / (
        1.0 + math.sqrt(1.0 + 4.0 * number_density * delta)
    )
    return 2.0 * (
        math.log(site_fraction) - 0.5 * site_fraction + 0.5
    )


def test_pure_2b_association_matches_independent_closed_form() -> None:
    model = _pure_2b_model()
    result = model.evaluate(
        temperature=300.0 * unit_registry.kelvin,
        molar_density=10000.0 * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=(1.0,),
    )

    assert result.association == pytest.approx(
        _pure_2b_association(300.0, 10000.0), rel=2.0e-11, abs=2.0e-12
    )
    assert result.association < 0.0
    assert result.residual_helmholtz == pytest.approx(
        result.hard_chain + result.dispersion + result.association
    )


def test_association_pressure_uses_total_density_derivative() -> None:
    model = _pure_2b_model()
    temperature = 300.0
    density = 10000.0
    step = density * 2.0e-4

    def residual_at(rho: float) -> float:
        return model.evaluate(
            temperature=temperature * unit_registry.kelvin,
            molar_density=rho * unit_registry.mole / unit_registry.meter**3,
            mole_fractions=(1.0,),
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
        mole_fractions=(1.0,),
    )
    expected = (1.0 + density * derivative) * density * GAS_CONSTANT * temperature

    assert result.pressure.to("pascal").magnitude == pytest.approx(
        expected, rel=3.0e-7
    )


def test_association_has_the_correct_zero_density_limit() -> None:
    result = _pure_2b_model().evaluate(
        temperature=300.0 * unit_registry.kelvin,
        molar_density=1.0e-8 * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=(1.0,),
    )

    assert abs(result.association) < 1.0e-8

import math
from itertools import combinations

import pytest

from epcsaft import EPCSAFT, ParameterBundle, unit_registry
from epcsaft.records import (
    ComponentRecord,
    ConstantCorrelation,
    ModelParameterRecord,
    PairParameterRecord,
    SingleParameterRecord,
    SourceRecord,
    ValidityDomain,
)


AVOGADRO = 6.02214076e23
BOLTZMANN = 1.380649e-23
GAS_CONSTANT = 8.31446261815324
PROVENANCE = {
    "source_id": "permittivity-fixture",
    "locator": "independent equation fixture",
    "domain_id": "permittivity-domain",
}


def _single(record_id: str, component_id: str, family: str, value: object):
    return SingleParameterRecord(record_id, component_id, family, value, **PROVENANCE)


def _base_component(component_id: str, charge: int, molar_mass: float):
    return (
        _single(
            f"{component_id}-molar-mass",
            component_id,
            "molar_mass",
            molar_mass * unit_registry.kilogram / unit_registry.mole,
        ),
        _single(f"{component_id}-m", component_id, "segment_count", 1),
        _single(f"{component_id}-sigma", component_id, "segment_diameter", 3 * unit_registry.angstrom),
        _single(f"{component_id}-epsilon", component_id, "dispersion_energy_over_k", 200 * unit_registry.kelvin),
        _single(f"{component_id}-charge", component_id, "charge_number", charge),
        _single(f"{component_id}-solvation", component_id, "solvation_factor", 1),
    )


def _parameters(
    kind: str,
    *,
    extra_models: tuple[ModelParameterRecord, ...] = (),
    include_second_solvent: bool = False,
) -> ParameterBundle:
    component_ids = (
        ("water", "methanol", "sodium", "chloride")
        if include_second_solvent
        else ("water", "sodium", "chloride")
    )
    singles = [
        *_base_component("water", 0, 0.01801528),
        *_base_component("sodium", 1, 0.02298977),
        *_base_component("chloride", -1, 0.03545),
        _single("sodium-born", "sodium", "born_diameter", 2.8232 * unit_registry.angstrom),
        _single("chloride-born", "chloride", "born_diameter", 2.7560 * unit_registry.angstrom),
    ]
    if include_second_solvent:
        singles.extend(_base_component("methanol", 0, 0.03204))
    models = [
        ModelParameterRecord(
            "ionic-region-permittivity",
            "ionic_region_relative_permittivity",
            8,
            **PROVENANCE,
        )
    ]
    correlations = []
    if kind in {"figiel-composition", "zuber-ion-specific"}:
        correlations.append(
            ConstantCorrelation(
                "water-relative-permittivity",
                "water",
                "relative_permittivity",
                78.09,
                **PROVENANCE,
            )
        )
        if include_second_solvent:
            correlations.append(
                ConstantCorrelation(
                    "methanol-relative-permittivity",
                    "methanol",
                    "relative_permittivity",
                    32.63,
                    **PROVENANCE,
                )
            )
    if kind == "figiel-composition":
        models.append(
            ModelParameterRecord(
                "figiel-suppression",
                "dielectric_ion_suppression_coefficient",
                7.01,
                **PROVENANCE,
            )
        )
    elif kind == "schreckenberg-density":
        singles.extend(
            (
                _single(
                    "water-schreckenberg-volume",
                    "water",
                    "schreckenberg_dielectric_volume",
                    0.3777 * unit_registry.decimeter**3 / unit_registry.mole,
                ),
                _single(
                    "water-schreckenberg-temperature",
                    "water",
                    "schreckenberg_dielectric_temperature",
                    1403 * unit_registry.kelvin,
                ),
            )
        )
    elif kind == "zuber-ion-specific":
        singles.extend(
            (
                _single("sodium-zuber-alpha", "sodium", "zuber_ion_suppression_coefficient", 3.62),
                _single("chloride-zuber-alpha", "chloride", "zuber_ion_suppression_coefficient", 6.80),
            )
        )
    elif kind == "rueben-perturbation":
        rueben = {
            "water": (5.199, 0.0, 0.1276),
            "sodium": (5.407, 0.0, 0.215),
            "chloride": (2.563, 0.0, 0.138),
        }
        for component_id, (dipole, polarizability, correlation) in rueben.items():
            singles.extend(
                (
                    _single(
                        f"{component_id}-rueben-dipole",
                        component_id,
                        "rueben_dipole_scaling",
                        dipole * unit_registry.debye**2,
                    ),
                    _single(
                        f"{component_id}-rueben-polarizability",
                        component_id,
                        "rueben_polarizability_scaling",
                        polarizability * unit_registry.angstrom**3,
                    ),
                    _single(
                        f"{component_id}-rueben-correlation",
                        component_id,
                        "rueben_correlation_integral_parameter",
                        correlation,
                    ),
                )
            )
    else:
        raise AssertionError(kind)
    models.extend(extra_models)

    pairs = tuple(
        PairParameterRecord(f"{a}-{b}-kij", a, b, "k_ij", 0, **PROVENANCE)
        for a, b in combinations(component_ids, 2)
    )
    return ParameterBundle.from_records(
        bundle_id="permittivity-model-fixture",
        bundle_version=1,
        purpose="package-test-fixture",
        sources=(SourceRecord("permittivity-fixture", "Primary-source equations", "test fixture"),),
        domains=(ValidityDomain("permittivity-domain", "unknown"),),
        components=tuple(ComponentRecord(component_id) for component_id in component_ids),
        singles=tuple(singles),
        pairs=pairs,
        correlations=tuple(correlations),
        models=tuple(models),
    ).select(component_ids)


def _model(kind: str) -> EPCSAFT:
    return EPCSAFT(_parameters(kind))


def _rueben_oracle(temperature: float, molar_density: float, fractions: tuple[float, ...]) -> float:
    dipole = (5.199, 5.407, 2.563)
    polarizability = (0.0, 0.0, 0.0)
    correlation = (0.1276, 0.215, 0.138)
    number_density_per_angstrom_cubed = molar_density * AVOGADRO * 1.0e-30
    beta = 1.0 / (BOLTZMANN * temperature)
    y_effective = 4.0 * math.pi / 9.0 * math.fsum(
        number_density_per_angstrom_cubed
        * fraction
        * (beta * dipole_scaling * 1.0e-19 + 3.0 * polarizability_scaling)
        for fraction, dipole_scaling, polarizability_scaling in zip(
            fractions, dipole, polarizability, strict=True
        )
    )
    integral = 1.0 + math.fsum(
        fraction * parameter for fraction, parameter in zip(fractions, correlation, strict=True)
    ) * (math.exp(-y_effective) - 1.0)
    return 1.0 + 3.0 * y_effective * (
        1.0
        + y_effective
        + ((17.0 / 16.0) * integral - 1.0) * y_effective**2
    )


@pytest.mark.parametrize(
    ("kind", "expected"),
    (
        ("figiel-composition", 78.09 / (1.0 + 7.01 * 0.1)),
        (
            "schreckenberg-density",
            1.0 + 55_000.0 * 0.9 * 0.3777e-3 * (1403.0 / 298.15 - 1.0),
        ),
        ("zuber-ion-specific", 78.09 / (1.0 + 3.62 * 0.05 + 6.80 * 0.05)),
        ("rueben-perturbation", _rueben_oracle(298.15, 55_000.0, (0.9, 0.05, 0.05))),
    ),
)
def test_bulk_permittivity_matches_independent_source_equation(kind: str, expected: float) -> None:
    model = _model(kind)
    result = model.evaluate(
        temperature=298.15 * unit_registry.kelvin,
        molar_density=55_000 * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=(0.9, 0.05, 0.05),
    )

    assert model.permittivity_model == kind
    assert result.bulk_relative_permittivity == pytest.approx(expected, rel=2.0e-12)


@pytest.mark.parametrize("kind", ("schreckenberg-density", "rueben-perturbation"))
def test_density_models_approach_vacuum_permittivity_at_zero_density(kind: str) -> None:
    model = _model(kind)
    result = model.evaluate(
        temperature=298.15 * unit_registry.kelvin,
        molar_density=1.0e-9 * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=(1.0 - 2.0e-12, 1.0e-12, 1.0e-12),
    )

    assert result.bulk_relative_permittivity == pytest.approx(1.0, abs=2.0e-12)


@pytest.mark.parametrize(
    "kind",
    (
        "figiel-composition",
        "schreckenberg-density",
        "zuber-ion-specific",
        "rueben-perturbation",
    ),
)
def test_permittivity_density_derivatives_remain_in_common_cppad_helmholtz_owner(
    kind: str,
) -> None:
    model = _model(kind)
    temperature = 298.15
    density = 40_000.0
    step = density * 2.0e-4

    def residual_at(rho: float) -> float:
        return model.evaluate(
            temperature=temperature * unit_registry.kelvin,
            molar_density=rho * unit_registry.mole / unit_registry.meter**3,
            mole_fractions=(0.9, 0.05, 0.05),
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
        mole_fractions=(0.9, 0.05, 0.05),
    )

    assert result.pressure.to("pascal").magnitude == pytest.approx(
        (1.0 + density * derivative) * density * GAS_CONSTANT * temperature,
        rel=8.0e-7,
    )


@pytest.mark.parametrize(
    "kind",
    (
        "figiel-composition",
        "schreckenberg-density",
        "zuber-ion-specific",
        "rueben-perturbation",
    ),
)
def test_permittivity_model_closes_a_single_pressure_root_with_exact_slope(kind: str) -> None:
    model = _model(kind)
    fractions = (0.826112, 0.086944, 0.086944)
    result = model.evaluate(
        temperature=298.15 * unit_registry.kelvin,
        pressure=2508.0 * unit_registry.pascal,
        mole_fractions=fractions,
    )

    assert result.density_diagnostics is not None
    assert result.density_diagnostics.branch == "single"
    assert result.density_diagnostics.stable is True
    assert abs(result.density_diagnostics.pressure_residual.magnitude) < 1.0e-5
    density = result.molar_density.to("mole / meter ** 3").magnitude
    step = density * 1.0e-4

    def pressure_at(rho: float) -> float:
        return (
            model.evaluate(
                temperature=298.15 * unit_registry.kelvin,
                molar_density=rho * unit_registry.mole / unit_registry.meter**3,
                mole_fractions=fractions,
            )
            .pressure.to("pascal")
            .magnitude
        )

    finite_difference = (pressure_at(density + step) - pressure_at(density - step)) / (
        2.0 * step
    )
    assert result.density_diagnostics.pressure_density_derivative.magnitude == pytest.approx(
        finite_difference,
        rel=8.0e-7,
    )


def test_mixed_permittivity_topologies_fail_during_model_resolution() -> None:
    mixed = _parameters(
        "schreckenberg-density",
        extra_models=(
            ModelParameterRecord(
                "figiel-suppression",
                "dielectric_ion_suppression_coefficient",
                7.01,
                **PROVENANCE,
            ),
        ),
    )

    with pytest.raises(ValueError, match="missing-or-mixed-topology"):
        EPCSAFT(mixed)


def test_zuber_mixed_solvent_selection_requires_source_volume_data() -> None:
    with pytest.raises(ValueError, match="zuber-mixed-solvent-volume-data-required"):
        EPCSAFT(_parameters("zuber-ion-specific", include_second_solvent=True))


def test_neutral_model_reports_permittivity_as_not_applicable() -> None:
    model = EPCSAFT(
        ParameterBundle.from_catalog("gross-2001-propane", version=1).select(("propane",))
    )
    result = model.evaluate(
        temperature=250.0 * unit_registry.kelvin,
        molar_density=1000.0 * unit_registry.mole / unit_registry.meter**3,
        mole_fractions=(1.0,),
    )

    assert model.permittivity_model == "none"
    assert result.bulk_relative_permittivity is None

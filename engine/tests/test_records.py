from dataclasses import FrozenInstanceError
from decimal import Decimal
from types import SimpleNamespace

import pytest

import epcsaft
from epcsaft import ParameterError, unit_registry
from epcsaft.records import (
    AssociationParameterRecord,
    ComponentRecord,
    ConstantCorrelation,
    ConstantPlusSumOfExponentialsCorrelation,
    ExponentialTerm,
    ModelParameterRecord,
    PairParameterRecord,
    SingleParameterRecord,
    SiteRecord,
    SourceRecord,
    ValidityDomain,
)
from pint import UnitRegistry


PROVENANCE = {
    "source_id": "figiel-yu-held-2025",
    "locator": "Table 2, water column",
    "domain_id": "figiel-reported-conditions",
}


def test_records_are_immutable_and_identifiers_are_exact() -> None:
    assert epcsaft.__all__ == (
        "ParameterBundle",
        "ParameterSet",
        "ParameterError",
        "unit_registry",
        "EPCSAFT",
        "EosResult",
        "EosError",
    )

    source = SourceRecord(
        source_id="gross-sadowski-2001",
        citation="Gross and Sadowski (2001)",
        use_basis="primary-paper factual parameters",
        doi="10.1021/ie0003887",
    )
    domain = ValidityDomain(
        domain_id="gross-methane-fit-range",
        kind="fit-range",
        temperature_min=97 * unit_registry.kelvin,
        temperature_max=300 * unit_registry.kelvin,
    )
    one_sided_domain = ValidityDomain(
        domain_id="sourced-upper-bound",
        kind="reported-conditions",
        pressure_max=2 * unit_registry.bar,
    )
    component = ComponentRecord("methane", name="methane", aliases=("CH4",))

    assert source.source_id == "gross-sadowski-2001"
    assert domain.temperature_min.magnitude == Decimal("97")
    assert one_sided_domain.pressure_max.magnitude == Decimal("200000")
    assert component.aliases == ("CH4",)
    with pytest.raises(FrozenInstanceError):
        component.name = "Methane"

    for bad_id in ("Methane", "methane_1", " methane", "methane ", ""):
        with pytest.raises(ParameterError) as exc_info:
            ComponentRecord(bad_id)
        assert exc_info.value.code == "invalid-id"
        assert exc_info.value.context == (("field", "component_id"), ("value", bad_id))


def test_field_native_units_and_exact_integers() -> None:
    foreign = UnitRegistry()
    mass = SingleParameterRecord(
        record_id="gross-methane-molar-mass",
        component_id="methane",
        family="molar_mass",
        value=16.043 * foreign.gram / foreign.mole,
        source_id="gross-sadowski-2001",
        locator="Table 2, methane row, M column",
        domain_id="gross-methane-fit-range",
    )
    charge = SingleParameterRecord(
        record_id="figiel-potassium-charge",
        component_id="potassium-cation",
        family="charge_number",
        value=1,
        **PROVENANCE,
    )
    site = SiteRecord(
        record_id="figiel-water-site-a",
        component_id="water",
        site_id="a",
        site_class="a",
        multiplicity=1,
        **PROVENANCE,
    )

    assert mass.value.units == unit_registry.kilogram / unit_registry.mole
    assert mass.value.magnitude == Decimal("0.016043")
    assert mass.value._REGISTRY is unit_registry
    assert charge.value == 1
    assert site.multiplicity == 1

    for invalid_charge in (True, 1.0, Decimal("1")):
        with pytest.raises(ParameterError, match="invalid-value"):
            SingleParameterRecord(
                record_id="figiel-invalid-charge",
                component_id="potassium-cation",
                family="charge_number",
                value=invalid_charge,
                **PROVENANCE,
            )
    for invalid_multiplicity in (True, 0, 1.0):
        with pytest.raises(ParameterError, match="invalid-value"):
            SiteRecord(
                record_id="figiel-invalid-site",
                component_id="water",
                site_id="a",
                site_class="a",
                multiplicity=invalid_multiplicity,
                **PROVENANCE,
            )
    with pytest.raises(ParameterError, match="invalid-unit"):
        SingleParameterRecord(
            record_id="gross-unitless-mass",
            component_id="methane",
            family="molar_mass",
            value=0.016043,
            source_id="gross-sadowski-2001",
            locator="Table 2, methane row, M column",
            domain_id="gross-methane-fit-range",
        )

    fake_quantity = SimpleNamespace(magnitude=1, units="angstrom")
    with pytest.raises(ParameterError, match="invalid-unit"):
        SingleParameterRecord(
            record_id="figiel-fake-born-diameter",
            component_id="potassium-cation",
            family="born_diameter",
            value=fake_quantity,
            **PROVENANCE,
        )
    with pytest.raises(ParameterError, match="invalid-unit"):
        ExponentialTerm(
            amplitude=fake_quantity,
            exponent_coefficient=1 / foreign.kelvin,
        )


def test_sparse_closed_records_correlations_and_error_contract() -> None:
    foreign = UnitRegistry()
    pair = PairParameterRecord(
        record_id="figiel-potassium-iodide-kij",
        component_id_a="potassium-cation",
        component_id_b="iodide-anion",
        family="k_ij",
        value=0,
        **PROVENANCE,
    )
    association = AssociationParameterRecord(
        record_id="figiel-water-association-energy",
        component_id_a="water",
        site_id_a="a",
        component_id_b="water",
        site_id_b="b",
        family="association_energy_over_k",
        value=2425.7 * foreign.kelvin,
        **PROVENANCE,
    )
    constant = ConstantCorrelation(
        record_id="figiel-water-relative-permittivity",
        component_id="water",
        family="relative_permittivity",
        constant=78.09,
        **PROVENANCE,
    )
    water_diameter = ConstantPlusSumOfExponentialsCorrelation(
        record_id="figiel-water-segment-diameter",
        component_id="water",
        family="segment_diameter",
        constant=2.7927 * foreign.angstrom,
        terms=(
            ExponentialTerm(10.11 * foreign.angstrom, -0.01775 / foreign.kelvin),
            ExponentialTerm(-1.417 * foreign.angstrom, -0.01146 / foreign.kelvin),
        ),
        **PROVENANCE,
    )

    assert pair.component_id_a == "iodide-anion"
    assert pair.component_id_b == "potassium-cation"
    assert pair.value == Decimal("0")
    assert association.value.units == unit_registry.kelvin
    assert constant.form == "constant"
    assert constant.constant == Decimal("78.09")
    assert water_diameter.form == "constant-plus-sum-of-exponentials"
    assert isinstance(water_diameter.terms, tuple)
    assert water_diameter.terms[0].amplitude.units == unit_registry.angstrom
    assert water_diameter.terms[0].exponent_coefficient.units == 1 / unit_registry.kelvin
    assert not hasattr(ComponentRecord("methane"), "charge_number")

    with pytest.raises(ParameterError, match="invalid-value"):
        SingleParameterRecord(
            record_id="unsupported-family",
            component_id="water",
            family="electrolyte_model",
            value=1,
            **PROVENANCE,
        )
    error = ParameterError(
        "invalid-value",
        (("field", "family"), ("value", "electrolyte_model")),
    )
    assert isinstance(error, ValueError)
    assert error.code == "invalid-value"
    assert error.context == (("field", "family"), ("value", "electrolyte_model"))
    assert str(error) == "invalid-value (field=family, value=electrolyte_model)"


def test_model_parameter_records_are_closed_finite_and_source_backed() -> None:
    ionic_permittivity = ModelParameterRecord(
        record_id="figiel-ionic-region-permittivity",
        family="ionic_region_relative_permittivity",
        value=8,
        **PROVENANCE,
    )
    suppression = ModelParameterRecord(
        record_id="figiel-dielectric-suppression",
        family="dielectric_ion_suppression_coefficient",
        value=Decimal("7.01"),
        **PROVENANCE,
    )

    assert ionic_permittivity.value == Decimal("8")
    assert suppression.value == Decimal("7.01")

    for family, value in (
        ("electrolyte_model", 1),
        ("ionic_region_relative_permittivity", Decimal("NaN")),
    ):
        with pytest.raises(ParameterError, match="invalid-value"):
            ModelParameterRecord(
                record_id="invalid-model-parameter",
                family=family,
                value=value,
                **PROVENANCE,
            )

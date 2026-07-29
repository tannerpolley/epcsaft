from epcsaft.parameters import _ParameterBundle
import csv
from dataclasses import FrozenInstanceError
from decimal import Decimal
import json
from pathlib import Path
import shutil

import pytest
from pint import UnitRegistry

import epcsaft
from epcsaft import BundlePurpose, ParameterError, Parameters, unit_registry
from epcsaft._canonical import _bundle_json
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


PROVENANCE = {
    "source_id": "source-a",
    "locator": "Table 1, test cell",
    "domain_id": "domain-a",
}
COMPONENT_COLUMNS = ("component_id", "name", "cas", "formula", "aliases")
MODEL_COLUMNS = (
    "record_id",
    "family",
    "value",
    "unit",
    "source_id",
    "locator",
    "domain_id",
)


def _bundle_records(registry: UnitRegistry = unit_registry, *, reverse: bool = False):
    source = SourceRecord("source-a", "Source A", "package-test evidence")
    domain = ValidityDomain("domain-a", "unknown")
    components = [
        ComponentRecord("methane", name="methane", formula="CH4", aliases=("CH4",)),
        ComponentRecord("ethane", name="ethane", formula="C2H6"),
    ]
    singles = [
        SingleParameterRecord(
            "methane-molar-mass",
            "methane",
            "molar_mass",
            Decimal("0.016") * registry.kilogram / registry.mole,
            **PROVENANCE,
        ),
        SingleParameterRecord(
            "ethane-segment-count",
            "ethane",
            "segment_count",
            Decimal("1.6"),
            **PROVENANCE,
        ),
    ]
    pair = PairParameterRecord(
        "methane-ethane-kij",
        "methane",
        "ethane",
        "k_ij",
        Decimal("-0"),
        **PROVENANCE,
    )
    if reverse:
        components.reverse()
        singles.reverse()
    return source, domain, components, singles, pair


def _from_records(
    registry: UnitRegistry = unit_registry,
    *,
    reverse: bool = False,
    ionic_permittivity: object = Decimal("8"),
    purpose: BundlePurpose = "package-test-fixture",
):
    source, domain, components, singles, pair = _bundle_records(registry, reverse=reverse)
    models = [
        ModelParameterRecord(
            "figiel-ionic-region-permittivity",
            "ionic_region_relative_permittivity",
            ionic_permittivity,
            **PROVENANCE,
        ),
        ModelParameterRecord(
            "figiel-dielectric-suppression",
            "ion_fraction_suppression_coefficient",
            Decimal("7.01"),
            **PROVENANCE,
        ),
    ]
    if reverse:
        models.reverse()
    return _ParameterBundle.from_records(
        bundle_id="test-bundle",
        bundle_version=1,
        purpose=purpose,
        sources=(source,),
        domains=(domain,),
        components=components,
        singles=singles,
        pairs=(pair,),
        models=models,
    )


def _advanced_records():
    sites = (
        SiteRecord("water-site-a", "water", "a", "a", 1, **PROVENANCE),
        SiteRecord("water-site-b", "water", "b", "b", 1, **PROVENANCE),
    )
    association = AssociationParameterRecord(
        "water-association-energy",
        "water",
        "a",
        "water",
        "b",
        "association_energy_over_k",
        1000 * unit_registry.kelvin,
        **PROVENANCE,
    )
    correlations = (
        ConstantCorrelation(
            "water-relative-permittivity",
            "water",
            "relative_permittivity",
            Decimal("78.09"),
            **PROVENANCE,
        ),
        ConstantPlusSumOfExponentialsCorrelation(
            "water-segment-diameter",
            "water",
            "segment_diameter",
            3 * unit_registry.angstrom,
            (ExponentialTerm(1 * unit_registry.angstrom, -0.01 / unit_registry.kelvin),),
            **PROVENANCE,
        ),
    )
    return (
        ComponentRecord("water", name="water", formula="H2O"),
        sites,
        association,
        correlations,
    )


def _parsed_records_bundle():
    source, domain, components, singles, pair = _bundle_records()
    water, sites, association, correlations = _advanced_records()
    return _ParameterBundle.from_records(
        bundle_id="test-bundle",
        bundle_version=1,
        purpose="package-test-fixture",
        sources=(source,),
        domains=(domain,),
        components=(*components, water),
        singles=singles,
        pairs=(pair,),
        sites=sites,
        associations=(association,),
        correlations=correlations,
        models=(
            ModelParameterRecord(
                "figiel-dielectric-suppression",
                "ion_fraction_suppression_coefficient",
                Decimal("7.01"),
                **PROVENANCE,
            ),
            ModelParameterRecord(
                "figiel-ionic-region-permittivity",
                "ionic_region_relative_permittivity",
                Decimal("8"),
                **PROVENANCE,
            ),
        ),
    )


def _write_bundle(
    path: Path,
    *,
    component_columns: tuple[str, ...] = COMPONENT_COLUMNS,
    model_columns: tuple[str, ...] = MODEL_COLUMNS,
    duplicate_pair: bool = False,
    extra_file: bool = False,
    bad_correlation: bool = False,
) -> None:
    path.mkdir()
    (path / "bundle.toml").write_text(
        """schema_version = 1
bundle_id = "test-bundle"
bundle_version = 1
purpose = "package-test-fixture"
files = ["components.csv", "single.csv", "pair.csv", "sites.csv", "association.csv", "model.csv", "correlations.toml"]

[[sources]]
source_id = "source-a"
citation = "Source A"
use_basis = "package-test evidence"

[[domains]]
domain_id = "domain-a"
kind = "unknown"
""",
        encoding="utf-8",
    )
    with (path / "components.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(component_columns)
        writer.writerow(("methane", "methane", "", "CH4", json.dumps(["CH4"])))
        writer.writerow(("ethane", "ethane", "", "C2H6", ""))
        writer.writerow(("water", "water", "", "H2O", ""))
    with (path / "single.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            (
                "record_id",
                "component_id",
                "family",
                "value",
                "unit",
                "source_id",
                "locator",
                "domain_id",
            )
        )
        writer.writerow(
            (
                "ethane-segment-count",
                "ethane",
                "segment_count",
                "1.6",
                "",
                "source-a",
                "Table 1, test cell",
                "domain-a",
            )
        )
        writer.writerow(
            (
                "methane-molar-mass",
                "methane",
                "molar_mass",
                "16",
                "gram / mole",
                "source-a",
                "Table 1, test cell",
                "domain-a",
            )
        )
    with (path / "pair.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            (
                "record_id",
                "component_id_a",
                "component_id_b",
                "family",
                "value",
                "unit",
                "source_id",
                "locator",
                "domain_id",
            )
        )
        writer.writerow(
            (
                "methane-ethane-kij",
                "methane",
                "ethane",
                "k_ij",
                "0",
                "",
                "source-a",
                "Table 1, test cell",
                "domain-a",
            )
        )
        writer.writerow(())
        if duplicate_pair:
            writer.writerow(
                (
                    "duplicate-kij",
                    "ethane",
                    "methane",
                    "k_ij",
                    "0",
                    "",
                    "source-a",
                    "Table 1, duplicate cell",
                    "domain-a",
                )
            )
    with (path / "sites.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            (
                "record_id",
                "component_id",
                "site_id",
                "site_class",
                "multiplicity",
                "source_id",
                "locator",
                "domain_id",
            )
        )
        for site_id in ("a", "b"):
            writer.writerow(
                (
                    f"water-site-{site_id}",
                    "water",
                    site_id,
                    site_id,
                    "1",
                    "source-a",
                    "Table 1, test cell",
                    "domain-a",
                )
            )
    with (path / "association.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            (
                "record_id",
                "component_id_a",
                "site_id_a",
                "component_id_b",
                "site_id_b",
                "family",
                "value",
                "unit",
                "source_id",
                "locator",
                "domain_id",
            )
        )
        writer.writerow(
            (
                "water-association-energy",
                "water",
                "a",
                "water",
                "b",
                "association_energy_over_k",
                "1000",
                "kelvin",
                "source-a",
                "Table 1, test cell",
                "domain-a",
            )
        )
    with (path / "model.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(model_columns)
        writer.writerow(
            (
                "figiel-dielectric-suppression",
                "ion_fraction_suppression_coefficient",
                "7.01",
                "",
                "source-a",
                "Table 1, test cell",
                "domain-a",
            )
        )
        writer.writerow(
            (
                "figiel-ionic-region-permittivity",
                "ionic_region_relative_permittivity",
                "8",
                "",
                "source-a",
                "Table 1, test cell",
                "domain-a",
            )
        )
    correlation_form = "expression" if bad_correlation else "constant"
    (path / "correlations.toml").write_text(
        f'''[[correlations]]
record_id = "water-relative-permittivity"
component_id = "water"
family = "relative_permittivity"
form = "{correlation_form}"
constant = 78.09
output_unit = "dimensionless"
source_id = "source-a"
locator = "Table 1, test cell"
domain_id = "domain-a"

[[correlations]]
record_id = "water-segment-diameter"
component_id = "water"
family = "segment_diameter"
form = "constant-plus-sum-of-exponentials"
constant = 3
output_unit = "angstrom"
input_unit = "kelvin"
coefficients = [{{ amplitude = 1, exponent_coefficient = -0.01 }}]
source_id = "source-a"
locator = "Table 1, test cell"
domain_id = "domain-a"
''',
        encoding="utf-8",
    )
    if extra_file:
        (path / "notes.txt").write_text("not in inventory", encoding="utf-8")


def _error_code(call) -> str:
    with pytest.raises(ParameterError) as exc_info:
        call()
    return exc_info.value.code


def test_record_lifecycle_order_and_canonical_fingerprints() -> None:
    foreign = UnitRegistry()
    bundle = _from_records()
    equivalent = _from_records(foreign, reverse=True)
    changed_model = _from_records(ionic_permittivity=Decimal("9"))

    assert isinstance(bundle, _ParameterBundle)
    assert tuple(component.component_id for component in bundle.components) == (
        "ethane",
        "methane",
    )
    assert tuple(record.record_id for record in bundle.records) == (
        "ethane-segment-count",
        "figiel-dielectric-suppression",
        "figiel-ionic-region-permittivity",
        "methane-ethane-kij",
        "methane-molar-mass",
    )
    assert all(
        isinstance(value, tuple)
        for value in (bundle.sources, bundle.domains, bundle.components, bundle.records)
    )
    assert bundle.fingerprint == equivalent.fingerprint
    assert bundle.fingerprint != changed_model.fingerprint
    assert bundle.fingerprint.startswith("sha256:")
    assert len(bundle.fingerprint) == 71

    selected = bundle.select(("ethane", "methane"))
    reverse_selected = bundle.select(("methane", "ethane"))
    assert isinstance(selected, Parameters)
    assert selected.component_ids == ("ethane", "methane")
    assert (
        tuple(component.component_id for component in selected.components) == selected.component_ids
    )
    assert selected.bundle_fingerprint == bundle.fingerprint
    assert selected.bundle_purpose == bundle.purpose
    assert selected.fingerprint != reverse_selected.fingerprint
    with pytest.raises(FrozenInstanceError):
        selected.component_ids = ("methane",)

    golden = _ParameterBundle.from_records(
        bundle_id="golden-bundle",
        bundle_version=1,
        purpose="package-test-fixture",
        sources=(SourceRecord("source-a", "Cafe\u0301\r\nLine", "test"),),
        domains=(ValidityDomain("unknown-domain", "unknown"),),
        components=(ComponentRecord("methane", name="Cafe\u0301\r\nLine", aliases=("CH4",)),),
    )
    expected_json = (
        '{"bundle":{"bundle_id":"golden-bundle","bundle_version":1,'
        '"components":[{"aliases":["CH4"],"component_id":"methane",'
        '"name":"Café\\nLine"}],"domains":[{"domain_id":"unknown-domain",'
        '"kind":"unknown"}],"purpose":"package-test-fixture","records":[],'
        '"sources":[{"citation":"Café\\nLine","source_id":"source-a",'
        '"use_basis":"test"}]},"tag":"epcsaft-parameters-v1"}'
    )
    assert _bundle_json(golden) == expected_json
    assert golden.fingerprint == (
        "sha256:2135d1f56ab314e599a56fa34f1c7beff8490c0b8042b961eab0bd174b73233f"
    )


def test_one_strict_path_parser_preserves_sparse_zero_and_matches_records(
    tmp_path: Path,
) -> None:
    bundle_path = tmp_path / "bundle"
    _write_bundle(bundle_path)
    parsed = _ParameterBundle.from_path(bundle_path)
    expected = _parsed_records_bundle()

    assert parsed.fingerprint == expected.fingerprint
    pair = next(record for record in parsed.records if record.record_id == "methane-ethane-kij")
    assert pair.component_id_a == "ethane"
    assert pair.component_id_b == "methane"
    assert pair.value == Decimal("0")
    assert sum(getattr(record, "family", None) == "k_ij" for record in parsed.records) == 1
    assert {type(record) for record in parsed.records} >= {
        SiteRecord,
        AssociationParameterRecord,
        ConstantCorrelation,
        ConstantPlusSumOfExponentialsCorrelation,
        ModelParameterRecord,
    }

    selected = parsed.select(("methane",))
    assert any(isinstance(record, ModelParameterRecord) for record in selected.records)

    bad_columns = tmp_path / "bad-columns"
    _write_bundle(bad_columns, component_columns=COMPONENT_COLUMNS + ("extra",))
    assert _error_code(lambda: _ParameterBundle.from_path(bad_columns)) == "invalid-schema"

    bad_model_columns = tmp_path / "bad-model-columns"
    _write_bundle(bad_model_columns, model_columns=MODEL_COLUMNS + ("extra",))
    assert _error_code(lambda: _ParameterBundle.from_path(bad_model_columns)) == "invalid-schema"

    bad_inventory = tmp_path / "bad-inventory"
    _write_bundle(bad_inventory, extra_file=True)
    assert _error_code(lambda: _ParameterBundle.from_path(bad_inventory)) == "invalid-schema"

    duplicate_pair = tmp_path / "duplicate-pair"
    _write_bundle(duplicate_pair, duplicate_pair=True)
    assert _error_code(lambda: _ParameterBundle.from_path(duplicate_pair)) == "duplicate-record"

    bad_correlation = tmp_path / "bad-correlation"
    _write_bundle(bad_correlation, bad_correlation=True)
    assert _error_code(lambda: _ParameterBundle.from_path(bad_correlation)) == "invalid-correlation"

    catalog_mismatch = tmp_path / "catalog-mismatch"
    _write_bundle(catalog_mismatch)
    assert (
        _error_code(
            lambda: _ParameterBundle._from_catalog_path("requested-bundle", 1, catalog_mismatch)
        )
        == "invalid-schema"
    )


def test_all_stable_errors_are_package_generated(tmp_path: Path) -> None:
    source, domain, components, singles, pair = _bundle_records()
    bundle = _from_records()
    bad_schema = tmp_path / "bad-schema"
    _write_bundle(bad_schema, component_columns=COMPONENT_COLUMNS + ("extra",))

    def invalid_unit_record():
        return SingleParameterRecord(
            "invalid-unit-record",
            "methane",
            "molar_mass",
            16,
            **PROVENANCE,
        )

    def invalid_correlation():
        return ConstantPlusSumOfExponentialsCorrelation(
            "invalid-correlation-record",
            "methane",
            "segment_diameter",
            3 * unit_registry.angstrom,
            (),
            **PROVENANCE,
        )

    missing_source = SingleParameterRecord(
        "missing-source-record",
        "methane",
        "segment_count",
        1,
        source_id="missing-source",
        locator="Table 1, test cell",
        domain_id="domain-a",
    )
    missing_domain = SingleParameterRecord(
        "missing-domain-record",
        "methane",
        "segment_count",
        1,
        source_id="source-a",
        locator="Table 1, test cell",
        domain_id="missing-domain",
    )
    missing_model_source = ModelParameterRecord(
        "missing-model-source",
        "ionic_region_relative_permittivity",
        8,
        source_id="missing-source",
        locator="Figure 1",
        domain_id="domain-a",
    )
    missing_model_domain = ModelParameterRecord(
        "missing-model-domain",
        "ion_fraction_suppression_coefficient",
        Decimal("7.01"),
        source_id="source-a",
        locator="Eq. 11",
        domain_id="missing-domain",
    )
    assert _error_code(lambda: bundle.select((["not-a-string"],))) == "invalid-selection"
    assert (
        _error_code(
            lambda: _ParameterBundle.from_records(
                bundle_id="test-bundle",
                bundle_version=1,
                purpose="package-test-fixture",
                sources=(source,),
                domains=(domain,),
                components=components,
                models=(missing_model_source,),
            )
        )
        == "invalid-source"
    )
    assert (
        _error_code(
            lambda: _ParameterBundle.from_records(
                bundle_id="test-bundle",
                bundle_version=1,
                purpose="package-test-fixture",
                sources=(source,),
                domains=(domain,),
                components=components,
                models=(missing_model_domain,),
            )
        )
        == "invalid-domain"
    )

    codes = {
        _error_code(lambda: _ParameterBundle.from_path(bad_schema)),
        _error_code(
            lambda: _ParameterBundle.from_records(
                bundle_id="Bad-ID",
                bundle_version=1,
                purpose="package-test-fixture",
                sources=(source,),
                domains=(domain,),
                components=components,
            )
        ),
        _error_code(
            lambda: _ParameterBundle.from_records(
                bundle_id="test-bundle",
                bundle_version=1,
                purpose="model-activation",
                sources=(source,),
                domains=(domain,),
                components=components,
            )
        ),
        _error_code(invalid_unit_record),
        _error_code(invalid_correlation),
        _error_code(
            lambda: _ParameterBundle.from_records(
                bundle_id="test-bundle",
                bundle_version=1,
                purpose="package-test-fixture",
                sources=(source,),
                domains=(domain,),
                components=components,
                singles=(missing_source,),
            )
        ),
        _error_code(
            lambda: _ParameterBundle.from_records(
                bundle_id="test-bundle",
                bundle_version=1,
                purpose="package-test-fixture",
                sources=(source,),
                domains=(domain,),
                components=components,
                singles=(missing_domain,),
            )
        ),
        _error_code(
            lambda: _ParameterBundle.from_records(
                bundle_id="test-bundle",
                bundle_version=1,
                purpose="package-test-fixture",
                sources=(source,),
                domains=(domain,),
                components=components,
                singles=singles,
                pairs=(pair, pair),
            )
        ),
        _error_code(lambda: _ParameterBundle.from_catalog("missing-catalog", version=1)),
        _error_code(lambda: bundle.select(("CH4",))),
    }
    assert codes == {
        "invalid-schema",
        "invalid-id",
        "invalid-value",
        "invalid-unit",
        "invalid-source",
        "invalid-domain",
        "invalid-correlation",
        "duplicate-record",
        "missing-record",
        "invalid-selection",
    }


def test_catalog_discovery_returns_sorted_immutable_bundles() -> None:
    catalog = _ParameterBundle.catalog()

    assert isinstance(catalog, tuple)
    assert tuple((bundle.bundle_id, bundle.bundle_version) for bundle in catalog) == tuple(
        sorted((bundle.bundle_id, bundle.bundle_version) for bundle in catalog)
    )
    assert {bundle.bundle_id for bundle in catalog} == {
        "ascani-2022-case-study-2-electrolyte-lle",
        "figiel-2025-aqueous-iodides",
            "figiel-2025-reference-electrolytes",
            "gross-2001-methane-ethane",
            "gross-2001-propane",
            "held-2008-water-self-ionization",
            "khudaida-2026-figure-2-electrolyte-lle",
        }
    assert all(bundle.purpose == "installed-catalog" for bundle in catalog)
    assert all(bundle.components and bundle.domains and bundle.sources for bundle in catalog)


def test_catalog_loading_rejects_a_non_catalog_bundle(tmp_path: Path) -> None:
    bundle = _from_records(purpose="user-provided")
    path = tmp_path / "test-bundle"
    bundle.to_path(path)

    with pytest.raises(ParameterError) as exc_info:
        _ParameterBundle._from_catalog_path("test-bundle", 1, path)

    assert exc_info.value.code == "invalid-schema"
    assert dict(exc_info.value.context) == {
        "field": "catalog-purpose",
        "value": "user-provided",
    }


def test_source_backed_catalog_and_disconnected_contract_fixture() -> None:
    gross = _ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1)
    assert gross.purpose == "installed-catalog"
    assert gross.fingerprint == (
        "sha256:49641a26a9f81b614454b83a755df7763ef5f762e13eaaea1600558dc2d54f0b"
    )
    assert {source.source_id: source.doi for source in gross.sources} == {
        "gross-sadowski-2001": "10.1021/ie0003887"
    }
    assert tuple(component.component_id for component in gross.components) == (
        "ethane",
        "methane",
    )
    singles = {
        (record.component_id, record.family): record
        for record in gross.records
        if isinstance(record, SingleParameterRecord)
    }
    expected = {
        ("methane", "molar_mass"): (Decimal("0.016043"), "kilogram / mole"),
        ("methane", "segment_count"): (Decimal("1.0000"), "dimensionless"),
        ("methane", "segment_diameter"): (Decimal("3.7039"), "angstrom"),
        ("methane", "dispersion_energy_over_k"): (Decimal("150.03"), "kelvin"),
        ("ethane", "molar_mass"): (Decimal("0.03007"), "kilogram / mole"),
        ("ethane", "segment_count"): (Decimal("1.6069"), "dimensionless"),
        ("ethane", "segment_diameter"): (Decimal("3.5206"), "angstrom"),
        ("ethane", "dispersion_energy_over_k"): (Decimal("191.42"), "kelvin"),
    }
    assert set(singles) == set(expected)
    for key, (magnitude, unit) in expected.items():
        value = singles[key].value
        if unit == "dimensionless":
            assert value == magnitude
        else:
            assert value.magnitude == magnitude
            assert str(value.units) == unit
        assert singles[key].source_id == "gross-sadowski-2001"
        assert singles[key].locator.startswith("Table 2")

    methane_domain = next(
        domain for domain in gross.domains if domain.domain_id == "methane-fit-range"
    )
    ethane_domain = next(
        domain for domain in gross.domains if domain.domain_id == "ethane-fit-range"
    )
    assert (
        methane_domain.temperature_min.magnitude,
        methane_domain.temperature_max.magnitude,
    ) == (
        Decimal("97"),
        Decimal("300"),
    )
    assert (
        ethane_domain.temperature_min.magnitude,
        ethane_domain.temperature_max.magnitude,
    ) == (
        Decimal("90"),
        Decimal("305"),
    )
    pair = next(record for record in gross.records if isinstance(record, PairParameterRecord))
    assert (pair.component_id_a, pair.component_id_b, pair.family, pair.value) == (
        "ethane",
        "methane",
        "k_ij",
        Decimal("0"),
    )
    assert pair.locator == (
        'p. 1254, paragraph beginning "Results for a ternary system of CO2-methane-ethane", '
        'statement "one can safely assume a kij equal to 0"'
    )
    assert all(
        record.value != Decimal("0.0003") for record in gross.records if hasattr(record, "value")
    )
    selected = gross.select(("ethane", "methane"))
    assert selected.component_ids == ("ethane", "methane")
    assert selected.fingerprint == (
        "sha256:99ae0ffd9f61aee8c2e0fb718960534a2d0700c1adf5b72c1984517eded21f5b"
    )

    fixture_path = Path(__file__).parent / "data" / "bundles" / "figiel-2025-contract" / "1"
    fixture = _ParameterBundle.from_path(fixture_path)
    assert fixture.purpose == "package-test-fixture"
    assert fixture.fingerprint == (
        "sha256:3cc8b454266c98d2225727f17944f57dd3e325548d8b621b93ea664803e8d22f"
    )
    assert {source.source_id: source.doi for source in fixture.sources} == {
        "figiel-yu-held-2025": "10.1021/acs.iecr.5c00475",
    }
    fixture_records = {record.record_id: record for record in fixture.records}
    assert fixture_records["potassium-charge-number"].value == 1
    assert fixture_records["potassium-born-diameter"].value.magnitude == Decimal("4.150")
    assert fixture_records["potassium-solvation-factor"].value == Decimal("1")
    assert fixture_records["potassium-iodide-kij"].value == Decimal("0")
    assert fixture_records["water-relative-permittivity"].constant == Decimal("78.09")
    diameter = fixture_records["water-segment-diameter-correlation"]
    assert diameter.constant.magnitude == Decimal("2.7927")
    assert tuple(
        (term.amplitude.magnitude, term.exponent_coefficient.magnitude) for term in diameter.terms
    ) == (
        (Decimal("10.11"), Decimal("-0.01775")),
        (Decimal("-1.417"), Decimal("-0.01146")),
    )
    assert fixture_records["water-association-energy"].value.magnitude == Decimal("2425.7")
    assert fixture_records["water-association-volume"].value == Decimal("0.04509")
    assert {
        (record.site_id, record.site_class, record.multiplicity)
        for record in fixture.records
        if isinstance(record, SiteRecord)
    } == {("a", "a", 1), ("b", "b", 1)}
    assert {
        record_id: (
            fixture_records[record_id].source_id,
            fixture_records[record_id].locator,
        )
        for record_id in (
            "potassium-born-diameter",
            "potassium-iodide-kij",
            "potassium-solvation-factor",
            "water-association-energy",
            "water-association-volume",
            "water-relative-permittivity",
            "water-segment-diameter-correlation",
        )
    } == {
        "potassium-born-diameter": (
            "figiel-yu-held-2025",
            "Table 3, K+ row, d_i^Born column",
        ),
        "potassium-iodide-kij": (
            "figiel-yu-held-2025",
            "Table 4, I- row, K+ column",
        ),
        "potassium-solvation-factor": (
            "figiel-yu-held-2025",
            "p. 9411, Model Parameters paragraph, parenthetical statement that the solvation parameter is set to 1 for all ions",
        ),
        "water-association-energy": (
            "figiel-yu-held-2025",
            "Table 2, water column, epsilon^AiBi/k_B row",
        ),
        "water-association-volume": (
            "figiel-yu-held-2025",
            "Table 2, water column, kappa^AiBi row",
        ),
        "water-relative-permittivity": (
            "figiel-yu-held-2025",
            "Table 1, water row, dielectric constant column",
        ),
        "water-segment-diameter-correlation": (
            "figiel-yu-held-2025",
            "Table 2 footnote b, water segment-diameter correlation",
        ),
    }

    assert (
        _error_code(lambda: _ParameterBundle.from_catalog("figiel-2025-contract", version=1))
        == "missing-record"
    )
    catalog = Path(epcsaft.__file__).parent / "data" / "bundles"
    assert sorted(
        path.relative_to(catalog).as_posix() for path in catalog.rglob("*") if path.is_file()
    ) == [
        "ascani-2022-case-study-2-electrolyte-lle/1/association.csv",
        "ascani-2022-case-study-2-electrolyte-lle/1/bundle.toml",
            "ascani-2022-case-study-2-electrolyte-lle/1/components.csv",
            "ascani-2022-case-study-2-electrolyte-lle/1/correlations.toml",
            "ascani-2022-case-study-2-electrolyte-lle/1/model.csv",
            "ascani-2022-case-study-2-electrolyte-lle/1/pair.csv",
        "ascani-2022-case-study-2-electrolyte-lle/1/single.csv",
        "ascani-2022-case-study-2-electrolyte-lle/1/sites.csv",
        "figiel-2025-aqueous-iodides/1/association.csv",
        "figiel-2025-aqueous-iodides/1/bundle.toml",
        "figiel-2025-aqueous-iodides/1/components.csv",
        "figiel-2025-aqueous-iodides/1/correlations.toml",
        "figiel-2025-aqueous-iodides/1/model.csv",
        "figiel-2025-aqueous-iodides/1/pair.csv",
        "figiel-2025-aqueous-iodides/1/single.csv",
        "figiel-2025-aqueous-iodides/1/sites.csv",
        "figiel-2025-reference-electrolytes/1/association.csv",
        "figiel-2025-reference-electrolytes/1/bundle.toml",
        "figiel-2025-reference-electrolytes/1/components.csv",
        "figiel-2025-reference-electrolytes/1/correlations.toml",
        "figiel-2025-reference-electrolytes/1/model.csv",
        "figiel-2025-reference-electrolytes/1/pair.csv",
        "figiel-2025-reference-electrolytes/1/single.csv",
        "figiel-2025-reference-electrolytes/1/sites.csv",
            "gross-2001-methane-ethane/1/bundle.toml",
            "gross-2001-methane-ethane/1/components.csv",
            "gross-2001-methane-ethane/1/model.csv",
            "gross-2001-methane-ethane/1/pair.csv",
        "gross-2001-methane-ethane/1/single.csv",
            "gross-2001-propane/1/bundle.toml",
            "gross-2001-propane/1/components.csv",
            "gross-2001-propane/1/model.csv",
            "gross-2001-propane/1/single.csv",
        "held-2008-water-self-ionization/1/association.csv",
        "held-2008-water-self-ionization/1/bundle.toml",
        "held-2008-water-self-ionization/1/components.csv",
        "held-2008-water-self-ionization/1/correlations.toml",
        "held-2008-water-self-ionization/1/model.csv",
        "held-2008-water-self-ionization/1/pair.csv",
        "held-2008-water-self-ionization/1/single.csv",
        "held-2008-water-self-ionization/1/sites.csv",
        "khudaida-2026-figure-2-electrolyte-lle/1/association.csv",
        "khudaida-2026-figure-2-electrolyte-lle/1/bundle.toml",
        "khudaida-2026-figure-2-electrolyte-lle/1/components.csv",
        "khudaida-2026-figure-2-electrolyte-lle/1/correlations.toml",
        "khudaida-2026-figure-2-electrolyte-lle/1/model.csv",
        "khudaida-2026-figure-2-electrolyte-lle/1/pair.csv",
        "khudaida-2026-figure-2-electrolyte-lle/1/single.csv",
            "khudaida-2026-figure-2-electrolyte-lle/1/sites.csv",
            "provenance_names.toml",
        ]


def test_gross_2001_propane_catalog_is_exact() -> None:
    bundle = _ParameterBundle.from_catalog("gross-2001-propane", version=1)
    parameters = bundle.select(("propane",))
    values = {
        record.family: record.value
        for record in parameters.records
        if isinstance(record, SingleParameterRecord)
    }

    assert values["molar_mass"].to("kilogram / mole").magnitude == Decimal("0.044096")
    assert values["segment_count"] == Decimal("2.0020")
    assert values["segment_diameter"].to("angstrom").magnitude == Decimal("3.6184")
    assert values["dispersion_energy_over_k"].to("kelvin").magnitude == Decimal("208.11")
    assert epcsaft.Mixture(parameters).component_ids == ("propane",)


def test_figiel_reference_electrolyte_catalog_is_exact_and_sparse() -> None:
    bundle = _ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1)
    assert bundle.purpose == "installed-catalog"
    assert {component.component_id for component in bundle.components} == {
        "bromide-anion",
        "chloride-anion",
        "ethanol",
        "lithium-cation",
        "methanol",
        "potassium-cation",
        "sodium-cation",
        "water",
    }
    figiel_domain = next(
        domain for domain in bundle.domains if domain.domain_id == "figiel-298-k-1-bar"
    )
    assert figiel_domain.ion_mole_fraction_max == Decimal("0.38")

    values = {
        (record.component_id, record.family): record.value
        for record in bundle.records
        if isinstance(record, SingleParameterRecord)
    }
    assert values[("water", "molar_mass")].to("gram / mole").magnitude == Decimal("18.0153")
    assert values[("water", "segment_count")] == Decimal("1.2047")
    assert values[("water", "dispersion_energy_over_k")].magnitude == Decimal("353.95")
    assert values[("water", "charge_number")] == 0
    assert values[("water", "solvation_factor")] == Decimal("1.5")
    assert values[("methanol", "molar_mass")].to("gram / mole").magnitude == Decimal("32.042")
    assert values[("ethanol", "molar_mass")].to("gram / mole").magnitude == Decimal("46.069")
    assert values[("ethanol", "segment_count")] == Decimal("2.3827")
    assert values[("ethanol", "dispersion_energy_over_k")].magnitude == Decimal("198.24")
    assert values[("ethanol", "solvation_factor")] == Decimal("1.6")
    assert values[("lithium-cation", "born_diameter")].magnitude == Decimal("2.784")
    assert values[("sodium-cation", "segment_diameter")].magnitude == Decimal("2.8232")
    assert values[("sodium-cation", "born_diameter")].magnitude == Decimal("3.445")
    assert values[("potassium-cation", "born_diameter")].magnitude == Decimal("4.150")
    assert values[("chloride-anion", "segment_diameter")].magnitude == Decimal("2.7560")
    assert values[("chloride-anion", "born_diameter")].magnitude == Decimal("4.100")
    assert values[("bromide-anion", "born_diameter")].magnitude == Decimal("4.480")

    pairs = {
        (record.component_id_a, record.component_id_b): record.value
        for record in bundle.records
        if isinstance(record, PairParameterRecord) and record.family == "k_ij"
    }
    assert pairs[("methanol", "water")] == Decimal("-0.0878")
    assert pairs[("ethanol", "water")] == Decimal("-0.0617")
    assert pairs[("ethanol", "sodium-cation")] == Decimal("0.05")
    assert pairs[("chloride-anion", "ethanol")] == Decimal("0.8")
    assert pairs[("lithium-cation", "water")] == Decimal("-0.4")
    assert pairs[("bromide-anion", "potassium-cation")] == Decimal("-0.35")
    assert pairs[("chloride-anion", "sodium-cation")] == Decimal("0.8")
    assert all(
        record.family == "k_ij"
        for record in bundle.records
        if isinstance(record, PairParameterRecord)
    )

    associations = {
        (
            record.component_id_a,
            record.site_id_a,
            record.component_id_b,
            record.site_id_b,
            record.family,
        ): record.value
        for record in bundle.records
        if isinstance(record, AssociationParameterRecord)
    }
    assert associations[("ethanol", "a", "ethanol", "b", "association_volume")] == Decimal(
        "0.03238"
    )
    assert associations[("ethanol", "b", "water", "a", "association_energy_over_k")].to(
        "kelvin"
    ).magnitude == Decimal("2539.55")
    assert associations[("ethanol", "a", "water", "b", "association_volume")] == Decimal(
        "0.03797842747120087"
    )

    for cation, anion in (
        ("lithium-cation", "chloride-anion"),
        ("sodium-cation", "chloride-anion"),
        ("potassium-cation", "chloride-anion"),
        ("lithium-cation", "bromide-anion"),
        ("sodium-cation", "bromide-anion"),
        ("potassium-cation", "bromide-anion"),
    ):
        selected = ("water", cation, anion)
        assert epcsaft.Mixture(bundle.select(selected)).component_ids == selected

    methanol_nabr = ("methanol", "sodium-cation", "bromide-anion")
    assert epcsaft.Mixture(bundle.select(methanol_nabr)).component_ids == methanol_nabr

    ethanol_nacl = ("water", "ethanol", "sodium-cation", "chloride-anion")
    assert epcsaft.Mixture(bundle.select(ethanol_nacl)).component_ids == ethanol_nacl


def test_figiel_aqueous_iodide_catalog_is_exact_and_source_complete() -> None:
    bundle = _ParameterBundle.from_catalog("figiel-2025-aqueous-iodides", version=1)
    assert bundle.fingerprint == (
        "sha256:d1fec8a0c7a44796fb23c0bebbdf9733bfd1776b568e8155b3f4eb120efe3a03"
    )
    source = next(source for source in bundle.sources if source.source_id == "figiel-yu-held-2025")
    assert source.doi == "10.1021/acs.iecr.5c00475"
    assert "ce80533925a91bc59d8d0d8056113c40611ca26c2edf04aced76986d50bd4bae" in (source.use_basis)
    assert {component.component_id for component in bundle.components} == {
        "water",
        "sodium-cation",
        "potassium-cation",
        "iodide-anion",
    }
    values = {
        (record.component_id, record.family): record.value
        for record in bundle.records
        if isinstance(record, SingleParameterRecord)
    }
    assert values[("iodide-anion", "segment_count")] == Decimal("1")
    assert values[("iodide-anion", "segment_diameter")].magnitude == Decimal("3.6672")
    assert values[("iodide-anion", "dispersion_energy_over_k")].magnitude == Decimal("200.00")
    assert values[("iodide-anion", "charge_number")] == -1
    assert values[("iodide-anion", "born_diameter")].magnitude == Decimal("4.985")

    pairs = {
        (record.component_id_a, record.component_id_b): record.value
        for record in bundle.records
        if isinstance(record, PairParameterRecord) and record.family == "k_ij"
    }
    assert pairs[("iodide-anion", "water")] == Decimal("-0.05")
    assert pairs[("iodide-anion", "sodium-cation")] == Decimal("0.45")
    assert pairs[("iodide-anion", "potassium-cation")] == Decimal("0")

    fingerprints = {
        "sodium-cation": "sha256:3f33a43a6ff8eb127bb107ed34d22a2c6368b127e08fc1f64456bab592d78784",
        "potassium-cation": "sha256:5cc92296330a7167964ca76c027ec0f7f1a6dbb75cea92c774d893edbc7727ad",
    }
    for cation, fingerprint in fingerprints.items():
        selected = ("water", cation, "iodide-anion")
        parameters = bundle.select(selected)
        assert parameters.fingerprint == fingerprint
        assert epcsaft.Mixture(parameters).component_ids == selected
    assert (
        _error_code(lambda: bundle.select(("water", "lithium-cation", "iodide-anion")))
        == "invalid-selection"
    )


def test_khudaida_figure_2_catalog_preserves_the_published_model_tuple() -> None:
    bundle = _ParameterBundle.from_catalog("khudaida-2026-figure-2-electrolyte-lle", version=1)
    assert bundle.purpose == "installed-catalog"
    assert {component.component_id for component in bundle.components} == {
        "water",
        "ethanol",
        "isobutanol",
        "sodium-cation",
        "chloride-anion",
    }

    values = {
        (record.component_id, record.family): record.value
        for record in bundle.records
        if isinstance(record, SingleParameterRecord)
    }
    assert values[("isobutanol", "segment_count")] == Decimal("2.8440")
    assert values[("isobutanol", "segment_diameter")].magnitude == Decimal("3.5561")
    assert values[("isobutanol", "dispersion_energy_over_k")].magnitude == Decimal("257.13")
    assert values[("sodium-cation", "born_diameter")].magnitude == Decimal("3.445")
    assert values[("chloride-anion", "born_diameter")].magnitude == Decimal("4.100")

    pair_records = [record for record in bundle.records if isinstance(record, PairParameterRecord)]
    assert len(pair_records) == 10
    pairs = {
        frozenset((record.component_id_a, record.component_id_b)): record.value
        for record in pair_records
    }
    assert pairs == {
        frozenset(("water", "ethanol")): Decimal("-0.06167"),
        frozenset(("water", "isobutanol")): Decimal("-0.0071"),
        frozenset(("water", "sodium-cation")): Decimal("-0.3"),
        frozenset(("water", "chloride-anion")): Decimal("-0.3"),
        frozenset(("ethanol", "isobutanol")): Decimal("0.00154"),
        frozenset(("ethanol", "sodium-cation")): Decimal("0.05"),
        frozenset(("ethanol", "chloride-anion")): Decimal("0.8"),
        frozenset(("isobutanol", "sodium-cation")): Decimal("0"),
        frozenset(("isobutanol", "chloride-anion")): Decimal("0"),
        frozenset(("sodium-cation", "chloride-anion")): Decimal("0.317"),
    }

    association_records = [
        record for record in bundle.records if isinstance(record, AssociationParameterRecord)
    ]
    assert len(association_records) == 18
    associations = {
        (
            record.component_id_a,
            record.site_id_a,
            record.component_id_b,
            record.site_id_b,
            record.family,
        ): record.value
        for record in association_records
    }
    assert len(associations) == 18
    assert associations[
        ("isobutanol", "b", "water", "a", "association_energy_over_k")
    ].magnitude == Decimal("2470.0")
    assert associations[("isobutanol", "b", "water", "a", "association_volume")] == Decimal(
        "0.01299623825352431"
    )
    assert associations[
        ("ethanol", "a", "isobutanol", "b", "association_energy_over_k")
    ].magnitude == Decimal("2583.85")
    assert associations[("ethanol", "a", "isobutanol", "b", "association_volume")] == Decimal(
        "0.01119849133244554"
    )

    model_values = {
        record.family: record.value
        for record in bundle.records
        if isinstance(record, ModelParameterRecord)
    }
    assert model_values["ion_fraction_suppression_coefficient"] == Decimal("7.01")
    assert model_values["ionic_region_relative_permittivity"] == Decimal("8")

    water_sigma = next(
        record
        for record in bundle.records
        if isinstance(record, ConstantPlusSumOfExponentialsCorrelation)
        and record.component_id == "water"
        and record.family == "segment_diameter"
    )
    assert water_sigma.constant.magnitude == Decimal("2.7927")
    assert tuple(
        (term.amplitude.magnitude, term.exponent_coefficient.magnitude)
        for term in water_sigma.terms
    ) == (
        (Decimal("10.11"), Decimal("-0.01775")),
        (Decimal("-1.417"), Decimal("-0.01146")),
    )

    figiel = _ParameterBundle.from_catalog("figiel-2025-reference-electrolytes", version=1)
    assert figiel.fingerprint == (
        "sha256:d56aa57f04f48f1881865fb85daa700b66d3f1bcb75e8a9d7149baedfc925352"
    )
    assert figiel.select(("water", "sodium-cation", "chloride-anion")).fingerprint == (
        "sha256:d01bc3e64f19293902c00d9cc52060fccc3945ce736f0cd3e2f756be226cf5bf"
    )


def test_ascani_case_study_2_catalog_preserves_explicit_source_families() -> None:
    bundle = _ParameterBundle.from_catalog("ascani-2022-case-study-2-electrolyte-lle", version=1)
    assert bundle.fingerprint == (
        "sha256:756c71c8d7bd8c66e499eda037644f273c911333e689ca21ce986b7f0102fd25"
    )
    assert bundle.select(
        ("water", "1-butanol", "sodium-cation", "potassium-cation", "chloride-anion")
    ).fingerprint == ("sha256:ce50a54621b13538e649e7478aeafcb01e5d6f0fb0f6d53f2328a10063b825c1")
    assert tuple(component.component_id for component in bundle.components) == (
        "1-butanol",
        "chloride-anion",
        "potassium-cation",
        "sodium-cation",
        "water",
    )

    singles = {
        (record.component_id, record.family): record.value
        for record in bundle.records
        if isinstance(record, SingleParameterRecord)
    }
    assert singles[("1-butanol", "segment_count")] == Decimal("2.7510")
    assert singles[("sodium-cation", "segment_diameter")].magnitude == Decimal("2.8232")
    assert singles[("potassium-cation", "segment_diameter")].magnitude == Decimal("3.3417")
    assert singles[("chloride-anion", "segment_diameter")].magnitude == Decimal("2.7560")
    assert all(family not in {"born_diameter", "solvation_factor"} for _, family in singles)
    assert singles[("sodium-cation", "relative_permittivity")] == Decimal("8")
    assert singles[("potassium-cation", "relative_permittivity")] == Decimal("8")
    assert singles[("chloride-anion", "relative_permittivity")] == Decimal("8")

    pairs = {
        (record.family, frozenset((record.component_id_a, record.component_id_b))): record.value
        for record in bundle.records
        if isinstance(record, PairParameterRecord)
    }
    assert len([key for key in pairs if key[0] == "k_ij"]) == 9
    assert ("k_ij", frozenset(("sodium-cation", "potassium-cation"))) not in pairs
    assert pairs[("k_ij", frozenset(("water", "1-butanol")))] == Decimal("-0.102")
    assert pairs[("k_ij", frozenset(("sodium-cation", "chloride-anion")))] == Decimal("0.317")
    assert pairs[("k_ij", frozenset(("potassium-cation", "chloride-anion")))] == Decimal("0.064")
    assert len([key for key in pairs if key[0] == "l_ij"]) == 10
    assert pairs[("l_ij", frozenset(("water", "1-butanol")))] == Decimal("-0.0044")
    assert all(
        value == 0
        for key, value in pairs.items()
        if key[0] == "l_ij" and value != Decimal("-0.0044")
    )

    associations = {
        (
            record.component_id_a,
            record.site_id_a,
            record.component_id_b,
            record.site_id_b,
            record.family,
        ): record.value
        for record in bundle.records
        if isinstance(record, AssociationParameterRecord)
    }
    assert len(associations) == 8
    assert associations[
        ("1-butanol", "b", "water", "a", "association_energy_over_k")
    ].magnitude == Decimal("2420.51662")
    assert associations[("1-butanol", "b", "water", "a", "association_volume")] == Decimal(
        "0.016946922015303192"
    )


def test_held_2008_water_self_ionization_catalog_is_exact_and_source_complete() -> None:
    bundle = _ParameterBundle.from_catalog("held-2008-water-self-ionization", version=1)
    assert bundle.fingerprint == (
        "sha256:ebb6350a0820ff88f10f6179ea75df909773d8dc871020fbf08c3fae0ec067bb"
    )
    source = next(
        source
        for source in bundle.sources
        if source.source_id == "held-cameretti-sadowski-2008"
    )
    assert source.doi == "10.1016/j.fluid.2008.06.010"
    selected = ("water", "hydronium-cation", "hydroxide-anion")
    parameters = bundle.select(selected)
    model = epcsaft.Mixture(parameters)

    assert parameters.component_ids == selected
    assert tuple(component.formula for component in parameters.components) == (
        "H2O",
        "H3O+",
        "OH-",
    )
    assert model._model.topology == "ionic-fully-dissociated-debye-huckel"

    singles = {
        (record.component_id, record.family): record.value
        for record in parameters.records
        if isinstance(record, SingleParameterRecord)
    }
    assert singles[("water", "segment_count")] == Decimal("1.204659")
    assert singles[("water", "dispersion_energy_over_k")].magnitude == Decimal("353.9449")
    assert singles[("hydronium-cation", "segment_diameter")].magnitude == Decimal("2.2740")
    assert singles[("hydronium-cation", "dispersion_energy_over_k")].magnitude == Decimal(
        "1616.4939"
    )
    assert singles[("hydroxide-anion", "segment_diameter")].magnitude == Decimal("1.6401")
    assert singles[("hydroxide-anion", "dispersion_energy_over_k")].magnitude == Decimal(
        "2444.7555"
    )

    model_records = {
        record.family: record.value
        for record in parameters.records
        if isinstance(record, ModelParameterRecord)
    }
    assert model_records == {
        "electrolyte_formulation": "fully-dissociated-debye-huckel",
        "relative_permittivity_formulation": "solvent-only",
    }


@pytest.mark.parametrize(
    ("filename", "record_id", "missing_family"),
    (
        (
            "model.csv",
            "held-relative-permittivity-formulation",
            "relative_permittivity_formulation",
        ),
        ("single.csv", "hydronium-segment-diameter", "segment_diameter"),
        ("association.csv", "water-association-energy", "association_energy_over_k"),
        ("pair.csv", "hydronium-hydroxide-kij", "k_ij"),
    ),
)
def test_held_2008_catalog_fails_closed_on_missing_source_family(
    tmp_path: Path,
    filename: str,
    record_id: str,
    missing_family: str,
) -> None:
    catalog = Path(epcsaft.__file__).parent / "data" / "bundles"
    source = catalog / "held-2008-water-self-ionization" / "1"
    mutated = tmp_path / "held-mutated"
    shutil.copytree(source, mutated)
    path = mutated / filename
    if filename == "model.csv":
        path.unlink()
        manifest = mutated / "bundle.toml"
        manifest.write_text(
            manifest.read_text(encoding="utf-8").replace(
                '"association.csv", "model.csv", "correlations.toml"',
                '"association.csv", "correlations.toml"',
            ),
            encoding="utf-8",
        )
    else:
        rows = list(csv.reader(path.read_text(encoding="utf-8").splitlines()))
        with path.open("w", encoding="utf-8", newline="") as stream:
            csv.writer(stream).writerows(row for row in rows if not row or row[0] != record_id)

    bundle = _ParameterBundle.from_path(mutated)
    with pytest.raises(ParameterError) as exc_info:
        epcsaft.Mixture(bundle.select(("water", "hydronium-cation", "hydroxide-anion")))

    assert exc_info.value.code == "missing-record"
    assert ("family", missing_family) in exc_info.value.context


@pytest.mark.parametrize(
    ("record_id", "mutated_value"),
    (
        ("water-hydronium-kij", "0.1"),
        ("hydronium-hydroxide-kij", "0.9"),
    ),
)
def test_held_2008_catalog_rejects_dispersion_rule_mutation(
    tmp_path: Path,
    record_id: str,
    mutated_value: str,
) -> None:
    catalog = Path(epcsaft.__file__).parent / "data" / "bundles"
    source = catalog / "held-2008-water-self-ionization" / "1"
    mutated = tmp_path / "held-mutated"
    shutil.copytree(source, mutated)
    path = mutated / "pair.csv"
    rows = list(csv.DictReader(path.read_text(encoding="utf-8").splitlines()))
    for row in rows:
        if row["record_id"] == record_id:
            row["value"] = mutated_value
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)

    bundle = _ParameterBundle.from_path(mutated)
    with pytest.raises(ParameterError) as exc_info:
        epcsaft.Mixture(bundle.select(("water", "hydronium-cation", "hydroxide-anion")))

    assert exc_info.value.code == "invalid-value"
    assert ("family", "k_ij") in exc_info.value.context


def test_electrolyte_formulation_identity_is_closed() -> None:
    with pytest.raises(ParameterError, match="invalid-value"):
        ModelParameterRecord(
            "unknown-electrolyte-formulation",
            "electrolyte_formulation",
            "unknown-formulation",
            **PROVENANCE,
        )


def test_validity_domain_rejects_reversed_pressure_bounds() -> None:
    with pytest.raises(ParameterError, match="reversed-bounds"):
        ValidityDomain(
            "reversed-pressure",
            "reported-conditions",
            pressure_min=2 * unit_registry.bar,
            pressure_max=1 * unit_registry.bar,
        )

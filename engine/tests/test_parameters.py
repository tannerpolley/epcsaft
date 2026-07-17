import csv
from dataclasses import FrozenInstanceError
from decimal import Decimal
import json
from pathlib import Path

import pytest
from pint import UnitRegistry

import epcsaft
from epcsaft import ParameterBundle, ParameterError, ParameterSet, unit_registry
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
):
    source, domain, components, singles, pair = _bundle_records(
        registry, reverse=reverse
    )
    models = [
        ModelParameterRecord(
            "figiel-ionic-region-permittivity",
            "ionic_region_relative_permittivity",
            ionic_permittivity,
            **PROVENANCE,
        ),
        ModelParameterRecord(
            "figiel-dielectric-suppression",
            "dielectric_ion_suppression_coefficient",
            Decimal("7.01"),
            **PROVENANCE,
        ),
    ]
    if reverse:
        models.reverse()
    return ParameterBundle.from_records(
        bundle_id="test-bundle",
        bundle_version=1,
        purpose="package-test-fixture",
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
    return ComponentRecord("water", name="water", formula="H2O"), sites, association, correlations


def _parsed_records_bundle():
    source, domain, components, singles, pair = _bundle_records()
    water, sites, association, correlations = _advanced_records()
    return ParameterBundle.from_records(
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
                "dielectric_ion_suppression_coefficient",
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
                "dielectric_ion_suppression_coefficient",
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
    assert epcsaft.__all__ == (
        "ParameterBundle",
        "ParameterSet",
        "ParameterError",
        "unit_registry",
        "EPCSAFT",
        "EosResult",
    )
    foreign = UnitRegistry()
    bundle = _from_records()
    equivalent = _from_records(foreign, reverse=True)
    changed_model = _from_records(ionic_permittivity=Decimal("9"))

    assert isinstance(bundle, ParameterBundle)
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
    assert isinstance(selected, ParameterSet)
    assert selected.component_ids == ("ethane", "methane")
    assert tuple(component.component_id for component in selected.components) == selected.component_ids
    assert selected.bundle_fingerprint == bundle.fingerprint
    assert selected.fingerprint != reverse_selected.fingerprint
    with pytest.raises(FrozenInstanceError):
        selected.component_ids = ("methane",)

    golden = ParameterBundle.from_records(
        bundle_id="golden-bundle",
        bundle_version=1,
        purpose="package-test-fixture",
        sources=(SourceRecord("source-a", "Cafe\u0301\r\nLine", "test"),),
        domains=(ValidityDomain("unknown-domain", "unknown"),),
        components=(
            ComponentRecord("methane", name="Cafe\u0301\r\nLine", aliases=("CH4",)),
        ),
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
    parsed = ParameterBundle.from_path(bundle_path)
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
    assert _error_code(lambda: ParameterBundle.from_path(bad_columns)) == "invalid-schema"

    bad_model_columns = tmp_path / "bad-model-columns"
    _write_bundle(bad_model_columns, model_columns=MODEL_COLUMNS + ("extra",))
    assert _error_code(lambda: ParameterBundle.from_path(bad_model_columns)) == "invalid-schema"

    bad_inventory = tmp_path / "bad-inventory"
    _write_bundle(bad_inventory, extra_file=True)
    assert _error_code(lambda: ParameterBundle.from_path(bad_inventory)) == "invalid-schema"

    duplicate_pair = tmp_path / "duplicate-pair"
    _write_bundle(duplicate_pair, duplicate_pair=True)
    assert _error_code(lambda: ParameterBundle.from_path(duplicate_pair)) == "duplicate-record"

    bad_correlation = tmp_path / "bad-correlation"
    _write_bundle(bad_correlation, bad_correlation=True)
    assert _error_code(lambda: ParameterBundle.from_path(bad_correlation)) == "invalid-correlation"

    catalog_mismatch = tmp_path / "catalog-mismatch"
    _write_bundle(catalog_mismatch)
    assert (
        _error_code(
            lambda: ParameterBundle._from_catalog_path(
                "requested-bundle", 1, catalog_mismatch
            )
        )
        == "invalid-schema"
    )


def test_all_stable_errors_are_package_generated(tmp_path: Path) -> None:
    source, domain, components, singles, pair = _bundle_records()
    bundle = _from_records()
    bad_schema = tmp_path / "bad-schema"
    _write_bundle(bad_schema, component_columns=COMPONENT_COLUMNS + ("extra",))

    invalid_unit_record = lambda: SingleParameterRecord(
        "invalid-unit-record",
        "methane",
        "molar_mass",
        16,
        **PROVENANCE,
    )
    invalid_correlation = lambda: ConstantPlusSumOfExponentialsCorrelation(
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
        "dielectric_ion_suppression_coefficient",
        Decimal("7.01"),
        source_id="source-a",
        locator="Eq. 11",
        domain_id="missing-domain",
    )
    assert _error_code(lambda: bundle.select((["not-a-string"],))) == "invalid-selection"
    assert (
        _error_code(
            lambda: ParameterBundle.from_records(
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
            lambda: ParameterBundle.from_records(
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
        _error_code(lambda: ParameterBundle.from_path(bad_schema)),
        _error_code(
            lambda: ParameterBundle.from_records(
                bundle_id="Bad-ID",
                bundle_version=1,
                purpose="package-test-fixture",
                sources=(source,),
                domains=(domain,),
                components=components,
            )
        ),
        _error_code(
            lambda: ParameterBundle.from_records(
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
            lambda: ParameterBundle.from_records(
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
            lambda: ParameterBundle.from_records(
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
            lambda: ParameterBundle.from_records(
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
        _error_code(
            lambda: ParameterBundle.from_catalog("missing-catalog", version=1)
        ),
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


def test_source_backed_catalog_and_disconnected_contract_fixture() -> None:
    gross = ParameterBundle.from_catalog("gross-2001-methane-ethane", version=1)
    assert gross.purpose == "installed-catalog"
    assert gross.fingerprint == (
        "sha256:3755d6d907178e8d8ce3fa6a25c7b7c0d670c9643bcb56ea16ad6555acb98c03"
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

    methane_domain = next(domain for domain in gross.domains if domain.domain_id == "methane-fit-range")
    ethane_domain = next(domain for domain in gross.domains if domain.domain_id == "ethane-fit-range")
    assert (methane_domain.temperature_min.magnitude, methane_domain.temperature_max.magnitude) == (
        Decimal("97"),
        Decimal("300"),
    )
    assert (ethane_domain.temperature_min.magnitude, ethane_domain.temperature_max.magnitude) == (
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
    assert all(record.value != Decimal("0.0003") for record in gross.records if hasattr(record, "value"))
    selected = gross.select(("ethane", "methane"))
    assert selected.component_ids == ("ethane", "methane")
    assert selected.fingerprint == (
        "sha256:cb0c63530f96d51d5a24f48a2d02ae6567d163abdc373099ce096f7ff0a896f5"
    )

    fixture_path = Path(__file__).parent / "data" / "bundles" / "figiel-2025-contract" / "1"
    fixture = ParameterBundle.from_path(fixture_path)
    assert fixture.purpose == "package-test-fixture"
    assert fixture.fingerprint == (
        "sha256:202ad0ca74959db80e834ad8c22a11fdb1526df622e14cbacee4cdb0a7b9d02c"
    )
    assert {source.source_id: source.doi for source in fixture.sources} == {
        "ascani-held-2021": "10.1002/zaac.202100032",
        "figiel-yu-held-2025": "10.1021/acs.iecr.5c00475",
    }
    fixture_records = {record.record_id: record for record in fixture.records}
    assert fixture_records["potassium-charge-number"].value == 1
    assert fixture_records["potassium-born-diameter"].value.magnitude == Decimal("4.150")
    assert fixture_records["potassium-solvation-factor"].value == Decimal("1")
    assert fixture_records["potassium-iodide-kij"].value == Decimal("0")
    assert not any(
        isinstance(record, PairParameterRecord)
        and {record.component_id_a, record.component_id_b}
        == {"hydrogen-cation", "iodide-anion"}
        and record.family == "k_ij"
        for record in fixture.records
    )
    assert fixture_records["water-one-butanol-lij"].value == Decimal("-0.0044")
    assert fixture_records["water-one-butanol-khb-ij"].value == Decimal("0.026")
    assert fixture_records["water-relative-permittivity"].constant == Decimal("78.09")
    diameter = fixture_records["water-segment-diameter-correlation"]
    assert diameter.constant.magnitude == Decimal("2.7927")
    assert tuple(
        (term.amplitude.magnitude, term.exponent_coefficient.magnitude)
        for term in diameter.terms
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
            "water-one-butanol-khb-ij",
            "water-one-butanol-lij",
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
        "water-one-butanol-khb-ij": (
            "ascani-held-2021",
            "Table 4, k_hb_ij row, 1-Butanol column",
        ),
        "water-one-butanol-lij": (
            "ascani-held-2021",
            "Table 4, l_ij row, 1-Butanol column",
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

    assert _error_code(
        lambda: ParameterBundle.from_catalog("figiel-2025-contract", version=1)
    ) == "missing-record"
    catalog = Path(epcsaft.__file__).parent / "data" / "bundles"
    assert sorted(
        path.relative_to(catalog).as_posix()
        for path in catalog.rglob("*")
        if path.is_file()
    ) == [
        "gross-2001-methane-ethane/1/bundle.toml",
        "gross-2001-methane-ethane/1/components.csv",
        "gross-2001-methane-ethane/1/pair.csv",
        "gross-2001-methane-ethane/1/single.csv",
    ]

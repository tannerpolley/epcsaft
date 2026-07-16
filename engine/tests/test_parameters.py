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


def _from_records(registry: UnitRegistry = unit_registry, *, reverse: bool = False):
    source, domain, components, singles, pair = _bundle_records(
        registry, reverse=reverse
    )
    return ParameterBundle.from_records(
        bundle_id="test-bundle",
        bundle_version=1,
        purpose="package-test-fixture",
        sources=(source,),
        domains=(domain,),
        components=components,
        singles=singles,
        pairs=(pair,),
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
    )


def _write_bundle(
    path: Path,
    *,
    component_columns: tuple[str, ...] = COMPONENT_COLUMNS,
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
files = ["components.csv", "single.csv", "pair.csv", "sites.csv", "association.csv", "correlations.toml"]

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
    )
    foreign = UnitRegistry()
    bundle = _from_records()
    equivalent = _from_records(foreign, reverse=True)

    assert isinstance(bundle, ParameterBundle)
    assert tuple(component.component_id for component in bundle.components) == (
        "ethane",
        "methane",
    )
    assert tuple(record.record_id for record in bundle.records) == (
        "ethane-segment-count",
        "methane-ethane-kij",
        "methane-molar-mass",
    )
    assert all(
        isinstance(value, tuple)
        for value in (bundle.sources, bundle.domains, bundle.components, bundle.records)
    )
    assert bundle.fingerprint == equivalent.fingerprint
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
    }

    bad_columns = tmp_path / "bad-columns"
    _write_bundle(bad_columns, component_columns=COMPONENT_COLUMNS + ("extra",))
    assert _error_code(lambda: ParameterBundle.from_path(bad_columns)) == "invalid-schema"

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
    assert _error_code(lambda: bundle.select((["not-a-string"],))) == "invalid-selection"

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

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

GROSS_SADOWSKI_2001_TABLE_2 = """\
methane|methane|16.043|1.0000|3.7039|150.03|97|300
ethane|ethane|30.07|1.6069|3.5206|191.42|90|305
propane|propane|44.096|2.0020|3.6184|208.11|85|523
butane|butane|58.123|2.3316|3.7086|222.88|135|573
pentane|pentane|72.146|2.6896|3.7729|231.20|143|469
hexane|hexane|86.177|3.0576|3.7983|236.77|177|503
heptane|heptane|100.203|3.4831|3.8049|238.40|182|623
octane|octane|114.231|3.8176|3.8373|242.78|216|569
nonane|nonane|128.25|4.2079|3.8448|244.51|219|595
decane|decane|142.285|4.6627|3.8384|243.87|243|617
undecane|undecane|156.312|4.9082|3.8893|248.82|247|639
dodecane|dodecane|170.338|5.3060|3.8959|249.21|263|658
tridecane|tridecane|184.365|5.6877|3.9143|249.78|267|675
tetradecane|tetradecane|198.392|5.9002|3.9396|254.21|279|693
pentadecane|pentadecane|212.419|6.2855|3.9531|254.14|283|708
hexadecane|hexadecane|226.446|6.6485|3.9552|254.70|291|723
heptadecane|heptadecane|240.473|6.9809|3.9675|255.65|295|736
octadecane|octadecane|254.5|7.3271|3.9668|256.20|301|747
nonadecane|nonadecane|268.527|7.7175|3.9721|256.00|305|758
eicosane|eicosane|282.553|7.9849|3.9869|257.75|309|775
isobutane|isobutane|58.123|2.2616|3.7574|216.53|113|407
isopentane|isopentane|72.15|2.5620|3.8296|230.75|113|460
neopentane|neopentane|72.15|2.3543|3.9550|225.69|256|433
2-methylpentane|2-methylpentane|86.177|2.9317|3.8535|235.58|119|498
2-2-dimethylbutane|2,2-dimethylbutane|86.177|2.6008|4.0042|243.51|174|488
2-3-dimethylbutane|2,3-dimethylbutane|86.177|2.6853|3.9545|246.07|145|500
3-methylpentane|3-methylpentane|86.177|2.8852|3.8605|240.48|110|504
2-methylhexane|2-methylhexane|100.204|3.3478|3.8612|237.42|154|530
cyclopentane|cyclopentane|70.13|2.3655|3.7114|265.83|193|503
cyclohexane|cyclohexane|84.147|2.5303|3.8499|278.11|279|553
methylcyclopentane|methylcyclopentane|84.156|2.6130|3.8253|265.12|183|532
methylcyclohexane|methylcyclohexane|98.182|2.6637|3.9993|282.33|203|572
ethylcyclopentane|ethylcyclopentane|98.182|2.9062|3.8873|270.50|134|569
ethylcyclohexane|ethylcyclohexane|112.215|2.8256|4.1039|294.04|263|609
ethylene|ethylene|28.05|1.5930|3.4450|176.47|104|400
propylene|propylene|42.081|1.9597|3.5356|207.19|87|364
1-butene|1-butene|56.107|2.2864|3.6431|222.00|87|419
1-pentene|1-pentene|70.134|2.6006|3.7399|231.99|108|464
1-hexene|1-hexene|84.616|2.9853|3.7753|236.81|133|504
1-octene|1-octene|112.215|3.7424|3.8133|243.02|171|567
cyclopentene|cyclopentene|68.114|2.2934|3.6668|267.76|223|393
benzene|benzene|78.114|2.4653|3.6478|287.35|278|562
toluene|toluene|92.141|2.8149|3.7169|285.69|178|594
ethylbenzene|ethylbenzene|106.167|3.0799|3.7974|287.35|178|617
m-xylene|m-xylene|106.167|3.1861|3.7563|283.98|225|619
o-xylene|o-xylene|106.167|3.1362|3.7600|291.05|248|630
p-xylene|p-xylene|106.167|3.1723|3.7781|283.77|286|616
n-propyl-benzene|n-propyl benzene|120.194|3.3438|3.8438|288.13|173|638
tetralin|tetralin|132.205|3.3131|3.8750|325.07|237|720
n-butylbenzene|n-butylbenzene|134.221|3.7662|3.8727|283.07|185|660
biphenyl|biphenyl|154.211|3.8877|3.8151|327.42|342|773
carbon-monoxide|carbon monoxide|28.01|1.3097|3.2507|92.15|80|133
nitrogen|nitrogen|28.01|1.2053|3.3130|90.96|63|126
argon|argon|39.948|0.9285|3.4784|122.23|84|151
carbon-dioxide|carbon dioxide|44.01|2.0729|2.7852|169.21|216|304
sulfur-dioxide|sulfur dioxide|64.065|2.8611|2.6826|205.35|198|431
chlorine|chlorine|70.905|1.5514|3.3672|265.67|172|417
carbon-disulfide|carbon disulfide|76.143|1.6919|3.6172|334.82|200|552
methyl-chloride|methyl chloride|50.488|1.9297|3.2293|240.56|175|416
chloroethane|chloroethane|64.514|2.2638|3.4160|245.43|134|460
2-chloropropane|2-chloropropane|78.541|2.4151|3.6184|251.47|156|489
1-chlorobutane|1-chlorobutane|92.568|2.8585|3.6424|258.66|150|537
chlorobenzene|chlorobenzene|112.558|2.6485|3.7533|315.04|228|632
bromobenzene|bromobenzene|157.01|2.6456|3.8360|334.37|242|670
dimethyl-ether|dimethyl ether|46.069|2.3071|3.2528|211.06|200|400
diethyl-ether|diethyl ether|74.123|2.9686|3.5147|220.09|240|467
methyl-n-propyl-ether|methyl-n-propyl ether|74.123|3.0087|3.4569|222.73|220|476
methyl-methanoate|methyl methanoate|60.053|2.6784|3.0875|242.63|174|487
ethyl-methanoate|ethyl methanoate|74.079|2.8876|3.3109|246.47|193|508
n-propyl-methanoate|n-propyl methanoate|88.106|3.2088|3.4168|246.46|270|538
ethyl-ethanoate|ethyl ethanoate|88.106|3.5375|3.3079|230.80|189|523
n-propyl-ethanoate|n-propyl ethanoate|102.133|3.7861|3.4227|235.76|290|549
n-butyl-ethanoate|n-butyl ethanoate|116.16|3.9808|3.5427|242.52|199|579
isopropyl-ethanoate|isopropyl ethanoate|102.133|3.6084|3.4818|231.91|199|532
methyl-propanoate|methyl propanoate|88.106|3.4793|3.3142|234.96|185|530
ethyl-propanoate|ethyl propanoate|102.133|3.8371|3.4031|232.78|199|546
n-propyl-propanoate|n-propyl propanoate|116.16|4.1155|3.4875|235.60|220|568
methyl-butanoate|methyl butanoate|102.133|3.6758|3.4437|240.62|187|554
"""


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
        "crespo-2021-deea",
        "figiel-2025-aqueous-iodides",
        "figiel-2025-reference-electrolytes",
        "figiel-yu-held-2025",
        "gross-2001-methane-ethane",
        "gross-2001-propane",
        "gross-sadowski-2001",
        "held-2008-water-self-ionization",
    }
    assert all(bundle.purpose == "installed-catalog" for bundle in catalog)
    assert all(bundle.components and bundle.domains and bundle.sources for bundle in catalog)


def test_gross_sadowski_2001_catalog_is_source_complete() -> None:
    bundle = _ParameterBundle.from_catalog("gross-sadowski-2001", version=1)
    expected_table = {
        component_id: (
            name,
            Decimal(molar_mass),
            Decimal(segment_count),
            Decimal(segment_diameter),
            Decimal(dispersion_energy),
            Decimal(temperature_min),
            Decimal(temperature_max),
        )
        for (
            component_id,
            name,
            molar_mass,
            segment_count,
            segment_diameter,
            dispersion_energy,
            temperature_min,
            temperature_max,
        ) in (line.split("|") for line in GROSS_SADOWSKI_2001_TABLE_2.splitlines())
    }
    assert len(expected_table) == 78
    components = {component.component_id: component for component in bundle.components}
    assert set(components) == set(expected_table)
    assert len(bundle.components) == 78
    singles = {
        (record.component_id, record.family): record
        for record in bundle.records
        if isinstance(record, SingleParameterRecord)
    }
    assert set(singles) == {
        (component_id, family)
        for component_id in expected_table
        for family in {
            "molar_mass",
            "segment_count",
            "segment_diameter",
            "dispersion_energy_over_k",
        }
    }
    domains = {domain.domain_id: domain for domain in bundle.domains}
    assert {domain_id for domain_id in domains if domain_id.endswith("-fit-range")} == {
        f"{component_id}-fit-range" for component_id in expected_table
    }
    for component_id, (
        name,
        molar_mass,
        segment_count,
        segment_diameter,
        dispersion_energy,
        temperature_min,
        temperature_max,
    ) in expected_table.items():
        assert components[component_id].name == name
        expected_values = {
            "molar_mass": (molar_mass / Decimal("1000"), "kilogram / mole"),
            "segment_count": (segment_count, "dimensionless"),
            "segment_diameter": (segment_diameter, "angstrom"),
            "dispersion_energy_over_k": (dispersion_energy, "kelvin"),
        }
        for family, (magnitude, unit) in expected_values.items():
            record = singles[(component_id, family)]
            if unit == "dimensionless":
                assert record.value == magnitude
            else:
                assert record.value.to(unit).magnitude == magnitude
            assert record.source_id == "gross-sadowski-2001"
            assert record.locator.startswith(f"Table 2, {name} row")
            assert record.domain_id == f"{component_id}-fit-range"
        domain = domains[f"{component_id}-fit-range"]
        assert domain.temperature_min.to("kelvin").magnitude == temperature_min
        assert domain.temperature_max.to("kelvin").magnitude == temperature_max

    pairs = {
        frozenset((record.component_id_a, record.component_id_b)): record.value
        for record in bundle.records
        if isinstance(record, PairParameterRecord)
    }
    assert pairs == {
        frozenset(("methane", "butane")): Decimal("0.022"),
        frozenset(("methane", "pentane")): Decimal("0.024"),
        frozenset(("methane", "hexane")): Decimal("0.021"),
        frozenset(("methane", "heptane")): Decimal("0.016"),
        frozenset(("methane", "decane")): Decimal("0.056"),
        frozenset(("methane", "isobutane")): Decimal("0.028"),
        frozenset(("methane", "cyclohexane")): Decimal("0.045"),
        frozenset(("methane", "benzene")): Decimal("0.037"),
        frozenset(("methane", "toluene")): Decimal("0.052"),
        frozenset(("methane", "m-xylene")): Decimal("0.045"),
        frozenset(("methane", "tetralin")): Decimal("0.075"),
        frozenset(("nitrogen", "ethylene")): Decimal("0.075"),
        frozenset(("nitrogen", "hexane")): Decimal("0.119"),
        frozenset(("ethylene", "heptane")): Decimal("0.018"),
        frozenset(("ethylene", "propylene")): Decimal("-0.001"),
        frozenset(("propane", "butane")): Decimal("0.003"),
        frozenset(("propylene", "isobutane")): Decimal("-0.006"),
        frozenset(("pentane", "heptane")): Decimal("0.011"),
        frozenset(("carbon-dioxide", "methane")): Decimal("0.065"),
        frozenset(("carbon-dioxide", "propane")): Decimal("0.109"),
        frozenset(("carbon-dioxide", "butane")): Decimal("0.120"),
        frozenset(("carbon-dioxide", "pentane")): Decimal("0.143"),
        frozenset(("carbon-dioxide", "heptane")): Decimal("0.129"),
        frozenset(("carbon-dioxide", "methylcyclohexane")): Decimal("0.144"),
        frozenset(("methane", "ethane")): Decimal("0"),
        frozenset(("carbon-dioxide", "ethane")): Decimal("0.095"),
        frozenset(("carbon-dioxide", "decane")): Decimal("0.128"),
        frozenset(("carbon-dioxide", "cyclohexane")): Decimal("0.13"),
        frozenset(("1-chlorobutane", "hexane")): Decimal("0.017"),
    }
    assert {source.source_id: source.doi for source in bundle.sources} == {
        "gross-sadowski-2001": "10.1021/ie0003887"
    }
    assert bundle.fingerprint == (
        "sha256:c038241875a0d4c55218589b3f73453eeef927bb831bda703a9edebed7e8c19e"
    )

    pure = Parameters.from_catalog(
        "gross-sadowski-2001",
        components=("eicosane",),
        version=1,
    )
    assert pure.fingerprint == (
        "sha256:3d65b05330ea05fe2180a8efc7da51dc0254970325404b9c64d605135846f6cd"
    )
    pure_state = epcsaft.Mixture(pure).state(
        T=400 * unit_registry.kelvin,
        rho=100 * unit_registry.mole / unit_registry.meter**3,
        x=(1.0,),
    )
    assert pure_state.ares().check("[energy] / [substance]")

    forward = Parameters.from_catalog(
        "gross-sadowski-2001",
        components=("carbon-dioxide", "cyclohexane"),
        version=1,
    )
    reverse = Parameters.from_catalog(
        "gross-sadowski-2001",
        components=("cyclohexane", "carbon-dioxide"),
        version=1,
    )
    assert forward.fingerprint == (
        "sha256:cc7445a5965e17091df613bf980733f6718cfa247d4ba95c6c97c8a3f1935217"
    )
    assert reverse.fingerprint == (
        "sha256:2de8e86ea87c2da0e261c0fe0d51596eb98830c87d7b2801bb2de579b7e59f4c"
    )
    forward_state = epcsaft.Mixture(forward).state(
        T=298.15 * unit_registry.kelvin,
        rho=100 * unit_registry.mole / unit_registry.meter**3,
        x=(0.4, 0.6),
    )
    reverse_state = epcsaft.Mixture(reverse).state(
        T=298.15 * unit_registry.kelvin,
        rho=100 * unit_registry.mole / unit_registry.meter**3,
        x=(0.6, 0.4),
    )
    assert forward_state.ares().to("joule / mole").magnitude == pytest.approx(
        reverse_state.ares().to("joule / mole").magnitude,
        rel=1e-13,
        abs=1e-12,
    )


def test_figiel_yu_held_2025_catalog_is_source_complete() -> None:
    bundle = _ParameterBundle.from_catalog("figiel-yu-held-2025", version=1)
    expected_components = {
        "water",
        "methanol",
        "ethanol",
        "hydrogen-cation",
        "lithium-cation",
        "sodium-cation",
        "potassium-cation",
        "vanadyl-cation",
        "vanadium-iii-cation",
        "chloride-anion",
        "bromide-anion",
        "iodide-anion",
        "sulfate-anion",
    }
    components = {component.component_id: component for component in bundle.components}
    assert set(components) == expected_components
    assert len(bundle.components) == 13
    observed_formulas = {
        component_id: components[component_id].formula
        for component_id in expected_components
    }
    assert observed_formulas == {
        "water": "H2O",
        "methanol": "CH4O",
        "ethanol": "C2H6O",
        "hydrogen-cation": "H+",
        "lithium-cation": "Li+",
        "sodium-cation": "Na+",
        "potassium-cation": "K+",
        "vanadyl-cation": "VO2+",
        "vanadium-iii-cation": "V3+",
        "chloride-anion": "Cl-",
        "bromide-anion": "Br-",
        "iodide-anion": "I-",
        "sulfate-anion": "SO4^2-",
    }

    singles = {
        (record.component_id, record.family): record
        for record in bundle.records
        if isinstance(record, SingleParameterRecord)
    }
    expected_ions = {
        "hydrogen-cation": ("3.4654", "500.00", "1.218", 1),
        "lithium-cation": ("2.8449", "360.00", "2.784", 1),
        "sodium-cation": ("2.8232", "230.00", "3.445", 1),
        "potassium-cation": ("3.3417", "200.00", "4.150", 1),
        "vanadyl-cation": ("3.4638", "1488.28", "2.930", 2),
        "vanadium-iii-cation": ("3.1401", "1415.66", "2.940", 3),
        "chloride-anion": ("2.7560", "170.00", "4.100", -1),
        "bromide-anion": ("3.0707", "190.00", "4.480", -1),
        "iodide-anion": ("3.6672", "200.00", "4.985", -1),
        "sulfate-anion": ("2.6491", "80.00", "4.953", -2),
    }
    for component_id, (diameter, energy, born, charge) in expected_ions.items():
        assert singles[(component_id, "segment_count")].value == Decimal("1")
        assert singles[(component_id, "segment_diameter")].value.to(
            "angstrom"
        ).magnitude == Decimal(diameter)
        assert singles[(component_id, "dispersion_energy_over_k")].value.to(
            "kelvin"
        ).magnitude == Decimal(energy)
        assert singles[(component_id, "born_diameter")].value.to(
            "angstrom"
        ).magnitude == Decimal(born)
        assert singles[(component_id, "charge_number")].value == charge
        assert singles[(component_id, "solvation_factor")].value == Decimal("1")
        assert all(
            singles[(component_id, family)].source_id == "figiel-yu-held-2025"
            for family in {
                "segment_count",
                "segment_diameter",
                "dispersion_energy_over_k",
                "born_diameter",
                "charge_number",
                "solvation_factor",
            }
        )

    expected_solvents = {
        "water": ("18.0153", "1.2047", "353.95", "1.5"),
        "methanol": ("32.042", "1.5255", "188.90", "1.4"),
        "ethanol": ("46.069", "2.3827", "198.24", "1.6"),
    }
    for component_id, (mass, segments, energy, solvation) in expected_solvents.items():
        assert singles[(component_id, "molar_mass")].value.to(
            "gram / mole"
        ).magnitude == Decimal(mass)
        assert singles[(component_id, "segment_count")].value == Decimal(segments)
        assert singles[(component_id, "dispersion_energy_over_k")].value.to(
            "kelvin"
        ).magnitude == Decimal(energy)
        assert singles[(component_id, "charge_number")].value == 0
        assert singles[(component_id, "solvation_factor")].value == Decimal(solvation)

    pair_records = {
        frozenset((record.component_id_a, record.component_id_b)): record
        for record in bundle.records
        if isinstance(record, PairParameterRecord)
    }
    pairs = {key: record.value for key, record in pair_records.items()}
    assert pairs == {
        frozenset(("hydrogen-cation", "chloride-anion")): Decimal("-0.9"),
        frozenset(("lithium-cation", "chloride-anion")): Decimal("0.8"),
        frozenset(("sodium-cation", "chloride-anion")): Decimal("0.8"),
        frozenset(("potassium-cation", "chloride-anion")): Decimal("0"),
        frozenset(("vanadium-iii-cation", "chloride-anion")): Decimal("1"),
        frozenset(("hydrogen-cation", "bromide-anion")): Decimal("-0.7"),
        frozenset(("lithium-cation", "bromide-anion")): Decimal("0.5"),
        frozenset(("sodium-cation", "bromide-anion")): Decimal("0.65"),
        frozenset(("potassium-cation", "bromide-anion")): Decimal("-0.35"),
        frozenset(("sodium-cation", "iodide-anion")): Decimal("0.45"),
        frozenset(("potassium-cation", "iodide-anion")): Decimal("0"),
        frozenset(("lithium-cation", "sulfate-anion")): Decimal("0.8"),
        frozenset(("sodium-cation", "sulfate-anion")): Decimal("0.7"),
        frozenset(("potassium-cation", "sulfate-anion")): Decimal("-0.6"),
        frozenset(("vanadyl-cation", "sulfate-anion")): Decimal("0"),
        frozenset(("water", "chloride-anion")): Decimal("-0.3"),
        frozenset(("methanol", "chloride-anion")): Decimal("0.5"),
        frozenset(("ethanol", "chloride-anion")): Decimal("0.8"),
        frozenset(("water", "bromide-anion")): Decimal("-0.3"),
        frozenset(("methanol", "bromide-anion")): Decimal("0.15"),
        frozenset(("ethanol", "bromide-anion")): Decimal("0"),
        frozenset(("water", "iodide-anion")): Decimal("-0.05"),
        frozenset(("methanol", "iodide-anion")): Decimal("0.37"),
        frozenset(("ethanol", "iodide-anion")): Decimal("0.18"),
        frozenset(("water", "sulfate-anion")): Decimal("-0.1"),
        frozenset(("water", "hydrogen-cation")): Decimal("0"),
        frozenset(("methanol", "hydrogen-cation")): Decimal("-0.3"),
        frozenset(("ethanol", "hydrogen-cation")): Decimal("-0.6"),
        frozenset(("water", "lithium-cation")): Decimal("-0.4"),
        frozenset(("methanol", "lithium-cation")): Decimal("-0.9"),
        frozenset(("ethanol", "lithium-cation")): Decimal("-0.8"),
        frozenset(("water", "sodium-cation")): Decimal("-0.3"),
        frozenset(("methanol", "sodium-cation")): Decimal("-0.25"),
        frozenset(("ethanol", "sodium-cation")): Decimal("0.05"),
        frozenset(("water", "potassium-cation")): Decimal("-0.1"),
        frozenset(("methanol", "potassium-cation")): Decimal("0.32"),
        frozenset(("ethanol", "potassium-cation")): Decimal("0.53"),
        frozenset(("water", "vanadyl-cation")): Decimal("-0.08"),
        frozenset(("water", "vanadium-iii-cation")): Decimal("-0.3"),
        frozenset(("water", "methanol")): Decimal("-0.0878"),
        frozenset(("water", "ethanol")): Decimal("-0.0617"),
    }
    assert len(pairs) == 41
    assert all(record.source_id == "figiel-yu-held-2025" for record in pair_records.values())
    assert all(
        record.locator.startswith(("Table 4,", "Table 5,"))
        for record in pair_records.values()
    )
    for blank_pair in (
        ("hydrogen-cation", "iodide-anion"),
        ("vanadyl-cation", "chloride-anion"),
        ("vanadium-iii-cation", "sulfate-anion"),
        ("methanol", "sulfate-anion"),
    ):
        assert frozenset(blank_pair) not in pairs

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
    assert len(associations) == 18
    assert associations[
        ("methanol", "b", "water", "a", "association_energy_over_k")
    ].to("kelvin").magnitude == Decimal("2662.6")
    assert associations[
        ("ethanol", "b", "methanol", "a", "association_volume")
    ] == Decimal("0.03374752517428578")
    association_records = [
        record
        for record in bundle.records
        if isinstance(record, AssociationParameterRecord)
    ]
    assert sum(record.source_id == "figiel-yu-held-2025" for record in association_records) == 6
    assert sum(record.source_id == "gross-sadowski-2002" for record in association_records) == 12

    correlations = {
        (record.component_id, record.family): record
        for record in bundle.records
        if isinstance(
            record,
            (ConstantCorrelation, ConstantPlusSumOfExponentialsCorrelation),
        )
    }
    assert set(correlations) == {
        (component_id, family)
        for component_id in ("water", "methanol", "ethanol")
        for family in ("segment_diameter", "relative_permittivity")
    }
    assert correlations[("methanol", "segment_diameter")].constant.to(
        "angstrom"
    ).magnitude == Decimal("3.2300")
    assert correlations[("ethanol", "segment_diameter")].constant.to(
        "angstrom"
    ).magnitude == Decimal("3.1771")
    assert correlations[("water", "relative_permittivity")].constant == Decimal(
        "78.09"
    )
    assert correlations[("methanol", "relative_permittivity")].constant == Decimal(
        "33.05"
    )
    assert correlations[("ethanol", "relative_permittivity")].constant == Decimal(
        "24.88"
    )
    water_sigma = correlations[("water", "segment_diameter")]
    assert isinstance(water_sigma, ConstantPlusSumOfExponentialsCorrelation)
    assert water_sigma.constant.to("angstrom").magnitude == Decimal("2.7927")
    assert tuple(
        (
            term.amplitude.to("angstrom").magnitude,
            term.exponent_coefficient.to("1 / kelvin").magnitude,
        )
        for term in water_sigma.terms
    ) == (
        (Decimal("10.11"), Decimal("-0.01775")),
        (Decimal("-1.417"), Decimal("-0.01146")),
    )

    models = {
        record.family: record.value
        for record in bundle.records
        if isinstance(record, ModelParameterRecord)
    }
    assert models == {
        "ionic_region_relative_permittivity": Decimal("8"),
        "ion_fraction_suppression_coefficient": Decimal("7.01"),
        "relative_permittivity_formulation": "ion-fraction-suppression",
    }

    domain = next(domain for domain in bundle.domains if domain.domain_id == "figiel-298-k-1-bar")
    assert domain.temperature_min.to("kelvin").magnitude == Decimal("298.15")
    assert domain.temperature_max.to("kelvin").magnitude == Decimal("298.15")
    assert domain.pressure_min.to("bar").magnitude == Decimal("1")
    assert domain.pressure_max.to("bar").magnitude == Decimal("1")
    assert domain.ion_mole_fraction_max == Decimal("0.38")
    assert {source.source_id: source.doi for source in bundle.sources}[
        "figiel-yu-held-2025"
    ] == "10.1021/acs.iecr.5c00475"
    sources = {source.source_id: source for source in bundle.sources}
    assert (
        "a3c940895c530f72f47f22f8c0f4796b5ad37918a4edd274f85497c6dd81ad1f"
        in sources["figiel-yu-held-2025"].use_basis
    )
    assert (
        "005b38ed566ec3c09b87e1ca3a9dd6eeafc9ba75e1a30b9322291d770bb93895"
        in sources["figiel-yu-held-2025"].use_basis
    )
    assert sources["gross-sadowski-2002"].doi == "10.1021/ie010954d"
    assert bundle.fingerprint == (
        "sha256:0c5f7b2ddeae9ca25357f269e1a2b135ace74eae7eb28c1c262fa1419ca487a6"
    )

    selections = {
        ("water", "sodium-cation", "chloride-anion"): (
            "sha256:10279776fc2a7515bbb79404a2da7356994e11861ef28a55d67abb00509408e3",
            (0.90, 0.05, 0.05),
        ),
        ("water", "vanadyl-cation", "sulfate-anion"): (
            "sha256:56a36f633eb5b47cd188eae864e16c26352b12fe28b5ddf2540e065195bd89ec",
            (0.96, 0.02, 0.02),
        ),
        ("water", "vanadium-iii-cation", "chloride-anion"): (
            "sha256:47f73ac9ec45616269024ec2fb1ffe82f4ba8f7b83dc984612651917d410aad6",
            (0.96, 0.01, 0.03),
        ),
    }
    for component_ids, (fingerprint, mole_fractions) in selections.items():
        parameters = bundle.select(component_ids)
        assert parameters.fingerprint == fingerprint
        state = epcsaft.Mixture(parameters).state(
            T=298.15 * unit_registry.kelvin,
            P=1 * unit_registry.bar,
            x=mole_fractions,
            phase="liquid",
        )
        assert state.ares().check("[energy] / [substance]")

    forward = epcsaft.Mixture(
        bundle.select(("water", "vanadium-iii-cation", "chloride-anion"))
    ).state(
        T=298.15 * unit_registry.kelvin,
        rho=1000 * unit_registry.mole / unit_registry.meter**3,
        x=(0.96, 0.01, 0.03),
    )
    reverse = epcsaft.Mixture(
        bundle.select(("chloride-anion", "vanadium-iii-cation", "water"))
    ).state(
        T=298.15 * unit_registry.kelvin,
        rho=1000 * unit_registry.mole / unit_registry.meter**3,
        x=(0.03, 0.01, 0.96),
    )
    assert forward.ares().to("joule / mole").magnitude == pytest.approx(
        reverse.ares().to("joule / mole").magnitude,
        rel=1e-13,
        abs=1e-12,
    )
    with pytest.raises(ParameterError) as exc_info:
        epcsaft.Mixture(
            bundle.select(("water", "hydrogen-cation", "iodide-anion"))
        )
    assert exc_info.value.code == "missing-record"


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
        "crespo-2021-deea/1/association.csv",
        "crespo-2021-deea/1/bundle.toml",
        "crespo-2021-deea/1/components.csv",
        "crespo-2021-deea/1/model.csv",
        "crespo-2021-deea/1/single.csv",
        "crespo-2021-deea/1/sites.csv",
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
        "figiel-yu-held-2025/1/association.csv",
        "figiel-yu-held-2025/1/bundle.toml",
        "figiel-yu-held-2025/1/components.csv",
        "figiel-yu-held-2025/1/correlations.toml",
        "figiel-yu-held-2025/1/model.csv",
        "figiel-yu-held-2025/1/pair.csv",
        "figiel-yu-held-2025/1/single.csv",
        "figiel-yu-held-2025/1/sites.csv",
        "gross-2001-methane-ethane/1/bundle.toml",
        "gross-2001-methane-ethane/1/components.csv",
        "gross-2001-methane-ethane/1/model.csv",
        "gross-2001-methane-ethane/1/pair.csv",
        "gross-2001-methane-ethane/1/single.csv",
        "gross-2001-propane/1/bundle.toml",
        "gross-2001-propane/1/components.csv",
        "gross-2001-propane/1/model.csv",
        "gross-2001-propane/1/single.csv",
        "gross-sadowski-2001/1/bundle.toml",
        "gross-sadowski-2001/1/components.csv",
        "gross-sadowski-2001/1/model.csv",
        "gross-sadowski-2001/1/pair.csv",
        "gross-sadowski-2001/1/single.csv",
        "held-2008-water-self-ionization/1/association.csv",
        "held-2008-water-self-ionization/1/bundle.toml",
        "held-2008-water-self-ionization/1/components.csv",
        "held-2008-water-self-ionization/1/correlations.toml",
        "held-2008-water-self-ionization/1/model.csv",
        "held-2008-water-self-ionization/1/pair.csv",
        "held-2008-water-self-ionization/1/single.csv",
        "held-2008-water-self-ionization/1/sites.csv",
        "provenance_names.toml",
    ]


def test_crespo_2021_deea_catalog_is_exact_and_source_complete() -> None:
    bundle = _ParameterBundle.from_catalog("crespo-2021-deea", version=1)
    assert bundle.fingerprint == (
        "sha256:e04e4799b7f3729e989e031feaf6e2a87c97f031bcd538d1bb09781f714ce65c"
    )
    assert bundle.purpose == "installed-catalog"
    assert tuple((source.source_id, source.doi, source.use_basis) for source in bundle.sources) == (
        (
            "crespo-et-al-2021",
            "10.1016/j.fluid.2021.113277",
            (
                "Table 3 neutral DEEA PC-SAFT parameters and 2B association topology; "
                "kappa_HB is used in the cited PC-SAFT sigma^3*kappa association-strength "
                "convention; c_i is excluded as the PIL volume correction"
            ),
        ),
    )
    assert tuple((domain.domain_id, domain.kind) for domain in bundle.domains) == (
        ("crespo-deea-fit-domain-unreported", "unknown"),
    )
    assert tuple(
        (component.component_id, component.name, component.aliases)
        for component in bundle.components
    ) == (
        (
            "diethylethanolamine",
            "N,N-diethylethanolamine",
            ("2-(diethylamino)ethanol", "DEEA"),
        ),
    )

    parameters = epcsaft.Parameters.from_catalog(
        "crespo-2021-deea",
        components=("diethylethanolamine",),
        version=1,
    )
    assert parameters.bundle_fingerprint == bundle.fingerprint
    assert epcsaft.Mixture(parameters)._model.topology == "associating"
    assert {record.record_id for record in parameters.records} == {
        "crespo-relative-permittivity-formulation",
        "deea-association-energy",
        "deea-association-volume",
        "deea-dispersion-energy",
        "deea-molar-mass",
        "deea-segment-count",
        "deea-segment-diameter",
        "deea-site-a",
        "deea-site-b",
    }
    assert {
        (record.site_id, record.site_class, record.multiplicity)
        for record in parameters.records
        if isinstance(record, SiteRecord)
    } == {("a", "a", 1), ("b", "b", 1)}

    singles = {
        record.family: record
        for record in parameters.records
        if isinstance(record, SingleParameterRecord)
    }
    assert singles["molar_mass"].value.to("gram / mole").magnitude == Decimal("117.19")
    assert singles["segment_count"].value == Decimal("3.844")
    assert singles["segment_diameter"].value.to("angstrom").magnitude == Decimal("3.6406")
    assert singles["dispersion_energy_over_k"].value.to("kelvin").magnitude == Decimal("248.55")

    associations = {
        record.family: record
        for record in parameters.records
        if isinstance(record, AssociationParameterRecord)
    }
    assert {
        (
            record.component_id_a,
            record.site_id_a,
            record.component_id_b,
            record.site_id_b,
        )
        for record in associations.values()
    } == {("diethylethanolamine", "a", "diethylethanolamine", "b")}
    assert associations["association_energy_over_k"].value.to("kelvin").magnitude == Decimal(
        "2033.2"
    )
    assert associations["association_volume"].value == Decimal("0.02412")
    assert {
        record.family: record.value
        for record in parameters.records
        if isinstance(record, ModelParameterRecord)
    } == {"relative_permittivity_formulation": "none"}
    assert all(record.source_id == "crespo-et-al-2021" for record in parameters.records)
    assert all(
        record.domain_id == "crespo-deea-fit-domain-unreported" for record in parameters.records
    )


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
        source for source in bundle.sources if source.source_id == "held-cameretti-sadowski-2008"
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

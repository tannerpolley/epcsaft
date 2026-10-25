"""Validate the frozen Figiel 2025 source/target ledger without package imports."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from decimal import Decimal
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LEDGER = ROOT / "inputs/source" / "figiel-2025-regression-target-ledger.csv"
DEFAULT_METADATA = ROOT / "inputs/source" / "figiel-2025-regression-target-ledger.yaml"
DEFAULT_HAMER = (
    ROOT.parents[2]
    / "data/packets/hamer-wu-1972-aqueous-alkali-halides/1/observations"
    / "hamer-wu-1972-aqueous-alkali-halides.csv"
)
DEFAULT_PARAMETERS = ROOT / "inputs/source" / "figiel-2025-parameter-provenance.csv"

EXPECTED_LEDGER_SHA256 = "4d6ce439f732d44c6cee08d9d800cef6f844eb297a30ebf935daa1c65403c5e7"
EXPECTED_HAMER_SHA256 = "2f63e13f06a5b0f4e8bca2980b6a8d9d7fb0f839153c43e3a71952daf9796595"
EXPECTED_SI_SHA256 = "005b38ed566ec3c09b87e1ca3a9dd6eeafc9ba75e1a30b9322291d770bb93895"
EXPECTED_PARAMETER_SHA256 = "932e8baa90fcefbaa8c3a8730cdeadd83a4c01f0a3b109f4e4cd0319aee9312b"
EXPECTED_COLUMNS = [
    "target_id",
    "regression_order",
    "role",
    "quantity",
    "value",
    "unit",
    "species",
    "solvent",
    "T_K",
    "P_bar",
    "composition_basis",
    "composition_value",
    "basis",
    "source_classification",
    "source_id",
    "source_locator",
    "source_row",
    "underlying_reference",
    "extraction",
    "source_sha256",
    "uncertainty_value",
    "uncertainty_unit",
    "uncertainty_meaning",
    "tracer_status",
]
EXPECTED_CLASSIFICATIONS = {
    "direct": 191,
    "table-derived": 47,
    "digitized": 90,
    "constructed": 43,
    "missing-primary-source": 36,
}
EXPECTED_QUANTITIES = {
    "stoichiometric_mean_ionic_activity_coefficient": 200,
    "single_ion_solvation_gibbs_energy": 90,
    "single_ion_transfer_gibbs_energy": 43,
    "relative_dielectric_constant": 29,
    "relative_dielectric_constant_ratio": 18,
    "solution_mass_density": 22,
    "osmotic_coefficient": 5,
}
EXPECTED_S5_AVERAGES = {
    "H+": Decimal("-1094.5"),
    "Li+": Decimal("-486.2"),
    "Na+": Decimal("-381.1"),
    "K+": Decimal("-309.1"),
    "Cl-": Decimal("-314.9"),
    "Br-": Decimal("-290.9"),
    "I-": Decimal("-252.6"),
    "SO4^2-": Decimal("-1062.1"),
    "VO^2+": Decimal("-1895.0"),
    "V^3+": Decimal("-4202.1"),
}
EXPECTED_S5_UNDERLYING_COUNTS = {
    "H+": 1,
    "Li+": 5,
    "Na+": 5,
    "K+": 5,
    "Cl-": 6,
    "Br-": 6,
    "I-": 6,
    "SO4^2-": 1,
    "VO^2+": 1,
    "V^3+": 1,
}
PARAMETER_COLUMNS = [
    "cell_id",
    "source_location",
    "parameter_family",
    "component_i",
    "component_j",
    "parameter",
    "reported_value",
    "unit",
    "temperature_K",
    "provenance_classification",
    "recovery_target",
    "staged_order",
    "source_reference",
    "required_recovery_dataset",
    "notes",
    "source_artifact_sha256",
]
EXPECTED_PARAMETER_PROVENANCE = {
    "fitted-in-Figiel-2025": 57,
    "inherited-from-cited-source": 33,
    "fixed-or-assumed": 13,
    "blank-or-unpublished": 19,
}
EXPECTED_PARAMETER_STAGES = {
    "A": 5,
    "B": 1,
    "C": 3,
    "D": 11,
    "E": 10,
    "F": 27,
    "G": 48,
    "H": 17,
}
MAIN_MARKDOWN_SHA256 = "ce80533925a91bc59d8d0d8056113c40611ca26c2edf04aced76986d50bd4bae"
EXPECTED_FITTED_ZERO_CELLS = {
    "t4-K+-Cl-",
    "t4-K+-I-",
    "t4-VO^2+-SO4^2-",
    "t5-Br--ethanol",
    "t5-H+-water",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path, expected_columns: list[str] | None = None) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if expected_columns is not None and reader.fieldnames != expected_columns:
            raise ValueError(f"unexpected columns in {path.name}: {reader.fieldnames}")
        return list(reader)


def check_hash_token(value: str) -> None:
    parts = value.split("+")
    if not parts or any(len(part) != 64 or any(char not in "0123456789abcdef" for char in part) for part in parts):
        raise ValueError(f"invalid source SHA-256 token: {value}")


def check_hamer_rows(rows: list[dict[str, str]], hamer_path: Path) -> None:
    if sha256(hamer_path) != EXPECTED_HAMER_SHA256:
        raise ValueError("retained Hamer-Wu source hash changed")
    hamer = read_csv(hamer_path)
    ledger = [row for row in rows if row["source_id"] == "hamer-wu-1972-aqueous-alkali-halides"]
    if len(hamer) != 164 or len(ledger) != len(hamer):
        raise ValueError("expected all 164 Hamer-Wu aqueous MIAC rows")
    observed = [
        (row["species"], row["composition_value"], row["value"])
        for row in ledger
    ]
    expected = [
        (row["salt"], row["molality_mol_kg"], row["gamma_pm_m"])
        for row in hamer
    ]
    if observed != expected:
        raise ValueError("Figiel ledger does not reproduce the retained Hamer-Wu packet")
    if any(row["tracer_status"] != "later_capability_gate" for row in ledger):
        raise ValueError("Hamer-Wu MIAC rows must remain outside the first Born tracer")


def check_s5(rows: list[dict[str, str]]) -> None:
    reported_targets = [row for row in rows if row["role"] == "reported_average_target"]
    reported_later = [row for row in rows if row["role"] == "reported_average_later_gate"]
    reported = reported_targets + reported_later
    underlying = [row for row in rows if row["role"] == "underlying_literature_value"]
    if len(reported_targets) != 5 or len(reported_later) != 5 or len(underlying) != 37:
        raise ValueError("unexpected SI Table S5 row count")
    observed_averages = {row["species"]: Decimal(row["value"]) for row in reported}
    if observed_averages != EXPECTED_S5_AVERAGES:
        raise ValueError("SI Table S5 reported averages changed")
    if Counter(row["species"] for row in underlying) != EXPECTED_S5_UNDERLYING_COUNTS:
        raise ValueError("SI Table S5 underlying literature columns changed")
    if any(row["source_sha256"] != EXPECTED_SI_SHA256 for row in reported + underlying):
        raise ValueError("SI Table S5 rows are not bound to the official SI")
    if any(row["source_classification"] != "table-derived" for row in reported + underlying):
        raise ValueError("SI Table S5 classifications changed")
    first_ions = {"Li+", "Na+", "K+", "Cl-", "Br-"}
    if {row["species"] for row in reported_targets} != first_ions:
        raise ValueError("only current-catalog SI averages may be Regression residual targets")
    first_support = [row for row in reported + underlying if row["species"] in first_ions]
    if len(first_support) != 32 or any(
        row["tracer_status"] != "first_tracer_support" for row in first_support
    ):
        raise ValueError("first tracer must retain exactly 32 current-ion SI S5 rows")
    if any(
        row["tracer_status"] != "later_capability_gate"
        for row in reported + underlying
        if row["species"] not in first_ions
    ):
        raise ValueError("expanded-ion SI S5 rows must remain later capability gates")


def check_constructed_rows(rows: list[dict[str, str]]) -> None:
    transfers = [row for row in rows if row["role"] == "ion_organic_solvation_construction_input"]
    constructed = [row for row in rows if row["role"] == "constructed_ion_solvation_target"]
    if len(transfers) != 43 or len(constructed) != 43:
        raise ValueError("expected paired Figure 6 transfer and constructed solvation rows")
    transfer_by_state = {
        (row["species"], row["solvent"], row["composition_value"]): Decimal(row["value"])
        for row in transfers
    }
    for row in constructed:
        state = (row["species"], row["solvent"], row["composition_value"])
        expected = EXPECTED_S5_AVERAGES[row["species"]] + transfer_by_state[state]
        if Decimal(row["value"]) != expected:
            raise ValueError(f"equation-19 construction changed at {state}")
        if row["source_classification"] != "constructed":
            raise ValueError("constructed solvation target classification changed")


def check_si_s4(rows: list[dict[str, str]]) -> None:
    density = [row for row in rows if row["quantity"] == "solution_mass_density"]
    osmotic = [row for row in rows if row["quantity"] == "osmotic_coefficient"]
    if len(density) != 22 or len(osmotic) != 5:
        raise ValueError("unexpected SI Table S4 target count")
    expected_density_temperatures = Counter({"298.15": 11, "303.15": 11})
    if Counter(row["T_K"] for row in density) != expected_density_temperatures:
        raise ValueError("SI Table S4 temperature partition changed")
    if Counter(row["species"] for row in density) != {"VCl3": 12, "VOSO4": 10}:
        raise ValueError("SI Table S4 salt partition changed")
    if any(row["source_classification"] != "direct" for row in density + osmotic):
        raise ValueError("SI Table S4 direct-data classification changed")
    if any(row["source_sha256"] != EXPECTED_SI_SHA256 for row in density + osmotic):
        raise ValueError("SI Table S4 rows are not bound to the official SI")


def check_parameter_provenance(parameter_path: Path, metadata: dict[str, object]) -> dict[str, int]:
    parameter_hash = sha256(parameter_path)
    if parameter_hash != EXPECTED_PARAMETER_SHA256:
        raise ValueError("Figiel parameter-provenance hash changed")
    rows = read_csv(parameter_path, PARAMETER_COLUMNS)
    if len(rows) != 122 or len({row["cell_id"] for row in rows}) != 122:
        raise ValueError("expected 122 unique parameter-provenance cells")
    provenance = Counter(row["provenance_classification"] for row in rows)
    stages = Counter(row["staged_order"] for row in rows)
    targets = Counter(row["recovery_target"] for row in rows)
    if provenance != EXPECTED_PARAMETER_PROVENANCE:
        raise ValueError(f"parameter provenance changed: {provenance}")
    if stages != EXPECTED_PARAMETER_STAGES:
        raise ValueError(f"parameter staging changed: {stages}")
    if targets != {"true": 90, "false": 32}:
        raise ValueError(f"parameter recovery-target partition changed: {targets}")
    if any(row["source_artifact_sha256"] != MAIN_MARKDOWN_SHA256 for row in rows):
        raise ValueError("parameter cell is not bound to the supplied main-paper artifact")

    for row in rows:
        value = row["reported_value"]
        provenance_value = row["provenance_classification"]
        target = row["recovery_target"]
        if not value and (provenance_value != "blank-or-unpublished" or target != "false"):
            raise ValueError("blank parameter cell was inferred or promoted")
        if value == "0" and (
            provenance_value != "fitted-in-Figiel-2025" or target != "true"
        ):
            raise ValueError("explicit published fitted zero lost its recovery-target identity")
        if provenance_value in {"fitted-in-Figiel-2025", "inherited-from-cited-source"}:
            if target != "true" or not value:
                raise ValueError("source-backed recovery parameter lost its target identity")
        if provenance_value in {"fixed-or-assumed", "blank-or-unpublished"} and target != "false":
            raise ValueError("fixed or unpublished parameter became a recovery target")
    if {row["cell_id"] for row in rows if row["reported_value"] == "0"} != (
        EXPECTED_FITTED_ZERO_CELLS
    ):
        raise ValueError("explicit published fitted-zero cell set changed")

    stage_a = [row for row in rows if row["staged_order"] == "A"]
    if {
        (row["component_i"], row["parameter"], row["provenance_classification"])
        for row in stage_a
    } != {
        (ion, "born_diameter", "fitted-in-Figiel-2025")
        for ion in ("Li+", "Na+", "K+", "Cl-", "Br-")
    }:
        raise ValueError("first tracer widened beyond current-catalog Born diameters")
    stage_b = [row for row in rows if row["staged_order"] == "B"]
    if len(stage_b) != 1 or stage_b[0]["reported_value"] != "7.01":
        raise ValueError("dielectric suppression stage changed")
    table4 = [row for row in rows if row["source_location"].startswith("main paper Table 4")]
    table5 = [row for row in rows if row["source_location"].startswith("main paper Table 5")]
    if len(table4) != 24 or len(table5) != 36:
        raise ValueError("Tables 4-5 cell coverage changed")

    retained = metadata["retained_files"]["parameter_provenance"]
    audit = metadata["parameter_provenance_audit"]
    if (
        retained["rows"] != len(rows)
        or retained["sha256"] != parameter_hash
        or audit["rows"] != len(rows)
        or audit["sha256"] != parameter_hash
    ):
        raise ValueError("metadata does not bind the parameter-provenance ledger")
    return {
        "rows": len(rows),
        "recovery_targets": targets["true"],
        "non_targets": targets["false"],
    }


def check(
    ledger_path: Path = DEFAULT_LEDGER,
    metadata_path: Path = DEFAULT_METADATA,
    hamer_path: Path = DEFAULT_HAMER,
    parameter_path: Path = DEFAULT_PARAMETERS,
) -> dict[str, object]:
    ledger_hash = sha256(ledger_path)
    if ledger_hash != EXPECTED_LEDGER_SHA256:
        raise ValueError("Figiel target ledger hash changed")
    rows = read_csv(ledger_path, EXPECTED_COLUMNS)
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    if len(rows) != 407 or len({row["target_id"] for row in rows}) != 407:
        raise ValueError("expected 407 unique Figiel target rows")
    if metadata["citation"]["main_doi"] != "10.1021/acs.iecr.5c00475":
        raise ValueError("unexpected main-paper DOI")
    if metadata["citation"]["supporting_information_doi"] != "10.1021/acs.iecr.5c00475.s001":
        raise ValueError("unexpected SI DOI")
    migration = metadata["migration_binding"]
    if migration != {
        "decision": "D-023",
        "gate_commit": "9ceb1dae796b44fad7789762efcca9c11491e9c2",
        "gate_tree": "6757a261bffb9486a25f45194bde870605094bb0",
    }:
        raise ValueError("unexpected Migration D-023 gate binding")
    official_si = metadata["official_sources"]["official_si_pdf"]
    if (
        official_si["figshare_article_id"] != 28888514
        or official_si["figshare_file_id"] != 54030209
        or official_si["filename"] != "ie5c00475_si_001.pdf"
        or official_si["size_bytes"] != 396182
        or official_si["md5"] != "4c3d3964d6e61b908c3ce3b17bb4e6b0"
        or official_si["sha256"] != EXPECTED_SI_SHA256
    ):
        raise ValueError("official SI identity or byte binding changed")
    retained = metadata["retained_files"]["target_ledger"]
    if retained["rows"] != len(rows) or retained["sha256"] != ledger_hash:
        raise ValueError("metadata does not bind the retained target ledger")

    classifications = Counter(row["source_classification"] for row in rows)
    quantities = Counter(row["quantity"] for row in rows)
    tracers = Counter(row["tracer_status"] for row in rows)
    if classifications != EXPECTED_CLASSIFICATIONS:
        raise ValueError(f"source classifications changed: {classifications}")
    if quantities != EXPECTED_QUANTITIES:
        raise ValueError(f"quantity counts changed: {quantities}")
    if tracers != {"first_tracer_support": 32, "later_capability_gate": 375}:
        raise ValueError(f"tracer partition changed: {tracers}")

    for row in rows:
        value = Decimal(row["value"])
        temperature = Decimal(row["T_K"])
        pressure = Decimal(row["P_bar"])
        if not value.is_finite() or not temperature.is_finite() or not pressure.is_finite():
            raise ValueError(f"nonfinite retained value in {row['target_id']}")
        if temperature <= 0 or pressure <= 0:
            raise ValueError(f"invalid state in {row['target_id']}")
        if row["uncertainty_value"] or row["uncertainty_unit"]:
            raise ValueError("the packet must not invent pointwise uncertainty")
        if not row["uncertainty_meaning"]:
            raise ValueError("every value must state its uncertainty meaning")
        check_hash_token(row["source_sha256"])

    check_hamer_rows(rows, hamer_path)
    check_s5(rows)
    check_constructed_rows(rows)
    check_si_s4(rows)
    parameter_report = check_parameter_provenance(parameter_path, metadata)

    dielectric = [row for row in rows if row["regression_order"] == "e"]
    if Counter(row["role"] for row in dielectric) != {
        "dielectric_correlation_fit": 36,
        "dielectric_validation": 11,
    }:
        raise ValueError("dielectric fit/validation roles changed")
    if any(row["solvent"] != "ethanol" for row in dielectric if row["role"] == "dielectric_validation"):
        raise ValueError("only ethanol may carry the dielectric validation role")

    comparison = metadata["comparison_contract"]
    if any(
        comparison[name] is not None
        for name in (
            "model_acceptance_threshold",
            "uncertainty_multiplier",
            "fitted_parameter_values",
        )
    ):
        raise ValueError("source packet must not create model thresholds or fitted parameters")
    if len(metadata["missing_primary_source_requests"]) != 19:
        raise ValueError("missing-primary-source request list changed")
    if metadata["source_readiness"]["status"] != "ready_for_provider_regression_design":
        raise ValueError("unexpected source-readiness status")
    markdown = metadata["official_sources"]["retained_lab_main_markdown"]
    if (
        markdown["size_bytes"] != 71826
        or markdown["sha256"] != MAIN_MARKDOWN_SHA256
        or len(markdown["execution_locators"]) != 2
    ):
        raise ValueError("main-paper read-only source identity changed")
    if set(metadata["staged_parameter_recovery"]) != set("ABCDEFGH"):
        raise ValueError("staged recovery order changed")
    proof = metadata["eventual_recovery_proof_contract"]
    if proof["separate_statuses"] != ["solver", "numerical", "physical", "predictive"]:
        raise ValueError("recovery status separation changed")

    return {
        "status": "source_packet_ready",
        "main_doi": metadata["citation"]["main_doi"],
        "si_doi": metadata["citation"]["supporting_information_doi"],
        "si_locator": metadata["official_sources"]["official_si_pdf"]["url"],
        "official_si_sha256": EXPECTED_SI_SHA256,
        "ledger_rows": len(rows),
        "classifications": dict(sorted(classifications.items())),
        "tracer_partition": dict(sorted(tracers.items())),
        "missing_primary_source_requests": len(metadata["missing_primary_source_requests"]),
        "ledger_sha256": ledger_hash,
        "parameter_provenance_sha256": sha256(parameter_path),
        "parameter_provenance": parameter_report,
        "metadata_sha256": sha256(metadata_path),
        "hamer_sha256": sha256(hamer_path),
        "comparison_threshold": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate the frozen Figiel 2025 source packet.")
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--hamer", type=Path, default=DEFAULT_HAMER)
    parser.add_argument("--parameters", type=Path, default=DEFAULT_PARAMETERS)
    args = parser.parse_args()
    print(
        json.dumps(
            check(args.ledger, args.metadata, args.hamer, args.parameters),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

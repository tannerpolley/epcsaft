from __future__ import annotations

import argparse
import csv
import hashlib
import json
from decimal import Decimal
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TARGET = ROOT / "inputs/source" / "glos-2004-propane-saturation.csv"
DEFAULT_METADATA = ROOT / "inputs/source" / "glos-2004-propane-saturation.yaml"
DEFAULT_SOURCE_result = ROOT / "inputs/source" / "glos-2004-propane-saturation-source.csv"
DEFAULT_AUXILIARY_CSV = ROOT / "inputs/source" / "nist-srd69-propane-reference-eos.csv"
DEFAULT_AUXILIARY_METADATA = ROOT / "inputs/source" / "nist-srd69-propane-reference-eos.yaml"

EXPECTED_SOURCE_result_SHA256 = (
    "ed5eb703ccd3e6bb4c4cfa82ecd58c58f9da0c93ab07a204dee94d8b0ae8d081"
)
EXPECTED_TARGET_SHA256 = "ccd1cfa15ec44432b06cbf22316d168c61b282631c9b1e1591e497b8d48b5676"
EXPECTED_AUXILIARY_SHA256 = (
    "4a05958854ed5fed28ce4f2b3dc1e98d2d5f833c3eb2f1290ce6428492f0051c"
)
EXPECTED_OFFICIAL_SOURCE_HASHES = {
    "nist_thermoml_json": "322495c5a01c003e83376e5bad544c3abced330d5054ff0411a7a00b70a963c9",
    "nist_thermoml_xml": "1b2e47d4cafff0f21cf7779d8d01b522bc2fa8d885ce4d6ebc04c151e0504829",
}
SOURCE_COLUMNS = ["dataset_number", "T_K", "value", "expanded_uncertainty"]
TARGET_COLUMNS = [
    "row_id",
    "component_id",
    "role",
    "T_K",
    "p_sat_Pa",
    "p_sat_expanded_uncertainty_Pa",
    "rho_sat_liq_kg_m3",
    "rho_sat_liq_expanded_uncertainty_kg_m3",
    "rho_sat_vap_kg_m3",
    "rho_sat_vap_expanded_uncertainty_kg_m3",
]
EXPECTED_SOURCE_TEMPERATURES = {
    1: tuple(range(110, 341, 10)),
    2: (90, 95, 100, *range(110, 341, 10)),
    3: tuple(range(230, 341, 10)),
}
EXPECTED_TRAINING = (150, 210, 270, 330)
EXPECTED_STRESS = (110, 340)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_hash(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def read_csv(path: Path, columns: list[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != columns:
            raise ValueError(f"unexpected columns in {path.name}: {reader.fieldnames}")
        return list(reader)


def check_contract_hash(metadata: dict[str, object], name: str) -> None:
    contract = metadata[name]
    if contract["sha256"] != canonical_hash(contract["payload"]):
        raise ValueError(f"{name} payload hash does not match metadata")


def check_source_result(
    source_result_path: Path,
    metadata: dict[str, object],
) -> dict[int, dict[int, tuple[Decimal, Decimal]]]:
    result_hash = sha256(source_result_path)
    verification = metadata["source_verification_contract"]["payload"]
    if result_hash != EXPECTED_SOURCE_result_SHA256:
        raise ValueError("retained Glos source result hash changed")
    if verification["mandatory_source_result_sha256"] != result_hash:
        raise ValueError("source verification contract does not bind the retained result")
    if metadata["retained_files"]["source_result"]["sha256"] != result_hash:
        raise ValueError("retained source result hash does not match metadata")
    if verification["required_source_rows"] != 63:
        raise ValueError("unexpected required source row count")

    rows = read_csv(source_result_path, SOURCE_COLUMNS)
    if len(rows) != verification["required_source_rows"]:
        raise ValueError("retained source result row count changed")

    values: dict[int, dict[int, tuple[Decimal, Decimal]]] = {1: {}, 2: {}, 3: {}}
    observed_temperatures: dict[int, list[int]] = {1: [], 2: [], 3: []}
    for row in rows:
        dataset = int(row["dataset_number"])
        temperature = int(row["T_K"])
        if dataset not in values or temperature in values[dataset]:
            raise ValueError("unexpected or duplicate ThermoML dataset state")
        value = Decimal(row["value"])
        uncertainty = Decimal(row["expanded_uncertainty"])
        if not value.is_finite() or not uncertainty.is_finite() or value <= 0 or uncertainty <= 0:
            raise ValueError("source values and expanded uncertainties must be positive and finite")
        values[dataset][temperature] = (value, uncertainty)
        observed_temperatures[dataset].append(temperature)

    for dataset, expected in EXPECTED_SOURCE_TEMPERATURES.items():
        if tuple(observed_temperatures[dataset]) != expected:
            raise ValueError(f"ThermoML dataset {dataset} temperature sequence changed")
    return values


def check_target(
    target_path: Path,
    metadata: dict[str, object],
    source: dict[int, dict[int, tuple[Decimal, Decimal]]],
) -> tuple[list[dict[str, str]], dict[str, list[str]]]:
    target_hash = sha256(target_path)
    if target_hash != EXPECTED_TARGET_SHA256:
        raise ValueError("retained Glos target hash changed")
    if metadata["retained_files"]["target_csv"]["sha256"] != target_hash:
        raise ValueError("retained target hash does not match metadata")

    rows = read_csv(target_path, TARGET_COLUMNS)
    if len(rows) != 24:
        raise ValueError("expected 24 matched vapor-pressure target rows")
    expected_temperatures = EXPECTED_SOURCE_TEMPERATURES[1]
    role_ids: dict[str, list[str]] = {"training": [], "held_out": [], "stress": []}
    pressures: list[Decimal] = []
    liquid_densities: list[Decimal] = []
    vapor_densities: list[Decimal] = []
    for row, temperature in zip(rows, expected_temperatures, strict=True):
        row_id = f"glos2004-propane-sat-{temperature}-k"
        expected_role = (
            "training"
            if temperature in EXPECTED_TRAINING
            else "stress"
            if temperature in EXPECTED_STRESS
            else "held_out"
        )
        if (
            row["row_id"] != row_id
            or row["component_id"] != "propane"
            or row["role"] != expected_role
            or row["T_K"] != str(temperature)
        ):
            raise ValueError(f"target identity, role, or temperature changed at {temperature} K")

        pressure, pressure_uncertainty = source[1][temperature]
        liquid_density, liquid_uncertainty = source[2][temperature]
        expected = (
            pressure * Decimal(1000),
            pressure_uncertainty * Decimal(1000),
            liquid_density,
            liquid_uncertainty,
        )
        observed = tuple(
            Decimal(row[name])
            for name in (
                "p_sat_Pa",
                "p_sat_expanded_uncertainty_Pa",
                "rho_sat_liq_kg_m3",
                "rho_sat_liq_expanded_uncertainty_kg_m3",
            )
        )
        if observed != expected:
            raise ValueError(f"target does not reproduce source result at {temperature} K")

        if temperature in source[3]:
            vapor_density, vapor_uncertainty = source[3][temperature]
            if (
                Decimal(row["rho_sat_vap_kg_m3"]),
                Decimal(row["rho_sat_vap_expanded_uncertainty_kg_m3"]),
            ) != (vapor_density, vapor_uncertainty):
                raise ValueError(f"vapor target does not reproduce source result at {temperature} K")
            vapor_densities.append(vapor_density)
        elif row["rho_sat_vap_kg_m3"] or row["rho_sat_vap_expanded_uncertainty_kg_m3"]:
            raise ValueError(f"vapor density was invented at {temperature} K")

        pressures.append(observed[0])
        liquid_densities.append(observed[2])
        role_ids[expected_role].append(row_id)

    if any(right <= left for left, right in zip(pressures[:-1], pressures[1:], strict=True)):
        raise ValueError("saturation pressure must increase across the target grid")
    if any(right >= left for left, right in zip(liquid_densities[:-1], liquid_densities[1:], strict=True)):
        raise ValueError("saturated-liquid density must decrease across the target grid")
    if any(right <= left for left, right in zip(vapor_densities[:-1], vapor_densities[1:], strict=True)):
        raise ValueError("saturated-vapor density must increase where reported")
    return rows, role_ids


def check_auxiliary(auxiliary_csv_path: Path, auxiliary_metadata_path: Path) -> dict[str, object]:
    auxiliary_metadata = json.loads(auxiliary_metadata_path.read_text(encoding="utf-8"))
    auxiliary_hash = sha256(auxiliary_csv_path)
    if auxiliary_hash != EXPECTED_AUXILIARY_SHA256:
        raise ValueError("auxiliary NIST reference-EOS CSV hash changed")
    if auxiliary_metadata["retained_files"]["csv"]["sha256"] != auxiliary_hash:
        raise ValueError("auxiliary NIST CSV hash does not match metadata")
    if auxiliary_metadata["dataset_id"] != "nist-webbook-propane-reference-eos-100-360-k-auxiliary-v1":
        raise ValueError("unexpected auxiliary NIST dataset identity")
    if "fit_target_contract" in auxiliary_metadata:
        raise ValueError("auxiliary NIST data must not define a fit target")
    if "not direct experimental" not in auxiliary_metadata["source_audit"]["source_kind"]:
        raise ValueError("auxiliary NIST source kind is not explicit")
    rows = read_csv(
        auxiliary_csv_path,
        [
            "row_id",
            "component_id",
            "role",
            "T_K",
            "p_sat_Pa",
            "rho_sat_liq_kg_m3",
            "source_pressure_relative_uncertainty",
            "source_liquid_density_relative_uncertainty",
        ],
    )
    if not rows or any(row["role"] != "reference_only" for row in rows):
        raise ValueError("all auxiliary NIST rows must be reference_only")
    check_contract_hash(auxiliary_metadata, "comparison_contract")
    comparison = auxiliary_metadata["comparison_contract"]["payload"]
    if any(
        comparison[name] is not None
        for name in (
            "pressure_comparison_allowance",
            "liquid_density_comparison_allowance",
            "model_accuracy_floor",
        )
    ):
        raise ValueError("auxiliary NIST evidence must not define a model allowance")
    return {"dataset_id": auxiliary_metadata["dataset_id"], "rows": len(rows), "role": "reference_only"}


def check(
    target_path: Path,
    metadata_path: Path,
    source_result_path: Path,
    auxiliary_csv_path: Path,
    auxiliary_metadata_path: Path,
) -> dict[str, object]:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata["dataset_id"] != "glos-2004-experimental-propane-saturation-110-340-k-v1":
        raise ValueError("unexpected experimental target identity")
    if metadata["citation"]["doi"] != "10.1016/j.jct.2004.07.017":
        raise ValueError("unexpected primary-paper DOI")
    if not metadata["citation"]["table_locator"].startswith("Table 2"):
        raise ValueError("unexpected primary-paper table locator")
    source_hashes = {
        name: metadata["sources"][name]["sha256"] for name in EXPECTED_OFFICIAL_SOURCE_HASHES
    }
    if source_hashes != EXPECTED_OFFICIAL_SOURCE_HASHES:
        raise ValueError("official ThermoML source hashes changed")
    for name in ("fit_target_contract", "source_verification_contract", "comparison_contract"):
        check_contract_hash(metadata, name)

    comparison = metadata["comparison_contract"]["payload"]
    if any(value is not None for name, value in comparison.items() if name != "contract_version" and name != "interpretation"):
        raise ValueError("experimental source packet must not define a model comparison allowance")

    source = check_source_result(source_result_path, metadata)
    rows, role_ids = check_target(target_path, metadata, source)
    target = metadata["fit_target_contract"]["payload"]
    for role in ("training", "held_out", "stress"):
        if target[f"{role}_row_ids"] != role_ids[role]:
            raise ValueError(f"{role} partition does not match retained target rows")
    if target["fit_parameter_names"] != [
        "segment_count",
        "segment_diameter_angstrom",
        "dispersion_energy_over_k_kelvin",
    ] or target["fit_parameter_units"] != ["1", "angstrom", "K"]:
        raise ValueError("fit coordinate identity changed")

    auxiliary = check_auxiliary(auxiliary_csv_path, auxiliary_metadata_path)
    return {
        "status": "accepted",
        "dataset_id": metadata["dataset_id"],
        "source_locator": f"doi:{metadata['citation']['doi']} — {metadata['citation']['table_locator']}",
        "source_result_verified": True,
        "rows": len(rows),
        "partitions": {role: len(row_ids) for role, row_ids in role_ids.items()},
        "source_result_sha256": sha256(source_result_path),
        "target_sha256": sha256(target_path),
        "metadata_sha256": sha256(metadata_path),
        "official_source_sha256": source_hashes,
        "fit_target_sha256": metadata["fit_target_contract"]["sha256"],
        "source_verification_sha256": metadata["source_verification_contract"]["sha256"],
        "comparison_contract_sha256": metadata["comparison_contract"]["sha256"],
        "auxiliary_reference_eos": auxiliary,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate the frozen propane source contract.")
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--source-result", type=Path, default=DEFAULT_SOURCE_result)
    parser.add_argument("--auxiliary-csv", type=Path, default=DEFAULT_AUXILIARY_CSV)
    parser.add_argument("--auxiliary-metadata", type=Path, default=DEFAULT_AUXILIARY_METADATA)
    args = parser.parse_args()
    print(
        json.dumps(
            check(
                args.target,
                args.metadata,
                args.source_result,
                args.auxiliary_csv,
                args.auxiliary_metadata,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()

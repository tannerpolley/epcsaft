from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CSV = ROOT / "inputs/source" / "may-2015-methane-ethane-vle.csv"
DEFAULT_METADATA = ROOT / "inputs/source" / "may-2015-methane-ethane-vle.yaml"
EXPECTED_COLUMNS = [
    "row_id",
    "pathway",
    "T_K",
    "P_Pa",
    "x_methane",
    "x_ethane",
    "u_x_methane",
    "uc_x_methane",
    "y_methane",
    "y_ethane",
    "u_y_methane",
    "uc_y_methane",
    "temperature_standard_uncertainty_K",
    "pressure_standard_uncertainty_Pa",
    "x_comparison_allowance",
    "y_comparison_allowance",
]
EXPECTED_SOURCE_HASHES = {
    "nist_thermoml_json": "77630e90db70bb6aabfdfa520f61f14cee5076ece0265754140a25f771659662",
    "nist_thermoml_xml": "311e35b53e27bc050e17c4146a466e087c24c7a15624c29ababc4b7897d7871a",
    "publisher_article_pdf": "53fd1bdd55dc6807ec76cf88626438d8dfceb3ec09149d4405ea36cfbe6b842a",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_hash(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def check(csv_path: Path, metadata_path: Path) -> dict[str, object]:
    with csv_path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != EXPECTED_COLUMNS:
            raise ValueError(f"unexpected CSV columns: {reader.fieldnames}")
        rows = list(reader)
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    if metadata["citation"]["doi"] != "10.1021/acs.jced.5b00610":
        raise ValueError("unexpected source DOI")
    source_hashes = {name: source["sha256"] for name, source in metadata["sources"].items()}
    if source_hashes != EXPECTED_SOURCE_HASHES:
        raise ValueError("official source download hashes changed")
    csv_hash = sha256(csv_path)
    if metadata["retained_files"]["csv"]["sha256"] != csv_hash:
        raise ValueError("retained CSV hash does not match metadata")
    tolerance = metadata["tolerance_contract"]
    if tolerance["sha256"] != canonical_hash(tolerance["payload"]):
        raise ValueError("tolerance payload hash does not match metadata")

    if len(rows) != 17:
        raise ValueError(f"expected 17 audited rows, found {len(rows)}")
    expected_ids = [f"may2015-ch4-c2h6-{index:03d}" for index in range(1, 18)]
    if [row["row_id"] for row in rows] != expected_ids:
        raise ValueError("row IDs are not unique and ordered")
    isothermal_temperatures: list[float] = []
    multiplier = tolerance["payload"]["combined_standard_uncertainty_multiplier"]
    if multiplier != 3.0 or tolerance["payload"]["model_accuracy_floor"] is not None:
        raise ValueError("unexpected pre-model tolerance policy")
    for row in rows:
        values = {name: float(row[name]) for name in EXPECTED_COLUMNS[2:]}
        if not all(math.isfinite(value) for value in values.values()):
            raise ValueError(f"nonfinite value in {row['row_id']}")
        if values["T_K"] <= 0.0 or values["P_Pa"] <= 0.0:
            raise ValueError(f"nonpositive state value in {row['row_id']}")
        if not (0.0 < values["x_methane"] < 1.0 and 0.0 < values["y_methane"] < 1.0):
            raise ValueError(f"invalid phase composition in {row['row_id']}")
        if not math.isclose(values["x_methane"] + values["x_ethane"], 1.0, abs_tol=1e-12):
            raise ValueError(f"liquid composition does not normalize in {row['row_id']}")
        if not math.isclose(values["y_methane"] + values["y_ethane"], 1.0, abs_tol=1e-12):
            raise ValueError(f"vapor composition does not normalize in {row['row_id']}")
        uncertainty_fields = (
            "u_x_methane",
            "uc_x_methane",
            "u_y_methane",
            "uc_y_methane",
            "temperature_standard_uncertainty_K",
            "pressure_standard_uncertainty_Pa",
        )
        if any(values[name] < 0.0 for name in uncertainty_fields):
            raise ValueError(f"negative uncertainty in {row['row_id']}")
        if not math.isclose(
            values["x_comparison_allowance"],
            multiplier * values["uc_x_methane"],
            abs_tol=1e-12,
        ):
            raise ValueError(f"invalid liquid allowance in {row['row_id']}")
        if not math.isclose(
            values["y_comparison_allowance"],
            multiplier * values["uc_y_methane"],
            abs_tol=1e-12,
        ):
            raise ValueError(f"invalid vapor allowance in {row['row_id']}")
        if row["pathway"] == "isothermal":
            isothermal_temperatures.append(values["T_K"])
        elif row["pathway"] != "isochoric":
            raise ValueError(f"unexpected pathway in {row['row_id']}")
    if len(isothermal_temperatures) < 3 or max(isothermal_temperatures) - min(isothermal_temperatures) > 0.02:
        raise ValueError("no audited isotherm with at least three coexistence rows")

    return {
        "status": "accepted",
        "doi": metadata["citation"]["doi"],
        "table_locator": metadata["citation"]["table_locator"],
        "rows": len(rows),
        "csv_sha256": csv_hash,
        "metadata_sha256": sha256(metadata_path),
        "tolerance_sha256": tolerance["sha256"],
        "source_sha256": source_hashes,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate the frozen May 2015 source contract.")
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    args = parser.parse_args()
    print(json.dumps(check(args.csv, args.metadata), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

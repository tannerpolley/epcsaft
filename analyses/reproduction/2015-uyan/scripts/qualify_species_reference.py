"""Retain the historical pure-CO2-vapor declaration and its actual table differences."""
import csv
import hashlib
import json
from pathlib import Path

import epcsaft
from epcsaft.equilibrium import SpeciesReference
import generate_data as inputs
from compare_starts import ABSOLUTE, RELATIVE, TABLE

OUTPUT = inputs.CAMPAIGN / "results/data/native/species-reference"


def write_table(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    rows = inputs.read_rows(None)
    rows = [r for r in rows if 0.18 <= float(r["mdea_molality_mol_per_kg_water"]) * inputs.MOLAR_MASSES[1]
            / (1 + float(r["mdea_molality_mol_per_kg_water"]) * inputs.MOLAR_MASSES[1]) <= 0.20]
    records = [inputs.evaluate_row((r, 3e5 if float(r["temperature_k"]) >= 380 else 1e5),
                                   {"co2": SpeciesReference.PureGas}) for r in rows]
    for field, filename in (("_attempts", "attempts.csv"), ("_rows", "rows.csv"), ("_references", "references.csv")):
        diagnostics = [item for record in records for item in record.pop(field)]
        if diagnostics:
            write_table(OUTPUT / filename, diagnostics)
    write_table(OUTPUT / TABLE, records)
    historical = {r["source_row"]: r for r in csv.DictReader((inputs.CAMPAIGN / "results/data" / TABLE).open())}
    comparisons = []
    for record in records:
        if record["grepe_status"] != "evaluated":
            comparisons.append(dict(source_row=record["source_row"], quantity="availability", actual=0,
                                    historical=1, error=1, tolerance=0, met=False))
            continue
        for quantity in (*RELATIVE, *ABSOLUTE, "maximum_original_row"):
            actual = float(record[quantity])
            old = 0 if quantity == "maximum_original_row" else float(historical[record["source_row"]][quantity])
            error = abs(actual / old - 1) if quantity in RELATIVE else abs(actual - old)
            tolerance = 1e-10 if quantity == "maximum_original_row" else 1e-7
            comparisons.append(dict(source_row=record["source_row"], quantity=quantity, actual=actual,
                                    historical=old, error=error, tolerance=tolerance, met=error <= tolerance))
    write_table(OUTPUT / "comparisons.csv", comparisons)
    evaluated = [r for r in records if r["grepe_status"] == "evaluated"]
    wheel = next((inputs.ROOT / "build/environment-wheel").glob("epcsaft-*.whl"))
    summary = {
        "source_row_count": len(records), "evaluated_row_count": len(evaluated), "failed_row_count": len(records)-len(evaluated),
        "maximum_original_row": max((float(r["maximum_original_row"]) for r in evaluated), default=None),
        "comparisons": len(comparisons), "misses": sum(not c["met"] for c in comparisons),
        "misses_by_quantity": {q: sum(not c["met"] for c in comparisons if c["quantity"] == q)
                               for q in (*RELATIVE, *ABSOLUTE, "maximum_original_row", "availability")},
        "maximum_relative_error": max((c["error"] for c in comparisons if c["quantity"] in RELATIVE), default=None),
        "maximum_composition_error": max((c["error"] for c in comparisons if c["quantity"] in ABSOLUTE), default=None),
        "liquid_denser_than_vapor_count": sum(float(r["liquid_density_mol_m3"]) > float(r["vapor_density_mol_m3"]) for r in evaluated),
        "input_sha256": {str(path.relative_to(inputs.ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in (inputs.ROOT / "analyses/reproduction/2015-uyan/scripts/generate_data.py",
                                      inputs.CAMPAIGN / "inputs/source/ermatchkov-2006-tables-1-3.csv",
                                      inputs.CAMPAIGN / "inputs/source/uyan-2015-source-ledger.json")},
        "wheel_sha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
        "installed_engine": epcsaft.__file__, "base_commit": "92382dae",
        "historical_commit": "566d57de5ec6b2bd36d1e1f0460b62a837e70e0e", "historical_wheel_sha256": None,
        "declaration": {"solvent": "water", "species_reference": {"co2": "PureGas"},
                        "basis": "RawMoleFractionInfiniteDilution", "source_pressure_Pa": 1e5},
        "meaning": "Numerical diagnostic of the e12a5ff8 generator declaration against an older table with different convention labels; no historical-equivalence or physical-validation claim.",
    }
    (OUTPUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

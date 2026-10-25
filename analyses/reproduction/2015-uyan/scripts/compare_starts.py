"""Compare native declared-start output with the retained 30-row execution record."""
import csv
import json
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "results/data"
NATIVE = DATA / "native"
TABLE = "uyan-2015-ermatchkov-2006-grepe-plotted-data.csv"
RELATIVE = ("grepe_total_pressure_kpa", "grepe_co2_partial_pressure_kpa",
            "liquid_density_mol_m3", "vapor_density_mol_m3")
ABSOLUTE = ("vapor_water_mole_fraction", "vapor_mdea_mole_fraction", "vapor_co2_mole_fraction")


def main():
    old = {row["source_row"]: row for row in csv.DictReader((DATA / TABLE).open())}
    new = {row["source_row"]: row for row in csv.DictReader((NATIVE / TABLE).open())}
    comparisons = []
    for identity, reference in old.items():
        current = new.get(identity)
        if current is None or current["grepe_status"] != "evaluated":
            comparisons.append(dict(source_row=identity, quantity="availability", old=1, new=0,
                                    error=1, tolerance=0, met=False))
            continue
        for quantity in (*RELATIVE, *ABSOLUTE, "selected_pressure_seed_kpa"):
            a, b = float(reference[quantity]), float(current[quantity])
            error = abs(b / a - 1) if quantity in RELATIVE else abs(b - a)
            tolerance = 0. if quantity == "selected_pressure_seed_kpa" else 1e-7
            comparisons.append(dict(source_row=identity, quantity=quantity, old=a, new=b,
                                    error=error, tolerance=tolerance, met=error <= tolerance))
        raw = float(current["maximum_original_row"])
        comparisons.append(dict(source_row=identity, quantity="maximum_original_row", old="", new=raw,
                                error=raw, tolerance=1e-10, met=raw <= 1e-10))
    with (NATIVE / "comparisons.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparisons[0]))
        writer.writeheader()
        writer.writerows(comparisons)
    summary = {"historical_commit": "566d57de5ec6b2bd36d1e1f0460b62a837e70e0e",
               "historical_wheel_sha256": None, "historical_rows": len(old), "native_rows": len(new),
               "missing_rows": sorted(set(old) - set(new)), "unexpected_rows": sorted(set(new) - set(old)),
               "comparisons": len(comparisons), "misses": sum(not r["met"] for r in comparisons),
               "maximum_relative_error": max((r["error"] for r in comparisons if r["quantity"] in RELATIVE), default=None),
               "maximum_composition_error": max((r["error"] for r in comparisons if r["quantity"] in ABSOLUTE), default=None),
               "meaning": "Numerical comparison of the declared historical reconstruction, not physical validation or a global-stability test."}
    (NATIVE / "comparison-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

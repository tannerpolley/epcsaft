"""Compare retained reference-convention lanes without evaluating the Engine."""

import csv
import json
from pathlib import Path

CAMPAIGN = Path(__file__).resolve().parents[1]
OLD = CAMPAIGN / "results/data"
NEW = OLD / "native/source-ledger"
QUANTITIES = {
    "grepe_total_pressure_kpa": "kPa",
    "grepe_co2_partial_pressure_kpa": "kPa",
    "grepe_co2_partial_pressure_relative_error": "1",
    "liquid_density_mol_m3": "mol/m3",
    "vapor_density_mol_m3": "mol/m3",
    "vapor_water_mole_fraction": "1",
    "vapor_mdea_mole_fraction": "1",
    "vapor_co2_mole_fraction": "1",
}


def main():
    name = "uyan-2015-ermatchkov-2006-grepe-plotted-data.csv"
    old = {r["source_row"]: r for r in csv.DictReader((OLD / name).open())}
    new = {r["source_row"]: r for r in csv.DictReader((NEW / name).open())}
    assert old.keys() == new.keys()
    comparisons = []
    for key, before in old.items():
        after = new[key]
        for field in (
            "temperature_k",
            "mdea_molality_mol_per_kg_water",
            "co2_molality_mol_per_kg_water",
        ):
            assert before[field] == after[field], (key, field)
        for quantity, unit in QUANTITIES.items():
            a = float(before[quantity]) if before.get(quantity) else None
            b = float(after[quantity]) if after.get(quantity) else None
            tolerance = None if quantity.endswith("relative_error") else 1e-7
            relative = quantity.endswith("pressure_kpa") or "density" in quantity
            error = (
                (abs(b / a - 1) if relative else abs(b - a))
                if a is not None and b is not None
                else None
            )
            comparisons.append(
                {
                    "source_row": key,
                    "temperature_k": before["temperature_k"],
                    "mdea_molality_mol_per_kg_water": before[
                        "mdea_molality_mol_per_kg_water"
                    ],
                    "co2_molality_mol_per_kg_water": before[
                        "co2_molality_mol_per_kg_water"
                    ],
                    "quantity": quantity,
                    "unit": unit,
                    "historical": a,
                    "source_ledger": b,
                    "change": b - a if a is not None and b is not None else None,
                    "historical_equivalence_error": error
                    if tolerance is not None
                    else None,
                    "historical_equivalence_tolerance": tolerance,
                    "historical_equivalence_met": (
                        error <= tolerance if error is not None else False
                    )
                    if tolerance is not None
                    else None,
                    "claim_change": "Difference between historical reconstruction and current source-ledger calculation",
                    "native_cause": after.get("grepe_failure", ""),
                }
            )
    with (NEW / "reference-convention-comparison.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(comparisons[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(comparisons)
    summary_name = "uyan-2015-ermatchkov-2006-grepe-summary.json"
    a, b = (json.loads((folder / summary_name).read_text()) for folder in (OLD, NEW))
    summary = {
        "historical_record_commit": "566d57de5ec6b2bd36d1e1f0460b62a837e70e0e",
        "historical_wheel_sha256": None,
        "wheel_sha256": b["wheel_sha256"],
        "conditions": "313.15, 353.15 and 393.15 K; CO2-free solvent MDEA mass fraction 0.18-0.20; 1 kg water basis",
        "reference_pressure_pa": 1e5,
        "row_tolerance": 1e-10,
        "source_reference_convention": b["source_reference_convention"],
        "fixed_start_anchor": b["fixed_start_anchor"],
        "fixed_start_anchor_rows": b["fixed_start_anchor_rows"],
        "pressure_only_replay": b["pressure_only_replay"],
        "coverage_interpretation": "Four additional roots use the owner-approved fixed t1-r08 composition and packing guesses at the original ordered pressure starts. Coverage recovery is an initialization change; equations, feeds, work limits and tolerances are unchanged.",
        "historical_generator_default": "pure-water liquid / pure-CO2 vapor / other solutes infinite dilution",
        "historical_table_lane_labels": sorted(
            {
                (r["reaction_constant_lane"], r["source_reference_lane"])
                for r in old.values()
            }
        ),
        "interpretation": "Difference between historical reconstruction and current source-ledger calculation. Historical wheel identity is unavailable and historical lane labels conflict, so reference-convention effects cannot be isolated from historical Engine differences.",
        "reported_quantity_comparisons": len(comparisons),
        "historical_equivalence_comparisons": sum(
            r["historical_equivalence_met"] is not None for r in comparisons
        ),
        "historical_equivalence_misses": sum(
            r["historical_equivalence_met"] is False for r in comparisons
        ),
        "changed_values": sum(
            r["change"] is not None and r["change"] != 0 for r in comparisons
        ),
        "unavailable_values": sum(r["source_ledger"] is None for r in comparisons),
        "historical_roots": a["evaluated_row_count"],
        "source_ledger_roots": b["evaluated_row_count"],
        "failed_rows": [
            key for key, r in new.items() if r["grepe_status"] != "evaluated"
        ],
        "historical_ard_percent": a["grepe_ard_percent"],
        "source_ledger_ard_percent": b["grepe_ard_percent"],
        "ard_comparison_basis": "The same 30 source rows when both lanes return 30 roots; otherwise each lane's evaluated rows with unequal coverage",
        "historical_ard_by_temperature_percent": a["grepe_ard_by_temperature_percent"],
        "source_ledger_ard_by_temperature_percent": b[
            "grepe_ard_by_temperature_percent"
        ],
        "maximum_original_row": b["maximum_original_row"],
        "limit": "Numerical verification, not physical validation",
    }
    (NEW / "reference-convention-summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

"""Read retained results to present reference identities and the historical diagnostic."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/data"
NEW = DATA / "native/species-reference"
OUTPUT = ROOT / "results/rendered/species-reference"
TABLE = "uyan-2015-ermatchkov-2006-grepe-plotted-data.csv"


def read(path: Path, delimiter: str = ",") -> list[dict]:
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def table(headers: list[str], rows: list[list]) -> str:
    return "\n".join(["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |",
                      *("| " + " | ".join(map(str, r)) + " |" for r in rows)]) + "\n"


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    old = json.loads((DATA / "uyan-2015-ermatchkov-2006-grepe-summary.json").read_text())
    source = json.loads((DATA / "native/source-ledger/uyan-2015-ermatchkov-2006-grepe-summary.json").read_text())
    new = json.loads((NEW / "summary.json").read_text())
    checks = json.loads((NEW / "qualification-summary.json").read_text())
    comparisons = read(NEW / "comparisons.csv")
    historical_tolerance = next(float(r["tolerance"]) for r in comparisons if r["quantity"] != "maximum_original_row")
    conditions = json.loads((DATA / "native/source-ledger/reference-convention-summary.json").read_text())
    (OUTPUT / "conditions.md").write_text(
        conditions["conditions"] + f". Source pressure is {conditions['reference_pressure_pa']/1e3:g} kPa. "
        f"The original-row limit is {conditions['row_tolerance']:g}; historical pressure/density and composition comparisons retain {historical_tolerance:g}.\n")
    (OUTPUT / "source-ledger.md").write_text(table(
        ["Retained quantity", "Historical table", "Source-ledger calculation"],
        [["Accepted roots", old["evaluated_row_count"], source["evaluated_row_count"]],
         ["CO₂ partial-pressure average relative deviation (%)", f"{old['grepe_ard_percent']:.6f}", f"{source['grepe_ard_percent']:.6f}"],
         *[[f"Same deviation at {t} K (%)", f"{old['grepe_ard_by_temperature_percent'][t]:.6f}",
             f"{source['grepe_ard_by_temperature_percent'][t]:.6f}"] for t in source["grepe_ard_by_temperature_percent"]],
         ["Maximum original-row residual", "Unidentified wheel", f"{source['maximum_original_row']:.9g}"]]))
    (OUTPUT / "diagnostic.md").write_text(table(
        ["Pure-CO₂-vapor diagnostic quantity", "Retained value"],
        [["Source rows / accepted roots / failures", f"{new['source_row_count']} / {new['evaluated_row_count']} / {new['failed_row_count']}"],
         ["Maximum original-row residual", f"{new['maximum_original_row']:.9g}"],
         ["Historical comparisons / misses", f"{new['comparisons']} / {new['misses']}"],
         ["Maximum pressure or density relative difference", f"{new['maximum_relative_error']:.9g}"],
         ["Maximum vapor mole-fraction absolute difference", f"{new['maximum_composition_error']:.9g}"],
         *[[q+" misses", n] for q, n in new["misses_by_quantity"].items()]]))
    (OUTPUT / "qualification.md").write_text(table(
        ["Numerical verification", "Retained evidence"], [[k,v] for k,v in checks["reported_quantities"].items()]))
    probe = json.loads((DATA / "native/source-ledger/failed-start-probe.json").read_text())
    (OUTPUT / "probe.md").write_text(table(
        ["Source row", "Probe total pressure (kPa)", "Original-row maximum", "Iterations"],
        [[r["source_row"], f"{r['pressure_pa']/1e3:.9f}", f"{r['maximum_original_row']:.9g}", r["iterations"]] for r in probe["attempts"]]))
    series = [("Historical table", read(DATA / TABLE), "o", "#0072B2"),
              ("Source-ledger reference", read(DATA / "native/source-ledger" / TABLE), "x", "#D55E00"),
              ("Pure-CO₂-vapor reference", read(NEW / TABLE), "+", "#009E73")]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), layout="constrained")
    temperatures = sorted({float(r["temperature_k"]) for r in series[0][1]})
    for ax, t in zip(axes, temperatures, strict=True):
        for label, rows, marker, color in series:
            rows = [r for r in rows if float(r["temperature_k"]) == t and r["grepe_status"] == "evaluated"]
            ax.scatter([float(r["co2_molality_mol_per_kg_water"]) for r in rows],
                       [float(r["grepe_co2_partial_pressure_kpa"]) for r in rows],
                       label=label, marker=marker, color=color, s=30)
        ax.set(title=f"{t:g} K", xlabel="CO₂ feed (mol/kg water)", yscale="log")
        ax.grid(alpha=.2)
    axes[0].set_ylabel("Predicted CO₂ partial pressure (kPa)")
    axes[0].legend(fontsize=8)
    fig.savefig(OUTPUT / "reference-comparison.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()

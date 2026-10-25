from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt

plt.switch_backend("Agg")

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "results/co2-partial-pressure-validation.csv"
OUTPUT = ROOT / "results/co2-partial-pressure-validation"


def main() -> None:
    with TABLE.open(newline="") as handle:
        rows = [
            row
            for row in csv.DictReader(handle)
            if row["prediction_status"] == "evaluated"
        ]
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.4), sharey=True)
    colors = {"Hilliard2008": "#0072B2", "Jou1995": "#D55E00"}
    names = {"Hilliard2008": "Hilliard", "Jou1995": "Jou"}
    for axis, temperature in zip(axes, (40.0, 60.0, 80.0), strict=True):
        selected = [row for row in rows if float(row["temperature_c"]) == temperature]
        for source in sorted({row["source_key"] for row in selected}):
            group = sorted(
                (row for row in selected if row["source_key"] == source),
                key=lambda row: float(row["co2_loading_mol_per_mol_mea"]),
            )
            x = [float(row["co2_loading_mol_per_mol_mea"]) for row in group]
            axis.scatter(
                x,
                [float(row["observed_pco2_kpa"]) for row in group],
                facecolors="none",
                edgecolors=colors[source],
                label=f"{names[source]} observations",
            )
            axis.plot(
                x,
                [float(row["predicted_pco2_kpa"]) for row in group],
                color=colors[source],
                label=f"Engine prediction ({names[source]} states)",
            )
        axis.set(yscale="log", title=f"{temperature:.0f} °C")
        axis.set_xlabel(r"$\alpha_{\mathrm{CO_2}}$ (mol/mol MEA)")
        axis.grid(True, which="both", alpha=0.2)
    axes[0].set_ylabel(r"$p_{\mathrm{CO_2}}$ (kPa)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncols=4, frameon=False)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.87))
    for suffix in ("pdf", "png", "svg"):
        fig.savefig(OUTPUT.with_suffix(f".{suffix}"), dpi=220)
    plt.close(fig)
    svg = OUTPUT.with_suffix(".svg")
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")


if __name__ == "__main__":
    main()

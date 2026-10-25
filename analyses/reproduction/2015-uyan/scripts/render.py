"""Display retained reference-convention predictions at source observation rows."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt

CAMPAIGN = Path(__file__).resolve().parents[1]
HISTORICAL = CAMPAIGN / "results/data/uyan-2015-ermatchkov-2006-grepe-plotted-data.csv"
CURRENT = (
    CAMPAIGN
    / "results/data/native/source-ledger/uyan-2015-ermatchkov-2006-grepe-plotted-data.csv"
)
OUTPUT = CAMPAIGN / "results/rendered/source-ledger"


def main() -> None:
    old, new = (list(csv.DictReader(path.open())) for path in (HISTORICAL, CURRENT))
    temperatures = sorted({float(row["temperature_k"]) for row in old})
    plt.rcParams.update(
        {"font.family": "serif", "axes.spines.top": False, "axes.spines.right": False}
    )
    figure, axes = plt.subplots(1, len(temperatures), figsize=(12.6, 4.4), sharey=True)
    for axis, temperature in zip(axes, temperatures, strict=True):
        observed = [r for r in old if float(r["temperature_k"]) == temperature]
        axis.errorbar(
            [float(r["co2_loading_mol_per_mol_mdea"]) for r in observed],
            [float(r["co2_partial_pressure_kpa"]) for r in observed],
            yerr=[float(r["co2_partial_pressure_uncertainty_kpa"]) for r in observed],
            fmt="o",
            color="#111111",
            markerfacecolor="white",
            ms=5,
            label="Ermatchkov observations",
        )
        for values, color, marker, label in (
            (old, "#D55E00", "x", "Historical reconstruction"),
            (new, "#0072B2", "s", "Current source-ledger calculation"),
        ):
            selected = [
                r
                for r in values
                if float(r["temperature_k"]) == temperature
                and r["grepe_status"] == "evaluated"
            ]
            axis.scatter(
                [float(r["co2_loading_mol_per_mol_mdea"]) for r in selected],
                [float(r["grepe_co2_partial_pressure_kpa"]) for r in selected],
                color=color,
                marker=marker,
                s=25,
                label=label,
            )
        axis.set_yscale("log")
        axis.set_title(f"{temperature:.2f} K")
        axis.set_xlabel(r"CO$_2$ loading / mol mol$_{\mathrm{MDEA}}^{-1}$")
        axis.grid(True, which="major", color="#D8D8D8", linewidth=0.7)
    axes[0].set_ylabel(r"CO$_2$ partial pressure / kPa")
    handles, labels = axes[-1].get_legend_handles_labels()
    figure.legend(handles, labels, loc="outside lower center", ncol=3, frameon=False)
    figure.suptitle(
        "Difference between historical reconstruction and current source-ledger calculation\n"
        "Finite reacting liquid + incipient neutral vapor; discrete source rows"
    )
    figure.subplots_adjust(bottom=0.23, top=0.80, left=0.07, right=0.99, wspace=0.10)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT / "co2-partial-pressure.png", dpi=180, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()

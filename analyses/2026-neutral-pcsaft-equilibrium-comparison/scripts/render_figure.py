#!/usr/bin/env python3
"""Render the retained comparison table without rerunning any model."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt

ANALYSIS = Path(__file__).resolve().parents[1]
TABLE = ANALYSIS / "results/values.csv"
FIGURES = ANALYSIS / "figures"


def main() -> None:
    with TABLE.open() as handle:
        rows = list(csv.DictReader(handle))
    names = [row["implementation"] for row in rows]
    colors = ["#1b6ca8", "#d95f02", "#7570b3", "#1b9e77", "#e7298a"]
    quantities = [
        ("pressure_Pa", "Saturation pressure [kPa]", 1.0e-3, "{:.6f}"),
        ("rho_liquid_mol_m3", "Liquid molar density [kmol m$^{-3}$]", 1.0e-3, "{:.6f}"),
        ("rho_vapor_mol_m3", "Vapor molar density [mol m$^{-3}$]", 1.0, "{:.6f}"),
    ]

    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.4), constrained_layout=True)
    for axis, (field, label, scale, value_format) in zip(axes, quantities):
        values = [float(row[field]) * scale for row in rows]
        axis.scatter(range(len(rows)), values, c=colors, s=75, zorder=3)
        for position, value in enumerate(values):
            axis.annotate(
                value_format.format(value),
                (position, value),
                xytext=(0, 9),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7,
                rotation=70,
            )
        axis.set_xticks(range(len(rows)), names, rotation=35, ha="right")
        axis.set_ylabel(label)
        axis.ticklabel_format(axis="y", style="plain", useOffset=False)
        axis.grid(axis="y", color="#d9d9d9", linewidth=0.7)
        axis.margins(x=0.12)
    fig.suptitle("Neutral propane PC-SAFT saturation at 300 K", fontsize=14)
    fig.text(
        0.5,
        0.01,
        "Same Gross–Sadowski 2001 parameters; points are independent implementations, not an experimental oracle",
        ha="center",
        fontsize=9,
    )
    fig.savefig(FIGURES / "neutral-propane-pcsaft-comparison.png", dpi=220)
    fig.savefig(FIGURES / "neutral-propane-pcsaft-comparison.svg")
    print(FIGURES / "neutral-propane-pcsaft-comparison.png")


if __name__ == "__main__":
    main()

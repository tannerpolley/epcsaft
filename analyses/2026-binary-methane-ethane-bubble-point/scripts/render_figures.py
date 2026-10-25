#!/usr/bin/env python3
"""Render binary result and performance figures from the retained CSV."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

ANALYSIS = Path(__file__).resolve().parents[1]
TABLE = ANALYSIS / "results/values.csv"
FIGURES = ANALYSIS / "figures"


def main() -> None:
    with TABLE.open() as handle:
        rows = list(csv.DictReader(handle))
    names = [row["implementation"] for row in rows]
    colors = ["#1b6ca8", "#d95f02", "#1b9e77"]

    result_fields = [
        ("pressure_Pa", "Bubble pressure [MPa]", 1.0e-6, "{:.9f}"),
        ("vapor_y_methane", "Incipient vapor $y_{CH4}$", 1.0, "{:.9f}"),
        ("rho_liquid_mol_m3", "Liquid density [kmol m$^{-3}$]", 1.0e-3, "{:.6f}"),
        ("rho_vapor_mol_m3", "Vapor density [kmol m$^{-3}$]", 1.0e-3, "{:.6f}"),
    ]
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.0), constrained_layout=True)
    for axis, (field, label, scale, fmt) in zip(axes.flat, result_fields):
        values = [float(row[field]) * scale for row in rows]
        axis.scatter(range(len(rows)), values, c=colors, s=80, zorder=3)
        for position, value in enumerate(values):
            axis.annotate(fmt.format(value), (position, value), xytext=(0, 8), textcoords="offset points", ha="center", va="bottom", fontsize=8, rotation=45)
        axis.set_xticks(range(len(rows)), names, rotation=25, ha="right")
        axis.set_ylabel(label)
        axis.yaxis.set_major_formatter(FormatStrFormatter("%.6f"))
        axis.grid(axis="y", color="#d9d9d9", linewidth=0.7)
        axis.margins(x=0.15)
    fig.suptitle("Methane/ethane PC-SAFT bubble point at 180 K", fontsize=14)
    fig.savefig(FIGURES / "binary-methane-ethane-bubble-point.png", dpi=220)
    fig.savefig(FIGURES / "binary-methane-ethane-bubble-point.svg")

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), constrained_layout=True)
    positions = range(len(rows))
    for axis, median_field, label in zip(axes, ("solve_median_ms", "setup_ms"), ("Warm bubble solve [ms]", "Model setup [ms]")):
        values = [float(row[median_field]) for row in rows]
        if median_field == "solve_median_ms":
            lows = [float(row["solve_p10_ms"]) for row in rows]
            highs = [float(row["solve_p90_ms"]) for row in rows]
            axis.errorbar(positions, values, yerr=[[v - lo for v, lo in zip(values, lows)], [hi - v for v, hi in zip(values, highs)]], fmt="o", color="#1b6ca8", ecolor="#777777", capsize=4)
        else:
            axis.scatter(positions, values, c=colors, s=80, zorder=3)
        axis.set_xticks(list(positions), names, rotation=25, ha="right")
        axis.set_ylabel(label)
        axis.set_yscale("log")
        axis.grid(axis="y", color="#d9d9d9", linewidth=0.7)
        axis.margins(x=0.15)
    fig.suptitle("Runtime statistics: warm solves and setup", fontsize=14)
    fig.text(0.5, 0.01, "5 warm-ups, 25 timed calls; p10–p90 whiskers on solve panel; imports excluded", ha="center", fontsize=9)
    fig.savefig(FIGURES / "binary-methane-ethane-runtime.png", dpi=220)
    fig.savefig(FIGURES / "binary-methane-ethane-runtime.svg")
    print(FIGURES)


if __name__ == "__main__":
    main()

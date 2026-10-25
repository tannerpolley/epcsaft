"""Render TP-flash value and runtime figures from the retained CSV."""

from __future__ import annotations

import csv
import math
import statistics
from pathlib import Path

import matplotlib.pyplot as plt

ANALYSIS = Path(__file__).resolve().parents[1]
TABLE = ANALYSIS / "results/values.csv"
FIGURES = ANALYSIS / "figures"


def finite(row: dict[str, str], field: str) -> float | None:
    value = row.get(field, "").strip()
    if not value:
        return None
    parsed = float(value)
    return parsed if math.isfinite(parsed) else None


def main() -> None:
    with TABLE.open() as handle:
        rows = list(csv.DictReader(handle))
    names = [row["implementation"] for row in rows]
    colors = ["#1b6ca8", "#9467bd", "#d95f02", "#7570b3", "#1b9e77"]
    FIGURES.mkdir(parents=True, exist_ok=True)

    fields = [
        ("beta_vapor", "Vapor fraction $\\beta_V$", 1.0, "{:.9f}"),
        ("x_methane", "Liquid $x_{CH4}$", 1.0, "{:.9f}"),
        ("y_methane", "Vapor $y_{CH4}$", 1.0, "{:.9f}"),
        ("rho_liquid_mol_m3", "Liquid density [kmol m$^{-3}$]", 1.0e-3, "{:.6f}"),
        ("rho_vapor_mol_m3", "Vapor density [kmol m$^{-3}$]", 1.0e-3, "{:.6f}"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(14.0, 8.0))
    for axis, (field, label, scale, fmt) in zip(axes.flat, fields):
        values = [finite(row, field) for row in rows]
        for i, value in enumerate(values):
            if value is None:
                continue
            scaled = value * scale
            axis.scatter(i, scaled, c=colors[i], s=70, zorder=3)
        axis.set_xticks(range(len(rows)), names, rotation=25, ha="right")
        axis.set_ylabel(label)
        axis.grid(axis="y", color="#d9d9d9", linewidth=0.7)
        axis.margins(x=0.15)
    axes.flat[-1].axis("off")
    fig.suptitle("Neutral methane/ethane PC-SAFT TP flash at 180 K and 1 MPa", fontsize=14)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.90, bottom=0.20, wspace=0.35, hspace=0.45)
    fig.savefig(FIGURES / "binary-methane-ethane-tp-flash-values.png", dpi=220)
    fig.savefig(FIGURES / "binary-methane-ethane-tp-flash-values.svg")

    fig, axis = plt.subplots(figsize=(10.5, 4.8))
    measured = [finite(row, "median_ms") for row in rows]
    positions = list(range(len(rows)))
    plotted = [value if value is not None else 0.1 for value in measured]
    bars = axis.bar(positions, plotted, color=colors, alpha=0.9)
    for bar, value in zip(bars, measured):
        if value is None:
            bar.set_alpha(0.2)
            axis.text(bar.get_x() + bar.get_width() / 2, 0.14, "not measured\\n(resource cap)", ha="center", va="bottom", fontsize=8, rotation=90)
        else:
            axis.text(bar.get_x() + bar.get_width() / 2, value * 1.15, f"{value:.3f}", ha="center", va="bottom", fontsize=8)
    axis.set_xticks(positions, names, rotation=25, ha="right")
    axis.set_ylabel("Warm TP-flash median [ms]")
    axis.set_yscale("log")
    axis.grid(axis="y", color="#d9d9d9", linewidth=0.7)
    axis.set_title("Runtime: 3 warm-ups, 15 timed calls; imports excluded")
    fig.subplots_adjust(left=0.10, right=0.98, top=0.88, bottom=0.24)
    fig.savefig(FIGURES / "binary-methane-ethane-tp-flash-runtime.png", dpi=220)
    fig.savefig(FIGURES / "binary-methane-ethane-tp-flash-runtime.svg")
    print(FIGURES)


def performance() -> None:
    """Discrete case medians and observed p10-p90 timing spread, not a fitted curve."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True, layout="constrained")
    labels = ["180 K / 1 MPa", "179 K", "181 K", "0.8 MPa", "1.2 MPa",
              "$z_{CH4}=0.4$", "$z_{CH4}=0.6$", "10 mol feed", "Reversed phases"]
    for version, color, marker, offset in (("before", "#D55E00", "x", -.12),
                                           ("after", "#0072B2", "o", .12)):
        with (ANALYSIS / f"results/performance-{version}.csv").open() as handle:
            samples = list(csv.DictReader(handle))
        for axis, mode in zip(axes, ("cold", "previous")):
            for case in range(len(labels)):
                values = [float(row["wall_ms"]) for row in samples
                          if int(row["case"]) == case and row["mode"] == mode]
                median = statistics.median(values)
                deciles = statistics.quantiles(values, n=10, method="inclusive")
                axis.errorbar(median, case + offset,
                    xerr=[[median-deciles[0]], [deciles[-1]-median]], fmt=marker,
                    color=color, capsize=3, label=version.capitalize() if case == 0 else None)
            axis.set_xlabel("Solve wall time [ms]")
            axis.grid(axis="x", alpha=.25)
            axis.set_xlim(left=0)
            axis.legend()
    axes[0].set_yticks(range(len(labels)), labels)
    axes[0].invert_yaxis()
    axes[0].set_title("No phase guesses; automatic initialization included")
    axes[1].set_title("Previous-result seed (same or nearby state)")
    fig.suptitle("Neutral methane/ethane TP flash: same general exact-Hessian driver\n"
                 "15 samples; warm process; points = median, intervals = p10–p90")
    fig.savefig(FIGURES / "performance-sprint.png", dpi=180)


if __name__ == "__main__":
    import sys
    performance() if "--performance" in sys.argv else main()

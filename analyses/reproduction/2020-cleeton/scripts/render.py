"""Render the native Figure 16 replay from retained CSVs (no Engine evaluation)."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt

CAMPAIGN = Path(__file__).resolve().parents[1]
DIGITIZED = CAMPAIGN / "inputs" / "digitized" / "cleeton-2020-figure-16-curves.csv"
GATE = CAMPAIGN / "results" / "data" / "cleeton-2020-figure-16-gate.csv"
DICKO = CAMPAIGN / "inputs" / "source" / "dicko-2010-table-1.csv"
OUT = CAMPAIGN / "results" / "rendered"
COLORS = {0.1: "#0072B2", 0.3: "#009E73", 0.5: "#D55E00", 0.8: "#CC79A7"}
MARKERS = {0.1: "o", 0.3: "s", 0.5: "^", 0.8: "D"}
LANES = {"main-0.88sigma": ("x", "native, d = 0.88σ (main)"), "diagnostic-sigma": ("+", "native, d = σ (diagnostic)")}


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    curves, gate, dicko = read(DIGITIZED), read(GATE), read(DICKO)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    for ax, (panel, gas, other, column) in zip(axes, (
            ("left", "co2", "H₂S", "figure_16_left_series_h2s_loading"),
            ("right", "h2s", "CO₂", "figure_16_right_series_co2_loading"))):
        for series in COLORS:
            pts = [(float(r["loading"]), float(r["pressure_pa"])) for r in curves
                   if r["panel"] == panel and float(r["series_loading"]) == series]
            ax.plot(*zip(*pts), color=COLORS[series], lw=1.2, label=f"Cleeton curve, {other} loading {series}")
            data = [(float(r[f"{gas}_loading"]), 1e6 * float(r[f"{gas}_pressure_mpa"])) for r in dicko
                    if r[column] and float(r[column]) == series]
            ax.plot(*zip(*data), MARKERS[series], mfc="none", color=COLORS[series], ms=7)
            for lane, (marker, _) in LANES.items():
                native = [(float(g["loading"]), float(g[f"p_{gas}_pa"])) for g in gate
                          if g["lane"] == lane and g["reading"] == "caption" and g["panel"] == panel
                          and float(g["series_loading"]) == series
                          and g["success"] == "True"]
                if native:
                    ax.plot(*zip(*native), marker, color=COLORS[series], ms=8, mew=1.6)
        ax.set_yscale("log")
        ax.set_xlabel(f"{gas.upper().replace('2', '₂')} loading / mol mol⁻¹ MDEA")
        ax.set_ylabel(f"p({gas.upper().replace('2', '₂')}) / Pa")
        ax.grid(alpha=0.25, which="both")
    handles = [plt.Line2D([], [], color="0.3", lw=1.2, label="Cleeton 2020 model curve (digitized)"),
               plt.Line2D([], [], color="0.3", marker="o", mfc="none", ls="", label="Dicko 2010 Table 1"),
               *[plt.Line2D([], [], color="0.3", marker=m, ls="", mew=1.6, label=t) for m, t in LANES.values()]]
    axes[0].legend(handles=handles, fontsize=8, loc="upper left")
    fig.suptitle("Cleeton 2020 Figure 16 replay, 323.15 K, 50 wt% MDEA (colour = second-gas loading 0.1/0.3/0.5/0.8)")
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"cleeton-2020-figure-16-native.{suffix}", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()

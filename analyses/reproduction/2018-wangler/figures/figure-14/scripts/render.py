from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

OUTPUT = Path(__file__).resolve().parents[1] / "output"
DATA = OUTPUT / "wangler-2018-table-10-native.csv"
LANES = {
    "printed-S0": ("#555555", "D", r"Printed Table 6 $k_{ij}$, S0"),
    "fitted-S0": ("#1261A0", "s", r"Figure 10 $k_{ij}$, S0 (self-dispersion excluded)"),
    "fitted-S1": ("#B24C14", "^", r"Figure 10 $k_{ij}$, S1 (self-dispersion kept)"),
}


def main() -> None:
    data = pd.read_csv(DATA).sort_values("co2_loading_mol_per_mol_mdea")
    plt.rcParams.update({"font.family": "serif", "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(8.6, 4.4), constrained_layout=True)
    measured = data[data["lane"] == "fitted-S0"]
    ax.errorbar(measured["co2_loading_mol_per_mol_mdea"], measured["experimental_total_pressure_bar"],
                yerr=measured["pressure_standard_uncertainty_bar"], color="#111111", marker="o",
                markerfacecolor="white", linestyle="none", markeredgewidth=1.4, capsize=2,
                label="Wangler Table 10", zorder=3)
    for lane, (color, marker, label) in LANES.items():
        group = data[data["lane"] == lane]
        ax.plot(group["co2_loading_mol_per_mol_mdea"], group["native_pressure_bar"], color=color,
                marker=marker, linewidth=1.8, label=label)
        if lane.startswith("fitted"):
            low = data[data["lane"] == f"{lane}-minus"]["native_pressure_bar"].to_numpy()
            high = data[data["lane"] == f"{lane}-plus"]["native_pressure_bar"].to_numpy()
            ax.fill_between(group["co2_loading_mol_per_mol_mdea"], low, high, color=color, alpha=0.25, linewidth=0)
    ax.set_yscale("log")
    ax.set_title("Wangler 2018 Table 10: 60 wt% MDEA + CO$_2$ at 303 K (predicted)")
    ax.set_xlabel(r"CO$_2$ loading / mol mol$_{\mathrm{MDEA}}^{-1}$")
    ax.set_ylabel("Total pressure / bar")
    ax.grid(True, which="both", color="#D8D8D8", linewidth=0.7)
    ax.legend(frameon=False, fontsize=8.0, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    for suffix, kwargs in (("svg", {}), ("png", {"dpi": 200}), ("pdf", {})):
        fig.savefig(OUTPUT / f"wangler-2018-table-10-native.{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)


if __name__ == "__main__":
    main()

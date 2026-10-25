from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

OUTPUT = Path(__file__).resolve().parents[1] / "output"
DATA = OUTPUT / "wangler-2018-mdeah-osmotic-plotted-data.csv"
TITLES = {"S0": "S0: ion self-dispersion excluded", "S1": "S1: ion self-dispersion kept (l.146)"}


def main() -> None:
    data = pd.read_csv(DATA)
    plt.rcParams.update({"font.family": "serif", "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.2), constrained_layout=True, sharex=True, sharey=True)
    colors = {273.0: "#1261A0", 313.0: "#B33A3A"}
    markers = {273.0: "o", 313.0: "s"}
    for ax, (lane, lane_data) in zip(axes, data.groupby("lane", sort=True), strict=True):
        for temperature, group in lane_data.groupby("temperature_k", sort=True):
            group = group.sort_values("molality_mol_per_kg")
            color = colors[float(temperature)]
            ax.errorbar(group["molality_mol_per_kg"], group["experimental_osmotic_coefficient"], yerr=0.005,
                        color=color, marker=markers[float(temperature)], markerfacecolor="white", linestyle="none",
                        markeredgewidth=1.4, capsize=2, label=rf"Tables 8–9, {temperature:.0f} K", zorder=3)
            ax.plot(group["molality_mol_per_kg"], group["fitted_model"], color=color, linewidth=2.0,
                    label=rf"Fitted $k_{{ij}}$, {temperature:.0f} K")
        ax.set_title(TITLES[lane])
        ax.set_xlabel(r"MDEAHCl molality / mol kg$_{\mathrm{water}}^{-1}$")
        ax.grid(True, color="#D8D8D8", linewidth=0.7)
    fig.suptitle(r"Wangler 2018 Figure 10: water–MDEAH$^+$ $k_{ij}$ fitted per isotherm")
    axes[0].set_ylabel(r"Osmotic coefficient, $\phi$ / 1")
    axes[0].set_ylim(0.80, 0.90)
    axes[0].legend(frameon=False, fontsize=8.0)
    for suffix, kwargs in (("svg", {}), ("png", {"dpi": 200}), ("pdf", {})):
        fig.savefig(OUTPUT / f"wangler-2018-mdeah-osmotic-fit.{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)


if __name__ == "__main__":
    main()

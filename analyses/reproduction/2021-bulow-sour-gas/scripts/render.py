from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CAMPAIGN = Path(__file__).resolve().parents[1]
DATA = CAMPAIGN / "results" / "data"
RENDERED = CAMPAIGN / "results" / "rendered"


def _rows(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _binary_panel(
    axis: plt.Axes,
    rows: list[dict[str, str]],
    system: str,
    pressure_scale: float,
    pressure_unit: str,
    fitted_rows: list[dict[str, str]] | None = None,
) -> None:
    grouped: dict[float, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row["system"] == system:
            grouped[float(row["temperature_k"])].append(row)
    fitted_grouped: dict[float, list[dict[str, str]]] = defaultdict(list)
    for row in fitted_rows or []:
        if row["system"] == system:
            fitted_grouped[float(row["temperature_k"])].append(row)
    colors = plt.cm.viridis_r(  # type: ignore[attr-defined]
        [index / max(1, len(grouped) - 1) for index in range(len(grouped))]
    )
    for color, (temperature, values) in zip(
        colors, sorted(grouped.items()), strict=True
    ):
        values.sort(key=lambda row: float(row["sulfolane_liquid_mole_fraction"]))
        x = [float(row["sulfolane_liquid_mole_fraction"]) for row in values]
        observed = [float(row["experimental_pressure"]) for row in values]
        axis.scatter(x, observed, s=22, color=color, marker="o", facecolors="none")
        evaluated = [row for row in values if row.get("predicted_pressure_pa")]
        if evaluated:
            axis.plot(
                [float(row["sulfolane_liquid_mole_fraction"]) for row in evaluated],
                [
                    float(row["predicted_pressure_pa"]) / pressure_scale
                    for row in evaluated
                ],
                color=color,
                linewidth=1.5,
                label=f"{temperature:g} K",
            )
        fitted = [
            row
            for row in sorted(
                fitted_grouped.get(temperature, []),
                key=lambda item: float(item["sulfolane_liquid_mole_fraction"]),
            )
            if row.get("predicted_pressure_pa")
        ]
        if fitted:
            axis.plot(
                [float(row["sulfolane_liquid_mole_fraction"]) for row in fitted],
                [
                    float(row["predicted_pressure_pa"]) / pressure_scale
                    for row in fitted
                ],
                color=color,
                linewidth=1.5,
                linestyle="--",
            )
    axis.set_xlabel(r"liquid $x_{\mathrm{sulfolane}}$")
    axis.set_ylabel(f"pressure / {pressure_unit}")
    axis.grid(alpha=0.2)
    axis.legend(frameon=False, fontsize=8, ncol=2, title="temperature")


def main() -> None:
    binary = _rows("bulow-2021-binary-plotted-data.csv")
    fitted = _rows("bulow-2021-sulfolane-co2-fitted-data.csv")
    reactive_path = DATA / "bulow-2021-reactive-mdea-co2-plotted-data.csv"
    reactive = _rows(reactive_path.name) if reactive_path.exists() else []
    with plt.rc_context(
        {
            "figure.dpi": 150,
            "font.size": 9,
            "font.family": "serif",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    ):
        figure, axes = plt.subplots(1, 3, figsize=(11.2, 3.45), constrained_layout=True)
        _binary_panel(
            axes[0], binary, "sulfolane-co2", 1.0e5, "bar", fitted_rows=fitted
        )
        axes[0].set_title("(a) Sulfolane + CO₂")
        _binary_panel(axes[1], binary, "sulfolane-water", 1.0e3, "kPa")
        axes[1].set_title("(b) Sulfolane + water")

        evaluated = [row for row in reactive if row.get("predicted_total_pressure_pa")]
        if evaluated:
            evaluated.sort(key=lambda row: float(row["temperature_k"]))
            temperature = [float(row["temperature_k"]) for row in evaluated]
            observed = [
                float(row["experimental_total_pressure_pa"]) / 1.0e3
                for row in evaluated
            ]
            predicted = [
                float(row["predicted_total_pressure_pa"]) / 1.0e3 for row in evaluated
            ]
            axes[2].plot(temperature, observed, "o-", label="Kuranov experiment")
            axes[2].plot(temperature, predicted, "s--", label="Bülow model replay")
        else:
            axes[2].text(
                0.5,
                0.5,
                "No reactive states evaluated",
                ha="center",
                va="center",
                transform=axes[2].transAxes,
            )
        axes[2].set_xlabel("temperature / K")
        axes[2].set_ylabel("total pressure / kPa")
        axes[2].set_title("(c) 19 wt% MDEA + CO₂")
        axes[2].grid(alpha=0.2)
        axes[2].legend(frameon=False, fontsize=8)

        figure.suptitle(
            "Bülow et al. (2021): declared-topology installed-engine replay",
            fontsize=11,
        )
        figure.legend(
            handles=(
                Line2D(
                    [],
                    [],
                    color="black",
                    marker="o",
                    markerfacecolor="none",
                    linestyle="none",
                    label="digitized binary data",
                ),
                Line2D(
                    [],
                    [],
                    color="black",
                    linestyle="-",
                    label=r"paper $k_{ij}=0.055$",
                ),
                Line2D(
                    [],
                    [],
                    color="black",
                    linestyle="--",
                    label=r"generic fit $k_{ij}=-0.05$",
                ),
            ),
            frameon=False,
            fontsize=8,
            loc="outside lower center",
            ncol=3,
        )
        RENDERED.mkdir(parents=True, exist_ok=True)
        stem = RENDERED / "bulow-2021-sour-gas-reproduction"
        for suffix in ("png", "svg", "pdf"):
            figure.savefig(stem.with_suffix(f".{suffix}"), bbox_inches="tight")
        plt.close(figure)


if __name__ == "__main__":
    main()

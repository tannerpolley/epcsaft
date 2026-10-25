"""Render the retained Figiel 2025 MIAC predictions without recomputing them."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


SALTS = ("LiCl", "NaCl", "KCl", "LiBr", "NaBr", "KBr")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--output-base", type=Path, required=True)
    arguments = parser.parse_args()

    with arguments.predictions.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    plt.rcParams.update(
        {
            "font.family": "serif",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
        }
    )
    figure, axes = plt.subplots(2, 3, figsize=(11.2, 7.2))
    for axis, salt in zip(axes.flat, SALTS, strict=True):
        selected = [row for row in rows if row["salt"] == salt]
        molality = [float(row["molality_mol_kg"]) for row in selected]
        literature = [float(row["gamma_pm_m_literature"]) for row in selected]
        model = [float(row["gamma_pm_m_model"]) for row in selected]
        axis.scatter(
            molality,
            literature,
            s=25,
            facecolors="white",
            edgecolors="#171717",
            linewidths=0.9,
            label="Hamer–Wu (1972)",
            zorder=3,
        )
        axis.plot(
            molality,
            model,
            color="#9C2F1F",
            linewidth=1.8,
            marker="o",
            markersize=2.4,
            label="Provider (fresh installed artifact)",
            zorder=2,
        )
        axis.set_xscale("log")
        axis.set_title(salt)
        axis.grid(True, which="major", color="#D7D7D7", linewidth=0.6)
        axis.tick_params(direction="out")

    figure.supxlabel(r"Formula-unit molality, $m$ / mol kg$^{-1}$", y=0.025)
    figure.supylabel(r"Mean ionic activity coefficient, $\gamma_{\pm}^{m}$", x=0.015)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.915),
        ncols=2,
        frameon=False,
    )
    figure.suptitle("Aqueous alkali-halide MIACs at 298.15 K and 1 bar", y=0.975)
    figure.subplots_adjust(
        left=0.08,
        right=0.985,
        bottom=0.10,
        top=0.82,
        wspace=0.18,
        hspace=0.30,
    )
    for suffix, options in ((".png", {"dpi": 220}), (".svg", {}), (".pdf", {})):
        figure.savefig(arguments.output_base.with_suffix(suffix), bbox_inches="tight", **options)
    svg_path = arguments.output_base.with_suffix(".svg")
    normalized_svg = "\n".join(
        line.rstrip()
        for line in svg_path.read_text(encoding="utf-8").splitlines()
    )
    svg_path.write_text(f"{normalized_svg}\n", encoding="utf-8")
    plt.close(figure)


if __name__ == "__main__":
    main()

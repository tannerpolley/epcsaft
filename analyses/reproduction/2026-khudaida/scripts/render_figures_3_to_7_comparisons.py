"""Render retained exploratory comparisons to Khudaida Figures 3-7."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from generate_figures_3_to_7_comparison_data import comparison_inputs

CAMPAIGN_ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN_RESULTS = (
    CAMPAIGN_ROOT / "results" / "data" / "khudaida-2026-figures-3-to-7-held2.csv"
)
FIGURES = {
    3: {"temperature_K": 303.15, "salt_wt_percent": 5, "feed_count": 8, "table": 3},
    4: {"temperature_K": 313.15, "salt_wt_percent": 5, "feed_count": 7, "table": 3},
    5: {"temperature_K": 293.15, "salt_wt_percent": 10, "feed_count": 5, "table": 4},
    6: {"temperature_K": 303.15, "salt_wt_percent": 10, "feed_count": 5, "table": 4},
    7: {"temperature_K": 313.15, "salt_wt_percent": 10, "feed_count": 6, "table": 4},
}
HEIGHT = math.sqrt(3.0) / 2.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ternary_xy(frame: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    x = frame["x_isobutanol_salt_free"] + 0.5 * frame["x_ethanol_salt_free"]
    y = HEIGHT * frame["x_ethanol_salt_free"]
    return x, y


def draw_grid(axis: plt.Axes) -> None:
    triangle = [(0.0, 0.0), (1.0, 0.0), (0.5, HEIGHT), (0.0, 0.0)]
    axis.plot(*zip(*triangle, strict=True), color="#202020", linewidth=1.2)
    for fraction in (0.2, 0.4, 0.6, 0.8):
        grid_lines = (
            ((fraction, 0.0, 1.0 - fraction), (fraction, 1.0 - fraction, 0.0)),
            ((0.0, fraction, 1.0 - fraction), (1.0 - fraction, fraction, 0.0)),
            ((0.0, 1.0 - fraction, fraction), (1.0 - fraction, 0.0, fraction)),
        )
        for left, right in grid_lines:
            x = [point[2] + 0.5 * point[1] for point in (left, right)]
            y = [HEIGHT * point[1] for point in (left, right)]
            axis.plot(x, y, color="#d4d4d4", linewidth=0.55, zorder=0)


def tie_lines(axis: plt.Axes, data: pd.DataFrame, dataset: str, **style: object) -> None:
    selected = data[data["dataset"].eq(dataset)]
    for point_id in sorted(selected["point_id"].unique()):
        pair = selected[selected["point_id"].eq(point_id)]
        if len(pair) == 2:
            axis.plot(*ternary_xy(pair), **style)


def scatter(
    axis: plt.Axes, data: pd.DataFrame, dataset: str, label: str, **style: object
) -> None:
    selected = data[data["dataset"].eq(dataset)]
    axis.scatter(*ternary_xy(selected), label=label, **style)


def render(figure_id: int) -> None:
    spec = FIGURES[figure_id]
    figure_root = CAMPAIGN_ROOT / "figures" / f"figure-{figure_id:02d}"
    plotted = figure_root / "output" / f"khudaida-figure-{figure_id:02d}-active-comparison.csv"
    output_stem = plotted.with_suffix("")
    provenance = figure_root / "output" / "comparison-provenance.json"
    data = pd.read_csv(plotted)
    active = data[data["dataset"].eq("active_epcsaft")]
    accepted_count = int(active["point_id"].nunique())

    with plt.rc_context(
        {
            "font.family": "serif",
            "font.size": 10.5,
            "axes.titleweight": "semibold",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    ):
        figure, axis = plt.subplots(figsize=(8.2, 6.7), constrained_layout=True)
        draw_grid(axis)
        tie_lines(
            axis,
            data,
            "experimental",
            color="#6f6f6f",
            linestyle="--",
            linewidth=0.8,
            alpha=0.75,
            zorder=1,
        )
        tie_lines(
            axis,
            data,
            "active_epcsaft",
            color="#1764ab",
            linestyle="-",
            linewidth=1.05,
            alpha=0.85,
            zorder=2,
        )
        scatter(
            axis,
            data,
            "experimental",
            f"Experiment (Table {spec['table']})",
            marker="o",
            s=40,
            facecolors="none",
            edgecolors="#292929",
            linewidths=1.1,
            zorder=4,
        )
        scatter(
            axis,
            data,
            "paper_epcsaft",
            f"Published ePC-SAFT (Figure {figure_id})",
            marker="s",
            s=42,
            facecolors="none",
            edgecolors="#c43b32",
            linewidths=1.2,
            zorder=5,
        )
        if not active.empty:
            scatter(
                axis,
                data,
                "active_epcsaft",
                "Exploratory active ePC-SAFT + HELD2",
                marker="D",
                s=34,
                color="#1764ab",
                edgecolors="white",
                linewidths=0.45,
                zorder=6,
            )
        scatter(
            axis,
            data,
            "feed",
            f"Figure {figure_id} feeds",
            marker="^",
            s=43,
            color="#178353",
            edgecolors="#124c35",
            linewidths=0.5,
            zorder=5,
        )
        active_organic = active[active["phase"].eq("organic")]
        x, y = ternary_xy(active_organic)
        for point_id, x_value, y_value in zip(
            active_organic["point_id"], x, y, strict=True
        ):
            axis.annotate(
                str(int(point_id)),
                (x_value, y_value),
                xytext=(3, 4),
                textcoords="offset points",
                color="#174e7d",
                fontsize=8,
            )

        axis.text(-0.025, -0.035, r"Water, $x_\mathrm{w}$", ha="left", va="top")
        axis.text(1.025, -0.035, r"Isobutanol, $x_\mathrm{iBuOH}$", ha="right", va="top")
        axis.text(0.5, HEIGHT + 0.035, r"Ethanol, $x_\mathrm{EtOH}$", ha="center")
        axis.set_title(
            f"Khudaida Figure {figure_id}: exploratory installed-artifact phase splits",
            pad=18,
        )
        axis.text(
            0.03,
            0.78,
            (
                "Salt-free solvent basis\n"
                f"{spec['temperature_K']:.2f} K, 100 kPa, "
                f"{spec['salt_wt_percent']} wt% NaCl\n"
                f"{accepted_count}/{spec['feed_count']} native splits accepted\n"
                "Figure 2 bundle extrapolated"
            ),
            transform=axis.transAxes,
            color="#404040",
            fontsize=9,
            linespacing=1.35,
        )
        axis.set_aspect("equal")
        axis.set_xlim(-0.06, 1.06)
        axis.set_ylim(-0.07, HEIGHT + 0.08)
        axis.set_xticks([])
        axis.set_yticks([])
        axis.spines[:].set_visible(False)
        axis.legend(loc="upper right", frameon=False, fontsize=9, handletextpad=0.7)

        svg_output = output_stem.with_suffix(".svg")
        figure.savefig(
            svg_output,
            metadata={"Title": axis.get_title(), "Date": "2026-07-31"},
        )
        svg_output.write_text(
            "\n".join(
                line.rstrip()
                for line in svg_output.read_text(encoding="utf-8").splitlines()
            )
            + "\n",
            encoding="utf-8",
        )
        figure.savefig(
            output_stem.with_suffix(".png"),
            dpi=220,
            metadata={"Title": axis.get_title(), "Author": "ePC-SAFT analysis"},
        )
        figure.savefig(
            output_stem.with_suffix(".pdf"),
            metadata={
                "Title": axis.get_title(),
                "Author": "ePC-SAFT analysis",
                "Creator": "Matplotlib",
                "CreationDate": None,
                "ModDate": None,
            },
        )
        plt.close(figure)
    outputs = {
        suffix: sha256(output_stem.with_suffix(suffix))
        for suffix in (".svg", ".png", ".pdf")
    }
    provenance.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "campaign_results": {
                    "path": str(CAMPAIGN_RESULTS.relative_to(CAMPAIGN_ROOT)),
                    "sha256": sha256(CAMPAIGN_RESULTS),
                },
                "plotted_data": {
                    "path": str(plotted.relative_to(figure_root)),
                    "sha256": sha256(plotted),
                },
                "source_inputs": comparison_inputs(figure_id),
                "coordinate_normalization": {
                    "method": "divide each ternary row by its retained input sum",
                    "maximum_absolute_correction": float(
                        data["input_coordinate_sum"].sub(1.0).abs().max()
                    ),
                },
                "outputs": outputs,
                "source_figure": f"source/figure-{figure_id:02d}.png",
                "accepted_native_splits": accepted_count,
                "feed_count": spec["feed_count"],
                "parameter_domain_status": "exploratory_parameter_extrapolation",
                "claim_scope": (
                    f"{accepted_count} accepted installed-artifact phase splits compared "
                    f"with retained Khudaida Figure {figure_id} and Table "
                    f"{spec['table']} digitizations under an explicit Figure 2 "
                    "parameter-bundle extrapolation assumption."
                ),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figures", type=int, nargs="+", choices=tuple(FIGURES), default=tuple(FIGURES))
    args = parser.parse_args()
    for figure_id in dict.fromkeys(args.figures):
        render(figure_id)


if __name__ == "__main__":
    main()

"""Render the retained active-model comparison to Khudaida Figure 2."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

FIGURE_ROOT = Path(__file__).resolve().parents[1]
PLOTTED = FIGURE_ROOT / "output" / "khudaida-figure-02-active-comparison.csv"
OUTPUT_STEM = FIGURE_ROOT / "output" / "khudaida-figure-02-active-comparison"
CAMPAIGN_RESULTS = (
    FIGURE_ROOT.parents[1] / "results" / "data" / "khudaida-2026-figure-2-held2.csv"
)
PROVENANCE = FIGURE_ROOT / "output" / "comparison-provenance.json"
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


def tie_lines(
    axis: plt.Axes, data: pd.DataFrame, dataset: str, **style: object
) -> None:
    selected = data[data["dataset"].eq(dataset)]
    for point_id in range(1, 9):
        pair = selected[selected["point_id"].eq(point_id)]
        if len(pair) != 2:
            continue
        x, y = ternary_xy(pair)
        axis.plot(x, y, **style)


def scatter(
    axis: plt.Axes,
    data: pd.DataFrame,
    dataset: str,
    label: str,
    **style: object,
) -> None:
    selected = data[data["dataset"].eq(dataset)]
    x, y = ternary_xy(selected)
    axis.scatter(x, y, label=label, **style)


def main() -> None:
    data = pd.read_csv(PLOTTED)
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
            "Experiment (Table 3)",
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
            "Published ePC-SAFT (Figure 2)",
            marker="s",
            s=42,
            facecolors="none",
            edgecolors="#c43b32",
            linewidths=1.2,
            zorder=5,
        )
        scatter(
            axis,
            data,
            "active_epcsaft",
            "Active installed ePC-SAFT + HELD2",
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
            "Figure 2 feeds",
            marker="^",
            s=43,
            color="#178353",
            edgecolors="#124c35",
            linewidths=0.5,
            zorder=5,
        )

        active_organic = data[
            data["dataset"].eq("active_epcsaft") & data["phase"].eq("organic")
        ]
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
        axis.text(
            1.025, -0.035, r"Isobutanol, $x_\mathrm{iBuOH}$", ha="right", va="top"
        )
        axis.text(0.5, HEIGHT + 0.035, r"Ethanol, $x_\mathrm{EtOH}$", ha="center")
        axis.set_title(
            "Khudaida Figure 2: active installed-artifact phase split",
            pad=18,
        )
        axis.text(
            0.03,
            0.80,
            "Salt-free solvent basis\n293.15 K, 100 kPa, 5 wt% NaCl",
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
        axis.legend(
            loc="upper right",
            frameon=False,
            fontsize=9,
            handletextpad=0.7,
        )

        svg_output = OUTPUT_STEM.with_suffix(".svg")
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
            OUTPUT_STEM.with_suffix(".png"),
            dpi=220,
            metadata={"Title": axis.get_title(), "Author": "ePC-SAFT analysis"},
        )
        figure.savefig(
            OUTPUT_STEM.with_suffix(".pdf"),
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
        suffix: sha256(OUTPUT_STEM.with_suffix(suffix))
        for suffix in (".svg", ".png", ".pdf")
    }
    feed_mapping = {
        str(int(row.feed_point_id)): int(row.point_id)
        for row in data[data["dataset"].eq("feed")][
            ["feed_point_id", "point_id"]
        ].itertuples(index=False)
    }
    if len(feed_mapping) != 8:
        raise AssertionError("Figure 2 comparison requires eight feed identities")
    PROVENANCE.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "campaign_results": {
                    "path": str(CAMPAIGN_RESULTS.relative_to(FIGURE_ROOT.parents[1])),
                    "sha256": sha256(CAMPAIGN_RESULTS),
                },
                "plotted_data": {
                    "path": str(PLOTTED.relative_to(FIGURE_ROOT)),
                    "sha256": sha256(PLOTTED),
                },
                "coordinate_normalization": {
                    "method": "divide each ternary row by its retained input sum",
                    "maximum_absolute_correction": float(
                        data["input_coordinate_sum"].sub(1.0).abs().max()
                    ),
                },
                "indexing": {
                    "point_id_role": "experimental Table 3 tie-line identity",
                    "feed_point_id_role": "Figure 2 solver-feed identity",
                    "feed_to_experimental_tie_line": feed_mapping,
                },
                "outputs": outputs,
                "source_figure": "source/figure-02.png",
                "claim_scope": (
                    "Eight accepted installed-artifact phase splits compared with "
                    "the retained Khudaida Figure 2 and Table 3 digitizations."
                ),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()

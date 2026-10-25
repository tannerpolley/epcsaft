from __future__ import annotations

import csv
from datetime import UTC, datetime
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt

FIGURE = Path(__file__).resolve().parents[1] / "results/analysis/neutral-mea-water/figures"
OUTPUT = FIGURE / "output"


def main() -> None:
    matplotlib.rcParams["svg.hashsalt"] = "issue-114-neutral-mea-water"
    with (OUTPUT / "topology-selection.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        rows = list(csv.DictReader(stream))
    figure, axis = plt.subplots(figsize=(6.4, 3.4), layout="constrained")
    for index, row in enumerate(rows):
        selected = row["selected"] == "true"
        if not row["grouped_mean_score"]:
            axis.text(
                0.98,
                index,
                "phase failure — no full grouped score",
                ha="right",
                va="center",
                color="#666666",
                transform=axis.get_yaxis_transform(),
            )
            continue
        axis.errorbar(
            float(row["grouped_mean_score"]),
            index,
            xerr=float(row["grouped_standard_error"]),
            marker="*" if selected else "o",
            markersize=11 if selected else 7,
            color="#0072B2" if row["eligible"] == "true" else "#666666",
            capsize=4,
            linestyle="none",
            label="selected" if selected else None,
        )
    axis.set_xscale("log")
    axis.set_yticks(range(len(rows)), [row["scheme"] for row in rows])
    axis.set_xlabel("Held-pressure grouped mean scaled square residual")
    axis.set_ylabel("MEA association topology")
    axis.set_title("2B remains the sole topology inside the one-SE threshold")
    axis.set_ylim(-0.4, 2.5)
    axis.grid(axis="x", alpha=0.25)
    svg = OUTPUT / "topology-selection.svg"
    figure.savefig(svg, dpi=180, metadata={"Date": "2026-08-30"})
    svg.write_text(
        "\n".join(
            line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines()
        )
        + "\n",
        encoding="utf-8",
    )
    figure.savefig(OUTPUT / "topology-selection.png", dpi=180)
    figure.savefig(
        OUTPUT / "topology-selection.pdf",
        dpi=180,
        metadata={"CreationDate": datetime(2026, 8, 30, tzinfo=UTC)},
    )
    plt.close(figure)


if __name__ == "__main__":
    main()

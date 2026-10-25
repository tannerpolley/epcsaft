from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1] / "results/analysis/pure-components/figures/source-domain-coverage/scripts"
OUTPUT = HERE.parent / "output"
DATA = OUTPUT / "source-domain-coverage.csv"


def main() -> None:
    with DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    labels = [row["family_id"] for row in rows]
    colors = {"fixed": "#0072B2", "candidate": "#D55E00"}
    mpl.rcParams.update({"svg.hashsalt": "issue-113-source-domains", "font.size": 9})
    figure, axis = plt.subplots(figsize=(8.0, 3.8), constrained_layout=True)
    for index, row in enumerate(rows):
        lower = float(row["temperature_min_k"])
        upper = float(row["temperature_max_k"])
        color = colors[row["classification"]]
        if lower == upper:
            axis.scatter(lower, index, color=color, marker="o", s=42, zorder=3)
        else:
            axis.plot(
                (lower, upper),
                (index, index),
                color=color,
                linewidth=5,
                solid_capstyle="round",
            )
            axis.scatter((lower, upper), (index, index), color=color, s=20, zorder=3)
    axis.axvspan(
        293.15, 393.15, color="#999999", alpha=0.12, label="active campaign window"
    )
    axis.set_yticks(range(len(labels)), labels)
    axis.set_xlabel("Source-supported temperature (K)")
    axis.set_title(
        "Pure-component source domains; intervals are not extrapolation claims"
    )
    axis.grid(axis="x", alpha=0.25)
    fixed = mpl.lines.Line2D(
        [], [], color=colors["fixed"], linewidth=5, label="fixed source family"
    )
    candidate = mpl.lines.Line2D(
        [], [], color=colors["candidate"], linewidth=5, label="candidate source family"
    )
    active = mpl.patches.Patch(
        color="#999999", alpha=0.2, label="active campaign window"
    )
    axis.legend(handles=(fixed, candidate, active), loc="lower right", frameon=False)
    svg = OUTPUT / "source-domain-coverage.svg"
    figure.savefig(svg, metadata={"Date": None})
    svg.write_text(
        "\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n",
        encoding="utf-8",
    )
    figure.savefig(
        OUTPUT / "source-domain-coverage.png", dpi=180, metadata={"Date": None}
    )
    figure.savefig(
        OUTPUT / "source-domain-coverage.pdf",
        metadata={"CreationDate": None, "ModDate": None},
    )
    plt.close(figure)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Render retained May 2015 source/model flash rows without package imports."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--output-stem", required=True, type=Path)
    return parser


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 17:
        raise ValueError("expected all 17 May coexistence rows")
    return rows


def _optional_float(value: str) -> float:
    return float(value) if value else math.nan


def _normalize_svg(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8")


def render(csv_path: Path, output_stem: Path) -> None:
    rows = _read_rows(csv_path)
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 9.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.22,
            "grid.linewidth": 0.6,
            "svg.hashsalt": "consumer-slice-2-may-2015-methane-ethane-flash-v1",
        }
    )
    fig, (phase_ax, residual_ax) = plt.subplots(1, 2, figsize=(9.0, 4.4))
    source_x = [float(row["source_x_methane"]) for row in rows]
    source_y = [float(row["source_y_methane"]) for row in rows]
    uc_x = [float(row["uc_x_methane"]) for row in rows]
    uc_y = [float(row["uc_y_methane"]) for row in rows]
    model_x = [_optional_float(row["model_x_methane"]) for row in rows]
    model_y = [_optional_float(row["model_y_methane"]) for row in rows]
    phase_ax.errorbar(
        source_x,
        source_y,
        xerr=uc_x,
        yerr=uc_y,
        fmt="o",
        ms=4.8,
        mfc="white",
        mec="#222222",
        ecolor="#777777",
        elinewidth=0.7,
        capsize=1.8,
        label=r"May et al. Table 5 ($u_c$)",
        zorder=3,
    )
    solved = [i for i, value in enumerate(model_x) if math.isfinite(value)]
    phase_ax.scatter(
        [model_x[i] for i in solved],
        [model_y[i] for i in solved],
        marker="^",
        s=34,
        color="#0072B2",
        label="Flash (fresh installed artifacts)",
        zorder=4,
    )
    for i in solved:
        phase_ax.plot(
            [source_x[i], model_x[i]],
            [source_y[i], model_y[i]],
            color="#BBBBBB",
            lw=0.7,
            zorder=1,
        )
    rejected = [i for i, value in enumerate(model_x) if not math.isfinite(value)]
    phase_ax.scatter(
        [source_x[i] for i in rejected],
        [source_y[i] for i in rejected],
        marker="X",
        s=62,
        color="#CC3311",
        label="Local state rejected",
        zorder=5,
    )
    phase_ax.set_title("Methane liquid/vapor coexistence")
    phase_ax.set_xlabel(r"Liquid methane fraction, $x_{CH_4}$")
    phase_ax.set_ylabel(r"Vapor methane fraction, $y_{CH_4}$")
    phase_ax.legend(frameon=False, fontsize=8)

    indices = list(range(1, len(rows) + 1))
    x_norm = [_optional_float(row["x_normalized_error"]) for row in rows]
    y_norm = [_optional_float(row["y_normalized_error"]) for row in rows]
    residual_ax.axhspan(-1.0, 1.0, color="#009E73", alpha=0.09, label=r"Frozen $3u_c$ allowance")
    residual_ax.axhline(1.0, color="#228833", ls="--", lw=0.9)
    residual_ax.axhline(-1.0, color="#228833", ls="--", lw=0.9)
    residual_ax.plot(indices, x_norm, "o-", ms=4.2, lw=1.0, color="#0072B2", label=r"$(x_{model}-x_{src})/(3u_{c,x})$")
    residual_ax.plot(indices, y_norm, "s-", ms=4.0, lw=1.0, color="#D55E00", label=r"$(y_{model}-y_{src})/(3u_{c,y})$")
    for i in rejected:
        residual_ax.scatter(i + 1, 0.0, marker="X", s=62, color="#CC3311", zorder=5)
        residual_ax.annotate("local rejection", (i + 1, 0.0), xytext=(4, 8), textcoords="offset points", fontsize=7.5, rotation=30)
    residual_ax.set_title("Composition error / frozen allowance")
    residual_ax.set_xlabel("Table 5 row")
    residual_ax.set_ylabel("Signed normalized composition error")
    residual_ax.set_xticks(indices)
    residual_ax.legend(frameon=False, fontsize=7.5, loc="lower left")

    admitted = sum(row["row_admission"] == "PASS" for row in rows)
    misses = sum(row["composition_agreement_status"] == "FAIL" for row in rows)
    rejected_count = len(rows) - len(solved)
    fig.suptitle("May et al. (2015) methane/ethane local two-phase flash", fontsize=12)
    fig.tight_layout(rect=(0.0, 0.09, 1.0, 0.92), w_pad=1.6)
    fig.text(
        0.5,
        0.02,
        f"{admitted}/17 admitted rows; {misses} solved model/data misses; {rejected_count} local package rejection. No globality certificate.",
        ha="center",
        fontsize=8.5,
    )
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    svg_path = output_stem.with_suffix(".svg")
    fig.savefig(svg_path, metadata={"Date": None})
    _normalize_svg(svg_path)
    fig.savefig(output_stem.with_suffix(".png"), dpi=220, metadata={"Software": "Matplotlib"})
    fig.savefig(
        output_stem.with_suffix(".pdf"),
        metadata={"CreationDate": None, "ModDate": None, "Title": "May 2015 methane ethane flash validation"},
    )
    plt.close(fig)


def main() -> int:
    args = build_parser().parse_args()
    render(args.csv.resolve(), args.output_stem.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

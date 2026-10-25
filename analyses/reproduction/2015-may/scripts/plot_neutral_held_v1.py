#!/usr/bin/env python3
"""Render neutral-held-v1 retained CSV evidence without package imports."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases-csv", required=True, type=Path)
    parser.add_argument("--surface-csv", required=True, type=Path)
    parser.add_argument("--record", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _normalize_svg(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8")


def _save_bundle(fig: plt.Figure, stem: Path, title: str) -> dict[str, dict[str, str]]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    svg = stem.with_suffix(".svg")
    png = stem.with_suffix(".png")
    pdf = stem.with_suffix(".pdf")
    fig.savefig(svg, metadata={"Date": None})
    _normalize_svg(svg)
    fig.savefig(png, dpi=220, metadata={"Software": "Matplotlib"})
    fig.savefig(pdf, metadata={"CreationDate": None, "ModDate": None, "Title": title})
    plt.close(fig)
    return {
        suffix: {"path": str(path), "sha256": sha256_file(path)}
        for suffix, path in (("svg", svg), ("png", png), ("pdf", pdf))
    }


def render(cases_csv: Path, surface_csv: Path, record_path: Path, output_dir: Path) -> None:
    cases = _read_csv(cases_csv)
    surfaces = _read_csv(surface_csv)
    may = [row for row in cases if row["case_role"] == "may_coexistence_midpoint"]
    if len(cases) != 18 or len(may) != 17:
        raise ValueError("expected 18 cases including all 17 May coexistence rows")
    representative = [row for row in surfaces if row["audit_case_id"] == "may2015-row-012-midpoint"]
    if len(representative) != 3003:
        raise ValueError("representative surface must contain 1001 points for three pressure-state branches")

    mpl.rcParams.update({
        "font.family": "serif", "font.size": 9.5, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.22,
        "grid.linewidth": 0.6, "svg.hashsalt": "neutral-held-v1-validation",
    })
    fig, (phase_ax, residual_ax) = plt.subplots(1, 2, figsize=(9.2, 4.4))
    source_x = [float(row["source_x_methane"]) for row in may]
    source_y = [float(row["source_y_methane"]) for row in may]
    source_ux = [float(row["x_comparison_allowance"]) / 3.0 for row in may]
    source_uy = [float(row["y_comparison_allowance"]) / 3.0 for row in may]
    phase_ax.errorbar(
        source_x, source_y, xerr=source_ux, yerr=source_uy, fmt="o", ms=4.8,
        mfc="white", mec="#222222", ecolor="#777777", elinewidth=0.7,
        capsize=1.8, label=r"May et al. Table 5 ($u_c$)", zorder=3,
    )
    solved = [row for row in may if row["returned_phase_count"] == "2"]
    if solved:
        phase_ax.scatter(
            [float(row["phase_1_x_methane"]) for row in solved],
            [float(row["phase_2_x_methane"]) for row in solved],
            marker="^", s=34, color="#0072B2", label="HELD (fresh installed artifacts)", zorder=4,
        )
    else:
        phase_ax.text(
            0.52, 0.06, "No two-phase HELD state returned",
            transform=phase_ax.transAxes, fontsize=8.5, color="#7A1E12",
        )
    rejected = [row for row in may if row["composition_agreement_status"] != "PASS"]
    phase_ax.scatter(
        [float(row["source_x_methane"]) for row in rejected],
        [float(row["source_y_methane"]) for row in rejected],
        marker="X", s=62, color="#CC3311", label="Predictive non-admission", zorder=5,
    )
    phase_ax.set_title("Methane/ethane coexistence")
    phase_ax.set_xlabel(r"Liquid methane fraction, $x_{CH_4}$")
    phase_ax.set_ylabel(r"Vapor methane fraction, $y_{CH_4}$")
    phase_ax.legend(frameon=False, fontsize=7.7)

    indices = list(range(1, 18))
    x_norm = [float(row["x_signed_error"]) / float(row["x_comparison_allowance"]) if row["x_signed_error"] else math.nan for row in may]
    y_norm = [float(row["y_signed_error"]) / float(row["y_comparison_allowance"]) if row["y_signed_error"] else math.nan for row in may]
    residual_ax.axhspan(-1.0, 1.0, color="#009E73", alpha=0.09, label=r"Frozen $3u_c$ allowance")
    residual_ax.axhline(1.0, color="#228833", ls="--", lw=0.9)
    residual_ax.axhline(-1.0, color="#228833", ls="--", lw=0.9)
    residual_ax.plot(indices, x_norm, "o-", ms=4.2, lw=1.0, color="#0072B2", label=r"Liquid error / $3u_c$")
    residual_ax.plot(indices, y_norm, "s-", ms=4.0, lw=1.0, color="#D55E00", label=r"Vapor error / $3u_c$")
    if not any(math.isfinite(value) for value in x_norm + y_norm):
        residual_ax.text(
            0.5, 0.53, "No two-phase composition\ncomparison was evaluable",
            transform=residual_ax.transAxes, ha="center", va="center",
            fontsize=9.5, color="#7A1E12",
        )
    residual_ax.set_title("Experimental composition comparison")
    residual_ax.set_xlabel("Table 5 row")
    residual_ax.set_ylabel("Signed normalized error")
    residual_ax.set_xticks(indices)
    residual_ax.legend(frameon=False, fontsize=7.5, loc="lower left")
    admitted = sum(row["composition_agreement_status"] == "PASS" for row in may)
    fig.suptitle("Neutral HELD v1: May et al. (2015) validation", fontsize=12)
    fig.tight_layout(rect=(0.0, 0.08, 1.0, 0.92), w_pad=1.6)
    fig.text(0.5, 0.02, f"{admitted}/17 rows meet the frozen composition contract. Predictive agreement is separate from HELD solver status.", ha="center", fontsize=8.3)
    coexistence = _save_bundle(fig, output_dir / "neutral-held-v1-coexistence", "Neutral HELD v1 coexistence validation")

    surface_fig, surface_ax = plt.subplots(figsize=(7.2, 4.8))
    colors = {"single": "#555555", "vapor": "#D55E00", "liquid": "#0072B2"}
    markers = {"single": ".", "vapor": "x", "liquid": "+"}
    for branch in ("single", "vapor", "liquid"):
        branch_rows = [row for row in representative if row["branch"] == branch and row["status"] == "accepted"]
        surface_ax.plot(
            [float(row["x_methane"]) for row in branch_rows],
            [float(row["g_bar"]) for row in branch_rows],
            linestyle="none", marker=markers[branch], ms=2.5, alpha=0.65,
            color=colors[branch], label=f"Public {branch} pressure branch",
        )
    support = [row for row in representative if row["status"] == "accepted" and row["support_line_g_bar"]]
    support_by_x: dict[float, float] = {}
    for row in support:
        support_by_x[float(row["x_methane"])] = float(row["support_line_g_bar"])
    surface_ax.plot(sorted(support_by_x), [support_by_x[x] for x in sorted(support_by_x)], color="#000000", lw=1.2, label="HELD one-phase tangent")
    surface_ax.set_title("Representative finite sampled Gibbs audit: May row 012")
    surface_ax.set_xlabel(r"Methane mole fraction, $x_{CH_4}$")
    surface_ax.set_ylabel(r"Dimensionless molar Gibbs value, $\bar{g}$")
    surface_ax.legend(frameon=False, fontsize=8, ncol=2)
    surface_fig.tight_layout()
    sampled = _save_bundle(surface_fig, output_dir / "neutral-held-v1-sampled-gibbs", "Neutral HELD v1 sampled Gibbs audit")

    record = json.loads(record_path.read_text(encoding="utf-8"))
    record["outputs"]["plots"] = {"coexistence": coexistence, "sampled_gibbs": sampled}
    record["plot_environment"] = {
        "python": sys.version,
        "python_executable": sys.executable,
        "matplotlib": mpl.__version__,
        "platform": platform.platform(),
    }
    record["plot_command"] = [str(Path(__file__).resolve()), "--cases-csv", str(cases_csv), "--surface-csv", str(surface_csv), "--record", str(record_path), "--output-dir", str(output_dir)]
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def main() -> int:
    args = build_parser().parse_args()
    render(args.cases_csv.resolve(), args.surface_csv.resolve(), args.record.resolve(), args.output_dir.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

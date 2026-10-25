"""Render retained Figiel result CSVs without importing Provider code."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from .figure_specs import FIGURE_8_PANELS, FIGURE_9_PANELS
from .paths import PAPER_ROOT


ORANGE = "#d9891b"
GREEN = "#228b22"
GRAY = "0.45"
LIGHT_GRAY = "0.82"
BLUE = "#1f4da8"
FIGURE_9_LINE_COLORS = {"0.4": GREEN, "0.8": "black"}


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "font.family": "DejaVu Serif",
            "axes.linewidth": 1.0,
            "axes.grid": False,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.frameon": False,
            "mathtext.default": "regular",
        }
    )


def load_rows(number: int) -> list[dict[str, str]]:
    path = PAPER_ROOT / "figures" / f"figure-{number:02d}" / "results" / "data" / f"figure-{number:02d}-series.csv"
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def select(rows: list[dict[str, str]], *, source_type: str | None = None, category: str | None = None, series: str | None = None) -> list[dict[str, str]]:
    return [
        row
        for row in rows
        if (source_type is None or row["source_type"] == source_type)
        and (category is None or row["category"] == category)
        and (series is None or row["series_id"] == series)
    ]


def xy(rows: list[dict[str, str]]) -> tuple[list[float], list[float]]:
    ordered = sorted(rows, key=lambda row: (int(row["display_order"]), float(row["x_value"])))
    return [float(row["x_value"]) for row in ordered], [float(row["y_value"]) for row in ordered]


def save(figure: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    figure.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    figure.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)


def render_figure_1(rows: list[dict[str, str]], output: Path) -> None:
    x, y = xy(rows)
    figure, axis = plt.subplots(figsize=(4.2, 3.2))
    axis.plot(x, y, color="black", linewidth=1.8)
    axis.set_xlim(0.0, 2.0)
    axis.set_ylim(0.0, 100.0)
    axis.set_xticks([0.0, 1.0, 2.0])
    axis.set_xticklabels([r"$r_i$", r"$r_i + \Delta r_i$", ""])
    axis.set_yticks([0.0, 25.0, 50.0, 75.0, 100.0])
    axis.set_ylabel(r"$\epsilon_r$ / -")
    axis.text(0.28, 12.5, r"$\epsilon_{r,\mathrm{ion}}$", fontsize=10)
    axis.text(1.42, 80.5, r"$\epsilon_{r,\mathrm{bulk}}$", fontsize=10)
    save(figure, output)


def render_figure_2(rows: list[dict[str, str]], output: Path) -> None:
    configure_style()
    styles = {
        "NaCl": (GRAY, "o", LIGHT_GRAY),
        "NaBr": (GREEN, "s", GREEN),
        "LiCl": (ORANGE, "v", ORANGE),
    }
    figure, axis = plt.subplots(figsize=(4.2, 3.2))
    for category, (color, marker, face) in styles.items():
        lit = select(rows, source_type="literature", category=category)
        x, y = xy(lit)
        axis.scatter(x, y, s=34, marker=marker, facecolor=face, edgecolor=color, linewidth=0.9, label=category)
        model = select(rows, source_type="model", category=category)
        x, y = xy(model)
        axis.plot(x, y, color="black", linewidth=1.5)
    axis.set_xlim(0.0, 0.30)
    axis.set_ylim(20.0, 85.0)
    axis.set_xlabel(r"$x_{\mathrm{ion}}$ / mol mol$^{-1}$")
    axis.set_ylabel(r"$\epsilon_r$ / -")
    axis.legend(loc="best", fontsize=8)
    save(figure, output)


def render_figure_3(rows: list[dict[str, str]], output: Path) -> None:
    configure_style()
    styles = {"water": (GRAY, "o", LIGHT_GRAY), "methanol": (BLUE, "^", BLUE), "ethanol": (GREEN, "s", GREEN)}
    figure, axis = plt.subplots(figsize=(4.2, 3.2))
    for category, (color, marker, face) in styles.items():
        source_type = "literature" if category == "water" else "digitized"
        x, y = xy(select(rows, source_type=source_type, category=category))
        axis.scatter(x, y, s=30, marker=marker, facecolor=face, edgecolor=color, linewidth=0.8, label=category.capitalize())
        x, y = xy(select(rows, source_type="model", category=category))
        axis.plot(x, y, color="black", linewidth=1.5)
    axis.set_xlim(0.0, 0.20)
    axis.set_ylim(0.0, 1.0)
    axis.set_xlabel(r"$x_{\mathrm{ion}}$ / mol mol$^{-1}$")
    axis.set_ylabel(r"$\epsilon_r / \epsilon_{r,\mathrm{solvent}}$ / -")
    axis.legend(loc="best", fontsize=8)
    save(figure, output)


def render_figure_4(rows: list[dict[str, str]], output: Path) -> None:
    configure_style()
    ions = ["Li", "Na", "K", "Cl", "Br", "I"]
    labels = [r"$Li^+$", r"$Na^+$", r"$K^+$", r"$Cl^-$", r"$Br^-$", r"$I^-$"]
    lit = {row["category"]: float(row["y_value"]) for row in select(rows, source_type="literature")}
    model = {row["category"]: float(row["y_value"]) for row in select(rows, source_type="model")}
    x = np.arange(len(ions), dtype=float)
    width = 0.30
    figure, axis = plt.subplots(figsize=(8.0, 4.2))
    axis.bar(x - width / 2, [lit[ion] for ion in ions], width=width, color=LIGHT_GRAY, edgecolor="black", linewidth=0.8, label="Literature")
    axis.bar(x + width / 2, [model[ion] for ion in ions], width=width, color=BLUE, edgecolor="black", linewidth=0.8, label="Provider (fresh installed artifact)")
    axis.set_xticks(x)
    axis.set_xticklabels(labels)
    axis.set_ylim(0.0, 800.0)
    axis.set_ylabel(r"$-\Delta G_i^{solv,\infty,x}$ / kJ mol$^{-1}$")
    axis.legend(loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.02))
    save(figure, output)


def render_figure_5(rows: list[dict[str, str]], output: Path) -> None:
    configure_style()
    styles = {"Li": (ORANGE, "^"), "Na": (GREEN, "s"), "K": (GRAY, "o")}
    panels = (("a)", ("LiCl", "NaCl", "KCl"), "Cl$^-$ salts in water"), ("b)", ("LiBr", "NaBr", "KBr"), "Br$^-$ salts in water"))
    figure, axes = plt.subplots(1, 2, figsize=(8.2, 3.9), sharey=True)
    for axis, (panel, salts, title) in zip(axes, panels, strict=True):
        for salt in salts:
            cation = salt[:2] if salt.startswith("Li") or salt.startswith("Na") else "K"
            color, marker = styles[cation]
            x, y = xy(select(rows, source_type="literature", category=salt))
            axis.scatter(x, y, marker=marker, s=26, facecolor="none", edgecolor=color, linewidth=1.0, label=salt)
            x, y = xy(select(rows, source_type="model", category=salt))
            axis.plot(x, y, color=color, linewidth=1.8)
        axis.set_title(title, fontsize=10)
        axis.set_xlim(0.0, 6.0)
        axis.set_ylim(0.4, 4.6)
        axis.set_xlabel(r"$\bar{m}_{salt}$ / mol kg$^{-1}$")
        axis.set_ylabel(r"$\gamma_{\pm}^{m,*}$ / -")
        axis.legend(loc="upper left", fontsize=8)
    figure.suptitle("Molality-based mean ionic activity coefficients of alkali halides in water\nat 298.15 K and 1 bar.", fontsize=11, y=0.99)
    figure.subplots_adjust(left=0.08, right=0.98, bottom=0.16, top=0.84, wspace=0.16)
    save(figure, output)


def render_figure_6(rows: list[dict[str, str]], output: Path) -> None:
    configure_style()
    panels = (("a)", "K_methanol", r"$x_{MeOH}$ / -", (0.0, 12.5)), ("b)", "Br_methanol", r"$x_{MeOH}$ / -", (0.0, 12.0)), ("c)", "Na_ethanol", r"$x_{EtOH}$ / -", (0.0, 20.0)), ("d)", "Cl_ethanol", r"$x_{EtOH}$ / -", (0.0, 30.0)))
    figure, axes = plt.subplots(2, 2, figsize=(7.6, 6.3))
    for axis, (panel, category, xlabel, ylim) in zip(axes.flat, panels, strict=True):
        x, y = xy(select(rows, source_type="model", category=category))
        axis.plot(x, y, color="black", linewidth=1.5, label="Provider (fresh installed artifact)")
        x, y = xy(select(rows, source_type="literature", category=category))
        axis.scatter(x, y, s=24, facecolor=LIGHT_GRAY, edgecolor=GRAY, linewidth=0.8, label="Literature")
        axis.set_xlim(0.0, 1.0)
        axis.set_ylim(*ylim)
        axis.set_xlabel(xlabel)
        axis.set_ylabel(r"$\Delta G_i^{trans,\infty,x}$ / kJ mol$^{-1}$")
        axis.text(0.03, 0.94, panel, transform=axis.transAxes)
    handles = [Line2D([0], [0], color="black", linewidth=1.5, label="Provider (fresh installed artifact)"), Line2D([0], [0], marker="o", linestyle="None", markerfacecolor=LIGHT_GRAY, markeredgecolor=GRAY, label="Literature")]
    figure.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.945), fontsize=9)
    figure.suptitle("Gibbs energies of transfer at infinite dilution from water to water + organic solvent systems\nat 298.15 K and 1 bar.", fontsize=11, y=0.995)
    figure.subplots_adjust(left=0.10, right=0.98, bottom=0.14, top=0.87, wspace=0.24, hspace=0.28)
    save(figure, output)


def render_figure_7(rows: list[dict[str, str]], output: Path) -> None:
    configure_style()
    figure, axis = plt.subplots(figsize=(3.4, 2.8))
    x, y = xy(select(rows, source_type="literature", category="NaBr_methanol"))
    axis.scatter(x, y, s=24, facecolor=LIGHT_GRAY, edgecolor=GRAY, linewidth=0.8, label="Literature")
    x, y = xy(select(rows, source_type="model", category="NaBr_methanol"))
    axis.plot(x, y, color="black", linewidth=1.5, label="Provider (fresh installed artifact)")
    axis.set_xlim(0.0, max(x))
    axis.set_ylim(0.0, 1.0)
    axis.set_xlabel(r"$\bar{m}_{NaBr}$ / mol kg$^{-1}$")
    axis.set_ylabel(r"$\gamma_{\pm}^{m,*}$ / -")
    axis.legend(loc="upper right", fontsize=8)
    save(figure, output)


def render_figure_8(rows: list[dict[str, str]], output: Path) -> None:
    configure_style()
    figure, axes = plt.subplots(1, 3, figsize=(11.2, 3.8))
    for axis, (panel, salt, x_max, y_max) in zip(axes, FIGURE_8_PANELS, strict=True):
        for solvent, marker, data_color, model_color in (("methanol", "o", GRAY, GRAY), ("ethanol", "s", GREEN, "black")):
            category = f"{salt}_{solvent}"
            x, y = xy(select(rows, source_type="literature", category=category))
            axis.scatter(x, y, s=24, marker=marker, facecolor="none" if solvent == "methanol" else GREEN, edgecolor=data_color, linewidth=0.9, label=f"{solvent.capitalize()} data")
            x, y = xy(select(rows, source_type="model", category=category))
            axis.plot(x, y, color=model_color, linewidth=1.5, label=f"{solvent.capitalize()} Provider (fresh installed artifact)")
        axis.set_xlim(0.0, x_max)
        axis.set_ylim(0.0, y_max)
        axis.set_title(salt, fontsize=10)
        axis.set_xlabel(r"$\bar{m}_{salt}$ / mol kg$^{-1}$")
        axis.set_ylabel(r"$\gamma_{\pm}^{m,*}$ / -")
    axes[1].legend(loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=2, fontsize=8)
    figure.suptitle(r"Molality-based salt mean ionic activity coefficients $\gamma_{\pm}^{m,*}$ in methanol and ethanol" + "\nat 298.15 K and 1 bar.", fontsize=11, y=0.995)
    figure.subplots_adjust(left=0.07, right=0.99, bottom=0.17, top=0.80, wspace=0.26)
    save(figure, output)


def render_figure_9(rows: list[dict[str, str]], output: Path) -> None:
    configure_style()
    colors = {"0.4": GREEN, "0.8": ORANGE}
    figure, axes = plt.subplots(2, 2, figsize=(9.3, 8.1))
    for axis, (panel, salt, system, organic, _organic_column, x_max) in zip(axes.flat, FIGURE_9_PANELS, strict=True):
        for target in ("0.8", "0.4"):
            category = f"{salt}_{system}_w{target}"
            x, y = xy(select(rows, source_type="model", category=category))
            axis.plot(x, y, color=FIGURE_9_LINE_COLORS[target], linewidth=2.2, zorder=5, label=fr"Provider (fresh installed artifact), $w_{{org}}^{{salt-free}}={target}$")
            x, y = xy(select(rows, source_type="literature", category=category))
            axis.scatter(x, y, s=24, marker="s" if target == "0.4" else "^", facecolor="none", edgecolor=colors[target], linewidth=1.0, zorder=6, label=fr"$w_{{org}}^{{salt-free}}={target}$")
        axis.set_xlim(0.0, x_max)
        axis.set_ylim(0.0, 1.125)
        axis.set_title(f"{salt} in aqueous {organic}", fontsize=10, pad=8)
        axis.set_xlabel(r"$\bar{m}_{salt}$ / mol kg$^{-1}$")
        axis.set_ylabel(r"$\gamma_{\pm}^{m,*}$ of salt / -")
    handles = [
        Line2D([0], [0], color="black", linewidth=2.2, label=r"Provider (fresh installed artifact), $w_{org}^{salt-free}=0.8$"),
        Line2D([0], [0], marker="^", linestyle="None", markerfacecolor="none", markeredgecolor=ORANGE, color=ORANGE, label=r"$w_{org}^{salt-free}=0.8$"),
        Line2D([0], [0], color=GREEN, linewidth=2.2, label=r"Provider (fresh installed artifact), $w_{org}^{salt-free}=0.4$"),
        Line2D([0], [0], marker="s", linestyle="None", markerfacecolor="none", markeredgecolor=GREEN, color=GREEN, label=r"$w_{org}^{salt-free}=0.4$"),
    ]
    figure.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.965), fontsize=9)
    figure.suptitle(r"Molality-based salt mean ionic activity coefficients $\gamma_{\pm}^{m,*}$ in aqueous-organic solvent mixtures", fontsize=12, y=0.995)
    figure.subplots_adjust(left=0.08, right=0.98, bottom=0.08, top=0.86, wspace=0.20, hspace=0.28)
    save(figure, output)


RENDERERS = {1: render_figure_1, 2: render_figure_2, 3: render_figure_3, 4: render_figure_4, 5: render_figure_5, 6: render_figure_6, 7: render_figure_7, 8: render_figure_8, 9: render_figure_9}


def render_figure(number: int) -> None:
    configure_style()
    output = PAPER_ROOT / "figures" / f"figure-{number:02d}" / "results" / "rendered" / f"figure-{number:02d}"
    RENDERERS[number](load_rows(number), output)

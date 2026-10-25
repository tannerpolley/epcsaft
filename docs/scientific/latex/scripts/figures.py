"""Draw the existing manuscript figures and the promoted issue-501 galleries from retained values.

Usage, from the repository root: python docs/scientific/latex/scripts/figures.py
Writes vector PDFs, the T10 LaTeX table and source records in figures/generated/.
"""

import csv
import gzip
import hashlib
import io
import json
import math
import re
import subprocess
import textwrap
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[4]
A = ROOT / "analyses"
OUT = Path(__file__).resolve().parents[1] / "figures" / "generated"
BLUE, ORANGE, GREEN, PURPLE, GOLD, SKY = (
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#E69F00",
    "#56B4E9",
)
plt.rcParams.update(
    {
        "font.family": "STIXGeneral",
        "mathtext.fontset": "stix",
        "font.size": 9,
        "axes.labelsize": 9,
        "legend.fontsize": 7.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.frameon": False,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "pdf.fonttype": 42,
    }
)


def rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def save_existing(fig, name):
    """Regenerate the four earlier PDFs without changing their recorded creation dates."""
    path = (OUT / name).relative_to(ROOT)
    original = subprocess.check_output(["git", "show", f"be78aed23:{path}"], cwd=ROOT)
    date = (
        re.search(rb"/CreationDate \(D:([^)]*)\)", original)[1]
        .decode()
        .replace("'", "")
    )
    buffer = io.BytesIO()
    fig.savefig(
        buffer,
        format="pdf",
        metadata={"CreationDate": datetime.strptime(date, "%Y%m%d%H%M%S%z")},
    )
    # Matplotlib 3.11.1 formats a negative fixed UTC offset incorrectly; preserve the source bytes.
    stamp = re.search(rb"/CreationDate \(D:[^)]*\)", original)[0]
    regenerated = re.sub(
        rb"/CreationDate \(D:[^)]*\)", lambda _: stamp, buffer.getvalue()
    )
    assert len(regenerated) == len(buffer.getvalue())
    if name != "chloride-calibration.pdf":
        assert regenerated == original, f"{name}: earlier figure changed"
    (OUT / name).write_bytes(regenerated)


def derivatives():
    """F2: (a) error/criterion of pressure-chart actions by order; (b) solved-state derivative differences."""
    actions = rows(
        A / "2026-shared-derivative-qualification/results/pressure-actions.csv"
    )
    fig, (a, b) = plt.subplots(
        1, 2, figsize=(6.6, 2.9), gridspec_kw={"width_ratios": [1, 1.5]}
    )
    selections = sorted({r["selection"] for r in actions})
    marks = {
        "residual_enthalpy": "o",
        "residual_entropy": "s",
        "residual_isobaric_heat_capacity": "^",
    }
    names = {
        "residual_enthalpy": "$h^{\\mathrm{res}}$",
        "residual_entropy": "$s^{\\mathrm{res}}$",
        "residual_isobaric_heat_capacity": "$c_p^{\\mathrm{res}}$",
    }
    for obs, mark in marks.items():
        sel = [r for r in actions if r["observable"] == obs]
        x = [
            int(r["order"]) + 0.07 * (selections.index(r["selection"]) - 1.5)
            for r in sel
        ]
        a.scatter(
            x,
            [float(r["scaled_error"]) / float(r["tolerance"]) for r in sel],
            marker=mark,
            color=BLUE,
            s=16,
            label=names[obs] + " error",
        )
        a.scatter(
            [v + 0.035 for v in x],
            [float(r["uncertainty_over_criterion"]) for r in sel],
            marker=mark,
            facecolors="none",
            edgecolors=ORANGE,
            s=16,
            linewidths=0.8,
        )
    a.scatter(
        [],
        [],
        marker="o",
        facecolors="none",
        edgecolors=ORANGE,
        s=16,
        label="reference uncertainty",
    )
    a.axhline(1.0, color="black", lw=0.8)
    a.set(
        yscale="log",
        ylim=(1e-9, 3),
        xticks=[1, 2, 3],
        xlabel="Action order",
        ylabel="Value / criterion",
    )
    a.legend(loc="lower right", handletextpad=0.2)
    a.set_title("(a) Property actions, four electrolyte records", fontsize=8.5)

    sensitivities = rows(
        A / "2026-head-to-head-comparisons/results/main-8298c602/h2h-sensitivities.csv"
    )
    for column, label, color, marker, offset in [
        ("clapeyron_vs_engine_relative", "Clapeyron.jl", BLUE, "o", -0.18),
        ("feos_vs_engine_relative", "FeOs", GREEN, "x", 0.0),
        ("engine_exact_vs_fd_relative", "Central differences", ORANGE, "+", 0.18),
    ]:
        selected = [r for r in sensitivities if r[column]]
        b.scatter(
            [abs(float(r[column])) for r in selected],
            [(0 if r["case"] == "C2-0.5" else 1) + offset for r in selected],
            color=color,
            marker=marker,
            s=20,
            label=label,
        )
    b.set(
        xscale="log",
        xlim=(1e-15, 1e-6),
        ylim=(-0.5, 1.75),
        yticks=[0, 1],
        yticklabels=["Bubble point", "Liquid--liquid"],
        xlabel="Absolute relative difference",
    )
    b.legend(loc="upper left", handletextpad=0.2)
    b.set_title("(b) First-order solved-state derivatives", fontsize=8.5)
    fig.tight_layout()
    save_existing(fig, "derivatives.pdf")


def electrolyte_lle():
    """F3: Clapeyron.jl minus Engine phase mole fractions for the eight Khudaida feeds."""
    data = rows(A / "2026-head-to-head-comparisons/results/h2h-khudaida-lle.csv")
    species = ["water", "ethanol", "isobutanol", "sodium-cation", "chloride-anion"]
    label = {
        "water": "water",
        "ethanol": "ethanol",
        "isobutanol": "isobutanol",
        "sodium-cation": "Na$^+$",
        "chloride-anion": "Cl$^-$",
    }
    colors = dict(zip(species, [BLUE, ORANGE, GREEN, PURPLE, GOLD]))
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.6, 2.7))
    for ax, variant, cols, title in (
        (
            a,
            "constant-permittivity",
            [
                ("clapeyron_default_init_minus_engine", "o", 0.0),
                ("clapeyron_engine_guesses_minus_engine", "s", 0.18),
            ],
            "(a) Matched model, constant $\\epsilon_r$",
        ),
        (
            b,
            "matched",
            [("clapeyron_engine_guesses_minus_engine", "s", 0.18)],
            "(b) Composition-dependent $\\epsilon_r$",
        ),
    ):
        for col, mark, shift in cols:
            for s in species:
                for phase, fill in (("water-rich", True), ("organic", False)):
                    sel = [
                        r
                        for r in data
                        if r["variant"] == variant
                        and r["quantity"] == s
                        and r["phase"] == phase
                    ]
                    x = [
                        int(r["point"]) + shift + (0.06 if phase == "organic" else 0)
                        for r in sel
                    ]
                    y = [abs(float(r[col])) for r in sel]
                    ax.scatter(
                        x,
                        y,
                        marker=mark,
                        s=12,
                        linewidths=0.8,
                        facecolors=colors[s] if fill else "none",
                        edgecolors=colors[s],
                    )
        ax.set(yscale="log", xlabel="Feed", xticks=range(1, 9))
        ax.set_title(title, fontsize=8.5)
    a.set(ylabel="|Clapeyron.jl $-$ package| mole fraction", ylim=(1e-13, 1e-7))
    b.set(ylim=(1e-6, 3e-2))
    for s in species:
        b.scatter([], [], marker="s", color=colors[s], s=12, label=label[s])
    a.scatter(
        [], [], marker="o", color="gray", s=12, label="Clapeyron.jl default start"
    )
    a.scatter([], [], marker="s", color="gray", s=12, label="from package phases")
    a.scatter(
        [],
        [],
        marker="s",
        facecolors="none",
        edgecolors="gray",
        s=12,
        label="open: organic liquid",
    )
    a.legend(loc="lower left", fontsize=6.5)
    b.legend(loc="lower right", ncol=2, fontsize=6.5)
    fig.tight_layout()
    save_existing(fig, "electrolyte-lle.pdf")


def acetic_acid_ph(axes=None, data=None):
    """F4: aqueous pH of acetic acid + water + toluene at 298.15 K against acid mole fraction."""
    if data is None:
        # Keep the existing PDF's result identity; the new gallery passes promoted data explicitly.
        paths = ["sweep-aq-literal.csv", "comparison.csv"]
        data = [
            list(
                csv.DictReader(
                    io.StringIO(
                        subprocess.check_output(
                            [
                                "git",
                                "show",
                                "12d0ba46b2d4db5731d3c827bea36411f9e5f5f6:"
                                + "analyses/2026-ascani-carboxylic-acid-lle/results/"
                                + p,
                            ],
                            cwd=ROOT,
                            text=True,
                        )
                    )
                )
            )
            for p in paths
        ]
    sweep, comp = data
    sweep = sorted(sweep, key=lambda r: float(r["x_acid_total_aq"]))
    own_figure = axes is None
    if own_figure:
        fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.7))
    a, b = axes
    a.plot(
        [float(r["x_acid_total_aq"]) for r in sweep],
        [float(r["pH_aq"]) for r in sweep],
        color=BLUE,
        lw=1.4,
        label="this package, published parameters",
    )
    xs = [float(r["x_acetic_acid_aq"]) for r in comp]
    a.scatter(
        xs,
        [float(r["pH_esi_model"]) for r in comp],
        marker="x",
        color=GREEN,
        s=22,
        label="published ePC-SAFT model",
    )
    a.scatter(
        xs,
        [float(r["pH_exp"]) for r in comp],
        marker="o",
        facecolors="none",
        edgecolors=ORANGE,
        s=18,
        label="experiment",
    )
    a.set(xlabel="Acetic acid mole fraction, aqueous liquid", ylabel="Aqueous pH")
    a.legend(loc="upper right")
    a.set_title("(a) pH against acid content", fontsize=8.5)
    b.axhspan(-0.005, 0.005, color=SKY, alpha=0.25, lw=0, label="$\\pm0.005$ criterion")
    b.scatter(
        xs,
        [float(r["pH_aq-literal_minus_esi"]) for r in comp],
        marker="o",
        color=BLUE,
        s=16,
        label="package $-$ published model",
    )
    b.axhline(0, color="black", lw=0.6)
    b.set(xlabel="Acetic acid mole fraction, aqueous liquid", ylabel="pH difference")
    b.legend(loc="lower left")
    b.set_title("(b) Difference from the published model", fontsize=8.5)
    if own_figure:
        fig.tight_layout()
        save_existing(fig, "acetic-acid-ph.pdf")


def chloride_calibration():
    """F5: signed log10(model/evaluated) of LiCl, NaCl and KCl activity and osmotic coefficients at 298.15 K."""
    fit = [
        r
        for r in rows(
            A / "2026-figiel-divalent-chlorides/results/stage1-minimax-rows.csv"
        )
        if r["set"] == "T" and r["start"] == "0"
    ]
    source = rows(
        A / "2026-figiel-divalent-chlorides/results/figiel-1to1-comparison.csv"
    )
    colors = {"LiCl": BLUE, "NaCl": ORANGE, "KCl": GREEN}
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.7), sharey=True)
    for ax, quantity, col, title in (
        (axes[0], "gamma", "log10_gamma_ratio", "(e) $\\gamma_\\pm$"),
        (axes[1], "osmotic", "log10_osmotic_ratio", "(f) $\\Phi$"),
    ):
        ax.axhspan(-0.02, 0.02, color=SKY, alpha=0.25, lw=0)
        for salt, color in colors.items():
            sel = [r for r in fit if r["group"] == f"T {salt} {quantity}"]
            ax.scatter(
                [float(r["abscissa"]) for r in sel],
                [float(r["error"]) for r in sel],
                marker="o",
                color=color,
                s=12,
                label=f"{salt}, fitted",
            )
            src = [r for r in source if r["salt"] == salt]
            ax.scatter(
                [float(r["molality_mol_kg"]) for r in src],
                [float(r[col]) for r in src],
                marker="x",
                color=color,
                s=14,
                linewidths=0.8,
                label=f"{salt}, published",
            )
        ax.axhline(0, color="black", lw=0.6)
        ax.set(xscale="log", xlabel="Molality / (mol kg$^{-1}$)")
        ax.set_title(title, fontsize=8.5)
    axes[0].set(ylabel="$\\log_{10}$(model / evaluated)")
    axes[1].legend(loc="lower left", ncol=2, fontsize=6.5)
    fig.tight_layout()
    save_existing(fig, "chloride-calibration.pdf")


# Investigator-selected identities, issue 501, owner comment 6373385229.
PROMOTED = {
    "cost": (
        "2026-release-computational-cost",
        "409e8e9b52052b320c645281832449ed90bb856a",
    ),
    "gross01": ("reproduction/2001-gross", "e75c60ad"),
    "gross02": ("reproduction/2002-gross", "e75c60ad"),
    "held08": ("reproduction/2008-held", "e75c60ad"),
    "opensource": ("2026-open-source-equilibrium-comparison", "e75c60ad"),
    "settler": ("2026-khudaida-idaes-settler", "e75c60ad"),
    "dipole": ("2026-polar-dipole", "83f95810"),
    "propane_fit": ("2026-pure-saturation-regression", "83f95810"),
    "dielectric": ("2026-dielectric-alternatives", "83f95810"),
    "wangler": ("reproduction/2018-wangler", "83f95810"),
    "cleeton": ("reproduction/2020-cleeton", "83f95810"),
    "glos": ("reproduction/2004-glos", "5912b74b"),
    "quadrupole": ("2026-polar-quadrupole", "b11d69d6"),
    "may": ("reproduction/2015-may", "0c44f995"),
    "schick": ("reproduction/2023-schick", "04f68bbf"),
    "density": ("2026-density-root-benchmark", "0aad8f08"),
    "discovery": ("2026-tpd-phase-discovery", "32044004"),
    "haslam": ("reproduction/2020-haslam", "676427c9"),
    "held14": ("reproduction/2014-held", "2ecfa3a9"),
    "refit": ("2026-chloride-brine-refit", "6e470be2"),
    "mixed": ("2026-mixed-chloride-water-activity", "c3817d27"),
    "chlorides": ("2026-held-2014-chlorides", "0aad8f08"),
    "esteso": ("reproduction/1989-esteso", "662dea9f"),
    "figiel": ("reproduction/2025-figiel", "9bfc702f"),
    "bulow20": ("reproduction/2020-bulow", "50792d0d"),
    "bulow21": ("reproduction/2021-bulow", "50792d0d"),
    "held12": ("reproduction/2012-held", "55cd00d0"),
    "khudaida": ("reproduction/2026-khudaida", "bb2356be"),
    "acid": ("2026-ascani-carboxylic-acid-lle", "2ecfa3a9"),
    "ascani": ("reproduction/2023-ascani", "37de3ed5"),
    "chase": ("reproduction/1998-chase", "3f6b115a"),
    "belov": ("reproduction/2023-belov", "5839e6d2"),
    "algorithms": ("2026-algorithm-illustrations", "d9469473"),
    "hda": ("2026-hda-two-flash", "816320bc"),
}
CAL, PRED, REPRO, NV = (
    "calibration",
    "prediction",
    "published-model reproduction",
    "numerical verification",
)
COLORS = [BLUE, ORANGE, GREEN, PURPLE, GOLD, SKY, "#333333", "#777777"]
MARKS = ["o", "s", "^", "v", "D", "P", "<", ">"]
TABLE, FILES, FIGURES, PLOTTED = [], {}, {}, []
USED = set()
PDFS = {}


def retained(study, path, absolute=False):
    """Read Git objects so a later leaf revision cannot change a promoted figure."""
    leaf, commit = PROMOTED[study]
    path = path if absolute else f"analyses/{leaf}/{path}"
    key = f"{commit}:{path}"
    if key not in FILES:
        data = subprocess.check_output(["git", "show", key], cwd=ROOT)
        FILES[key] = {
            "commit": subprocess.check_output(
                ["git", "rev-parse", commit], cwd=ROOT, text=True
            ).strip(),
            "path": path,
            "sha256": hashlib.sha256(data).hexdigest(),
            "data": (gzip.decompress(data) if path.endswith(".gz") else data).decode(),
        }
    USED.add(key)
    return FILES[key]["data"]


def csv_at(study, path, absolute=False):
    return list(csv.DictReader(io.StringIO(retained(study, path, absolute))))


def json_at(study, path):
    return json.loads(retained(study, path))


def values(data, column):
    return [float(r[column]) for r in data]


def pairs(data, model, reference):
    return [
        (float(r[model]), float(r[reference]))
        for r in data
        if r[model] and r[reference] and math.isfinite(float(r[model]))
    ]


def metric(operands, name="AARD / %"):
    assert operands and all(math.isfinite(a) and math.isfinite(b) for a, b in operands)
    errors = [a - b for a, b in operands]
    if name == "AARD / %":
        return 100 * sum(abs(a / b - 1) for a, b in operands) / len(operands)
    if name.startswith("RMS"):
        return math.sqrt(sum(e * e for e in errors) / len(errors))
    if name.startswith("MAE"):
        return sum(abs(e) for e in errors) / len(errors)
    if name.startswith("max |log10"):
        return max(abs(math.log10(a / b)) for a, b in operands)
    if name.startswith("max |ln"):
        return max(abs(math.log(a) - math.log(b)) for a, b in operands)
    if name.startswith("max"):
        return max(abs(e) for e in errors)
    raise ValueError(name)


def panel(ax, label, title, evidence, xlabel, ylabel):
    ax.set_title("")
    ax.set_label(label)
    width = 36 if ax.get_subplotspec().get_gridspec().ncols == 3 else 65
    number = ax.figure.axes.index(ax) + (1 if label == "F4a" else 0)
    if label.startswith("F18"):
        number += 2
    if label.startswith("F8c-dimethyl"):
        number += 4
    visible = f"({chr(97 + number)})"
    heading = textwrap.fill(
        f"{visible}  {title}", width, break_long_words=False, break_on_hyphens=False
    )
    evidence = textwrap.fill(
        evidence, width, break_long_words=False, break_on_hyphens=False
    )
    ax.set_title(f"{heading}\n{evidence}", fontsize=7.6, loc="left", pad=5)
    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.tick_params(labelsize=7)
    ax.xaxis.label.set_size(8)
    ax.yaxis.label.set_size(8)


def canvas(count, columns=2):
    rows_count = math.ceil(count / columns)
    fig, axes = plt.subplots(
        rows_count, columns, figsize=(7.1, 2.65 * rows_count), squeeze=False
    )
    flat = list(axes.flat)
    for ax in flat[count:]:
        ax.remove()
    return fig, flat[:count]


def comparison(
    ax,
    data,
    x,
    model,
    reference,
    label="",
    color=BLUE,
    mark="o",
    continuous=False,
    uncertainty=None,
    scale=1,
):
    """Only retained sweeps have segments; observation-coordinate model values are crosses."""
    data = sorted(data, key=lambda r: float(r[x]))
    xs = values(data, x)
    ref = [float(r[reference]) / scale for r in data]
    mod = [float(r[model]) / scale if r[model] else math.nan for r in data]
    ax.errorbar(
        xs,
        ref,
        yerr=[float(r[uncertainty]) / scale for r in data] if uncertainty else None,
        fmt=mark,
        ms=2.8,
        mfc="none",
        mec=color,
        ecolor=color,
        elinewidth=0.5,
        mew=0.7,
        label=label or "Reference",
        alpha=0.85,
    )
    ax.plot(
        xs,
        mod,
        "-" if continuous else "x",
        color=color,
        lw=1,
        ms=3,
        mew=0.7,
        label="Model" if not label else None,
    )
    for r in data:
        PLOTTED.append(
            {
                "panel": ax.get_label(),
                "series": label,
                "x": float(r[x]),
                "reference": float(r[reference]) / scale,
                "model": float(r[model]) / scale if r[model] else None,
            }
        )


def record(
    label,
    system,
    quantity,
    temperatures,
    pressures,
    n,
    parameter,
    fitted,
    evidence,
    operands,
    metric_name="AARD / %",
    note="",
    fluid_class="",
    property_class="",
):
    TABLE.append(
        {
            "panel": label,
            "system": system,
            "property": quantity,
            "T_K": sorted(set(temperatures)),
            "P_Pa": sorted(set(pressures)),
            "points": n,
            "comparisons": len(operands),
            "parameter_source": parameter,
            "fitted": fitted,
            "evidence": evidence,
            "metric": metric_name,
            "value": metric(operands, metric_name),
            "operands": operands,
            "note": note,
            "fluid_class": fluid_class,
            "property_class": property_class,
        }
    )


def finish(fig, name, notes):
    bottom = 0.09 if fig.legends else 0.025
    fig.tight_layout(pad=0.8, h_pad=1.2, w_pad=0.9, rect=(0, bottom, 1, 1))
    footer = "Open symbols: measurements or references; crosses: discrete model evaluations; lines: continuous model calculations."
    if name.startswith("F15-"):
        footer = "Top: measurements and discrete model values. Middle: relative deviations. Bottom: chloride calibration errors."
    elif name.startswith("F6-"):
        footer = "Points: median paired block ratios; intervals: full observed range; vertical line: equal time."
    elif name.startswith("F18-"):
        footer = (
            "Open symbols: model bubble pressure at each measured state; line: parity."
        )
    fig.text(0.5, 0.005, footer, ha="center", va="bottom", fontsize=6.5)
    if name not in PDFS:
        PDFS[name] = PdfPages(
            OUT / name, metadata={"CreationDate": None, "ModDate": None}
        )
    PDFS[name].savefig(fig)
    plt.close(fig)
    panels = [r for r in TABLE if "figure" not in r]
    for r in panels:
        r.update(figure=name, sources=sorted(USED))
    info = FIGURES.setdefault(name, {"panels": [], "sources": [], "notes": ""})
    info["panels"].extend(r["panel"] for r in panels)
    info["sources"] = sorted(set(info["sources"]) | USED)
    if notes not in info["notes"]:
        info["notes"] += (" " if info["notes"] else "") + notes
    USED.clear()


def computational_cost():
    """F6: value-matched operation cost and marginal actions on solved states."""
    data = csv_at("cost", "results/ratios.csv")
    operations = [
        ("property-bundle", "Fixed-density bundle"),
        ("pressure-state", "Density closure"),
        ("C1-300", "Propane saturation"),
        ("C2-0.5", "Ethanol–butane bubble"),
        ("C4-298.15", "Water–pentanol split"),
        ("explicit-discovery", "Methane–ethane discovery"),
    ]
    actions = [
        ("C2-0.5/actions-on-solved-state", "Ethanol–butane bubble"),
        ("C4-298.15/actions-on-solved-state", "Water–pentanol split"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.1))
    for ax, cases, comparisons in [
        (
            axes[0],
            operations,
            [("FeOs", BLUE, "o", -0.12), ("Clapeyron", ORANGE, "x", 0.12)],
        ),
        (axes[1], actions, [("package central differences", GREEN, "o", 0)]),
    ]:
        for comparator, color, marker, offset in comparisons:
            selected = [
                next(
                    r
                    for r in data
                    if r["operation"] == operation and r["comparator"] == comparator
                )
                for operation, _ in cases
            ]
            medians = values(selected, "package_over_comparator")
            ax.errorbar(
                medians,
                [i + offset for i in range(len(cases))],
                xerr=[
                    [
                        m - float(r["minimum_paired_ratio"])
                        for m, r in zip(medians, selected)
                    ],
                    [
                        float(r["maximum_paired_ratio"]) - m
                        for m, r in zip(medians, selected)
                    ],
                ],
                fmt=marker,
                color=color,
                ms=3.5,
                capsize=2,
                lw=0.8,
                label="Clapeyron.jl" if comparator == "Clapeyron" else comparator,
            )
        ax.axvline(1, color="black", lw=0.8)
        ax.set(
            xscale="log",
            yticks=range(len(cases)),
            yticklabels=[label for _, label in cases],
            ylim=(len(cases) - 0.5, -0.6),
        )
        ax.set_xlabel(
            "Package / comparator time"
            if ax is axes[0]
            else "Exact action / central-difference time"
        )
    axes[0].legend(loc="upper right", fontsize=7)
    axes[0].set_title("(a) Value-matched operations", loc="left", fontsize=8.5)
    axes[1].set_title(
        "(b) Actions on solved states\nCentral solve excluded", loc="left", fontsize=8.5
    )
    finish(
        fig,
        "F6-computational-cost.pdf",
        "Warm serial measurements on one host; ranges are observed spread, not confidence intervals. Density closure against Clapeyron.jl is a tie. Actions exclude the central solve and imply no total solve speedup. Electrolyte timings are package-only and excluded.",
    )
    FIGURES["F6-computational-cost.pdf"]["panels"] = ["F6a", "F6b"]


def pure_fluids():
    fig, axes = canvas(3)
    data = csv_at("glos", "results/propane-saturation.csv")
    for ax, quantity, phase, label, ylabel, scale in zip(
        axes,
        ["pressure", "liquid_density", "vapor_density"],
        ["", "L", "V"],
        ["F8a-P", "F8a-L", "F8a-V"],
        [
            "Saturation pressure / MPa",
            r"Liquid density / kg m$^{-3}$",
            r"Vapor density / kg m$^{-3}$",
        ],
        [1e6, 1, 1],
    ):
        sel = [r for r in data if r["quantity"] == quantity and r["phase"] == phase]
        assert sel, (quantity, phase)
        panel(ax, label, "Propane / Glos", PRED, "Temperature / K", ylabel)
        comparison(
            ax,
            sel,
            "T_K",
            "predicted",
            "observed",
            uncertainty="uncertainty",
            scale=scale,
        )
        if quantity == "pressure":
            ax.set_yscale("log")
        ax.legend(fontsize=6)
        record(
            label,
            "Propane",
            ylabel.split(" / ")[0],
            values(sel, "T_K"),
            values(sel, "observed") if quantity == "pressure" else [],
            len(sel),
            "Gross-Sadowski 2001 Table 2",
            "no",
            PRED,
            pairs(sel, "predicted", "observed"),
            note="Expanded measurement uncertainties; no additional fit. Saturation conditions; density-row pressures are not separately retained.",
            fluid_class="Pure fluids",
            property_class="Saturation" if quantity == "pressure" else "Density",
        )
    finish(
        fig,
        "F8-pure-fluids.pdf",
        "Glos propane predictions retain expanded measurement uncertainties.",
    )
    gross_pure()
    dipole_pure()
    fig, axes = canvas(2)
    q = csv_at("quadrupole", "results/engine-rows.csv")
    for ax, species in zip(axes, ["carbon-dioxide", "benzene"]):
        sel = [
            dict(r, T=r["case"].split()[-2])
            for r in q
            if r["kind"] == "p_sat" and r["case"].startswith(species)
        ]
        label = "F8d-CO2" if species == "carbon-dioxide" else "F8d-Bz"
        panel(
            ax,
            label,
            species.replace("-", " "),
            REPRO,
            "Temperature / K",
            "Saturation pressure / MPa",
        )
        comparison(
            ax,
            sel,
            "T",
            "engine",
            "reference",
            label="Reference correlation",
            scale=1e6,
        )
        if species == "carbon-dioxide":
            nonpolar = [
                dict(r, T=r["case"].split()[-2])
                for r in q
                if r["kind"] == "p_sat_nonpolar_gs2001"
                and r["case"].startswith(species)
            ]
            ax.plot(
                values(nonpolar, "T"),
                [float(r["engine"]) / 1e6 for r in nonpolar],
                "+",
                color=ORANGE,
                ms=4,
                label="Nonpolar parameters",
            )
        ax.legend(fontsize=6)
        record(
            label,
            species,
            "Saturation pressure",
            values(sel, "T"),
            values(sel, "reference"),
            len(sel),
            "Gross 2005 Table A1",
            "original rows unknown",
            REPRO,
            pairs(sel, "engine", "reference"),
            note="Evaluated correlation, not raw measurements; benzene DIPPR coefficient attribution unverified.",
            fluid_class="Pure fluids",
            property_class="Saturation",
        )
    finish(
        fig,
        "F8-pure-fluids.pdf",
        "Glos: all observations and expanded uncertainties. Quadrupole: reference correlations; separate nonpolar parameters.",
    )


def neutral_mixtures():
    fig, axes = canvas(2)
    data = csv_at("may", "results/comparison.csv")
    for ax, phase, observed in zip(
        axes, ["liquid", "vapor"], ["x_methane", "y_methane"]
    ):
        label = "F9a-" + phase
        panel(
            ax,
            label,
            "Methane/ethane",
            PRED,
            "Pressure / MPa",
            f"Methane fraction, {phase}",
        )
        sel = [dict(r, P_MPa=float(r["P_Pa"]) / 1e6) for r in data]
        comparison(
            ax,
            sel,
            "P_MPa",
            "model_" + observed,
            observed,
            uncertainty="uc_" + observed,
        )
        missing = [r for r in data if not r["model_" + observed]]
        for r in missing:
            ax.annotate(
                "undetermined",
                (float(r["P_Pa"]) / 1e6, float(r[observed])),
                xytext=(-45, -15),
                textcoords="offset points",
                fontsize=6,
                arrowprops={"arrowstyle": "-", "lw": 0.5},
            )
        ax.legend(fontsize=6)
        record(
            label,
            "Methane/ethane",
            f"{phase} methane mole fraction",
            values(data, "T_K"),
            values(data, "P_Pa"),
            len(data),
            "Gross-Sadowski 2001; kij = 0",
            "no",
            PRED,
            pairs(data, "model_" + observed, observed),
            "MAE / mole fraction",
            note="Points span the retained temperatures; no isothermal line is implied. Combined composition uncertainty only; failed comparisons retained.",
            fluid_class="Neutral mixtures",
            property_class="VLE",
        )
    finish(
        fig,
        "F9-neutral-mixtures.pdf",
        "May methane/ethane predictions retain the undetermined state.",
    )
    dipole_binary()
    fig, axes = canvas(1, 1)
    sel = csv_at("schick", "results/data/comparison.csv")
    ax = axes[0]
    panel(ax, "F9d", "CO2/THF", REPRO, "Pressure / bar", "Liquid CO2 mole fraction")
    for j, temp in enumerate(sorted({r["temperature_k"] for r in sel})):
        part = [r for r in sel if r["temperature_k"] == temp]
        comparison(
            ax,
            part,
            "pressure_bar",
            "engine_x",
            "experimental_liquid_co2_mole_fraction",
            f"{float(temp):g} K",
            COLORS[j],
            MARKS[j],
        )
        ax.plot(
            values(part, "pressure_bar"),
            values(part, "paper_model_liquid_co2_mole_fraction"),
            "+",
            color=COLORS[j],
            ms=4,
        )
    ax.plot([], [], "+", color="black", label="Published model")
    ax.legend(fontsize=6, ncol=2)
    record(
        "F9d",
        "CO2/THF",
        "Liquid CO2 mole fraction",
        values(sel, "temperature_k"),
        [p * 1e5 for p in values(sel, "pressure_bar")],
        len(sel),
        "Schick 2023 Tables 1,6/S3",
        "inherited source fits",
        REPRO,
        pairs(sel, "engine_x", "experimental_liquid_co2_mole_fraction"),
        note="Metric against measurements; printed model values shown separately. Source-model non-reproduction has an undetermined cause; measurement uncertainty absent.",
        fluid_class="Neutral mixtures",
        property_class="VLE",
    )
    finish(
        fig,
        "F9-neutral-mixtures.pdf",
        "May: retain all source states, including the undetermined state. Schick: plus symbols are the published model; discrepancies remain.",
    )


def aqueous_electrolytes():
    fig, axes = canvas(2)
    data = csv_at("haslam", "results/data/haslam-2020-table8.csv")
    cond = json_at("haslam", "results/data/haslam-2020-table8.json")
    for ax, obs, label, ylabel in zip(
        axes,
        ["gamma_pm_m", "osmotic_coefficient"],
        ["F11a-gamma", "F11a-Phi"],
        [r"$\gamma_\pm$", r"$\Phi$"],
    ):
        sel = [r for r in data if r["observable"] == obs]
        assert sel, {r["observable"] for r in data}
        panel(
            ax,
            label,
            "Aqueous 1:1 salts",
            CAL,
            r"Molality / mol kg$^{-1}$ water",
            ylabel,
        )
        for j, salt in enumerate(sorted({r["salt"] for r in sel if r["model_value"]})):
            comparison(
                ax,
                [r for r in sel if r["salt"] == salt],
                "molality_mol_kg",
                "model_value",
                "literature_value",
                salt,
                COLORS[j],
                MARKS[j],
            )
        ax.set_xscale("log")
        ax.legend(fontsize=5.7, ncol=4, loc="upper left", bbox_to_anchor=(0, 0.9))
        absent = sum(not r["model_value"] for r in sel)
        ax.text(
            0.02,
            0.96,
            f"LiI: {absent} values unavailable",
            transform=ax.transAxes,
            va="top",
            fontsize=6,
        )
        record(
            label,
            "Li/Na/K chlorides, bromides, iodides",
            obs,
            [cond["temperature_K"]],
            [cond["pressure_Pa"]],
            len(sel),
            "Figiel 2025 inherited parameters",
            "partly; exact overlap incomplete",
            CAL,
            pairs(sel, "model_value", "literature_value"),
            note="Hamer-Wu evaluated coefficients, substantial calibration overlap; LiI unavailable; packet pressure extrapolation unbounded.",
            fluid_class="Aqueous electrolytes",
            property_class="Activity / osmotic",
        )
    finish(
        fig,
        "F11-aqueous-electrolytes.pdf",
        "Haslam: Hamer-Wu evaluated coefficients retain calibration overlap and unavailable LiI values.",
    )
    held_2008()
    fig, axes = canvas(2)
    line = [
        r
        for r in csv_at("held14", "results/figure-line-comparison.csv")
        if r["figure"] in ("figure-1", "figure-2")
    ]
    sweep = csv_at("held14", "results/figure-curves.csv")
    cold = csv_at("held14", "results/figure-2-cold-markers.csv")
    identity = json_at("held14", "results/run-identity.json")
    pressure = float(re.search(r"([\d.]+) bar", identity["conditions"])[1]) * 1e5
    for ax, temp in zip(axes, sorted({r["temperature_K"] for r in line}, reverse=True)):
        sel = [r for r in line if r["temperature_K"] == temp]
        label = "F11c-" + str(round(float(temp)))
        panel(
            ax,
            label,
            f"Held {float(temp):g} K",
            REPRO,
            r"Molality / mol kg$^{-1}$ water",
            r"$\Phi$",
        )
        for j, series in enumerate(sorted({r["series"] for r in sel})):
            part = [r for r in sel if r["series"] == series]
            ax.plot(
                values(part, "molality_mol_kg"),
                values(part, "osmotic_paper_line"),
                "o",
                color=COLORS[j],
                mfc="none",
                ms=1.6,
                mew=0.4,
            )
            trace = sorted(
                [r for r in sweep if r["series"] == series],
                key=lambda r: float(r["molality_mol_kg"]),
            )
            ax.plot(
                values(trace, "molality_mol_kg"),
                values(trace, "osmotic_engine"),
                color=COLORS[j],
                ls="-" if j % 2 == 0 else "--",
                lw=0.9,
                label=series,
            )
        if float(temp) == min(values(line, "temperature_K")):
            for j, salt in enumerate(sorted({r["salt"] for r in cold})):
                part = [r for r in cold if r["salt"] == salt]
                ax.scatter(
                    values(part, "molality_mol_kg"),
                    values(part, "osmotic_digitized"),
                    marker="s",
                    s=10,
                    facecolors="none",
                    edgecolors="black",
                    linewidths=0.5,
                    label=f"{salt} measured",
                )
        ax.legend(fontsize=5.5, ncol=2)
        record(
            label,
            "NaCl/KCl/KBr in water",
            "Osmotic coefficient, model-line reproduction",
            [float(temp)],
            [pressure],
            len(sel),
            "Held 2014 strategies 1,2; 0.88 sigma",
            "source calibration",
            REPRO,
            pairs(sel, "osmotic_engine", "osmotic_paper_line"),
            "max absolute Phi difference",
            note="Digitized model-line samples, not independent observations; cold measured symbols shown separately. Inferred screening convention and corrected water diameter.",
            fluid_class="Aqueous electrolytes",
            property_class="Activity / osmotic",
        )
    finish(
        fig,
        "F11-aqueous-electrolytes.pdf",
        "Haslam selects Hamer-Wu evaluated data, not Haslam SAFT-gamma predictions. Held panels are bounded halide model-line reproductions with cold experimental symbols.",
    )


def mixed_brines():
    fig, axes = canvas(7)
    train = csv_at("refit", "results/stage-a-calcium/training-residuals-F.csv")
    hold = [
        r
        for r in csv_at("refit", "results/stage-a-calcium/nominal-residuals-F.csv")
        if r["role"] == "mixture" and r["source"] != "rard-miller-1987"
    ]
    cond = json_at("refit", "inputs/fit-declaration.json")["conditions"]
    limits = csv_at("refit", "results/source-agreement-limits.csv")
    limits = [r for r in limits if r["realization_id"] == "r0000"]
    train = [
        dict(r, I=float(r["m_NaCl"]) + float(r["m_KCl"]) + 3 * float(r["m_divalent"]))
        for r in train
    ]
    for ax, sel, label, evidence, prop, x in [
        (axes[0], train, "F12a-fit", CAL, "Phi", "I"),
        (
            axes[1],
            [r for r in hold if r["quantity"] == "Phi"],
            "F12a-Phi",
            PRED,
            "Phi",
            "I_mol_kg",
        ),
        (
            axes[2],
            [r for r in hold if r["quantity"] == "aw"],
            "F12a-aw",
            PRED,
            "aw",
            "I_mol_kg",
        ),
    ]:
        panel(
            ax,
            label,
            "Joint chloride refit",
            evidence,
            r"Ionic strength / mol kg$^{-1}$ water",
            r"$\Delta\Phi$ (model − measured)" if prop == "Phi" else r"$\Delta a_w$ (model − measured)",
        )
        for j, source in enumerate(sorted({r["source"] for r in sel})):
            part = [r for r in sel if r["source"] == source]
            ax.plot(values(part, x), [float(r["model"]) - float(r["observed"]) for r in part],
                    MARKS[j], color=COLORS[j], ms=3, label=source.split("-")[0].capitalize())
        ax.axhline(0, color="black", lw=0.5)
        for lo, hi, bucket in [(0, 3, "I<=3"), (3, max(values(sel, x)), "I>3")]:
            bound = 0.01 if prop == "Phi" else float(next(r for r in limits if r["I_bucket"] == bucket)["limit_RMS"])
            ax.fill_between([lo, hi], -bound, bound, color="grey", alpha=0.15,
                            label="RMS comparison limit" if lo == 0 else None)
        ax.legend(fontsize=5.7, ncol=2)
        record(
            label,
            "Na/Mg/Ca/K chloride-water",
            prop,
            [cond["temperature_K"]],
            [cond["pressure_Pa"]],
            len(sel),
            "Joint F refit, inherited Figiel form",
            "yes" if evidence == CAL else "no",
            evidence,
            pairs(sel, "model", "observed"),
            "RMS / dimensionless",
            note="Separate property populations. Held-out source limit and the failing realization retained; no Na/Mg/Ca combined feed.",
            fluid_class="Mixed brines",
            property_class="Activity / osmotic",
        )
    baseline = [
        r
        for r in csv_at("mixed", "results/nominal-residuals.csv")
        if r["role"] == "mixture"
    ]
    cond = json_at("mixed", "inputs/comparison-specification.json")["conditions"]
    for ax, prop in zip(axes[3:5], ["Phi", "aw"]):
        label = "F12b-" + prop
        panel(
            ax,
            label,
            "Unchanged parameters",
            PRED,
            r"Ionic strength / mol kg$^{-1}$ water",
            r"$\Delta\Phi$ (model − measured)" if prop == "Phi" else r"$\Delta a_w$ (model − measured)",
        )
        sel = [r for r in baseline if r["quantity"] == prop]
        sources = sorted({r["source"] for r in sel})
        for j, source in enumerate(sources):
            for declaration, mark in [("Held", "x"), ("Figiel", "+")]:
                part = [
                    r
                    for r in sel
                    if r["source"] == source and r["declaration"] == declaration
                ]
                ax.plot(
                    values(part, "I_mol_kg"),
                    [float(r["model"]) - float(r["observed"]) for r in part],
                    mark,
                    ms=3,
                    color=COLORS[j],
                    label=source.split("-")[0].capitalize() + " " + declaration,
                )
        bound = 0.01 if prop == "Phi" else 0.0013
        ax.axhspan(-bound, bound, color="grey", alpha=0.15, label="Original RMS limit")
        ax.axhline(0, color="black", lw=0.5)
        ax.legend(fontsize=5.4, ncol=2)
        record(
            label,
            "NaCl/MgCl2; NaCl/CaCl2 in water",
            prop,
            [cond["temperature_K"]],
            [cond["pressure_Pa"]],
            len(sel) // 2,
            "Unchanged Held / Figiel declarations",
            "no",
            PRED,
            pairs(sel, "model", "observed"),
            "RMS / dimensionless",
            note="Metric pools both declarations, not extra observations. Both declarations fail every nominal source/strength bucket; failures are not extrapolated to the full EOS.",
            fluid_class="Mixed brines",
            property_class="Activity / osmotic",
        )
    sigma = csv_at("chlorides", "results/sigma-replay/evaluated-comparison.csv")
    shorter = csv_at("chlorides", "results/dh-088/evaluated-comparison.csv")
    cc = csv_at("chlorides", "results/sigma-replay/model-curves.csv")
    for ax, prop in zip(axes[5:], ["gamma", "osmotic"]):
        label = "F12c-" + prop
        panel(
            ax,
            label,
            "Divalent and 1:1 limits",
            REPRO,
            r"Molality / mol kg$^{-1}$ water",
            r"$\gamma_\pm$" if prop == "gamma" else r"$\Phi$",
        )
        operands = []
        n = 0
        for j, salt in enumerate(["LiCl", "NaCl", "MgCl2", "CaCl2"]):
            for data, mark, suffix in [
                (sigma, "x", r"$\sigma$"),
                (shorter, "+", r"$0.88\sigma$"),
            ]:
                mode = "held2014" if data is sigma else "held2014_dh088"
                part = [r for r in data if r["mode"] == mode and r["salt"] == salt]
                assert part
                if data is sigma:
                    n += len(part)
                    ax.plot(
                        values(part, "molality_mol_kg"),
                        values(part, prop + "_evaluated"),
                        "o",
                        color=COLORS[j],
                        mfc="none",
                        ms=2,
                        mew=0.6,
                    )
                ax.plot(
                    values(part, "molality_mol_kg"),
                    values(part, prop + "_engine"),
                    mark,
                    color=COLORS[j],
                    ms=2.7,
                    label=salt + " " + suffix,
                )
                operands.extend(pairs(part, prop + "_engine", prop + "_evaluated"))
        ax.legend(fontsize=5.3, ncol=2)
        record(
            label,
            "LiCl/NaCl/MgCl2/CaCl2-water",
            prop,
            values(cc, "temperature_K"),
            values(cc, "pressure_Pa"),
            n,
            "Held 2014, screening sigma / 0.88 sigma",
            "source calibration; row overlap unknown",
            REPRO,
            operands,
            "max |log10(model/reference)| / decade",
            note="Paired screening-convention comparisons. Maximum includes both conventions. Divalent and shorter-screening failures remain.",
            fluid_class="Mixed brines",
            property_class="Activity / osmotic",
        )
    fail = csv_at(
        "refit",
        "results/stage-a-calcium/held-out/F-optimum-1/dinane-mounir-failures.csv",
    )
    axes[2].text(
        0.03,
        0.08,
        "One held-out RMS limit exceeded",
        transform=axes[2].transAxes,
        fontsize=6,
    )
    assert any(float(r["RMS"]) > float(r["limit_RMS"]) for r in fail)
    finish(
        fig,
        "F12-mixed-brines.pdf",
        "Training, held-out and unchanged declarations are separate. One retained realization exceeds the stated held-out RMS limit. Divalent panels compare the two screening conventions.",
    )


def categorical(
    ax, data, model, reference, labels, color=BLUE, offset=0, label="Model"
):
    xs = [i + offset for i in range(len(data))]
    ax.plot(
        xs,
        values(data, reference),
        "o",
        color="black",
        mfc="none",
        ms=3,
        mew=0.7,
        label="Reference",
    )
    ax.plot(xs, values(data, model), "x", color=color, ms=3.5, mew=0.8, label=label)
    ax.set_xticks(range(len(data)), labels, rotation=60, ha="right", fontsize=5.5)


def mixed_solvents():
    fig, axes = canvas(12, 3)
    esteso = csv_at("esteso", "results/comparisons.csv")
    figiel = csv_at("figiel", "results/figure-09-comparison.csv")
    cond12 = json_at("held12", "inputs/source/comparison-policy.json")
    # Esteso does not report experimental pressure; the retained model assumes one bar.
    panel(
        axes[0],
        "F13a",
        "NaCl / ethanol-water",
        PRED,
        r"Molality / mol kg$^{-1}$ solvent",
        r"$\gamma_\pm$",
    )
    for j, mass in enumerate(
        sorted({r["ethanol_mass_percent"] for r in esteso}, key=float)
    ):
        comparison(
            axes[0],
            [r for r in esteso if r["ethanol_mass_percent"] == mass],
            "molality_mol_kg",
            "predicted_gamma",
            "gamma_pm_m",
            f"{mass} wt%",
            COLORS[j],
            MARKS[j],
        )
    axes[0].set_xscale("log")
    axes[0].legend(fontsize=5.5, ncol=2)
    record(
        "F13a",
        "NaCl / ethanol-water",
        "Mean ionic activity coefficient",
        values(esteso, "temperature_K"),
        [],
        len(esteso),
        "Figiel fixed mixed-solvent parameters",
        "no new mixture fit",
        PRED,
        pairs(esteso, "predicted_gamma", "gamma_pm_m"),
        note="Experimental pressure unreported; model assumes one bar. Voltage-derived coefficients; development independence unestablished. Same observations recur in F13b.",
        fluid_class="Mixed solvents",
        property_class="Activity / osmotic",
    )
    for ax, system in zip(axes[1:5], sorted({r["system"] for r in figiel})):
        sel = [r for r in figiel if r["system"] == system and r["included"] == "True"]
        suffix = (
            system.removeprefix("water-")
            .replace("methanol", "M")
            .replace("ethanol", "E")
            .replace("nabr", "Br")
            .replace("nacl", "Cl")
        )
        label = "F13b-" + suffix
        panel(
            ax,
            label,
            system.removeprefix("water-").capitalize().replace("-", "/").replace("nacl", "NaCl").replace("nabr", "NaBr"),
            PRED,
            r"Molality / mol kg$^{-1}$ solvent",
            r"$\gamma_\pm$",
        )
        for j, fraction in enumerate(
            sorted({r["x_alcohol_salt_free"] for r in sel}, key=float)
        ):
            comparison(
                ax,
                [r for r in sel if r["x_alcohol_salt_free"] == fraction],
                "molality_mol_kg",
                "value",
                "source_coefficient",
                f"{float(fraction):.2f}",
                COLORS[j % len(COLORS)],
                MARKS[j % len(MARKS)],
            )
        ax.set_xscale("log")
        ax.legend(
            fontsize=5,
            ncol=3,
            title="Salt-free alcohol mole fraction",
            title_fontsize=5.2,
        )
        record(
            label,
            system,
            "Mean ionic activity coefficient",
            values(sel, "T_K"),
            values(sel, "P_Pa"),
            len(sel),
            "Figiel 2025 fixed parameters",
            "no new mixture fit; NaBr pure-solvent activity inputs fitted",
            PRED,
            pairs(sel, "value", "source_coefficient"),
            "RMS / dimensionless",
            note="Fixed-parameter mixture prediction with incomplete primary-source lineage; not literal paper-curve parity or independent model-development evidence. The four populations exclude 57 pure-solvent endpoints in total. All 123 ethanol/NaCl observations repeat Esteso.",
            fluid_class="Mixed solvents",
            property_class="Activity / osmotic",
        )
    identity = json_at("bulow20", "results/run-identity.json")
    hydration = csv_at("bulow20", "results/hydration.csv")
    transfer = csv_at("bulow20", "results/transfer.csv")
    ax = axes[5]
    panel(
        ax,
        "F13c-hydration",
        "Bülow I, hydration",
        REPRO,
        "Ion",
        r"$\Delta G$ / kJ mol$^{-1}$",
    )
    categorical(
        ax,
        hydration,
        "engine_advanced",
        "paper_advanced",
        [r["ion"] for r in hydration],
    )
    ax.plot(
        range(len(hydration)),
        values(hydration, "fawcett"),
        "+",
        color=ORANGE,
        ms=4,
        label="Fawcett",
    )
    ax.legend(fontsize=5.3)
    record(
        "F13c-hydration",
        "Seven hydrated ions",
        "Single-ion hydration Gibbs energy",
        [identity["temperature_K"]],
        [identity["pressure_Pa"]],
        len(hydration),
        "Bülow I advanced parameters",
        "no new fit",
        REPRO,
        pairs(hydration, "engine_advanced", "paper_advanced"),
        "MAE / kJ mol-1",
        note="Metric against paper advanced model; Fawcett separately marked. Convention-dependent single-ion quantity.",
        fluid_class="Solvation",
        property_class="Solvation / transfer",
    )
    ax = axes[6]
    panel(
        ax,
        "F13c-transfer-I",
        "Bülow I, transfer",
        REPRO,
        "Ion / destination",
        r"$\Delta G$ / kJ mol$^{-1}$",
    )
    categorical(
        ax,
        transfer,
        "engine_advanced",
        "paper_advanced",
        [r["ion"] + "/" + r["solvent"][:2] for r in transfer],
    )
    ax.errorbar(
        range(len(transfer)),
        values(transfer, "marcus"),
        yerr=[
            [float(r["marcus"]) - float(r["marcus_low"]) for r in transfer],
            [float(r["marcus_high"]) - float(r["marcus"]) for r in transfer],
        ],
        fmt="+",
        color=ORANGE,
        ms=3.5,
        elinewidth=0.6,
        label="Marcus range",
    )
    ax.legend(fontsize=5.3)
    record(
        "F13c-transfer-I",
        "Water to methanol/ethanol, seven ions",
        "Single-ion transfer Gibbs energy",
        [identity["temperature_K"]],
        [identity["pressure_Pa"]],
        len(transfer),
        "Bülow I advanced parameters",
        "no new fit",
        REPRO,
        pairs(transfer, "engine_advanced", "paper_advanced"),
        "MAE / kJ mol-1",
        note="Metric against paper model. Marcus ranges separately plotted; single-ion convention retained.",
        fluid_class="Solvation",
        property_class="Solvation / transfer",
    )
    transfer2 = csv_at("bulow21", "results/fig8b-transfer-gibbs.csv")
    id2 = json_at("bulow21", "results/run-identity.json")
    t2, p2 = map(
        float, re.search(r"([\d.]+) K, ([\de+.]+) Pa", id2["conditions"]).groups()
    )
    for ax, kij in zip(axes[7:9], ["prediction", "correlation"]):
        label = "F13c-II-" + kij
        panel(
            ax,
            label,
            "Bülow II / " + kij,
            REPRO,
            "Ion / destination",
            r"$\Delta G$ / kJ mol$^{-1}$",
        )
        sel = [
            r
            for r in transfer2
            if r["water_sigma"] == "correlation"
            and r["kij_set"] == kij
            and r["diameter"] == "sigma"
        ]
        categorical(
            ax,
            sel,
            "engine_kJ_mol",
            "paper_model_kJ_mol",
            [r["ion"] + "/" + r["solvent"][:2] for r in sel],
            label=r"Model, $\sigma$",
        )
        unified = [
            r
            for r in transfer2
            if r["water_sigma"] == "correlation"
            and r["kij_set"] == kij
            and r["diameter"] == "unified-0.88"
        ]
        ax.plot(
            range(len(unified)),
            values(unified, "engine_kJ_mol"),
            "+",
            color=GREEN,
            ms=4,
            label=r"Model, unified $0.88\sigma$",
        )
        ax.errorbar(
            range(len(sel)),
            values(sel, "experiment_kJ_mol"),
            yerr=[
                [
                    float(r["experiment_kJ_mol"]) - float(r["experiment_lower_kJ_mol"])
                    for r in sel
                ],
                [
                    float(r["experiment_upper_kJ_mol"]) - float(r["experiment_kJ_mol"])
                    for r in sel
                ],
            ],
            fmt="s",
            mfc="none",
            color=ORANGE,
            ms=2.5,
            elinewidth=0.5,
            label="Reported range",
        )
        ax.legend(fontsize=4.8, ncol=2)
        record(
            label,
            "Water to methanol/ethanol, five ions",
            "Single-ion transfer Gibbs energy",
            [t2],
            [p2],
            len(sel),
            "Bülow II zero-k / Table 8 interactions",
            "source solubility fit for Table 8",
            REPRO,
            pairs(sel, "engine_kJ_mol", "paper_model_kJ_mol"),
            "MAE / kJ mol-1",
            note="Metric uses sigma, paper-model reference. Unified 0.88 sigma shown separately. Experimental scale unspecified; no convention-free measurement claim.",
            fluid_class="Solvation",
            property_class="Solvation / transfer",
        )
    held = csv_at("held12", "results/comparisons.csv")
    held_pairs = sorted({(r["solvent"], r["salt"]) for r in held})
    held_handles = {}
    for ax, quantity in zip(axes[9:], ["density", "ionic", "osmotic"]):
        sel = [
            r for r in held if r["quantity"] == quantity and r["group"] != "Shared #730"
        ]
        label = "F13d-" + quantity
        evidence = CAL if quantity == "density" else REPRO
        ylabel = (
            r"Density / kg m$^{-3}$"
            if quantity == "density"
            else (r"$\gamma_\pm$" if quantity == "ionic" else r"$\Phi$")
        )
        panel(
            ax,
            label,
            "Held alcohol/salt",
            evidence,
            r"Molality / mol kg$^{-1}$ solvent",
            ylabel,
        )
        groups = sorted(
            {(r["solvent"], r["salt"], r["series"], r["evidence_label"]) for r in sel}
        )
        operands = []
        for j, group in enumerate(groups):
            part = [
                r
                for r in sel
                if (r["solvent"], r["salt"], r["series"], r["evidence_label"]) == group
            ]
            style = held_pairs.index(group[:2])
            pair_label = (
                group[0].replace("methanol", "MeOH").replace("ethanol", "EtOH")
                + "/"
                + group[1]
            )
            comparison(
                ax,
                part,
                "molality",
                "computed_value",
                "source_value",
                pair_label,
                COLORS[style % len(COLORS)],
                MARKS[style % len(MARKS)],
            )
            held_handles[pair_label] = ax.lines[-1]
            if group[-1] == (
                "physical comparison"
                if quantity == "density"
                else "numerical verification"
            ):
                operands.extend(pairs(part, "computed_value", "source_value"))
        missing = sum(not r["computed_value"] for r in sel)
        if missing:
            ax.text(
                0.04,
                0.08,
                f"{missing} model values unavailable",
                transform=ax.transAxes,
                fontsize=5.5,
            )
        record(
            label,
            "Alcohol / lithium, sodium and ammonium salts",
            quantity,
            [cond12["temperature_K"]],
            [cond12["pressure_Pa"]],
            len(sel),
            "Held 2012 solvent-specific parameters",
            "source calibration",
            evidence,
            operands,
            "MAE / kg m-3" if quantity == "density" else "MAE / dimensionless",
            note=(
                "Metric against measured density only."
                if quantity == "density"
                else "Metric against digitized paper-model references only."
            )
            + " All selected source values remain plotted, including unavailable model values. No independent validation; source-correction limits retained.",
            fluid_class="Mixed solvents",
            property_class="Density" if quantity == "density" else "Activity / osmotic",
        )
    fig.legend(
        list(held_handles.values()),
        list(held_handles),
        loc="lower center",
        bbox_to_anchor=(0.5, 0.025),
        ncol=5,
        fontsize=5.5,
        title="Held alcohol/salt series (bottom row)",
        title_fontsize=6,
    )
    finish(
        fig,
        "F13-mixed-solvents-solvation.pdf",
        "Esteso and the ethanol/NaCl Figiel subset are duplicate observations, not extra evidence. Figiel has incomplete primary lineage. Bülow panels distinguish paper-model reproduction, reference ranges, and diameter conventions. Held selects density, ionic activity and osmotic comparisons; pressure and solvent activity are not added. Calibration and source-correction limits remain.",
    )


def reacting_lle():
    import numpy as np

    fig, axes = canvas(12, 3)
    raw = list(
        csv.reader(io.StringIO(retained("khudaida", "results/figure-2-comparison.csv")))
    )
    data = [
        {
            "tie": int(r[0]),
            "phase": r[1],
            "model": list(map(float, r[2:6])),
            "reference": list(map(float, r[6:10])),
        }
        for r in raw[3:]
    ]
    conditions = csv_at("khudaida", "inputs/source/table-3-4-experimental-tielines.csv")
    conditions = [
        r for r in conditions
        if r["source_table"] == "Table 3"
        and float(r["temperature_K"]) == 293.15
        and float(r["salt_wtfrac"]) == 0.05
    ]
    uncertainty = csv_at("khudaida", "inputs/source/table-3-standard-uncertainties.csv")
    u_salt = next(
        float(r["standard_uncertainty_mole_fraction"])
        for r in uncertainty
        if r["component"] == "nacl"
    )
    ax = axes[0]
    panel(
        ax,
        "F14a-ties",
        "Khudaida / salt-free projection",
        PRED,
        "Salt-free water mole fraction",
        "Salt-free isobutanol mole fraction",
    )
    ax.plot([0, 1, 0], [0, 0, 1], color="0.6", lw=0.5)
    for tie in sorted({r["tie"] for r in data}):
        part = [r for r in data if r["tie"] == tie]
        for key, style, color in [("reference", "o--", "black"), ("model", "x-", BLUE)]:
            # Projection only: NaCl is removed and the three solvent coordinates renormalized.
            ax.plot(
                [r[key][0] / sum(r[key][:3]) for r in part],
                [r[key][2] / sum(r[key][:3]) for r in part],
                style,
                ms=2.8,
                lw=0.6,
                color=color,
                mfc="none",
            )
    ax.set(xlim=(0, 1), ylim=(0, 0.7))
    ax.plot(
        [], [], "o--", color="black", mfc="none", ms=3, lw=0.6, label="Measured ties"
    )
    ax.plot([], [], "x-", color=BLUE, ms=3, lw=0.6, label="Model ties")
    ax.legend(fontsize=5.5)
    record(
        "F14a-ties",
        "Water/ethanol/isobutanol/NaCl",
        "Four formula-basis phase fractions",
        values(conditions, "temperature_K"),
        [],
        len(data) // 2,
        "Khudaida 2026 literature interactions",
        "no",
        PRED,
        [(a, b) for r in data for a, b in zip(r["model"], r["reference"])],
        "MAE / mole fraction",
        note="Metric uses all four unprojected formula coordinates. Segments join coexisting phases, not a continuous binodal. Model pressure is recorded with the declaration; source table has no pressure column.",
        fluid_class="Electrolyte LLE",
        property_class="LLE",
    )
    ax = axes[1]
    panel(
        ax,
        "F14a-salt",
        "Khudaida / NaCl formula fraction",
        PRED,
        "Tie line",
        "NaCl formula mole fraction",
    )
    for phase, mark, color, offset in [
        ("aqueous", "o", BLUE, -0.08),
        ("organic", "s", ORANGE, 0.08),
    ]:
        part = [r for r in data if r["phase"] == phase]
        xs = [r["tie"] + offset for r in part]
        ax.errorbar(
            xs,
            [r["reference"][3] for r in part],
            yerr=u_salt,
            fmt=mark,
            color=color,
            mfc="none",
            ms=3,
            elinewidth=0.6,
            label=phase,
        )
        ax.plot(xs, [r["model"][3] for r in part], "x", color=color, ms=3)
    ax.legend(fontsize=5.5)
    record(
        "F14a-salt",
        "Water/ethanol/isobutanol/NaCl",
        "NaCl formula phase fraction",
        values(conditions, "temperature_K"),
        [],
        len(data),
        "Khudaida 2026 literature interactions",
        "no",
        PRED,
        [(r["model"][3], r["reference"][3]) for r in data],
        "MAE / mole fraction",
        note="One-standard-uncertainty bars from Table 3. Salt stays visible beside the solvent projection.",
        fluid_class="Electrolyte LLE",
        property_class="LLE",
    )
    ax = axes[2]
    calc = csv_at("held14", "results/figure-6a.csv")
    lines = [
        r
        for r in csv_at("held14", "inputs/paper-lines.csv")
        if r["figure"] == "figure-6a"
    ]
    markers = [
        r
        for r in csv_at("held14", "inputs/paper-markers.csv")
        if r["figure"] == "figure-6a"
    ]
    id14 = json_at("held14", "results/run-identity.json")
    p14 = float(re.search(r"([\d.]+) bar", id14["conditions"])[1]) * 1e5
    t14 = float(re.search(r"; ([\d.]+) K", id14["conditions"])[1])
    panel(
        ax,
        "F14b",
        "1-butanol/NH4Cl/water",
        REPRO,
        "Aqueous NH4Cl / wt%",
        "1-butanol / wt%",
    )
    operands = []
    for series, column, color in [
        ("1-butanol phase", "w_butanol_organic", BLUE),
        ("aqueous phase", "w_butanol_aqueous", ORANGE),
    ]:
        trace = sorted(
            [r for r in lines if r["series"] == series], key=lambda r: float(r["x"])
        )
        ax.plot(
            values(trace, "x"),
            values(trace, "y"),
            "o",
            ms=1.5,
            color=color,
            mfc="none",
            mew=0.4,
        )
        ax.plot(
            values(calc, "w_nh4cl_aqueous"),
            values(calc, column),
            "-",
            color=color,
            lw=1,
            label=series,
        )
        part = [
            r
            for r in calc
            if min(values(trace, "x"))
            <= float(r["w_nh4cl_aqueous"])
            <= max(values(trace, "x"))
        ]
        operands.extend(
            (
                float(r[column]),
                float(
                    np.interp(
                        float(r["w_nh4cl_aqueous"]),
                        values(trace, "x"),
                        values(trace, "y"),
                    )
                ),
            )
            for r in part
        )
    ax.plot(
        values(markers, "x"),
        values(markers, "y"),
        "s",
        mfc="none",
        color="black",
        ms=3,
        label="Measured symbols",
    )
    ax.legend(fontsize=5.5)
    record(
        "F14b",
        "1-butanol/NH4Cl/water",
        "Butanol phase mass percentage",
        [t14],
        [p14],
        len(calc),
        "Held 2014 fitted ion-butanol corrections",
        "source calibration",
        REPRO,
        operands,
        "max absolute wt%-point difference",
        note="Metric interpolates digitized model lines inside each branch domain only; no extrapolation. Experimental symbols separately shown.",
        fluid_class="Electrolyte LLE",
        property_class="LLE",
    )
    acid = csv_at("acid", "results/comparison.csv")
    sweep = csv_at("acid", "results/sweep-aq-literal.csv")
    ca = json_at("acid", "results/provenance.json")
    for path in ["results/comparison.csv", "results/sweep-aq-literal.csv"]:
        # The unchanged original renderer is reused only after checking its input identity.
        assert (A / PROMOTED["acid"][0] / path).read_text() == retained("acid", path)
    acetic_acid_ph(axes[3:5], (sweep, acid))
    for ax, label, title, quantity, evidence, operands, name in [
        (
            axes[3],
            "F14c-pH",
            "Acetic acid / water-toluene",
            "Aqueous pH",
            PRED,
            pairs(acid, "pH_aq-literal", "pH_exp"),
            "MAE / pH unit",
        ),
        (
            axes[4],
            "F14c-difference",
            "Source-model disagreement",
            "Aqueous pH minus published model",
            REPRO,
            pairs(acid, "pH_aq-literal", "pH_esi_model"),
            "max absolute pH difference",
        ),
    ]:
        panel(
            ax,
            label,
            title,
            evidence,
            "Aqueous analytical acid fraction",
            "Aqueous pH" if evidence == PRED else "pH difference",
        )
        ax.legend(fontsize=4.8)
        record(
            label,
            "Acetic acid/water/toluene",
            quantity,
            [ca["T_K"]],
            [ca["P_Pa"]],
            len(acid),
            "Ascani 2022 published parameters",
            "no pH fit",
            evidence,
            operands,
            name,
            note="Failed published-model comparison retained; pH error is not relative hydrogen-ion activity error.",
            fluid_class="Reacting liquids",
            property_class="pH / partition",
        )
    ax = axes[5]
    panel(
        ax,
        "F14c-K",
        "Acetic acid partition",
        CAL,
        "Aqueous analytical acid fraction",
        r"$K=x_{acid,aq}/x_{acid,org}$",
    )
    sweep = sorted(sweep, key=lambda r: float(r["x_acid_total_aq"]))
    ax.plot(
        values(sweep, "x_acid_total_aq"),
        values(sweep, "K_aq_org_molecular"),
        color=BLUE,
        lw=1,
        label="Model sweep",
    )
    comparison(ax, acid, "x_acetic_acid_aq", "K_aq-literal", "K_exp")
    ax.plot(
        values(acid, "x_acetic_acid_aq"),
        values(acid, "K_esi_model"),
        "+",
        color=GREEN,
        ms=3,
        label="Published model",
    )
    ax.legend(fontsize=5)
    record(
        "F14c-K",
        "Acetic acid/water/toluene",
        "Molecular acid partition coefficient",
        [ca["T_K"]],
        [ca["P_Pa"]],
        len(acid),
        "Ascani 2022 acid-organic fit",
        "calibration family; exact rows unknown",
        CAL,
        pairs(acid, "K_aq-literal", "K_exp"),
        note="No partition row is demonstrated held out; published-model misses retained.",
        fluid_class="Reacting liquids",
        property_class="pH / partition",
    )
    ties = [
        r
        for r in csv_at("ascani", "results/reactive-tie-lines.csv")
        if r["rule"] == "printed"
    ]
    refs = csv_at("ascani", "inputs/data/senina-2022-table-3-reactive-tie-lines.csv")
    ah = [
        r
        for r in csv_at("ascani", "results/homogeneous-ce.csv")
        if r["rule"] == "printed"
    ]
    ar = csv_at("ascani", "inputs/data/senina-2022-table-2-homogeneous-ce.csv")
    ac = json_at("ascani", "results/provenance.json")
    ax = axes[6]
    panel(
        ax,
        "F14d-ties",
        "Esterification tie lines",
        REPRO,
        "Water mole fraction",
        "Pentyl acetate mole fraction",
    )
    operands = []
    for r in ties:
        ref = next(v for v in refs if v["tie_line"] == r["tie_line"])
        for d, style, color in [(ref, "o--", "black"), (r, "x-", BLUE)]:
            ax.plot(
                [float(d[p + "_x_water"]) for p in ["aq", "org"]],
                [float(d[p + "_x_pentyl_acetate"]) for p in ["aq", "org"]],
                style,
                color=color,
                ms=3,
                mfc="none",
                lw=0.6,
            )
        operands.extend((float(r[k]), float(ref[k])) for k in ref if k != "tie_line")
    ax.plot([], [], "o--", color="black", mfc="none", ms=3, label="Measured ties")
    ax.plot([], [], "x-", color=BLUE, ms=3, label="Model")
    ax.legend(fontsize=5.5)
    record(
        "F14d-ties",
        "Acetic acid/1-pentanol/pentyl acetate/water",
        "Four-component two-liquid compositions",
        [ac["temperatures_k"]["1"]],
        [ac["pressure_pa"]],
        len(ties),
        "Ascani 2023; printed association rule",
        "partly; reaction constant fitted",
        REPRO,
        operands,
        "RMS / mole fraction",
        note="Projection shows water/ester coordinates; metric includes all four components and both phases. Source 1.013-bar wording differs from retained one-bar model.",
        fluid_class="Reacting liquids",
        property_class="LLE",
    )
    ax = axes[7]
    panel(
        ax,
        "F14d-homogeneous",
        "Homogeneous reference states",
        REPRO,
        "Measured point",
        "Acetic acid mole fraction",
    )
    operands = []
    for r in ah:
        ref = next(v for v in ar if v["point"] == r["point"])

        def overall(species, r=r):
            return sum(
                float(r[f"phase{i}_amount"]) * float(r[f"phase{i}_x_{species}"])
                for i in [1, 2]
                if r[f"phase{i}_amount"]
            ) / sum(
                float(r[f"phase{i}_amount"]) for i in [1, 2] if r[f"phase{i}_amount"]
            )

        ax.plot(
            int(r["point"]),
            float(ref["x_acetic_acid"]),
            "o",
            color="black",
            mfc="none",
            ms=3,
        )
        ax.plot(int(r["point"]), overall("acetic_acid"), "x", color=BLUE, ms=3)
        if int(r["phases"]) > 1:
            ax.plot(
                int(r["point"]),
                overall("acetic_acid"),
                "s",
                color=ORANGE,
                mfc="none",
                ms=6,
            )
        operands.extend(
            (overall(k.removeprefix("x_")), float(ref[k])) for k in ref if k != "point"
        )
    ax.plot([], [], "s", color=ORANGE, mfc="none", label="Predicted split")
    ax.legend(fontsize=5.5)
    record(
        "F14d-homogeneous",
        "Acetic acid/1-pentanol/pentyl acetate/water",
        "Overall four-component composition",
        [ac["temperatures_k"]["1"]],
        [ac["pressure_pa"]],
        len(ah),
        "Ascani 2023; Senina-derived reaction constant",
        "partly; unnamed source row",
        REPRO,
        operands,
        "RMS / mole fraction",
        note="Overall fractions reconstructed by phase amounts. Homogeneous measured states predicted to split remain marked; invariant feed agreement is not a reaction-extent test.",
        fluid_class="Reacting liquids",
        property_class="Reaction amounts",
    )
    chase = csv_at("chase", "results/comparison.csv")
    ax = axes[8]
    cp = [r for r in chase if r["quantity"] == "Cp_J_mol_K"]
    panel(
        ax,
        "F14e-Cp",
        "Chase / printed properties",
        NV,
        "Temperature / K",
        r"$C_p$ / J mol$^{-1}$ K$^{-1}$",
    )
    for j, species in enumerate(sorted({r["case"] for r in cp})):
        comparison(
            ax,
            [r for r in cp if r["case"] == species],
            "T_K",
            "value",
            "reference",
            species,
            COLORS[j],
        )
    ax.legend(fontsize=5.2, ncol=2)
    record(
        "F14e-Cp",
        "CO/CO2/H2/H2O/N2O4/NO2",
        "Printed ideal-gas heat capacity",
        values(cp, "T_K"),
        values(cp, "P_Pa"),
        len(cp),
        "Chase/NIST Shomate coefficients",
        "no",
        NV,
        pairs(cp, "value", "reference"),
        "max absolute J mol-1 K-1 difference",
        note="Printed-property numerical verification, not physical caloric validation.",
        fluid_class="Ideal reacting gas",
        property_class="Calorics",
    )
    for ax, case, col, scale, xlabel in [
        (axes[9], "water-gas-shift", "T_K", 1, "Temperature / K"),
        (axes[10], "n2o4-no2", "P_Pa", 1e5, "Pressure / bar"),
    ]:
        sel = [
            dict(r, x_plot=float(r[col]) / scale)
            for r in chase
            if r["case"] == case and r["quantity"].endswith("_amount_mol")
        ]
        label = "F14e-" + case
        title = {"water-gas-shift": "Water-gas shift", "n2o4-no2": r"N$_2$O$_4$/NO$_2$"}[case]
        panel(ax, label, title, NV, xlabel, "Species amount / mol")
        for j, species in enumerate(sorted({r["quantity"] for r in sel})):
            comparison(
                ax,
                [r for r in sel if r["quantity"] == species],
                "x_plot",
                "value",
                "reference",
                species.removesuffix("_amount_mol"),
                COLORS[j],
                MARKS[j],
            )
        ax.legend(fontsize=5.2, ncol=2)
        record(
            label,
            case,
            "Ideal-gas equilibrium amounts",
            values(sel, "T_K"),
            values(sel, "P_Pa"),
            len(sel),
            "Chase/NIST source thermochemistry",
            "no",
            NV,
            pairs(sel, "value", "reference"),
            "max absolute amount difference / mol",
            note="Analytic/source-derived numerical references, not measured compositions.",
            fluid_class="Ideal reacting gas",
            property_class="Reaction amounts",
        )
    belov = [
        r
        for r in csv_at("belov", "results/comparison.csv")
        if r["quantity"].endswith("_amount_mol")
    ]
    bc = json_at("belov", "inputs/cases/test-01/source/source.json")
    ax = axes[11]
    panel(
        ax, "F14e-Belov", "Eight-gas restriction", NV, "Species", "Species amount / mol"
    )
    categorical(
        ax,
        belov,
        "value",
        "reference",
        [r["quantity"].removesuffix("_amount_mol") for r in belov],
    )
    ax.set_yscale("log")
    ax.legend(fontsize=5.5)
    record(
        "F14e-Belov",
        "O/O2/O3/C/C2/CO/CO2/C2O",
        "Ideal-gas equilibrium amounts",
        [bc["temperature_k"]],
        [bc["pressure_pa"]],
        len(belov),
        "Belov-Aristova fixed Gibbs encoding",
        "no",
        NV,
        pairs(belov, "value", "reference"),
        "max |ln(model/reference)|",
        note="Independent high-precision gas-only numerical reference; ozone retained. Excludes graphite, diamond, heterogeneous equilibrium and caloric extrapolation.",
        fluid_class="Ideal reacting gas",
        property_class="Reaction amounts",
    )
    finish(
        fig,
        "F14-lle-reacting.pdf",
        "Khudaida includes all formula-basis NaCl coordinates beside the solvent projection. Held is a calibrated paper-model reproduction. Acetic-acid pH retains source-model misses; partitioning is calibration-related. Esterification keeps incorrect phase counts. Chase and Belov are numerical verification only.",
    )


def algorithm_figures():
    fig, axes = canvas(4)
    base = "results/main-8298c602-pressure-residual/"
    cases = json_at("density", base + "cases.json")
    records = [
        json_at("density", base + f"{c['id']}-state-{s['index']}.json")
        for c in cases["cases"]
        for s in c["states"]
    ]
    for ax, label, chosen, title in [
        (axes[0], "F2-propane", ["A-propane-vapor"], "Propane, 300 K"),
        (
            axes[1],
            "F2-electrolyte",
            ["B-C3-1-mol-kg", "B-C3-3-mol-kg", "B-C3-6-mol-kg"],
            "Aqueous NaCl, 298.15 K",
        ),
    ]:
        selected = [
            r for r in records if r["case"] in chosen and r["state"]["index"] == 0
        ]
        panel(ax, label, title, NV, r"Density / mol m$^{-3}$", "Pressure / MPa")
        operands = []
        for j, r in enumerate(selected):
            scan = r["scans"][0]
            for quantity, mark in [("roots", "o"), ("spinodals", "x")]:
                pts = scan[quantity]
                ax.plot(
                    [p["molar_density"] for p in pts],
                    [p["pressure"] / 1e6 for p in pts],
                    mark,
                    color=COLORS[j],
                    mfc="none",
                    ms=4,
                    label=(
                        r["case"].removeprefix("B-C3-") + " "
                        if len(selected) > 1
                        else ""
                    )
                    + quantity,
                )
            operands.extend(
                (r["scans"][0]["roots"][i]["molar_density"] / p["molar_density"], 1)
                for i, p in enumerate(r["scans"][1]["roots"])
            )
        ax.legend(fontsize=5.6)
        record(
            label,
            title,
            "Root/spinodal locations, nested-grid root comparison",
            [r["state"]["T_K"] for r in selected],
            [r["state"]["P_Pa"] for r in selected],
            sum(
                len(r["scans"][0][q]) for r in selected for q in ["roots", "spinodals"]
            ),
            "Declared propane and aqueous NaCl parameters",
            "no",
            NV,
            operands,
            "max absolute relative density difference",
            note="Discrete crossings only; no interpolated P-rho curve or complete-root/global-stability claim.",
            fluid_class="Numerical methods",
            property_class="Density",
        )
    pressure_errors, density_errors, temperatures, pressures = [], [], [], []
    for r in records:
        state = r["state"]
        temperatures.append(state["T_K"])
        pressures.append(state["P_Pa"])
        for v in r["returns"]:
            pressure_errors.append(
                abs(v["pressure_Pa"] - state["P_Pa"])
                / (v["rho_mol_m3"] * cases["gas_constant_J_mol_K"] * state["T_K"])
            )
            if v["mode"] == "pressure":
                density_errors.append(
                    abs(v["rho_mol_m3"] / state["oracle_rho_mol_m3"] - 1)
                )
    for ax, label, errors, title, xlabel in [
        (
            axes[2],
            "F2-pressure",
            pressure_errors,
            "Pressure closure",
            r"$|P_{calc}-P|/(\rho RT)$",
        ),
        (
            axes[3],
            "F2-density",
            density_errors,
            "Density agreement",
            "Absolute relative density difference",
        ),
    ]:
        panel(ax, label, title, NV, xlabel, "Empirical cumulative fraction")
        positive = sorted(x for x in errors if x > 0)
        zero_count = len(errors) - len(positive)
        # Zero residuals are counted explicitly, never moved to an arbitrary log floor.
        ax.plot(
            positive,
            [(i + zero_count + 1) / len(errors) for i in range(len(positive))],
            "o",
            color=BLUE,
            ms=2,
            mfc="none",
            mew=0.5,
        )
        ax.text(
            0.04,
            0.93,
            f"{zero_count} requests with zero difference (not shown)",
            transform=ax.transAxes,
            fontsize=6,
        )
        ax.set_xscale("log")
        record(
            label,
            "All declared density benchmark sequences",
            title,
            temperatures,
            pressures,
            len(errors),
            "Declared propane and aqueous NaCl parameters",
            "no",
            NV,
            [(x, 0) for x in errors],
            "max absolute dimensionless difference",
            note="All per-state JSON returns reduced directly. Binary64 zeros counted separately; no artificial positive floor.",
            fluid_class="Numerical methods",
            property_class="Density",
        )
    finish(
        fig,
        "F2-density-roots.pdf",
        "Root and spinodal markers from the first retained states; distributions use every per-state density-benchmark records. Shared EOS evaluation is numerical verification, not physical validation.",
    )

    fig, axes = canvas(1, 1)
    ax = axes[0]
    study = json_at("algorithms", "inputs/study.json")
    trial = csv_at("algorithms", "results/tpd.csv")
    search = csv_at("algorithms", "results/search.csv")
    earlier = csv_at("discovery", "results/main-8298c602/feos.csv")
    # Same source case; the earlier promotion supplies finite search points, not a continuous path.
    previous = [
        float(r["value"])
        for r in earlier
        if re.fullmatch(r"case-5/FeOs/explicit-test/trial-[012]/quantities", r["case"])
        and r["quantity"] == "tpd_over_rt"
    ]
    assert sorted(previous) == sorted(values(search, "tpd_over_rt"))
    panel(
        ax,
        "F3",
        "Water/hexane, 300 K, 0.1 MPa",
        NV,
        "Water mole fraction (logit scale)",
        r"Tangent-plane distance / $RT$",
    )
    ax.plot(
        values(trial, "x_water"),
        values(trial, "tpd_over_rt"),
        color=BLUE,
        lw=1.2,
        label="Tangent-plane distance",
    )
    ax.plot(
        values(search, "x_water"),
        values(search, "tpd_over_rt"),
        "x",
        color=ORANGE,
        ms=6,
        label="Native liquid trial search",
    )
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xscale("logit")
    ax.legend(fontsize=7)
    # Direct evaluation at the retained search coordinates is present in the swept table.
    operands = []
    for r in search:
        t = min(trial, key=lambda t: abs(float(t["x_water"]) - float(r["x_water"])))
        assert abs(float(t["x_water"]) - float(r["x_water"])) < 1e-14
        operands.append((float(t["tpd_over_rt"]), float(r["tpd_over_rt"])))
    record(
        "F3",
        "Water/hexane",
        "Tangent-plane distance / RT",
        [study["tpd_temperature_K"]],
        [study["tpd_pressure_Pa"]],
        len(trial),
        "Declared water/hexane parameters",
        "no",
        NV,
        operands,
        "max absolute TPD/RT difference",
        note="Composition sweep, not optimizer iterates. Includes native search compositions by construction. No global minimum certificate.",
        fluid_class="Numerical methods",
        property_class="Stability",
    )
    finish(
        fig,
        "F3-tpd-path.pdf",
        "Numerical composition-path illustration; search points agree with the separately retained discovery values. No exhaustive search or independent discovery claim.",
    )

    fig, axes = canvas(1, 1)
    ax = axes[0]
    study = json_at("algorithms", "inputs/study.json")
    envelope = csv_at("algorithms", "results/envelope.csv")
    trace = json_at("algorithms", "results/trace.json")
    # These pre-existing reference points are explicitly incorporated into the promoted algorithm study.
    ref = csv_at("algorithms", study["envelope_reference"], absolute=True)
    cross = csv_at("algorithms", "results/crosschecks.csv")
    panel(
        ax,
        "F4a",
        "Equal-feed methane/ethane phase envelope",
        NV,
        "Temperature / K",
        "Pressure / MPa",
    )
    for bubble, color, style, label in [
        ("True", BLUE, "-", "Bubble"),
        ("False", ORANGE, "--", "Dew"),
    ]:
        part = [r for r in envelope if r["bubble"] == bubble]
        ax.plot(
            values(part, "T_K"),
            [p / 1e6 for p in values(part, "P_Pa")],
            style,
            color=color,
            label=label,
        )
    ax.plot(
        values(ref, "T_K"),
        [p / 1e6 for p in values(ref, "P_Pa")],
        "o",
        mfc="none",
        color="black",
        ms=4,
        label="FeOs",
    )
    for k in trace["critical"]:
        ax.plot(
            k["T_K"],
            k["P_Pa"] / 1e6,
            "*",
            color="black",
            ms=8,
            label="Critical candidate",
        )
    ax.legend(fontsize=7)
    # The retained crosschecks store relative interpolation differences, not new EOS calculations.
    operands = [
        (float(r[c]), 0)
        for r in cross
        for c in ["P_Pa", "rho_feed", "rho_incipient", "x0_incipient"]
    ]
    record(
        "F4a",
        "Equal-feed methane/ethane",
        "Bubble/dew pressure and critical candidate",
        values(envelope, "T_K"),
        values(envelope, "P_Pa"),
        len(envelope),
        study["envelope_packet"].removeprefix("data/packets/"),
        "no",
        NV,
        operands,
        "max relative interpolation difference",
        note="Metric includes pressure, both densities and incipient composition at retained FeOs points. Critical candidate gives no global-stability certificate.",
        fluid_class="Numerical methods",
        property_class="VLE",
    )
    finish(
        fig,
        "F4-phase-envelope.pdf",
        "Traced native continuation, fixed feed composition; retained FeOs points and one critical candidate. Numerical illustration only.",
    )


def gross_pure():
    data = csv_at("gross01", "results/pure-saturation.csv")
    assert len(data) == 78 and all(r["numerical_tolerance_met"] == "True" for r in data)
    fig, axes = canvas(1, 1)
    ax = axes[0]
    panel(
        ax,
        "F8b-catalogue",
        "78 calculated compounds; no measured references",
        NV,
        "Temperature / K",
        "Calculated saturation pressure / MPa",
    )
    ax.plot(
        values(data, "temperature_k"),
        [p / 1e6 for p in values(data, "pressure_pa")],
        "x",
        ms=4,
        color=BLUE,
    )
    ax.set_yscale("log")
    record(
        "F8b-catalogue",
        "78 pure compounds",
        "Saturation calculation catalogue",
        values(data, "temperature_k"),
        values(data, "pressure_pa"),
        len(data),
        "Gross-Sadowski 2001 Table 2",
        "original fits; no references here",
        NV,
        [(float(r["phase_pressure_relative_closure"]), 0) for r in data],
        "max relative phase-pressure closure",
        note="One reference-free calculation per compound; breadth is not validation of 78 compounds.",
        fluid_class="Pure fluids",
        property_class="Saturation",
    )
    finish(
        fig,
        "F8-pure-fluids.pdf",
        "The 78-compound catalogue has no measured references or accuracy metric.",
    )
    comparisons = csv_at("gross01", "results/figure-comparisons.csv")
    assert (
        len(comparisons) == 122
        and sum(r["numerical_tolerance_met"] == "True" for r in comparisons) == 120
    )
    data = [r for r in comparisons if r["kind"] == "density"]
    fig, axes = canvas(4)
    for ax, compound in zip(axes, sorted({r["component"] for r in data})):
        selected = [r for r in data if r["component"] == compound]
        label = "F8b-" + compound.replace(" ", "-")
        panel(
            ax,
            label,
            compound + " liquid",
            REPRO + " / " + PRED,
            "Temperature / K",
            r"Liquid density / kg m$^{-3}$",
        )
        for source, mark, color in [("model", "s", BLUE), ("experiment", "o", ORANGE)]:
            part = [
                dict(
                    r,
                    qualified=r["predicted"]
                    if r["numerical_tolerance_met"] == "True"
                    else "",
                )
                for r in selected
                if r["reference"] == source
            ]
            comparison(
                ax,
                part,
                "temperature_k",
                "qualified",
                "observed",
                source,
                color,
                mark,
                uncertainty="x_uncertainty",
            )
            valid = [r for r in part if r["qualified"]]
            record(
                label + "-" + source,
                compound,
                "Liquid density",
                values(part, "temperature_k"),
                values(part, "pressure_pa"),
                len(part),
                "Gross-Sadowski 2001 Table 2",
                "original fit overlap unknown",
                REPRO if source == "model" else PRED,
                pairs(valid, "qualified", "observed"),
                note="Sparse figure samples; uncertainty bars are reading intervals, not measurement uncertainty. "
                "Cold propylene and diethyl ether are numerically excluded, with finite returned values retained separately.",
                fluid_class="Pure fluids",
                property_class="Density",
            )
        excluded = [r for r in selected if r["numerical_tolerance_met"] != "True"]
        if excluded:
            ax.plot(
                values(excluded, "temperature_k"),
                values(excluded, "predicted"),
                "D",
                mfc="none",
                color="black",
                ms=5,
                label="Pressure closure failed",
            )
        ax.legend(fontsize=5.5, ncol=2)
    finish(
        fig,
        "F8-pure-fluids.pdf",
        "Sparse Gross liquid-density comparisons preserve both cold numerical exclusions.",
    )


def dipole_pure():
    data = csv_at("dipole", "results/pure-rows.csv")
    notes = "Dipolar pure fluids: evaluated reference equations and DMSO measurements, published parameters and local refits. No isolated fixed-parameter dipole effect is inferred."
    for i, compound in enumerate(sorted({r["component_id"] for r in data})):
        if i % 2 == 0:
            fig, axes = canvas(4)
            csv_at("dipole", "results/pure-rows.csv")
        selected = [r for r in data if r["component_id"] == compound]
        for j, (quantity, unit, scale) in enumerate(
            [
                ("p", "Saturation pressure / MPa", 1e6),
                ("rho", r"Liquid density / kmol m$^{-3}$", 1000),
            ]
        ):
            ax = axes[2 * (i % 2) + j]
            label = f"F8c-{compound}-{quantity}"
            panel(
                ax,
                label,
                compound.replace("-", " "),
                REPRO + " / calibration / prediction",
                "Temperature / K",
                unit,
            )
            suffix = "Pa" if quantity == "p" else "mol_m3"
            reference, model = f"{quantity}_ref_{suffix}", f"{quantity}_model_{suffix}"
            ref = [
                r for r in selected if r["model"] == "polar/published" and r[reference]
            ]
            unc = "u_p_sat_rel" if quantity == "p" else "u_rho_liq_rel"
            ax.errorbar(
                values(ref, "T_K"),
                [x / scale for x in values(ref, reference)],
                yerr=[
                    float(r[reference]) * float(r[unc]) / scale if r[unc] else 0
                    for r in ref
                ],
                fmt="o",
                mfc="none",
                color="black",
                ms=2.5,
                lw=0.5,
                label="DMSO measurements"
                if "sulfoxide" in compound
                else "Reference equation",
            )
            for k, name in enumerate(
                [
                    "polar/published",
                    "nonpolar/published",
                    "polar/refit",
                    "nonpolar/refit",
                ]
            ):
                part = [r for r in selected if r["model"] == name and r[reference]]
                if not part:
                    continue
                ax.plot(
                    values(part, "T_K"),
                    [x / scale for x in values(part, model)],
                    ["x", "+", "1", "2"][k],
                    ms=4,
                    color=COLORS[k],
                    label=name,
                )
                for role in (
                    ["calibration", "validation"] if "refit" in name else [None]
                ):
                    cohort = [r for r in part if role is None or r["role"] == role]
                    evidence = (
                        CAL
                        if role == "calibration"
                        else PRED
                        if role == "validation"
                        else REPRO
                    )
                    record(
                        label
                        + "-"
                        + name.replace("/", "-")
                        + ("-" + role if role else ""),
                        compound,
                        unit.split(" / ")[0],
                        values(cohort, "T_K"),
                        [float(r["p_ref_Pa"]) for r in cohort if r["p_ref_Pa"]],
                        len(cohort),
                        name,
                        "local alternating-row fit"
                        if role
                        else "original fit rows unavailable",
                        evidence,
                        pairs(cohort, model, reference),
                        note="Published parameters and local refits are distinct. Alternating withholding applies only to local refits; "
                        "reference equations may share experimental ancestry with original fits. DMSO measurements postdate published parameters. "
                        "Stated source bounds are not predictive uncertainty; failed published pressure criteria remain in the comparison.",
                        fluid_class="Pure fluids",
                        property_class="Saturation" if quantity == "p" else "Density",
                    )
                for r in part:
                    PLOTTED.append(
                        {
                            "panel": label,
                            "series": name,
                            "role": r["role"]
                            if "refit" in name
                            else "original overlap unknown",
                            "x": float(r["T_K"]),
                            "reference": float(r[reference]) / scale,
                            "model": float(r[model]) / scale,
                        }
                    )
            if quantity == "p":
                ax.set_yscale("log")
            ax.legend(fontsize=6.5, ncol=2)
        if i % 2 == 1:
            finish(fig, "F8-pure-fluids.pdf", notes)


def dipole_binary():
    data = csv_at("dipole", "results/binary-xy-rows.csv")
    fig, axes = canvas(6)
    for i, temp in enumerate(sorted(set(values(data, "T_K")))):
        for j, (model, reference, ylabel, metric_name) in enumerate(
            [
                ("P_model_kPa", "P_data_kPa", "Bubble pressure / kPa", "AARD / %"),
                (
                    "y_model",
                    "y_data",
                    "Vapor pentane mole fraction",
                    "MAE / mole fraction",
                ),
            ]
        ):
            ax = axes[2 * i + j]
            label = f"F9b-{temp:g}-" + ("P" if j == 0 else "y")
            evidence = CAL if i == 0 else PRED
            panel(
                ax,
                label,
                f"Acetone/n-pentane, {temp:g} K",
                evidence + "; Margules: fit",
                "Liquid pentane mole fraction",
                ylabel,
            )
            selected = [r for r in data if float(r["T_K"]) == temp]
            ref = [
                r
                for r in selected
                if r["model"] == "polar/GV2006" and r["variant"] == "k_ij calibrated"
            ]
            ax.plot(
                values(ref, "x_pentane"),
                values(ref, reference),
                "o",
                mfc="none",
                color="black",
                ms=3,
                label="Campbell observations",
            )
            for k, (name, variant, title) in enumerate(
                [
                    ("polar/GV2006", "k_ij calibrated", "Polar, fitted interaction"),
                    (
                        "nonpolar/refit",
                        "k_ij calibrated",
                        "Nonpolar, fitted interaction",
                    ),
                    ("polar/GV2006", "k_ij=0", "Polar, zero interaction"),
                    (
                        "Margules/Campbell1986",
                        "four-suffix eq 8 (Barker P-x fit) per isotherm",
                        "Margules, isotherm fit",
                    ),
                ]
            ):
                part = [
                    r
                    for r in selected
                    if r["model"] == name and r["variant"] == variant
                ]
                assert part
                ax.plot(
                    values(part, "x_pentane"),
                    values(part, model),
                    ["x", "+", "1", "2"][k],
                    color=COLORS[k],
                    ms=4,
                    label=title,
                )
                cls = CAL if k == 3 else PRED if k == 2 else evidence
                record(
                    label + f"-{k}",
                    "Acetone/n-pentane",
                    ylabel,
                    [temp],
                    [p * 1000 for p in values(part, "P_data_kPa")],
                    len(part),
                    title,
                    "per-isotherm fit"
                    if k == 3
                    else "372.7 K interaction fit"
                    if k < 2
                    else "no binary fit",
                    cls,
                    pairs(part, model, reference),
                    metric_name,
                    note="Published polar parameters, nonpolar pure refit, calibrated interaction and Margules fits are distinct. "
                    "Nineteen warmer paired observations are withheld from the interaction fit; Margules remains in-sample at every temperature. "
                    "Source uncertainty and global stability unqualified.",
                    fluid_class="Neutral mixtures",
                    property_class="VLE",
                )
                for r in part:
                    PLOTTED.append(
                        {
                            "panel": label,
                            "series": title,
                            "x": float(r["x_pentane"]),
                            "reference": float(r[reference]),
                            "model": float(r[model]),
                        }
                    )
            ax.legend(fontsize=5.1, ncol=2)
    finish(
        fig,
        "F9-neutral-mixtures.pdf",
        "Acetone/n-pentane: 372.7 K interaction calibration and 397.7/422.6 K withheld prediction; Margules fits each isotherm independently.",
    )


def gross_binaries():
    data = [
        r
        for r in csv_at("gross01", "results/figure-comparisons.csv")
        if r["kind"] != "density"
    ]
    csv_at("gross01", "inputs/source/figures-manifest.csv")
    systems = {
        "07": "Methane/n-butane",
        "08": "Ethane/n-decane",
        "09": "Ethane/diethyl ether",
        "10": "Propane/benzene",
        "12": "CO2/n-decane",
        "13": "CO2/cyclohexane",
        "15": "Hexane/1-chlorobutane",
    }
    fig, axes = canvas(7)
    for ax, number in zip(axes, sorted({r["figure"] for r in data})):
        selected = [r for r in data if r["figure"] == number]
        label = "F9e-" + number
        isobar = selected[0]["kind"] == "temperature"
        title = systems[number]
        panel(
            ax,
            label,
            title,
            REPRO + " / " + (CAL if number in {"12", "13", "15"} else PRED),
            "Liquid " + selected[0]["component"] + " mole fraction",
            "Temperature / K" if isobar else "Pressure / bar",
        )
        for j, temp in enumerate(
            sorted({r["temperature_c"] for r in selected}, key=float)
        ):
            for role, mark in [("model", "s"), ("experiment", "o")]:
                part = [
                    dict(
                        r,
                        plotted_reference=float(r["observed"])
                        + (273.15 if isobar else 0),
                        plotted_model=float(r["predicted"]) + (273.15 if isobar else 0),
                    )
                    for r in selected
                    if r["temperature_c"] == temp and r["reference"] == role
                ]
                if not part:
                    continue
                comparison(
                    ax,
                    part,
                    "x",
                    "plotted_model",
                    "plotted_reference",
                    ("Source model" if role == "model" else "Measurements")
                    + ("" if isobar else f", {temp} °C"),
                    COLORS[j],
                    mark,
                )
                record(
                    label + "-" + temp + "-" + role,
                    title,
                    "Temperature" if isobar else "Bubble pressure",
                    values(part, "temperature_k"),
                    values(part, "pressure_pa"),
                    len(part),
                    "Gross-Sadowski 2001 figure-specific interaction",
                    "caption-specific; original fit overlap unknown",
                    REPRO
                    if role == "model"
                    else CAL
                    if number in {"12", "13", "15"}
                    else PRED,
                    pairs(part, "plotted_model", "plotted_reference"),
                    note="Sparse digitized model curves and measured symbols are separate populations. Hydrocarbon zero-interaction predictions "
                    "differ from caption-specific carbon-dioxide/chlorobutane interactions. No original aggregate K-value statistic reproduced; "
                    "reading intervals are not measurement uncertainties.",
                    fluid_class="Neutral mixtures",
                    property_class="VLE",
                )
        ax.legend(fontsize=5.2, ncol=2)
    finish(
        fig,
        "F9-neutral-mixtures.pdf",
        "Gross binary comparisons use sparse source points with figure-specific interactions; source-model reproduction and physical comparison remain separate.",
    )


def associating():
    data = csv_at("gross02", "results/data/state-comparisons.csv")
    pure = [r for r in data if r["figure"] == "figure-01"]
    fig, axes = canvas(6)
    for i, compound in enumerate(sorted({r["names"] for r in pure})):
        for j, phase in enumerate(["liquid", "vapor"]):
            ax, label = axes[2 * i + j], f"F10a-{compound}-{phase}"
            selected = [
                r for r in pure if r["names"] == compound and r["series"] == phase
            ]
            panel(
                ax,
                label,
                compound + ", " + phase,
                REPRO + " / calibration",
                "Temperature / K",
                r"Density / kg m$^{-3}$",
            )
            for role, mark, color in [
                ("paper_pc_saft_curve", "s", BLUE),
                ("experimental_marker", "o", ORANGE),
            ]:
                part = [r for r in selected if r["dataset"] == role]
                if not part:
                    continue
                comparison(
                    ax,
                    part,
                    "T",
                    "value",
                    "reference",
                    "Source model" if role.startswith("paper") else "Measured markers",
                    color,
                    mark,
                )
                record(
                    label + "-" + role,
                    compound + " " + phase,
                    "Saturated density",
                    values(part, "T"),
                    values(part, "P_solved_Pa"),
                    len(part),
                    "Gross-Sadowski 2002 associating parameters",
                    "source calibration; exact overlap incomplete",
                    REPRO if role.startswith("paper") else CAL,
                    pairs(part, "value", "reference"),
                    "MAE / kg m-3" if role.startswith("paper") else "AARD / %",
                    note="Liquid and vapor statistics are separate. Obscured source-model coordinates are omitted, while measured markers retain their role. "
                    "The methanol marker at 505 K and 180 kg/m3 is vapor; original fit-population averages remain unavailable.",
                    fluid_class="Associating fluids",
                    property_class="Density",
                )
            ax.legend(fontsize=5.7)
    finish(
        fig,
        "F10-associating.pdf",
        "Associating pure-fluid densities: liquid and vapor remain separate; measured markers used for source calibration are labelled.",
    )
    # Facet mixed quantities and phases instead of pooling liquid-liquid and vapor-liquid statistics.
    for numbers, group_label in [
        ([f"figure-{i:02}" for i in range(2, 10)], "F10a"),
        (["figure-10"], "F10b"),
    ]:
        csv_at("gross02", "results/data/state-comparisons.csv")
        selected = [r for r in data if r["figure"] in numbers]
        keys = sorted({(r["figure"], r["names"], r["quantity"]) for r in selected})
        fig, axes = canvas(len(keys))
        for ax, (number, names, quantity) in zip(axes, keys):
            part = [
                r
                for r in selected
                if (r["figure"], r["names"], r["quantity"]) == (number, names, quantity)
            ]
            label = (
                group_label
                + "-"
                + number[-2:]
                + "-"
                + names.replace("/", "-")
                + "-"
                + quantity
            )
            xlabel = (
                "Temperature / K"
                if quantity == "x"
                else ("Water" if names == "water/1-pentanol" else names.split("/")[0].capitalize()) + " mole fraction"
            )
            ylabel = (
                "First-component liquid mole fraction"
                if quantity == "x"
                else "Pressure / bar"
                if quantity == "P_bar"
                else "Temperature / K"
            )
            panel(ax, label, names, REPRO + " / calibration", xlabel, ylabel)
            curves = sorted({(r["series"], r["mode"]) for r in part})
            for j, (series, mode) in enumerate(curves):
                for dataset in sorted({r["dataset"] for r in part}):
                    cohort = [
                        r
                        for r in part
                        if (r["series"], r["mode"], r["dataset"])
                        == (series, mode, dataset)
                    ]
                    if not cohort:
                        continue
                    is_model = dataset.startswith("paper_")
                    comparison(
                        ax,
                        cohort,
                        "T" if quantity == "x" else "x",
                        "value",
                        "reference",
                        f"{series.replace('_', ' ')} / "
                        + ("model" if is_model else "measured"),
                        COLORS[j % len(COLORS)],
                        "s" if is_model else "o",
                    )
                    record(
                        label + f"-{j}-" + dataset,
                        names,
                        ylabel,
                        values(cohort, "T_solved_K"),
                        values(cohort, "P_solved_Pa"),
                        1 if mode == "vlle" else len(cohort),
                        "Gross-Sadowski 2002 figure-specific interactions",
                        "source fits; supplements separate",
                        REPRO if is_model else CAL,
                        pairs(cohort, "value", "reference"),
                        "MAE / mole fraction"
                        if quantity == "x"
                        else "MAE / bar"
                        if quantity == "P_bar"
                        else "MAE / K",
                        note="Source model/correlation, experimental markers and supplemental compilations are distinct. "
                        "Calibration overlap prevents independent-validation claims. VLLE plotting coordinates repeat one coexistence state, "
                        "not four independent states; above-coexistence bubble discrepancies remain visible. Supplemental liquid-only water/pentanol states above coexistence can be metastable to vapor; primary lineage of these compilations is incomplete. No global stability certificate.",
                        fluid_class="Associating fluids",
                        property_class="LLE" if quantity == "x" else "VLE",
                    )
            if group_label == "F10b" and quantity == "T_K":
                vlle = next(r for r in part if r["mode"] == "vlle")
                xs = json.loads(vlle["compositions"])[::2]
                ax.plot(
                    xs,
                    [float(vlle["T_solved_K"])] * 3,
                    "*",
                    ms=6,
                    color="black",
                    label="One three-phase state",
                )
                ax.text(
                    0.02,
                    0.02,
                    f"Three-phase T = {float(vlle['T_solved_K']):.3f} K",
                    transform=ax.transAxes,
                    fontsize=6,
                )
                PLOTTED.append(
                    {
                        "panel": label,
                        "series": "one VLLE state",
                        "liquid_water_fractions": xs[:2],
                        "vapor_water_fraction": xs[2],
                        "T_K": float(vlle["T_solved_K"]),
                    }
                )
            handles, labels = ax.get_legend_handles_labels()
            labels = [
                s.replace("lle ", "LLE ").replace("1 pentanol", "1-pentanol")
                .replace(" rich", "-rich").replace(" lean", "-lean")
                if s.startswith("lle ") else s for s in labels
            ]
            ax.legend(handles, labels, fontsize=4.7, ncol=2)
        finish(
            fig,
            "F10-associating.pdf",
            "Associating binary and water/1-pentanol comparisons retain calibration overlap, supplemental-source limits, and one three-phase coexistence state.",
        )


def held_2008():
    data = csv_at("held08", "results/figure-comparisons.csv.gz")
    keys = sorted(
        {
            (r["figure"], r["panel"], r["salt"], r["temperature_K"], r["quantity"])
            for r in data
        }
    )
    assert len(keys) == 45 and len({r["salt"] for r in data}) == 19
    units = {
        "density_difference": r"Density difference / kg m$^{-3}$",
        "vapor_pressure": "Vapor pressure / bar",
        "water_activity_coefficient": r"Water coefficient $a_w/x_w$",
        "mean_ionic_activity_coefficient_molal": r"Molal ionic coefficient $\gamma_\pm$",
    }
    for start in range(0, len(keys), 9):
        fig, axes = canvas(min(9, len(keys) - start), 3)
        # Reading every page's exact source through the pinned reader retains its membership.
        csv_at("held08", "results/figure-comparisons.csv.gz")
        for ax, key in zip(axes, keys[start : start + 9]):
            part = sorted(
                [
                    r
                    for r in data
                    if tuple(
                        r[k]
                        for k in [
                            "figure",
                            "panel",
                            "salt",
                            "temperature_K",
                            "quantity",
                        ]
                    )
                    == key
                ],
                key=lambda r: float(r["molality"]),
            )
            number, _, salt, temp, quantity = key
            label = f"F11b-{start + list(axes).index(ax) + 1:02}"
            panel(
                ax,
                label,
                f"{salt}, {float(temp):g} K",
                REPRO + "; partial",
                r"Molality / mol kg$^{-1}$ water",
                units[quantity],
            )
            for role, marker, color in [("line", "s", BLUE), ("markers", "o", ORANGE)]:
                points = [
                    r for r in part if r["source_kind"] == role and r["qa"] == "pass"
                ]
                if not points:
                    continue
                comparison(
                    ax,
                    points,
                    "molality",
                    "native_value",
                    "value",
                    "Digitized source-model curve" if role == "line" else "Measurements",
                    color,
                    marker,
                    uncertainty="value_half_width",
                )
                record(
                    label + "-" + role,
                    salt + " in water",
                    quantity,
                    [float(temp)],
                    [p * 1e5 for p in values(points, "native_value")]
                    if quantity == "vapor_pressure"
                    else [1e5],
                    len(points),
                    "Held-Cameretti-Sadowski 2008; printed screening sigma",
                    "withheld RbCl (source Figure 8)"
                    if number == "8"
                    else "mostly original calibration",
                    REPRO if role == "line" else PRED if number == "8" else CAL,
                    pairs(points, "native_value", "value"),
                    "MAE / kg m-3"
                    if quantity == "density_difference"
                    else "MAE / bar"
                    if quantity == "vapor_pressure"
                    else "MAE / dimensionless",
                    note="Partial reproduction only; correlated model-ink samples are not independent observations or aggregate experimental accuracy. "
                    "Reading bars are digitization intervals. Unresolved source candidates remain excluded; disjoint intervals remain disagreements. "
                    "Fixed pressure assumed 1 bar; off-temperature permittivity conditional on the stated supplement. Original Table 4 population/weights unavailable.",
                    fluid_class="Aqueous electrolytes",
                    property_class="Density"
                    if quantity == "density_difference"
                    else "VLE"
                    if quantity == "vapor_pressure"
                    else "Activity / osmotic",
                )
            unresolved = [r for r in part if r["qa"] != "pass"]
            # Keep unresolved source locations visible without presenting them as admitted references.
            ax.plot(
                values(unresolved, "molality"),
                values(unresolved, "value"),
                "+",
                color="#777777",
                ms=3,
            )
            ax.legend(fontsize=5.1, loc="best")
            PLOTTED.extend(
                {
                    "panel": label,
                    "series": r["source_kind"],
                    "x": float(r["molality"]),
                    "reference": float(r["value"]),
                    "model": float(r["native_value"]),
                    "source_admitted": r["qa"] == "pass",
                    "outcome": r["outcome"],
                    "reference_half_width": float(r["value_half_width"]),
                    "native_interval_samples": json.loads(r["native_ordinates"]),
                }
                for r in part
            )
        finish(
            fig,
            "F11-aqueous-electrolytes.pdf",
            "Held 2008: all 19 salts and 45 series, partial reproduction. Unresolved source candidates and model-curve disagreements stay visible; no aggregate experimental accuracy follows.",
        )


def dielectric_panels():
    data = csv_at("dielectric", "results/comparison.csv")
    fig, axes = canvas(3)
    for ax, prop, role, label, ylabel in [
        (axes[0], "miac", "fit", "F11d-fit", r"Molal ionic coefficient $\gamma_\pm$"),
        (axes[1], "miac", "held", "F11d-held", r"Molal ionic coefficient $\gamma_\pm$"),
        (
            axes[2],
            "relative_permittivity",
            "reference",
            "F11d-epsilon",
            "Relative permittivity / 1",
        ),
    ]:
        panel(
            ax,
            label,
            "NaCl/water, 298.15 K, 1 bar",
            CAL if role == "fit" else PRED if role == "held" else REPRO,
            r"Molality / mol kg$^{-1}$ water",
            ylabel,
        )
        selected = [
            r
            for r in data
            if r["property"] == prop
            and (
                role == "reference"
                or (float(r["condition_value"]) <= 2) == (role == "fit")
            )
        ]
        for j, variant in enumerate(["F0", "F1", "D1", "P1"]):
            part = [r for r in selected if r["variant"] == variant]
            if j == 0:
                ax.plot(
                    values(part, "condition_value"),
                    values(part, "observed"),
                    "o",
                    mfc="none",
                    color="black",
                    ms=3,
                    label="Evaluated reference"
                    if prop == "miac"
                    else "Source-derived samples",
                )
            ax.plot(
                values(part, "condition_value"),
                values(part, "model"),
                ["x", "+", "1", "2"][j],
                color=COLORS[j],
                ms=4,
                label=variant,
            )
            evidence = (
                REPRO
                if prop == "relative_permittivity"
                else CAL
                if variant == "F0" or role == "fit"
                else PRED
            )
            record(
                label + "-" + variant,
                "NaCl/water",
                prop,
                [298.15],
                [1e5],
                len(part),
                variant + ": inherited parameters; one shared interaction",
                "all activity rows" if variant == "F0" else "activity molality <= 2",
                evidence,
                pairs(part, "model", "observed"),
                note="F0 used all activity rows; F1/D1/P1 train on 21 lower-molality rows, with eight higher rows withheld. "
                "Inherited Figiel/Rueben fits retain overlap. Seventeen permittivity samples lack verified primary-table and pointwise uncertainty. "
                "Hybrid-parameter sensitivity, not validation of a superior form; no measured-density comparison.",
                fluid_class="Aqueous electrolytes",
                property_class="Permittivity"
                if prop == "relative_permittivity"
                else "Activity / osmotic",
            )
            for r in part:
                PLOTTED.append(
                    {
                        "panel": label,
                        "series": variant,
                        "x": float(r["condition_value"]),
                        "reference": float(r["observed"]),
                        "model": float(r["model"]),
                    }
                )
        ax.legend(fontsize=5.5, ncol=2)
    finish(
        fig,
        "F11-aqueous-electrolytes.pdf",
        "NaCl permittivity forms: low-molality calibration and conditional higher-molality prediction; F0 has fitting overlap throughout. Permittivity primary table and pointwise uncertainty are unverified.",
    )


def propane_fit():
    data = csv_at("propane_fit", "results/propane-replay-2026-10-01.csv")
    train = set(
        values(
            csv_at("propane_fit", "inputs/source/propane-baseline.csv"), "temperature_k"
        )
    )
    assert train == {210, 270, 330}
    for r in data:
        r["selected_role"] = (
            "training"
            if float(r["T_K"]) in train
            else "stress"
            if r["role"] == "stress"
            else "held_out"
        )
    assert sum(r["selected_role"] == "held_out" for r in data) == 19
    assert (
        next(r for r in data if float(r["T_K"]) == 150)["selected_role"] == "held_out"
    )
    fig, axes = canvas(4)
    for j, (q, suffix, scale, ylabel) in enumerate(
        [
            ("p", "Pa", 1e6, "Saturation pressure / MPa"),
            ("rho", "kg_m3", 1, r"Liquid density / kg m$^{-3}$"),
        ]
    ):
        for residual in [False, True]:
            ax = axes[j + (2 if residual else 0)]
            label = "F15a-" + q + ("-deviation" if residual else "")
            panel(
                ax,
                label,
                "Propane: train 210/270/330 K; 19 withheld",
                CAL + " / " + PRED,
                "Temperature / K",
                "Relative deviation / %" if residual else ylabel,
            )
            for k, role in enumerate(["training", "held_out", "stress"]):
                part = [r for r in data if r["selected_role"] == role]
                model, reference, unc = (
                    f"{q}_fitted_{suffix}",
                    f"{q}_source_{suffix}",
                    f"{q}_expanded_uncertainty_{suffix}",
                )
                if residual:
                    ax.plot(
                        values(part, "T_K"),
                        [100 * (float(r[f"{q}_packet_{suffix}"]) / float(r[reference]) - 1) for r in part],
                        "+", color="#777777", ms=3,
                    )
                    ax.plot(
                        values(part, "T_K"),
                        [
                            100 * (float(r[model]) / float(r[reference]) - 1)
                            for r in part
                        ],
                        MARKS[k],
                        color=COLORS[k],
                        ms=3.5,
                        label=role.replace("_", " "),
                    )
                    ax.errorbar(
                        values(part, "T_K"),
                        [0] * len(part),
                        yerr=[100 * float(r[unc]) / float(r[reference]) for r in part],
                        fmt="none",
                        ecolor=COLORS[k],
                        lw=0.5,
                    )
                else:
                    comparison(
                        ax,
                        part,
                        "T_K",
                        model,
                        reference,
                        role.replace("_", " "),
                        COLORS[k],
                        MARKS[k],
                        uncertainty=unc,
                        scale=scale,
                    )
                    ax.plot(
                        values(part, "T_K"),
                        [float(r[f"{q}_packet_{suffix}"]) / scale for r in part],
                        "+",
                        color="#777777",
                        ms=3,
                    )
                record(
                    label + "-" + role,
                    "Propane",
                    ylabel.split(" / ")[0],
                    values(part, "T_K"),
                    values(part, "p_source_Pa"),
                    len(part),
                    "Three-state pure-parameter fit; published packet shown separately",
                    "210/270/330 K"
                    if role == "training"
                    else "no; withheld"
                    if role == "held_out"
                    else "no; stress",
                    CAL if role == "training" else PRED,
                    pairs(part, model, reference),
                    note="Training membership comes from the actual three-row fit input; 150 K is withheld, yielding 19 held-out states. "
                    "110/340 K are stress cases. Expanded measurement uncertainties are shown; fit weights are unit log-ratio weights. "
                    "Conditional covariance is not predictive uncertainty. Cold-pressure and withheld-density misses remain visible.",
                    fluid_class="Pure fluids",
                    property_class="Saturation" if q == "p" else "Density",
                )
            if not residual and q == "p":
                r150 = next(r for r in data if float(r["T_K"]) == 150)
                ax.annotate(
                    "150 K: withheld",
                    (150, float(r150["p_fitted_Pa"]) / 1e6),
                    xytext=(178, 1e-5),
                    fontsize=6,
                    arrowprops={"arrowstyle": "-", "lw": 0.5},
                )
            if residual:
                ax.axhline(0, color="black", lw=0.5)
            elif q == "p":
                ax.set_yscale("log")
            ax.plot([], [], "+", color="#777777", label="Published parameters")
            ax.legend(fontsize=5.5, ncol=2)
    finish(
        fig,
        "F15-propane-fit.pdf",
        "Worked propane fit: 210/270/330 K training, 150 K withheld, 19 held-out states and two stress states. Physical-property comparisons precede deviation diagnostics.",
    )


def reactive_amines():
    fig, axes = canvas(6)
    osmotic = csv_at(
        "wangler",
        "figures/figure-10/output/wangler-2018-mdeah-osmotic-plotted-data.csv",
    )
    pressures = csv_at(
        "wangler", "figures/figure-14/output/wangler-2018-table-10-native.csv"
    )
    for j, temp in enumerate(sorted(set(values(osmotic, "temperature_k")))):
        ax, label = axes[j], f"F16a-osmotic-{temp:g}"
        panel(
            ax,
            label,
            f"Wangler MDEA/HCl, {temp:g} K",
            CAL,
            r"Molality / mol kg$^{-1}$",
            "Osmotic coefficient / 1",
        )
        for k, lane in enumerate(["S0", "S1"]):
            part = [
                r
                for r in osmotic
                if r["lane"] == lane and float(r["temperature_k"]) == temp
            ]
            comparison(
                ax,
                part,
                "molality_mol_per_kg",
                "fitted_model",
                "experimental_osmotic_coefficient",
                lane,
                COLORS[k],
                MARKS[k],
            )
            record(
                label + "-" + lane,
                "MDEA/HCl/water",
                "Osmotic coefficient",
                [temp],
                [],
                len(part),
                "Wangler MDEA-water interaction fit, " + lane,
                "yes",
                CAL,
                pairs(part, "fitted_model", "experimental_osmotic_coefficient"),
                note="Ten calibration observations in total; S0/S1 repeat the same observations under separate declared conventions. "
                "The subsequent pressure prediction fails; calibration does not qualify absorption prediction.",
                fluid_class="Reactive amines",
                property_class="Activity / osmotic",
            )
        ax.legend(fontsize=6)
    ax = axes[2]
    panel(
        ax,
        "F16a-pressure",
        "Wangler: failed pressure transfer, 303 K",
        PRED,
        "CO2 loading / mol per mol MDEA",
        "Total pressure / bar",
    )
    for k, lane in enumerate(["fitted-S0", "fitted-S1", "printed-S0"]):
        part = [r for r in pressures if r["lane"] == lane]
        if k == 0:
            ax.errorbar(
                values(part, "co2_loading_mol_per_mol_mdea"),
                values(part, "experimental_total_pressure_bar"),
                yerr=values(part, "pressure_standard_uncertainty_bar"),
                fmt="o",
                color="black",
                mfc="none",
                ms=4,
                lw=0.5,
                label="Withheld observations",
            )
        ax.plot(
            values(part, "co2_loading_mol_per_mol_mdea"),
            values(part, "native_pressure_bar"),
            ["x", "+", "1"][k],
            color=COLORS[k],
            ms=5,
            label=lane + (" diagnostic" if k == 2 else ""),
        )
        record(
            "F16a-pressure-" + lane,
            "CO2 in 60 wt% MDEA",
            "Total pressure",
            [303],
            [p * 1e5 for p in values(part, "experimental_total_pressure_bar")],
            len(part),
            lane,
            "no pressure fit; osmotic calibration"
            if k < 2
            else "printed-parameter diagnostic",
            PRED if k < 2 else REPRO,
            pairs(part, "native_pressure_bar", "experimental_total_pressure_bar"),
            note="Failed withheld pressure prediction under the declared reading. Root cause unresolved despite investigated declaration/convention probes; "
            "no universal MDEA or MEA conclusion. Printed S0 is a separate diagnostic. Pressure bars are source standard uncertainties.",
            fluid_class="Reactive amines",
            property_class="VLE",
        )
        for r in part:
            PLOTTED.append(
                {
                    "panel": "F16a-pressure",
                    "series": lane,
                    "x": float(r["co2_loading_mol_per_mol_mdea"]),
                    "reference": float(r["experimental_total_pressure_bar"]),
                    "model": float(r["native_pressure_bar"]),
                }
            )
    ax.set_yscale("log")
    ax.legend(fontsize=5.5)
    data = csv_at("cleeton", "results/data/cleeton-2020-dicko-rows.csv")
    curves = [
        r
        for r in csv_at("cleeton", "results/data/cleeton-2020-figure-16-gate.csv")
        if r["lane"] == "main-0.88sigma" and r["reading"] == "caption"
    ]
    assert (
        len(curves) == 30 and sum(r["within_tolerance"] == "False" for r in curves) == 8
    )
    for ax, gas, key in [(axes[3], "co2", "left"), (axes[4], "h2s", "right")]:
        label = "F16b-" + gas
        panel(
            ax,
            label,
            f"Cleeton {gas.upper()}, 323.15 K, 50 wt% MDEA",
            REPRO,
            gas.upper() + " loading / mol per mol MDEA",
            "Partial pressure / kPa",
        )
        part = [
            r
            for r in data
            if r["lane"] == "main-0.88sigma"
            and r[
                f"figure_16_{key}_series_"
                + ("h2s" if gas == "co2" else "co2")
                + "_loading"
            ]
        ]
        other = "h2s" if gas == "co2" else "co2"
        group_col = f"figure_16_{key}_series_{other}_loading"
        for j, loading in enumerate(sorted({r[group_col] for r in part}, key=float)):
            cohort = [
                dict(
                    r,
                    measured_kPa=float(r[f"{gas}_pressure_mpa"]) * 1000,
                    calculated_kPa=float(r[f"p_{gas}_pa"]) / 1000,
                )
                for r in part
                if r[group_col] == loading
            ]
            comparison(
                ax,
                cohort,
                gas + "_loading",
                "calculated_kPa",
                "measured_kPa",
                f"{other.upper()} loading {loading}",
                COLORS[j],
                MARKS[j],
            )
        record(
            label + "-measured",
            "CO2/H2S in MDEA",
            gas.upper() + " partial pressure",
            [323.15],
            [float(r[f"{gas}_pressure_mpa"]) * 1e6 for r in part],
            len(part),
            "Cleeton printed-model replay; 0.88 sigma",
            "inherited pure/binary fits; no Dicko refit",
            REPRO,
            [
                (float(r[f"p_{gas}_pa"]), float(r[f"{gas}_pressure_mpa"]) * 1e6)
                for r in part
            ],
            note="Bounded physical comparison, with large pressure discrepancies. Methane is omitted although the source states contain it. "
            "Source code/supporting information unavailable; inherited calibration and incomplete component populations prevent an absorption-accuracy claim.",
            fluid_class="Reactive amines",
            property_class="VLE",
        )
        ax.set_yscale("log")
        ax.legend(fontsize=5.2, ncol=2)
    ax = axes[5]
    panel(
        ax,
        "F16b-curves",
        "Cleeton: 8 of 30 source-curve checks miss",
        REPRO,
        "Acid-gas loading / mol per mol MDEA",
        r"$\ln(P_{model}/P_{source\ curve})$",
    )
    for j, gas in enumerate(["co2", "h2s"]):
        part = [r for r in curves if r["gas"] == gas]
        errors = [
            math.log(float(r[f"p_{gas}_pa"]) / float(r["figure_pressure_pa"]))
            for r in part
        ]
        ax.errorbar(
            values(part, "loading"),
            errors,
            yerr=values(part, "tolerance"),
            fmt=MARKS[j],
            color=COLORS[j],
            ms=3,
            lw=0.5,
            mfc="none",
            label=gas.upper() + "; bars: own reading criterion",
        )
        record(
            "F16b-curves-" + gas,
            "CO2/H2S in MDEA",
            gas.upper() + " source-curve pressure",
            [323.15],
            values(part, "figure_pressure_pa"),
            len(part),
            "Cleeton caption reading; 0.88 sigma",
            "inherited fits",
            REPRO,
            pairs(part, f"p_{gas}_pa", "figure_pressure_pa"),
            "max |ln pressure ratio|",
            note="Thirty primary caption-reading comparisons retain eight misses against their own slope-aware reading criteria. "
            "Alternative group-loading interpretation and unreproduced source terms remain limitations; not accurate absorption prediction.",
            fluid_class="Reactive amines",
            property_class="VLE",
        )
        PLOTTED.extend(
            {
                "panel": "F16b-curves",
                "series": gas,
                "x": float(r["loading"]),
                "reference": float(r["figure_pressure_pa"]),
                "model": float(r[f"p_{gas}_pa"]),
                "log_reading_criterion": float(r["tolerance"]),
                "within_criterion": r["within_tolerance"] == "True",
            }
            for r in part
        )
    ax.axhline(0, color="black", lw=0.5)
    ax.legend(fontsize=5.3)
    finish(
        fig,
        "F16-reactive-amines.pdf",
        "Wangler osmotic calibration followed by failed pressure transfer is a limit under the declared MDEA reading, not an MEA result. Cleeton is a bounded published-model replay with eight curve misses and large physical discrepancies; methane omission remains explicit.",
    )


def settler():
    data = json_at("settler", "results/live/summary.json")
    identity = json_at("settler", "results/live/identity.json")
    assert data["identity"]["record_sha256"] == identity["record_sha256"]
    encoded = retained("settler", "results/live/record.netstrings")
    assert hashlib.sha256(encoded.encode()).hexdigest() == identity["record_sha256"]
    tokens = []
    while encoded:
        count, _, body = encoded.partition(":")
        size = int(count)
        assert body[size] == ","
        tokens.append(body[:size])
        encoded = body[size + 1 :]
    temperature = float.fromhex(tokens[tokens.index("temperature") + 1])
    pressure = float.fromhex(tokens[tokens.index("pressure") + 1])
    assert (temperature, pressure) == (293.15, 100000)
    runs, baseline = data["optimization"], data["baseline"]
    for run in runs:
        assert math.isclose(run["objective"], run["si"] + run["ss"], rel_tol=1e-14)
        assert math.isclose(
            run["recovery"],
            run["flows"][6] / (run["flows"][1] + run["flows"][6]),
            rel_tol=1e-13,
        )
    fig, axes = canvas(2)
    ax = axes[0]
    panel(
        ax,
        "F17a-flows",
        "293.15 K, 100 kPa: local two-liquid branch",
        NV,
        "Component",
        r"Outlet molar flow / mol s$^{-1}$",
    )
    species = ["Water", "Ethanol", "Isobutanol", "Na+", "Cl−"]
    for j, (name, run) in enumerate(
        [("Baseline", baseline), ("Local optimum", runs[0])]
    ):
        for phase, color, marker in [(0, BLUE, "o"), (1, ORANGE, "s")]:
            ax.plot(
                [i + (j - 0.5) * 0.15 for i in range(5)],
                run["flows"][phase * 5 : phase * 5 + 5],
                marker,
                color=color,
                mfc="none" if j == 0 else color,
                ms=4,
                label=name + (", aqueous" if phase == 0 else ", organic"),
            )
            PLOTTED.append(
                {
                    "panel": "F17a-flows",
                    "series": name,
                    "phase": "aqueous" if phase == 0 else "organic",
                    "components": species,
                    "flows_mol_s": run["flows"][phase * 5 : phase * 5 + 5],
                }
            )
    ax.set(xticks=range(5), xticklabels=species, yscale="log")
    ax.legend(fontsize=5.6, ncol=2)
    record(
        "F17a-flows",
        "Water/ethanol/isobutanol/NaCl settler",
        "Outlet flows; material balance",
        [293.15],
        [1e5],
        10,
        "Khudaida source parameters; declared feed map",
        "no process fit",
        NV,
        [(r["material"], 0) for r in runs],
        "max material-balance residual / mol s-1",
        note="Fixed temperature/pressure and declared local two-liquid branch; no measured settler data, physical validation, global optimum or reactive qualification. "
        "Three starts share one selected feed box, 0.85 to 1.15 for both controls.",
        fluid_class="Electrolyte LLE",
        property_class="LLE",
    )
    ax = axes[1]
    panel(
        ax,
        "F17b-local",
        "Three starts: numerical local optimum",
        NV,
        "Start control multiplier",
        "Objective / 1",
    )
    for start, run in zip([0.9, 1, 1.1], runs):
        ax.plot(start, run["si"] + run["ss"], "x", color=BLUE, ms=6)
        PLOTTED.append(
            {
                "panel": "F17b-local",
                "start": start,
                "objective": run["objective"],
                "recovery": run["recovery"],
                "stationarity": run["stationarity"],
                "si": run["si"],
                "ss": run["ss"],
            }
        )
    maximum = max(r["stationarity"] for r in runs)
    assert math.isclose(maximum, 3.665966266010231e-12, rel_tol=1e-14)
    ax.text(
        0.03,
        0.95,
        f"Baseline objective: {baseline['objective']:.6f}\nLocal objective: {runs[0]['objective']:.6f}\nEthanol recovery: {runs[0]['recovery']:.6f}\nMax stored stationarity: {maximum:.4e}",
        transform=ax.transAxes,
        va="top",
        fontsize=6.7,
    )
    ax.ticklabel_format(axis="y", useOffset=False)
    ax.set_ylim(1.9755, 1.9762)
    record(
        "F17b-local",
        "Water/ethanol/isobutanol/NaCl settler",
        "Stored stationarity norm",
        [293.15],
        [1e5],
        len(runs),
        "Khudaida source parameters; declared local branch",
        "no process fit",
        NV,
        [(r["stationarity"], 0) for r in runs],
        "max stored stationarity norm",
        note="Stored stationarity norm for the selected numerical solution. Objective reduced from retained controls and recovery from retained ethanol flows. "
        "No independent re-evaluation of the full nonlinear program or global optimum certificate.",
        fluid_class="Electrolyte LLE",
        property_class="LLE",
    )
    finish(
        fig,
        "F17-settler.pdf",
        "Numerical settler use at 293.15 K and 100 kPa on a declared local two-liquid branch, controls bounded 0.85–1.15. Stored stationarity is shown for the selected run; no physical-validation or global-optimum claim.",
    )


def hda_hydrogen_pairs():
    data = [
        r
        for r in csv_at("hda", "results/hydrogen_pair_rows.csv")
        if r["partition"] != "stress"
    ]
    fit = json_at("hda", "results/hydrogen_pair_fit.json")
    fig, axes = canvas(2)
    for ax, pair, letter in zip(axes, ["hydrogen/benzene", "hydrogen/methane"], "ab"):
        label = f"F18{letter}-{pair.split('/')[1]}"
        k = fit[pair]["k_hat"]
        panel(
            ax,
            label,
            f"{pair.capitalize()} bubble pressure, fitted $k_{{ij}}$ = {k:.3f}",
            f"{CAL} (training) and {PRED} (withheld)",
            "Measured pressure / MPa",
            "Model pressure / MPa",
        )
        subset = [r for r in data if r["pair"] == pair]
        top = max(float(r["P_obs_Pa"]) for r in subset) / 1e6 * 1.08
        ax.plot([0, top], [0, top], color="black", lw=0.6)
        for part, color, mark, fitted in [
            ("training", BLUE, "o", "fitted"),
            ("holdout", ORANGE, "s", "withheld"),
        ]:
            rows_part = [r for r in subset if r["partition"] == part]
            operands = pairs(rows_part, "khat_P_Pa", "P_obs_Pa")
            ax.plot(
                [b / 1e6 for a, b in operands],
                [a / 1e6 for a, b in operands],
                mark,
                color=color,
                mfc="none",
                ms=3,
                mew=0.7,
                label=f"{fitted.capitalize()} rows ({len(operands)})",
            )
            PLOTTED.extend(
                {"panel": label, "series": part, "reference": b, "model": a}
                for a, b in operands
            )
            record(
                f"{label}-{part}",
                pair.capitalize(),
                "Bubble pressure",
                [float(r["T_K"]) for r in rows_part],
                [float(r["P_obs_Pa"]) for r in rows_part],
                len(rows_part),
                "One fitted hydrogen pair coefficient",
                "yes" if part == "training" else "no (withheld)",
                CAL if part == "training" else PRED,
                operands,
                note="Fixed-composition bubble pressure at the measured liquid composition. "
                + (
                    "The methane pair keeps a large discrepancy after fitting; hot HDA states lie outside the binary source domains."
                    if pair == "hydrogen/methane"
                    else "Hot HDA states lie outside the binary source domain."
                ),
                fluid_class="Neutral mixtures",
                property_class="VLE",
            )
        ax.set(xlim=(0, top), ylim=(0, top))
        ax.legend(fontsize=6, loc="upper left")
    finish(
        fig,
        "F18-hda-hydrogen-pairs.pdf",
        "Hydrogen pair coefficients for a two-flash hydrodealkylation example: bubble pressures at fitted and withheld measured states. Discrepancy allowances are not statistical uncertainty; hot process states are outside the binary source domains.",
    )


def discrete_equilibrium():
    base = "results/issue-224-sliceB/reference/"
    compared = csv_at("opensource", base + "comparisons.csv")
    ref = csv_at("opensource", base + "envelope-points.csv")
    critical = csv_at("opensource", base + "critical-binary.csv")
    pure = csv_at("opensource", base + "critical-pure.csv")
    fig, axes = canvas(3)
    ax = axes[0]
    panel(
        ax,
        "F4b",
        "Methane/ethane: nine discrete bubble/dew states",
        NV,
        "Temperature / K",
        "Pressure / MPa",
    )
    for j, bubble in enumerate(["1", "0"]):
        part = [r for r in ref if r["bubble"] == bubble]
        mode = "bubble" if bubble == "1" else "dew"
        ax.plot(
            values(part, "T_K"),
            [p / 1e6 for p in values(part, "P_Pa")],
            MARKS[j],
            mfc="none",
            ms=5,
            color=COLORS[j],
            label="FeOs " + mode,
        )
        for r in part:
            # The retained scalar keys preserve both distinct retrograde dew points near 268 K.
            prefix = f"{mode}/T={r['T_K']}/"
            p = next(
                c
                for c in compared
                if c["falsifier"] == "E2" and c["quantity"] == prefix + "P_Pa"
            )
            ax.plot(
                float(r["T_K"]), float(p["engine"]) / 1e6, "x", color=COLORS[j], ms=3
            )
            PLOTTED.append(
                {
                    "panel": "F4b",
                    "series": mode,
                    "x": float(r["T_K"]),
                    "reference": float(r["P_Pa"]),
                    "model": float(p["engine"]),
                }
            )
    part = [r for r in compared if r["falsifier"] == "E2"]
    record(
        "F4b",
        "Equal-feed methane/ethane",
        "Nine discrete bubble/dew states; four scalar quantities each",
        values(ref, "T_K"),
        values(ref, "P_Pa"),
        len(ref),
        "Gross 2001 methane/ethane, zero interaction",
        "no",
        NV,
        [
            (abs(float(r["engine"]) - float(r["reference"])) / float(r["scale"]), 0)
            for r in part
        ],
        "max scaled difference",
        note="Discrete numerical comparisons only: pressure, two densities and incipient composition. Nine states, not 36 physical observations. "
        "No continuous envelope or global-stability claim follows from these reference points.",
        fluid_class="Numerical methods",
        property_class="VLE",
    )
    ax.plot([], [], "x", color="black", label="Engine, same states")
    ax.legend(fontsize=5.5)
    ax = axes[1]
    panel(
        ax,
        "F4c",
        "Binary critical points: source-specific compositions",
        NV,
        "Methane mole fraction",
        "Critical temperature / K",
    )
    for j, source in enumerate(["feos", "teqp"]):
        part = [
            dict(
                r,
                engine=next(
                    c["engine"]
                    for c in compared
                    if c["falsifier"] == "C3"
                    and c["quantity"] == f"{source}/z0={r['z0']}/T_K"
                ),
            )
            for r in critical
            if r["source"] == source
        ]
        comparison(
            ax,
            part,
            "z0",
            "engine",
            "T_K",
            "FeOs fixed composition" if source == "feos" else "teqp polished locus",
            COLORS[j],
            MARKS[j],
        )
    for group, selected, temperatures, pressures, label in [
        ("C3", critical, values(critical, "T_K"), values(critical, "P_Pa"), "F4c"),
        ("C1", pure, values(pure, "T_K"), values(pure, "P_Pa"), "F4d"),
    ]:
        part = [r for r in compared if r["falsifier"] == group]
        record(
            label,
            "Methane/ethane" if group == "C3" else "Propane",
            "Critical temperature, density and pressure",
            temperatures,
            pressures,
            len(selected),
            "Gross 2001 pure/binary declarations",
            "no",
            NV,
            [
                (abs(float(r["engine"]) - float(r["reference"])) / float(r["scale"]), 0)
                for r in part
            ],
            "max scaled difference",
            note="FeOs uses fixed binary composition; teqp uses its own polished-locus compositions, which differ. "
            "Critical temperature, density and pressure are compared at each declared composition. Numerical verification, no physical validation or global-stability certificate.",
            fluid_class="Numerical methods",
            property_class="Critical point",
        )
    ax.legend(fontsize=5.5)
    ax = axes[2]
    panel(
        ax,
        "F4d",
        "Propane critical-point scalar comparisons",
        NV,
        "Critical quantity",
        "Absolute scaled difference",
    )
    for j, source in enumerate(["feos", "teqp"]):
        part = [
            r
            for r in compared
            if r["falsifier"] == "C1" and r["quantity"].startswith(source + "/")
        ]
        ax.plot(
            [i + (j - 0.5) * 0.12 for i in range(3)],
            [
                abs(float(r["engine"]) - float(r["reference"])) / float(r["scale"])
                for r in part
            ],
            MARKS[j],
            mfc="none",
            color=COLORS[j],
            ms=4,
            label=source,
        )
        PLOTTED.extend(
            {
                "panel": "F4d",
                "series": source,
                "quantity": r["quantity"].split("/")[-1],
                "model": float(r["engine"]),
                "reference": float(r["reference"]),
                "scale": float(r["scale"]),
            }
            for r in part
        )
    ax.set(
        xticks=[0, 1, 2],
        xticklabels=["Temperature", "Density", "Pressure"],
        yscale="log",
    )
    ax.legend(fontsize=6)
    finish(
        fig,
        "F4-phase-envelope.pdf",
        "Added discrete FeOs/teqp critical and FeOs bubble/dew comparisons are numerical verification at source-specific compositions. Only the separately retained native continuation supports lines; these points imply no continuous envelope or global stability.",
    )


def coverage_map():
    from matplotlib.patches import Patch, Rectangle

    classes = [
        "Pure fluids",
        "Neutral mixtures",
        "Associating fluids",
        "Reactive amines",
        "Aqueous electrolytes",
        "Mixed brines",
        "Mixed solvents",
        "Solvation",
        "Electrolyte LLE",
        "Reacting liquids",
        "Ideal reacting gas",
    ]
    props = [
        "Saturation",
        "Density",
        "VLE",
        "Activity / osmotic",
        "Solvation / transfer",
        "LLE",
        "pH / partition",
        "Reaction amounts",
        "Calorics",
        "Permittivity",
    ]
    evidence = [CAL, PRED, REPRO, NV]
    shades = dict(zip(evidence, [ORANGE, BLUE, GREEN, PURPLE]))
    short = {CAL: "C", PRED: "P", REPRO: "R", NV: "N"}
    fig, ax = plt.subplots(figsize=(7.1, 5.1))
    cells = []
    for i, cls in enumerate(classes):
        for j, prop in enumerate(props):
            entries = [
                r
                for r in TABLE
                if r["fluid_class"] == cls and r["property_class"] == prop
            ]
            present = [e for e in evidence if any(r["evidence"] == e for r in entries)]
            for k, e in enumerate(present):
                ax.add_patch(
                    Rectangle(
                        (j - 0.5 + k / len(present), i - 0.5),
                        1 / len(present),
                        1,
                        facecolor=shades[e],
                        alpha=0.55,
                        edgecolor="white",
                        lw=0.5,
                    )
                )
            if present:
                ax.text(
                    j,
                    i,
                    "/".join(short[e] for e in present),
                    ha="center",
                    va="center",
                    fontsize=8,
                )
            cells.append(
                {
                    "fluid_class": cls,
                    "property": prop,
                    "evidence": present,
                    "panels": [r["panel"] for r in entries],
                }
            )
    ax.set(
        xlim=(-0.5, len(props) - 0.5),
        ylim=(len(classes) - 0.5, -0.5),
        xticks=range(len(props)),
        yticks=range(len(classes)),
    )
    ax.set_xticklabels(props, rotation=35, ha="right", fontsize=7)
    ax.set_yticklabels(classes, fontsize=8)
    ax.grid(False)
    ax.set_title("Evidence class by fluid class and property", fontsize=10)
    ax.legend(
        handles=[
            Patch(color=shades[e], alpha=0.55, label=short[e] + ": " + e)
            for e in evidence
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, 1.04),
        ncol=2,
        fontsize=7,
    )
    fig.tight_layout()
    fig.savefig(
        OUT / "F7-coverage.pdf", metadata={"CreationDate": None, "ModDate": None}
    )
    plt.close(fig)
    (OUT / "F7-coverage.json").write_text(json.dumps(cells, indent=2) + "\n")
    FIGURES["F7-coverage.pdf"] = {
        "panels": [],
        "sources": [],
        "notes": "Derived only from T10 panel evidence. Split cells retain multiple classes. Blank means no selected panel, not absence of Engine capability. Numerical methods are excluded from fluid-class coverage.",
    }


def tex_escape(text):
    return "".join(
        {
            "&": r"\&",
            "%": r"\%",
            "_": r"\_\allowbreak{}",
            "/": r"/\allowbreak{}",
            "-": r"-\allowbreak{}",
            "#": r"\#",
            "$": r"\$",
            "{": r"\{",
            "}": r"\}",
            "\\": r"\textbackslash{}",
        }.get(c, c)
        for c in str(text)
    )


def number_range(numbers, scale=1):
    if not numbers:
        return "unreported"
    lo, hi = min(numbers) / scale, max(numbers) / scale
    return f"{lo:.6g}" if lo == hi else f"{lo:.6g} to {hi:.6g}"


def write_catalog():
    for pdf in PDFS.values():
        pdf.close()
    assert len({r["panel"] for r in TABLE}) == len(TABLE)
    assert all(
        r["points"] > 0 and r["comparisons"] > 0 and r["sources"] for r in TABLE
    ), [r["panel"] for r in TABLE if not r["sources"]]
    # Population checks from the promoted reviews prevent calibration/holdout and availability drift.
    by_panel = {r["panel"]: r for r in TABLE}
    assert by_panel["F12a-fit"]["points"] == 443
    assert by_panel["F12a-Phi"]["points"] + by_panel["F12a-aw"]["points"] == 72
    assert sum(r["points"] for r in TABLE if r["panel"].startswith("F13b-")) == 409
    assert all(r["comparisons"] == 184 for r in TABLE if r["panel"].startswith("F11a-"))
    assert by_panel["F14a-ties"]["comparisons"] == 64
    (OUT / "T10-master-table.json").write_text(
        json.dumps(
            [{k: v for k, v in r.items() if k != "sources"} for r in TABLE],
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )
    columns = [
        "panel",
        "system",
        "property",
        "T_K",
        "P_MPa",
        "points",
        "parameter_source",
        "fitted",
        "metric",
        "value",
        "evidence",
    ]
    with (OUT / "T10-master-table.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        for r in TABLE:
            writer.writerow(
                {
                    **{k: r[k] for k in columns if k not in ["T_K", "P_MPa"]},
                    "T_K": number_range(r["T_K"]),
                    "P_MPa": number_range(r["P_Pa"], 1e6),
                }
            )
    lines = [
        r"% Generated by scripts/figures.py. Requires longtable, pdflscape and array.",
        r"\begin{landscape}",
        r"\begingroup\scriptsize\setlength{\tabcolsep}{2pt}",
        r"\newcolumntype{L}[1]{>{\raggedright\arraybackslash}p{#1}}",
        r"\begin{longtable}{L{1.9cm}L{2.5cm}L{2.1cm}L{1.6cm}L{1.6cm}rL{2.4cm}L{1.9cm}L{2.2cm}rL{1cm}}",
        r"\caption{Selected scientific comparisons. Counts refer to source rows or retained numerical samples; repeated model evaluations do not add observations.}\label{tab:T10}\\",
        r"\toprule Panel & System & Property & $T$/K & $P$/MPa & $N$ & Parameter source & Fitted to data & Metric & Value & Class \\ \midrule\endfirsthead",
        r"\toprule Panel & System & Property & $T$/K & $P$/MPa & $N$ & Parameter source & Fitted to data & Metric & Value & Class \\ \midrule\endhead",
    ]
    codes = {CAL: "CAL", PRED: "PRED", REPRO: "REPRO", NV: "NV"}
    for r in TABLE:
        cells = [
            r["panel"],
            r["system"],
            r["property"],
            number_range(r["T_K"]),
            number_range(r["P_Pa"], 1e6),
            r["points"],
            r["parameter_source"],
            r["fitted"],
            r["metric"],
            f"{r['value']:.6g}",
            codes[r["evidence"]],
        ]
        lines.append(" & ".join(tex_escape(v) for v in cells) + r" \\")
    lines += [
        r"\bottomrule\end{longtable}",
        (
            r"CAL: calibration; PRED: prediction; REPRO: published-model reproduction; NV: numerical verification. "
            r"Model-line sample counts are not measured sample counts. Unreported conditions are left explicit. "
            r"The JSON companion records exact operands, comparison counts, selection limits and comparison limits."
        ),
        r"\endgroup\end{landscape}",
        "",
    ]
    (OUT / "T10-master-table.tex").write_text("\n".join(lines))
    provenance = {
        key: {k: v for k, v in f.items() if k != "data"} for key, f in FILES.items()
    }
    (OUT / "sources.json").write_text(json.dumps(provenance, indent=2) + "\n")
    (OUT / "figure-notes.json").write_text(
        json.dumps(
            {
                name: {k: v for k, v in info.items() if k != "sources"}
                for name, info in FIGURES.items()
            },
            indent=2,
        )
        + "\n"
    )
    (OUT / "plotted-comparisons.json").write_text(
        json.dumps(PLOTTED, indent=2, allow_nan=False) + "\n"
    )
    readme = [
        "# Paper figures from promoted retained results",
        "",
        "Regenerate from the repository root with `nice -n 10 .venv/bin/python docs/scientific/latex/scripts/figures.py`.",
        "",
        (
            "The script is the single figure owner. It reads CSV and JSON Git objects at the investigator-promoted commits from "
            "[the first owner promotion](https://github.com/tannerpolley/ePC-SAFT/issues/501#issuecomment-6373385229) and "
            "[the second owner promotion](https://github.com/tannerpolley/ePC-SAFT/issues/501#issuecomment-6376503363) and "
            "[the cost promotion](https://github.com/tannerpolley/ePC-SAFT/issues/681#issuecomment-6388971415). "
            "It does not import the Engine or calculate a new thermodynamic state. Every new PDF is vector output at 7.1-inch two-column width; dense galleries use multiple pages. "
            "The four earlier PDFs retain their original drawing functions and dimensions."
        ),
        "",
        (
            "`T10-master-table.tex` requires `longtable`, `pdflscape` and `array`; it is generated for later inclusion and does not change manuscript text. "
            "`T10-master-table.json` retains full-precision metrics, comparison populations, conditions and limits. "
            "A visual panel may have separate table rows for parameter variants, reference roles and training/withheld populations; those rows do not add independent observations. "
            "`sources.json` gives exact promoted commits, paths and SHA-256 values. "
            "`figure-notes.json` gives panel membership and caption limits; source membership is recorded here. "
            "`plotted-comparisons.json` retains ordinary paired series; categorical, projected and algorithm series are read directly from the exact source files."
        ),
        "",
        (
            "AARD is the mean absolute relative difference times 100; RMS is the square root of the mean squared direct difference; "
            "MAE is the mean absolute direct difference. Maxima state their quantity and unit. "
            "Missing values stay outside arithmetic and inside coverage counts. No unavailable value becomes a zero. "
            "Metrics against digitized model lines describe reproduction, not measured accuracy. "
            "In mixed-convention or mixed-declaration panels the metric includes both, as identified in each table note."
        ),
        "",
    ]
    for name, info in FIGURES.items():
        readme += [
            f"## {name}",
            "",
            "Panels: " + (", ".join(info["panels"]) or "fluid class by property") + ".",
            "",
            info["notes"],
            "",
        ]
        grouped = {}
        for key in info["sources"]:
            f = FILES[key]
            grouped.setdefault((f["commit"], str(Path(f["path"]).parent)), []).append(
                Path(f["path"]).name
            )
        for (commit, directory), names in grouped.items():
            if len(names) > 30:
                readme.append(
                    f"- `{directory}/*-state-*.json` at `{commit}` (exact files in `sources.json`)."
                )
            else:
                readme.append(
                    f"- `{directory}/` at `{commit}`: "
                    + ", ".join(f"`{n}`" for n in names)
                    + "."
                )
        for r in TABLE:
            if r["figure"] == name:
                readme.append(f"- {r['panel']}: **{r['evidence']}**; {r['note']}")
        readme.append("")
    readme += [
        "## Existing figures",
        "",
        (
            "- `derivatives.pdf`: `2026-shared-derivative-qualification/results/pressure-actions.csv` and "
            "`2026-head-to-head-comparisons/results/main-8298c602/h2h-sensitivities.csv`; numerical verification."
        ),
        "- `electrolyte-lle.pdf`: `2026-head-to-head-comparisons/results/h2h-khudaida-lle.csv`; numerical verification.",
        "- `acetic-acid-ph.pdf`: unchanged original at `12d0ba46b2d4db5731d3c827bea36411f9e5f5f6`; its drawing function is reused with the promoted inputs in F14.",
        (
            "- `chloride-calibration.pdf`: `2026-figiel-divalent-chlorides/results/stage1-minimax-rows.csv` and "
            "`figiel-1to1-comparison.csv`; calibration. Existing figures retain their earlier authority and are not newly promoted studies."
        ),
        "",
        (
            "Hubach and Rezaee are promoted only for the later limits section and are not added as physical-validation panels. "
            "The ten studies in the second owner promotion are included at their exact promoted commits. "
            "HDA two-flash is not promoted and is excluded. F4 separates retained continuation from discrete independent-code comparisons."
        ),
        "",
    ]
    (OUT / "README.md").write_text("\n".join(readme))
    print(
        f"Generated {len(FIGURES)} new figures and {len(TABLE)} T10 rows from {len(FILES)} pinned retained files."
    )


def promoted_figures():
    pure_fluids()
    neutral_mixtures()
    gross_binaries()
    associating()
    aqueous_electrolytes()
    dielectric_panels()
    propane_fit()
    reactive_amines()
    settler()
    hda_hydrogen_pairs()
    mixed_brines()
    mixed_solvents()
    reacting_lle()
    algorithm_figures()
    discrete_equilibrium()
    coverage_map()
    computational_cost()
    write_catalog()


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    derivatives()
    electrolyte_lle()
    acetic_acid_ph()
    chloride_calibration()
    assert math.isclose(
        max(
            abs(float(r["error"]))
            for r in rows(
                A / "2026-figiel-divalent-chlorides/results/stage1-minimax-rows.csv"
            )
            if r["set"] == "T"
        ),
        0.019370073426250307,
        rel_tol=1e-12,
    )

    promoted_figures()

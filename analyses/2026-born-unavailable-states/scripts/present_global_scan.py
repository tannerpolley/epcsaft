"""Plot retained fixed-target chart points; never calculate thermodynamics."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update(
    {"font.size": 9, "axes.titlesize": 11, "axes.labelsize": 9, "figure.titlesize": 13}
)
HERE = Path(__file__).resolve().parents[1]
fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
plotted = []
for case in (1, 2):
    sources = [
        p / "attempts.csv"
        for p in (HERE / "results/global-scan").iterdir()
        if p.name == f"case{case}" or p.name.startswith(f"case{case}-completion")
    ]
    rows = [
        dict(r, source=str(source.relative_to(HERE)))
        for source in sources
        for r in csv.DictReader(source.open())
    ]
    for column, coordinate in enumerate((7, 5)):
        ax = axes[case - 1, column]
        for extra, color, label in (
            (False, "#0072B2", "Other retained / preceding starts"),
            (True, "#D55E00", "Prior corrected-chart start"),
        ):
            series = [
                r
                for r in rows
                if int(r["coordinate"]) == coordinate
                and r["chart_criterion_met"] == "True"
                and (r["seed"] == "retained_modified_chart") == extra
            ]
            if not series:
                continue
            ax.scatter(
                [float(r["ln_x_target"]) for r in series],
                [float(r["carbon_residual"]) for r in series],
                s=9,
                color=color,
                marker="x" if extra else "o",
                label=label,
            )
            plotted.extend(
                (
                    case,
                    coordinate,
                    r["attempt"],
                    r["ln_x_target"],
                    r["carbon_residual"],
                    r["seed"],
                    r["source"],
                )
                for r in series
            )
        available = [
            r
            for r in rows
            if int(r["coordinate"]) == coordinate and r["chart_criterion_met"] == "True"
        ]
        if not available:
            best = {}
            for r in rows:
                if int(r["coordinate"]) != coordinate or not r["chart_row_max"]:
                    continue
                t = float(r["ln_x_target"])
                if t not in best or float(r["chart_row_max"]) < float(
                    best[t]["chart_row_max"]
                ):
                    best[t] = r
            ax.scatter(
                list(best),
                [float(r["chart_row_max"]) for r in best.values()],
                s=9,
                color="0.35",
                marker="x",
            )
            ax.axhline(1e-10, color="black", linewidth=0.7)
            ax.set_yscale("log")
            ax.text(
                0.04,
                0.06,
                "No corrected chart points; showing achieved 19-row maxima",
                transform=ax.transAxes,
                fontsize=7,
            )
            plotted.extend(
                (
                    case,
                    coordinate,
                    r["attempt"],
                    r["ln_x_target"],
                    r["chart_row_max"],
                    "achieved_19_row_max",
                    r["source"],
                )
                for r in best.values()
            )
        failed = [
            r
            for r in rows
            if int(r["coordinate"]) == coordinate and r["chart_criterion_met"] != "True"
        ]
        ax.plot(
            [float(r["ln_x_target"]) for r in failed],
            [0.025] * len(failed),
            "|",
            color="0.65",
            markersize=4,
            transform=ax.get_xaxis_transform(),
            label="Failed corrections (coordinate only)",
        )
        if available:
            ax.axhline(0, color="black", linewidth=0.7)
            ax.set_yscale("symlog", linthresh=1e-5)
            ax.set_ylim(min(float(r["carbon_residual"]) for r in available) * 1.3, 2e-5)
        ax.set_xlabel(
            "ln(liquid hydronium mole fraction)"
            if coordinate == 7
            else "ln(liquid bicarbonate mole fraction)"
        )
        ax.set_title(
            f"Case {case}: " + ("c_dielectric = 0" if case == 1 else "epsilon_ion = 2")
        )
        ax.grid(alpha=0.2)
        ax.set_ylabel(
            "Carbon balance row 16\n(dimensionless)"
            if available
            else "Best achieved 19-row maximum\n(dimensionless)"
        )
axes[0, 0].legend(loc="best", fontsize=7)
axes[1, 0].legend(loc="upper left", fontsize=7)
fig.suptitle("Fixed-target chart samples; failed corrections leave gaps", fontsize=13)
fig.savefig(HERE / "figures/global-scan.png", dpi=170)
with (HERE / "results/global-scan/plotted-values.csv").open("w") as file:
    writer = csv.writer(file, lineterminator="\n")
    writer.writerow(
        ("case", "coordinate", "attempt", "ln_x_target", "value", "series", "source")
    )
    writer.writerows(plotted)

"""Render retained continuation values; no thermodynamic calculation runs here."""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--results", type=Path, default=HERE / "results/native-parameter")
parser.add_argument("--output-dir", type=Path)
args = parser.parse_args()
out = args.output_dir or args.results
out.mkdir(parents=True, exist_ok=True)
fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
with (out / "plotted-values.csv").open("w") as file:
    writer = csv.writer(file, lineterminator="\n")
    writer.writerow(
        ("case", "method", "hydronium_mole_fraction", "parameter", "row_max")
    )
    for case, ax in enumerate(axes, 1):
        native = []
        for path in sorted(args.results.glob(f"*/case{case}-states.csv")):
            residuals = {}
            for row in csv.DictReader(path.with_name(f"case{case}-rows.csv").open()):
                if row["accepted"] == "True":
                    value = float(row["value"])
                    residuals[value] = max(residuals.get(value, 0), abs(float(row["residual"])))
            for row in csv.DictReader(path.open()):
                if row["phase"] == "0" and row["species"] == "hydronium-cation":
                    value = float(row["value"])
                    native.append((float(row["mole_fraction"]), value, residuals[value]))
        native = sorted(set(native))
        assert native and max(r[2] for r in native) <= 1e-10, (case, "missing or unclosed native path")
        for H, value, r in native:
            writer.writerow((case, "native parameter continuation", H, value, r))
        ax.plot(
            [r[0] for r in native],
            [r[1] for r in native],
            "-o",
            color="#0072B2",
            markersize=3,
            label="Native parameter continuation",
        )
        chart_series = []
        for folder, style, label in (
            ("turning", "-", "Historical hydronium chart"),
            ("turning-refined", "--", "Historical refined chart"),
        ):
            path = HERE / "results" / folder / f"case{case}-path.csv"
            if not path.exists():
                continue
            rows = [
                r
                for r in csv.DictReader(path.open())
                if float(r["chart_row_max"]) <= 1e-10
            ]
            chart_series.append((rows, style))
            for r in rows:
                writer.writerow(
                    (case, folder, r["hydronium_target"], r["parameter"], r["row_max"])
                )
            ax.plot(
                [float(r["hydronium_target"]) for r in rows],
                [float(r["parameter"]) for r in rows],
                style,
                color="#D55E00",
                label=label,
            )
        if case == 2:
            inset = ax.inset_axes([0.12, 0.38, 0.5, 0.3])
            for rows, style in chart_series:
                inset.plot([float(r["hydronium_target"]) for r in rows],
                           [float(r["parameter"]) for r in rows], style, color="#D55E00")
            inset.axhline(2, color="black", linestyle=":")
            inset.set_xscale("log")
            inset.set_xlim(0.005, 0.06)
            inset.set_ylim(1.99, 2.065)
            inset.set_xticks([0.01, 0.03], labels=["0.01", "0.03"])
            inset.tick_params(axis="x", which="minor", labelbottom=False)
            inset.set_title("Two lower turns", fontsize=8)
            inset.tick_params(labelsize=7)
            inset.grid(alpha=0.2)
        ax.axhline(
            0 if case == 1 else 2,
            color="black",
            linestyle=":",
            label="Requested target",
        )
        ax.set_xscale("log")
        ax.set_xlabel("Liquid hydronium mole fraction (dimensionless)")
        ax.set_ylabel(
            "c_dielectric (dimensionless)"
            if case == 1
            else "Ionic-region relative permittivity (dimensionless)"
        )
        ax.set_title(f"Case {case}: sampled solution branches")
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
fig.savefig(out / "continuation.png", dpi=180)

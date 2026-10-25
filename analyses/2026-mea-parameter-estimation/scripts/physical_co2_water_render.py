from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1] / "results/analysis/physical-co2-water"


def main() -> None:
    plt.rcParams["svg.hashsalt"] = "issue-115-co2-water"
    with (HERE / "predictions.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    temperatures = sorted({float(row["temperature_k"]) for row in rows})
    colors = dict(
        zip(
            temperatures,
            plt.get_cmap("viridis").resampled(len(temperatures)).colors,
            strict=True,
        )
    )

    figure, axes = plt.subplots(1, 3, figsize=(11.2, 3.6), constrained_layout=True)
    parity, composition, pressure = axes
    observed = [float(row["observed_pressure_pa"]) / 1.0e6 for row in rows]
    baseline = [float(row["baseline_pressure_pa"]) / 1.0e6 for row in rows]
    constant = [float(row["constant_k_pressure_pa"]) / 1.0e6 for row in rows]
    candidate = [float(row["temperature_k_pressure_pa"]) / 1.0e6 for row in rows]
    limits = (
        min(observed + baseline + constant + candidate),
        max(observed + baseline + constant + candidate),
    )
    parity.plot(limits, limits, color="0.3", linewidth=1.0, label="parity")
    parity.scatter(observed, baseline, marker="x", s=24, color="0.6", label="$k=0$")
    parity.scatter(
        observed,
        constant,
        marker="o",
        s=24,
        facecolors="none",
        edgecolors="#0072B2",
        label="constant $k$",
    )
    parity.scatter(
        observed, candidate, marker=".", s=22, color="#D55E00", label="T candidate"
    )
    parity.set(
        xscale="log",
        yscale="log",
        xlabel="Observed pressure / MPa",
        ylabel="Predicted pressure / MPa",
        title="a  Pressure parity",
    )
    parity.legend(frameon=False, fontsize=7)

    for temperature_k in temperatures:
        group = [row for row in rows if float(row["temperature_k"]) == temperature_k]
        x = [float(row["liquid_co2_mole_fraction"]) for row in group]
        residual = [
            math.log(
                float(row["constant_k_pressure_pa"])
                / float(row["observed_pressure_pa"])
            )
            for row in group
        ]
        composition.plot(
            x,
            residual,
            marker="o",
            linewidth=1.0,
            markersize=3.5,
            color=colors[temperature_k],
            label=f"{temperature_k:.2f} K",
        )
        pressure.scatter(
            [float(row["observed_pressure_pa"]) / 1.0e6 for row in group],
            residual,
            marker="o",
            s=18,
            color=colors[temperature_k],
        )
    for axis in (composition, pressure):
        axis.axhline(0.0, color="0.3", linewidth=0.8)
        axis.set_ylabel("log($P_{calc}/P_{obs}$)")
    composition.set(
        xscale="log",
        xlabel="Liquid CO$_2$ mole fraction",
        title="b  Constant-k residual by composition",
    )
    composition.legend(frameon=False, fontsize=7)
    pressure.set(
        xscale="log",
        xlabel="Observed pressure / MPa",
        title="c  Constant-k residual by pressure",
    )

    figure.savefig(
        HERE / "pressure-evidence.png",
        dpi=220,
        metadata={"Software": "Matplotlib"},
    )
    svg = HERE / "pressure-evidence.svg"
    figure.savefig(
        svg,
        dpi=220,
        metadata={"Date": None},
    )
    svg.write_text(
        "\n".join(
            line.rstrip() for line in svg.read_text(encoding="utf-8").splitlines()
        )
        + "\n",
        encoding="utf-8",
    )
    plt.close(figure)


if __name__ == "__main__":
    main()

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt

OUTPUT = Path(__file__).resolve().parents[1] / "results/analysis/figures/output"
plt.rcParams["svg.hashsalt"] = "mea-parameter-estimation"


def rows(name: str) -> list[dict[str, str]]:
    with (OUTPUT / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def save(figure: plt.Figure, stem: str) -> None:
    for suffix in ("svg", "png", "pdf"):
        path = OUTPUT / f"{stem}.{suffix}"
        metadata = (
            {"Date": None}
            if suffix == "svg"
            else ({"CreationDate": None, "ModDate": None} if suffix == "pdf" else {})
        )
        figure.savefig(
            path,
            bbox_inches="tight",
            dpi=180,
            metadata=metadata,
        )
        if suffix == "svg":
            path.write_text(
                "\n".join(line.rstrip() for line in path.read_text().splitlines())
                + "\n"
            )
    plt.close(figure)


def neutral_selection() -> None:
    data = rows("neutral-association-selection.csv")
    figure, axis = plt.subplots(figsize=(6.2, 3.8))
    for index, item in enumerate(data):
        selected = item["selected"] == "True"
        axis.errorbar(
            index,
            float(item["grouped_mean_score"]),
            yerr=float(item["grouped_standard_error"]),
            fmt="o" if selected else "x",
            color="#0072B2" if selected else "#D55E00",
            markersize=8,
            capsize=4,
            label="selected" if selected else ("not selected" if index == 1 else None),
        )
    axis.set_xticks(range(len(data)), [item["scheme"] for item in data])
    axis.set_yscale("log")
    axis.set_xlabel("MEA association family")
    axis.set_ylabel("Grouped diagnostic score (dimensionless)")
    axis.set_title("2B is the only low-score, phase-successful family")
    axis.grid(axis="y", alpha=0.25)
    axis.legend(frameon=False)
    save(figure, "neutral-association-selection")


def kiepe_parity() -> None:
    data = rows("kiepe-pressure-parity.csv")
    figure, axis = plt.subplots(figsize=(5.2, 4.4))
    temperatures = sorted({float(item["temperature_k"]) for item in data})
    markers = ("o", "s", "^", "D")
    colors = ("#0072B2", "#D55E00", "#009E73", "#CC79A7")
    for marker, color, temperature in zip(markers, colors, temperatures, strict=True):
        selected = [
            item for item in data if float(item["temperature_k"]) == temperature
        ]
        axis.scatter(
            [float(item["observed_pressure_pa"]) for item in selected],
            [float(item["predicted_pressure_pa"]) for item in selected],
            marker=marker,
            facecolors="none",
            edgecolors=color,
            label=f"{temperature:.2f} K",
        )
    values = [
        float(item[key])
        for item in data
        for key in ("observed_pressure_pa", "predicted_pressure_pa")
    ]
    lower, upper = min(values), max(values)
    axis.plot(
        [lower, upper], [lower, upper], "--", color="black", linewidth=1, label="parity"
    )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("Observed pressure (Pa)")
    axis.set_ylabel("Predicted pressure (Pa)")
    axis.set_title("Historical CO2-water diagnostic (pre-correction m=2.079)")
    axis.grid(which="both", alpha=0.2)
    axis.legend(frameon=False, fontsize=8)
    save(figure, "kiepe-pressure-parity")


if __name__ == "__main__":
    neutral_selection()
    kiepe_parity()

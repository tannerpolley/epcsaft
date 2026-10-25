from __future__ import annotations

import csv
from collections import defaultdict
from itertools import pairwise
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

ANALYSIS = Path(__file__).resolve().parents[1]
RESULTS = ANALYSIS / "results"
FIGURES = ANALYSIS / "figures"

COLORS = {
    "blue": "#1769aa",
    "teal": "#00897b",
    "amber": "#d97706",
    "red": "#c62828",
    "gray": "#5f6368",
}


def rows(name: str) -> list[dict[str, str]]:
    with (RESULTS / name).open(newline="") as stream:
        return list(csv.DictReader(stream))


def save(figure: plt.Figure, name: str) -> None:
    figure.savefig(
        FIGURES / f"{name}.png", dpi=180, bbox_inches="tight", metadata={"Date": None}
    )
    figure.savefig(
        FIGURES / f"{name}.svg", bbox_inches="tight", metadata={"Date": None}
    )
    plt.close(figure)


def style() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "svg.hashsalt": "issue-10-neutral-adapter",
        }
    )


def domain_envelope() -> None:
    grouped: dict[float, list[dict[str, str]]] = defaultdict(list)
    for row in rows("domain.csv"):
        if row["source"] == "conditioning-grid":
            grouped[float(row["temperature_k"])].append(row)

    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.1))
    palette = [
        COLORS["blue"], COLORS["teal"], COLORS["amber"],
        COLORS["red"], COLORS["gray"],
    ]
    for color, temperature in zip(palette, sorted(grouped), strict=True):
        data = sorted(
            grouped[temperature], key=lambda item: float(item["density_mol_per_m3"])
        )
        density = [float(item["density_mol_per_m3"]) for item in data]
        pressure = [float(item["pressure_pa"]) / 1000 for item in data]
        slope = [float(item["pressure_slope_over_rt"]) for item in data]
        axes[0].plot(
            density, pressure, marker="o", ms=3.5, color=color,
            label=f"{temperature:g} K",
        )
        axes[1].plot(density, slope, marker="o", ms=3.5, color=color)

    axes[0].set(
        title="Sampled pressure branches",
        xlabel="Molar density (mol m$^{-3}$)",
        ylabel="Pressure (kPa)",
    )
    axes[0].legend(ncol=2)
    axes[1].axhline(
        0.8, color=COLORS["red"], linestyle="--", label="admission limit"
    )
    axes[1].set(
        title="Density conditioning",
        xlabel="Molar density (mol m$^{-3}$)",
        ylabel=r"$(\partial P/\partial \rho)_T/(RT)$",
        ylim=(0.78, 1.01),
    )
    axes[1].legend()
    figure.suptitle(
        "Sampled neutral propane pressure branch: 45 retained grid states",
        fontweight="bold",
    )
    figure.tight_layout()
    save(figure, "domain-envelope")


def derivative_performance() -> None:
    derivative_rows = [
        row
        for row in rows("derivatives.csv")
        if row["property"] == "pressure" and row["step"]
    ]
    by_step: dict[float, dict[str, float]] = defaultdict(
        lambda: {"gradient": 0.0, "hessian": 0.0}
    )
    for row in derivative_rows:
        step = float(row["step"])
        by_step[step]["gradient"] = max(
            by_step[step]["gradient"],
            float(row["gradient_error"]) / float(row["gradient_target"]),
        )
        by_step[step]["hessian"] = max(
            by_step[step]["hessian"],
            float(row["hessian_error"]) / float(row["hessian_target"]),
        )

    samples: dict[tuple[str, int], list[float]] = defaultdict(list)
    for row in rows("cost.csv"):
        if row["kind"] == "warm_batch":
            samples[(row["path"], int(row["actual_derivative_order"]))].append(
                1000 * float(row["elapsed_s"])
            )
    medians = {
        key: sorted(values)[len(values) // 2] for key, values in samples.items()
    }

    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.1))
    steps = sorted(by_step, reverse=True)
    axes[0].loglog(
        steps, [by_step[s]["gradient"] for s in steps], "o-",
        color=COLORS["blue"], label="gradient",
    )
    axes[0].loglog(
        steps, [by_step[s]["hessian"] for s in steps], "s-",
        color=COLORS["amber"], label="Hessian",
    )
    axes[0].axhline(
        1.0, color=COLORS["red"], linestyle="--", label="acceptance ratio = 1"
    )
    axes[0].invert_xaxis()
    axes[0].set(
        title="Worst pressure finite-difference ratio",
        xlabel="Dimensionless centered-difference step",
        ylabel="error / acceptance target",
    )
    axes[0].legend()

    orders = [0, 1, 2]
    width = 0.34
    axes[1].bar(
        [item - width / 2 for item in orders],
        [medians[("native", item)] for item in orders],
        width, color=COLORS["teal"], label="native",
    )
    axes[1].bar(
        [item + width / 2 for item in orders],
        [medians[("asl", item)] for item in orders],
        width, color=COLORS["blue"], label="ASL",
    )
    axes[1].set_xticks(orders, ["value", "gradient", "Hessian"])
    axes[1].set(title="Warm evaluation cost", ylabel="ms per 100 calls")
    axes[1].legend()
    figure.suptitle(
        "Exact derivatives reach Ipopt with bounded adapter overhead",
        fontweight="bold",
    )
    figure.tight_layout()
    save(figure, "derivative-performance")


def process_results() -> None:
    process = rows("process.csv")
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.1))
    colors = [
        COLORS["blue"] if float(row["cap"]) == 200000 else COLORS["amber"]
        for row in process
    ]
    positions = list(range(len(process)))
    labels = [
        f"{float(row['cap']) / 1000:g} kPa\nstart {float(row['start_rho']):g}"
        for row in process
    ]
    for position, row, color in zip(positions, process, colors, strict=True):
        start = float(row["start_rho"])
        solved = float(row["rho"])
        axes[0].plot([position, position], [start, solved], color=color, alpha=0.65)
        axes[0].scatter(position, start, marker="x", s=60, color=color)
        axes[0].scatter(position, solved, marker="o", s=48, color=color)
    axes[0].set_xticks(positions, labels)
    axes[0].set(
        title="Three starts converge to one optimum",
        ylabel="Molar density (mol m$^{-3}$)",
    )
    axes[0].text(
        0.03, 0.95,
        "x start   o optimum\nblue 200 kPa   amber 150 kPa",
        transform=axes[0].transAxes, va="top",
    )

    stationarity = [float(row["stationarity_inf"]) for row in process]
    complementarity = [float(row["complementarity_inf"]) for row in process]
    width = 0.36
    axes[1].bar(
        [item - width / 2 for item in positions], stationarity, width,
        color=COLORS["teal"], label="stationarity",
    )
    axes[1].bar(
        [item + width / 2 for item in positions], complementarity, width,
        color=COLORS["blue"], label="complementarity",
    )
    axes[1].axhline(
        1e-7, color=COLORS["red"], linestyle="--", label="KKT criterion"
    )
    axes[1].set_yscale("log")
    axes[1].set_xticks(positions, labels)
    axes[1].set(title="Selected KKT residuals", ylabel="Infinity norm")
    axes[1].legend()
    figure.suptitle(
        "Isothermal throughput optimization at 400 K", fontweight="bold"
    )
    figure.tight_layout()
    save(figure, "process-results")


def box(
    axis: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    title: str,
    body: str,
    color: str,
) -> None:
    axis.add_patch(
        Rectangle(xy, width, height, facecolor="white", edgecolor=color, linewidth=2)
    )
    axis.text(
        xy[0] + width / 2, xy[1] + height * 0.68, title,
        ha="center", va="center", fontsize=11, fontweight="bold", color=color,
    )
    axis.text(
        xy[0] + width / 2, xy[1] + height * 0.32, body,
        ha="center", va="center", fontsize=9,
    )


def arrow(
    axis: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    label: str = "",
) -> None:
    axis.add_patch(
        FancyArrowPatch(
            start, end, arrowstyle="-|>", mutation_scale=14,
            linewidth=1.7, color=COLORS["gray"],
        )
    )
    if label:
        axis.text(
            (start[0] + end[0]) / 2,
            (start[1] + end[1]) / 2 + 0.045,
            label, ha="center", fontsize=9,
        )


def flowsheet() -> None:
    figure, axis = plt.subplots(figsize=(11, 4.3))
    axis.set(xlim=(0, 1), ylim=(0, 1))
    axis.axis("off")
    box(
        axis, (0.05, 0.36), 0.19, 0.28,
        "Inlet port", "$F_{in}$, 400 K, $P_{in}$", COLORS["blue"],
    )
    box(
        axis, (0.38, 0.27), 0.24, 0.46,
        "Neutral throughput unit",
        "$F_{in}=F_{out}$\n$P_{in}=P_{out}$\n$F_{out}=\\rho Q$",
        COLORS["teal"],
    )
    box(
        axis, (0.76, 0.36), 0.19, 0.28,
        "Outlet port", "$F_{out}$, 400 K, $P_{out}$", COLORS["blue"],
    )
    arrow(axis, (0.24, 0.5), (0.38, 0.5), "feed")
    arrow(axis, (0.62, 0.5), (0.76, 0.5), "product")
    box(
        axis, (0.07, 0.02), 0.25, 0.14,
        "Optimization",
        "$P_{out}\\leq P_{cap}$\nminimize $0.5q^2-f_{out}$",
        COLORS["red"],
    )
    box(
        axis, (0.38, 0.02), 0.24, 0.14,
        "EOS closure", "$P_{out}=P_{ePC-SAFT}(400\\,K,\\rho)$", COLORS["amber"],
    )
    arrow(axis, (0.32, 0.09), (0.42, 0.27))
    arrow(axis, (0.50, 0.16), (0.50, 0.27))
    axis.text(
        0.5, 0.91, "Issue #10 process flowsheet",
        ha="center", fontsize=14, fontweight="bold",
    )
    axis.text(
        0.5, 0.84,
        "Steady, isothermal, zero pressure drop; no holdup, reaction, phase "
        "equilibrium, heat duty, or compressor work",
        ha="center", fontsize=9.5, color=COLORS["gray"],
    )
    save(figure, "flowsheet")


def adapter_path() -> None:
    figure, axis = plt.subplots(figsize=(12, 5.2))
    axis.set(xlim=(0, 1), ylim=(0, 1))
    axis.axis("off")
    items = [
        ((0.03, 0.58), 0.17, "IDAES / Pyomo", "unit variables\nand constraints", COLORS["blue"]),
        ((0.23, 0.58), 0.14, "NL file", "symbolic rows\nand scaling", COLORS["gray"]),
        ((0.40, 0.58), 0.15, "Ipopt + ASL", "requests F / G / H", COLORS["amber"]),
        ((0.58, 0.58), 0.17, "libepcsaft_asl", "validates arguments\nand maps arrays", COLORS["teal"]),
        ((0.78, 0.58), 0.19, "Engine Properties", "one Helmholtz tape\nexact contractions", COLORS["blue"]),
    ]
    for (x, y), width, title, body, color in items:
        box(axis, (x, y), width, 0.25, title, body, color)
    for left, right in pairwise(items):
        arrow(axis, (left[0][0] + left[1], 0.705), (right[0][0], 0.705))
    arrow(
        axis, (0.875, 0.53), (0.475, 0.29),
        "value + exact gradient + packed Hessian",
    )
    arrow(
        axis, (0.475, 0.25), (0.115, 0.53),
        "solver residuals and KKT steps",
    )
    axis.text(
        0.5, 0.94, "How the compiled adapter transports exact curvature",
        ha="center", fontsize=14, fontweight="bold",
    )
    axis.text(
        0.5, 0.08,
        "The adapter maps arguments and errors only. Thermodynamic equations "
        "and differentiation remain in the Engine.",
        ha="center", fontsize=10, color=COLORS["gray"],
    )
    save(figure, "adapter-path")


def main() -> None:
    FIGURES.mkdir(exist_ok=True)
    style()
    domain_envelope()
    derivative_performance()
    process_results()
    flowsheet()
    adapter_path()


if __name__ == "__main__":
    main()

"""Summarize and plot retained values without importing the Engine."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

CASE = Path(__file__).resolve().parents[1]
DATA = CASE / "results/data"
RENDERED = CASE / "results/rendered"


def main():
    RENDERED.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(DATA / "state-comparisons.csv")
    default = pd.read_csv(DATA / "default-start-comparisons.csv")
    grouped = data.groupby(["figure", "dataset", "names", "mode", "quantity"], sort=False)
    summaries = []
    for key, group in grouped:
        good = group[group.success]
        default_group = default.loc[group.index]
        is_curve = key[1] in ("paper_model", "paper_pc_saft_curve")
        summaries.append(dict(zip(("figure", "dataset", "components", "mode", "quantity"), key),
                              declared=len(group), default_returned=int(default_group.success.sum()),
                              default_failed=int((~default_group.success).sum()), returned=len(good), failed=len(group)-len(good),
                              within_curve_tolerance=int(good.within_curve_tolerance.eq(True).sum()) if is_curve else "",
                              mean_absolute_deviation=good.deviation.abs().mean(), maximum_absolute_deviation=good.deviation.abs().max(),
                              aad_percent=good.absolute_relative_percent.mean(), maximum_residual=good.residual.max()))
    assert sum(s["default_returned"] for s in summaries) == int(default.success.sum()), "summary disagrees with default-start comparisons"
    assert sum(s["returned"] for s in summaries) == int(data.success.sum()), "summary disagrees with state comparisons"
    pd.DataFrame(summaries).to_csv(DATA / "comparison-summary.csv", index=False)
    source = list(csv.reader((CASE / "results/tables/table-001/source/table-001.csv").open()))[1:]
    params = pd.read_csv(CASE / "inputs/parameters/pure/any-solvent.csv").set_index("component")
    pure = []
    for r in source:
        if not r[2]:
            continue
        predicted = data[(data.figure == "table-001") & (data.names == r[0])].iloc[0]
        liquid = data[(data.names == r[0]) & (data.dataset == "experimental_marker") & (data.series == "liquid")]
        good = liquid[liquid.success]
        densities = json.loads(predicted.densities_mol_m3)
        mass = params.loc["H2O" if r[0] == "water" else r[0], "MW"]
        pure.append(dict(component=r[0], association="2B", published_pressure_aad_percent=float(r[7]),
                         published_liquid_aad_percent=float(r[8]), source_temperature_range_K=r[9],
                         marker_liquid_declared=len(liquid), marker_liquid_returned=len(good),
                         marker_density_aad_percent=good.liquid_density_aad_percent.mean(),
                         marker_volume_aad_percent=good.liquid_volume_aad_percent.mean(),
                         prediction_temperature_K=predicted.T_solved_K, prediction_pressure_Pa=predicted.P_solved_Pa,
                         prediction_liquid_density_kg_m3=densities[0]*mass,
                         original_fit_aad_recomputed=False))
    pd.DataFrame(pure).to_csv(DATA / "pure-component-summary.csv", index=False)
    plotted = data[data.figure != "table-001"].copy()
    plotted.to_csv(DATA / "plotted-states.csv", index=False)
    fig, axes = plt.subplots(5, 2, figsize=(12, 19), constrained_layout=True)
    titles = ["Pure alcohol densities", "Methanol / isobutane", "1-propanol / ethylbenzene", "1-pentanol / benzene",
              "Propanol / benzene", "1-butanol / butane", "Ethanol / butane", "Methanol / cyclohexane",
              "Methanol / 1-octanol", "Water / 1-pentanol"]
    for number, (ax, title) in enumerate(zip(axes.flat, titles), 1):
        group = plotted[plotted.figure == f"figure-{number:02}"]
        for role, color, marker in (("paper", "#333333", "o"), ("experiment", "#D55E00", "+"), ("supplemental", "#CC79A7", "."), ("correlation", "#009E73", "s"), ("Engine", "#0072B2", "x")):
            selected = (group[group.dataset.isin(("paper_model", "paper_pc_saft_curve"))] if role == "paper"
                        else group[group.dataset.isin(("experiment", "experimental_marker", "experiment_paper"))] if role == "experiment"
                        else group[group.dataset == "experiment_supplemental"] if role == "supplemental"
                        else group[group.dataset == "paper_correlation"] if role == "correlation" else group[group.success])
            if selected.empty:
                continue
            if number == 1:
                x, y = selected["T"], selected.value if role == "Engine" else selected.reference
            else:
                x = selected.x.copy()
                y = selected.value.copy() if role == "Engine" else selected.reference.copy()
                lle = selected.quantity == "x"
                x.loc[lle] = selected.loc[lle, "value"] if role == "Engine" else selected.loc[lle, "reference"]
                y.loc[lle] = selected.loc[lle, "T"]
            ax.scatter(x, y, s=17, color=color, marker=marker, label=role, linewidths=.8,
                       facecolors="none" if marker == "o" else color)
        original = group[group.source.str.startswith("Gross and Sadowski 2002 Figure") & group.dataset.isin(("experiment_paper", "paper_correlation"))]
        if number == 8:
            ax.errorbar(original.reference, original["T"], xerr=original.tolerance, yerr=original.input_uncertainty,
                        fmt="none", ecolor="#D55E00", alpha=.35, linewidth=.6)
        if number == 10:
            ax.errorbar(original.x, original.reference, xerr=original.input_uncertainty, yerr=original.tolerance,
                        fmt="none", ecolor="#D55E00", alpha=.35, linewidth=.6)
        ax.set_title(f"Figure {number}: {title}")
        ax.set_xlabel("Temperature / K" if number == 1 else "Mole fraction of first declared component")
        ax.set_ylabel("Mass density / kg m⁻³" if number == 1 else "Pressure / bar" if number in (2, 4, 5, 6, 7) else "Temperature / K")
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.savefig(RENDERED / "model-and-experiment.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()

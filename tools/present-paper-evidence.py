#!/usr/bin/env python3
"""Format the five issue-831 notebooks from retained values; never load the Engine."""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "analyses"
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]


def table(leaf, name, frame, caption):
    out = leaf / "figures"
    out.mkdir(exist_ok=True)
    (out / f"notebook-{name}.md").write_text(
        caption + "\n\n" + frame.astype(object).where(frame.notna(), None).to_markdown(
            index=False, floatfmt=".6g", missingval="not applicable") + "\n"
    )


def figure(leaf, fig):
    fig.savefig(leaf / "figures/notebook-comparison.png", dpi=180)
    plt.close(fig)


def polar():
    leaf = ROOT / "2026-polar-dipole"
    rows = pd.read_csv(leaf / "results/pure-rows.csv")
    summary = pd.read_csv(leaf / "results/pure-summary.csv")
    table(leaf, "pure", summary[["component_id", "model", "AAD_psat_all", "AAD_rho_all",
                                 "AAD_psat_validation", "AAD_rho_validation", "failures"]],
          "Pure-fluid average absolute relative deviations (%), from `results/pure-summary.csv`. "
          "The validation columns mean withheld rows only for local refits; published models retain unknown original fit overlap.")
    binary = pd.read_csv(leaf / "results/binary-xy-summary.csv")
    table(leaf, "binary", binary,
          "Paired Campbell comparisons from `results/binary-xy-summary.csv`. Pressure deviations are %, "
          "composition deviations are mole fraction; validation means withheld from this study's binary interaction fit.")
    domain = rows[rows.model == "polar/published"].groupby("component_id").agg(
        rows=("T_K", "size"), minimum_T_K=("T_K", "min"), maximum_T_K=("T_K", "max"))
    table(leaf, "domain", domain.reset_index(), "Source-state counts and temperature domain from `results/pure-rows.csv`.")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    for model, color, marker in zip(["polar/published", "nonpolar/refit"], COLORS, ["s", "^"]):
        group = rows[(rows.component_id == "acetone") & (rows.model == model)]
        axes[0].plot(group.T_K, group.rho_model_mol_m3, marker, color=color, label=model)
    ref = rows[(rows.component_id == "acetone") & (rows.model == "polar/published")]
    axes[0].errorbar(ref.T_K, ref.rho_ref_mol_m3,
                     yerr=ref.rho_ref_mol_m3 * ref.u_rho_liq_rel, fmt="o", mfc="none", color="black",
                     label="reference EOS; stated bounds")
    xy = pd.read_csv(leaf / "results/binary-xy-rows.csv")
    held = xy[(xy.role == "validation") & (xy.variant == "k_ij calibrated")]
    for model, color, marker in zip(["polar/GV2006", "nonpolar/refit"], COLORS, ["s", "^"]):
        group = held[held.model == model]
        assert not group.empty, model
        axes[1].plot(group.P_data_kPa, group.P_model_kPa, marker, color=color, label=model)
    low, high = held.P_data_kPa.min(), held.P_data_kPa.max()
    axes[1].plot([low, high], [low, high], color="black", lw=1, label="equality")
    axes[0].set(xlabel="Temperature / K", ylabel="Acetone liquid density / mol m⁻³")
    axes[1].set(xlabel="Campbell measured pressure / kPa", ylabel="Held-out bubble pressure / kPa")
    for ax in axes:
        ax.legend(fontsize=8)
        ax.grid(alpha=.2)
    figure(leaf, fig)


def propane():
    leaf = ROOT / "2026-pure-saturation-regression"
    rows = pd.read_csv(leaf / "results/propane-replay-2026-10-01.csv")
    training = pd.read_csv(leaf / "inputs/source/propane-baseline.csv")
    replay = json.loads((leaf / "results/propane-replay-2026-10-01.json").read_text())
    # The modern fit owns training membership; retain the packet's stress declarations.
    rows["replay_role"] = "held_out"
    rows.loc[rows.role == "stress", "replay_role"] = "stress"
    rows.loc[rows.T_K.isin(training.temperature_k), "replay_role"] = "training"
    assert set(rows.loc[rows.replay_role == "training", "T_K"]) == set(training.temperature_k)
    assert (rows.replay_role == "held_out").sum() == (rows.role == "held_out").sum() + 1
    table(leaf, "roles", rows[["row_id", "T_K", "role", "replay_role"]],
          "Role correction from the actual training CSV. `role` is the unchanged historical packet label; "
          "`replay_role` governs every panel and metric below.")
    metrics = []
    for role, group in rows.groupby("replay_role", sort=False):
        for model in ("fitted", "packet"):
            assert (group[f"{model}_status"] == "Available").all()
            metrics.append([role, model, len(group), group.T_K.min(), group.T_K.max(),
                            100 * group[f"p_{model}_relative_error"].abs().mean(),
                            100 * group[f"rho_{model}_relative_error"].abs().mean()])
    table(leaf, "metrics", pd.DataFrame(metrics, columns=["Role", "Model", "States", "Minimum T / K",
          "Maximum T / K", "Pressure AARD / %", "Liquid density AARD / %"]),
          "Average absolute relative deviation (AARD), computed from `results/propane-replay-2026-10-01.csv`; "
          "unavailable predictions would stop presentation, rather than disappear from the mean.")
    fit = replay["fitted"]
    table(leaf, "parameters", pd.DataFrame({"Coordinate": fit["identities"], "Unit": fit["units"],
          "Fitted": fit["values"], "Packet": fit["packet_values"],
          "Conditional standard deviation": replay["covariance"]["physical_std_scaled_by_variance_factor"]}),
          "Selected first-start parameters and conditional standard deviations from the dated JSON; "
          "these are local fit diagnostics, not prediction intervals.")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    for role, color, marker in zip(["training", "held_out", "stress"], COLORS, ["s", "o", "^"]):
        group = rows[rows.replay_role == role]
        for ax, quantity in zip(axes, ["p", "rho"]):
            ax.plot(group.T_K, 100 * group[f"{quantity}_fitted_relative_error"], marker,
                    color=color, label=role)
            ax.errorbar(group.T_K, [0] * len(group),
                        yerr=100 * group[f"{quantity}_expanded_uncertainty_" + ("Pa" if quantity == "p" else "kg_m3")]
                        / group[f"{quantity}_source_" + ("Pa" if quantity == "p" else "kg_m3")],
                        fmt="none", ecolor=color, alpha=.4)
    for ax, name in zip(axes, ["Saturation pressure", "Saturated liquid density"]):
        ax.set(xlabel="Temperature / K", ylabel=f"{name} relative deviation / %")
        ax.axhline(0, color="black", lw=.7)
        ax.legend(fontsize=8)
        ax.grid(alpha=.2)
    figure(leaf, fig)


def dielectric():
    leaf = ROOT / "2026-dielectric-alternatives"
    rows = pd.read_csv(leaf / "results/comparison.csv")
    fits = pd.read_csv(leaf / "results/refit.csv")
    summary = json.loads((leaf / "results/summary.json").read_text())
    (leaf / "figures/notebook-conditions.md").write_text(summary["question"] + "\n")
    assert (rows.status == "ok").all()
    values = []
    for variant, group in rows.groupby("variant"):
        activity = group[group.property == "miac"]
        eps = group[group.property == "relative_permittivity"]
        values.append([variant, group.k_water_ion.iloc[0], len(activity[activity.condition_value <= 2]),
                       100 * activity[activity.condition_value <= 2].relative_error.abs().mean(),
                       len(activity[activity.condition_value > 2]),
                       100 * activity[activity.condition_value > 2].relative_error.abs().mean(),
                       100 * eps.relative_error.abs().mean()])
    frame = pd.DataFrame(values, columns=["Variant", "Shared k", "Low-m rows", "Low-m MARD / %",
                                         "High-m rows", "High-m MARD / %", "Permittivity MARD / %"])
    # One decimal is the source of the reported rounding defect.
    frame["High-m MARD / %"] = frame["High-m MARD / %"].map(lambda x: f"{x:.1f}")
    table(leaf, "metrics", frame,
          "Mean absolute relative deviations (MARD), from `results/comparison.csv`. High-m values are rounded "
          "to one decimal before display. The low/high division is the modern-fit calibration cutoff; "
          "F0's published interaction used both groups.")
    table(leaf, "fits", fits[["form", "calibration_rows", "k_water_ion", "sigma_k", "rank", "variance_factor", "at_bound"]],
          "Conditional native-fit diagnostics from `results/refit.csv`; sigma covers the shared interaction only.")
    properties = []
    for variant in frame.Variant:
        metrics = summary["metrics"]
        solvation = next(value for key, value in metrics.items()
                         if key.startswith(f"{variant}/solvation_gibbs_kJ_mol/"))
        permittivity = next(value for key, value in metrics.items()
                            if key.startswith(f"{variant}/relative_permittivity/"))
        high_m = next(value for key, value in metrics.items()
                      if key.startswith(f"{variant}/miac/") and key.endswith("_m_gt_2"))
        density = metrics[f"{variant}/liquid_molar_density_mol_m3/difference_from_same_fit_F"]
        sigma = high_m["mard_sigma_from_k"]
        properties.append([variant, 100 * permittivity["max_ard"], solvation["mean_abs_error"],
                           100 * density["max_ard"], None if sigma is None else 100 * sigma])
    table(leaf, "properties", pd.DataFrame(properties, columns=["Variant", "Maximum permittivity deviation / %",
          "Mean absolute solvation difference / kJ mol⁻¹", "Maximum density difference from same-fit F / %",
          "High-m activity MARD sigma / percentage points"]),
          "Other retained properties and interaction-only propagation from `results/summary.json:metrics`. "
          "Density differences use a model reference, not observations; zero for F is its comparison with itself. "
          "Not applicable means no fitted-interaction propagation for a fixed-parameter variant.")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    for variant, color, marker in zip(["F1", "D1", "P1"], COLORS, ["s", "^", "x"]):
        group = rows[rows.variant == variant]
        for ax, prop in zip(axes, ["miac", "relative_permittivity"]):
            selected = group[group.property == prop]
            ax.plot(selected.condition_value, selected.model, marker, color=color, label=variant)
    for ax, prop in zip(axes, ["miac", "relative_permittivity"]):
        ref = rows[(rows.variant == "F1") & (rows.property == prop)]
        ax.plot(ref.condition_value, ref.observed, "o", mfc="none", color="black", label="reference")
        ax.legend(fontsize=8)
        ax.grid(alpha=.2)
    axes[0].set(xlabel="NaCl molality / mol kg⁻¹ water", ylabel="Mean ionic activity coefficient / 1")
    axes[0].axvline(2, color="black", ls=":", lw=1)
    axes[1].set(xlabel="Ion mole fraction / 1", ylabel="Relative permittivity / 1")
    figure(leaf, fig)


def wangler():
    leaf = ROOT / "reproduction/2018-wangler"
    osmotic = json.loads((leaf / "figures/figure-10/output/wangler-2018-mdeah-osmotic-summary.json").read_text())
    pressure = json.loads((leaf / "figures/figure-14/output/wangler-2018-table-10-native-summary.json").read_text())
    rows = pd.read_csv(leaf / "figures/figure-14/output/wangler-2018-table-10-native.csv")
    osmotic_rows = pd.read_csv(leaf / "figures/figure-10/output/wangler-2018-mdeah-osmotic-plotted-data.csv")
    domain = osmotic_rows[osmotic_rows.lane == "S0"].groupby("temperature_k").agg(
        observations=("molality_mol_per_kg", "size"),
        minimum_molality=("molality_mol_per_kg", "min"), maximum_molality=("molality_mol_per_kg", "max"))
    table(leaf, "domain", domain.reset_index(),
          "Osmotic calibration observations shared by both lanes, from the Figure 10 plotted-data CSV; "
          "temperature is in K, molality in mol/kg water, at the pressure specified in the osmotic summary limits.")
    values = []
    for lane in ("S0", "S1"):
        fit, pred = osmotic["lanes"][lane], pressure["lanes"][f"fitted-{lane}"]
        values.append([lane, fit["ard_percent"], fit["F1"]["ard_limit_percent"], fit["k_ij_303_k"],
                       fit["k_ij_303_k_allowance"], pred["native_ard_percent"],
                       pred["solved_row_count"], 100 * pred["F3"]["maximum_relative_pressure_change"]])
    table(leaf, "metrics", pd.DataFrame(values, columns=["Lane", "Osmotic calibration AARD / %",
          "Calibration limit / %", "Transferred k", "k allowance", "Pressure prediction AARD / %",
          "Solved pressure states", "Maximum pressure change from k allowance / %"]),
          "Calibration, transfer and failed prediction from the Figure 10 and Table 10 native summaries. "
          "The pressure change is the retained allowance probe, not a full uncertainty interval.")
    selected = rows[rows.lane.isin(["fitted-S0", "fitted-S1", "printed-S0"])]
    assert selected.success.all()
    table(leaf, "pressure", selected[["lane", "temperature_k", "co2_loading_mol_per_mol_mdea",
          "experimental_total_pressure_bar", "pressure_standard_uncertainty_bar", "native_pressure_bar",
          "reaction_affinity_max", "matching_k_ij_303_k"]],
          "Table 10 pressure states from the retained native CSV. Row-matching interactions are diagnostic fits "
          "and are never used for the pressure predictions.")
    fig, ax = plt.subplots(figsize=(8, 4.5), layout="constrained")
    measured = rows[rows.lane == "fitted-S0"]
    ax.errorbar(measured.co2_loading_mol_per_mol_mdea, measured.experimental_total_pressure_bar,
                yerr=measured.pressure_standard_uncertainty_bar, fmt="o", mfc="none", color="black", label="withheld Table 10")
    for lane, color, marker in zip(["fitted-S0", "fitted-S1", "printed-S0"], COLORS, ["s", "^", "x"]):
        group = rows[rows.lane == lane]
        ax.plot(group.co2_loading_mol_per_mol_mdea, group.native_pressure_bar, marker, color=color, label=lane)
    ax.set(yscale="log", xlabel="CO₂ loading / mol mol⁻¹ MDEA", ylabel="Total pressure / bar")
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    ax.grid(alpha=.2, which="both")
    figure(leaf, fig)


def cleeton():
    leaf = ROOT / "reproduction/2020-cleeton"
    summary = json.loads((leaf / "results/data/cleeton-2020-summary.json").read_text())
    values = []
    for lane, entry in summary["lanes"].items():
        gate = entry["figure_16_gate"]
        values.append([lane, entry["solved"], gate["within_tolerance"], gate["points"], gate["max_abs_ln_ratio"],
                       entry["dicko_aad_percent"]["co2_figure_16_left_12"],
                       entry["dicko_aad_percent"]["h2s_figure_16_right_18"],
                       entry["dicko_symmetric_deviation_percent"]["co2_figure_16_left_12"],
                       entry["dicko_symmetric_deviation_percent"]["h2s_figure_16_right_18"]])
    table(leaf, "metrics", pd.DataFrame(values, columns=["Lane", "Solved Dicko states", "Curve checks within limit",
          "Curve checks", "Maximum absolute ln ratio", "CO₂ AARD / %", "H₂S AARD / %",
          "CO₂ symmetric deviation / %", "H₂S symmetric deviation / %"]),
          "Caption-reading comparisons and separate physical deviations from `results/data/cleeton-2020-summary.json`. "
          "The gas subsets are those defined by Figure 16, rather than all Dicko states.")
    table(leaf, "ladder", pd.DataFrame(summary["ladder"]).T.reset_index(names="Comparison"),
          "Component and binary comparison ladder from the same summary; figure subsets differ from the source's full fitted datasets.")


if __name__ == "__main__":
    polar()
    propane()
    dielectric()
    wangler()
    cleeton()

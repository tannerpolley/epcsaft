"""Prepare the Wangler notebook from retained CSV/JSON only; never load a model."""

import json
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

LEAF = Path(__file__).resolve().parents[1]
OUT = LEAF / "figures"
COLORS = ["#0072B2", "#D55E00", "#009E73"]


def table(name, frame, caption):
    (OUT / f"notebook-{name}.md").write_text(
        caption
        + "\n\n"
        + frame.to_html(
            index=False,
            float_format=lambda v: f"{v:.8g}",
            na_rep="unavailable",
            border=0,
        )
        + "\n"
    )


def main():
    osmotic = json.loads(
        (OUT / "figure-10/output/wangler-2018-mdeah-osmotic-summary.json").read_text()
    )
    pressure = json.loads(
        (OUT / "figure-14/output/wangler-2018-table-10-native-summary.json").read_text()
    )
    rows = pd.read_csv(OUT / "figure-14/output/wangler-2018-table-10-native.csv")
    osmotic_rows = pd.read_csv(
        OUT / "figure-10/output/wangler-2018-mdeah-osmotic-plotted-data.csv"
    )
    neutral = pd.read_csv(LEAF / "results/neutral-figure-1.csv")
    cross = pd.read_csv(LEAF / "results/cross-code.csv")
    terms = pd.read_csv(LEAF / "results/cross-code-terms.csv")
    pure = pd.read_csv(LEAF / "results/pure-mdea.csv")
    obs = pd.read_csv(LEAF / "results/pure-observations.csv")
    sweep = pd.read_csv(LEAF / "results/convention-sweep.csv")
    baseline = pd.read_csv(LEAF / "results/reference-pressure-baseline.csv")
    diagnostics = json.loads((LEAF / "results/neutral-diagnostics.json").read_text())
    table(
        "domain",
        osmotic_rows[osmotic_rows.lane == "S0"]
        .groupby("temperature_k")
        .agg(
            observations=("molality_mol_per_kg", "size"),
            minimum_molality=("molality_mol_per_kg", "min"),
            maximum_molality=("molality_mol_per_kg", "max"),
        )
        .reset_index(),
        "Calibration observations from Tables 8–9; temperature in K, molality in mol/kg water, at 1 bar.",
    )
    values = []
    for lane in ("S0", "S1"):
        fit, pred = osmotic["lanes"][lane], pressure["lanes"][f"fitted-{lane}"]
        values.append(
            [
                lane,
                fit["ard_percent"],
                fit["k_ij_303_k"],
                fit["k_ij_303_k_allowance"],
                pred["native_ard_percent"],
                pred["solved_row_count"],
                100 * pred["F3"]["maximum_relative_pressure_change"],
            ]
        )
    table(
        "metrics",
        pd.DataFrame(
            values,
            columns=[
                "Convention",
                "Osmotic AARD / %",
                "Transferred k",
                "k allowance",
                "Pressure AARD / %",
                "Solved states",
                "Maximum pressure change from k allowance / %",
            ],
        ),
        "Native osmotic calibration and withheld Table 10 prediction. The allowance is measurement propagation, not total model uncertainty.",
    )
    selected = rows[rows.lane.isin(["fitted-S0", "fitted-S1"])]
    assert len(selected) == 6 and selected.success.all()
    table(
        "pressure",
        selected[
            [
                "lane",
                "co2_loading_mol_per_mol_mdea",
                "experimental_total_pressure_bar",
                "pressure_standard_uncertainty_bar",
                "native_pressure_bar",
                "maximum_original_row",
                "reaction_affinity_max",
            ]
        ],
        "Table 10 at 303 K, apparent solvent mass fraction w(MDEA) = 0.6. Pressures and standard pressure uncertainties are in bar; residuals are dimensionless.",
    )
    failures = rows[~rows.success]
    table(
        "pressure-failures",
        failures[
            [
                "lane",
                "co2_loading_mol_per_mol_mdea",
                "pressure_seed_bar",
                "solver_status",
                "failure",
            ]
        ],
        "Unsolved printed-interaction diagnostic states. No iterate pressure is retained as a prediction.",
    )
    changes = []
    for lane in ("S0", "S1"):
        b = baseline[baseline.lane == lane]
        changes.append(
            [
                lane,
                100 * (b.pressure_bar / b.observed_bar - 1).abs().mean(),
                pressure["lanes"][f"fitted-{lane}"]["native_ard_percent"],
            ]
        )
    table(
        "reference",
        pd.DataFrame(
            changes,
            columns=[
                "Convention",
                "One-bar standard AARD / %",
                "Same-pressure standard AARD / %",
            ],
        ),
        "Isolated reference-pressure correction, Eqs. 20–21. Fitted parameters and species standards are held fixed.",
    )
    table(
        "neutral",
        neutral,
        "Neutral source-curve readouts from Wangler p. 20 Figure 1 at liquid x(MDEA) = 0.2. The ±0.01 bar allowance describes graphical reading, not experimental uncertainty.",
    )
    table(
        "cross",
        cross,
        "Independent bubble solutions at 353 K and liquid x(MDEA) = 0.2; units are stated per quantity.",
    )
    errors = pd.DataFrame(
        {
            "Term": terms.term,
            "Difference in a/(RT)": terms.FeOs_a_over_RT - terms.Engine_a_over_RT,
            "Difference in mu_water/(RT)": terms.FeOs_mu_water_over_RT
            - terms.Engine_mu_water_over_RT,
            "Difference in mu_MDEA/(RT)": terms.FeOs_mu_mdea_over_RT
            - terms.Engine_mu_mdea_over_RT,
        }
    )
    table(
        "terms",
        errors,
        "FeOs minus Engine at the same Engine liquid density. Engine term chemical potentials use derivative verification; FeOs term derivatives and Engine totals are exact native derivatives.",
    )
    table(
        "pure",
        pure[
            [
                "T_K",
                "Engine_P_Pa",
                "FeOs_P_Pa",
                "Engine_rho_kg_m3",
                "Kim_Riedel_P_Pa",
                "Kim_Antoine_P_Pa",
            ]
        ],
        "Pure MDEA saturation using Wangler Table 4. Correlations are from Kim p. 2523 Table 5; pressures in Pa and density in kg/m³.",
    )
    table(
        "pure-observations",
        obs,
        "Kim p. 2522 Table 3 neat-MDEA measurements; pressure in Pa, density in kg/m³. Relative errors use each measured pressure as denominator.",
    )
    table(
        "conventions",
        sweep.drop(columns=["failure"]),
        "One-at-a-time neutral declarations at x(MDEA) = 0.2. All pressure columns and maximum deviations are in bar. Empty values denote an unsolved control; see the retained failure column in results/convention-sweep.csv.",
    )
    (OUT / "notebook-diagnosis-numbers.md").write_text(
        f"The independent bubble pressures differ by {diagnostics['bubble_relative_difference']:.3g} relative. "
        f"The largest term chemical-potential difference is {diagnostics['maximum_term_mu_difference']:.3g} in $\\mu/(RT)$. "
        f"The Engine bubble residual is {diagnostics['Engine_bubble_residual']:.3g}; "
        f"FeOs liquid–vapor pressure closure is {diagnostics['FeOs_pressure_closure_Pa']:.3g} Pa. "
        f"The derivative-step change is at most {diagnostics['term_step_change']:.3g}, and the term sum differs from the exact Engine total by at most {diagnostics['term_sum_error']:.3g}.\n"
    )
    at353 = pure[pure.T_K == 353].iloc[0]
    (OUT / "notebook-pure-numbers.md").write_text(
        f"At 353 K, the Engine gives {at353.Engine_P_Pa:.5g} Pa, bracketed by Kim's Antoine and Riedel values "
        f"of {at353.Kim_Antoine_P_Pa:.5g} and {at353.Kim_Riedel_P_Pa:.5g} Pa. "
        f"Against the seven neat measurements, the AARD is {diagnostics['Kim_neat_measurement_AARD_percent']:.6f}%. "
        f"Pahlavanzadeh Table 1's printed pressure correlation gives {at353.pahlav_Table1_literal_P_Pa:.6g} Pa at that temperature; "
        f"its printed density correlation gives {at353.pahlav_Table1_literal_rho_kg_m3:.6g} kg/m³ when assigned MDEA's molar mass.\n"
    )
    fig, ax = plt.subplots(figsize=(7.4, 4.6), layout="constrained")
    ax.errorbar(
        neutral.T_K,
        neutral.source_curve_bar,
        yerr=neutral.graphical_allowance_bar,
        fmt="o",
        mfc="none",
        color="black",
        capsize=3,
        label="digitized from Wangler et al. 2018, Fig. 1 (±0.01 bar)",
    )
    ax.plot(
        neutral.T_K,
        neutral.Engine_P_bar,
        "s",
        color=COLORS[0],
        label="Engine, printed parameters",
    )
    literal = sweep[sweep["mode"] == "eq12-literal"].iloc[0]
    ax.plot(
        neutral.T_K,
        [literal[f"P_{int(t)}_K_bar"] for t in neutral.T_K],
        "^",
        color=COLORS[1],
        label="Engine, literal Eq. 12 radical",
    )
    ax.set(
        xlabel="Temperature / K",
        ylabel="Total bubble pressure / bar",
        title="Neutral water–MDEA, liquid x(MDEA) = 0.2",
    )
    ax.legend(fontsize=9)
    ax.grid(alpha=0.2)
    fig.savefig(OUT / "notebook-neutral.png", dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.4, 4.6), layout="constrained")
    measured = selected[selected.lane == "fitted-S0"]
    ax.errorbar(
        measured.co2_loading_mol_per_mol_mdea,
        measured.experimental_total_pressure_bar,
        yerr=measured.pressure_standard_uncertainty_bar,
        fmt="o",
        mfc="none",
        color="black",
        capsize=3,
        label="Table 10 observations",
    )
    for lane, color, marker in zip(("S0", "S1"), COLORS, ("s", "^")):
        group = selected[selected.lane == f"fitted-{lane}"]
        ax.plot(
            group.co2_loading_mol_per_mol_mdea,
            group.native_pressure_bar,
            marker,
            color=color,
            label=f"Osmotic fit, {lane}",
        )
    ax.set(
        yscale="log",
        xlabel="CO₂ loading / mol per initial mol MDEA",
        ylabel="Total pressure / bar",
    )
    ax.legend()
    ax.grid(alpha=0.2, which="both")
    fig.savefig(OUT / "notebook-comparison.png", dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7.4, 4.6), layout="constrained")
    ax.errorbar(
        obs.T_K,
        obs.observed_P_Pa,
        yerr=300.0,
        fmt="o",
        mfc="none",
        color="black",
        capsize=3,
        label="Kim Table 3 (±300 Pa)",
    )
    for col, color, marker, label in [
        ("Engine_P_Pa", COLORS[0], "s", "Engine"),
        ("Kim_Riedel_P_Pa", COLORS[1], "^", "Kim Riedel"),
        ("Kim_Antoine_P_Pa", COLORS[2], "x", "Kim Antoine"),
    ]:
        ax.plot(obs.T_K, obs[col], marker, color=color, label=label)
    ax.set(xlabel="Temperature / K", ylabel="Pure-MDEA saturation pressure / Pa")
    ax.legend()
    ax.grid(alpha=0.2)
    fig.savefig(OUT / "notebook-pure.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()

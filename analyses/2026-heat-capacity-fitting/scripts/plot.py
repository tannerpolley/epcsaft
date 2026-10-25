"""Show the retained, discrete reference and fitted predictions without running thermodynamics."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SET = Path(__file__).resolve().parents[1]
with (SET / "results/predictions.csv").open() as stream:
    rows = [{k: v if k in ("case", "quantity") else float(v) for k, v in row.items()}
            for row in csv.DictReader(stream)]
fig, axes = plt.subplots(1, 3, figsize=(11, 3.7), layout="constrained")
for ax, quantity, unit in zip(axes, ("rho_mol_m3", "cp_J_mol_K", "P_sat_Pa"),
    ("Density (mol m⁻³)", "Total cp (J mol⁻¹ K⁻¹)", "Saturation pressure (Pa)")):
    selected = [row for row in rows if row["quantity"] == quantity and row["case"] == "with_cp"]
    T = [row["T_K"] for row in selected]
    ax.plot(T, [row["reference"] for row in selected], "o", color="black", markersize=5, label="IAPWS-95 reference")
    ax.plot(T, [row["before"] for row in selected], "+", color="#009E73", markersize=7, label="Source parameters")
    ax.plot(T, [row["after"] for row in selected], "x", color="#0072B2", markersize=6, label="Fit with cp")
    companion = [row for row in rows if row["quantity"] == quantity and row["case"] == "without_cp"]
    ax.plot(T, [row["after"] for row in companion], "s", color="#D55E00", markersize=4, fillstyle="none", label="Fit without cp")
    ax.set(xlabel="Temperature (K)", ylabel=unit)
    if quantity == "P_sat_Pa":
        ax.set_yscale("log")
    ax.grid(alpha=0.2)
axes[0].legend(fontsize=8)
(SET / "figures").mkdir(exist_ok=True)
fig.savefig(SET / "figures/water-comparison.png", dpi=180)

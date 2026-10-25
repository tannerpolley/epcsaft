"""Plot retained chloride activities without recalculating the EOS."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
curves = list(csv.DictReader((HERE / "results/model-curves.csv").open()))
observed = list(csv.DictReader((HERE / "inputs/evaluated.csv").open()))
fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
for ax, salt in zip(axes.flat, ("NaCl", "LiCl", "MgCl2", "CaCl2")):
    for mode, label, color, dash in (("packet3", "Packet 3: DH = 0.88σ", "#0072B2", "-"),
                                     ("packet4", "Packet 4: DH = σ", "#D55E00", "--")):
        rows = [r for r in curves if r["mode"] == mode and r["salt"] == salt]
        ax.plot([float(r["molality_mol_kg"]) for r in rows], [float(r["gamma_pm_m"]) for r in rows],
                label=label, color=color, linestyle=dash)
    rows = [r for r in observed if r["salt"] == salt]
    label = "Hamer–Wu evaluated" if salt in ("NaCl", "LiCl") else "Issue #209 supplied anchor*"
    ax.scatter([float(r["molality_mol_kg"]) for r in rows], [float(r["gamma_evaluated"]) for r in rows],
               s=18, color="black", marker="o", label=label, zorder=3)
    ax.set(xscale="log", yscale="log", xlabel="Salt molality (mol/kg water)",
           ylabel="Mean ionic activity coefficient γ± (molal)", title=salt)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2)
fig.suptitle("298.15 K, 1 bar: correcting the DH length does not qualify a new divalent fit\n"
             "*Goldberg–Nuttall source unavailable; Mg/Ca anchors are not independently verified", fontsize=11)
fig.savefig(HERE / "figures/chloride-activities.png", dpi=160)

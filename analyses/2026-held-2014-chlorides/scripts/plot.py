"""Plot retained model predictions and evaluated coefficients; no Engine calls."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
curves = list(csv.DictReader((HERE / "results/model-curves.csv").open()))
observed = list(csv.DictReader((HERE / "inputs/evaluated.csv").open()))
fig, axes = plt.subplots(2, 4, figsize=(13, 6), layout="constrained")
for column, salt in enumerate(("LiCl", "NaCl", "MgCl2", "CaCl2")):
    for row, (prediction, reference, label) in enumerate((
        ("gamma_pm_m", "gamma_evaluated", r"$\gamma_\pm$ (molal)"),
        ("osmotic_coefficient", "osmotic_evaluated", r"$\Phi$"),
    )):
        ax = axes[row, column]
        for mode, color, style, title in (("packet4", "#0072B2", "--", "Held 2008 packet 4"),
                                         ("held2014", "#D55E00", "-", r"Held 2014, $\epsilon_r=78.09$")):
            rows = [r for r in curves if r["mode"] == mode and r["salt"] == salt]
            ax.plot([float(r["molality_mol_kg"]) for r in rows],
                    [float(r[prediction]) for r in rows], style, color=color, label=title)
        rows = [r for r in observed if r["salt"] == salt]
        ax.scatter([float(r["molality_mol_kg"]) for r in rows], [float(r[reference]) for r in rows],
                   marker="o", s=13, facecolor="none", edgecolor="black", label="Evaluated source")
        ax.set(xlim=(0, 4.5), ylabel=label)
        if row == 0:
            ax.set_title(salt)
        else:
            ax.set_xlabel("Salt molality (mol/kg water)")
        ax.grid(alpha=.18)
axes[0, 0].legend(fontsize=8)
fig.savefig(HERE / "figures/chloride-coefficients.png", dpi=180)
plt.close(fig)

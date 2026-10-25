"""Plot retained row errors for the issue-398 closest-approach comparison."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
fig, axes = plt.subplots(2, 4, figsize=(13, 6), layout="constrained", sharex=True, sharey=True)
for column, salt in enumerate(("LiCl", "NaCl", "MgCl2", "CaCl2")):
    for row, (quantity, label) in enumerate((("log10_gamma_ratio", r"$\gamma_\pm$"),
                                           ("log10_osmotic_ratio", r"$\Phi$"))):
        ax = axes[row, column]
        for folder, mode, color, marker, title in (
            ("sigma-replay", "held2014", "#0072B2", "o", r"DH $\sigma$"),
            ("dh-088", "held2014_dh088", "#D55E00", "x", r"DH $0.88\sigma$ (inferred)"),
        ):
            rows = [r for r in csv.DictReader((HERE / "results" / folder / "evaluated-comparison.csv").open())
                    if r["mode"] == mode and r["salt"] == salt]
            ax.scatter([float(r["molality_mol_kg"]) for r in rows], [float(r[quantity]) for r in rows],
                       color=color, marker=marker, s=16, label=title)
        ax.axhspan(-0.02, 0.02, color="grey", alpha=0.15, label="Issue #398 target")
        ax.axhline(0, color="grey", linewidth=0.6)
        ax.set(xlim=(0, 4.6), ylabel=label + " error (decade)")
        if row == 0:
            ax.set_title(salt)
        else:
            ax.set_xlabel("Salt molality (mol/kg water)")
        ax.grid(alpha=0.15)
axes[0, 0].legend(fontsize=8)
fig.savefig(HERE / "figures/closest-approach-errors.png", dpi=180)
plt.close(fig)

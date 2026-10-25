"""Figures from retained results only (no model evaluation)."""
from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / "results/figures"
OUT.mkdir(exist_ok=True)
rows = lambda name: list(csv.DictReader((HERE / "results" / name).open()))
SIGMA, UNIFIED, PAPER, EXP = "#1f77b4", "#d95f02", "#7f7f7f", "black"
CONVENTIONS = (("sigma", SIGMA, "Engine, σ"), ("unified-0.88", UNIFIED, "Engine, 0.88σ unified"))

# Fig. 4 / S1: solvent activity coefficient (formula salt axis, ion-alcohol k_ij = 0), both diameter conventions
data = [r for r in rows("fig4-solvent-activity.csv") if r["basis"] == "formula" and r["kij_set"] == "prediction"]
titles = {"fig4-licl": "LiCl in methanol", "fig4-libr": "LiBr in methanol", "figS1-licl-ethanol": "LiCl in ethanol"}
fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), constrained_layout=True)
for ax, (panel, title) in zip(axes, titles.items()):
    for model, style in (("advanced", "-"), ("revised-2014", "--")):
        sel = lambda d: sorted((float(r["x_salt"]), float(r["gamma_species_basis"]), float(r["digitized_gamma"])) for r in data
                               if r["panel"] == panel and r["model"] == model and r["kind"] == "paper-model" and r["diameter"] == d)
        line = sel("sigma")
        ax.plot([p[0] for p in line], [p[2] for p in line], style, color=PAPER, lw=4, alpha=0.45, label=f"paper, {model}")
        for diameter, color, label in CONVENTIONS:
            line = sel(diameter)
            ax.plot([p[0] for p in line], [p[1] for p in line], style, color=color, lw=1.2, label=f"{label}, {model}")
    exp = [(float(r["x_salt"]), float(r["digitized_gamma"])) for r in data
           if r["panel"] == panel and r["model"] == "advanced" and r["kind"] == "experiment" and r["diameter"] == "sigma"]
    ax.plot(*zip(*exp), "o", color=EXP, ms=4, label="experiment (paper)")
    ax.set_title(title)
    ax.set_xlabel(r"$x_\mathrm{salt} = n_\mathrm{salt}/(n_\mathrm{salt}+n_\mathrm{solvent})$")
axes[0].set_ylabel(r"$\gamma_\mathrm{solvent}$ (species basis)")
axes[2].legend(fontsize=6.5)
fig.savefig(OUT / "fig4-solvent-activity.png", dpi=150)

# Fig. 8b: transfer Gibbs energy water -> alcohol, both diameter conventions
data = [r for r in rows("fig8b-transfer-gibbs.csv") if r["water_sigma"] == "correlation"]
fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), constrained_layout=True, sharey=True)
w = 0.2
for ax, kij in zip(axes, ("prediction", "correlation")):
    by = {d: [r for r in data if r["kij_set"] == kij and r["diameter"] == d] for d in ("sigma", "unified-0.88")}
    sel = by["sigma"]
    x = range(len(sel))
    value = [float(r["experiment_kJ_mol"]) for r in sel]
    yerr = [[v - float(r["experiment_lower_kJ_mol"]) for v, r in zip(value, sel)],
            [float(r["experiment_upper_kJ_mol"]) - v for v, r in zip(value, sel)]]
    ax.bar([i - 1.5 * w for i in x], value, w, color="#d9d9d9", edgecolor="k", yerr=yerr, capsize=2, label="experiment (paper)")
    ax.bar([i - 0.5 * w for i in x], [float(r["paper_model_kJ_mol"]) for r in sel], w, color=PAPER, label="paper model")
    for k, (diameter, color, label) in enumerate(CONVENTIONS):
        ax.bar([i + (0.5 + k) * w for i in x], [float(r["engine_kJ_mol"]) for r in by[diameter]], w, color=color, label=label)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(list(x), [f"{r['ion']}\n{r['solvent'][:4]}" for r in sel], fontsize=8)
    ax.set_title("ion–alcohol $k_{ij}=0$ (prediction)" if kij == "prediction" else "Table 8 $k_{ij}$ (correlation)")
axes[0].set_ylabel(r"$\Delta^\mathrm{Tr}G_i$ water → alcohol / kJ mol$^{-1}$")
axes[1].legend(fontsize=7)
fig.savefig(OUT / "fig8b-transfer-gibbs.png", dpi=150)

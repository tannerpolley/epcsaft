"""Render retained derivative refinement; never execute the thermodynamic model."""

import csv
import hashlib
from pathlib import Path

import matplotlib.pyplot as plt

campaign = Path(__file__).resolve().parents[1]
source = campaign / "results/live/ladders.csv"
with source.open() as stream:
    rows = list(csv.DictReader(stream))
fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True, layout="constrained")
for order, ax in enumerate(axes, 1):
    for case, color, marker in (("anchor", "#0072B2", "o"), ("optimum", "#D55E00", "x")):
        selected = [r for r in rows if r["case"] == case and int(r["order"]) == order]
        steps = sorted({float(r["step"]) for r in selected}, reverse=True)
        ratios = [max(float(r["error"]) / float(r["limit"]) for r in selected if float(r["step"]) == h) for h in steps]
        ax.loglog(steps, ratios, marker=marker, color=color, label=case)
    ax.axhline(1, color="0.25", linestyle="--", linewidth=1, label="verification threshold")
    ax.invert_xaxis()
    ax.set(xlabel="Dimensionless perturbation h", title="Gradient" if order == 1 else "Hessian-vector")
    ax.grid(which="major", alpha=0.2)
axes[0].set_ylabel("Worst error / verification threshold")
axes[1].legend(fontsize=8)
fig.suptitle("Refinement resolves exact equilibrium derivatives")
fig.savefig(campaign / "figures/live-derivatives.png", dpi=180, metadata={
    "Description": "Worst of ten outputs and eight directions; numerical refinement, not experimental uncertainty.",
    "SourceSHA256": hashlib.sha256(source.read_bytes()).hexdigest(),
})
plt.close(fig)

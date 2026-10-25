"""Plot retained dimensionless first-derivative comparison errors."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

study = Path(__file__).resolve().parents[1]
rows = list(csv.DictReader((study / "results/derivatives.csv").open()))
cases = list(dict.fromkeys(row["comparison"].rsplit("/", 3)[0] for row in rows))
fig, ax = plt.subplots(figsize=(8, 4.5), layout="constrained")
for rung, offset, marker, color in (("fine", -0.08, "o", "#0072B2"), ("coarse", 0.08, "x", "#D55E00")):
    worst = [max(float(row["absolute_error"]) / float(row["allowed_error"]) for row in rows
                 if row["comparison"].rsplit("/", 3)[0] == case and row["comparison"].endswith("/" + rung))
             for case in cases]
    ax.scatter([i + offset for i in range(len(cases))], worst, marker=marker, color=color,
               label="0.0025 optimizer units" if rung == "fine" else "0.005 optimizer units", s=55)
ax.axhline(1, color="black", linewidth=1, label="Required limit")
ax.set_yscale("log")
ax.set_ylim(1e-6, 2)
ax.set_xticks(range(len(cases)), [case.replace("reaction-", "").replace("/", "\n").replace("_", " ") for case in cases])
ax.set_ylabel("Largest absolute error / allowed error (dimensionless)")
ax.set_title("Exact first columns satisfy both central-difference rungs")
ax.legend(loc="lower left", fontsize=9)
ax.spines[["top", "right"]].set_visible(False)
(study / "figures").mkdir(exist_ok=True)
fig.savefig(study / "figures/derivative-errors.png", dpi=180)

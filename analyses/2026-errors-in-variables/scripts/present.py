"""Plot the signed means and intervals retained by the native driver."""
import csv, pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
leaf = Path(__file__).resolve().parents[1]; (leaf / "figures").mkdir(exist_ok=True)
fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), layout="constrained")
for ax, case in zip(axes, ("mea", "ionic")):
    rows = list(csv.DictReader((leaf / "results" / case / "summary.csv").open())); assert len(rows) == 2
    for stem in ("summary", "decision", "numerical-components"): (leaf / "figures" / f"{case}-{stem}.html").write_text(pd.read_csv(leaf / "results" / case / f"{stem}.csv").T.to_html(header=False, float_format="{:.9g}".format))
    for i, row in enumerate(rows):
        ax.plot([float(row["signed_lower"]), float(row["signed_upper"])], [i, i], color=("#0072B2", "#D55E00")[i])
        ax.scatter(float(row["signed_mean"]), i, marker=("o", "x")[i], color=("#0072B2", "#D55E00")[i])
    ax.axvline(0, color="0.4", linestyle=":"); ax.set(title=case.upper(), xlabel="Signed mean parameter error\n/ coordinate scale", yticks=[0, 1], yticklabels=["Fixed inputs", "Input and output errors"], ylim=(-0.5, 1.5))
fig.savefig(leaf / "figures/signed-mean-intervals.png", dpi=180)
pd.read_csv(leaf / "results/requalification.csv").replace({"holds": {1: "holds", 0: "NOT MET"}}).fillna("").to_html(leaf / "figures/requalification.html", index=False, float_format="{:.9g}".format)

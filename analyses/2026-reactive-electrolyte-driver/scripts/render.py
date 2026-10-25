"""Render retained same-wheel exact-response timings; never run the Engine."""

import csv
import math
from pathlib import Path
from statistics import median, quantiles

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
with (ROOT / "results/timings.csv").open() as stream:
    rows = list(csv.DictReader(stream))
methods = [
    ("exact_three_directions_legacy", "Three separate requests", "#0072B2"),
    ("exact_three_directions", "One shared-center batch", "#D55E00"),
]
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(6.2, 3.8), layout="constrained")
for position, (method, _, color) in enumerate(methods):
    samples = [float(row["wall_ms"]) for row in rows if row["method"] == method]
    assert len(samples) >= 2 and all(math.isfinite(value) for value in samples)
    mid = median(samples)
    percentiles = quantiles(samples, n=10, method="inclusive")
    offsets = [position + 0.035 * (i - (len(samples) - 1) / 2) for i in range(len(samples))]
    ax.scatter(offsets, samples, color=color, s=28,
               label="Individual samples" if position == 0 else None)
    ax.errorbar(position, mid, yerr=[[mid - percentiles[0]], [percentiles[-1] - mid]],
                fmt="x", color="#222222", capsize=6, markersize=8,
                label="Median and p10–p90" if position == 0 else None)
    ax.annotate(f"{mid:.3f} ms", (position, max(samples)), xytext=(0, 9),
                textcoords="offset points", ha="center")
ax.set(xticks=[0, 1], xticklabels=[item[1] for item in methods],
       ylabel="Wall time (ms)", title="Same three exact response directions", xlim=(-0.45, 1.45))
ax.set_ylim(0, ax.get_ylim()[1] * 1.18)
ax.grid(axis="y", alpha=0.2)
ax.legend(loc="lower left", frameon=False)
(ROOT / "figures").mkdir(exist_ok=True)
fig.savefig(ROOT / "figures/first-response-runtime.png", dpi=160)
plt.close(fig)

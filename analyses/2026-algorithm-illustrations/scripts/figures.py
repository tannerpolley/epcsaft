"""Paper figures from retained native values; no model calls."""

import json
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parents[1]


def save(fig, name):
    fig.text(0.5, 0.01, "Numerical illustration of the algorithms; no physical validation", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    for suffix in ("pdf", "png"):
        fig.savefig(HERE / "figures" / f"{name}.{suffix}", dpi=240)


e = pd.read_csv(HERE / "results/envelope.csv")
s = json.loads((HERE / "inputs/study.json").read_text())
ref = pd.read_csv(HERE.parents[1] / s["envelope_reference"])
k = json.loads((HERE / "results/trace.json").read_text())["critical"][0]
fig, ax = plt.subplots(figsize=(6, 4))
for bubble, color, style, label in zip((True, False), ("#0072B2", "#D55E00"), ("-", "--"), ("Bubble", "Dew")):
    part = e[e.bubble == bubble]
    ax.plot(part.T_K, part.P_Pa / 1e6, style, color=color, label=label)
ax.scatter(
    ref.T_K, ref.P_Pa / 1e6, marker="o", facecolors="none", edgecolors="black", label="Retained FeOs points", zorder=4
)
ax.plot(k["T_K"], k["P_Pa"] / 1e6, "*", color="black", ms=10, label="Critical candidate")
ax.set(xlabel="Temperature / K", ylabel="Pressure / MPa", title="Methane/ethane: feed mole fractions (0.5, 0.5)")
ax.legend(fontsize=8)
save(fig, "F4-phase-envelope")
t, p = pd.read_csv(HERE / "results/tpd.csv"), pd.read_csv(HERE / "results/search.csv")
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(t.x_water, t.tpd_over_rt, color="#0072B2", label="Native supplied liquid trials")
ax.scatter(p.x_water, p.tpd_over_rt, marker="x", color="#D55E00", label="Native liquid trial search", zorder=4)
ax.axhline(0, color="black", lw=0.7)
ax.set(
    xscale="logit",
    xlabel="Water mole fraction (logit scale)",
    ylabel=r"Tangent-plane distance / $RT$",
    title="Water/hexane: 300 K, 0.1 MPa",
)
ax.legend(fontsize=8)
save(fig, "F3-tpd-path")

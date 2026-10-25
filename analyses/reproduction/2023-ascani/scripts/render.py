"""Plot the retained system-1 comparison in Senina's transformed variables (alpha_1 = x_acid + x_water,
alpha_2 = x_alcohol + x_water) and the acid content of the water-rich liquid. Reads results only."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
X = ("x_acetic_acid", "x_1_pentanol", "x_pentyl_acetate", "x_water")
EXP, MODEL = "#555555", "#1f6fb4"


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def alpha(row, prefix=""):
    x = [float(row[prefix + k]) for k in X]
    return x[0] + x[3], x[1] + x[3]


def main():
    hom = [r for r in rows(ROOT / "results/homogeneous-ce.csv") if r["rule"] == "printed"]
    ties = [r for r in rows(ROOT / "results/reactive-tie-lines.csv") if r["rule"] == "printed"]
    exp_hom = rows(ROOT / "inputs/data/senina-2022-table-2-homogeneous-ce.csv")
    exp_ties = rows(ROOT / "inputs/data/senina-2022-table-3-reactive-tie-lines.csv")
    fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4.4), gridspec_kw={"width_ratios": [1.25, 1]})
    for e in exp_ties:
        (a1, a2), (b1, b2) = alpha(e, "aq_"), alpha(e, "org_")
        left.plot([a1, b1], [a2, b2], "--o", color=EXP, mfc="white", ms=4, lw=0.8)
    for t in ties:
        (a1, a2), (b1, b2) = alpha(t, "aq_"), alpha(t, "org_")
        left.plot([a1, b1], [a2, b2], "-s", color=MODEL, ms=3.5, lw=0.9)
    left.plot(*zip(*(alpha(r) for r in exp_hom)), "o", color=EXP, mfc="white", ms=6, label="experiment")
    for r in hom:  # overall composition; two liquids where the model splits
        a = (float(r["phase1_x_acetic_acid"]) + float(r["phase1_x_water"]), float(r["phase1_x_1_pentanol"]) + float(r["phase1_x_water"])) \
            if r["phases"] == "1" else None
        if a:
            left.plot(*a, "s", color=MODEL, ms=3.5)
        else:
            left.plot(*alpha(r, "phase1_"), "s", color=MODEL, mfc="white", ms=4)
            left.plot(*alpha(r, "phase2_"), "s", color=MODEL, mfc="white", ms=4)
    left.plot([], [], "--o", color=EXP, mfc="white", label="experimental tie-line")
    left.plot([], [], "-s", color=MODEL, label="PC-SAFT, Engine")
    left.plot([], [], "s", color=MODEL, mfc="white", label="PC-SAFT two liquids at a measured one-liquid point")
    left.set(xlabel=r"$\alpha_1 = x_\mathrm{acid} + x_\mathrm{water}$", ylabel=r"$\alpha_2 = x_\mathrm{alcohol} + x_\mathrm{water}$",
             xlim=(0, 1.02), ylim=(0, 1.02), title="(a) chemical equilibrium surface, 318.15 K, 1 bar")
    fig.legend(*left.get_legend_handles_labels(), frameon=False, fontsize=8, loc="lower center", ncol=4)
    exp_by = {e["tie_line"]: e for e in exp_ties}
    xe = [float(exp_by[t["tie_line"]]["aq_x_acetic_acid"]) for t in ties]
    xm = [float(t["aq_x_acetic_acid"]) for t in ties]
    right.plot([0, 0.32], [0, 0.32], color="#bbbbbb", lw=0.8)
    right.plot(xe, xm, "s", color=MODEL)
    right.set(xlabel=r"experimental $x_\mathrm{acid}$ in the water-rich liquid",
              ylabel=r"PC-SAFT $x_\mathrm{acid}$ in the water-rich liquid", xlim=(0, 0.32), ylim=(0, 0.32),
              title="(b) acid in the water-rich liquid")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    (ROOT / "figures").mkdir(exist_ok=True)
    fig.savefig(ROOT / "figures/system-1-reacting-lle.png", dpi=160)


if __name__ == "__main__":
    main()

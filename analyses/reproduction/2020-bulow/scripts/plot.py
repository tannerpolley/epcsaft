"""Render retained results only; no Engine calls."""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
M = {"methanol": 0.03204186, "ethanol": 0.04606844}


def read(path):
    with (HERE / path).open() as f:
        return list(csv.DictReader(f))


def miac():
    curves, points, han = read("results/miac-curves.csv"), read("results/miac-digitized-points.csv"), read("results/miac-han-pan.csv")
    panels = [("figure-07/libr-methanol", "LiBr", "methanol", 0.2), ("figure-07/nabr-ethanol", "NaBr", "ethanol", 0.010),
              ("figure-07/licl-methanol", "LiCl", "methanol", 0.15), ("figure-07/nabr-methanol", "NaBr", "methanol", 0.05)]
    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    for ax, (panel, salt, solvent, xmax) in zip(axes.flat, panels):
        digit = read(f"figures/{panel}.csv")
        for variant, color in (("advanced", "tab:green"), ("revised", "black")):
            c = [r for r in curves if r["variant"] == variant and r["salt"] == salt and r["solvent"] == solvent]
            ax.plot([float(r["x_salt"]) for r in c], [float(r["gamma_pm_x"]) for r in c], color=color, lw=1.6,
                    label=f"Engine {variant}, $\\gamma^*_\\pm$")
            if variant == "advanced":
                ax.plot([float(r["x_salt"]) for r in c], [float(r["gamma_pm_x_times_factor"]) for r in c], color=color,
                        lw=1.2, ls="--", label="Engine advanced $\\times(1+\\nu m M_s)$")
            d = sorted((float(r["x"]), float(r["y"])) for r in digit if r["series"] == variant)
            ax.plot(*zip(*d), color=color, lw=0, marker=".", ms=1.5, alpha=0.5, label=f"corrigendum {variant} (digitized)")
        p = [r for r in points if r["panel"] == panel]
        ax.plot([float(r["x_salt"]) for r in p], [float(r["gamma_pm_x_digitized"]) for r in p], "o", mfc="none", color="tab:red",
                label="plotted data (digitized)")
        h = [r for r in han if r["solvent"] == solvent and salt == "NaBr"]
        if h:
            m = [float(r["molality_mol_kg"]) for r in h]
            ax.plot([mi * M[solvent] / (1 + mi * M[solvent]) for mi in m],
                    [float(r["gamma_pm_m_han_pan"]) * (1 + 2 * mi * M[solvent]) for r, mi in zip(h, m)], "x", color="tab:blue",
                    label="Han & Pan 1993 (tabulated)")
        ax.set(xlim=(0, xmax), ylim=(0, None), xlabel=r"$x_\mathrm{salt}$ / mol mol$^{-1}$", ylabel=r"$\gamma^*_\pm$",
               title=f"{salt} in {solvent}, 298.15 K")
    axes.flat[0].legend(fontsize=6.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(HERE / "figures/miac-comparison.png", dpi=150)


def infinite_dilution():
    hyd, tr = read("results/hydration.csv"), read("results/transfer.csv")
    fig, (a, b) = plt.subplots(1, 2, figsize=(9, 4))
    for key, marker in (("advanced", "o"), ("revised", "s")):
        a.plot([float(r[f"paper_{key}"]) for r in hyd], [float(r[f"engine_{key}"]) for r in hyd], marker, mfc="none", label=key)
        b.plot([float(r[f"paper_{key}"]) for r in tr], [float(r[f"engine_{key}"]) for r in tr], marker, mfc="none", label=key)
    for ax, lim, title in ((a, (-900, 50), "hydration, Fig. 2"), (b, (-12, 32), "water → alcohol transfer, Fig. 4")):
        ax.plot(lim, lim, color="0.6", lw=0.8)
        ax.set(xlim=lim, ylim=lim, xlabel="paper (digitized) / kJ mol$^{-1}$", ylabel="Engine / kJ mol$^{-1}$", title=title)
        ax.legend()
    fig.tight_layout()
    fig.savefig(HERE / "figures/infinite-dilution-parity.png", dpi=150)


if __name__ == "__main__":
    miac()
    infinite_dilution()

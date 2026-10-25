"""Plot retained Held 2014 comparisons; reads results/ and inputs/ only (no Engine)."""
import csv
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from model import HERE


def read(path):
    return list(csv.DictReader((HERE / path).open()))


def by(rows, *keys):
    out = defaultdict(list)
    for r in rows:
        out[tuple(r[k] for k in keys)].append(r)
    return out


def osmotic_figure():
    curves = by(read("results/figure-curves.csv"), "figure", "series")
    lines = by(read("inputs/paper-lines.csv"), "figure", "series")
    hw = by(read("inputs/hamer-wu-1972-osmotic.csv"), "salt")
    cold = by(read("results/figure-2-cold-markers.csv"), "salt")
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    panels = [("figure-1", "Fig. 1, 298.15 K", {"NaCl": "NaCl", "KBr": "KBr"}),
              ("figure-2", "Fig. 2, strategy 2", {"NaCl": "NaCl", "KCl": "KCl"}),
              ("figure-3b", "Fig. 3b, strategy 2, 298.15 K", {"LiAc": "LiAc", "NaAc": "NaAc", "KAc": "KAc"})]
    for ax, (figure, title, salts) in zip(axes.flat, panels):
        for i, (key, rows) in enumerate(sorted((k, v) for k, v in curves.items() if k[0] == figure)):
            color = f"C{i}"
            m = [float(r["molality_mol_kg"]) for r in rows if float(r["molality_mol_kg"]) >= 0.05]
            ax.plot(m, [float(r["osmotic_engine"]) for r in rows if float(r["molality_mol_kg"]) >= 0.05], color=color,
                    label=f"{key[1]}: Engine, closest approach 0.88 σ")
            ax.plot(m, [float(r["osmotic_engine_sigma_closest_approach"]) for r in rows
                        if float(r["molality_mol_kg"]) >= 0.05], color=color, ls=":", lw=1)
            paper = sorted((float(r["x"]), float(r["y"])) for r in lines[key])
            ax.plot(*zip(*paper), ls="none", marker=".", ms=2, color="k")
        for salt in salts:
            points = [(float(r["molality_mol_kg"]), float(r["osmotic_coefficient"])) for r in hw[(salt,)]
                      if 0.05 <= float(r["molality_mol_kg"]) <= 4]
            ax.plot(*zip(*points), ls="none", marker="o", mfc="none", color="0.4", ms=4)
        if figure == "figure-2":
            for salt, rows in cold.items():
                ax.plot([float(r["molality_mol_kg"]) for r in rows], [float(r["osmotic_digitized"]) for r in rows],
                        ls="none", marker="s", mfc="none", color="0.6", ms=4)
        ax.set(title=title, xlabel="molality / mol kg$^{-1}$", ylabel="osmotic coefficient Φ")
        ax.legend(fontsize=6)
    ax = axes[1, 1]
    table = read("results/table-5-osmotic-ard.csv")
    for dh, marker, label in (("0.88", "o", "closest approach 0.88 σ"), ("1.0", "x", "closest approach σ")):
        rows = [r for r in table if r["dh_factor"] == dh and r["strategy"] == "2"]
        ax.plot([float(r["ard_printed_percent"]) for r in rows], [float(r["ard_engine_percent"]) for r in rows],
                ls="none", marker=marker, label=label)
    ax.plot([0, 5], [0, 5], color="k", lw=0.8)
    ax.set(title="Table 5, strategy 2, Φ ARD", xlabel="printed ARD / %", ylabel="Engine vs Hamer–Wu ARD / %",
           xlim=(0, 5), ylim=(0, 5))
    ax.legend(fontsize=7)
    fig.suptitle("Engine (lines) against digitized Held 2014 model lines (black dots); open circles Hamer–Wu 1972, "
                 "open squares 273.15 K figure data; dotted: σ closest approach", fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / "figures/osmotic-coefficients.png", dpi=150)


def ternary_figure():
    lines = by(read("inputs/paper-lines.csv"), "figure", "series")
    markers = by(read("inputs/paper-markers.csv"), "figure")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    rows = by(read("results/figure-4a.csv"), "solute")
    ax = axes[0]
    for i, solute in enumerate(("glycine", "alanine")):
        r = rows[(solute,)]
        ax.plot([float(x["m_amino_acid_mol_kg"]) for x in r], [float(x["osmotic_python_contraction"]) for x in r], color=f"C{i}",
                ls="--", label=f"{solute}: Python contraction (no native Φ)")
        paper = sorted((float(p["x"]), float(p["y"])) for p in lines[("figure-4a", f"{solute} strategy 2")])
        ax.plot(*zip(*paper), ls="none", marker=".", ms=2, color="k")
    exp = markers[("figure-4a",)]
    ax.plot([float(p["x"]) for p in exp], [float(p["y"]) for p in exp], ls="none", marker="s", mfc="none", color="0.5")
    ax.set(title="Fig. 4a, 1 mol/kg KCl, 298.15 K", xlabel="amino-acid molality / mol kg$^{-1}$", ylabel="Φ")
    ax.legend(fontsize=7)
    ax = axes[1]
    rows = by(read("results/figure-5.csv"), "salt")
    for i, salt in enumerate(("NaCl", "NaI")):
        r = rows[(salt,)]
        cut = 3.6 if salt == "NaCl" else 1.2
        r = [x for x in r if float(x["m_salt_mol_kg"]) <= cut]
        ax.plot([float(x["m_salt_mol_kg"]) for x in r], [float(x["m_benzene_mmol_kg"]) for x in r], color=f"C{i}",
                label=f"{salt}: Engine")
        paper = sorted((float(p["x"]), float(p["y"])) for p in lines[("figure-5", f"{salt} strategy 2")])
        ax.plot(*zip(*paper), ls="none", marker=".", ms=2, color="k")
    exp = markers[("figure-5",)]
    ax.plot([float(p["x"]) for p in exp], [float(p["y"]) for p in exp], ls="none", marker="o", mfc="none", color="0.5")
    ax.set(title="Fig. 5, benzene in aqueous salt, 298.15 K", xlabel="salt molality / mol kg$^{-1}$",
           ylabel="benzene molality / mmol kg$^{-1}$")
    ax.legend(fontsize=7)
    ax = axes[2]
    r = read("results/figure-6a.csv")
    w = [float(x["w_nh4cl_aqueous"]) for x in r]
    ax.plot(w, [float(x["w_butanol_organic"]) for x in r], color="C0", label="1-butanol phase: Engine")
    ax.plot(w, [float(x["w_butanol_aqueous"]) for x in r], color="C1", label="aqueous phase: Engine")
    for series in ("1-butanol phase", "aqueous phase"):
        paper = sorted((float(p["x"]), float(p["y"])) for p in lines[("figure-6a", series)])
        ax.plot(*zip(*paper), ls="none", marker=".", ms=2, color="k")
    exp = markers[("figure-6a",)]
    ax.plot([float(p["x"]) for p in exp], [float(p["y"]) for p in exp], ls="none", marker="o", mfc="none", color="0.5")
    ax.set(title="Fig. 6a, 1-butanol + NH$_4$Cl + water, 298.15 K", xlabel="NH$_4$Cl in aqueous phase / wt%",
           ylabel="1-butanol / wt%")
    ax.legend(fontsize=7)
    fig.suptitle("Engine (lines), digitized Held 2014 model lines (black dots), digitized figure data (open symbols)",
                 fontsize=8)
    fig.tight_layout()
    fig.savefig(HERE / "figures/ternary-systems.png", dpi=150)


if __name__ == "__main__":
    (HERE / "figures").mkdir(exist_ok=True)
    osmotic_figure()
    ternary_figure()

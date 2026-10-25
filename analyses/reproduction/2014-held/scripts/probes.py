"""Retained diagnostic probes for the unresolved Fig. 3a and Fig. 5 NaCl discrepancies (all native).

Fig. 3a: strategy-1 acetate gamma+- against Hamer-Wu 1972, to compare with Held and Sadowski 2009
Table 3 classical deviations (a gamma+- check; it does not verify a Phi curve).
Fig. 5 NaCl: aqueous benzene molality relative to the traced published line at four salt molalities
for the Table 4 Cl- -benzene k_ij, two probe values, and all ion-benzene k_ij = 0.
"""
import csv

import model
from run import Salt, write
from ternary import lle, mixture, MW
from model import HERE


def acetate_gamma():
    observed = [r for r in csv.DictReader((HERE / "inputs/hamer-wu-1972-osmotic.csv").open())
                if r["salt"] in ("LiAc", "NaAc", "KAc") and float(r["molality_mol_kg"]) <= 5.0]
    salts = {name: Salt(1, name[:-2] + "+", "Ac-", 298.15) for name in ("LiAc", "NaAc", "KAc")}
    return [dict(salt=r["salt"], strategy=1, molality_mol_kg=float(r["molality_mol_kg"]),
                 gamma_engine=salts[r["salt"]].reference.activity(float(r["molality_mol_kg"])).mean_ionic_molality_scale,
                 gamma_hamer_wu=float(r["gamma_pm_m"])) for r in observed]


def benzene_kij(lines):
    table4 = dict(model.TABLE4)
    rows = []
    for label, change in (("Table 4 (Cl- 0.35)", {}), ("Cl- 0.45", {("Cl-", "benzene"): (0.45, 0.0)}),
                          ("Cl- 0.50", {("Cl-", "benzene"): (0.5, 0.0)}),
                          ("all ion-benzene 0", {("Na+", "benzene"): (0.0, 0.0), ("Cl-", "benzene"): (0.0, 0.0)})):
        model.TABLE4.update(change)
        mix, comp = mixture("Na+", "Cl-", "benzene")
        model.TABLE4.clear(); model.TABLE4.update(table4)
        for m in (0.5, 1.5, 2.5, 3.5):
            n_w = 1 / MW
            total = n_w + 2 * m
            r = lle(mix, comp, (n_w, 0.5, m, m), [n_w / total * 0.999, 4e-4, m / total, m / total],
                    [3e-3, 0.997], comp[:2], None)
            x = r.mole_fractions[:4]
            m_aq, benzene = x[2] / (x[0] * MW), 1000 * x[1] / (x[0] * MW)
            paper = sorted(lines)
            line = next(y0 + (y1 - y0) * (m_aq - x0) / (x1 - x0)
                        for (x0, y0), (x1, y1) in zip(paper, paper[1:]) if x0 <= m_aq <= x1)
            rows.append(dict(probe=label, m_salt_feed_mol_kg=m, m_salt_mol_kg=m_aq, m_benzene_mmol_kg=benzene,
                             m_benzene_paper_line_mmol_kg=line, relative_difference=benzene / line - 1))
    return rows


def main():
    lines = [(float(r["x"]), float(r["y"])) for r in csv.DictReader((HERE / "inputs/paper-lines.csv").open())
             if r["figure"] == "figure-5" and r["series"] == "NaCl strategy 2"]
    write("probe-figure-3a-gamma.csv", acetate_gamma())
    write("probe-figure-5-kij.csv", benzene_kij(lines))


if __name__ == "__main__":
    main()

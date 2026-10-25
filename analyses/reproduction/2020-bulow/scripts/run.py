"""Bülow, Ascani & Held 2020 (ePC-SAFT advanced I) through the Engine: Figs. 2-7 with the corrigendum.

The Engine computes every fugacity, density and activity natively; Python declares states, reads results
and compares them. Tolerances (TOL) were fixed before the first comparison run.
"""
import csv
import hashlib
import json
import math
import sys
import time
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro, R, elementary_charge

import epcsaft

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "scripts"))
import parameters as P  # noqa: E402

T, PRESSURE = 298.15, 1.0e5
RT = R * T / 1000  # kJ/mol
# CODATA 2018, as in engine/native/core/physical_constants.hpp; SciPy now ships CODATA 2022 (8.8541878188e-12),
# whose 6.8e-10 relative difference would exceed the closed-form tolerance by itself.
EPSILON_0 = 8.8541878128e-12
SOLVENT_M = {s: P.COMPONENTS[s][1] / 1000 for s in ("water", "methanol", "ethanol")}  # kg/mol
ION = {"Li+": "lithium-cation", "Na+": "sodium-cation", "K+": "potassium-cation", "F-": "fluoride-anion",
       "Cl-": "chloride-anion", "Br-": "bromide-anion", "I-": "iodide-anion"}
PARTNER = {"Li+": "Cl-", "Na+": "Cl-", "K+": "Cl-", "F-": "Na+", "Cl-": "Na+", "Br-": "Na+", "I-": "Na+"}
SALTS = {"LiBr": ("Li+", "Br-"), "NaBr": ("Na+", "Br-"), "LiCl": ("Li+", "Cl-")}
# MIAC panels of the corrigendum (digitized curve file, salt, solvent, x_salt upper end of the plotted axis)
PANELS = {"figure-07/libr-methanol": ("LiBr", "methanol", 0.2), "figure-07/nabr-ethanol": ("NaBr", "ethanol", 0.010),
          "figure-07/licl-methanol": ("LiCl", "methanol", 0.15), "figure-07/nabr-methanol": ("NaBr", "methanol", 0.05),
          "figure-06/libr-methanol-left": ("LiBr", "methanol", 0.2)}
TOL = {  # stated before running; kJ/mol unless noted
    "figure-2 hydration": 10.0,           # bar heights on a ~1000 kJ/mol axis
    "figure-3 Born": 5.0,                 # packet low/high read spread is about +-4.5
    "figure-4 transfer": 1.5,             # ~0.15 kJ/mol per pixel; bar edges and the 0-line
    "figure-5 Born transfer": 1.5,
    "curve pixel distance": 3.0,          # px, orthogonal distance from a digitized curve pixel to the Engine curve
    "Born closed form relative": 1e-10,   # numerical identity at infinite dilution
    "partner independence": 1e-3,         # kJ/mol, ion ln phi_inf from two different counter-ions
}


def write(name, rows):
    with (HERE / "results" / name).open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def read(path):
    with (HERE / path).open() as f:
        return list(csv.DictReader(f))


def document(variant):
    if variant != "electrostatics-off":
        return P.document(variant)
    # Diagnostic only: one equal, very large permittivity makes Born cancel in gamma and DH vanish, leaving
    # hard-chain + dispersion + association; it splits the Fig. 6 right-panel contributions.
    doc = P.document("dh-mixed-permittivity")
    for component in doc["components"]:
        for c in component["coefficients"]:
            if c["family"] == "relative_permittivity":
                c["value"]["magnitude"] = 1e8
    return doc


@lru_cache(maxsize=None)
def mixture(variant, components):
    return epcsaft.Mixture(epcsaft.Parameters.from_mapping(document(variant), components=list(components)))


@lru_cache(maxsize=None)
def reference(variant, solvent, cation, anion):
    model = mixture(variant, (solvent, ION[cation], ION[anion]))
    return model.reference(T, PRESSURE, solvent_fractions=[1, 0, 0], ions=(ION[cation], ION[anion]), stoichiometry=(1, 1))


def ln_phi_inf(variant, solvent, ion, partner=None):
    partner = partner or PARTNER[ion]
    cation, anion = (ion, partner) if ion.endswith("+") else (partner, ion)
    ref = reference(variant, solvent, cation, anion)
    return ref.log_fugacity[0 if ion == cation else 1], ref.convergence_error


def born_closed_form(ion, solvent):
    """Eq. (2) insertion limit: -(e^2/(4 pi eps0 kT)) (1 - 1/eps_s) z^2/d_Born, per mole in kJ/mol."""
    d = P.PACKING * P.COMPONENTS[ION[ion]][3] * 1e-10
    eps = P.COMPONENTS[solvent][5]
    return -Avogadro * elementary_charge**2 / (4 * math.pi * EPSILON_0 * d) * (1 - 1 / eps) / 1000


def infinite_dilution():
    paper2 = {r["x"]: r for r in read("figures/figure-02/source/water-comparisons.csv")}
    paper3 = {r[""]: r for r in read("figures/figure-03/source/water-contributions.csv")}
    paper4 = {s: {r[""]: r for r in read(f"figures/figure-04/source/water-{s}-comparison.csv")} for s in ("methanol", "ethanol")}
    rows, hyd, transfer, checks = [], [], [], []
    variants = ("advanced", "revised", "dh-mixed-permittivity", "born-diameter-sigma", "sodium-water-held2014", "advanced-dh-sigma")
    lnphi = {}
    for variant in variants:
        for solvent in SOLVENT_M:
            for ion in ION:
                value, conv = ln_phi_inf(variant, solvent, ion)
                lnphi[variant, solvent, ion] = value
                rows.append(dict(variant=variant, solvent=solvent, ion=ion, counter_ion=PARTNER[ion], ln_phi_inf=value,
                                 solvation_gibbs_kJ_mol=RT * value, reference_terminal_change=conv))
    write("infinite-dilution.csv", rows)
    for ion in ION:
        born = RT * (lnphi["advanced", "water", ion] - lnphi["dh-mixed-permittivity", "water", ion])
        hyd.append(dict(ion=ion,
                        engine_advanced=RT * lnphi["advanced", "water", ion], paper_advanced=float(paper2["advanced"][ion]),
                        engine_revised=RT * lnphi["revised", "water", ion], paper_revised=float(paper2["revised"][ion]),
                        fawcett=float(paper2["data"][ion]), saft_vre=float(paper2["SAFT-VR"][ion]),
                        engine_born_diameter_sigma=RT * lnphi["born-diameter-sigma", "water", ion],
                        engine_sodium_water_held2014=RT * lnphi["sodium-water-held2014", "water", ion],
                        engine_born=born, paper_born=float(paper3["born avg"][ion]),
                        paper_born_low=float(paper3["born low"][ion]), paper_born_high=float(paper3["born hi"][ion]),
                        engine_non_born=RT * lnphi["advanced", "water", ion] - born,
                        paper_hc_disp_assoc=sum(float(paper3[f"{k} avg"][ion]) for k in ("hc", "disp", "assoc"))))
        for solvent in ("methanol", "ethanol"):
            data = paper4[solvent]
            b = RT * ((lnphi["advanced", solvent, ion] - lnphi["dh-mixed-permittivity", solvent, ion])
                      - (lnphi["advanced", "water", ion] - lnphi["dh-mixed-permittivity", "water", ion]))
            transfer.append(dict(solvent=solvent, ion=ion,
                                 engine_advanced=RT * (lnphi["advanced", solvent, ion] - lnphi["advanced", "water", ion]),
                                 paper_advanced=float(data["advanced"][ion]),
                                 engine_revised=RT * (lnphi["revised", solvent, ion] - lnphi["revised", "water", ion]),
                                 paper_revised=float(data["revised"][ion]),
                                 marcus_low=float(data["data low"][ion]), marcus=float(data["data median"][ion]),
                                 marcus_high=float(data["data high"][ion]),
                                 engine_born_diameter_sigma=RT * (lnphi["born-diameter-sigma", solvent, ion] - lnphi["born-diameter-sigma", "water", ion]),
                                 engine_born=b))
        for solvent in SOLVENT_M:  # Born closed form at infinite dilution: identity, not a source comparison
            b = RT * (lnphi["advanced", solvent, ion] - lnphi["dh-mixed-permittivity", solvent, ion])
            checks.append(dict(check="Born closed form", solvent=solvent, ion=ion, engine=b, reference=born_closed_form(ion, solvent),
                               relative_error=abs(b / born_closed_form(ion, solvent) - 1), tolerance=TOL["Born closed form relative"]))
    for ion, other in (("Na+", "I-"), ("Cl-", "Li+"), ("Br-", "K+"), ("K+", "Br-")):
        for solvent in SOLVENT_M:
            a, b = ln_phi_inf("advanced", solvent, ion)[0], ln_phi_inf("advanced", solvent, ion, other)[0]
            checks.append(dict(check="counter-ion independence", solvent=solvent, ion=f"{ion} ({PARTNER[ion]} vs {other})",
                               engine=RT * a, reference=RT * b, relative_error=abs(RT * (a - b)), tolerance=TOL["partner independence"]))
    write("hydration.csv", hyd)
    write("transfer.csv", transfer)
    paper5 = {s: {r["x"]: r for r in read(f"figures/figure-05/source/water-{s}-contributions.csv")} for s in ("methanol", "ethanol")}
    contributions = []
    for label, data in paper5.items():
        for ion in ("Na+", "Cl-", "I-"):
            engine = {s: next(r["engine_born"] for r in transfer if r["solvent"] == s and r["ion"] == ion) for s in ("methanol", "ethanol")}
            contributions.append(dict(packet_file_label=label, ion=ion, paper_born=float(data["Born"][ion]),
                                      paper_hc_disp_assoc=sum(float(data[k][ion]) for k in ("hc", "disp", "assoc")),
                                      engine_born_methanol=engine["methanol"], engine_born_ethanol=engine["ethanol"]))
    write("transfer-born-contributions.csv", contributions)
    return checks


def x_to_m(x, solvent):
    return x / ((1 - x) * SOLVENT_M[solvent])


def activity(variant, salt, solvent, m):
    cation, anion = SALTS[salt]
    ref = reference(variant, solvent, cation, anion)
    return ref.activity(m)


def plotted_factor(m, solvent):
    """1 + nu m M_s: the corrigendum's advanced curves equal the Engine gamma_x times this factor (README)."""
    return 1 + 2 * m * SOLVENT_M[solvent]


def curve(variant, salt, solvent, x_max, n=240):
    xs = np.unique(np.concatenate([np.geomspace(1e-6, x_max, n), np.linspace(x_max / n, x_max, n)]))
    out = []
    for x in xs:
        m = x_to_m(x, solvent)
        a = activity(variant, salt, solvent, m)
        out.append(dict(variant=variant, salt=salt, solvent=solvent, x_salt=x, molality_mol_kg=m,
                        gamma_pm_x=a.mean_ionic_mole_fraction_scale, gamma_pm_m=a.mean_ionic_molality_scale,
                        gamma_pm_x_times_factor=a.mean_ionic_mole_fraction_scale * plotted_factor(m, solvent),
                        osmotic_coefficient=a.osmotic_coefficient if isinstance(a.osmotic_coefficient, float) else None))
    return out


def pixel_distance(dx, dy, curve_x, curve_y, calibration):
    px, py = calibration["px_per_x"], calibration["px_per_y"]
    return float(np.min(np.hypot((curve_x - dx) * px, (curve_y - dy) * py)))


def activity_coefficients():
    calibration = json.loads((HERE / "figures/digitized-calibration.json").read_text())["panels"]
    curves, verification, points = [], [], []
    cache = {}
    common = ("advanced", "revised", "advanced-dh-sigma", "revised-dh-sigma", "born-diameter-sigma")
    sets = {("LiBr", "methanol"): common + ("born-constant-permittivity", "dh-mixed-permittivity", "electrostatics-off")}
    for panel, (salt, solvent, x_max) in PANELS.items():
        for variant in sets.get((salt, solvent), common):
            if (variant, salt, solvent) not in cache:
                cache[variant, salt, solvent] = curve(variant, salt, solvent, x_max)
                curves += cache[variant, salt, solvent]
        digitized = read(f"figures/{panel}.csv")
        cal = calibration[panel]
        discs = [(float(r["pixel_x"]), float(r["pixel_y"])) for r in digitized if r["series"] == "points"]
        for series in sorted({r["series"] for r in digitized} - {"points"}):
            for engine_variant in (series, f"{series}-dh-sigma"):
                if (engine_variant, salt, solvent) not in cache:
                    continue
                c = cache[engine_variant, salt, solvent]
                for basis in ("gamma_pm_x", "gamma_pm_x_times_factor"):
                    cx, cy = np.array([r["x_salt"] for r in c]), np.array([r[basis] for r in c])
                    dist = [(pixel_distance(float(r["x"]), float(r["y"]), cx, cy, cal), r) for r in digitized if r["series"] == series]
                    worst = max(dist, key=lambda d: d[0])
                    dy = [abs(float(r["y"]) - float(np.interp(float(r["x"]), cx, cy))) for _, r in dist]
                    # Screened: past the steep dilute drop (5% of the axis) and > 10 px from every detected data disc.
                    screened = [d for d, r in dist if float(r["x"]) > 0.05 * x_max and min(
                        math.hypot(float(r["pixel_x"]) - qx, float(r["pixel_y"]) - qy) for qx, qy in discs) > 10]
                    verification.append(dict(panel=panel, series=series, engine_variant=engine_variant, engine_basis=basis,
                                             pixels=len(dist), median_px=float(np.median([d for d, _ in dist])), max_px=worst[0],
                                             max_px_screened=max(screened), screened_pixels=len(screened),
                                             worst_x=float(worst[1]["x"]), max_abs_gamma_at_same_x=max(dy),
                                             within_tolerance_fraction=sum(d <= TOL["curve pixel distance"] for d, _ in dist) / len(dist),
                                             tolerance_px=TOL["curve pixel distance"]))
        green = sorted((float(r["x"]), float(r["y"])) for r in digitized if r["series"] == "advanced")
        for r in (r for r in digitized if r["series"] == "points" and float(r["x"]) > 0.002 * x_max):
            x = float(r["x"])
            m = x_to_m(x, solvent)
            row = dict(panel=panel, salt=salt, solvent=solvent, x_salt=x, molality_mol_kg=m, gamma_pm_x_digitized=float(r["y"]),
                       paper_advanced_gamma_x=float(np.interp(x, *zip(*green))), note=r["note"])
            for variant in common[:4]:
                row[f"engine_{variant}"] = activity(variant, salt, solvent, m).mean_ionic_mole_fraction_scale
            row["engine_advanced_times_factor"] = row["engine_advanced"] * plotted_factor(m, solvent)
            points.append(row)
    # Fig. 6 right: Born contribution to ln gamma*, the difference between models with and without Born.
    ln = {v: np.log([r["gamma_pm_x"] for r in cache[v, "LiBr", "methanol"]])
          for v in ("advanced", "dh-mixed-permittivity", "electrostatics-off")}
    cx = np.array([r["x_salt"] for r in cache["advanced", "LiBr", "methanol"]])
    contributions = {"born-contribution": ln["advanced"] - ln["dh-mixed-permittivity"],
                     "debye-huckel-contribution": ln["dh-mixed-permittivity"] - ln["electrostatics-off"],
                     "advanced": ln["advanced"]}
    write("figure-06-contributions.csv", [dict(x_salt=x, ln_gamma_x=a, born=b, debye_huckel=d, hard_chain_dispersion_association=s)
                                          for x, a, b, d, s in zip(cx, ln["advanced"], contributions["born-contribution"],
                                                                   contributions["debye-huckel-contribution"], ln["electrostatics-off"])])
    right = read("figures/figure-06/libr-methanol-right.csv")
    cal = calibration["figure-06/libr-methanol-right"]
    for series, cy in contributions.items():
        dist = [(pixel_distance(float(r["x"]), float(r["y"]), cx, cy, cal), r) for r in right if r["series"] == series]
        verification.append(dict(panel="figure-06/libr-methanol-right", series=series + " (ln)", engine_variant="advanced",
                                 engine_basis="ln gamma_pm_x", pixels=len(dist),
                                 max_px_screened=max(d for d, r in dist if float(r["x"]) > 0.01),
                                 screened_pixels=sum(float(r["x"]) > 0.01 for _, r in dist),
                                 median_px=float(np.median([d for d, _ in dist])), max_px=max(d for d, _ in dist),
                                 worst_x=float(max(dist, key=lambda d: d[0])[1]["x"]),
                                 max_abs_gamma_at_same_x=max(abs(float(r["y"]) - float(np.interp(float(r["x"]), cx, cy))) for _, r in dist),
                                 within_tolerance_fraction=sum(d <= TOL["curve pixel distance"] for d, _ in dist) / len(dist),
                                 tolerance_px=TOL["curve pixel distance"]))
    write("miac-curves.csv", curves)
    write("miac-curve-verification.csv", verification)
    write("miac-digitized-points.csv", points)
    tabulated = []
    for solvent, path in (("ethanol", "figures/figure-07/source/ethanol-nabr.csv"), ("methanol", "figures/figure-07/source/methanol-nabr.csv")):
        for r in read(path):
            m = float(r["molality_mol_kg"] if "molality_mol_kg" in r else r["molality"])
            row = dict(salt="NaBr", solvent=solvent, molality_mol_kg=m, gamma_pm_m_han_pan=float(r["miac_m"]))
            for variant in common[:4]:
                row[f"engine_{variant}"] = activity(variant, "NaBr", solvent, m).mean_ionic_molality_scale
            tabulated.append(row)
    write("miac-han-pan.csv", tabulated)
    osmotic = []
    for r in read("inputs/source/zafarani-moattar-1998-osmotic.csv"):
        m = float(r["molality_mol_kg"])
        row = dict(salt=r["salt"], solvent="methanol", molality_mol_kg=m, osmotic_measured=float(r["osmotic_coefficient"]))
        for variant in common[:4]:
            phi = activity(variant, r["salt"], "methanol", m).osmotic_coefficient
            row[f"engine_{variant}"] = phi if isinstance(phi, float) else None
        osmotic.append(row)
    write("osmotic-zafarani-moattar.csv", osmotic)


def main():
    start = time.perf_counter()
    checks = infinite_dilution()
    write("checks.csv", checks)
    activity_coefficients()
    wall = time.perf_counter() - start
    inputs = sorted(p for p in (HERE / "inputs").rglob("*") if p.is_file()) + sorted((HERE / "figures").rglob("*.csv"))
    (HERE / "results/run-identity.json").write_text(json.dumps({
        "engine": list(epcsaft.provenance()), "calculation_wall_seconds": wall, "temperature_K": T, "pressure_Pa": PRESSURE,
        "scripts_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((HERE / "scripts").glob("*.py"))},
        "inputs_sha256": {p.relative_to(HERE).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        "conditions": "298.15 K, 1e5 Pa; salt-free solvent molality; gamma_pm_x on the ion mole-fraction scale",
        "tolerances": TOL}, indent=1) + "\n")
    print(f"done in {wall:.1f} s")


def declarations():
    return [
        mixture(v, (s, ION[c], ION[a]))
        for v in (*P.VARIANTS, 'electrostatics-off')
        for s in SOLVENT_M
        for c, a in {*SALTS.values(), *((i, p) if i.endswith('+') else (p, i) for i, p in PARTNER.items())}
    ]


if __name__ == "__main__":
    main()

"""Figiel-Yu-Held 2025 1:1 chloride baselines and design probes from the installed Engine.

Design evidence for issue #209; no coefficient is fitted here.
"""
from pathlib import Path
import copy
import csv
import hashlib
import json
import math
import time

from epcsaft import Mixture, Parameters, provenance
from scipy.constants import Avogadro, Boltzmann, elementary_charge, epsilon_0

HERE = Path(__file__).resolve().parents[1]
OUTPUT = HERE / "results/current-main"
REPO = HERE.parents[1]
T, P, MW = 298.15, 100000.0, 0.01801528
PACKET = REPO / "data/packets/figiel-yu-held-2025/4/parameters/parameters.json"
HELD2014 = REPO / "data/packets/held-2014-strategy2/1/parameters/parameters.json"
EVALUATED = [REPO / "analyses/2026-held-2014-chlorides/inputs/evaluated.csv", HERE / "inputs/kcl-hamer-wu.csv"]
CATIONS = {"LiCl": "lithium-cation", "NaCl": "sodium-cation", "KCl": "potassium-cation",
           "MgCl2": "magnesium-cation", "CaCl2": "calcium-cation"}


def write(name, rows):
    with (OUTPUT / name).open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def reference(document, cation, nu):
    model = Mixture(Parameters.from_mapping(document, components=["water", cation, "chloride-anion"]))
    return model, model.reference(T, P, solvent_fractions=[1, 0, 0], ions=(cation, "chloride-anion"),
                                  stoichiometry=(1, nu))


def coefficient(document, identity):
    for owner in document["components"] + document["pairs"]:
        for c in owner["coefficients"]:
            if c["identity"] == identity:
                return c["value"]


def divalent_document(figiel, held, cation, born_diameter):
    """Figiel basis plus one Held 2014 cation (sigma, u, segment); k_ij = 0 start, DH = packing = 0.88 sigma."""
    document = copy.deepcopy(figiel)
    source = next(c for c in held["components"] if c["component_id"] == cation)
    values = {c["family"]: c["value"]["magnitude"] for c in source["coefficients"]}
    sigma = values["segment_diameter"]
    ion = copy.deepcopy(next(c for c in document["components"] if c["component_id"] == "sodium-cation"))
    ion.update(component_id=cation, name=cation.replace("-", " "), formula=cation, aliases=[])
    ion["fixed"]["charge_number"]["value"]["magnitude"] = 2
    for c in ion["coefficients"]:
        c["identity"] = c["identity"].replace("sodium-cation", cation)
        c["provenance"] = {"source_id": "figiel-yu-held-2025", "locator": "issue #209 design probe", "domain_id": "figiel-298-k-1-bar"}
    new = {"segment_diameter": sigma, "debye_huckel_diameter": .88 * sigma, "packing_diameter": .88 * sigma,
           "dispersion_energy_over_k": values["dispersion_energy_over_k"], "born_diameter": born_diameter}
    for c in ion["coefficients"]:
        if c["family"] in new:
            c["value"]["magnitude"] = new[c["family"]]
    document["components"].append(ion)
    for other in ("water", "chloride-anion"):
        pair = copy.deepcopy(next(p for p in document["pairs"] if {p["component_id_a"], p["component_id_b"]} == {"sodium-cation", other}))
        pair["component_id_a"], pair["component_id_b"] = cation, other
        for c in pair["coefficients"]:
            c["identity"] = f"pair/{cation}/{other}/k_ij"
            c["value"]["magnitude"] = 0.0
        document["pairs"].append(pair)
    return document


def baselines(figiel, evaluated):
    rows, limits = [], []
    for salt in ("LiCl", "NaCl", "KCl"):
        model, ref = reference(figiel, CATIONS[salt], 1)
        # Salt-free water relative permittivity (Figiel Table 1) for the Debye limiting law.
        dielectric = next(c["constant"]["value"]["magnitude"] for c in figiel["correlations"]
                          if c["correlation_id"] == "component/water/relative_permittivity/constant")
        pure = model.state(T, P=P, x=[1, 0, 0], phase="liquid")
        for observed in (r for r in evaluated if r["salt"] == salt):
            m = float(observed["molality_mol_kg"])
            a = ref.activity(m)
            gm, ph = float(observed["gamma_evaluated"]), float(observed["osmotic_evaluated"])
            eg, ep = math.log10(a.mean_ionic_molality_scale / gm), math.log10(a.osmotic_coefficient / ph)
            rows.append(dict(salt=salt, molality_mol_kg=m, gamma_engine=a.mean_ionic_molality_scale, gamma_evaluated=gm,
                             log10_gamma_ratio=eg, osmotic_engine=a.osmotic_coefficient, osmotic_evaluated=ph,
                             log10_osmotic_ratio=ep, gamma_pass=abs(eg) <= .02, osmotic_pass=abs(ep) <= .02,
                             source_key=observed["source_key"], source_table=observed["source_table"]))
            if m == .001:
                # Debye limiting law at infinite dilution.
                bjerrum = elementary_charge**2 / (4 * math.pi * epsilon_0 * dielectric * Boltzmann * T)
                a_ln = math.sqrt(2 * math.pi * pure.molar_density * MW * Avogadro) * bjerrum**1.5
                d = a_ln * math.sqrt(m)
                limits.append(dict(salt=salt, molality_mol_kg=m, gamma_engine=a.mean_ionic_molality_scale,
                                   gamma_debye_limiting=math.exp(-d), osmotic_engine=a.osmotic_coefficient,
                                   osmotic_debye_limiting=1 - d / 3))
    write("figiel-1to1-comparison.csv", rows)
    write("figiel-1to1-debye-limit.csv", limits)
    summary = []
    for salt in ("LiCl", "NaCl", "KCl"):
        data = [r for r in rows if r["salt"] == salt]
        g = max(data, key=lambda r: abs(r["log10_gamma_ratio"]))
        f = max(data, key=lambda r: abs(r["log10_osmotic_ratio"]))
        summary.append(dict(salt=salt, n_rows=len(data), max_abs_gamma_decade=abs(g["log10_gamma_ratio"]),
                            gamma_worst_molality=g["molality_mol_kg"], max_abs_osmotic_decade=abs(f["log10_osmotic_ratio"]),
                            osmotic_worst_molality=f["molality_mol_kg"],
                            gamma_fail_rows=sum(not r["gamma_pass"] for r in data),
                            osmotic_fail_rows=sum(not r["osmotic_pass"] for r in data)))
    write("figiel-1to1-maximum-errors.csv", summary)
    return summary


def divalent_probes(figiel, held):
    """Born-diameter sensitivity and finite-difference step agreement at the Held 2014 start, k_ij = 0."""
    born, steps = [], []
    for salt, nu in (("MgCl2", 2), ("CaCl2", 2)):
        cation = CATIONS[salt]
        for d_born in (3.0, 3.5, 4.0):
            _, ref = reference(divalent_document(figiel, held, cation, d_born), cation, 2)
            for m in (.001, .1, 1.0, 4.5):
                a = ref.activity(m)
                born.append(dict(salt=salt, born_diameter_angstrom=d_born, molality_mol_kg=m,
                                 log10_gamma=math.log10(a.mean_ionic_molality_scale),
                                 log10_osmotic=math.log10(a.osmotic_coefficient)))
        base = divalent_document(figiel, held, cation, 3.5)
        for identity in (f"pair/{cation}/water/k_ij", f"pair/{cation}/chloride-anion/k_ij", f"component/{cation}/segment_diameter"):
            for h in (1e-3, 1e-4):
                values = {}
                for sign in (1, -1):
                    document = copy.deepcopy(base)
                    v = coefficient(document, identity)
                    v["magnitude"] += sign * h
                    if identity.endswith("segment_diameter"):
                        for fam in ("debye_huckel_diameter", "packing_diameter"):
                            coefficient(document, identity.replace("segment_diameter", fam))["magnitude"] = .88 * v["magnitude"]
                    _, ref = reference(document, cation, 2)
                    values[sign] = [ref.activity(m) for m in (.01, 1.0, 4.5)]
                for i, m in enumerate((.01, 1.0, 4.5)):
                    up, down = values[1][i], values[-1][i]
                    steps.append(dict(salt=salt, coordinate=identity, step=h, molality_mol_kg=m,
                                      dlog10gamma=(math.log10(up.mean_ionic_molality_scale) - math.log10(down.mean_ionic_molality_scale)) / (2 * h),
                                      dlog10phi=(math.log10(up.osmotic_coefficient) - math.log10(down.osmotic_coefficient)) / (2 * h)))
    write("born-diameter-sensitivity.csv", born)
    write("finite-difference-steps.csv", steps)


def born_form_checks(figiel):
    """Which Born form reproduces Figiel's own results: Table S5 modelled x-process solvation
    Gibbs energies (RT ln phi_i^inf) and the digitized Figure 5 LiBr model curve."""
    s5 = {"lithium-cation": -487.0, "sodium-cation": -381.9, "potassium-cation": -308.9,
          "chloride-anion": -315.0, "bromide-anion": -290.7}  # Figiel 2025 SI Table S5, modelled column, kJ/mol
    curve = list(csv.DictReader((REPO / "analyses/reproduction/2025-figiel/figures/figure-05/source/libr-model-digitized.csv").open()))
    rows = []
    for c_shell, c_dielectric, form in ((1, 1, "corrected Eq. 7 (Engine packet)"), (1, 0, "shell only"),
                                        (0, 0, "direct Born at d_Born (printed Eq. 7 reduces to this)")):
        document = copy.deepcopy(figiel)
        family = next(f for f in document["model_families"] if f["kind"] == "electrolyte")
        family["c_shell"], family["c_dielectric"] = c_shell, c_dielectric
        energies = {}
        for cation, anion in (("lithium-cation", "bromide-anion"), ("sodium-cation", "chloride-anion"), ("potassium-cation", "chloride-anion")):
            model = Mixture(Parameters.from_mapping(document, components=["water", cation, anion]))
            ref = model.reference(T, P, solvent_fractions=[1, 0, 0], ions=(cation, anion), stoichiometry=(1, 1))
            energies.update(zip((cation, anion), (8.314462618e-3 * T * v for v in ref.log_fugacity)))
            if cation == "lithium-cation":
                points = []
                for p in curve:
                    m, g, u = float(p["molality_mol_kg"]), float(p["gamma_pm_m"]), float(p["digitization_uncertainty"])
                    engine = ref.activity(m).mean_ionic_molality_scale
                    points.append(dict(molality_mol_kg=m, gamma_engine=engine, gamma_digitized=g, digitization_uncertainty=u,
                                       log10_ratio=math.log10(engine / g), difference_in_uncertainty_units=(engine - g) / u))
                curve_error = max(abs(r["log10_ratio"]) for r in points)
                if c_shell == c_dielectric == 1:
                    write("libr-curve-comparison.csv", points)
        rows.append(dict(form=form, c_shell=c_shell, c_dielectric=c_dielectric,
                         max_abs_solvation_error_kj_mol=max(abs(energies[i] - s5[i]) for i in s5),
                         libr_curve_max_abs_decade=curve_error, libr_curve_points=len(curve),
                         **{f"dG_{i}_kj_mol": energies[i] for i in s5}))
    write("born-form-discrimination.csv", rows)


def run():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    figiel, held = json.loads(PACKET.read_text()), json.loads(HELD2014.read_text())
    evaluated = [r for path in EVALUATED for r in csv.DictReader(path.open())]
    start = time.perf_counter()
    for row in baselines(figiel, evaluated):
        print(row)
    elapsed_baseline = time.perf_counter() - start
    divalent_probes(figiel, held)
    born_form_checks(figiel)
    (OUTPUT / "run-identity.json").write_text(json.dumps({
        "engine_provenance": list(provenance()),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "packet_sha256": hashlib.sha256(PACKET.read_bytes()).hexdigest(),
        "baseline_seconds_78_rows": round(elapsed_baseline, 2),
        "conditions": "298.15 K, 100000 Pa, fully dissociated formula-unit molality/kg water",
    }, indent=2) + "\n")


if __name__ == "__main__":
    run()

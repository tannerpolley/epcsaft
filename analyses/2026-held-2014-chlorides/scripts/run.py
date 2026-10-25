"""Qualify binary chlorides with native properties from the installed Engine."""
from pathlib import Path
import argparse
import copy
import csv
import hashlib
import json
import math

import epcsaft
from epcsaft import Mixture, Parameters
from scipy.constants import Avogadro, Boltzmann, elementary_charge, epsilon_0

HERE = Path(__file__).resolve().parents[1]
OUTPUT = HERE / "results/sigma-replay"
T, P, MW = 298.15, 100000.0, 0.01801528
SALTS = {"LiCl": ("lithium-cation", 1), "NaCl": ("sodium-cation", 1),
         "MgCl2": ("magnesium-cation", 2), "CaCl2": ("calcium-cation", 2),
         "SrCl2": ("strontium-cation", 2)}


def write(name, rows):
    with (OUTPUT / name).open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def check_energy_coordinates(document):
    """Binary cross-energy invariance; this is a probe, not a parameter fit."""
    rows = []
    for salt in ("MgCl2", "CaCl2"):
        cation, nu_minus = SALTS[salt]
        altered = copy.deepcopy(document)
        component = next(c for c in altered["components"] if c["component_id"] == cation)
        energy = next(c for c in component["coefficients"] if c["family"] == "dispersion_energy_over_k")
        energy["value"]["magnitude"] *= 1.21
        for pair in altered["pairs"]:
            if cation in (pair["component_id_a"], pair["component_id_b"]):
                coefficient = pair["coefficients"][0]["value"]
                coefficient["magnitude"] = 1 - (1 - coefficient["magnitude"]) / 1.1
        models = [Mixture(Parameters.from_mapping(d, components=["water", cation, "chloride-anion"]))
                  for d in (document, altered)]
        refs = [model.reference(T, P, solvent_fractions=[1, 0, 0],
                                ions=(cation, "chloride-anion"), stoichiometry=(1, nu_minus)) for model in models]
        for m in (.001, .4, 1.5, 2.25, 4.5):
            original, changed = [ref.activity(m) for ref in refs]
            eg = abs(math.log10(changed.mean_ionic_molality_scale / original.mean_ionic_molality_scale))
            ep = abs(changed.osmotic_coefficient - original.osmotic_coefficient)
            rows.append(dict(salt=salt, molality_mol_kg=m, ion_energy_factor=1.21,
                             cross_energy_preserving_kij_map="1-(1-kij)/sqrt(1.21)",
                             original_gamma=original.mean_ionic_molality_scale,
                             changed_gamma=changed.mean_ionic_molality_scale,
                             original_phi=original.osmotic_coefficient, changed_phi=changed.osmotic_coefficient,
                             gamma_difference_decade=eg, phi_absolute_difference=ep))
            assert eg < 1e-8 and ep < 1e-6  # Far below the 0.02-decade qualification target.
    write("energy-coordinate-invariance.csv", rows)


PACKETS = {"packet4": "data/packets/held-cameretti-sadowski-2008/4/parameters/parameters.json",
           "held2014": "data/packets/held-2014-strategy2/1/parameters/parameters.json"}


def run(dh_088=False):
    global OUTPUT
    OUTPUT = HERE / "results/current-main" / ("dh-088" if dh_088 else "sigma-replay")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    documents = {key: json.loads((HERE.parents[1] / path).read_text()) for key, path in PACKETS.items()}
    # Bounded source ambiguity check, not a new packet or a fitted coefficient.
    # Floriano-Nascimento Eq.(3), Table I at 298 K and 0.1 MPa (not 298.15 K).
    epsilon_nearest = 78.85 + 14.1113 * math.log((341.5902 + .1) / (341.5902 + 10))
    sensitivity = copy.deepcopy(documents["held2014"])
    sensitivity["correlations"][0]["constant"]["value"]["magnitude"] = epsilon_nearest
    documents["held2014_nearest_dielectric"] = sensitivity
    if dh_088:
        variant = copy.deepcopy(documents["held2014"])
        for component in variant["components"]:
            coefficients = {c["family"]: c for c in component["coefficients"]}
            if "debye_huckel_diameter" in coefficients:
                coefficients["debye_huckel_diameter"]["value"]["magnitude"] = (
                    0.88 * coefficients["segment_diameter"]["value"]["magnitude"])
                coefficients["debye_huckel_diameter"]["provenance"]["locator"] = (
                    "Inferred 0.88 sigma convention: issue #387, PRs #395, #397, #399; issue #398")
        documents = {"held2014_dh088": variant}
    else:
        check_energy_coordinates(documents["held2014"])
    evaluated = list(csv.DictReader((HERE / "inputs/evaluated.csv").open()))
    curves, comparisons, limits = [], [], []
    for mode, document in documents.items():
        for salt, (cation, nu_minus) in SALTS.items():
            if salt == "SrCl2" and mode != "packet4":
                continue  # Held 2014 Tables 2-3 give no Sr parameters.
            if mode.endswith("dielectric") and salt not in ("MgCl2", "CaCl2"):
                continue
            model = Mixture(Parameters.from_mapping(document, components=["water", cation, "chloride-anion"]))
            reference = model.reference(T, P, solvent_fractions=[1, 0, 0],
                                        ions=(cation, "chloride-anion"), stoichiometry=(1, nu_minus))
            selected = [r for r in evaluated if r["salt"] == salt]
            observed_by_m = {float(r["molality_mol_kg"]): r for r in selected}
            molalities = set(observed_by_m)
            if not mode.endswith("dielectric"):
                molalities |= {0.001 * 4500 ** (i / 80) for i in range(81)}
            pure = model.state(T, P=P, x=[1, 0, 0], phase="liquid")
            for m in sorted(molalities):
                activity = reference.activity(m)
                nw = 1 / MW
                total = nw + (1 + nu_minus) * m
                x = [nw / total, m / total, nu_minus * m / total]
                state = model.state(T, P=P, x=x, phase="liquid")
                logs = state.log_fugacity_coefficient
                # Independent contraction checks the molal conversion and stoichiometry.
                log_gamma = sum(n * (logs[i] - reference.log_fugacity[i - 1])
                                for i, n in [(1, 1), (2, nu_minus)]) / (1 + nu_minus) + math.log(x[0])
                # Native Phi uses prescribed P; remove density-closure residuals from the check.
                pressure_correction = ((state.pressure-P)/state.molar_density
                                       - (pure.pressure-P)/pure.molar_density) / (epcsaft.GAS_CONSTANT_J_PER_MOL_K*T)
                phi_formula = -math.fsum((math.log(x[0]), logs[0],
                                          -reference.pure_solvent_log_fugacity, -pressure_correction)) / ((1 + nu_minus) * m * MW)
                gamma, phi = activity.mean_ionic_molality_scale, activity.osmotic_coefficient
                row = dict(mode=mode, salt=salt, temperature_K=T, pressure_Pa=P, molality_mol_kg=m,
                           gamma_pm_m=gamma, osmotic_coefficient=phi,
                           molar_density_mol_m3=state.molar_density, pressure_residual_Pa=state.pressure-P,
                           pressure_slope=state.pressure_density_derivative, charge_residual=nu_minus*x[1]-x[2],
                           reference_molality=reference.molality, reference_terminal_change_log_phi=reference.convergence_error,
                           gamma_formula_error=abs(math.log(gamma)-log_gamma),
                           osmotic_formula_error=abs(phi-phi_formula))
                assert abs(state.pressure-P) < .01 and state.pressure_density_derivative > 0
                assert row["gamma_formula_error"] < 1e-12 and row["osmotic_formula_error"] < 1e-10
                assert gamma > 0 and math.isfinite(phi)
                curves.append(row)
                if m in observed_by_m:
                    observed = observed_by_m[m]
                    gm, ph = float(observed["gamma_evaluated"]), float(observed["osmotic_evaluated"])
                    eg = math.log10(gamma / gm)
                    ep = math.log10(phi / ph) if phi > 0 else None
                    comparisons.append(dict(mode=mode, salt=salt, molality_mol_kg=m,
                                            gamma_engine=gamma, gamma_evaluated=gm, log10_gamma_ratio=eg,
                                            osmotic_engine=phi, osmotic_evaluated=ph, log10_osmotic_ratio=ep,
                                            gamma_pass=abs(eg) <= .02, osmotic_pass=ep is not None and abs(ep) <= .02,
                                            source_key=observed["source_key"], attachment_key=observed["attachment_key"],
                                            source_table=observed["source_table"], source_row=observed["source_row"],
                                            saturation_status=observed["saturation_status"]))
                if m == .001:
                    # SI Debye limiting law with pure-water model density, formula-unit molality.
                    dielectric = document["correlations"][0]["constant"]["value"]["magnitude"]
                    bjerrum = elementary_charge**2 / (4 * math.pi * epsilon_0 * dielectric * Boltzmann * T)
                    a_ln = math.sqrt(2 * math.pi * pure.molar_density * MW * Avogadro) * bjerrum**1.5
                    ionic_strength = .5 * (nu_minus**2 + nu_minus) * m
                    decrement = nu_minus * a_ln * math.sqrt(ionic_strength)
                    assert 0 < gamma < 1 and 0 < phi < 1
                    limits.append(dict(mode=mode, salt=salt, molality_mol_kg=m, ionic_strength_mol_kg=ionic_strength,
                                       relative_permittivity=dielectric, pure_water_density_kg_m3=pure.molar_density*MW,
                                       debye_A_natural_log=a_ln, gamma_engine=gamma,
                                       gamma_debye_limiting=math.exp(-decrement), osmotic_engine=phi,
                                       osmotic_debye_limiting=1-decrement/3,
                                       gamma_finite_size_remainder_log=math.log(gamma)+decrement,
                                       osmotic_finite_size_remainder=phi-(1-decrement/3)))
            print(mode, salt, len(molalities), "states", flush=True)
    write("model-curves.csv", curves)
    write("evaluated-comparison.csv", comparisons)
    write("debye-limit.csv", limits)
    summaries = []
    for mode in documents:
        for salt in SALTS:
            data = [r for r in comparisons if r["mode"] == mode and r["salt"] == salt]
            if not data:
                continue
            g = max(data, key=lambda r: abs(r["log10_gamma_ratio"]))
            f = max(data, key=lambda r: abs(r["log10_osmotic_ratio"]))
            summaries.append(dict(mode=mode, salt=salt, n_rows=len(data),
                                  max_abs_gamma_decade=abs(g["log10_gamma_ratio"]), gamma_worst_molality=g["molality_mol_kg"],
                                  max_abs_osmotic_decade=abs(f["log10_osmotic_ratio"]), osmotic_worst_molality=f["molality_mol_kg"],
                                  gamma_fail_rows=sum(not r["gamma_pass"] for r in data),
                                  osmotic_fail_rows=sum(not r["osmotic_pass"] for r in data),
                                  gamma_ard_percent=100*sum(abs(r["gamma_engine"]/r["gamma_evaluated"]-1) for r in data)/len(data),
                                  osmotic_ard_percent=100*sum(abs(r["osmotic_engine"]/r["osmotic_evaluated"]-1) for r in data)/len(data),
                                  passes=all(r["gamma_pass"] and r["osmotic_pass"] for r in data)))
    if not dh_088:
        # #594 returns the extrapolated infinite-dilution reference, not its finite rung.
        expected = {"MgCl2": ("0.085605", 2.25, "0.045990", 1.5),
                    "CaCl2": ("0.089351", 2.25, "0.045722", 1.5)}
        for row in summaries:
            if row["mode"] != "held2014":
                continue
            if row["salt"] in expected:
                actual = (f'{row["max_abs_gamma_decade"]:.6f}', row["gamma_worst_molality"],
                          f'{row["max_abs_osmotic_decade"]:.6f}', row["osmotic_worst_molality"])
                assert actual == expected[row["salt"]], (actual, expected[row["salt"]])
            else:
                assert row["passes"], row
    write("maximum-errors.csv", summaries)
    repo = HERE.parents[1]
    input_files = sorted(p for p in (HERE / "inputs").rglob("*") if p.is_file()) + [repo / path for path in PACKETS.values()]
    (OUTPUT / "run-identity.json").write_text(json.dumps({
        "engine_provenance": list(epcsaft.provenance()),
        "dh_factor": 0.88 if dh_088 else 1.0,
        "packing_factor": 0.88,
        "convention": "inferred: #387 and PRs #395, #397, #399" if dh_088 else "sigma record",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": {(p.relative_to(HERE) if p.is_relative_to(HERE) else p.relative_to(repo)).as_posix():
                         hashlib.sha256(p.read_bytes()).hexdigest() for p in input_files},
        "conditions": "298.15 K, 100000 Pa, fully dissociated formula-unit molality/kg water",
        "limits": "Held2014 dielectric 78.09 is a declared supplement; nearest-temperature Floriano sensitivity is not an exact author replay. No refit or ternary validation."
    }, indent=2)+"\n")
    for row in summaries:
        print(row)


def declarations():
    return [
        Mixture(
            Parameters.from_json(
                HERE.parents[1] / f'data/packets/{folder}/parameters/parameters.json',
                components=('water', ion, 'chloride-anion'),
            )
        )
        for folder in ('held-cameretti-sadowski-2008/4', 'held-2014-strategy2/1')
        for ion, _ in SALTS.values()
        if not (folder == 'held-2014-strategy2/1' and ion == 'strontium-cation')
    ]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dh-088", action="store_true", help="Inferred DH 0.88 sigma; keep packing fixed")
    run(parser.parse_args().dh_088)

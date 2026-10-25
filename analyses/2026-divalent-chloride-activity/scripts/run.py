"""Held 2008 diameter diagnosis through the installed Engine wheel."""
from pathlib import Path
import copy
import csv
import hashlib
import json
import math

from epcsaft import Mixture, Parameters, provenance

HERE = Path(__file__).resolve().parents[1]
OUTPUT = HERE / "results/current-main"
T, P, MW = 298.15, 100000.0, 0.01801528
SALTS = {"LiCl": ("lithium-cation", 1), "NaCl": ("sodium-cation", 1),
         "MgCl2": ("magnesium-cation", 2), "CaCl2": ("calcium-cation", 2)}


def write(name, rows):
    with (OUTPUT / name).open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    before = json.loads((HERE.parents[1] / "data/packets/held-cameretti-sadowski-2008/3/parameters/parameters.json").read_text())
    after = json.loads((HERE.parents[1] / "data/packets/held-cameretti-sadowski-2008/4/parameters/parameters.json").read_text())
    documents = {"packet3": before, "packet4": after}
    # Rejected probes distinguish the two independent lengths and dispersion.
    for mode in ("packing_sigma", "both_sigma", "opposite_dispersion"):
        doc = copy.deepcopy(after if mode == "both_sigma" else before)
        if mode == "opposite_dispersion":
            for pair in doc["pairs"]:
                if "water" not in (pair["component_id_a"], pair["component_id_b"]):
                    pair["coefficients"][0]["value"]["magnitude"] = 0.0
        else:
            for component in doc["components"]:
                coefficients = {c["family"]: c for c in component["coefficients"]}
                if "packing_diameter" in coefficients:
                    coefficients["packing_diameter"]["value"]["magnitude"] = coefficients["segment_diameter"]["value"]["magnitude"]
        documents[mode] = doc
    evaluated = list(csv.DictReader((HERE / "inputs/evaluated.csv").open()))
    anchors = {"LiCl": 0.02, "NaCl": 3.0, "MgCl2": 0.4, "CaCl2": 1.5}
    comparisons, curves, probes, limits, conditioning = [], [], [], [], []
    for mode, doc in documents.items():
        for salt, (cation, nu_minus) in SALTS.items():
            model = Mixture(Parameters.from_mapping(doc, components=["water", cation, "chloride-anion"]))
            reference = model.reference(T, P, solvent_fractions=[1, 0, 0],
                                        ions=(cation, "chloride-anion"), stoichiometry=(1, nu_minus))
            selected = [r for r in evaluated if r["salt"] == salt]
            molalities = {anchors[salt]} if mode not in ("packet3", "packet4") else {
                0.001 * (4500 ** (i / 80)) for i in range(81)
            } | {float(r["molality_mol_kg"]) for r in selected}
            for m in sorted(molalities):
                activity = reference.activity(m)
                solvent_moles = 1 / MW
                total = solvent_moles + m * (1 + nu_minus)
                x = [solvent_moles / total, m / total, nu_minus * m / total]
                xw = x[0]
                state = model.state(T, P=P, x=x, phase="liquid")
                logs = state.log_fugacity_coefficient
                manual_log = sum(n * (logs[i] - reference.log_fugacity[i - 1])
                                 for i, n in [(1, 1), (2, nu_minus)]) / (1 + nu_minus) + math.log(xw)
                # #214/#556 use prescribed pressure; separately closed ln(phi) carries root residuals.
                manual_phi = reference.solvent_activity(x=x).osmotic_coefficient
                gamma = activity.mean_ionic_molality_scale
                row = dict(mode=mode, salt=salt, temperature_K=T, pressure_Pa=P, molality_mol_kg=m,
                           gamma_pm_m=gamma, gamma_pm_x=activity.mean_ionic_mole_fraction_scale,
                           osmotic_coefficient=activity.osmotic_coefficient, molar_density_mol_m3=state.molar_density,
                           pressure_residual_Pa=state.pressure-P,
                           pressure_slope=state.pressure_density_derivative,
                           charge_residual=nu_minus*x[1]-x[2],
                           reference_molality=reference.molality,
                           reference_terminal_change_log_phi=reference.convergence_error,
                           gamma_formula_error=abs(math.log(gamma)-manual_log),
                           osmotic_formula_error=abs(activity.osmotic_coefficient-manual_phi))
                assert abs(state.pressure-P) < 0.01 and state.pressure_density_derivative > 0
                assert row["gamma_formula_error"] < 1e-12 and row["osmotic_formula_error"] < 1e-10, (mode, salt, m, row["gamma_formula_error"], row["osmotic_formula_error"])
                if mode not in ("packet3", "packet4"):
                    probes.append(row)
                else:
                    curves.append(row)
                    for observed in selected:
                        if float(observed["molality_mol_kg"]) == m:
                            gm = float(observed["gamma_evaluated"])
                            phi = observed["osmotic_evaluated"]
                            comparisons.append(dict(mode=mode, salt=salt, molality_mol_kg=m,
                                                    gamma_engine=gamma, gamma_evaluated=gm,
                                                    log10_gamma_ratio=math.log10(gamma/gm),
                                                    osmotic_engine=activity.osmotic_coefficient,
                                                    osmotic_evaluated=phi,
                                                    source=observed["source"], source_verified=observed["source_verified"]))
            if mode in ("packet3", "packet4"):
                pure = model.state(T, P=P, x=[1, 0, 0], phase="liquid")
                # Same-pressure chemical potentials remove log(P_state/P_pure)
                # amplification by 1/m from independently terminated density roots.
                for m in (1e-3, 1e-6, 1e-8, 1e-10):
                    nw = 1 / MW
                    total = nw + (1 + nu_minus) * m
                    x = [nw / total, m / total, nu_minus * m / total]
                    finite = model.state(T, P=P, x=x, phase="liquid")
                    log_as = (math.log(x[0]) + finite.residual_chemical_potential_over_rt[0]
                              - pure.residual_chemical_potential_over_rt[0]
                              + math.log(finite.molar_density / pure.molar_density))
                    conditioning.append(dict(mode=mode, salt=salt, molality_mol_kg=m,
                                             osmotic_native=reference.activity(m).osmotic_coefficient,
                                             osmotic_equal_pressure_probe=-log_as / ((1 + nu_minus) * m * MW),
                                             finite_pressure_residual_Pa=finite.pressure-P,
                                             pure_pressure_residual_Pa=pure.pressure-P))
                for m in (1e-8, 1e-10, 1e-12, 1e-14):
                    a = reference.activity(m)
                    limits.append(dict(mode=mode, salt=salt, molality_mol_kg=m,
                                       gamma_pm_m=a.mean_ionic_molality_scale,
                                       osmotic_coefficient=a.osmotic_coefficient))
            print(mode, salt, "computed", len(molalities), flush=True)
    write("model-curves.csv", curves)
    write("evaluated-comparison.csv", comparisons)
    write("rejected-probes.csv", probes)
    write("dilution-limits.csv", limits)
    write("solvent-activity-conditioning.csv", conditioning)
    summary = []
    for mode in ("packet3", "packet4"):
        for salt in SALTS:
            selected = [r for r in comparisons if r["mode"] == mode and r["salt"] == salt]
            for max_m in (4.5, 3.2 if salt == "NaCl" else 4.5):
                if max_m == 4.5 and summary and summary[-1].get("mode") == mode and summary[-1].get("salt") == salt:
                    continue
                data = [r for r in selected if r["molality_mol_kg"] <= max_m]
                summary.append(dict(mode=mode, salt=salt, maximum_molality_mol_kg=max_m, n_points=len(data),
                                    gamma_aad=sum(abs(r["gamma_engine"]-r["gamma_evaluated"]) for r in data)/len(data),
                                    gamma_ard_percent=100*sum(abs(r["gamma_engine"]/r["gamma_evaluated"]-1) for r in data)/len(data),
                                    max_abs_decade=max(abs(r["log10_gamma_ratio"]) for r in data),
                                    basis="Hamer-Wu evaluated grid" if salt in ("NaCl", "LiCl") else "one unverified source anchor supplied in issue209"))
    write("aggregate-comparison.csv", summary)
    (OUTPUT / "run-identity.json").write_text(json.dumps({
        "engine_provenance": list(provenance()),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "conditions": "298.15 K, 100000 Pa, fully dissociated chloride in pure water",
        "limits": "No original Lobo fitting rows or verified Goldberg-Nuttall grid; no exact reproduction of Held Table4 AAD; no refit."
    }, indent=2)+"\n")


def declarations():
    return [
        Mixture(
            Parameters.from_json(
                HERE.parents[1] / f'data/packets/held-cameretti-sadowski-2008/{v}/parameters/parameters.json',
                components=('water', ion, 'chloride-anion'),
            )
        )
        for v in (3, 4)
        for ion, _ in SALTS.values()
    ]


if __name__ == "__main__":
    run()

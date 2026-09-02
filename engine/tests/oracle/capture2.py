"""Oracle: physical quantities from the installed engine, by chemistry.

No solver status, certificate, route label, or fingerprint is recorded. Cases
are named for the species and the phenomenon, never for an author, figure, or
downstream application. A case the current engine cannot evaluate is recorded
as unevaluable with its reason: that is a baseline to beat, not a failure.
"""
from __future__ import annotations
import csv, json, traceback
from pathlib import Path
import epcsaft
from epcsaft import equilibrium, unit_registry as u

REPO = Path("~/Workspaces/Engineering/ePC-SAFT-project")
PACKETS = REPO / "data/packets"
OUT = Path(__file__).resolve().parent

def pk(name, ver="2"):
    return json.loads((PACKETS/name/ver/"parameters/parameters.json").read_text())

def rows_of(name, pattern):
    f = next(iter(sorted((PACKETS/name).rglob(pattern))), None)
    return list(csv.DictReader(f.open())) if f else []

def r6(v): return None if v is None else round(float(v), 12)

# ---------- 1. pure saturation of light alkanes, on the measured grid ----------
def light_alkane_saturation():
    mapping, out = pk("gross-sadowski-2001"), []
    ref = rows_of("nist-srd69-pure-saturation-properties", "*.csv")
    for row in ref:
        species = (row.get("species") or "").strip().lower()
        if species not in {"methane", "ethane", "propane"}:
            continue
        t = float(row["T_K"])
        try:
            params = epcsaft.Parameters.from_mapping(mapping, components=[species])
            s = epcsaft.Pure(params).saturation(T=t*u.kelvin)
            out.append({"species": species, "temperature_k": r6(t),
                        "saturation_pressure_pa": r6(s.saturation_pressure_pa),
                        "vapor_molar_volume_m3_per_mol": r6(s.vapor.volume_m3/s.vapor.amount_mol),
                        "liquid_molar_volume_m3_per_mol": r6(s.liquid.volume_m3/s.liquid.amount_mol),
                        "measured_pressure_pa": r6(row.get("p_sat_Pa")),
                        "measured_liquid_density_kg_m3": r6(row.get("rho_sat_liq_kg_m3"))})
        except Exception as exc:
            out.append({"species": species, "temperature_k": r6(t),
                        "unevaluable": f"{type(exc).__name__}: {str(exc)[:90]}"})
    return out

# ---------- 2. mean ionic activity of aqueous alkali halides ----------
SALTS = {"LiCl": ("lithium-cation","chloride-anion"), "LiBr": ("lithium-cation","bromide-anion"),
         "KCl": ("potassium-cation","chloride-anion"), "KBr": ("potassium-cation","bromide-anion")}

def aqueous_alkali_halide_activity():
    mapping, out = pk("figiel-2025-reference-electrolytes"), []
    measured = rows_of("hamer-wu-1972-aqueous-alkali-halides", "*.csv")
    want = {0.001, 0.01, 0.1, 0.5, 1.0, 2.0}
    for salt, ions in SALTS.items():
        try:
            params = epcsaft.Parameters.from_mapping(mapping, components=["water", *ions])
            mix = epcsaft.Mixture(params)
        except Exception as exc:
            out.append({"salt": salt, "unevaluable": f"{type(exc).__name__}: {str(exc)[:90]}"}); continue
        for row in measured:
            if (row.get("salt") or "").strip() != salt: continue
            mol = float(row["molality_mol_kg"])
            if not any(abs(mol-w) < 1e-9 for w in want): continue
            try:
                ref = mix.reference(T=298.15*u.kelvin, P=1.0*u.bar, solvent_x=(1.0,0.0,0.0),
                                    phase="liquid", ion_component_ids=ions,
                                    formula_stoichiometry=(1,1),
                                    reference_convention="molality-infinite-dilution")
                lf = ref.infinite_dilution_log_fugacity_coefficient
                out.append({"salt": salt, "cation": ions[0], "anion": ions[1],
                            "temperature_k": 298.15, "pressure_pa": 100000.0,
                            "molality_mol_kg": r6(mol),
                            "infinite_dilution_log_fugacity_cation": r6(lf[0]),
                            "infinite_dilution_log_fugacity_anion": r6(lf[1]),
                            "pure_solvent_log_fugacity": r6(ref._pure_solvent_log_fugacity),
                            "convergence_error": r6(ref.convergence_error),
                            "measured_mean_ionic_activity_coefficient": r6(row.get("gamma_pm_m"))})
            except Exception as exc:
                out.append({"salt": salt, "molality_mol_kg": r6(mol),
                            "unevaluable": f"{type(exc).__name__}: {str(exc)[:90]}"})
    return out

# ---------- 3. reactive equilibrium: does today's engine reach it at all? ----------
def reactive_attempts():
    out = []
    for label, packet_name, comps in (
        ("aqueous_amine_carbon_dioxide_reactive", "mea-co2-h2o-nine-species-estimation", None),
        ("aqueous_lithium_chloride_reactive", "figiel-2025-reference-electrolytes",
         ["water", "lithium-cation", "chloride-anion"]),
    ):
        entry = {"system": label, "packet": packet_name}
        try:
            mapping = pk(packet_name)
            ids = [c["component_id"] for c in mapping.get("components", [])]
            entry["components_available"] = comps or ids
            params = epcsaft.Parameters.from_mapping(mapping, components=(comps or ids))
            epcsaft.Mixture(params)
            entry["model_constructed"] = True
            entry["reactive_problem_constructed"] = False
            entry["note"] = ("model builds; a reaction system and declared phases are required "
                             "to solve, which this capture does not assemble")
        except Exception as exc:
            entry["model_constructed"] = False
            entry["unevaluable"] = f"{type(exc).__name__}: {str(exc)[:140]}"
        out.append(entry)
    return out

def main():
    cases = {"light_alkane_saturation": light_alkane_saturation(),
             "aqueous_alkali_halide_activity": aqueous_alkali_halide_activity(),
             "reactive_equilibrium_reachability": reactive_attempts()}
    (OUT/"oracle.json").write_text(json.dumps({"cases": cases}, indent=2, sort_keys=True)+"\n")
    for k, v in cases.items():
        bad = sum(1 for r in v if "unevaluable" in r)
        print(f"  {k:38s} {len(v):3d} rows ({bad} unevaluable)")

if __name__ == "__main__":
    main()

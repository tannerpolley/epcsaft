"""Stage-A joint fit (#575): declare the frozen training rows and starts; the native Ceres fit owns every calculation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import subprocess
import sys
import time
from pathlib import Path

from epcsaft import Parameters, regression as r

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
DECLARATION = json.loads((HERE / "inputs/fit-declaration.json").read_text())
T, P = DECLARATION["conditions"]["temperature_K"], DECLARATION["conditions"]["pressure_Pa"]
COMPONENTS = tuple(DECLARATION["conditions"]["components"])
CATIONS = ("sodium-cation", "magnesium-cation", "calcium-cation")
SALTS = {"NaCl": ("sodium-cation", (1, 1)), "MgCl2": ("magnesium-cation", (1, 2)), "CaCl2": ("calcium-cation", (1, 2))}
BOUNDS = DECLARATION["bounds"]
RESULTS = HERE / "results" / DECLARATION["run"]
SAAD = ROOT / "data/packets/saad-padova-marcus-1980-nacl-kcl-cacl2-observations/1/observations/saad-padova-marcus-1980-tables-i-iii.csv"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def comparison():
    """The frozen #556 comparison code, loaded by path and checked against its declared hash."""
    path = ROOT / DECLARATION["held_out"]["code"]
    assert sha(path) == DECLARATION["held_out"]["code_sha256"], path
    spec = importlib.util.spec_from_file_location("compare", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules["compare"] = module
    spec.loader.exec_module(module)
    return module


def composition(water_mass: float, m_na: float, m_mg: float, m_ca: float, m_k: float = 0.0) -> dict[str, float]:
    amounts = (1 / water_mass, m_na, m_mg, m_ca, m_na + m_k + 2 * (m_mg + m_ca), m_k)
    return dict(zip(COMPONENTS, (n / sum(amounts) for n in amounts)))


def phi_nacl(hamer_wu: list[tuple[float, float]], m: float) -> float:
    """Source arithmetic only: three-point Lagrange interpolation of the Hamer-Wu NaCl osmotic coefficient."""
    near = sorted(hamer_wu, key=lambda point: abs(point[0] - m))[:3]
    return sum(y * math.prod((m - xk) / (xj - xk) for xk, _ in near if xk != xj) for xj, y in near)


def training(p: Parameters) -> list[dict]:
    """Rows in declaration order: source label, inputs, observed Phi and the native observation."""
    evaluated = ROOT / "analyses/2026-held-2014-chlorides/inputs/evaluated.csv"
    assert sha(evaluated) == DECLARATION["inputs_sha256"][str(evaluated.relative_to(ROOT))]
    data = list(csv.DictReader(evaluated.open()))
    rows = []
    for spec in DECLARATION["training"][:3]:
        ion, nu = SALTS[spec["salt"]]
        selected = [v for v in data if v["salt"] == spec["salt"]
                    and (spec["salt"] == "NaCl" or float(v["molality_mol_kg"]) <= 3.3)]
        assert len(selected) == spec["rows"], spec
        for v in selected:
            m = float(v["molality_mol_kg"])
            rows.append({"source": spec["source"], "salt": spec["salt"], "m_NaCl": m if ion == CATIONS[0] else 0.0,
                         "m_divalent": 0.0 if ion == CATIONS[0] else m, "observed": float(v["osmotic_evaluated"]),
                         "observation": r.electrolyte(p, "osmotic_coefficient", T, P, {"water": 1}, ions=(ion, "chloride-anion"),
                                                      stoichiometry=nu, molality=m, observed=float(v["osmotic_evaluated"]),
                                                      form="difference", scale=spec["scale"])})
    compare = comparison()
    nominal = compare.read(compare.INPUTS / "frozen-realizations.csv")[0]
    assert nominal["operation"] == "nominal"
    spec = DECLARATION["training"][3]
    mixed = [row for row in compare.realize("rard-miller-1987", nominal) if row["role"] == "mixture"]
    assert len(mixed) == spec["rows"]
    water_mass = p.molar_masses[0]
    for row in mixed:
        rows.append({"source": spec["source"], "salt": "NaCl-MgCl2", "m_NaCl": row["m_NaCl"], "m_divalent": row["m_divalent"],
                     "observed": row["observed"], "row_id": row["row_id"], "I_mol_kg": row["I_mol_kg"],
                     "observation": r.mixture_osmotic_coefficient(p, T, P, composition(water_mass, row["m_NaCl"], row["m_divalent"], 0.0),
                                                                  observed=row["observed"], scale=spec["scale"])})
    for spec in DECLARATION["training"][4:]:
        assert sha(SAAD) == DECLARATION["inputs_sha256"][str(SAAD.relative_to(ROOT))]
        hamer_wu = [(float(v["molality_mol_kg"]), float(v["osmotic_evaluated"])) for v in data if v["salt"] == "NaCl"]
        admitted = []
        for v in csv.DictReader(SAAD.open()):
            n, m_r, y, ionic = (float(v[k]) for k in ("n_KCl_per_NaCl", "reference_molality_mol_kg", "y_CaCl2_ionic_strength_fraction",
                                                      "ionic_strength_mol_kg"))
            observed = float(v["osmotic_coefficient_observed"])
            balance = 2 * m_r * phi_nacl(hamer_wu, m_r) / (ionic * (2 - y))
            if abs(observed - balance) > spec["admission_tolerance"]:
                continue
            m_a, m_ca = ionic * (1 - y), y * ionic / 3
            admitted.append({"source": spec["source"], "salt": spec["salt"], "row_id": f"{v['table']}/{v['row']}",
                             "m_NaCl": m_a / (1 + n), "m_KCl": n * m_a / (1 + n), "m_divalent": m_ca, "I_mol_kg": ionic, "observed": observed,
                             "observation": r.mixture_osmotic_coefficient(p, T, P, composition(water_mass, m_a / (1 + n), 0.0, m_ca, n * m_a / (1 + n)),
                                                                          observed=observed, scale=spec["scale"])})
        assert len(admitted) == spec["rows"], len(admitted)
        rows += admitted
    return rows


def coordinates(p: Parameters, names: list[str]) -> tuple[list, list[str]]:
    """One coordinate per cation and name; labels follow the native physical order (one per target)."""
    result, labels = [], []
    for ion in CATIONS:
        floor = BOUNDS["diameter_floor_angstrom"][ion]
        for name in names:
            if name.startswith("k_ij/"):
                other = name.split("/")[1]
                targets = [("k_ij", (ion, other), 1)]
                bounds = (BOUNDS["k_ij_numerical_lower"], BOUNDS["k_ij_upper"])
                labels.append("k_ij/" + "+".join(sorted((ion, other))))
            else:
                targets = [(name, ion, 1)] + ([("packing_diameter", ion, 0.88)] if name == "segment_diameter" else [])
                bounds = (floor, BOUNDS["diameter_numerical_upper_angstrom"])
                labels += [f"{family}/{ion}" for family, _, _ in targets]
            result.append(r.coordinate(p, targets, origin=0, scale=1, bounds=bounds))
    return result, labels


def summary(result, labels: list[str], rows: list[dict], cs: list) -> dict:
    by_source: dict[str, list[float]] = {}
    for row, value in zip(rows, result.predictions):
        by_source.setdefault(f"{row['source']}/{row['salt']}", []).append(value - row["observed"])
    covariance = result.covariance
    width = len(result.physical)
    sigma = None
    if covariance.status == "available" and covariance.variance_factor is not None:
        sigma = [math.sqrt(covariance.physical[i * width + i] * covariance.variance_factor) for i in range(width)]
    lower = [c.lower for c in cs]
    upper = [c.upper for c in cs]
    numerical = [i for i, (x, lo, hi) in enumerate(zip(result.optimizer, lower, upper))
                 if (hi == BOUNDS["diameter_numerical_upper_angstrom"] and abs(x - hi) <= 1e-6)
                 or (lo == BOUNDS["k_ij_numerical_lower"] and abs(x - lo) <= 1e-6)]
    physical_bound = [i for i, (x, lo, hi) in enumerate(zip(result.optimizer, lower, upper))
                      if i not in numerical and (abs(x - lo) <= 1e-6 or abs(x - hi) <= 1e-6)]
    return {
        "status": result.status, "message": result.message, "initial_cost": result.initial_cost,
        "final_cost": result.final_cost, "iterations": result.iterations, "trial_failures": result.trial_failures,
        "wall_seconds": sum(result.iteration_seconds), "labels": labels, "physical": list(result.physical),
        "optimizer": list(result.optimizer), "native_failures": sum(str(s) != "ObservationStatus.Available" for s in result.statuses),
        "training_rms": {k: math.sqrt(sum(v * v for v in e) / len(e)) for k, e in by_source.items()},
        "training_max": {k: max(map(abs, e)) for k, e in by_source.items()},
        "training_bias": {k: sum(e) / len(e) for k, e in by_source.items()},
        "covariance_status": covariance.status, "variance_factor": covariance.variance_factor,
        "rank": covariance.rank, "active_bounds": list(covariance.active_bounds),
        "numerical_box_coordinates": numerical, "physical_bound_coordinates": physical_bound,
        "standard_error": sigma,
    }


def optima(starts: list[dict]) -> list[dict]:
    """Distinct converged optima by cost; the declaration fixes the identity rule and the 10% band."""
    converged = sorted((s for s in starts if s["status"] == "converged"), key=lambda s: s["final_cost"])
    found: list[dict] = []
    for start in converged:
        same = next((o for o in found if max(abs(a - b) for a, b in zip(o["physical"], start["physical"])) <= 1e-3), None)
        if same:
            same["starts"].append(start["start"])
        else:
            found.append({**start, "starts": [start["start"]]})
    best = found[0]["final_cost"] if found else None
    for item in found:
        item["within_band"] = item["final_cost"] <= 1.10 * best
    return found


def declarations():
    """Construct every form's training rows, coordinates and starts without solving."""
    for spec in DECLARATION["forms"].values():
        p = Parameters.from_json(ROOT / f"data/packets/{spec['base_packet']}/parameters/parameters.json", components=COMPONENTS)
        training(p)
        cs, _ = coordinates(p, spec["coordinates_per_cation"])
        assert all(len(cs) == len([x for ion in CATIONS for x in start["values"][ion]]) for start in spec["starts"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("forms", nargs="*", default=list(DECLARATION["forms"]))
    args = parser.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    for form in args.forms:
        spec = DECLARATION["forms"][form]
        p = Parameters.from_json(ROOT / f"data/packets/{spec['base_packet']}/parameters/parameters.json", components=COMPONENTS)
        rows = training(p)
        cs, labels = coordinates(p, spec["coordinates_per_cation"])
        starts = [[x for ion in CATIONS for x in start["values"][ion]] for start in spec["starts"]]
        controls = DECLARATION["controls"]
        began = time.perf_counter()
        results = r.fit(p, cs, [row["observation"] for row in rows],
                        controls=r.FitControls(maximum_iterations=controls["maximum_iterations"],
                                               maximum_elapsed_time_seconds=controls["maximum_elapsed_time_seconds"],
                                               workers=controls["workers"]),
                        starts=starts, shared_iterations=controls["shared_iterations"], shared_seconds=controls["shared_seconds"])
        wall = time.perf_counter() - began
        records = [{"start": start["id"], **summary(result, labels, rows, cs)} for start, result in zip(spec["starts"], results)]
        found = optima(records)
        best = min(range(len(results)), key=lambda i: results[i].final_cost if records[i]["status"] == "converged" else math.inf)
        with (RESULTS / f"training-residuals-{form}.csv").open("w", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(["source", "salt", "row_id", "m_NaCl", "m_KCl", "m_divalent", "observed", "model", "residual", "status"])
            for row, value, status in zip(rows, results[best].predictions, results[best].statuses):
                writer.writerow([row["source"], row["salt"], row.get("row_id", ""), row["m_NaCl"], row.get("m_KCl", 0.0), row["m_divalent"],
                                 row["observed"], value, value - row["observed"], str(status).split(".")[-1]])
        record = {
            "form": form, "base_packet": spec["base_packet"], "declaration_sha256": sha(HERE / "inputs/fit-declaration.json"),
            "engine_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True, cwd=ROOT).strip(),
            "engine_dirty": bool(subprocess.check_output(["git", "status", "--porcelain", "engine"], text=True, cwd=ROOT).strip()),
            "rows": len(rows), "coordinates": len(cs), "wall_seconds": wall, "lowest_cost_start": records[best]["start"],
            "starts": records, "optima": [{k: o[k] for k in ("starts", "final_cost", "within_band", "physical", "training_rms",
                                                              "training_max", "standard_error", "covariance_status",
                                                              "numerical_box_coordinates", "physical_bound_coordinates")}
                                          for o in found],
        }
        (RESULTS / f"fit-{form}.json").write_text(json.dumps(record, indent=2) + "\n")
        print(form, f"wall {wall:.1f} s", "best", records[best]["start"], f"cost {records[best]['final_cost']:.6g}",
              json.dumps(records[best]["training_rms"]), f"optima {len(found)} in band {sum(o['within_band'] for o in found)}", flush=True)


if __name__ == "__main__":
    main()

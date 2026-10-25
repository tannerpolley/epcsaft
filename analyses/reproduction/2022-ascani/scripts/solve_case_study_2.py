"""Ascani 2022 Case Study 2: water/1-butanol/NaCl/KCl two-liquid split, ions in both liquids.

Python declares the Engine problems and reads results; every equation-of-state and
equilibrium calculation runs natively. Writes results/summary.json and results/phases.csv.
"""

from __future__ import annotations

import copy
import csv
import hashlib
import json
import math
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import epcsaft
from epcsaft import equilibrium as eq

CASE = Path(__file__).resolve().parents[1]
REPO = CASE.parents[2]
PACKET = REPO / "data/packets/ascani-2022-case-study-2-electrolyte-lle/3/parameters/parameters.json"
TRACER = CASE / "inputs/source/ascani-2022-case-study-2-tracer.csv"
RESULTS = CASE / "results"
C = ("water", "1-butanol", "sodium-cation", "potassium-cation", "chloride-anion")
Z = (0, 0, 1, 1, -1)
T, P = 298.15, 1.0e5
ROW_TOLERANCE = 1.0e-10
# Water/1-butanol (k_ij, l_ij) at 298.15 K. Nann 2013 Table 3: k_ij = 2.94e-4 T/K - 0.102.
# Ascani 2022 SI Table S3 prints -0.102 as k_ij(298.15 K), the intercept of that expression.
# l_ij = -0.0044 enters the SI Eq. S4 (1 - l_ij) form; Nann 2013 Eq. 2 writes (1 + l_ij), so the
# diagnostic variant flips its sign in the Engine's (1 - l_ij) form.
VARIANTS = {"nann-2013-expression": (2.94e-4 * T - 0.102, -0.0044), "si-table-s3-printed": (-0.102, -0.0044),
            "l-ij-sign-flipped": (2.94e-4 * T - 0.102, 0.0044)}
# Generic starts that do not use the paper's answer; the paper start checks start independence.
GENERIC = {"aqueous": (0.95, 0.02, 0.005, 0.005, 0.01), "organic": (0.50, 0.4996, 0.0001, 0.0001, 0.0002)}


def tracer() -> dict[str, str]:
    with TRACER.open(newline="") as stream:
        return next(csv.DictReader(stream))


def paper_phases(row: dict[str, str]) -> dict[str, tuple[float, ...]]:
    names = ("water", "butanol", "na", "k", "cl")
    return {phase: tuple(float(row[f"{phase}_x_{n}_explicit"]) for n in names) for phase in ("organic", "aqueous")}


def document(variant: str) -> dict:
    doc = json.loads(PACKET.read_text())
    values = dict(zip(("k_ij", "l_ij"), VARIANTS[variant]))
    for pair in doc["pairs"]:
        if {pair["component_id_a"], pair["component_id_b"]} == {"water", "1-butanol"}:
            for item in pair["coefficients"]:
                if item["family"] in values:
                    item["value"]["magnitude"] = values[item["family"]]
    return doc


def without_born(doc: dict) -> dict:
    doc = copy.deepcopy(doc)
    for family in doc["model_families"]:
        if family["kind"] == "electrolyte":
            family["choice"] = "fully-dissociated-debye-huckel"
            family.pop("c_shell")
            family.pop("c_dielectric")
    doc["model_coefficients"] = []
    for component in doc["components"]:
        component["coefficients"] = [i for i in component["coefficients"] if i["family"] not in ("born_diameter", "solvation_factor")]
    return doc


def log_fugacities(model, x, rho=None):
    state = model.state(T, rho=rho, x=x) if rho is not None else model.state(T, P=P, x=x, phase="liquid")
    ln_f = [math.log(xi) + lp + math.log(P / 1.0e5) for xi, lp in zip(x, state.log_fugacity_coefficient)]
    if len(x) == 2:
        return state, ln_f, {"water": ln_f[0], "1-butanol": ln_f[1]}
    # Mean ionic fugacity of a 1:1 salt: ln f_pm = (ln f_+ + ln f_-)/2 (Ascani 2022 Table 5).
    return state, ln_f, {"water": ln_f[0], "1-butanol": ln_f[1], "KCl": 0.5 * (ln_f[3] + ln_f[4]), "NaCl": 0.5 * (ln_f[2] + ln_f[4])}


def declare(variant: str, system: str, start: str):
    components = C if system == "quaternary" else C[:2]
    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(document(variant), components=components))
    row = tracer()
    if system == "quaternary":
        feed = tuple(float(row[f"feed_x_{n}_explicit"]) for n in ("water", "butanol", "na", "k", "cl"))
        guesses = GENERIC if start == "generic" else paper_phases(row)
        guesses = (guesses["aqueous"], guesses["organic"])
    else:
        feed, guesses = (0.8, 0.2), ((0.98, 0.02), (0.5, 0.5))
    problem = eq.Problem(
        phases=[eq.Phase(name, kind="liquid", support=components, composition_guess=g) for name, g in zip(("phase-0", "phase-1"), guesses)],
        T=T, P=P, feed=eq.MoleFractions(dict(zip(components, feed))), row_tolerance=ROW_TOLERANCE)
    return model, problem, feed


TASKS = [(v, "binary", "generic") for v in VARIANTS] + [(v, "quaternary", s) for v in VARIANTS for s in ("generic", "paper")]


def solve(task: tuple[str, str, str]) -> dict:
    variant, system, start = task
    model, problem, feed = declare(*task)
    n = len(feed)
    started = time.perf_counter()
    result = eq.solve_equilibrium(model, problem)
    seconds = time.perf_counter() - started
    record = {"variant": variant, "system": system, "start": start, "k_ij_water_butanol": VARIANTS[variant][0], "l_ij_water_butanol": VARIANTS[variant][1], "feed": feed,
              "success": bool(result.success), "requested_tolerance_met": result.requested_tolerance_met,
              "message": str(result.message), "solver_status": str(result.solver_status),
              "iterations": int(result.iterations), "evaluations": int(result.evaluations),
              "hessian_calls": int(result.hessian_calls), "seconds": seconds, "phases": []}
    if not result.success or len(result.amounts) != 2:
        return record
    phases = []
    for p in range(2):
        x = tuple(result.mole_fractions[p * n:(p + 1) * n])
        rho = float(result.molar_densities[p])
        state, ln_f, groups = log_fugacities(model, x, rho)
        phases.append({"amount": float(result.amounts[p]), "x": x, "rho": rho, "pressure": state.pressure,
                       "dp_drho": state.pressure_density_derivative, "permittivity": state.bulk_relative_permittivity,
                       "charge": sum(z * xi for z, xi in zip(Z, x)), "ln_f": ln_f, "ln_f_groups": groups})
    phases.sort(key=lambda phase: phase["x"][0], reverse=True)
    phases[0]["label"], phases[1]["label"] = "aqueous", "organic"
    record["phases"] = phases
    record["material_balance_max_abs"] = max(abs(sum(ph["amount"] * ph["x"][i] for ph in phases) - feed[i]) for i in range(n))
    record["charge_max_abs"] = max(abs(ph["charge"]) for ph in phases)
    record["pressure_max_rel"] = max(abs(ph["pressure"] - P) / P for ph in phases)
    record["ln_f_gap_max_abs"] = max(abs(phases[0]["ln_f_groups"][k] - phases[1]["ln_f_groups"][k]) for k in phases[0]["ln_f_groups"])
    record["classified_residual_max_abs"] = max(abs(float(v)) for v in result.classified_residuals)
    record["reduced_hessian_minimum"] = float(result.reduced_hessian_minimum)
    record["objective_g_over_rt"] = float(result.objective)
    if system == "quaternary":
        _, ln_f_feed, _ = log_fugacities(model, feed)
        g_feed = sum(z * f for z, f in zip(feed, ln_f_feed))
        g_split = sum(ph["amount"] * sum(xi * f for xi, f in zip(ph["x"], ph["ln_f"])) for ph in phases)
        record["delta_g_split_minus_feed_J_per_mol"] = 8.314462618 * T * (g_split - g_feed)
        try:
            options = eq.StabilityOptions()
            report = eq.stability_test(eq.compile_problem(model, problem), result, options)
            record["stability"] = {"completed": report.completed, "message": report.message, "trials": len(report.trials)}
        except epcsaft.ParameterError as error:
            record["stability"] = {"refused": str(error)}
    return record


def paper_point_diagnostics() -> dict:
    """The Engine model evaluated at the paper's printed phases (no solve)."""
    row, out = tracer(), {}
    phases = paper_phases(row)
    paper = {"water": float(row["ln_f_water_bar"]), "1-butanol": float(row["ln_f_butanol_bar"]),
             "KCl": float(row["ln_f_kcl_pair_bar"]), "NaCl": float(row["ln_f_nacl_pair_bar"])}
    doc = document("nann-2013-expression")
    full = epcsaft.Mixture(epcsaft.Parameters.from_mapping(doc, components=C))
    bare = epcsaft.Mixture(epcsaft.Parameters.from_mapping(without_born(doc), components=C))
    for name, x in phases.items():
        state, _, groups = log_fugacities(full, x)
        _, _, no_born = log_fugacities(bare, x)
        born = {k: groups[k] - no_born[k] for k in groups}
        out[name] = {"ln_f": groups, "born_part": born, "permittivity": state.bulk_relative_permittivity,
                     "born_scale_matching_paper": {k: (paper[k] - no_born[k]) / born[k] for k in ("KCl", "NaCl")}}
    # Each formula-component balance z_i = b x_i(org) + (1 - b) x_i(aq) alone fixes the organic mole share b.
    # Formula basis from the renormalized explicit columns: x_salt = x_cation / (1 - x_Cl).
    feed = [float(row[f"feed_x_{n}_explicit"]) for n in ("water", "butanol", "na", "k", "cl")]
    z, org, aq = ([v / (1 - x[4]) for v in x[:4]] for x in (feed, phases["organic"], phases["aqueous"]))
    out["printed_phase_organic_share_by_balance"] = {n: (z[i] - aq[i]) / (org[i] - aq[i]) for i, n in enumerate(("water", "1-butanol", "NaCl", "KCl"))}
    out["paper_ln_f"] = paper
    return out


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    started = time.perf_counter()
    with ProcessPoolExecutor(max_workers=min(len(TASKS), os.cpu_count() or 1)) as pool:
        records = list(pool.map(solve, TASKS))
    # Phase discovery is the Engine's own stability-driven search; record whether it admits this model.
    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(document("nann-2013-expression"), components=C))
    row = tracer()
    feed = {c: float(row[f"feed_x_{n}_explicit"]) for c, n in zip(C, ("water", "butanol", "na", "k", "cl"))}
    try:
        eq.solve_equilibrium(model, eq.Problem(phases=eq.Discover(), T=T, P=P, feed=eq.MoleFractions(feed), row_tolerance=ROW_TOLERANCE))
        discovery = "admitted"
    except epcsaft.ParameterError as error:
        discovery = str(error)
    summary = {
        "schema": "ascani-2022-case-study-2.v1",
        "temperature_K": T, "pressure_Pa": P, "components": C, "row_tolerance": ROW_TOLERANCE,
        "parameters": {"path": str(PACKET.relative_to(REPO)), "sha256": sha256(PACKET), "water_butanol_k_ij_l_ij": VARIANTS},
        "tracer_sha256": sha256(TRACER),
        "engine": epcsaft.provenance(),
        "git_commit": subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True).stdout.strip(),
        "cpu_count": os.cpu_count(),
        "phase_discovery": discovery,
        "solves": records,
        "paper_point_diagnostics": paper_point_diagnostics(),
        "wall_seconds": time.perf_counter() - started,
    }
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    with (RESULTS / "phases.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["variant", "system", "start", "phase", "amount", *[f"x_{c}" for c in C], "rho_mol_m3", "permittivity", "charge"])
        for r in records:
            for ph in r["phases"]:
                x = list(ph["x"]) + [""] * (len(C) - len(ph["x"]))
                writer.writerow([r["variant"], r["system"], r["start"], ph["label"], ph["amount"], *x, ph["rho"], ph["permittivity"], ph["charge"]])
    print(json.dumps({"wall_seconds": summary["wall_seconds"], "discovery": discovery,
                      "solves": [(r["variant"], r["system"], r["start"], r["success"], r["message"][:60]) for r in records]}, indent=1))


def declarations():
    return [declare(*task) for task in TASKS]


if __name__ == "__main__":
    main()

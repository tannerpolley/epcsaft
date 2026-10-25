"""Test the retained declared candidates with native ALG-STABILITY.

The shared Ipopt driver owns trial enumeration, TPD minimization, reconstruction,
uncertainty and finite donor splits. Every trial is retained, including failures.
Density-independent charged families retain the adopted explicit refusal.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path

import epcsaft
from epcsaft import equilibrium

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parents[1]
METHANE_ETHANE = ROOT / "data/packets/gross-2001-methane-ethane/3/parameters/parameters.json"
TERNARY = HERE / "inputs/phase-stability-neutral-ternary.json"
CHARGED = HERE / "inputs/phase-stability-charged-singleton.json"
KHUDAIDA = ROOT / "data/packets/khudaida-2026-figure-2-electrolyte-lle/6/parameters/parameters.json"
HELD2 = ROOT / "analyses/reproduction/2026-khudaida/results/data/khudaida-2026-figure-2-held2.csv"
KHUDAIDA_C = ("water", "ethanol", "isobutanol", "sodium-cation", "chloride-anion")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_rows(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    output = parser.parse_args().output_dir
    output.mkdir(parents=True, exist_ok=False)
    rows, trials, quantities, residuals, comparisons = [], [], [], [], []
    historical = list(csv.DictReader((HERE / "results/current-engine.csv").open()))
    historical_identity = json.loads((HERE / "results/summary.json").read_text())["identity"]
    load = lambda path, components=None: epcsaft.Mixture(epcsaft.Parameters.from_json(path, components=components))

    def retain_quantities(case, entry, values):
        quantities.extend(dict(case=case, entry=entry, quantity=k, value=v) for k, v in values.items())

    def record(case, label, model, T, P, feed, guesses, kind="liquid"):
        ids = model.component_ids
        problem = equilibrium.Problem(
            phases=[equilibrium.Phase(f"phase-{k}", kind=kind, support=ids, composition_guess=g)
                    for k, g in enumerate(guesses)],
            T=T, P=P, feed=equilibrium.Amounts(dict(zip(ids, feed, strict=True))), row_tolerance=1e-10)
        result = equilibrium.solve_equilibrium(model, problem)
        key = f"{case}/{label}"
        for r, raw, classified in zip(result.rows, result.residuals, result.classified_residuals, strict=True):
            residuals.append(dict(case=key, kind=str(r.kind), phase=r.phase, index=r.index,
                                  raw=raw, classified=classified, eqid=r.eqid))
        if not result.success or not result.requested_tolerance_met:
            write_rows(output / "declared-residuals.csv", residuals)
            raise RuntimeError(f"{key}: {result.message}")
        chemical = {(r.phase, r.index): v for r, v in zip(result.rows, result.residuals, strict=True)
                    if str(r.kind).endswith("Chemical")}
        neutral_groups = [[i] for i, name in enumerate(ids) if name not in ("sodium-cation", "chloride-anion")]
        if "sodium-cation" in ids:
            neutral_groups.append([ids.index("sodium-cation"), ids.index("chloride-anion")])
        chemical_gap = max((abs(sum(chemical[p, i]-chemical[0, i] for i in group))
                            for p in range(1, result.phase_count) for group in neutral_groups), default=0.0)
        material_defect = max(abs(v) for r, v in zip(result.rows, result.residuals, strict=True)
                              if str(r.kind).endswith("Moiety"))
        compiled = equilibrium.compile_problem(model, problem)
        report = None
        refusal = ""
        try:
            report = equilibrium.stability_test(compiled, result, equilibrium.StabilityOptions())
        except epcsaft.ParameterError as error:
            if str(error) != "permittivity model has no low-density limit; declare phases":
                raise
            refusal = str(error)
        n = len(ids)
        x = [list(result.mole_fractions[p*n:(p+1)*n]) for p in range(result.phase_count)]
        row = dict(case_id=case, declared_topology=label, T_K=T, P_Pa=P, components=json.dumps(ids),
                   feed_mol=json.dumps(feed), phase_count=result.phase_count, solver_success=result.success,
                   solver_message=result.message, g_over_rt_mol=result.objective,
                   classified_residual_max_abs=max(map(abs, result.classified_residuals)),
                   phase_mole_fractions=json.dumps(x), phase_amounts_mol=json.dumps(result.amounts),
                   molar_densities_mol_m3=json.dumps(result.molar_densities),
                   stability_verdict=refusal if report is None else report.message,
                   trial_count=0 if report is None else len(report.trials),
                   dmu_max_abs=chemical_gap, material_max_abs=material_defect,
                   coordinates=json.dumps(result.coordinates), multipliers=json.dumps(result.multipliers),
                   constraint_multipliers=json.dumps(result.constraint_multipliers),
                   stability_completed=False if report is None else report.completed,
                   candidate_coordinates="" if report is None else json.dumps(report.candidate_coordinates),
                   standard_over_rt="" if report is None else json.dumps(report.standard_over_rt),
                   plane_over_rt="" if report is None else json.dumps(report.plane_over_rt))
        if report is not None:
            retain_quantities(key, "work", report.work)
            for p, phase in enumerate(report.present_phases):
                retain_quantities(key, f"present-{p}", phase)
            for i, trial in enumerate(report.trials):
                q = trial.quantities
                for name in ("charge_residual", "normalization_residual"):
                    if name in q and abs(q[name]) > 1e-12:
                        raise RuntimeError(f"{key}/trial-{i}: {name}={q[name]}")
                trials.append(dict(case=key, trial=i, kind=str(trial.phase.kind), support=json.dumps(trial.phase.support),
                                   seed=trial.seed, evidence=trial.evidence, message=trial.message,
                                   composition=json.dumps(trial.composition), donor=trial.donor,
                                   resolved_descent=trial.resolved_descent, optimizer_success=trial.optimizer.success,
                                   optimizer_message=trial.optimizer.message, optimizer_status=trial.optimizer.solver_status,
                                   iterations=trial.optimizer.iterations, evaluations=trial.optimizer.evaluations,
                                   hessian_calls=trial.optimizer.hessian_calls,
                                   raw_coordinates=json.dumps(trial.raw_coordinates), raw_residuals=json.dumps(trial.raw_residuals),
                                   optimizer_classified_residuals=json.dumps(trial.optimizer.classified_residuals),
                                   split_coordinates=json.dumps(trial.split_coordinates), basis=json.dumps(trial.basis),
                                   row_transformation=json.dumps(trial.row_transformation),
                                   optimizer_state=json.dumps({name: getattr(trial.optimizer, name) for name in
                                       ("temperature", "pressure", "amounts", "mole_fractions", "molar_densities",
                                        "packing_fractions", "multipliers", "constraint_multipliers", "coordinates",
                                        "objective", "reduced_hessian_minimum", "requested_tolerance_met")}
                                       | {name: [] if (v := getattr(trial.optimizer, name)) is None else [v] for name in ("temperature", "pressure")}),
                                   convergence=json.dumps({name: str(getattr(trial.optimizer.convergence, name)) for name in
                                       ("residual", "floor", "amplification", "contraction", "iterations", "evaluations", "stop")})))
                retain_quantities(key, f"trial-{i}", q)
                for j, split in enumerate(trial.split_checks):
                    retain_quantities(key, f"trial-{i}/split-{j}", split)
            for branch in ("Liquid", "Vapor"):
                available = [t for t in report.trials if str(t.phase.kind).endswith(branch)
                             and "tpd_over_rt" in t.quantities and math.isfinite(t.quantities["tpd_over_rt"])]
                if available:
                    best = min(available, key=lambda t: t.quantities["tpd_over_rt"])
                    row[f"tpd_min_{branch.lower()}"] = best.quantities["tpd_over_rt"]
                    row[f"tpd_argmin_{branch.lower()}"] = json.dumps(best.composition)
        rows.append(row)
        old = next(r for r in historical if r["case_id"] == case and r["declared_topology"] == label)
        # Presentation-only conversion of the retained pressure-based Gibbs reference
        # to the native density-based reference (rho0=1 mol/m3), at pressure closure.
        old_native = sum(feed) * (float(old["g_rt"]) + math.log(P / (epcsaft._core.GAS_CONSTANT_J_PER_MOL_K * T)))
        comparisons.append(dict(case=key, T_K=T, P_Pa=P, before_g_over_rt_mol=old_native,
                                after_g_over_rt_mol=result.objective, delta_g_over_rt_mol=result.objective-old_native,
                                before_dmu_max_abs=old["dmu_max_abs"], before_material_max_abs=old["material_max_abs"],
                                after_classified_residual_max_abs=row["classified_residual_max_abs"],
                                after_dmu_max_abs=chemical_gap, after_material_max_abs=material_defect,
                                before_phase_mole_fractions=old["phase_mole_fractions"], after_phase_mole_fractions=row["phase_mole_fractions"],
                                before_phase_amount_fractions=old["phase_amounts"], after_phase_amounts_mol=row["phase_amounts_mol"],
                                before_tpd_min_liquid=old["tpd_min_liquid"], after_tpd_min_liquid=row.get("tpd_min_liquid", ""),
                                before_tpd_min_vapor=old["tpd_min_vapor"], after_tpd_min_vapor=row.get("tpd_min_vapor", ""),
                                before_argmin_liquid=old["tpd_argmin_liquid"], after_argmin_liquid=row.get("tpd_argmin_liquid", ""),
                                before_argmin_vapor=old["tpd_argmin_vapor"], after_argmin_vapor=row.get("tpd_argmin_vapor", ""),
                                before_grid_points_liquid=old["grid_points_liquid"], before_grid_points_vapor=old["grid_points_vapor"],
                                after_trial_count=row["trial_count"], after_verdict=row["stability_verdict"]))
        print(json.dumps(row), flush=True)

    feed7 = next(json.loads(r["explicit_feed_json"]) for r in csv.DictReader(HELD2.open()) if r["point_id"] == "7")
    record("NEUTRAL_SINGLETON", "1 vapor", load(METHANE_ETHANE, ("methane", "ethane")),
           300.0, 1e5, (0.5, 0.5), [(0.5, 0.5)], "vapor")
    record("CHARGED_SINGLETON", "1 liquid", load(CHARGED), 298.15, 1e5,
           (0.98, 0.01, 0.01), [(0.98, 0.01, 0.01)])
    z = (1/3, 1/3, 1/3)
    for kind in ("vapor", "liquid"):
        record("NEUTRAL_TERNARY_IMMISCIBLE", f"1 {kind}", load(TERNARY), 300.0, 1e5, (1., 1., 1.), [z], kind)
    record("NEUTRAL_TERNARY_IMMISCIBLE", "2 liquids (alpha-rich / rest)", load(TERNARY), 300.0, 1e5,
           (1., 1., 1.), [(0.90, 0.05, 0.05), (0.05, 0.475, 0.475)])
    guesses = [(0.90, 0.06, 0.04), (0.11, 0.83, 0.06), (0.14, 0.12, 0.74)]
    for P in (1e5, 1e7):
        record("NEUTRAL_TERNARY_IMMISCIBLE", f"3 liquids at {P/1e6:g} MPa", load(TERNARY),
               300.0, P, (1., 1., 1.), guesses)
    record("KHUDAIDA_ROW_7", "2 liquids", load(KHUDAIDA, KHUDAIDA_C), 293.15, 1e5, tuple(feed7),
           [(.8, .05, .05, .05, .05), (.2, .1, .6998, .0001, .0001)])
    for name, table in (("current-engine", rows), ("trials", trials), ("quantities", quantities),
                        ("declared-residuals", residuals), ("before-after", comparisons)):
        write_rows(output / f"{name}.csv", table)
    wheel = Path(json.loads(importlib.metadata.distribution("epcsaft").read_text("direct_url.json"))["url"].removeprefix("file://"))
    (output / "summary.json").write_text(json.dumps(dict(
        identity=dict(wheel_sha256=sha(wheel), historical_identity=historical_identity,
                      script_sha256=sha(Path(__file__)), inputs={str(p.relative_to(ROOT)): sha(p) for p in
                      (METHANE_ETHANE, TERNARY, CHARGED, KHUDAIDA, HELD2)},
                      outputs={p.name: sha(p) for p in output.glob("*.csv")}),
        method="equilibrium.stability_test at each native declared candidate; default native trials and budgets",
        row_tolerance=1e-10, trial_charge_tolerance=1e-12,
        claim_limits=["finite trial set, estimated uncertainty; no global stability certificate or physical validation",
                      "density-independent ionic families refused; no charged stability verdict inferred",
                      "historical and native trial seeds differ; explicit mole fractions replace formula-salt coordinates",
                      "historical Gibbs reference converted to native density reference for comparison; no identical-wheel replay",
                      "no discovery solve, new topology, solid phase or parameter adoption"]), indent=2) + "\n")
    return 0


def declarations():
    from epcsaft import Mixture, Parameters

    return [
        Mixture(Parameters.from_json(p, components=c))
        for p, c in (
            (METHANE_ETHANE, ('methane', 'ethane')),
            (CHARGED, None),
            (TERNARY, None),
            (KHUDAIDA, KHUDAIDA_C),
        )
    ]


if __name__ == "__main__":
    raise SystemExit(main())

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import epcsaft
from epcsaft import equilibrium
from idaes.core import FlowsheetBlock, UnitModelBlockData, declare_process_block_class
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.common.config import ConfigValue
from pyomo.environ import (ConcreteModel, Constraint, Expression, NonNegativeReals, Objective, Param, Var, check_optimal_termination, value, units as pyunits)
from pyomo.network import Port

CAMPAIGN = Path(__file__).resolve().parents[1]
REPO = CAMPAIGN.parents[1]
PACKET = CAMPAIGN / "inputs/data/khudaida-2026-figure-2-electrolyte-lle/7"
LOCK = CAMPAIGN / "inputs/packet-lock.json"
SOURCE = REPO / "analyses/reproduction/2026-khudaida/figures/figure-02/source/source-points.csv"
PARAMETERS = PACKET / "parameters/parameters.json"
COMPONENTS = ("water", "ethanol", "isobutanol", "sodium-cation", "chloride-anion")
PHASES = ("aqueous", "organic")
CHARGES = (0.0, 0.0, 0.0, 1.0, -1.0)
TEMPERATURE = 293.15
PRESSURE = 100000.0
ROW_TOLERANCE = 1.0e-10
GAP = 1.0e-3
SOLVER_OPTIONS = {"tol": 1.0e-9, "constr_viol_tol": 1.0e-10, "max_iter": 200}
CSV_FIELDS = (
    "point_id", "phase", "engine_phase_ordinal", "temperature_K", "pressure_Pa",
    "feed_formula_water", "feed_formula_ethanol", "feed_formula_isobutanol", "feed_formula_nacl",
    "inlet_flow_mol_s", "phase_amount_mol", "phase_fraction", "flow_mol_s",
    "x_water", "x_ethanol", "x_isobutanol", "x_sodium_cation", "x_chloride_anion",
    "density_mol_m3", "volume_m3", "state_pressure_Pa", "pressure_density_derivative",
    "engine_solver_status", "engine_message", "engine_g_rt",
    "engine_material_balance_max_abs", "engine_phase_charge_max_abs", "engine_pressure_max_relative", "engine_classified_residual_max_abs", "engine_reduced_hessian_minimum",
    "idaes_component_balance_max_abs", "idaes_charge_residual_abs", "idaes_phase_fraction_abs", "idaes_phase_flow_closure_abs", "idaes_composition_max_abs",
    "idaes_dof", "idaes_termination", "status",
)
RESIDUAL_FIELDS = ("point_id", "eqid", "kind", "phase", "index", "raw_residual", "classified_residual")

def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def atomic_json(path: Path, record: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)

def atomic_csv(path: Path, rows: list[dict[str, object]], fields: tuple[str, ...] = CSV_FIELDS) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)

def packet_fingerprint(root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix().encode()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(path.stat().st_size.to_bytes(8, "big"))
        digest.update(path.read_bytes())
    return digest.hexdigest()

def validate_references() -> dict[str, dict[str, str]]:
    records = json.loads((CAMPAIGN / "inputs/references.json").read_text(encoding="utf-8"))
    if not isinstance(records, dict) or not records:
        raise RuntimeError("references.json must contain path/hash records")
    validated = {}
    for name, record in records.items():
        relative = Path(str(record.get("path", "")))
        path = (REPO / relative).resolve()
        if relative.is_absolute() or ".." in relative.parts or not path.is_relative_to(REPO):
            raise RuntimeError(f"reference path escapes repository: {name}")
        expected = str(record.get("sha256", ""))
        if not path.is_file() or sha256(path) != expected:
            raise RuntimeError(f"reference path/hash mismatch: {name} ({path})")
        validated[name] = {"path": str(relative), "sha256": expected}
    return validated

def input_identity() -> dict[str, object]:
    references = validate_references()
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    for item in lock["files"]:
        path = PACKET / item["path"]
        if not path.is_file() or sha256(path) != item["sha256"]:
            raise RuntimeError(f"packet materialization changed: {path}")
    if packet_fingerprint(PACKET) != lock["packet_fingerprint"]:
        raise RuntimeError("packet fingerprint changed")
    return {
        "packet": {**lock, "computed_fingerprint": packet_fingerprint(PACKET)},
        "source_files": references,
    }

def solver_identity(solver) -> dict[str, object]:
    executable = Path(str(solver.executable())).resolve()
    version = subprocess.check_output((str(executable), "-v"), text=True, stderr=subprocess.STDOUT).strip()
    return {"executable": str(executable), "sha256": sha256(executable), "version": version, "options": SOLVER_OPTIONS}

def campaign_identity(solver) -> dict[str, object]:
    distribution = importlib.metadata.distribution("epcsaft")
    direct = json.loads(distribution.read_text("direct_url.json") or "{}")
    wheel_url = str(direct.get("url", ""))
    wheel = Path(wheel_url.removeprefix("file://")).resolve()
    if wheel.suffix != ".whl" or not wheel.is_file():
        raise RuntimeError("epcsaft must be installed from one local wheel")
    return {
        "source_commit": subprocess.check_output(("git", "log", "-1", "--format=%H", "--", "engine"), cwd=REPO, text=True).strip(),
        "source_dirty": bool(subprocess.check_output(("git", "status", "--porcelain", "--", "engine"), cwd=REPO, text=True).strip()),
        "engine_wheel": str(wheel), "engine_wheel_sha256": sha256(wheel),
        "epcsaft_import": str(Path(epcsaft.__file__).resolve()),
        "epcsaft_core_import": str(Path(epcsaft._core.__file__).resolve()),
        "epcsaft_version": importlib.metadata.version("epcsaft"),
        "idaes_version": importlib.metadata.version("idaes-pse"),
        "pyomo_version": importlib.metadata.version("pyomo"),
        "solver": solver_identity(solver), "inputs": input_identity(),
        "units": {"flows": "mol/s", "temperature": "K", "pressure": "Pa", "composition": "mol/mol", "density": "mol/m^3"},
    }

def feed_rows() -> dict[int, tuple[float, ...]]:
    result = {}
    with SOURCE.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row["dataset"] != "feed":
                continue
            formula = tuple(float(row[name]) for name in ("x_water", "x_ethanol", "x_isobutanol", "x_nacl"))
            explicit = (*formula[:3], formula[3], formula[3])
            total = math.fsum(explicit)
            result[int(row["point_id"])] = tuple(item / total for item in explicit)
    if sorted(result) != list(range(1, 9)):
        raise RuntimeError("source feed must contain exactly Figure 2 points 1 through 8")
    return result

def engine_problem(feed: tuple[float, ...]) -> equilibrium.Problem:
    values = dict(zip(COMPONENTS, feed, strict=True))
    return equilibrium.Problem(
        phases=[equilibrium.Phase("phase-0", kind="liquid", support=COMPONENTS, composition_guess=(0.80, 0.05, 0.05, 0.05, 0.05)), equilibrium.Phase("phase-1", kind="liquid", support=COMPONENTS, composition_guess=(0.20, 0.10, 0.6998, 0.0001, 0.0001))],
        T=TEMPERATURE, P=PRESSURE, feed=equilibrium.Amounts(values), row_tolerance=ROW_TOLERANCE,
    )

def declarations():
    return epcsaft.Mixture(epcsaft.Parameters.from_json(PARAMETERS, components=COMPONENTS)), [engine_problem(feed) for feed in feed_rows().values()]

def engine_row_evidence(result) -> list[dict[str, object]]:
    raw = list(result.residuals)
    classified = list(result.classified_residuals)
    assert len(result.rows) == len(raw) == len(classified)
    return [{"eqid": str(row.eqid), "kind": getattr(row.kind, "name", str(row.kind)), "phase": int(row.phase), "index": int(row.index), "raw_residual": raw[index] if index < len(raw) else None, "classified_residual": classified[index] if index < len(classified) else None} for index, row in enumerate(result.rows)]

def phase_data(model, result, feed: tuple[float, ...]) -> tuple[list[dict[str, object]], dict[str, object], list[str]]:
    n = len(COMPONENTS)
    failures = []
    metrics: dict[str, object] = {"engine_solver_status": str(result.solver_status), "engine_message": str(result.message), "engine_success": bool(result.success), "engine_requested_tolerance_met": result.requested_tolerance_met, "engine_phase_count": int(result.phase_count), "row_evidence": engine_row_evidence(result)}
    if result.phase_count != 2:
        failures.append(f"engine returned invalid phase count {result.phase_count}")
        return [], metrics, failures
    if len(result.mole_fractions) != n * 2 or len(result.amounts) != 2 or len(result.molar_densities) != 2:
        failures.append("engine returned incomplete phase arrays")
        return [], metrics, failures
    if not result.success or result.requested_tolerance_met is not True:
        failures.append(f"engine solve status={result.solver_status} success={result.success} tolerance={result.requested_tolerance_met}")
    phases = []
    for ordinal in range(2):
        x = tuple(float(v) for v in result.mole_fractions[ordinal * n:(ordinal + 1) * n])
        rho = float(result.molar_densities[ordinal])
        state = model.state(TEMPERATURE, rho=rho, x=x)
        phases.append({"ordinal": ordinal, "x": x, "amount": float(result.amounts[ordinal]), "rho": rho, "state": state, "volume": float(result.amounts[ordinal]) / rho})
        if not math.isfinite(state.pressure_density_derivative) or state.pressure_density_derivative <= 0.0:
            failures.append(f"phase {ordinal} has nonpositive pressure_density_derivative")
    labels = sorted(phases, key=lambda phase: phase["x"][0], reverse=True)
    if labels[0]["x"][0] == labels[1]["x"][0]:
        failures.append("phase water ordering is not unique")
        return [], metrics, failures
    labels[0]["label"] = "aqueous"
    labels[1]["label"] = "organic"
    if max(abs(a - b) for a, b in zip(phases[0]["x"], phases[1]["x"], strict=True)) < GAP:
        failures.append("declared liquid phases are not distinct")
    reconstructed = [math.fsum(phase["amount"] * phase["x"][i] for phase in phases) for i in range(n)]
    material = max(abs(a - b) for a, b in zip(reconstructed, feed, strict=True))
    charge = max(abs(phase["x"][3] - phase["x"][4]) for phase in phases)
    pressure = max(abs(phase["state"].pressure - PRESSURE) / PRESSURE for phase in phases)
    classified = max(abs(float(v)) for v in result.classified_residuals) if result.classified_residuals else math.inf
    if material > 1e-10 or charge > 1e-10 or pressure > 1e-8 or classified > ROW_TOLERANCE or not math.isfinite(result.reduced_hessian_minimum) or result.reduced_hessian_minimum <= 0.0:
        failures.append("engine certificate failed")
    metrics.update({"material": material, "charge": charge, "pressure": pressure, "classified": classified,
                    "hessian": float(result.reduced_hessian_minimum), "g_rt": float(result.objective)})
    return phases, metrics, failures


@declare_process_block_class("KhudaidaSettler")
class KhudaidaSettlerData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare("inlet_flows", ConfigValue(default=()))
    CONFIG.declare("phase_data", ConfigValue(default=()))

    def build(self):
        super().build()
        flow_units = pyunits.mol / pyunits.s
        self.component_order = COMPONENTS
        inlet = tuple(self.config.inlet_flows)
        data = {name: record for name, record in self.config.phase_data}
        self.flow_mol_in = Var(COMPONENTS, domain=NonNegativeReals, units=flow_units, initialize=dict(zip(COMPONENTS, inlet, strict=True)))
        for name, amount in zip(COMPONENTS, inlet, strict=True):
            self.flow_mol_in[name].fix(amount)
        self.temperature_in = Var(units=pyunits.K, initialize=TEMPERATURE)
        self.pressure_in = Var(units=pyunits.Pa, initialize=PRESSURE)
        self.temperature_in.fix(TEMPERATURE)
        self.pressure_in.fix(PRESSURE)
        self.inlet_flow_total = Expression(expr=sum(self.flow_mol_in.values()))
        self.accepted_phase_amount = Param(PHASES, initialize={p: data[p]["amount"] for p in PHASES})
        self.accepted_phase_mole_fraction = Param(PHASES, COMPONENTS, initialize={(p, COMPONENTS[i]): data[p]["x"][i] for p in PHASES for i in range(len(COMPONENTS))})
        self.mapping_factor = Param(PHASES, COMPONENTS, initialize={(p, COMPONENTS[i]): data[p]["amount"] * data[p]["x"][i] / math.fsum(inlet) for p in PHASES for i in range(len(COMPONENTS))})
        self.flow_mol_aqueous = Var(COMPONENTS, domain=NonNegativeReals, units=flow_units, initialize=0.5)
        self.flow_mol_organic = Var(COMPONENTS, domain=NonNegativeReals, units=flow_units, initialize=0.5)
        self.temperature_out = Var(PHASES, units=pyunits.K, initialize=TEMPERATURE)
        self.pressure_out = Var(PHASES, units=pyunits.Pa, initialize=PRESSURE)
        self.inlet, self.aqueous, self.organic = Port(), Port(), Port()
        for port, flow, temperature, pressure in (
            (self.inlet, self.flow_mol_in, self.temperature_in, self.pressure_in),
            (self.aqueous, self.flow_mol_aqueous, self.temperature_out["aqueous"], self.pressure_out["aqueous"]),
            (self.organic, self.flow_mol_organic, self.temperature_out["organic"], self.pressure_out["organic"]),
        ):
            port.add(flow, "flow_mol")
            port.add(temperature, "temperature")
            port.add(pressure, "pressure")
        self.aqueous_flow_mapping = Constraint(COMPONENTS, rule=lambda b, i: b.flow_mol_aqueous[i] == b.inlet_flow_total * b.mapping_factor["aqueous", i])
        self.organic_flow_mapping = Constraint(COMPONENTS, rule=lambda b, i: b.flow_mol_organic[i] == b.inlet_flow_total * b.mapping_factor["organic", i])
        self.temperature_mapping = Constraint(PHASES, rule=lambda b, p: b.temperature_out[p] == TEMPERATURE * pyunits.K)
        self.pressure_mapping = Constraint(PHASES, rule=lambda b, p: b.pressure_out[p] == PRESSURE * pyunits.Pa)
        self.component_balance_residual = Expression(COMPONENTS, rule=lambda b, i: b.flow_mol_in[i] - b.flow_mol_aqueous[i] - b.flow_mol_organic[i])
        self.charge_residual = Expression(PHASES, rule=lambda b, p: sum(CHARGES[i] * (b.flow_mol_aqueous[COMPONENTS[i]] if p == "aqueous" else b.flow_mol_organic[COMPONENTS[i]]) for i in range(len(COMPONENTS))))


def solve_settler(feed: tuple[float, ...], phases: list[dict[str, object]], solver) -> dict[str, object]:
    data = [(name, next(phase for phase in phases if phase["label"] == name)) for name in PHASES]
    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)
    model.fs.settler = KhudaidaSettler(dynamic=False, has_holdup=False, inlet_flows=feed, phase_data=data)
    model.objective = Objective(expr=0.0)
    dof = degrees_of_freedom(model)
    if dof != 0:
        raise RuntimeError(f"settler degrees of freedom is {dof}")
    result = solver.solve(model, tee=False)
    unit = model.fs.settler
    if not check_optimal_termination(result):
        raise RuntimeError(f"IDAES solver did not reach optimal termination: {result.solver.termination_condition}")
    flows = {
        "aqueous": [value(unit.flow_mol_aqueous[i]) for i in COMPONENTS],
        "organic": [value(unit.flow_mol_organic[i]) for i in COMPONENTS],
    }
    inlet_total = value(unit.inlet_flow_total)
    checks = {}
    checks["component_balance"] = max(abs(value(unit.component_balance_residual[i])) for i in COMPONENTS)
    checks["charge"] = max(abs(value(unit.charge_residual[p])) for p in PHASES)
    checks["phase_fraction"] = max(
        abs(math.fsum(flows[p]) / inlet_total - next(item["amount"] for item in phases if item["label"] == p))
        for p in PHASES
    )
    checks["phase_flow_closure"] = abs(math.fsum(math.fsum(flows[p]) for p in PHASES) - inlet_total)
    checks["composition"] = max(
        abs(flows[p][i] / math.fsum(flows[p]) - next(item for item in phases if item["label"] == p)["x"][i])
        for p in PHASES for i in range(len(COMPONENTS))
    )
    if max(checks.values()) > 1e-12:
        raise RuntimeError(f"IDAES mapping check failed: {checks}")
    return {"dof": dof, "termination": str(result.solver.termination_condition), "status": str(result.solver.status), "flows": flows, "inlet_total": inlet_total, **checks}


def row_records(point: int, feed: tuple[float, ...], phases: list[dict[str, object]], metrics: dict[str, object], stage: dict[str, object]) -> list[dict[str, object]]:
    rows = []
    formula_total = math.fsum((*feed[:3], feed[3]))
    formula = tuple(item / formula_total for item in (*feed[:3], feed[3]))
    for phase in phases:
        p = phase["label"]
        flows = stage["flows"][p]
        rows.append({
            "point_id": point, "phase": p, "engine_phase_ordinal": phase["ordinal"], "temperature_K": TEMPERATURE, "pressure_Pa": PRESSURE,
            "feed_formula_water": formula[0], "feed_formula_ethanol": formula[1], "feed_formula_isobutanol": formula[2], "feed_formula_nacl": formula[3],
            "inlet_flow_mol_s": value(stage["inlet_total"]), "phase_amount_mol": phase["amount"], "phase_fraction": phase["amount"], "flow_mol_s": math.fsum(flows),
            "x_water": phase["x"][0], "x_ethanol": phase["x"][1], "x_isobutanol": phase["x"][2], "x_sodium_cation": phase["x"][3], "x_chloride_anion": phase["x"][4],
            "density_mol_m3": phase["rho"], "volume_m3": phase["volume"], "state_pressure_Pa": phase["state"].pressure, "pressure_density_derivative": phase["state"].pressure_density_derivative,
            "engine_solver_status": metrics["engine_solver_status"], "engine_message": metrics["engine_message"], "engine_g_rt": metrics["g_rt"],
            "engine_material_balance_max_abs": metrics["material"], "engine_phase_charge_max_abs": metrics["charge"], "engine_pressure_max_relative": metrics["pressure"], "engine_classified_residual_max_abs": metrics["classified"], "engine_reduced_hessian_minimum": metrics["hessian"],
            "idaes_component_balance_max_abs": stage["component_balance"], "idaes_charge_residual_abs": stage["charge"], "idaes_phase_fraction_abs": stage["phase_fraction"], "idaes_phase_flow_closure_abs": stage["phase_flow_closure"], "idaes_composition_max_abs": stage["composition"], "idaes_dof": stage["dof"], "idaes_termination": stage["termination"], "status": "accepted",
        })
    return rows


def run(output: Path) -> int:
    if output.exists() and any(output.iterdir()):
        raise RuntimeError(f"output directory must be new or empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    completed: list[int] = []
    rows: list[dict[str, object]] = []
    failures: list[str] = []
    engine_points: list[dict[str, object]] = []
    residual_rows: list[dict[str, object]] = []
    summary: dict[str, object] = {
        "schema": "epcsaft.khudaida.idaes-settler.v1",
        "status": "running",
        "completed_point_ids": [],
        "unattempted_point_ids": list(range(1, 9)),
        "failures": failures,
        "claim_limits": ["declared two-liquid local branch only", "no phase-discovery/globality claim", "fixed IDAES mapping is not a live equilibrium callback"],
    }
    atomic_json(output / "summary.json", summary)
    atomic_csv(output / "stage-results.csv", rows)
    atomic_csv(output / "engine-residuals.csv", residual_rows, RESIDUAL_FIELDS)
    exit_code = 1
    stage = "identity"
    try:
        solver = get_solver(solver_options=SOLVER_OPTIONS)
        identity = campaign_identity(solver)
        atomic_json(output / "identity.json", identity)
        summary.update({"identity": identity})
        feeds = feed_rows()
        model = epcsaft.Mixture(epcsaft.Parameters.from_json(PARAMETERS, components=COMPONENTS))
        previous = None
        for point in range(1, 9):
            stage = f"engine-point-{point}"
            feed = feeds[point]
            result = equilibrium.solve_equilibrium(model, engine_problem(feed), previous)
            phases, metrics, point_failures = phase_data(model, result, feed)
            engine_points.append({"point_id": point, "success": metrics["engine_success"], "requested_tolerance_met": metrics["engine_requested_tolerance_met"], "phase_count": metrics["engine_phase_count"], "solver_status": metrics["engine_solver_status"], "message": metrics["engine_message"]})
            residual_rows.extend({"point_id": point, **item} for item in metrics["row_evidence"])
            atomic_csv(output / "engine-residuals.csv", residual_rows, RESIDUAL_FIELDS)
            if point_failures:
                raise RuntimeError(f"point {point}: {'; '.join(point_failures)}")
            stage = f"idaes-point-{point}"
            mapping = solve_settler(feed, phases, solver)
            rows.extend(row_records(point, feed, phases, metrics, mapping))
            atomic_csv(output / "stage-results.csv", rows)
            completed.append(point)
            previous = result
        summary["maxima"] = {
            field: max(float(row[field]) for row in rows)
            for field in (
                "engine_material_balance_max_abs", "engine_phase_charge_max_abs",
                "engine_pressure_max_relative", "engine_classified_residual_max_abs",
                "idaes_component_balance_max_abs", "idaes_charge_residual_abs",
                "idaes_phase_fraction_abs", "idaes_phase_flow_closure_abs",
                "idaes_composition_max_abs",
            )
        }
        summary["minimum_engine_reduced_hessian_minimum"] = min(float(row["engine_reduced_hessian_minimum"]) for row in rows)
        summary.update({
            "status": "accepted",
            "typed_cause": None,
            "accepted_count": len(completed),
            "stopping_stage": "complete",
        })
        exit_code = 0
    except KeyboardInterrupt:
        summary.update({"status": "interrupted", "typed_cause": "keyboard-interrupt", "stopping_stage": stage})
        exit_code = 130
    except Exception as error:
        failures.append(f"{stage}: {type(error).__name__}: {error}")
        summary.update({"status": "failed", "typed_cause": type(error).__name__, "error": str(error), "stopping_stage": stage})
    finally:
        summary.update({"completed_point_ids": completed, "unattempted_point_ids": [p for p in range(1, 9) if p not in completed], "engine_points": engine_points, "engine_residuals": {"path": "engine-residuals.csv", "sha256": sha256(output / "engine-residuals.csv"), "row_count": len(residual_rows), "max_raw_abs": max((abs(float(row["raw_residual"])) for row in residual_rows), default=0.0), "max_classified_abs": max((abs(float(row["classified_residual"])) for row in residual_rows), default=0.0)}, "runtime_s": time.perf_counter() - started, "line_budget_physical_python_lines": sum(1 for _ in Path(__file__).read_text(encoding="utf-8").splitlines()), "failures": failures})
        atomic_json(output / "summary.json", summary)
    return exit_code


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True, type=Path)
    return run(parser.parse_args().output_dir.resolve())


if __name__ == "__main__":
    sys.exit(main())

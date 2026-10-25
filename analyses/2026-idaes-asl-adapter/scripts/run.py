from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import re
import selectors
import signal
import statistics
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from decimal import Decimal, localcontext
from pathlib import Path

import epcsaft
import numpy as np
from idaes.core import FlowsheetBlock, UnitModelBlockData, declare_process_block_class
from idaes.core.solvers import get_solver
from pyomo.contrib.pynumero.interfaces.ampl_nlp import AslNLP
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Expression,
    ExternalFunction,
    Objective,
    Suffix,
    Var,
    value,
)
from pyomo.environ import (
    units as pyunits,
)
from pyomo.network import Port

CAMPAIGN = Path(__file__).resolve().parents[1]
REPO = CAMPAIGN.parents[1]
sys.path.insert(0, str(REPO / "tools"))
from engine_identity import asl_library, loaded_identity
SOURCE_PACKET = REPO / "data/packets/gross-2001-propane/3"
PACKET = CAMPAIGN / "inputs/data/gross-2001-propane/3"
PARAMETERS = PACKET / "parameters/parameters.json"
MATERIALIZATION = PACKET.parent / "materialization-v3.json"


def packet_parameter_strings() -> tuple[str, str, str, str]:
    document = json.loads(PARAMETERS.read_text(encoding="utf-8"))
    component = document["components"][0]
    coefficients = {
        item["family"]: item["value"]["magnitude"]
        for item in component["coefficients"]
    }
    return tuple(
        str(coefficients[family])
        for family in ("segment_count", "segment_diameter", "dispersion_energy_over_k")
    ) + (str(component["fixed"]["molar_mass"]["value"]["magnitude"]),)


PARAMETER_STRINGS = packet_parameter_strings()
TEMPERATURE_SCALE = 400.0
DENSITY_SCALE = 100.0
PRESSURE_SCALE = 200000.0
FLOW_SCALE = 0.1
VOLUME_FLOW_SCALE = 0.001
TEMPERATURES = (390.0, 400.0, 410.0)
GRID_TEMPERATURES = (385.0, 392.5, 400.0, 407.5, 415.0)
DENSITIES = (20.0, 60.0, 90.0)
DIRECTIONS = (
    (1.0, 0.0),
    (0.0, 1.0),
    (1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0)),
    (1.0 / math.sqrt(2.0), -1.0 / math.sqrt(2.0)),
)
STEPS = (1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6)
GAS_CONSTANT = 8.31446261815324
TERM_GRACE_S = 2.0
KILL_GRACE_S = 0.1
TERMINATION_BUDGET_S = TERM_GRACE_S + KILL_GRACE_S
PRIMAL_NAMES = (
    "fs.unit.f_in",
    "fs.unit.f_out",
    "fs.unit.p_in",
    "fs.unit.p_out",
    "fs.unit.r",
    "fs.unit.q",
)
CONSTRAINT_NAMES = (
    "fs.unit.flow_balance",
    "fs.unit.pressure_balance",
    "fs.unit.volume_balance",
    "eos_pressure",
)
NL_IO_OPTIONS = {
    "symbolic_solver_labels": True,
    "linear_presolve": False,
    "scale_model": False,
    "skip_trivial_constraints": False,
}
SOLVER_OPTIONS = {
    "linear_solver": "mumps",
    "hessian_approximation": "exact",
    "point_perturbation_radius": 0,
    "bound_relax_factor": 0,
    "nlp_scaling_method": "none",
    "tol": 1e-9,
    "constr_viol_tol": 1e-8,
    "dual_inf_tol": 1e-7,
    "compl_inf_tol": 1e-7,
    "acceptable_iter": 0,
    "max_iter": 200,
}
DERIVATIVE_TEST_OPTIONS = {
    "derivative_test": "second-order",
    "derivative_test_perturbation": 1e-6,
    "derivative_test_tol": 1e-4,
}
REQUIRED_IPOPT_OPTIONS = tuple(sorted(set(SOLVER_OPTIONS) | set(DERIVATIVE_TEST_OPTIONS)))


@declare_process_block_class("NeutralThroughput")
class NeutralThroughputData(UnitModelBlockData):
    def build(self):
        super().build()
        flow_units = pyunits.mol / pyunits.s
        pressure_units = pyunits.Pa
        density_units = pyunits.mol / pyunits.m**3
        volume_flow_units = pyunits.m**3 / pyunits.s
        self.f_in = Var(bounds=(0.0, 1.0), initialize=0.03)
        self.f_out = Var(bounds=(0.0, 1.0), initialize=0.03)
        self.p_in = Var(bounds=(1.0 / PRESSURE_SCALE, 1.0), initialize=0.5)
        self.p_out = Var(bounds=(1.0 / PRESSURE_SCALE, 1.0), initialize=0.5)
        self.r = Var(bounds=(0.1, 1.0), initialize=0.6)
        self.q = Var(bounds=(0.1, 1.0), initialize=0.5)
        self.temperature_in = Expression(expr=400.0 * pyunits.K)
        self.temperature_out = Expression(expr=400.0 * pyunits.K)
        self.flow_mol_in = Expression(expr=FLOW_SCALE * flow_units * self.f_in)
        self.flow_mol_out = Expression(expr=FLOW_SCALE * flow_units * self.f_out)
        self.pressure_in = Expression(expr=PRESSURE_SCALE * pressure_units * self.p_in)
        self.pressure_out = Expression(expr=PRESSURE_SCALE * pressure_units * self.p_out)
        self.molar_density_out = Expression(expr=DENSITY_SCALE * density_units * self.r)
        self.volumetric_flow_out = Expression(expr=VOLUME_FLOW_SCALE * volume_flow_units * self.q)
        self.inlet = Port()
        self.inlet.add(self.flow_mol_in, "flow_mol")
        self.inlet.add(self.temperature_in, "temperature")
        self.inlet.add(self.pressure_in, "pressure")
        self.outlet = Port()
        self.outlet.add(self.flow_mol_out, "flow_mol")
        self.outlet.add(self.temperature_out, "temperature")
        self.outlet.add(self.pressure_out, "pressure")
        self.flow_balance = Constraint(expr=self.f_in - self.f_out == 0.0)
        self.pressure_balance = Constraint(expr=self.p_in - self.p_out == 0.0)
        self.volume_balance = Constraint(expr=self.f_out - self.r * self.q == 0.0)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ipopt_preflight(path: Path) -> dict[str, object]:
    version = subprocess.check_output((str(path), "-v"), text=True).strip()
    options = subprocess.check_output((str(path), "-="), text=True, stderr=subprocess.STDOUT)
    supported = {
        line.split(maxsplit=1)[0]
        for line in options.splitlines()
        if line and not line[0].isspace()
    }
    missing = [name for name in REQUIRED_IPOPT_OPTIONS if name not in supported]
    if missing:
        raise RuntimeError(f"selected Ipopt is missing required options: {missing}")
    return {
        "sha256": sha256(path),
        "version": version,
        "asl_version": version.split("ASL(", 1)[1].split(")", 1)[0],
        "required_options": REQUIRED_IPOPT_OPTIONS,
        "supported_options": {name: name in supported for name in REQUIRED_IPOPT_OPTIONS},
        "all_required_options_supported": True,
    }


def materialized_packet_identity() -> dict[str, object]:
    record = json.loads(MATERIALIZATION.read_text(encoding="utf-8"))
    expected_destination = str(PACKET.relative_to(REPO))
    if record.get("destination") != expected_destination:
        raise RuntimeError("materialized packet destination changed")
    for item in record.get("files", []):
        path = PACKET / str(item["path"])
        if not path.is_file() or sha256(path) != item.get("sha256"):
            raise RuntimeError(f"materialized packet file changed: {path}")
    fingerprint = record.get("packet_fingerprint")
    if not isinstance(fingerprint, str) or len(fingerprint) != 64:
        raise RuntimeError("materialized packet fingerprint is unavailable")
    return {
        "packet_path": str(SOURCE_PACKET.relative_to(REPO)),
        "materialized_path": expected_destination,
        "packet_fingerprint": fingerprint,
        "materialization_sha256": sha256(MATERIALIZATION),
    }


def installed_identity() -> tuple[Path | None, str | None, Path]:
    identity = loaded_identity()
    return Path(identity["wheel"]) if identity["wheel"] else None, identity["wheel_sha256"], asl_library(identity)


def shared_library_identity(path: Path) -> dict[str, object]:
    dynamic = subprocess.check_output(("readelf", "-d", str(path)), text=True)
    exports = subprocess.check_output(("nm", "-D", "--defined-only", str(path)), text=True)
    needed = [line.split("[", 1)[1].split("]", 1)[0] for line in dynamic.splitlines() if "(NEEDED)" in line]
    runpaths = [line.split("[", 1)[1].split("]", 1)[0] for line in dynamic.splitlines() if "(RUNPATH)" in line or "(RPATH)" in line]
    symbols = [line.rsplit(None, 1)[-1] for line in exports.splitlines() if line.strip()]
    if symbols != ["funcadd_ASL"]:
        raise RuntimeError(f"unexpected ASL exports: {symbols}")
    # A repaired wheel's ASL finds its bundled Ipopt in epcsaft.libs/ (#576); no other runpath is allowed.
    if runpaths not in ([], ["$ORIGIN/../epcsaft.libs"]):
        raise RuntimeError(f"ASL library has an install-relative runpath: {runpaths}")
    return {"exports": symbols, "needed": needed, "runpaths": runpaths}


def campaign_identity(ipopt: Path) -> dict[str, object]:
    if not ipopt.is_file() or not os.access(ipopt, os.X_OK):
        raise RuntimeError(f"selected Ipopt executable is unavailable: {ipopt}")
    wheel, wheel_hash, library = installed_identity()
    ipopt_identity = ipopt_preflight(ipopt)
    return {
        "source_commit": subprocess.check_output(
            ("git", "rev-parse", "HEAD"), cwd=REPO, text=True
        ).strip(),
        "source_dirty": bool(subprocess.check_output(
            ("git", "status", "--porcelain"), cwd=REPO, text=True
        ).strip()),
        "wheel": str(wheel) if wheel else None,
        "wheel_sha256": wheel_hash,
        "prepared_environment": loaded_identity(),
        "asl_library": str(library),
        "asl_library_sha256": sha256(library),
        "asl_library_identity": shared_library_identity(library),
        "adapter_source_sha256": sha256(REPO / "engine/native/bindings/asl_adapter.cpp"),
        "ipopt": str(ipopt),
        "ipopt_sha256": ipopt_identity["sha256"],
        "ipopt_version": ipopt_identity["version"],
        "asl_version": ipopt_identity["asl_version"],
        "ipopt_preflight": ipopt_identity,
        "idaes_version": importlib.metadata.version("idaes-pse"),
        "pyomo_version": importlib.metadata.version("pyomo"),
        "parameter_strings": PARAMETER_STRINGS,
        **materialized_packet_identity(),
    }


def native_model():
    parameters = epcsaft.Parameters.from_json(PARAMETERS)
    return epcsaft._core.HelmholtzModel(parameters._native)


def asl_functions(library: Path):
    return (
        ExternalFunction(
            library=str(library),
            function="epcsaft_neutral_ar",
            units=pyunits.dimensionless,
            arg_units=[pyunits.K, pyunits.mol / pyunits.m**3, None, None, None, None],
        ),
        ExternalFunction(
            library=str(library),
            function="epcsaft_neutral_pressure",
            units=pyunits.Pa,
            arg_units=[pyunits.K, pyunits.mol / pyunits.m**3, None, None, None, None],
        ),
    )


def external_fgh(function, temperature: float, density: float, fgh: int, fixed=None):
    result = function.evaluate_fgh(
        [temperature, density, *PARAMETER_STRINGS], fgh=fgh, fixed=fixed
    )
    value_, gradient, hessian = result
    if fgh == 0:
        return float(value_)
    gradient_ = np.array([gradient[0], gradient[1]], dtype=float)
    if fgh == 1:
        return float(value_), gradient_
    return float(value_), gradient_, np.array(
        [[hessian[0], hessian[1]], [hessian[1], hessian[2]]], dtype=float
    )


def native_fgh(model, property_, temperature: float, density: float, order):
    result = model.neutral_property(temperature, density, property_, order)
    gradient = np.array(result.gradient, dtype=float) if result.gradient else None
    hessian = None
    if result.hessian:
        hessian = np.array(result.hessian, dtype=float).reshape(2, 2)
    return float(result.value), gradient, hessian


def scaled_gradient(gradient: np.ndarray, output_scale: float) -> np.ndarray:
    return gradient * np.array([TEMPERATURE_SCALE, DENSITY_SCALE]) / output_scale


def scaled_hessian(hessian: np.ndarray, output_scale: float) -> np.ndarray:
    scales = np.array([TEMPERATURE_SCALE, DENSITY_SCALE])
    return hessian * scales[:, None] * scales[None, :] / output_scale


def tolerance(reference, result, absolute: float, relative: float) -> float:
    reference_norm = float(np.linalg.norm(np.atleast_1d(reference), ord=np.inf))
    result_norm = float(np.linalg.norm(np.atleast_1d(result), ord=np.inf))
    return absolute + relative * max(reference_norm, result_norm)


def ladder_summary(
    estimates: list[np.ndarray], expected: np.ndarray, steps: tuple[float, ...],
    absolute: float, relative: float,
) -> dict[str, object]:
    errors = [float(np.linalg.norm(item - expected, ord=np.inf)) for item in estimates]
    targets = [tolerance(expected, item, absolute, relative) for item in estimates]
    adjacent = next(
        (index for index in range(len(steps) - 1)
         if errors[index] <= targets[index] and errors[index + 1] <= targets[index + 1]),
        None,
    )
    supports = []
    for index in range(len(steps) - 2):
        floor_target = 0.1 * min(targets[index:index + 3])
        floor = max(
            errors[index:index + 3]
            + [
                float(np.linalg.norm(estimates[index] - estimates[index + 1], ord=np.inf)),
                float(np.linalg.norm(estimates[index + 1] - estimates[index + 2], ord=np.inf)),
            ]
        )
        orders = []
        for left, right in ((index, index + 1), (index + 1, index + 2)):
            if errors[left] > 0.0 and errors[right] > 0.0:
                orders.append(math.log(errors[left] / errors[right])
                              / math.log(steps[left] / steps[right]))
        scale = max(1.0, float(np.linalg.norm(expected, ord=np.inf)))
        resolved = len(orders) == 2 and all(1.5 <= order <= 2.5 for order in orders)
        resolved = resolved and all(error > 100.0 * np.finfo(float).eps * scale
                                    for error in errors[index:index + 3])
        exact = (
            all(error <= floor_target for error in errors[index:index + 3])
            and float(np.linalg.norm(estimates[index] - estimates[index + 1], ord=np.inf))
                <= floor_target
            and float(np.linalg.norm(estimates[index + 1] - estimates[index + 2], ord=np.inf))
                <= floor_target
        )
        if floor <= floor_target and (resolved or exact):
            supports.append((index, floor, floor_target, "resolved" if resolved else "roundoff"))
    if adjacent is None or not supports:
        raise RuntimeError(
            f"derivative ladder failed: adjacent={adjacent}, floors={supports}, errors={errors}"
        )
    support, floor, floor_target, mode = supports[0]
    return {
        "errors": errors,
        "targets": targets,
        "adjacent_support": adjacent,
        "floor": floor,
        "floor_target": floor_target,
        "floor_support": support,
        "mode": mode,
    }


def unit_pressure_chain_checks(
    function, native_model, library: Path, rows: list[dict[str, object]]
) -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="epcsaft-unit-chain-") as temporary:
        check = ConcreteModel()
        check.tau = Var(bounds=(0.9, 1.1), initialize=1.0)
        check.r = Var(bounds=(0.05, 1.2), initialize=0.6)
        pressure = function(
            400.0 * check.tau * pyunits.K,
            100.0 * check.r * pyunits.mol / pyunits.m**3,
            *PARAMETER_STRINGS,
        )
        check.pa_row = Constraint(
            expr=(pressure - 100000.0 * pyunits.Pa) / (200000.0 * pyunits.Pa) == 0.0
        )
        check.bar_row = Constraint(
            expr=(pyunits.convert(pressure, to_units=pyunits.bar) - pyunits.bar)
            / (2.0 * pyunits.bar) == 0.0
        )
        check.objective = Objective(expr=0.0)
        path = Path(temporary) / "unit-chain.nl"
        write_test_nl(check, path)
        os.environ["AMPLFUNC"] = str(library)
        nlp = AslNLP(str(path))
        columns, row_names = named_file_rows(path)
        if nlp.n_primals() != 2 or nlp.n_constraints() != 2:
            raise RuntimeError("unit-chain check NL does not have two variables and two rows")
        if set(columns) != {"tau", "r"} or set(row_names[:2]) != {"pa_row", "bar_row"}:
            raise RuntimeError(f"unexpected unit-chain names: columns={columns}, rows={row_names}")
        positions = {name: index for index, name in enumerate(columns)}
        constraint_positions = {name: index for index, name in enumerate(row_names[:2])}
        maximums = {"value": 0.0, "gradient": 0.0, "hessian": 0.0}
        for temperature in TEMPERATURES:
            for density in DENSITIES:
                native_value, native_gradient, native_hessian = native_fgh(
                    native_model, epcsaft._core.NeutralProperty.Pressure,
                    temperature, density, epcsaft._core.DerivativeOrder.Hessian,
                )
                expected_value = native_value / PRESSURE_SCALE - 0.5
                expected_gradient = np.array([400.0 * native_gradient[0], 100.0 * native_gradient[1]]) / PRESSURE_SCALE
                scales = np.array([400.0, 100.0])
                expected_hessian = native_hessian * scales[:, None] * scales[None, :] / PRESSURE_SCALE
                nlp.set_primals(np.array([temperature / 400.0, density / 100.0]))
                lower = np.asarray(nlp.constraints_lb(), dtype=float)
                upper = np.asarray(nlp.constraints_ub(), dtype=float)
                if not np.allclose(lower, upper, atol=0.0, rtol=0.0):
                    raise RuntimeError("unit-chain rows are not equalities")
                residuals = np.asarray(nlp.evaluate_constraints(), dtype=float) - lower
                jacobian = np.asarray(nlp.evaluate_jacobian().toarray(), dtype=float)
                value_errors = {}
                gradient_errors = {}
                hessian_errors = {}
                value_targets = {}
                gradient_targets = {}
                hessian_targets = {}
                for row_name, row_index in constraint_positions.items():
                    value_errors[row_name] = abs(residuals[row_index] - expected_value)
                    jacobian_row = jacobian[row_index, [positions["tau"], positions["r"]]]
                    gradient_errors[row_name] = float(np.linalg.norm(
                        jacobian_row - expected_gradient,
                        ord=np.inf,
                    ))
                    value_targets[row_name] = tolerance(
                        np.array([expected_value]), np.array([residuals[row_index]]), 1e-11, 1e-9
                    )
                    gradient_targets[row_name] = tolerance(
                        expected_gradient, jacobian_row, 1e-11, 1e-9
                    )
                    duals = np.zeros(nlp.n_constraints())
                    duals[row_index] = 1.0
                    nlp.set_obj_factor(0.0)
                    nlp.set_duals(duals)
                    hessian = np.asarray(nlp.evaluate_hessian_lag().toarray(), dtype=float)
                    hessian_errors[row_name] = float(np.linalg.norm(hessian - expected_hessian, ord=np.inf))
                    hessian_targets[row_name] = tolerance(
                        expected_hessian, hessian, 1e-11, 1e-9
                    )
                errors = {
                    "value": max(value_errors.values()),
                    "gradient": max(gradient_errors.values()),
                    "hessian": max(hessian_errors.values()),
                }
                targets = {
                    "value": max(value_targets.values()),
                    "gradient": max(gradient_targets.values()),
                    "hessian": max(hessian_targets.values()),
                }
                if not all(math.isfinite(error) for error in errors.values()):
                    raise RuntimeError(f"nonfinite unit-chain check at T={temperature} rho={density}")
                if errors["value"] > targets["value"] or errors["gradient"] > targets["gradient"]:
                    raise RuntimeError(
                        f"unit-chain value/gradient mismatch at T={temperature} rho={density}: "
                        f"errors={errors}, targets={targets}"
                    )
                if errors["hessian"] > max(hessian_targets.values()):
                    raise RuntimeError(
                        f"unit-chain Hessian mismatch at T={temperature} rho={density}: "
                        f"error={errors['hessian']}, target={max(hessian_targets.values())}"
                    )
                maximums = {key: max(maximums[key], errors[key]) for key in maximums}
                rows.append({
                    "kind": "unit_pressure_chain", "temperature_k": temperature,
                    "density_mol_per_m3": density,
                    "pa_value_error": value_errors["pa_row"],
                    "bar_value_error": value_errors["bar_row"],
                    "gradient_error": errors["gradient"], "hessian_error": errors["hessian"],
                    "gradient_target": targets["gradient"], "hessian_target": targets["hessian"],
                    "passed": True,
                })
        return {
            "state_count": len(TEMPERATURES) * len(DENSITIES),
            "max_value_error": maximums["value"],
            "max_gradient_error": maximums["gradient"],
            "max_hessian_error": maximums["hessian"],
            "passed": True,
        }


def derivative_checks(
    model, functions, library: Path, rows: list[dict[str, object]],
    output_path: Path | None = None,
) -> dict[str, object]:
    properties = (
        ("a_res", epcsaft._core.NeutralProperty.ResidualHelmholtz, functions[0], 1.0),
        ("pressure", epcsaft._core.NeutralProperty.Pressure, functions[1], PRESSURE_SCALE),
    )
    max_transport = 0.0
    max_transport_gradient = 0.0
    max_transport_hessian = 0.0
    max_gradient_error = 0.0
    max_hessian_error = 0.0
    floor_max = 0.0
    ladders = 0
    for name, property_, function, output_scale in properties:
        for temperature in TEMPERATURES:
            for density in DENSITIES:
                native_value, native_gradient, native_hessian = native_fgh(
                    model, property_, temperature, density, epcsaft._core.DerivativeOrder.Hessian
                )
                asl_value, asl_gradient, asl_hessian = external_fgh(function, temperature, density, 2)
                native_scaled_value = native_value / output_scale
                asl_scaled_value = asl_value / output_scale
                native_scaled_gradient = scaled_gradient(native_gradient, output_scale)
                asl_scaled_gradient = scaled_gradient(asl_gradient, output_scale)
                native_scaled_hessian = scaled_hessian(native_hessian, output_scale)
                asl_scaled_hessian = scaled_hessian(asl_hessian, output_scale)
                transport_errors = {
                    "value": abs(native_scaled_value - asl_scaled_value),
                    "gradient": float(np.linalg.norm(
                        native_scaled_gradient - asl_scaled_gradient, ord=np.inf
                    )),
                    "hessian": float(np.linalg.norm(
                        native_scaled_hessian - asl_scaled_hessian, ord=np.inf
                    )),
                }
                transport_targets = {
                    name: tolerance(
                        native, asl, 1e-11, 1e-9
                    ) for name, native, asl in (
                        ("value", np.array([native_scaled_value]), np.array([asl_scaled_value])),
                        ("gradient", native_scaled_gradient, asl_scaled_gradient),
                        ("hessian", native_scaled_hessian, asl_scaled_hessian),
                    )
                }
                if (
                    not all(math.isfinite(error) for error in transport_errors.values())
                    or any(transport_errors[name] > transport_targets[name] for name in transport_errors)
                ):
                    raise RuntimeError(
                        f"native/ASL transport mismatch at {name} T={temperature} rho={density}: "
                        f"errors={transport_errors}, targets={transport_targets}"
                    )
                max_transport = max(max_transport, transport_errors["value"])
                max_transport_gradient = max(max_transport_gradient, transport_errors["gradient"])
                max_transport_hessian = max(max_transport_hessian, transport_errors["hessian"])
                rows.append({
                    "kind": "transport", "property": name,
                    "temperature_k": temperature, "density_mol_per_m3": density,
                    "value_error": transport_errors["value"],
                    "value_target": transport_targets["value"],
                    "gradient_error": transport_errors["gradient"],
                    "gradient_target": transport_targets["gradient"],
                    "hessian_error": transport_errors["hessian"],
                    "hessian_target": transport_targets["hessian"], "passed": True,
                })
                fixed_t = external_fgh(
                    function, temperature, density, 2,
                    fixed=[True, False, False, False, False, False],
                )
                fixed_rho = external_fgh(
                    function, temperature, density, 2,
                    fixed=[False, True, False, False, False, False],
                )
                if not np.allclose(fixed_t[1][1], asl_gradient[1], atol=1e-12, rtol=1e-12):
                    raise RuntimeError("fixed-temperature ASL derivative mask changed rho gradient")
                if not np.allclose(fixed_rho[1][0], asl_gradient[0], atol=1e-12, rtol=1e-12):
                    raise RuntimeError("fixed-density ASL derivative mask changed T gradient")
                if not np.isclose(fixed_t[2][1, 1], asl_hessian[1, 1], atol=1e-12, rtol=1e-12):
                    raise RuntimeError("fixed-temperature ASL derivative mask changed rho curvature")
                if not np.isclose(fixed_rho[2][0, 0], asl_hessian[0, 0], atol=1e-12, rtol=1e-12):
                    raise RuntimeError("fixed-density ASL derivative mask changed T curvature")
                for direction in DIRECTIONS:
                    direction_vector = np.array(direction)
                    gradient_errors: list[float] = []
                    hessian_errors: list[float] = []
                    gradient_estimates: list[np.ndarray] = []
                    hessian_estimates: list[np.ndarray] = []
                    for step in STEPS:
                        dt, drho = direction_vector * step * np.array([TEMPERATURE_SCALE, DENSITY_SCALE])
                        plus_value = external_fgh(function, temperature + dt, density + drho, 0) / output_scale
                        minus_value = external_fgh(function, temperature - dt, density - drho, 0) / output_scale
                        centered_gradient = (plus_value - minus_value) / (2.0 * step)
                        expected_gradient = float(np.dot(scaled_gradient(asl_gradient, output_scale), direction_vector))
                        gradient_error = abs(centered_gradient - expected_gradient)
                        plus_gradient = external_fgh(function, temperature + dt, density + drho, 1)[1]
                        minus_gradient = external_fgh(function, temperature - dt, density - drho, 1)[1]
                        centered_hessian = (
                            scaled_gradient(plus_gradient, output_scale)
                            - scaled_gradient(minus_gradient, output_scale)
                        ) / (2.0 * step)
                        expected_hessian = scaled_hessian(asl_hessian, output_scale) @ direction_vector
                        hessian_error = float(np.linalg.norm(centered_hessian - expected_hessian, ord=np.inf))
                        gradient_errors.append(gradient_error)
                        hessian_errors.append(hessian_error)
                        gradient_estimates.append(np.array([centered_gradient]))
                        hessian_estimates.append(centered_hessian)
                        gradient_target = tolerance(
                            np.array([expected_gradient]), np.array([centered_gradient]), 1e-8, 1e-6
                        )
                        hessian_target = tolerance(expected_hessian, centered_hessian, 1e-7, 1e-4)
                        passed = gradient_error <= gradient_target and hessian_error <= hessian_target
                        rows.append({
                            "property": name, "temperature_k": temperature, "density_mol_per_m3": density,
                            "direction": f"{direction[0]},{direction[1]}", "step": step,
                            "gradient_error": gradient_error, "gradient_target": gradient_target,
                            "hessian_error": hessian_error, "hessian_target": hessian_target,
                            "passed": passed,
                        })
                    max_gradient_error = max(max_gradient_error, max(gradient_errors))
                    max_hessian_error = max(max_hessian_error, max(hessian_errors))
                    gradient_summary = ladder_summary(
                        gradient_estimates, np.array([expected_gradient]), STEPS, 1e-8, 1e-6
                    )
                    hessian_summary = ladder_summary(
                        hessian_estimates, expected_hessian, STEPS, 1e-7, 1e-4
                    )
                    floor_max = max(floor_max, gradient_summary["floor"], hessian_summary["floor"])
                    ladders += 1
                    for row in rows[-len(STEPS):]:
                        row.update({
                            "gradient_floor": gradient_summary["floor"],
                            "gradient_floor_target": gradient_summary["floor_target"],
                            "gradient_floor_support": gradient_summary["floor_support"],
                            "gradient_mode": gradient_summary["mode"],
                            "hessian_floor": hessian_summary["floor"],
                            "hessian_floor_target": hessian_summary["floor_target"],
                            "hessian_floor_support": hessian_summary["floor_support"],
                            "hessian_mode": hessian_summary["mode"],
                        })
                    if output_path is not None:
                        write_rows(output_path, rows)
    unit_summary = unit_pressure_chain_checks(functions[1], model, library, rows)
    if output_path is not None:
        write_rows(output_path, rows)
    return {
        "transport_max_scaled_value_error": max_transport,
        "transport_max_scaled_gradient_error": max_transport_gradient,
        "transport_max_scaled_hessian_error": max_transport_hessian,
        "unit_pressure_chain": unit_summary,
        "transport_passed": True,
        "gradient_max_scaled_error": max_gradient_error,
        "hessian_max_scaled_error": max_hessian_error,
        "derivative_ladders_passed": ladders,
        "maximum_estimated_floor": floor_max,
    }


def throughput_model(
    library: Path, pressure_cap: float, start: tuple[float, float],
    parameter_strings: tuple[str, str, str, str] = PARAMETER_STRINGS,
):
    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)
    model.fs.unit = globals()["NeutralThroughput"](dynamic=False, has_holdup=False)
    unit = model.fs.unit
    unit.p_in.setub(pressure_cap / PRESSURE_SCALE)
    unit.p_out.setub(pressure_cap / PRESSURE_SCALE)
    density, volume_flow = start
    unit.r.set_value(density / DENSITY_SCALE)
    unit.q.set_value(volume_flow / VOLUME_FLOW_SCALE)
    unit.f_in.set_value(density * volume_flow / FLOW_SCALE)
    unit.f_out.set_value(density * volume_flow / FLOW_SCALE)
    unit.p_in.set_value(0.5)
    unit.p_out.set_value(0.5)
    model.pressure = ExternalFunction(
        library=str(library), function="epcsaft_neutral_pressure", units=pyunits.Pa,
        arg_units=[pyunits.K, pyunits.mol / pyunits.m**3, None, None, None, None],
    )
    model.eos_pressure = Constraint(
        expr=(model.pressure(
            TEMPERATURE_SCALE * pyunits.K,
            DENSITY_SCALE * unit.r * pyunits.mol / pyunits.m**3,
            *parameter_strings,
        ) - PRESSURE_SCALE * unit.p_out * pyunits.Pa) / (PRESSURE_SCALE * pyunits.Pa) == 0.0
    )
    model.objective = Objective(
        expr=0.5 * unit.q ** 2 - unit.f_out
    )
    return model


def process_solution(model, result) -> dict[str, object]:
    unit = model.fs.unit
    return {
        "termination": str(result.solver.termination_condition),
        "status": str(result.solver.status),
        "rho": DENSITY_SCALE * value(unit.r),
        "Q": VOLUME_FLOW_SCALE * value(unit.q),
        "F": FLOW_SCALE * value(unit.f_out),
        "P": PRESSURE_SCALE * value(unit.p_out),
        "flow_balance": value(unit.flow_balance.body),
        "volume_balance": value(unit.volume_balance.body),
        "pressure_balance": value(unit.pressure_balance.body),
        "eos_residual": value(model.eos_pressure.body),
        "objective": value(model.objective),
    }


def write_test_nl(model, path: Path) -> None:
    model.write(str(path), format="nl", io_options=NL_IO_OPTIONS)


def named_file_rows(path: Path) -> tuple[list[str], list[str]]:
    columns = path.with_suffix(".col").read_text(encoding="utf-8").splitlines()
    rows = path.with_suffix(".row").read_text(encoding="utf-8").splitlines()
    return columns, rows


def kkt_diagnostics(model, nl_path: Path, library: Path) -> dict[str, object]:
    os.environ["AMPLFUNC"] = str(library)
    nlp = AslNLP(str(nl_path))
    columns, rows = named_file_rows(nl_path)
    if set(columns) != set(PRIMAL_NAMES) or set(rows[:len(CONSTRAINT_NAMES)]) != set(CONSTRAINT_NAMES):
        raise RuntimeError(f"unexpected retained NL names: columns={columns}, rows={rows}")
    values = {
        "fs.unit.f_in": value(model.fs.unit.f_in),
        "fs.unit.f_out": value(model.fs.unit.f_out),
        "fs.unit.p_in": value(model.fs.unit.p_in),
        "fs.unit.p_out": value(model.fs.unit.p_out),
        "fs.unit.r": value(model.fs.unit.r),
        "fs.unit.q": value(model.fs.unit.q),
    }
    primals = np.array([values[name] for name in columns], dtype=float)
    constraint_duals = {
        "fs.unit.flow_balance": float(model.dual[model.fs.unit.flow_balance]),
        "fs.unit.pressure_balance": float(model.dual[model.fs.unit.pressure_balance]),
        "fs.unit.volume_balance": float(model.dual[model.fs.unit.volume_balance]),
        "eos_pressure": float(model.dual[model.eos_pressure]),
    }
    duals = -np.array([constraint_duals[name] for name in rows[:len(CONSTRAINT_NAMES)]], dtype=float)
    lower = np.array(nlp.primals_lb(), dtype=float)
    upper = np.array(nlp.primals_ub(), dtype=float)
    z_lower = np.array([
        float(model.ipopt_zL_out[getattr(model.fs.unit, name.rsplit(".", 1)[-1])])
        if model.ipopt_zL_out.get(getattr(model.fs.unit, name.rsplit(".", 1)[-1])) is not None else 0.0
        for name in columns
    ])
    z_upper = -np.array([
        float(model.ipopt_zU_out[getattr(model.fs.unit, name.rsplit(".", 1)[-1])])
        if model.ipopt_zU_out.get(getattr(model.fs.unit, name.rsplit(".", 1)[-1])) is not None else 0.0
        for name in columns
    ])
    nlp.set_primals(primals)
    nlp.set_obj_factor(1.0)
    nlp.set_duals(duals)
    stationarity = nlp.evaluate_grad_objective() + nlp.evaluate_jacobian().T @ duals
    stationarity = stationarity - z_lower + z_upper
    constraints = np.asarray(nlp.evaluate_constraints(), dtype=float)
    lower_violation = np.maximum(lower - primals, 0.0)
    upper_violation = np.maximum(primals - upper, 0.0)
    complementarity = np.concatenate((z_lower * (primals - lower), z_upper * (upper - primals)))
    result = {
        "stationarity_inf": float(np.linalg.norm(stationarity, ord=np.inf)),
        "primal_residual_inf": float(np.linalg.norm(constraints, ord=np.inf)),
        "bound_violation_inf": float(max(np.max(lower_violation), np.max(upper_violation))),
        "dual_sign_violation": float(max(np.max(-z_lower), np.max(-z_upper), 0.0)),
        "complementarity_inf": float(np.max(np.abs(complementarity))),
    }
    result["kkt_passed"] = all(
        result[name] <= limit for name, limit in {
            "stationarity_inf": 1e-7,
            "primal_residual_inf": 1e-8,
            "bound_violation_inf": 1e-8,
            "dual_sign_violation": 1e-8,
            "complementarity_inf": 1e-7,
        }.items()
    )
    return result


def prepare_process_nl(
    library: Path, pressure_cap: float, start: tuple[float, float],
    output_dir: Path, case_name: str,
) -> Path:
    case_dir = output_dir / "solver" / case_name
    case_dir.mkdir(parents=True, exist_ok=False)
    path = case_dir / f"{case_name}.nl"
    write_test_nl(throughput_model(library, pressure_cap, start), path)
    return path


def load_file_solution(model, result, nl_path: Path) -> None:
    columns, rows = named_file_rows(nl_path)
    solution = result.solution[0]
    for index, name in enumerate(columns):
        variable = model.find_component(name)
        values = solution.variable[f"v{index}"]
        variable.set_value(float(values["Value"]))
        for suffix_name in ("ipopt_zL_out", "ipopt_zU_out"):
            getattr(model, suffix_name)[variable] = float(values.get(suffix_name, 0.0))
    for index, name in enumerate(rows[:len(CONSTRAINT_NAMES)]):
        model.dual[model.find_component(name)] = float(
            solution.constraint[f"c{index}"].get("Dual", 0.0)
        )


def solve_process(
    library: Path, pressure_cap: float, start: tuple[float, float],
    ipopt: Path, output_dir: Path, nl_path: Path, case_name: str = "process",
    derivative_test: bool = False, entry_probe: bool = False,
) -> dict[str, object]:
    model = throughput_model(library, pressure_cap, start)
    model.dual = Suffix(direction=Suffix.IMPORT)
    model.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
    model.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
    if not nl_path.is_file():
        raise RuntimeError(f"prewritten solver NL is unavailable: {nl_path}")
    case_dir = output_dir / "solver" / case_name
    ipopt_log = case_dir / "ipopt.log"
    os.environ["AMPLFUNC"] = str(library)
    options = {
        **SOLVER_OPTIONS,
        "print_level": 5 if entry_probe else 0,
        "output_file": str(ipopt_log),
        "file_print_level": 5,
    }
    if derivative_test:
        options.update(DERIVATIVE_TEST_OPTIONS)
    solver = get_solver(solver_options=options)
    solver.set_executable(str(ipopt))
    started = time.perf_counter()
    nl_hash = sha256(nl_path)
    try:
        result = solver.solve(
            str(nl_path), tee=entry_probe, load_solutions=False,
            suffixes=["dual", "ipopt_zL_out", "ipopt_zU_out"],
        )
        elapsed = time.perf_counter() - started
        load_file_solution(model, result, nl_path)
        record = process_solution(model, result)
    except Exception as error:
        elapsed = time.perf_counter() - started
        record = {
            "termination": "error", "status": "failed", "error": repr(error),
        }
    record["elapsed_s"] = elapsed
    text = ipopt_log.read_text(errors="replace") if ipopt_log.is_file() else ""
    record["derivative_check"] = (
        "passed"
        if "No errors detected by derivative checker." in text
        else "failed"
        if derivative_test
        else "not-requested"
    )
    if derivative_test and record["derivative_check"] != "passed":
        record["status"] = "failed"
        record["error"] = "Ipopt derivative checker did not report a clean result"
    record["hessian_evaluations"] = next(
        (int(line.rsplit(None, 1)[-1]) for line in text.splitlines()
         if "Number of Lagrangian Hessian evaluations" in line), None
    )
    record["solver_command"] = str(getattr(getattr(solver, "_command", None), "cmd", ""))
    if sha256(nl_path) != nl_hash:
        raise RuntimeError("Ipopt changed the prewritten NL file")
    if record.get("termination") == "optimal":
        record.update(kkt_diagnostics(model, nl_path, library))
    record.update({
        "solver_nl": str(nl_path),
        "solver_row": str(nl_path.with_suffix(".row")),
        "solver_col": str(nl_path.with_suffix(".col")),
        "solver_log": str(ipopt_log),
        "solver_nl_sha256": nl_hash,
        "writer_options": NL_IO_OPTIONS,
    })
    for solution in case_dir.glob("*.sol"):
        solution.unlink()
    for duplicate_log in case_dir.glob("*_ipopt.log"):
        duplicate_log.unlink()
    return record


def lagrangian_hessian_checks(
    library: Path, nl_path: Path, rows: list[dict[str, object]],
    output_path: Path | None = None,
) -> dict[str, object]:
    os.environ["AMPLFUNC"] = str(library)
    nlp = AslNLP(str(nl_path))
    columns, row_names = named_file_rows(nl_path)
    if set(columns) != set(PRIMAL_NAMES) or set(row_names[:4]) != set(CONSTRAINT_NAMES):
        raise RuntimeError(f"unexpected whole-NLP names: columns={columns}, rows={row_names}")
    if nlp.n_primals() != 6 or nlp.n_constraints() != 4:
        raise RuntimeError("whole-NLP model does not have six variables and four rows")

    positions = {name: index for index, name in enumerate(columns)}
    starts = ((20.0, 2e-4), (60.0, 5e-4), (90.0, 9e-4))

    def start_point(start: tuple[float, float]) -> np.ndarray:
        rho, volume_flow = start
        point = np.array(nlp.init_primals(), dtype=float)
        point[positions["fs.unit.f_in"]] = rho * volume_flow / FLOW_SCALE
        point[positions["fs.unit.f_out"]] = rho * volume_flow / FLOW_SCALE
        point[positions["fs.unit.p_in"]] = 0.5
        point[positions["fs.unit.p_out"]] = 0.5
        point[positions["fs.unit.r"]] = rho / DENSITY_SCALE
        point[positions["fs.unit.q"]] = volume_flow / VOLUME_FLOW_SCALE
        return point

    x = start_point(starts[1])
    invalid = x.copy()
    invalid[positions["fs.unit.r"]] = -0.01
    nlp.set_primals(invalid)
    try:
        invalid_rejected = not np.all(np.isfinite(nlp.evaluate_constraints()))
    except (RuntimeError, TypeError, ValueError):
        invalid_rejected = True
    if not invalid_rejected:
        raise RuntimeError("whole-NLP invalid-domain recovery failed")
    with tempfile.TemporaryDirectory(prefix="epcsaft-changed-inputs-") as changed_temp:
        changed_model = throughput_model(
            library, PRESSURE_SCALE, (60.0, 5e-4), ("2.003", *PARAMETER_STRINGS[1:])
        )
        changed_path = Path(changed_temp) / "changed-inputs.nl"
        write_test_nl(changed_model, changed_path)
        os.environ["AMPLFUNC"] = str(library)
        changed_nlp = AslNLP(str(changed_path))
        changed_nlp.set_primals(changed_nlp.init_primals())
        try:
            changed_input_rejected = not np.all(np.isfinite(changed_nlp.evaluate_constraints()))
        except (RuntimeError, TypeError, ValueError):
            changed_input_rejected = True
    nlp.set_primals(x)
    valid_after_failure = bool(np.all(np.isfinite(nlp.evaluate_constraints())))
    if not changed_input_rejected:
        raise RuntimeError("whole-NLP immutable parameter change was accepted")
    if not valid_after_failure:
        raise RuntimeError("whole-NLP valid recovery after immutable-input failure failed")

    directions = [np.eye(nlp.n_primals())[index] for index in range(nlp.n_primals())]
    for sign in (1.0, -1.0):
        direction = np.zeros(nlp.n_primals())
        direction[positions["fs.unit.r"]] = 1.0 / math.sqrt(2.0)
        direction[positions["fs.unit.q"]] = sign / math.sqrt(2.0)
        directions.append(direction)

    def gradient(point: np.ndarray, objective_factor: float, duals: np.ndarray) -> np.ndarray:
        nlp.set_primals(point)
        nlp.set_obj_factor(objective_factor)
        nlp.set_duals(duals)
        return objective_factor * nlp.evaluate_grad_objective() + nlp.evaluate_jacobian().T @ duals

    multipliers = (
        (0.0, 0.0, 0.0, 0.0),
        (1.0, -2.0, 3.0, -4.0),
        (-1.0, 2.0, -3.0, 4.0),
    )
    maximum = 0.0
    maximum_floor = 0.0
    cases = 0
    matrix_cases = 0
    for start in starts:
        x = start_point(start)
        native = native_model()
        native_hessian = native_fgh(
            native, epcsaft._core.NeutralProperty.Pressure, TEMPERATURE_SCALE,
            start[0], epcsaft._core.DerivativeOrder.Hessian,
        )[2]
        for objective_factor in (0.0, 0.5, 1.0):
            for multiplier in multipliers:
                multiplier_by_name = dict(zip(CONSTRAINT_NAMES, multiplier))
                duals = np.array([multiplier_by_name[name] for name in row_names[:4]], dtype=float)
                nlp.set_primals(x)
                nlp.set_obj_factor(objective_factor)
                nlp.set_duals(duals)
                expected_matrix = np.zeros((nlp.n_primals(), nlp.n_primals()))
                r_index = positions["fs.unit.r"]
                q_index = positions["fs.unit.q"]
                lambda_flow = multiplier_by_name["fs.unit.volume_balance"]
                lambda_pressure = multiplier_by_name["eos_pressure"]
                expected_matrix[q_index, q_index] = objective_factor
                expected_matrix[r_index, q_index] = -lambda_flow
                expected_matrix[q_index, r_index] = -lambda_flow
                expected_matrix[r_index, r_index] = (
                    lambda_pressure * DENSITY_SCALE**2 * float(native_hessian[1, 1])
                    / PRESSURE_SCALE
                )
                actual_matrix = np.asarray(nlp.evaluate_hessian_lag().toarray(), dtype=float)
                matrix_error = float(np.linalg.norm(actual_matrix - expected_matrix, ord=np.inf))
                matrix_target = tolerance(expected_matrix, actual_matrix, 1e-7, 1e-4)
                if not np.all(np.isfinite(actual_matrix)) or matrix_error > matrix_target:
                    raise RuntimeError(
                        f"whole-NLP Hessian entries mismatch at start={start}, "
                        f"objective_factor={objective_factor}, multipliers={multiplier}: "
                        f"error={matrix_error}, target={matrix_target}"
                    )
                rows.append({
                    "kind": "matrix", "start_rho": start[0], "start_Q": start[1],
                    "objective_factor": objective_factor,
                    "multipliers": ",".join(str(item) for item in multiplier),
                    "matrix_error": matrix_error, "matrix_target": matrix_target,
                    "expected_H_rr": expected_matrix[r_index, r_index],
                    "expected_H_rq": expected_matrix[r_index, q_index],
                    "expected_H_qq": expected_matrix[q_index, q_index], "passed": True,
                })
                if output_path is not None:
                    write_rows(output_path, rows)
                    print(json.dumps({"e": "row", "p": "whole_nlp", "count": len(rows)}), flush=True)
                matrix_cases += 1
                for direction_index, direction in enumerate(directions):
                    nlp.set_primals(x)
                    nlp.set_obj_factor(objective_factor)
                    nlp.set_duals(duals)
                    expected = expected_matrix @ direction
                    estimates = []
                    for step in STEPS:
                        centered = (
                            gradient(x + step * direction, objective_factor, duals)
                            - gradient(x - step * direction, objective_factor, duals)
                        ) / (2.0 * step)
                        estimates.append(np.asarray(centered, dtype=float))
                    ladder = ladder_summary(estimates, expected, STEPS, 1e-7, 1e-4)
                    maximum = max(maximum, max(ladder["errors"]))
                    maximum_floor = max(maximum_floor, ladder["floor"])
                    cases += 1
                    for step, estimate, error, target in zip(
                        STEPS, estimates, ladder["errors"], ladder["targets"]
                    ):
                        rows.append({
                            "kind": "directional", "start_rho": start[0], "start_Q": start[1],
                            "objective_factor": objective_factor,
                            "multipliers": ",".join(str(item) for item in multiplier),
                            "direction": direction_index,
                            "step": step,
                            "error": error,
                            "target": target,
                            "passed": error <= target,
                            "floor": ladder["floor"],
                            "floor_target": ladder["floor_target"],
                            "floor_support": ladder["floor_support"],
                            "mode": ladder["mode"],
                        })
                    if output_path is not None:
                        write_rows(output_path, rows)
                        print(json.dumps({"e": "row", "p": "whole_nlp", "count": len(rows)}), flush=True)
    return {
        "maximum_scaled_lagrangian_hessian_error": maximum,
        "maximum_estimated_floor": maximum_floor,
        "cases": cases,
        "matrix_cases": matrix_cases,
        "directions_per_case": len(directions),
        "expected_step_rows": cases * len(STEPS),
        "invalid_domain_rejected": invalid_rejected,
        "changed_input_rejected": changed_input_rejected,
        "valid_after_failure": valid_after_failure,
    }


def pressure_cap_state(model, cap: float) -> tuple[float, float]:
    state = model.state_at_pressure(400.0, cap, [1.0], epcsaft._core.DensityBranch.Vapor)
    if not math.isfinite(state.pressure) or abs(state.pressure - cap) > 2e-4:
        raise RuntimeError("native pressure-cap density missed pressure target")
    return state.molar_density, state.pressure


def dispersion_reference(density: float) -> float:
    with localcontext() as context:
        context.prec = 50
        pi = Decimal("3.1415926535897932384626433832795028841971693993751")
        number = Decimal("6.02214076e23")
        c = number * Decimal("1e-30")
        sigma = Decimal("3.6184")
        epsilon = Decimal("208.11")
        temperature = Decimal(400)
        exponent = epsilon / temperature
        diameter = sigma * (Decimal(1) - Decimal("0.12") * (-Decimal(3) * exponent).exp())
        eta = pi * c * Decimal(str(density)) * diameter**3 / Decimal(6)
        a_coefficients = tuple(Decimal(str(item)) for item in (
            0.9105631445, 0.6361281449, 2.6861347891, -26.547362491,
            97.759208784, -159.59154087, 91.297774084,
        ))
        b_coefficients = tuple(Decimal(str(item)) for item in (
            0.7240946941, 2.2382791861, -4.0025849485, -21.003576815,
            26.855641363, 206.55133841, -355.60235612,
        ))
        integral_one = sum(coefficient * eta**index
                           for index, coefficient in enumerate(a_coefficients))
        integral_two = sum(coefficient * eta**index
                           for index, coefficient in enumerate(b_coefficients))
        c1 = Decimal(1) / (Decimal(1) + (Decimal(8) * eta - Decimal(2) * eta**2)
                           / (Decimal(1) - eta)**4)
        result = (-Decimal(2) * pi * c * Decimal(str(density)) * integral_one
                  * sigma**3 * exponent
                  - pi * c * Decimal(str(density)) * c1 * integral_two
                  * sigma**3 * exponent**2)
        return float(result)


def independent_reference_checks(
    model, output_path: Path | None = None,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    reference_document = json.loads(PARAMETERS.read_text(encoding="utf-8"))
    for coefficient in reference_document["components"][0]["coefficients"]:
        if coefficient["family"] == "segment_count":
            coefficient["value"]["magnitude"] = 1.0
    m_one = epcsaft.Parameters.from_mapping(reference_document)
    reference_model = epcsaft._core.HelmholtzModel(m_one._native)
    dispersion_rows = []
    for density in DENSITIES:
        actual = reference_model.state(400.0, density, [1.0]).dispersion
        expected = dispersion_reference(density)
        error = abs(actual - expected)
        target = 1e-12 + 1e-9 * max(abs(actual), abs(expected))
        passed = all(map(math.isfinite, (actual, expected, error, target))) and error <= target
        dispersion_rows.append({
            "check": "gross2001_m1_dispersion",
            "density_mol_per_m3": density,
            "actual": actual,
            "reference": expected,
            "error": error,
            "target": target,
            "passed": passed,
        })
        if output_path is not None:
            write_rows(output_path, dispersion_rows)
        if not passed:
            raise RuntimeError(f"independent dispersion reference mismatch: {dispersion_rows[-1]}")

    m = 2.002
    sigma = 3.6184
    epsilon_over_k = 208.11
    temperature = 400.0
    c = 6.02214076e23 * 1e-30
    exponent = epsilon_over_k / temperature
    diameter = sigma * (1.0 - 0.12 * math.exp(-3.0 * exponent))
    chain_factor_1 = (m - 1.0) / m
    chain_factor_2 = chain_factor_1 * (m - 2.0) / m
    a0 = 0.9105631445 + chain_factor_1 * -0.3084016918 + chain_factor_2 * -0.0906148351
    b0 = 0.7240946941 + chain_factor_1 * -0.5755498075 + chain_factor_2 * 0.0976883116
    k_eta = math.pi * c * m * diameter**3 / 6.0
    virial = ((1.5 * m + 2.5) * k_eta
              - 2.0 * math.pi * c * a0 * m * m * sigma**3 * exponent
              - math.pi * c * m * b0 * m * m * sigma**3 * exponent**2)
    if abs(virial) <= 1e-8:
        raise RuntimeError("second-virial reference is not a resolved signal")
    dilute_rows = []
    for density in (1.0, 0.1, 0.01, 0.001):
        pressure = model.neutral_property(
            temperature, density, epcsaft._core.NeutralProperty.Pressure,
            epcsaft._core.DerivativeOrder.Gradient,
        )
        q1 = pressure.value / (density * GAS_CONSTANT * temperature) - 1.0
        q2 = pressure.gradient[1] / (GAS_CONSTANT * temperature) - 1.0
        q1_reference = virial * density
        q2_reference = 2.0 * virial * density
        q1_error = abs(q1 - q1_reference)
        q2_error = abs(q2 - q2_reference)
        q1_target = 0.01 * abs(q1_reference) + 1e-12
        q2_target = 0.01 * abs(q2_reference) + 1e-12
        passed = q1_error <= q1_target and q2_error <= q2_target
        dilute_rows.append({
            "check": "propane_second_virial",
            "density_mol_per_m3": density,
            "q1": q1,
            "q1_reference": q1_reference,
            "q1_error": q1_error,
            "q1_target": q1_target,
            "q2": q2,
            "q2_reference": q2_reference,
            "q2_error": q2_error,
            "q2_target": q2_target,
            "passed": passed,
        })
        if output_path is not None:
            write_rows(output_path, dispersion_rows + dilute_rows)
        if not passed:
            raise RuntimeError(f"second-virial reference mismatch: {dilute_rows[-1]}")
    if max(abs(dilute_rows[-1][key]) for key in ("q1", "q2")) > 1e-5:
        raise RuntimeError("native dilute ideal limit did not resolve")
    return dispersion_rows, dilute_rows


def reference_summary(rows: list[dict[str, object]]) -> dict[str, object]:
    groups = {
        "dispersion": [row for row in rows if row["check"] == "gross2001_m1_dispersion"],
        "second_virial": [row for row in rows if row["check"] == "propane_second_virial"],
    }
    return {
        "row_count": len(rows),
        "all_passed": all(bool(row["passed"]) for row in rows),
        "groups": {
            name: {
                "row_count": len(group),
                "all_passed": all(bool(row["passed"]) for row in group),
                "maximum_abs_errors": {
                    field: max(
                        abs(float(row[field])) for row in group if row.get(field) not in (None, "")
                    )
                    if any(row.get(field) not in (None, "") for row in group)
                    else 0.0
                    for field in ("error", "q1_error", "q2_error")
                },
            }
            for name, group in groups.items()
        },
    }


class CampaignAttempt:
    def __init__(self, output_dir: Path) -> None:
        self.output_dir = output_dir
        self.attempt_id = f"{time.strftime('%Y%m%dT%H%M%S', time.gmtime())}-{os.getpid()}"

    def start(self, identity: dict[str, object]) -> None:
        if self.output_dir.exists() and any(self.output_dir.iterdir()):
            raise RuntimeError(f"campaign output directory is not empty: {self.output_dir}")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._write_summary({**identity, "status": "running", "stage": "started"})

    def checkpoint(self, stage: str, rows: dict[str, list[dict[str, object]]] | None = None,
                   **metadata: object) -> None:
        if rows:
            for name, values in rows.items():
                if values:
                    write_rows(self.output_dir / f"{name}.csv", values)
        self._write_summary({"status": "running", "stage": stage, **metadata})

    def fail(self, error: BaseException, stage: str = "failed") -> None:
        partial_rows = {}
        for path in self.output_dir.glob("*.csv"):
            with path.open(newline="", encoding="utf-8") as stream:
                partial_rows[path.stem] = max(sum(1 for _ in stream) - 1, 0)
        self._write_summary({
            "status": "failed attempt", "stage": stage,
            "error": repr(error), "error_type": type(error).__name__,
            "partial_rows": partial_rows,
        })

    def finish(self, summary: dict[str, object]) -> None:
        summary["attempt_id"] = self.attempt_id
        summary["lifecycle"] = "accepted"
        self._write_summary(summary)

    def _write_summary(self, payload: dict[str, object]) -> None:
        path = self.output_dir / "summary.json"
        summary = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
        if summary.get("attempt_id", self.attempt_id) != self.attempt_id:
            raise RuntimeError("campaign summary belongs to another attempt")
        summary.update({"attempt_id": self.attempt_id, **payload})
        temporary = self.output_dir / "summary.json.tmp"
        temporary.write_text(
            json.dumps(summary, indent=2, sort_keys=True, default=str) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, path)


def process_group_alive(pid: int) -> bool:
    try:
        os.killpg(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def terminate_process_group(process: subprocess.Popen[str]) -> dict[str, object]:
    signals = []
    try:
        os.killpg(process.pid, signal.SIGTERM)
        signals.append("SIGTERM")
    except ProcessLookupError:
        process.wait()
        return {"signals": signals, "group_terminated": not process_group_alive(process.pid)}
    term_deadline = time.monotonic() + TERM_GRACE_S
    try:
        process.wait(timeout=TERM_GRACE_S)
    except subprocess.TimeoutExpired:
        pass
    remaining = term_deadline - time.monotonic()
    if remaining > 0.0 and process_group_alive(process.pid):
        time.sleep(remaining)
    if process_group_alive(process.pid):
        try:
            os.killpg(process.pid, signal.SIGKILL)
            signals.append("SIGKILL")
        except ProcessLookupError:
            pass
        if process.poll() is None:
            process.wait()
        time.sleep(KILL_GRACE_S)
    return {"signals": signals, "group_terminated": not process_group_alive(process.pid)}


@contextmanager
def managed_child(command):
    with subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        env=os.environ.copy(), start_new_session=True,
    ) as process:
        try:
            yield process
        finally:
            if process.poll() is None or process_group_alive(process.pid):
                termination = terminate_process_group(process)
                if not termination["group_terminated"]:
                    raise RuntimeError(f"child process group {process.pid} survived cleanup")


def interrupt_campaign(signum: int, _frame) -> None:
    raise KeyboardInterrupt(f"campaign interrupted by {signal.Signals(signum).name}")


def collect_child_output(
    process: subprocess.Popen[str], timeout: float,
    watch=None, phase_deadline=None,
) -> dict[str, object]:
    selector = selectors.DefaultSelector()
    buffers = {"stdout": bytearray(), "stderr": bytearray()}
    for name, stream in (("stdout", process.stdout), ("stderr", process.stderr)):
        if stream is not None:
            os.set_blocking(stream.fileno(), False)
            selector.register(stream, selectors.EVENT_READ, name)
    deadline = time.monotonic() + timeout
    drain_deadline = None
    entry_observed = False
    timed_out = False
    timeout_scope = "total"
    while selector.get_map():
        phase_limit = phase_deadline() if phase_deadline is not None else None
        if phase_limit is not None and time.monotonic() >= phase_limit:
            timed_out = True
            timeout_scope = "phase"
            break
        if process.poll() is None:
            limit = deadline
            if phase_limit is not None and phase_limit < limit:
                limit = phase_limit
                timeout_scope = "phase"
            else:
                timeout_scope = "total"
            remaining = limit - time.monotonic()
            if remaining <= 0.0:
                timed_out = True
                break
        else:
            if drain_deadline is None:
                drain_deadline = time.monotonic() + 2.0
            remaining = drain_deadline - time.monotonic()
            if remaining <= 0.0:
                break
        events = selector.select(min(remaining, 0.25))
        for key, _ in events:
            active_phase_limit = phase_deadline() if phase_deadline is not None else None
            try:
                chunk = os.read(key.fileobj.fileno(), 64)
            except BlockingIOError:
                continue
            if not chunk:
                selector.unregister(key.fileobj)
                continue
            name = key.data
            buffers[name].extend(chunk)
            if name == "stdout" and watch is not None:
                entry_observed = bool(watch(bytes(buffers[name]).decode(errors="replace")))
                if entry_observed:
                    selector.close()
                    return {
                        "stdout": bytes(buffers["stdout"]).decode(errors="replace"),
                        "stderr": bytes(buffers["stderr"]).decode(errors="replace"),
                        "timed_out": False, "entry_observed": True,
                        "timeout_scope": timeout_scope,
                    }
            if active_phase_limit is not None and time.monotonic() >= active_phase_limit:
                timed_out = True
                timeout_scope = "phase"
                break
        if timed_out:
            break
    selector.close()
    return {
        "stdout": bytes(buffers["stdout"]).decode(errors="replace"),
        "stderr": bytes(buffers["stderr"]).decode(errors="replace"),
        "timed_out": timed_out, "entry_observed": entry_observed,
        "timeout_scope": timeout_scope,
    }


def solver_nl_from_output(output: str) -> Path | None:
    for line in output.splitlines():
        match = re.search(r"([^\s'\"]+\.nl)", line)
        if match is None:
            match = re.search(r"([^\s'\"]+)\s+nl(?:\s|$)", line)
        if match:
            path = Path(match.group(1) if match.group(1).endswith(".nl") else match.group(1) + ".nl")
            if path.is_file():
                return path.resolve()
    return None


def iteration_zero_line(output: str) -> str | None:
    for line in output.splitlines():
        if re.match(
            r"^\s*0\s+[-+0-9.]+(?:[eE][-+0-9]+)?"
            r"(?:\s+[-+0-9.]+(?:[eE][-+0-9]+)?){2,}", line
        ):
            return line.strip()
    return None


def live_ipopt_process(
    process_group: int, executable: Path, nl_path: Path | None = None
) -> dict[str, object] | None:
    expected = executable.resolve()
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        try:
            if os.getpgid(pid) != process_group:
                continue
            actual_executable = Path(os.readlink(entry / "exe")).resolve()
            command_line = (entry / "cmdline").read_bytes().replace(b"\0", b" ").decode(errors="replace").strip()
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        command_matches_nl = nl_path is None or str(nl_path) in command_line or (
            str(nl_path.with_suffix("")) + " nl" in command_line
        )
        if actual_executable == expected and command_matches_nl:
            return {
                "pid": pid, "process_group": process_group,
                "executable": str(actual_executable), "cmdline": command_line,
            }
    return None


def retain_solver_maps(
    output_dir: Path, case_name: str, nl_path: Path,
) -> tuple[Path, Path]:
    paths = [nl_path.with_suffix(suffix) for suffix in (".row", ".col")]
    if case_name == "process-0":
        return paths[0], paths[1]
    canonical = [
        output_dir / "solver" / "process-0" / f"process-0{suffix}"
        for suffix in (".row", ".col")
    ]
    for source, retained in zip(paths, canonical, strict=True):
        if not source.is_file() and retained.is_file():
            continue
        if sha256(source) != sha256(retained):
            raise RuntimeError(f"symbolic solver {source.suffix} changed between cases")
        source.unlink()
    return canonical[0], canonical[1]


def retained_iteration_line(case_dir: Path, fallback: str | None) -> str | None:
    path = case_dir / "ipopt.log"
    if path.is_file():
        return iteration_zero_line(path.read_text(errors="replace")) or fallback
    return fallback


def timed_process_solve(
    cap: float, start: tuple[float, float], derivative_test: bool, ipopt: Path,
    output_dir: Path, nl_path: Path, case_name: str, timeout: float = 30.0,
    forced_timeout: bool = False, deadline_scope: str = "solve",
) -> dict[str, object]:
    output_dir = output_dir.resolve()
    command = (
        sys.executable, "-u", str(Path(__file__).resolve()), "--solve-case",
        str(cap), str(start[0]), str(start[1]), "--derivative-test", str(derivative_test),
        "--ipopt", str(ipopt), "--output-dir", str(output_dir), "--case-name", case_name,
        "--nl-path", str(nl_path.resolve()),
    )
    if forced_timeout:
        command += ("--entry-probe",)
    started = time.perf_counter()
    (output_dir / "solver" / case_name).mkdir(parents=True, exist_ok=True)
    with managed_child(command) as process:
        case_dir = output_dir / "solver" / case_name

        def probe_entry(output: str) -> bool:
            return iteration_zero_line(output) is not None

        watch = probe_entry
        collected = collect_child_output(process, timeout, watch if forced_timeout else None)
        stdout = str(collected["stdout"])
        stderr = str(collected["stderr"])
        group = process.pid
        if not forced_timeout and not collected["timed_out"]:
            process.wait()
        if collected["entry_observed"] and forced_timeout:
            live = live_ipopt_process(group, ipopt)
            nl_path = solver_nl_from_output(live["cmdline"]) if live is not None else None
            iteration = iteration_zero_line(stdout)
            live = live_ipopt_process(group, ipopt, nl_path) if live is not None and nl_path is not None else None
            if live is not None and nl_path is not None:
                cancellation_deadline = time.monotonic()
                termination = terminate_process_group(process)
                drained = collect_child_output(process, 2.0)
                stdout += str(drained["stdout"])
                stderr += str(drained["stderr"])
                iteration = retained_iteration_line(case_dir, iteration)
                record = {
                    "termination": "forced-post-entry-timeout", "status": "failed",
                    "timed_out": True, "entry_observed": True,
                    "numeric_iteration_zero_line": iteration,
                    "solver_nl": str(nl_path),
                    "solver_row": str(nl_path.with_suffix(".row")),
                    "solver_col": str(nl_path.with_suffix(".col")),
                    "solver_nl_sha256": sha256(nl_path), "ipopt_process": live,
                    "cancellation_deadline_monotonic": cancellation_deadline,
                    "signals": termination["signals"],
                    "process_group_terminated": termination["group_terminated"],
                }
            else:
                termination = terminate_process_group(process)
                drained = collect_child_output(process, 2.0)
                stdout += str(drained["stdout"])
                stderr += str(drained["stderr"])
                record = {
                    "termination": "inconclusive-completed-before-cancellation", "status": "failed",
                    "timed_out": False, "entry_observed": True,
                    "numeric_iteration_zero_line": iteration,
                    "solver_nl": str(nl_path) if nl_path is not None else None,
                    "signals": termination["signals"],
                    "process_group_terminated": termination["group_terminated"],
                }
        elif forced_timeout:
            termination = terminate_process_group(process)
            drained = collect_child_output(process, 2.0)
            stdout += str(drained["stdout"])
            stderr += str(drained["stderr"])
            record = {
                "termination": "startup-timeout" if collected["timed_out"]
                else "inconclusive-completed-before-cancellation", "status": "failed",
                "timed_out": bool(collected["timed_out"]),
                "signals": termination["signals"],
                "process_group_terminated": termination["group_terminated"],
            }
        elif collected["timed_out"]:
            termination = terminate_process_group(process)
            drained = collect_child_output(process, 2.0)
            stdout += str(drained["stdout"])
            stderr += str(drained["stderr"])
            record = {
                "termination": "timeout", "status": "failed", "timed_out": True,
                "timeout_scope": deadline_scope,
                "signals": termination["signals"],
                "process_group_terminated": termination["group_terminated"],
            }
        else:
            output = stdout.strip().splitlines()
            if process.returncode != 0:
                record = {
                    "termination": "child-error", "status": "failed",
                    "error": f"solve child exited with {process.returncode}",
                }
            elif not output:
                record = {
                    "termination": "child-parse-error", "status": "failed",
                    "error": "solve child emitted no JSON record",
                }
            else:
                try:
                    record = json.loads(output[-1])
                except json.JSONDecodeError as error:
                    record = {
                        "termination": "child-parse-error", "status": "failed",
                        "error": repr(error),
                    }
            record["timed_out"] = False
            record["process_group_terminated"] = not process_group_alive(group)
        case_dir.mkdir(parents=True, exist_ok=True)
        stdout_path = case_dir / "parent.stdout.log"
        stderr_path = case_dir / "parent.stderr.log"
        stdout_path.write_text(stdout, encoding="utf-8")
        if stderr:
            stderr_path.write_text(stderr, encoding="utf-8")
        else:
            stderr_path.unlink(missing_ok=True)
        record["cap"] = cap
        record["start_rho"] = start[0]
        record["start_Q"] = start[1]
        record["outer_elapsed_s"] = time.perf_counter() - started
        record["returncode"] = process.returncode
        record["process_group"] = group
        record["stdout_tail"] = stdout[-2000:]
        record["stderr_tail"] = stderr[-2000:]
        record["parent_stdout_log"] = str(stdout_path)
        record["parent_stderr_log"] = str(stderr_path) if stderr else None
        if record.get("solver_nl"):
            row_path, col_path = retain_solver_maps(output_dir, case_name, nl_path)
            record["solver_row"] = str(row_path)
            record["solver_col"] = str(col_path)
        record.setdefault(
            "timeout_scope",
            "forced-post-entry" if record.get("termination") == "forced-post-entry-timeout"
            else deadline_scope if record.get("timed_out") else "none",
        )
        return record


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    fields = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def append_failure_row(path: Path, failure: dict[str, object]) -> None:
    rows: list[dict[str, object]] = []
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as stream:
            rows.extend(csv.DictReader(stream))
    rows.append(failure)
    write_rows(path, rows)


def recording_count(diagnostic) -> int:
    count = float(diagnostic.evaluate([]))
    if not math.isfinite(count) or count < 0.0 or count != round(count):
        raise RuntimeError(f"invalid ASL recording counter: {count}")
    return int(count)


def cost_child(mode: str, library: Path) -> dict[str, object]:
    states = tuple((temperature, density) for temperature in TEMPERATURES for density in DENSITIES)
    orders = (
        (0, epcsaft._core.DerivativeOrder.Value),
        (1, epcsaft._core.DerivativeOrder.Gradient),
        (2, epcsaft._core.DerivativeOrder.Hessian),
    )
    rows: list[dict[str, object]] = []
    started = time.perf_counter()

    def phase_event(state: str, name: str, limit: float | None = None) -> None:
        event = {"e": state, "p": name}
        if limit is not None:
            event["l"] = limit
            event["t"] = time.monotonic()
        print(json.dumps(event, separators=(",", ":")), flush=True)

    phase_event("start", "cold", 60.0)
    if mode == "native":
        evaluator = native_model()
        evaluator.reset_counters()
        def counter():
            return int(evaluator.counters().recordings)
    elif mode == "asl":
        evaluator = asl_functions(library)[1]
        diagnostic = ExternalFunction(
            library=str(library), function="epcsaft_neutral_recordings",
            units=pyunits.dimensionless, arg_units=[],
        )
        def counter():
            return recording_count(diagnostic)
    else:
        raise RuntimeError(f"unknown cost mode: {mode}")
    cold_before = counter()
    cold_call_started = time.perf_counter()
    if mode == "native":
        native_fgh(evaluator, epcsaft._core.NeutralProperty.Pressure,
                   states[0][0], states[0][1], epcsaft._core.DerivativeOrder.Value)
    else:
        external_fgh(evaluator, states[0][0], states[0][1], 0)
    cold_call_elapsed = time.perf_counter() - cold_call_started
    cold_after = counter()
    cold_elapsed = time.perf_counter() - started
    rows.append({
        "path": mode, "kind": "cold_value", "actual_derivative_order": 0,
        "batch": "cold", "call_count": 1, "elapsed_s": cold_call_elapsed,
        "recording_before": cold_before, "recording_after": cold_after,
        "recording_delta": cold_after - cold_before,
    })
    print(json.dumps({"e": "row", "r": rows[-1]}, separators=(",", ":")), flush=True)
    if cold_elapsed > 60.0 or cold_after - cold_before != 1:
        raise RuntimeError(f"{mode} cold setup/tape measurement failed: {rows[-1]}")
    phase_event("done", "cold")
    first_calls = {0: cold_call_elapsed}
    for order_number, order in orders[1:]:
        phase_name = "first_gradient" if order_number == 1 else "first_hessian"
        phase_event("start", phase_name, 10.0)
        before = counter()
        call_started = time.perf_counter()
        if mode == "native":
            native_fgh(evaluator, epcsaft._core.NeutralProperty.Pressure,
                       states[0][0], states[0][1], order)
        else:
            external_fgh(evaluator, states[0][0], states[0][1], order_number)
        elapsed = time.perf_counter() - call_started
        after = counter()
        first_calls[order_number] = elapsed
        rows.append({
            "path": mode, "kind": "first_call", "actual_derivative_order": order_number,
            "batch": "first", "call_count": 1, "elapsed_s": elapsed,
            "recording_before": before, "recording_after": after,
            "recording_delta": after - before,
        })
        print(json.dumps({"e": "row", "r": rows[-1]}, separators=(",", ":")), flush=True)
        if after - before != 0:
            raise RuntimeError(f"{mode} derivative order {order_number} re-recorded its prepared tape")
        phase_event("done", phase_name)
    warm: dict[str, object] = {}
    for order_number, order in orders:
        batches = []
        for batch in range(1, 4):
            phase_event("start", f"order_{order_number}_batch_{batch}", 10.0)
            before = counter()
            batch_started = time.perf_counter()
            for index in range(100):
                temperature, density = states[index % len(states)]
                if mode == "native":
                    native_fgh(evaluator, epcsaft._core.NeutralProperty.Pressure,
                               temperature, density, order)
                else:
                    external_fgh(evaluator, temperature, density, order_number)
            elapsed = time.perf_counter() - batch_started
            after = counter()
            row = {
                "path": mode, "kind": "warm_batch",
                "actual_derivative_order": order_number, "batch": batch,
                "call_count": 100, "elapsed_s": elapsed,
                "recording_before": before, "recording_after": after,
                "recording_delta": after - before,
            }
            rows.append(row)
            print(json.dumps({"e": "row", "r": row}, separators=(",", ":")), flush=True)
            if elapsed > 10.0:
                raise RuntimeError(
                    f"{mode} derivative order {order_number} batch {batch} exceeded 10 seconds"
                )
            if after - before != 0:
                raise RuntimeError(
                    f"{mode} derivative order {order_number} batch {batch} re-recorded its prepared tape"
                )
            batches.append(elapsed)
            phase_event("done", f"order_{order_number}_batch_{batch}")
        warm[str(order_number)] = {
            "median_s": statistics.median(batches),
            "range_s": [min(batches), max(batches)],
            "call_count": 100, "batch_count": len(batches),
        }
    return {
        "mode": mode, "fresh_process": True, "rows": rows,
        "summary": {
            "cold_elapsed_s": cold_elapsed, "first_call_s": first_calls,
            "warm": warm, "state_count": len(states), "states": states,
            "cold_recording_count": cold_after,
            "warm_recording_deltas": [
                row["recording_delta"] for row in rows
                if row["kind"] in ("first_call", "warm_batch")
            ],
            "recording_reuse_passed": all(
                row["recording_delta"] == 0
                for row in rows if row["kind"] in ("first_call", "warm_batch")
            ),
        },
    }


def cost_experiment(
    output_dir: Path, ipopt: Path, total_deadline: float,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    rows: list[dict[str, object]] = []
    summaries: dict[str, object] = {}
    for mode in ("native", "asl"):
        remaining = total_deadline - time.monotonic() - TERMINATION_BUDGET_S
        if remaining <= 0.0:
            failure = {
                "path": mode, "kind": "failure", "status": "failed",
                "phase": "startup", "failure": "post-build campaign deadline exceeded",
                "timeout_scope": "total", "phase_limit_s": None,
            }
            append_failure_row(output_dir / "cost.csv", failure)
            raise RuntimeError(f"cost child {mode} exceeded the post-build deadline: {failure}")
        command = (
            sys.executable, "-u", str(Path(__file__).resolve()), "--cost-mode", mode,
            "--output-dir", str(output_dir.resolve()), "--ipopt", str(ipopt),
        )
        child_started = time.perf_counter()
        phase = {"name": None, "deadline": None, "limit_s": None, "events": []}
        processed_lines = 0
        child_startup_elapsed = None
        mode_rows: list[dict[str, object]] = []

        def observe_progress(
            output: str, phase=phase, child_started=child_started,
            mode_rows=mode_rows,
        ) -> bool:
            nonlocal processed_lines, child_startup_elapsed
            for line in output.splitlines(keepends=True)[processed_lines:]:
                if not line.endswith(("\n", "\r")):
                    break
                processed_lines += 1
                try:
                    event = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if event.get("e") == "start":
                    phase["name"] = event["p"]
                    phase["limit_s"] = float(event["l"])
                    phase["deadline"] = float(event["t"]) + phase["limit_s"]
                    phase["events"].append(event)
                    if child_startup_elapsed is None:
                        child_startup_elapsed = time.perf_counter() - child_started
                elif event.get("e") == "done":
                    phase["events"].append(event)
                    phase["name"] = None
                    phase["deadline"] = None
                    phase["limit_s"] = None
                elif event.get("e") == "row":
                    row = dict(event["r"])
                    row["child_startup_elapsed_s"] = child_startup_elapsed
                    mode_rows.append(row)
                    write_rows(output_dir / "cost.csv", rows + mode_rows)
            return False

        def current_phase_deadline(phase=phase) -> float | None:
            return phase["deadline"]

        with managed_child(command) as process:
            collected = collect_child_output(
                process, remaining, watch=observe_progress,
                phase_deadline=current_phase_deadline,
            )
            child_total_elapsed = time.perf_counter() - child_started
            stdout = str(collected["stdout"])
            stderr = str(collected["stderr"])
            if collected["timed_out"]:
                termination = terminate_process_group(process)
                drained = collect_child_output(process, 2.0)
                stdout += str(drained["stdout"])
                stderr += str(drained["stderr"])
                failure = {
                    "path": mode, "kind": "failure", "status": "failed",
                    "phase": phase["name"] or "startup",
                    "failure": "cost child deadline exceeded",
                    "timeout_scope": collected["timeout_scope"],
                    "phase_limit_s": phase["limit_s"],
                    "child_startup_elapsed_s": child_startup_elapsed,
                    "child_total_elapsed_s": child_total_elapsed,
                    "signals": termination["signals"],
                    "process_group_terminated": termination["group_terminated"],
                    "returncode": process.returncode,
                    "phase_events": json.dumps(phase["events"], separators=(",", ":")),
                    "stdout_tail": stdout[-2000:], "stderr_tail": stderr[-2000:],
                }
                for row in mode_rows:
                    row["child_total_elapsed_s"] = child_total_elapsed
                write_rows(output_dir / "cost.csv", rows + mode_rows + [failure])
                raise RuntimeError(f"cost child {mode} timed out: {failure}")
            process.wait()
            output = stdout.strip().splitlines()
            if process.returncode != 0 or not output:
                failure = {
                    "path": mode, "kind": "failure", "status": "failed",
                    "phase": phase["name"] or "startup", "failure": "cost child failed",
                    "timeout_scope": "none", "phase_limit_s": phase["limit_s"],
                    "child_startup_elapsed_s": child_startup_elapsed,
                    "child_total_elapsed_s": child_total_elapsed,
                    "returncode": process.returncode,
                    "phase_events": json.dumps(phase["events"], separators=(",", ":")),
                    "stdout_tail": stdout[-2000:], "stderr_tail": stderr[-2000:],
                }
                for row in mode_rows:
                    row["child_total_elapsed_s"] = child_total_elapsed
                write_rows(output_dir / "cost.csv", rows + mode_rows + [failure])
                raise RuntimeError(f"cost child {mode} failed: {failure}")
            try:
                result = json.loads(output[-1])
            except json.JSONDecodeError as error:
                failure = {
                    "path": mode, "kind": "failure", "status": "failed",
                    "phase": phase["name"] or "result", "failure": repr(error),
                    "timeout_scope": "none", "phase_limit_s": phase["limit_s"],
                    "child_startup_elapsed_s": child_startup_elapsed,
                    "child_total_elapsed_s": child_total_elapsed,
                    "returncode": process.returncode,
                    "phase_events": json.dumps(phase["events"], separators=(",", ":")),
                    "stdout_tail": stdout[-2000:], "stderr_tail": stderr[-2000:],
                }
                for row in mode_rows:
                    row["child_total_elapsed_s"] = child_total_elapsed
                write_rows(output_dir / "cost.csv", rows + mode_rows + [failure])
                raise RuntimeError(f"cost child {mode} emitted invalid JSON: {failure}") from error
            streamed = [
                {key: value for key, value in row.items() if not key.startswith("child_")}
                for row in mode_rows
            ]
            if result["rows"] != streamed:
                failure = {
                    "path": mode, "kind": "failure", "status": "failed",
                    "phase": "result", "failure": "cost child stream/result rows disagree",
                    "timeout_scope": "none", "child_total_elapsed_s": child_total_elapsed,
                }
                write_rows(output_dir / "cost.csv", rows + mode_rows + [failure])
                raise RuntimeError(f"cost child {mode} stream/result rows disagree")
            for row in mode_rows:
                row["child_startup_elapsed_s"] = child_startup_elapsed
                row["child_total_elapsed_s"] = child_total_elapsed
            rows.extend(mode_rows)
            write_rows(output_dir / "cost.csv", rows)
            summaries[mode] = {
                "fresh_process": True, "child_mode": mode,
                "child_startup_elapsed_s": child_startup_elapsed,
                "child_total_elapsed_s": child_total_elapsed,
                "phase_count": len(phase["events"]) // 2,
                **result["summary"],
            }
    return rows, summaries


def numerical_child(
    stage: str, output_dir: Path, ipopt: Path, nl_path: Path | None,
) -> dict[str, object]:
    target = output_dir / ("lagrangian.csv" if stage == "whole-nlp" else "derivatives.csv")
    try:
        wheel, wheel_hash, library = installed_identity()
        print(json.dumps({"e": "start", "p": stage, "t": time.monotonic()}), flush=True)
        if stage == "whole-nlp":
            if nl_path is None or not nl_path.is_file():
                raise RuntimeError(f"whole-NLP solver file is unavailable: {nl_path}")
            rows: list[dict[str, object]] = []
            summary = lagrangian_hessian_checks(library, nl_path, rows, target)
            print(json.dumps({"e": "done", "p": stage}), flush=True)
            return {"status": "ok", "whole_nlp": summary}

        os.environ["AMPLFUNC"] = str(library)
        library_identity = shared_library_identity(library)
        ipopt_identity = ipopt_preflight(ipopt)
        model = native_model()
        functions = asl_functions(library)
        derivative_rows: list[dict[str, object]] = []
        try:
            derivative_summary = derivative_checks(
                model, functions, library, derivative_rows, target,
            )
        finally:
            if derivative_rows:
                write_rows(target, derivative_rows)

        domain_rows = []
        conditioning_states = (
            ((temperature, float(density), "conditioning-grid")
             for temperature in GRID_TEMPERATURES
             for density in np.linspace(5.0, 110.0, 9)),
            ((temperature, density, "derivative-state")
             for temperature in TEMPERATURES for density in DENSITIES),
        )
        for states in conditioning_states:
            for temperature, density, source in states:
                pressure = native_fgh(
                    model, epcsaft._core.NeutralProperty.Pressure, temperature, density,
                    epcsaft._core.DerivativeOrder.Gradient,
                )
                ratio = pressure[1][1] / (GAS_CONSTANT * temperature)
                passed = math.isfinite(pressure[0]) and math.isfinite(ratio) and ratio > 0.8
                domain_rows.append({
                    "temperature_k": temperature,
                    "density_mol_per_m3": density,
                    "pressure_pa": pressure[0],
                    "pressure_slope_pa_m3_per_mol": pressure[1][1],
                    "pressure_slope_over_rt": ratio,
                    "source": source,
                    "passed": passed,
                })
                write_rows(output_dir / "domain.csv", domain_rows)
                if not passed:
                    raise RuntimeError("adapter domain conditioning gate failed")
        pressure_at_low = model.neutral_property(
            400.0, 10.0, epcsaft._core.NeutralProperty.Pressure,
            epcsaft._core.DerivativeOrder.Value,
        ).value
        pressure_at_high = model.neutral_property(
            400.0, 100.0, epcsaft._core.NeutralProperty.Pressure,
            epcsaft._core.DerivativeOrder.Value,
        ).value
        if not pressure_at_low < PRESSURE_SCALE < pressure_at_high:
            raise RuntimeError("pressure bracket gate failed")
        pressure_span = pressure_at_high - pressure_at_low
        if pressure_span <= 0.0:
            raise RuntimeError("pressure nonconstancy gate failed")

        changed_value = functions[1].evaluate(
            [400.0, 60.0, "2.003", *PARAMETER_STRINGS[1:]],
        )
        immutable_input_rejected = math.isnan(changed_value)
        if not immutable_input_rejected:
            raise RuntimeError("changed immutable parameter tuple was accepted")
        invalid_domain_rejected = math.isnan(
            functions[1].evaluate([400.0, -1.0, *PARAMETER_STRINGS]),
        )
        if not invalid_domain_rejected:
            raise RuntimeError("invalid density was accepted")
        dispersion_rows, dilute_rows = independent_reference_checks(
            model, output_dir / "references.csv",
        )
        reference_rows = dispersion_rows + dilute_rows
        references = reference_summary(reference_rows)
        if not references["all_passed"]:
            raise RuntimeError("independent reference acceptance gate failed")
        rho_reference, pressure_reference = pressure_cap_state(model, PRESSURE_SCALE)
        changed_rho, changed_pressure = pressure_cap_state(model, 150000.0)
        print(json.dumps({"e": "done", "p": stage}), flush=True)
        return {
            "status": "ok",
            "wheel": str(wheel) if wheel else None, "wheel_sha256": wheel_hash,
            "prepared_environment": loaded_identity(),
            "library": str(library), "library_identity": library_identity,
            "library_sha256": sha256(library), "ipopt_identity": ipopt_identity,
            "derivatives": derivative_summary,
            "conditioning": {
                "state_count": len(domain_rows),
                "grid_state_count": 45,
                "derivative_state_count": 9,
                "pressure_min_pa": min(row["pressure_pa"] for row in domain_rows),
                "pressure_max_pa": max(row["pressure_pa"] for row in domain_rows),
                "pressure_slope_over_rt_min": min(
                    row["pressure_slope_over_rt"] for row in domain_rows
                ),
                "minimum_pressure_slope": min(
                    row["pressure_slope_pa_m3_per_mol"] for row in domain_rows
                ),
            },
            "pressure_bracket": {
                "low_pa": pressure_at_low, "high_pa": pressure_at_high,
                "span_pa": pressure_span,
            },
            "immutable_input_rejected": immutable_input_rejected,
            "invalid_domain_rejected": invalid_domain_rejected,
            "references": references,
            "reduced_reference": {"rho": rho_reference, "P": pressure_reference},
            "changed_cap_reference": {"rho": changed_rho, "P": changed_pressure},
        }
    except Exception as error:
        append_failure_row(target, {
            "kind": "failure", "status": "failed", "stage": stage,
            "failure": repr(error),
        })
        return {"status": "failed", "stage": stage, "error": repr(error)}


def run_numerical_stage(
    stage: str, output_dir: Path, ipopt: Path, nl_path: Path | None,
    numerical_limit: float, total_deadline: float,
) -> tuple[dict[str, object], float]:
    target = output_dir / ("lagrangian.csv" if stage == "whole-nlp" else "derivatives.csv")
    remaining = total_deadline - time.monotonic() - TERMINATION_BUDGET_S
    if remaining <= 0.0:
        failure = {
            "kind": "failure", "status": "failed", "stage": stage,
            "failure": "post-build campaign deadline exceeded before stage start",
            "timeout_scope": "total", "process_group_terminated": True,
        }
        append_failure_row(target, failure)
        raise RuntimeError(f"numerical stage {stage} exceeded the post-build deadline")
    command = [
        sys.executable, "-u", str(Path(__file__).resolve()),
        "--numerical-stage", stage, "--output-dir", str(output_dir.resolve()),
        "--ipopt", str(ipopt),
    ]
    if nl_path is not None:
        command.extend(("--nl-path", str(nl_path.resolve())))
    started = time.perf_counter()
    phase_started = None
    phase_deadline = None
    processed_lines = 0

    def observe_progress(output: str) -> bool:
        nonlocal processed_lines, phase_started, phase_deadline
        for line in output.splitlines(keepends=True)[processed_lines:]:
            if not line.endswith(("\n", "\r")):
                break
            processed_lines += 1
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                continue
            if event.get("e") == "start" and event.get("p") == stage:
                phase_started = float(event["t"])
                phase_deadline = phase_started + numerical_limit
            elif event.get("e") == "done" and event.get("p") == stage:
                phase_deadline = None
        return False

    with managed_child(command) as process:
        collected = collect_child_output(
            process, remaining, watch=observe_progress,
            phase_deadline=lambda: phase_deadline,
        )
        stdout = str(collected["stdout"])
        stderr = str(collected["stderr"])
        if collected["timed_out"]:
            termination = terminate_process_group(process)
            drained = collect_child_output(process, TERMINATION_BUDGET_S)
            stdout += str(drained["stdout"])
            stderr += str(drained["stderr"])
            failure = {
                "kind": "failure", "status": "failed", "stage": stage,
                "failure": "numerical stage deadline exceeded",
                "timeout_scope": collected["timeout_scope"],
                "numerical_limit_s": numerical_limit,
                "signals": termination["signals"],
                "process_group_terminated": termination["group_terminated"],
                "returncode": process.returncode,
                "stdout_tail": stdout[-2000:], "stderr_tail": stderr[-2000:],
            }
            append_failure_row(target, failure)
            raise RuntimeError(f"numerical stage {stage} timed out: {failure}")
        process.wait()
        output = stdout.strip().splitlines()
        if process.returncode != 0 or not output:
            failure = {
                "kind": "failure", "status": "failed", "stage": stage,
                "failure": "numerical child failed", "timeout_scope": "none",
                "returncode": process.returncode,
                "stdout_tail": stdout[-2000:], "stderr_tail": stderr[-2000:],
            }
            append_failure_row(target, failure)
            raise RuntimeError(f"numerical stage {stage} failed: {failure}")
        try:
            result = json.loads(output[-1])
        except json.JSONDecodeError as error:
            failure = {
                "kind": "failure", "status": "failed", "stage": stage,
                "failure": repr(error), "timeout_scope": "none",
                "returncode": process.returncode,
                "stdout_tail": stdout[-2000:], "stderr_tail": stderr[-2000:],
            }
            append_failure_row(target, failure)
            raise RuntimeError(f"numerical stage {stage} emitted invalid JSON") from error
        if result.get("status") != "ok":
            raise RuntimeError(f"numerical stage {stage} failed: {result}")
        elapsed = (
            time.monotonic() - phase_started
            if phase_started is not None
            else time.perf_counter() - started
        )
        return result, elapsed


def enforce_elapsed(started: float, limit: float, stage: str) -> None:
    elapsed = time.perf_counter() - started
    if elapsed > limit:
        raise RuntimeError(f"{stage} exceeded {limit:g} seconds: {elapsed:g}")


def run_campaign(args, attempt: CampaignAttempt) -> int:
    started = time.perf_counter()
    total_deadline = time.monotonic() + 600.0
    ipopt = args.ipopt.resolve()
    if not ipopt.is_file() or not os.access(ipopt, os.X_OK):
        raise RuntimeError(f"selected Ipopt executable is unavailable: {ipopt}")
    initial, numerical_elapsed = run_numerical_stage(
        "initial", args.output_dir, ipopt, None, 120.0, total_deadline,
    )
    wheel = Path(str(initial["wheel"])) if initial["wheel"] else None
    wheel_hash = initial["wheel_sha256"]
    library = Path(str(initial["library"]))
    library_identity = dict(initial["library_identity"])
    ipopt_identity = dict(initial["ipopt_identity"])
    derivative_summary = dict(initial["derivatives"])
    conditioning = dict(initial["conditioning"])
    pressure_bracket = dict(initial["pressure_bracket"])
    immutable_input_rejected = bool(initial["immutable_input_rejected"])
    invalid_domain_rejected = bool(initial["invalid_domain_rejected"])
    references_summary = dict(initial["references"])
    attempt.checkpoint(
        "numerical-initial", derivatives=derivative_summary,
        references=references_summary,
        conditioning_state_count=conditioning["state_count"],
    )

    def solve_budget() -> tuple[float, str]:
        remaining = max(
            0.0, total_deadline - time.monotonic() - TERMINATION_BUDGET_S,
        )
        return min(30.0, remaining), "solve" if remaining >= 30.0 else "total"

    process_cases = (
        (PRESSURE_SCALE, (20.0, 2e-4), "process-0"),
        (PRESSURE_SCALE, (60.0, 5e-4), "process-1"),
        (PRESSURE_SCALE, (90.0, 9e-4), "process-2"),
        (150000.0, (60.0, 5e-4), "process-3"),
    )
    process_nl = {
        name: prepare_process_nl(library, cap, start, args.output_dir, name)
        for cap, start, name in process_cases
    }
    first_solver_nl = process_nl["process-0"]
    whole_nlp, whole_nlp_elapsed = run_numerical_stage(
        "whole-nlp", args.output_dir, ipopt, first_solver_nl,
        max(0.0, 120.0 - numerical_elapsed), total_deadline,
    )
    numerical_elapsed += whole_nlp_elapsed
    hessian_summary = dict(whole_nlp["whole_nlp"])
    attempt.checkpoint(
        "whole-nlp", whole_nlp=hessian_summary,
        numerical_runtime_s=numerical_elapsed,
    )
    if numerical_elapsed > 120.0:
        raise RuntimeError(
            f"derivative/domain/whole-NLP checks exceeded 120 seconds: {numerical_elapsed:g}"
        )

    process_rows: list[dict[str, object]] = []
    for index, (cap, start, case_name) in enumerate(process_cases):
        timeout, timeout_scope = solve_budget()
        row = timed_process_solve(
            cap, start, index == 0, ipopt, args.output_dir, process_nl[case_name],
            case_name, timeout=timeout, deadline_scope=timeout_scope,
        )
        process_rows.append(row)
        attempt.checkpoint("process", {"process": process_rows})
        if row.get("status") != "ok" or row.get("termination") != "optimal":
            raise RuntimeError(f"required process solve failed: {row}")
        if float(row["outer_elapsed_s"]) > 30.0:
            raise RuntimeError(f"required process solve exceeded 30 seconds: {row}")

    reduced_reference = dict(initial["reduced_reference"])
    rho_reference = float(reduced_reference["rho"])
    pressure_reference = float(reduced_reference["P"])
    optimum = {
        "rho": rho_reference,
        "Q": rho_reference / DENSITY_SCALE * VOLUME_FLOW_SCALE,
        "F": rho_reference * (rho_reference / DENSITY_SCALE * VOLUME_FLOW_SCALE),
        "P": pressure_reference,
    }
    changed_cap_reference = dict(initial["changed_cap_reference"])
    changed_rho = float(changed_cap_reference["rho"])
    changed_pressure = float(changed_cap_reference["P"])
    timeout, timeout_scope = solve_budget()
    failure_nl = process_nl["process-1"]
    failure_probe = timed_process_solve(
        PRESSURE_SCALE, (60.0, 5e-4), False, ipopt, args.output_dir,
        failure_nl, "forced-post-entry-timeout", timeout=timeout, forced_timeout=True,
        deadline_scope=timeout_scope,
    )
    write_rows(args.output_dir / "failures.csv", [failure_probe])
    if (
        failure_probe.get("termination") != "forced-post-entry-timeout"
        or not failure_probe.get("entry_observed")
        or not failure_probe.get("process_group_terminated")
    ):
        raise RuntimeError(f"bounded failure probe did not prove post-entry timeout cleanup: {failure_probe}")
    attempt.checkpoint(
        "failure-probe", {"failures": [failure_probe], "process": process_rows},
        failure_probe=failure_probe,
    )

    cost_rows, cost_summary = cost_experiment(args.output_dir, ipopt, total_deadline)
    write_rows(args.output_dir / "cost.csv", cost_rows)
    attempt.checkpoint(
        "cost", {"cost": cost_rows, "failures": [failure_probe], "process": process_rows,
                 },
        cost=cost_summary,
    )
    enforce_elapsed(started, 600.0, "post-build campaign")
    slope_error = 2e-4 / float(conditioning["minimum_pressure_slope"])
    references = (
        (optimum, process_rows[:3]),
        ({"rho": changed_rho, "P": changed_pressure,
          "Q": changed_rho / DENSITY_SCALE * VOLUME_FLOW_SCALE,
          "F": changed_rho * (changed_rho / DENSITY_SCALE * VOLUME_FLOW_SCALE)},
         process_rows[3:]),
    )
    for reference, cases in references:
        reference_objective = 0.5 * (reference["Q"] / VOLUME_FLOW_SCALE) ** 2 - reference["F"] / FLOW_SCALE
        for row in cases:
            residuals = max(abs(row.get(name, math.inf)) for name in (
                "flow_balance", "volume_balance", "pressure_balance", "eos_residual",
            ))
            if (
                row.get("status") != "ok" or row.get("termination") != "optimal"
                or residuals > 1e-8 or not row.get("kkt_passed", False)
                or abs(row["rho"] / DENSITY_SCALE - reference["rho"] / DENSITY_SCALE)
                    > 1e-6 + slope_error / DENSITY_SCALE
                or abs(row["Q"] / VOLUME_FLOW_SCALE - reference["Q"] / VOLUME_FLOW_SCALE) > 1e-6
                or abs(row["F"] / FLOW_SCALE - reference["F"] / FLOW_SCALE) > 1e-6
                or abs(row["objective"] - reference_objective) > 1e-7
            ):
                raise RuntimeError(f"IDAES solution missed reduced reference: {row}")
    elapsed = time.perf_counter() - started
    write_rows(args.output_dir / "process.csv", process_rows)
    ipopt_version = str(ipopt_identity["version"])
    asl_version = str(ipopt_identity["asl_version"])
    compiler = subprocess.check_output(("c++", "--version"), text=True).splitlines()[0]
    cppad = next(
        (path.name.removeprefix("CppAD-") for path in (REPO / "engine/build/subprojects").glob("CppAD-*") if path.is_dir()),
        "unavailable",
    )
    summary = {
        "status": "bounded numerical prototype passed",
        "stage": "completed",
        "source_commit": subprocess.check_output(("git", "rev-parse", "HEAD"), cwd=REPO, text=True).strip(),
        "source_dirty": bool(subprocess.check_output(("git", "status", "--porcelain"), cwd=REPO, text=True).strip()),
        "wheel": str(wheel) if wheel else None, "wheel_sha256": wheel_hash,
        "prepared_environment": loaded_identity(),
        "asl_library": str(library), "asl_library_sha256": sha256(library),
        "asl_library_identity": library_identity,
        "adapter_source_sha256": sha256(REPO / "engine/native/bindings/asl_adapter.cpp"),
        "idaes_version": importlib.metadata.version("idaes-pse"),
        "pyomo_version": importlib.metadata.version("pyomo"),
        "ipopt": str(ipopt), "ipopt_sha256": ipopt_identity["sha256"],
        "ipopt_version": ipopt_version, "asl_version": asl_version,
        "ipopt_preflight": ipopt_identity,
        "derivative_check": {
            "requested": DERIVATIVE_TEST_OPTIONS["derivative_test"],
            "result": process_rows[0]["derivative_check"],
            "passed": process_rows[0]["derivative_check"] == "passed",
        },
        "cppad_version": cppad, "compiler": compiler,
        "build_mode": "Meson C++17 installed wheel" if wheel else "Meson C++17 prepared environment", "python": sys.version,
        "platform": platform.platform(), "parameter_strings": PARAMETER_STRINGS,
        "conditioning_grid": {
            "temperatures_k": GRID_TEMPERATURES, "density_count": 9,
            "state_count": conditioning["state_count"],
            "grid_state_count": conditioning["grid_state_count"],
            "derivative_state_count": conditioning["derivative_state_count"],
            "pressure_min_pa": conditioning["pressure_min_pa"],
            "pressure_max_pa": conditioning["pressure_max_pa"],
            "pressure_slope_over_rt_min": conditioning["pressure_slope_over_rt_min"],
        },
        "pressure_bracket": pressure_bracket,
        "immutable_input_rejected": immutable_input_rejected,
        "invalid_domain_rejected": invalid_domain_rejected,
        "reduced_reference": {
            **optimum,
            "objective": 0.5 * (optimum["Q"] / VOLUME_FLOW_SCALE) ** 2 - optimum["F"] / FLOW_SCALE,
            "pressure_error_pa": abs(optimum["P"] - PRESSURE_SCALE),
            "density_error_bound": slope_error,
        },
        "changed_cap_reference": {
            "rho": changed_rho, "pressure": changed_pressure,
            "pressure_error_pa": abs(changed_pressure - 150000.0),
        },
        "process_cases": len(process_rows), "process_accepted": True,
        "process_runtime_s": sum(float(row["outer_elapsed_s"]) for row in process_rows),
        "numerical_runtime_s": numerical_elapsed,
        "analysis_elapsed_s": elapsed,
        "references": references_summary,
        "derivatives": derivative_summary, "whole_nlp": hessian_summary,
        "whole_nlp_before_solves": True,
        "cost": cost_summary, "failure_probe": failure_probe,
        "attempt_limits": {
            "cold_load_admission_tape_s": 60.0, "warm_batch_s": 10.0,
            "derivative_domain_whole_nlp_s": 120.0, "solve_s": 30.0,
            "post_build_s": 600.0,
        },
        "solver_options": SOLVER_OPTIONS | DERIVATIVE_TEST_OPTIONS | {
            "print_level": 0, "file_print_level": 5,
        },
        "writer_options": NL_IO_OPTIONS,
        "whole_nlp_solver_files": {
            "nl": str(first_solver_nl),
            "row": str(first_solver_nl.with_suffix(".row")),
            "col": str(first_solver_nl.with_suffix(".col")),
            "nl_sha256": sha256(first_solver_nl),
        },
        "claim_limits": ["pure neutral PC-SAFT", "propane packet", "385-415 K and 5-110 mol/m3 sampled branch", "not global phase stability or application qualification"],
    }
    attempt.finish(summary)
    print(json.dumps({"status": summary["status"], "wheel_sha256": wheel_hash, "process_rows": len(process_rows), **derivative_summary, **hessian_summary}, sort_keys=True))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--solve-case", nargs=3, metavar=("CAP", "RHO", "Q"))
    parser.add_argument("--cost-mode", choices=("native", "asl"))
    parser.add_argument("--numerical-stage", choices=("initial", "whole-nlp"))
    parser.add_argument("--nl-path", type=Path)
    parser.add_argument("--derivative-test", default="False")
    parser.add_argument("--case-name", default="process")
    parser.add_argument("--entry-probe", action="store_true")
    parser.add_argument("--ipopt", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir = args.output_dir.resolve()
    if args.numerical_stage is not None:
        result = numerical_child(
            args.numerical_stage, args.output_dir, args.ipopt.resolve(), args.nl_path,
        )
        print(json.dumps(result, sort_keys=True))
        return 0
    if args.cost_mode is not None:
        _, _, library = installed_identity()
        print(json.dumps(cost_child(args.cost_mode, library), sort_keys=True))
        return 0
    if args.solve_case is not None:
        if args.nl_path is None:
            parser.error("--solve-case requires --nl-path")
        _, _, library = installed_identity()
        os.environ["AMPLFUNC"] = str(library)
        record = solve_process(
            library, float(args.solve_case[0]),
            (float(args.solve_case[1]), float(args.solve_case[2])),
            args.ipopt.resolve(), args.output_dir, args.nl_path.resolve(), args.case_name,
            args.derivative_test.lower() == "true", args.entry_probe,
        )
        print(json.dumps(record, sort_keys=True))
        return 0
    identity = campaign_identity(args.ipopt.resolve())
    attempt = CampaignAttempt(args.output_dir)
    attempt.start(identity)
    previous_handler = signal.signal(signal.SIGTERM, interrupt_campaign)
    try:
        return run_campaign(args, attempt)
    except BaseException as error:
        stage = "interrupted" if isinstance(error, KeyboardInterrupt) else "failed"
        attempt.fail(error, stage=stage)
        raise
    finally:
        signal.signal(signal.SIGTERM, previous_handler)


def declarations():
    return native_model()


if __name__ == "__main__":
    raise SystemExit(main())

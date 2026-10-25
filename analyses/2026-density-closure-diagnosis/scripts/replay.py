#!/usr/bin/env python3
"""Reproduce the bounded #89 density diagnosis on the qualified wheel.

The script intentionally keeps the numerical work finite.  It regenerates the
original control rows, checks four fixed alignments of the 256-cell scan on
those same controls, and compares the consumed pressure-density action with an
independent fixed-density pressure ladder.  It is an analysis replay, not a
replacement density solver and it never imports the Engine source tree.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import math
import time
import subprocess
from pathlib import Path
from typing import Any, Callable
from urllib.parse import unquote, urlparse

import epcsaft


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = Path(__file__).resolve()
WHEEL = ROOT / "build" / "environment-wheel" / "epcsaft-1.0.0-cp313-cp313-linux_x86_64.whl"
PACKET = ROOT / "data" / "packets" / "gross-2001-methane-ethane" / "2" / "parameters" / "parameters.json"
DECLARATION = ROOT / "analyses/2026-density-closure-diagnosis/inputs/parameters.json"
DECLARATION_SHA256 = "feb55fe933ce45d200554838e9542edc718f7220293d77ef2d533e26d4830e1e"
PACKET_SHA256 = "7dd1899712aacfbc3395d86349768b0d93dd5dfd75b2e0f0700640be618aa9c2"
SOURCE_COMMIT = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
MODEL_ID = "gross-2001-methane-ethane@2:ethane"
R = epcsaft.GAS_CONSTANT_J_PER_MOL_K

# These are the only states and alignments admitted by this replay.  The
# 16,384-cell reference is the already retained finite control reference; it
# is not a completeness claim.
CONTROLS = (
    ("one_root", 200.0, 1.0e8, "vapor"),
    ("two_root", 200.0, 217800.69494692178, "vapor"),
    ("near_coalescence", 308.8, 5147500.0, "vapor"),
    ("out_of_domain", 200.0, 1.0e10, "vapor"),
)
REFERENCE_INTERVALS = 16384
GRID_INTERVALS = 256
GRID_SHIFTS = (0.0, 0.25, 0.5, 0.75)
LOCAL_INTERVALS = 64
LOWER_FRACTION = 1.0e-8
UPPER_FRACTION = 0.999
REFERENCE_UPPER_FRACTION = 0.999999
LOG_ROOT_TOLERANCE = 1.0e-14
# Cancellation in the liquid pressure-difference check consumes a tighter root bracket.
LADDER_LOG_TOLERANCE = 5.0e-16
ACTION_RTOL = 1.0e-12
FIRST_ORDER = 5.0e-6
WALL_BUDGET_SECONDS = 60.0
DOMAIN_ROWS: list[dict[str, object]] = []


class Budget:
    def __init__(self) -> None:
        self.started = time.monotonic()

    def check(self, label: str) -> None:
        elapsed = time.monotonic() - self.started
        if elapsed > WALL_BUDGET_SECONDS:
            raise RuntimeError(f"bounded replay exceeded {WALL_BUDGET_SECONDS:g}s at {label}")

    @property
    def elapsed(self) -> float:
        return time.monotonic() - self.started


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def enum_name(value: object | None) -> str:
    if value is None:
        return ""
    return getattr(value, "name", str(value).rsplit(".", 1)[-1])


def json_field(value: object) -> str:
    return json.dumps(value, separators=(",", ":"), sort_keys=True)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    rows = [
        {
            **row,
            "replay_script": str(SCRIPT.relative_to(ROOT)),
            "replay_script_sha256": sha256(SCRIPT),
        }
        for row in rows
    ]
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def check_inputs() -> dict[str, object]:
    if not WHEEL.is_file() or not PACKET.is_file():
        raise RuntimeError(f"missing qualified input: wheel={WHEEL} packet={PACKET}")
    observed = {"wheel": sha256(WHEEL), "packet": sha256(PACKET)}
    if observed["packet"] != PACKET_SHA256:
        raise RuntimeError(f"parameter packet hash mismatch: {observed}")

    subprocess.run([str(ROOT / ".venv/bin/python"), str(ROOT / "tools/verify-installed-wheel.py"),
                    "--wheel", str(WHEEL)], check=True)
    if sha256(DECLARATION) != DECLARATION_SHA256:
        raise RuntimeError("neutral declaration hash mismatch")
    module_path = Path(epcsaft.__file__).resolve()
    if ROOT / ".venv" not in module_path.parents:
        raise RuntimeError(f"epcsaft is not imported from the qualified .venv: {module_path}")
    distribution = importlib.metadata.distribution("epcsaft")
    if distribution.version != "1.0.0":
        raise RuntimeError(f"unexpected installed epcsaft version: {distribution.version}")
    direct_url_text = distribution.read_text("direct_url.json")
    if direct_url_text is None:
        raise RuntimeError("installed epcsaft has no direct_url.json")
    direct_url = json.loads(direct_url_text)
    parsed = urlparse(direct_url.get("url", ""))
    installed_wheel = Path(unquote(parsed.path)).resolve() if parsed.scheme == "file" else None
    if installed_wheel != WHEEL.resolve():
        raise RuntimeError(f"installed wheel provenance mismatch: {installed_wheel} != {WHEEL}")
    return {
        "wheel": str(WHEEL.relative_to(ROOT)),
        "wheel_sha256": observed["wheel"],
        "parameter_packet": str(PACKET.relative_to(ROOT)),
        "parameter_packet_sha256": observed["packet"],
        "neutral_declaration": str(DECLARATION.relative_to(ROOT)),
        "neutral_declaration_sha256": DECLARATION_SHA256,
        "installed_distribution": distribution.version,
        "installed_module": str(module_path),
        "installed_wheel": str(installed_wheel),
    }


def finite_state(model: epcsaft.Pure, temperature: float, density: float, budget: Budget) -> Any | None:
    budget.check("fixed-density state")
    try:
        return model.state(temperature, rho=density)
    except Exception:
        return None


def pressure_residual(
    model: epcsaft.Pure,
    temperature: float,
    target_pressure: float,
    density: float,
    budget: Budget,
) -> tuple[float, Any] | None:
    state = finite_state(model, temperature, density, budget)
    if state is None:
        return None
    return state.pressure - target_pressure, state


def admitted_interval(
    model: epcsaft.Pure,
    temperature: float,
    target_pressure: float,
    budget: Budget,
) -> tuple[float, float, float]:
    ideal_density = target_pressure / (R * temperature)
    lower = ideal_density * LOWER_FRACTION
    anchor = max(1.0e-3, lower)
    diagnostic = model.density_diagnostics(
        temperature, target_pressure, nodes=(anchor, anchor * 2), adaptive=False,
        density_tolerance=1e-9, log_tolerance=LOG_ROOT_TOLERANCE,
        max_evaluations=100000, max_seconds=max(0.001, WALL_BUDGET_SECONDS - budget.elapsed),
        finite_anchor=anchor, finite_upper=model.molar_density_upper_bound(temperature) * 2,
        evaluation_order=epcsaft._core.DerivativeOrder.Hessian,
    )
    boundary = diagnostic.finite_domain
    if boundary is None or boundary.stop != "finite_to_unavailable":
        raise RuntimeError(f"finite-domain search unresolved: {diagnostic.stop}")
    DOMAIN_ROWS.append({"temperature_K": temperature, "target_pressure_Pa": target_pressure,
                        "last_finite_mol_m3": boundary.lower, "first_unavailable_mol_m3": boundary.upper,
                        "geometric_cap_mol_m3": diagnostic.geometric_cap,
                        "evaluation_kind": "Hessian", "cause": boundary.cause, "stop": boundary.stop})
    return lower, boundary.lower * UPPER_FRACTION, boundary.lower * REFERENCE_UPPER_FRACTION


def scan_roots(
    model: epcsaft.Pure, temperature: float, target_pressure: float,
    lower: float, upper: float, intervals: int, budget: Budget, *, shift: float = 0.0, adaptive: bool = False,
    log_tolerance: float = LOG_ROOT_TOLERANCE,
) -> tuple[list[dict[str, float]], list[dict[str, float]]]:
    budget.check("declared native density grid")
    log_lower, log_upper = math.log(lower), math.log(upper)
    nodes = [lower, *(math.exp(log_lower + (index + shift) / intervals * (log_upper - log_lower))
                      for index in range(1, intervals)), upper]
    diagnostic = model.density_diagnostics(
        temperature, target_pressure, nodes=nodes, log_tolerance=log_tolerance,
        density_tolerance=1e-9, max_evaluations=100000, finite_upper=0.0,
        evaluation_order=epcsaft._core.DerivativeOrder.PressureSlope,
        adaptive=adaptive, max_seconds=max(0.001, WALL_BUDGET_SECONDS - budget.elapsed),
    )
    if diagnostic.stop not in ("resolved_crossings", "no_sampled_crossing"):
        raise RuntimeError(f"declared density grid unresolved: {diagnostic.stop}")
    roots = [{"density": row.molar_density, "pressure": row.pressure,
              "pressure_residual": row.residual, "pressure_slope": row.pressure_density_derivative,
              "stable": row.pressure_density_derivative > 0.0,
              "bracket_lower": row.lower, "bracket_upper": row.upper,
              "log_bracket_width": math.log(row.upper / row.lower), "stop": row.stop}
             for row in diagnostic.pressure_roots]
    assert all(row["log_bracket_width"] <= log_tolerance or row["stop"] == "exact_numerical_zero" for row in roots)
    brackets = [{"lower_density": row.lower, "upper_density": row.upper,
                 "lower_residual_Pa": row.lower_residual, "upper_residual_Pa": row.upper_residual}
                for row in diagnostic.declared_brackets]
    return roots, brackets


def reference_roots(
    model: epcsaft.Pure,
    temperature: float,
    target_pressure: float,
    lower: float,
    upper: float,
    budget: Budget,
) -> tuple[list[dict[str, float]], list[dict[str, float]]]:
    return scan_roots(
        model, temperature, target_pressure, lower, upper,
        REFERENCE_INTERVALS, budget,
    )


def current_result(
    model: epcsaft.Pure,
    temperature: float,
    target_pressure: float,
    phase: str,
    budget: Budget,
) -> dict[str, object]:
    budget.check("native cold density closure")
    try:
        state = model.state(temperature, P=target_pressure, phase=phase)
    except Exception as error:
        return {
            "status": "density_unavailable",
            "density": "",
            "branch": "",
            "stable_root_count": "",
            "reduced_residual": "",
            "iterations": "",
            "message": str(error),
        }
    convergence = state.density_convergence
    return {
        "status": "available",
        "density": state.molar_density,
        "branch": enum_name(state.density_branch),
        "stable_root_count": state.stable_root_count,
        "reduced_residual": "" if convergence is None else convergence.residual,
        "iterations": "" if convergence is None else convergence.iterations,
        "message": "",
    }


def reference_role(index: int, roots: list[dict[str, float]]) -> str:
    if not roots:
        return ""
    if not roots[index]["stable"]:
        return "unstable"
    stable = [root for root in roots if root["stable"]]
    if len(stable) == 1:
        return "stable"
    return "stable_vapor" if roots[index] is stable[0] else "stable_liquid"


def controls_and_domains(
    model: epcsaft.Pure,
    budget: Budget,
) -> tuple[list[dict[str, object]], dict[str, tuple[float, float, float, list[dict[str, float]], list[dict[str, float]]]]]:
    rows: list[dict[str, object]] = []
    domains: dict[str, tuple[float, float, float, list[dict[str, float]], list[dict[str, float]]]] = {}
    for case, temperature, target_pressure, phase in CONTROLS:
        lower, native_upper, reference_upper = admitted_interval(
            model, temperature, target_pressure, budget,
        )
        roots, _ = reference_roots(
            model, temperature, target_pressure, lower, reference_upper, budget,
        )
        current_phase = phase
        if case == "two_root":
            current_phase = "vapor"
        current_by_phase = {current_phase: current_result(
            model, temperature, target_pressure, current_phase, budget,
        )}
        if case == "two_root":
            current_by_phase["liquid"] = current_result(
                model, temperature, target_pressure, "liquid", budget,
            )
        domains[case] = (lower, native_upper, reference_upper, roots, [])
        if roots:
            for index, root in enumerate(roots):
                # The retained near-coalescence control deliberately records
                # one cold vapor request for every reference root.  The native
                # cold scan therefore shows the same low root for both stable
                # reference branches.  The regular two-root control retains
                # its explicit liquid request for the high branch.
                row_phase = (
                    "liquid"
                    if case == "two_root" and reference_role(index, roots) == "stable_liquid"
                    else current_phase
                )
                current = current_by_phase.setdefault(row_phase, current_result(
                    model, temperature, target_pressure, row_phase, budget,
                ))
                rows.append({
                    "case": case,
                    "model_id": MODEL_ID,
                    "parameter_packet_sha256": PACKET_SHA256,
                    "temperature_K": temperature,
                    "target_pressure_Pa": target_pressure,
                    "composition": "[1.0]",
                    "native_domain_upper_mol_m3": native_upper,
                    "reference_domain_upper_mol_m3": reference_upper,
                    "reference_root_index": index,
                    "reference_root_role": reference_role(index, roots),
                    "reference_density_mol_m3": root["density"],
                    "reference_pressure_Pa": root["pressure"],
                    "reference_bracket_lower_mol_m3": root["bracket_lower"],
                    "reference_bracket_upper_mol_m3": root["bracket_upper"],
                    "reference_log_bracket_width": root["log_bracket_width"],
                    "reference_stop": root["stop"],
                    "reference_pressure_residual_Pa": root["pressure_residual"],
                    "reference_pressure_slope_Pa_m3_per_mol": root["pressure_slope"],
                    "reference_stable": str(root["stable"]).lower(),
                    "current_phase": row_phase,
                    "current_status": current["status"],
                    "current_message": current["message"],
                    "current_density_mol_m3": current["density"],
                    "current_branch": current["branch"],
                    "current_stable_root_count": current["stable_root_count"],
                    "current_reduced_residual": current["reduced_residual"],
                    "current_iterations": current["iterations"],
                })
        else:
            current = current_by_phase[current_phase]
            rows.append({
                "case": case,
                "model_id": MODEL_ID,
                "parameter_packet_sha256": PACKET_SHA256,
                "temperature_K": temperature,
                "target_pressure_Pa": target_pressure,
                "composition": "[1.0]",
                "native_domain_upper_mol_m3": native_upper,
                "reference_domain_upper_mol_m3": reference_upper,
                "reference_root_index": "",
                "reference_root_role": "",
                "reference_density_mol_m3": "",
                "reference_pressure_Pa": "",
                "reference_pressure_residual_Pa": "",
                "reference_pressure_slope_Pa_m3_per_mol": "",
                "reference_stable": "",
                "current_phase": current_phase,
                "current_status": current["status"],
                "current_message": current["message"],
                "current_density_mol_m3": current["density"],
                "current_branch": current["branch"],
                "current_stable_root_count": current["stable_root_count"],
                "current_reduced_residual": current["reduced_residual"],
                "current_iterations": current["iterations"],
            })
    return rows, domains


def grid_rows(
    model: epcsaft.Pure,
    domains: dict[str, tuple[float, float, float, list[dict[str, float]], list[dict[str, float]]]],
    budget: Budget,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case, temperature, target_pressure, phase in CONTROLS:
        lower, native_upper, reference_upper, references, _ = domains[case]
        current = current_result(model, temperature, target_pressure, phase, budget)
        for shift in GRID_SHIFTS:
            # Four shifts are frozen before this loop; no count search is added.
            # The endpoint must match density.cpp's 0.999 close-packing bound.
            roots, brackets = scan_roots(
                model, temperature, target_pressure, lower, native_upper,
                GRID_INTERVALS, budget, shift=shift,
            )
            adaptive_roots, adaptive_brackets = scan_roots(
                model, temperature, target_pressure, lower, native_upper,
                GRID_INTERVALS, budget, shift=shift, adaptive=True,
            )
            assert len(adaptive_brackets) == len(brackets)
            reference_crossings = len(references)
            stable_reference = sum(1 for root in references if root["stable"])
            gap_logs = [
                math.log(references[index + 1]["density"] / references[index]["density"])
                for index in range(len(references) - 1)
            ]
            grid_step = math.log(native_upper / lower) / GRID_INTERVALS
            if len(brackets) < reference_crossings:
                diagnosis = "missed_reference_crossing_aliasing"
            elif len(brackets) == reference_crossings:
                diagnosis = "all_reference_crossings_bracketed"
            else:
                diagnosis = "extra_sign_change_requires_review"
            rows.append({
                "case": case,
                "temperature_K": temperature,
                "target_pressure_Pa": target_pressure,
                "composition": "[1.0]",
                "grid_intervals": GRID_INTERVALS,
                "shift_fraction": shift,
                "endpoints_included": True,
                "admitted_lower_mol_m3": lower,
                "native_admitted_upper_mol_m3": native_upper,
                "reference_admitted_upper_mol_m3": reference_upper,
                "lower_endpoint_residual_Pa": pressure_residual(
                    model, temperature, target_pressure, lower, budget,
                )[0],
                "upper_endpoint_residual_Pa": pressure_residual(
                    model, temperature, target_pressure, native_upper, budget,
                )[0],
                "grid_log_step": grid_step,
                "reference_root_count": reference_crossings,
                "reference_stable_root_count": stable_reference,
                "reference_root_gaps_log": json_field(gap_logs),
                "reference_min_root_gap_log": min(gap_logs) if gap_logs else "",
                "grid_bracket_count": len(brackets),
                "grid_brackets": json_field(brackets),
                "grid_refined_root_count": len(roots),
                "adaptive_refined_root_count": len(adaptive_roots),
                "current_phase": phase,
                "current_status": current["status"],
                "current_message": current["message"],
                "current_stable_root_count": current["stable_root_count"],
                "aliasing_diagnosis": diagnosis,
            })
    return rows


def local_root(
    model: epcsaft.Pure,
    temperature: float,
    target_pressure: float,
    reference_density: float,
    lower_domain: float,
    upper_domain: float,
    budget: Budget,
) -> float:
    # A fixed local ladder around one independently identified stable branch;
    # it is deliberately separate from the cold native scan.
    window = 0.025
    lower = max(lower_domain, reference_density * math.exp(-window))
    upper = min(upper_domain, reference_density * math.exp(window))
    roots, _ = scan_roots(
        model, temperature, target_pressure, lower, upper,
        LOCAL_INTERVALS, budget, log_tolerance=LADDER_LOG_TOLERANCE,
    )
    if not roots:
        raise RuntimeError(
            f"local pressure ladder could not bracket branch near rho={reference_density:g}"
        )
    candidate = min(roots, key=lambda root: abs(math.log(root["density"] / reference_density)))
    if not candidate["stable"]:
        raise RuntimeError(
            f"local pressure ladder selected a nonstable root near rho={reference_density:g}"
        )
    return candidate["density"]


def action_rows(
    model: epcsaft.Pure,
    domains: dict[str, tuple[float, float, float, list[dict[str, float]], list[dict[str, float]]]],
    budget: Budget,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    direction = epcsaft.PressureActionDirection(0.0, 1.0, (0.0,), ())
    for case, temperature, target_pressure, default_phase in CONTROLS:
        if case == "out_of_domain":
            continue
        lower, _, reference_upper, references, _ = domains[case]
        for index, root in enumerate(references):
            if not root["stable"]:
                continue
            role = reference_role(index, references)
            phase = "liquid" if role == "stable_liquid" else "vapor"
            reference_density = root["density"]
            # The reference reciprocal is for the independent finite-difference
            # ladder. Native identity is evaluated at the density native code
            # actually returned, below.
            expected = 1.0 / root["pressure_slope"]
            # The independent pressure ladder uses a relative 1e-6 step.  It
            # is small enough to test the smooth local branch near the frozen
            # near-coalescence roots while remaining well above double
            # precision noise for these SI pressures.
            h = target_pressure * 1.0e-6
            rho_minus = local_root(model, temperature, target_pressure - h, reference_density, lower, reference_upper, budget)
            rho_plus = local_root(model, temperature, target_pressure + h, reference_density, lower, reference_upper, budget)
            rho_minus_half = local_root(model, temperature, target_pressure - h / 2.0, reference_density, lower, reference_upper, budget)
            rho_plus_half = local_root(model, temperature, target_pressure + h / 2.0, reference_density, lower, reference_upper, budget)
            coarse = (rho_plus - rho_minus) / (2.0 * h)
            fine = (rho_plus_half - rho_minus_half) / h
            coarse_error = abs(coarse - expected) / abs(expected)
            fine_error = abs(fine - expected) / abs(expected)
            fine_pass = fine_error <= FIRST_ORDER
            coarse_pass = coarse_error <= 4.0 * FIRST_ORDER

            budget.check("native pressure-density action")
            native_status = ""
            native_density: object = ""
            native_pressure_slope: object = ""
            native_action: object = ""
            native_message = ""
            try:
                batch = model.pressure_actions(
                    temperature, target_pressure, phase=phase, observables=(), directions=(direction,)
                )
                native_status = enum_name(batch.density_status)
                native_density = batch.molar_density if batch.molar_density is not None else ""
                native_action = batch.density_action[0] if batch.density_action else ""
                native_message = batch.density_message
                if native_density != "":
                    try:
                        native_state = model.state(temperature, rho=float(native_density))
                        native_pressure_slope = native_state.pressure_density_derivative
                    except Exception:
                        pass
            except Exception as error:
                native_status = "exception"
                native_message = str(error)

            target_match = (
                native_density != ""
                and abs(float(native_density) - reference_density)
                <= max(1.0e-8, 1.0e-10 * abs(reference_density))
            )
            native_expected: object = ""
            if native_pressure_slope != "":
                native_expected = 1.0 / float(native_pressure_slope)
            identity_error: object = ""
            identity_pass: object = ""
            if native_action != "" and native_expected != "":
                identity_error = abs(float(native_action) - float(native_expected)) / max(
                    abs(float(native_action)), abs(float(native_expected)), 1.0e-300,
                )
                identity_pass = identity_error <= ACTION_RTOL
            if target_match:
                current_branch_status = "target_branch_returned"
            else:
                current_branch_status = "missed_by_current_cold_scan"
            rows.append({
                "case": case,
                "stable_root_role": role,
                "reference_root_index": index,
                "temperature_K": temperature,
                "target_pressure_Pa": target_pressure,
                "branch": phase,
                "reference_density_mol_m3": reference_density,
                "reference_pressure_slope_Pa_m3_per_mol": root["pressure_slope"],
                "expected_d_rho_dP_mol_m3_per_Pa": expected,
                "ladder_step_Pa": h,
                "rho_minus_h_mol_m3": rho_minus,
                "rho_plus_h_mol_m3": rho_plus,
                "rho_minus_half_h_mol_m3": rho_minus_half,
                "rho_plus_half_h_mol_m3": rho_plus_half,
                "coarse_d_rho_dP_mol_m3_per_Pa": coarse,
                "fine_d_rho_dP_mol_m3_per_Pa": fine,
                "coarse_normalized_error": coarse_error,
                "fine_normalized_error": fine_error,
                "coarse_ladder_pass_4x5e-6": coarse_pass,
                "fine_ladder_pass_5e-6": fine_pass,
                "native_status": native_status,
                "native_density_mol_m3": native_density,
                "native_pressure_slope_Pa_m3_per_mol": native_pressure_slope,
                "native_expected_d_rho_dP_mol_m3_per_Pa": native_expected,
                "native_density_action_mol_m3_per_Pa": native_action,
                "native_message": native_message,
                "current_branch_status": current_branch_status,
                "identity_relative_error": identity_error,
                "identity_pass_rtol_1e-12": identity_pass,
            })
    return rows


def anchor_rows(model: epcsaft.Pure, budget: Budget) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    branches = (
        ("vapor", epcsaft.DensityBranch.Vapor),
        ("liquid", epcsaft.DensityBranch.Liquid),
    )
    for phase, branch in branches:
        budget.check("positive anchor")
        cold = model._model.state_at_pressure(280.0, 1.0e6, (1.0,), branch)
        anchored = model._model.state_at_pressure(
            280.0, 1.0e6, (1.0,), branch, cold.molar_density,
        )
        rows.append({
            "control": "positive_replay",
            "temperature_K": 280.0,
            "pressure_Pa": 1.0e6,
            "phase": phase,
            "anchor_density_mol_m3": cold.molar_density,
            "replayed_density_mol_m3": anchored.molar_density,
            "density_relative_error": abs(anchored.molar_density - cold.molar_density) / cold.molar_density,
            "status": "available",
            "branch": enum_name(anchored.density_branch),
            "message": "",
        })

    for control, anchor_density in (("negative_anchor", -1.0), ("finite_positive_outside_domain", 1.0e-6)):
        budget.check(control)
        try:
            model._model.state_at_pressure(
                200.0, 217800.69494692178, (1.0,), epcsaft.DensityBranch.Liquid, anchor_density,
            )
        except Exception as error:
            rows.append({
                "control": control,
                "temperature_K": 200.0,
                "pressure_Pa": 217800.69494692178,
                "phase": "liquid",
                "anchor_density_mol_m3": anchor_density,
                "replayed_density_mol_m3": "",
                "density_relative_error": "",
                "status": "rejected",
                "branch": "",
                "message": f"{type(error).__name__}: {error}",
            })
        else:
            raise RuntimeError(f"anchor control unexpectedly succeeded: {control}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output_dir = args.output_dir if args.output_dir.is_absolute() else ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    budget = Budget()
    provenance = check_inputs()
    model = epcsaft.Pure(epcsaft.Parameters.from_json(DECLARATION, components=("ethane",)))
    root_rows, domains = controls_and_domains(model, budget)
    write_csv(output_dir / "root-controls.csv", root_rows)
    previous = list(csv.DictReader((DECLARATION.parent / "root-controls-before.csv").open()))
    comparisons = []
    for old, new in zip(previous, root_rows, strict=True):
        assert old["case"] == new["case"] and str(new["reference_root_index"]) == old["reference_root_index"]
        for quantity in ("reference_density_mol_m3", "reference_pressure_Pa",
                         "reference_pressure_slope_Pa_m3_per_mol", "native_domain_upper_mol_m3",
                         "reference_domain_upper_mol_m3", "current_density_mol_m3", "current_stable_root_count"):
            before, after = old[quantity], new[quantity]
            comparisons.append({"case": new["case"], "root_index": new["reference_root_index"],
                                "quantity": quantity, "before": before, "after": after,
                                "difference": float(after) - float(before) if before != "" and after != "" else ""})
        if old["reference_density_mol_m3"]:
            assert abs(float(new["reference_density_mol_m3"]) / float(old["reference_density_mol_m3"]) - 1) <= 1e-7
            assert str(new["reference_stable"]) == old["reference_stable"]
    write_csv(output_dir / "before-after.csv", comparisons)
    write_csv(output_dir / "finite-domain.csv", DOMAIN_ROWS)
    grid = grid_rows(model, domains, budget)
    write_csv(output_dir / "scan-grids.csv", grid)
    actions = action_rows(model, domains, budget)
    write_csv(output_dir / "pressure-actions.csv", actions)
    assert [int(row["grid_bracket_count"]) for row in grid if row["case"] == "near_coalescence"] == [3, 1, 1, 3]
    assert [row["adaptive_refined_root_count"] for row in grid if row["case"] == "near_coalescence"] == [3, 3, 3, 3]
    assert [len(domains[name][3]) for name, *_ in CONTROLS] == [1, 3, 3, 0]
    for actual, expected in zip(domains["two_root"][3], (138.1005458522085, 3218.7215247759887, 17541.25129083695)):
        assert abs(actual["density"] - expected) <= 1e-9
    assert all(row["identity_pass_rtol_1e-12"] and row["fine_ladder_pass_5e-6"] and row["coarse_ladder_pass_4x5e-6"] for row in actions)
    anchors = anchor_rows(model, budget)
    write_csv(output_dir / "anchors.csv", anchors)
    assert max(float(row["density_relative_error"]) for row in anchors if row["status"] == "available") <= 1e-10
    coefficients = {row["family"]: row["value"]["magnitude"] for row in json.loads(PACKET.read_text())["components"][0]["coefficients"]}
    limits = []
    for T in (200.0, 308.8):
        diameter = coefficients["segment_diameter"] * (1 - .12 * math.exp(-3 * coefficients["dispersion_energy_over_k"] / T))
        geometric_reference = (math.pi / (3 * math.sqrt(2))) / (math.pi / 6 * 6.02214076e23 * 1e-30 * coefficients["segment_count"] * diameter**3)
        boundary = next(row for row in DOMAIN_ROWS if row["temperature_K"] == T)
        assert float(boundary["last_finite_mol_m3"]) <= geometric_reference * (1 + 1e-13)
        assert float(boundary["first_unavailable_mol_m3"]) >= geometric_reference * (1 - 1e-13)
        assert abs(float(boundary["geometric_cap_mol_m3"]) / geometric_reference - 1) <= 1e-13
        errors = []
        for rho in (1e-6, 5e-7):
            state = model.state(T, rho=rho)
            z_error = abs(state.pressure / (rho * R * T) - 1)
            slope_error = abs(state.pressure_density_derivative / (R * T) - 1)
            errors.append((z_error, slope_error))
            limits.append({"temperature_K": T, "density_mol_m3": rho, "pressure_over_rho_RT_error": z_error,
                           "slope_over_RT_error": slope_error, "independent_geometric_cap_mol_m3": geometric_reference})
        assert max(errors[0]) <= 1e-6 and all(b < a for a, b in zip(*errors))
    write_csv(output_dir / "dilute-limit.csv", limits)
    budget.check("replay completion")
    summary = {
        "schema_version": 1,
        "analysis": "2026-density-closure-diagnosis",
        "issue": 89,
        "source_commit": SOURCE_COMMIT,
        "provenance": provenance,
        "replay_script": str(SCRIPT.relative_to(ROOT)),
        "replay_script_sha256": sha256(SCRIPT),
        "frozen_probe": {
            "controls": [case for case, _, _, _ in CONTROLS],
            "grid_intervals": GRID_INTERVALS,
            "grid_shifts": list(GRID_SHIFTS),
            "reference_intervals": REFERENCE_INTERVALS,
            "native_upper_fraction": UPPER_FRACTION,
            "reference_upper_fraction": REFERENCE_UPPER_FRACTION,
            "local_ladder_intervals": LOCAL_INTERVALS,
            "wall_budget_seconds": WALL_BUDGET_SECONDS,
        },
        "criteria": {
            "root_log_bracket_tolerance": LOG_ROOT_TOLERANCE,
            "verification_ladder_log_bracket_tolerance": LADDER_LOG_TOLERANCE,
            "ladder_error_denominator": "abs(expected_d_rho_dP), not max(1, abs(expected_d_rho_dP))",
            "action_identity_rtol": ACTION_RTOL,
            "ladder_fine_rtol": FIRST_ORDER,
            "ladder_coarse_rtol": 4.0 * FIRST_ORDER,
            "physical_prediction_qualification": "not established by this diagnosis",
        },
        "claim_limits": [
            "The independent finite root ladder supports the retained controls only; it is not global completeness or phase discovery.",
            "A failed cold numerical closure is preserved as a solver/domain result and is not evidence of physical nonexistence.",
            "The four shifted 256-cell alignments measure declared-grid crossings separately from adaptive discoveries.",
        ],
        "historical_reference_commit": "9c0deb53",
        "historical_reference_sha256": {name: sha256(DECLARATION.parent / name) for name in
                                        ("root-controls-before.csv", "pressure-actions-before.csv")},
        "outputs": ["ladder-refinement.csv", "before-after.csv", "dilute-limit.csv", "finite-domain.csv", "root-controls.csv", "scan-grids.csv", "pressure-actions.csv", "anchors.csv"],
        "elapsed_seconds": budget.elapsed,
    }
    (output_dir / "replay-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"output_dir": str(output_dir), "elapsed_seconds": budget.elapsed, "rows": {
        "root_controls": len(root_rows), "scan_grids": len(grid),
        "pressure_actions": len(actions), "anchors": len(anchors),
    }}, sort_keys=True))


def declarations():
    from epcsaft import Parameters, Pure
    return Pure(Parameters.from_json(DECLARATION, components=('ethane',)))


if __name__ == "__main__":
    main()

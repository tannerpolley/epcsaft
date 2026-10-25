#!/usr/bin/env python3
"""Bounded #90 initialization and continuation diagnosis on the installed wheel."""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "analyses" / "2026-initialization-diagnosis" / "results"
WHEEL = ROOT / "build" / "environment-wheel" / "epcsaft-0.2.0.dev0-cp313-cp313-linux_x86_64.whl"
WHEEL_SHA256 = "9f73aeb466c54eff80c89a8776c90d49cb9c9c2d130f0eb15bf0311c1f76e606"
PACKET = ROOT / "data" / "packets" / "gross-2001-methane-ethane" / "2" / "parameters" / "parameters.json"
PACKET_SHA256 = "7dd1899712aacfbc3395d86349768b0d93dd5dfd75b2e0f0700640be618aa9c2"
SYNTHETIC_SOURCE = ROOT / "engine" / "tests" / "equilibrium" / "test_native_transport.py"

sys.path.insert(0, str(ROOT / "engine"))

from epcsaft import Mixture, Parameters
from epcsaft.equilibrium import (  # noqa: E402
    Amounts,
    Continuation,
    Free,
    Phase,
    Pinned,
    Problem,
    continue_equilibrium,
    solve_equilibrium,
)
from tests.equilibrium.test_native_transport import (  # noqa: E402
    correlation,
    reactive_model,
    reactive_problem,
)

TARGET_T = 180.0
TARGET_P = 1596859.1113026983
ROW_TOLERANCE = 1.0e-10
TARGET_T_ATOL = 1.0e-6
TARGET_P_RTOL = 1.0e-8
# These are the existing E4 branch-comparison criteria. The density window is
# the approved P2 local-envelope value, not a new acceptance tolerance.
BRANCH_PRESSURE_RTOL = 1.0e-8
BRANCH_COMPOSITION_ATOL = 1.0e-8
BRANCH_AMOUNT_ATOL_MOL = 1.0e-8
BRANCH_DENSITY_LOG_WINDOW = 0.009950330853168083
CONTINUATION = {"initial_step": 0.25, "minimum_step": 0.25, "max_steps": 4, "max_retries_per_step": 1}
LONG_CONTINUATION = {"initial_step": 0.25, "minimum_step": 0.25, "max_steps": 8, "max_retries_per_step": 1}

# Freeze the matrix here before executing any retained rows. A cold free-T start is
# diagnostic: a typed topology event is preserved and is never counted as physical nonexistence.
COLD_GUESSES = (170.0, 179.0, 180.0, 181.0, 190.0)
WARM_SEEDS = (179.0, 181.0)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def enum_name(value: object | None) -> str:
    if value is None:
        return ""
    return str(value).rsplit(".", 1)[-1]


def finite_max(values: list[float]) -> float | None:
    return max((abs(value) for value in values), default=None)


def convergence_fields(result: Any) -> dict[str, object]:
    convergence = result.convergence
    criterion = convergence.criterion
    return {
        "convergence_stop": enum_name(convergence.stop),
        "residual_inf": finite_max(result.residuals),
        "classified_residual_inf": finite_max(result.classified_residuals),
        "convergence_residual": convergence.residual,
        "criterion_name": None if criterion is None else criterion.name,
        "criterion_value": None if criterion is None else criterion.value,
        "floor": convergence.floor,
    }


def encoded(values: object) -> str:
    return json.dumps(list(values), separators=(",", ":"))


def result_fields(result: Any, phase_order: tuple[str, ...] | None = None) -> dict[str, object]:
    fields = {
        "success": result.success,
        "solver_status": result.solver_status,
        "message": result.message,
        "topology_event": enum_name(result.topology_event),
        "topology_phase": result.topology_phase,
        "topology_pair": None if result.topology_pair is None else list(result.topology_pair),
        "phase_order": None if phase_order is None else ",".join(phase_order),
        "amounts_mol": encoded(result.amounts),
        "molar_densities_mol_m3": encoded(result.molar_densities),
        "mole_fractions": encoded(result.mole_fractions),
        "temperature_K": result.temperature,
        "pressure_Pa": result.pressure,
        "iterations": result.iterations,
        "evaluations": result.evaluations,
        "hessian_calls": result.hessian_calls,
        "requested_tolerance_met": result.requested_tolerance_met,
        "phase_count": result.phase_count,
    }
    fields.update(convergence_fields(result))
    return fields


def bubble_boundary(model: Mixture, temperature: float, pressure: float | Free) -> Problem:
    return Problem(
        phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))],
        T=temperature,
        P=pressure,
        feed=Amounts({"methane": 0.5, "ethane": 0.5}),
    )


def free_temperature(model: Mixture, guess: float) -> Problem:
    return bubble_boundary(model, Free(guess), TARGET_P)


def accepted(result: Any) -> bool:
    return (
        result.success
        and result.requested_tolerance_met is True
        and enum_name(result.topology_event) == "None_"
        and result.convergence.residual <= ROW_TOLERANCE
    )


def comparison_fields(
    reference: Any | None,
    result: Any | None,
    phase_order: tuple[str, ...] | None,
    reason: str | None = None,
) -> dict[str, object]:
    fields: dict[str, object] = {
        "anchor_match": None,
        "anchor_pressure_rel_error": None,
        "anchor_composition_abs_error": None,
        "anchor_density_log_error": None,
        "anchor_amount_abs_error_mol": None,
        "anchor_compare_reason": reason,
    }
    if reason is not None or reference is None or result is None:
        return fields
    if not result.success or result.temperature is None or result.pressure is None:
        fields["anchor_compare_reason"] = "state_unavailable_after_typed_failure"
        return fields
    if phase_order is None or len(result.amounts) != len(reference.amounts):
        fields["anchor_compare_reason"] = "layout_unavailable"
        return fields
    if (
        len(result.molar_densities) != len(reference.molar_densities)
        or len(result.mole_fractions) != len(reference.mole_fractions)
    ):
        fields["anchor_compare_reason"] = "state_unavailable"
        return fields
    pressure_error = abs(result.pressure / reference.pressure - 1.0)
    composition_error = max(
        (abs(left - right) for left, right in zip(result.mole_fractions, reference.mole_fractions)),
        default=math.inf,
    )
    amount_error = max(
        (abs(left - right) for left, right in zip(result.amounts, reference.amounts)),
        default=math.inf,
    )
    density_error = max(
        (abs(math.log(left / right)) for left, right in zip(result.molar_densities, reference.molar_densities)
         if left > 0.0 and right > 0.0),
        default=math.inf,
    )
    fields.update({
        "anchor_pressure_rel_error": pressure_error,
        "anchor_composition_abs_error": composition_error,
        "anchor_density_log_error": density_error,
        "anchor_amount_abs_error_mol": amount_error,
        "anchor_match": (
            result.success
            and enum_name(result.topology_event) == "None_"
            and pressure_error <= BRANCH_PRESSURE_RTOL
            and composition_error <= BRANCH_COMPOSITION_ATOL
            and density_error <= BRANCH_DENSITY_LOG_WINDOW
            and amount_error <= BRANCH_AMOUNT_ATOL_MOL
        ),
        "anchor_compare_reason": "compared_to_fixed_180_anchor",
    })
    return fields


def direct_row(
    case_id: str,
    family: str,
    expected: str,
    result: Any,
    *,
    phase_order: tuple[str, ...] | None = None,
    anchor: Any | None = None,
    **inputs: object,
) -> dict[str, object]:
    target_temperature = inputs.get("target_temperature_K")
    target_pressure = inputs.get("target_pressure_Pa")
    target_ok = True
    if isinstance(target_temperature, (int, float)) and result.temperature is not None:
        target_ok = abs(result.temperature - target_temperature) <= TARGET_T_ATOL
    if isinstance(target_pressure, (int, float)) and result.pressure is not None:
        target_ok = target_ok and abs(result.pressure / target_pressure - 1.0) <= TARGET_P_RTOL
    row = {
        "case_id": case_id,
        "family": family,
        "expected": expected,
        **inputs,
        **result_fields(result, phase_order),
        "accepted_by_matrix": accepted(result) and target_ok,
    }
    row.update(comparison_fields(anchor, result, phase_order, "not_requested" if anchor is None else None))
    return row


def continuation_row(
    case_id: str,
    direction: str,
    expected: str,
    result: Any,
    *,
    phase_order: tuple[str, ...] | None = None,
    anchor: Any | None = None,
    **inputs: object,
) -> dict[str, object]:
    row: dict[str, object] = {
        "case_id": case_id,
        "family": "continuation",
        "expected": expected,
        "direction": direction,
        "phase_order": None if phase_order is None else ",".join(phase_order),
        **inputs,
        "success": result.success,
        "solver_status": "",
        "message": result.message,
        "topology_event": "",
        "topology_phase": result.event_phase,
        "topology_pair": None if result.event_pair is None else list(result.event_pair),
        "temperature_K": None,
        "pressure_Pa": None,
        "iterations": None,
        "evaluations": None,
        "hessian_calls": None,
        "requested_tolerance_met": None,
        "phase_count": None,
        "convergence_stop": "",
        "residual_inf": None,
        "classified_residual_inf": None,
        "convergence_residual": None,
        "criterion_name": None,
        "criterion_value": None,
        "floor": None,
        "accepted_by_matrix": result.success and enum_name(result.stop) == "reached_target"
        and all(accepted(state) for state in result.states),
        "continuation_stop": enum_name(result.stop),
        "cause_stop": enum_name(result.cause_stop),
        "cause_message": result.cause_message,
        "path_reached": result.reached_path,
        "path_attempted": result.attempted_path,
        "attempts": result.attempts,
        "accepted_steps": result.accepted_steps,
        "retries": result.retries,
        "state_count": len(result.states),
    }
    last = result.last_result
    if last is not None:
        row.update({f"last_{key}": value for key, value in result_fields(last, phase_order).items()})
    else:
        row.update({f"last_{key}": None for key in result_fields_placeholder()})
    terminal = result.states[-1] if result.success and result.states else None
    row.update(comparison_fields(
        anchor,
        terminal,
        phase_order,
        "not_requested" if anchor is None else ("no_target_state" if terminal is None else None),
    ))
    return row


def result_fields_placeholder() -> tuple[str, ...]:
    return (
        "success", "solver_status", "message", "topology_event", "topology_phase", "topology_pair",
        "phase_order", "amounts_mol", "molar_densities_mol_m3", "mole_fractions",
        "temperature_K", "pressure_Pa", "iterations", "evaluations", "hessian_calls",
        "requested_tolerance_met", "phase_count", "convergence_stop", "residual_inf",
        "classified_residual_inf", "convergence_residual", "criterion_name", "criterion_value", "floor",
        "accepted_by_matrix",
    )


def state_rows(
    case_id: str,
    result: Any,
    path_value: float,
    index: int,
    *,
    phase_order: tuple[str, ...] | None = None,
    anchor: Any | None = None,
) -> dict[str, object]:
    fields = result_fields(result, phase_order)
    if anchor is None:
        fields.update(comparison_fields(None, None, phase_order, "not_requested"))
    elif (
        result.temperature
        and result.pressure is not None
        and abs(result.temperature - TARGET_T) <= TARGET_T_ATOL
        and abs(result.pressure / anchor.pressure - 1.0) <= TARGET_P_RTOL
    ):
        fields.update(comparison_fields(anchor, result, phase_order))
    else:
        fields.update(comparison_fields(anchor, None, phase_order, "not_target_endpoint"))
    return {"case_id": case_id, "state_index": index, "path_value": path_value, **fields}


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    expected_hashes = {
        "wheel": WHEEL_SHA256,
        "packet": PACKET_SHA256,
    }
    observed_hashes = {"wheel": sha256(WHEEL), "packet": sha256(PACKET)}
    if observed_hashes != expected_hashes:
        raise SystemExit(f"input hash mismatch: expected {expected_hashes}, observed {observed_hashes}")
    distribution = importlib.metadata.distribution("epcsaft")
    if distribution.version != "0.2.0.dev0":
        raise SystemExit(f"unexpected installed epcsaft version: {distribution.version}")

    model = Mixture(Parameters.from_json(PACKET, components=("methane", "ethane")))
    fixed = solve_equilibrium(model, bubble_boundary(model, TARGET_T, Free(1.0e5)))
    if not accepted(fixed):
        raise SystemExit(f"fixed-temperature anchor did not satisfy acceptance: {result_fields(fixed)}")
    pressure = fixed.pressure
    if abs(pressure / TARGET_P - 1.0) > 1.0e-14:
        raise SystemExit(f"canonical pressure changed: {pressure}")

    rows: list[dict[str, object]] = []
    state_rows_all: list[dict[str, object]] = []
    for guess in COLD_GUESSES:
        result = solve_equilibrium(model, free_temperature(model, guess))
        rows.append(direct_row(
            f"free-T/cold/{guess:g}", "free_temperature_cold", "target branch or typed topology event",
            result, start_temperature_K=guess, target_temperature_K=TARGET_T, target_pressure_Pa=TARGET_P,
            start_mode="cold", units="T:K; P:Pa; feed:mol; residual:1",
            phase_order=("L", "V"), anchor=fixed,
        ))

    for seed_temperature in WARM_SEEDS:
        seed = solve_equilibrium(model, bubble_boundary(model, seed_temperature, Free(1.0e5)))
        result = solve_equilibrium(model, free_temperature(model, seed_temperature), seed)
        rows.append(direct_row(
            f"free-T/warm-from-fixed/{seed_temperature:g}", "free_temperature_warm", "accepted target branch",
            result, start_temperature_K=seed_temperature, target_temperature_K=TARGET_T, target_pressure_Pa=TARGET_P,
            start_mode="accepted fixed-T/P boundary result", units="T:K; P:Pa; feed:mol; residual:1",
            phase_order=("L", "V"), anchor=fixed,
        ))

    continuation_cases = (
        ("continuation/forward/179-to-180", "forward", 179.0, 180.0, CONTINUATION, "reached target"),
        ("continuation/reverse/181-to-180", "reverse", 181.0, 180.0, CONTINUATION, "reached target"),
        ("continuation/forward/179-to-181", "forward-long", 179.0, 181.0, LONG_CONTINUATION, "bounded branch loss"),
        ("continuation/reverse/181-to-179", "reverse-long", 181.0, 179.0, LONG_CONTINUATION, "bounded branch loss"),
    )
    for case_id, direction, start_t, target_t, options, expected in continuation_cases:
        result = continue_equilibrium(
            model,
            bubble_boundary(model, start_t, Free(1.0e5)),
            bubble_boundary(model, target_t, Free(1.0e5)),
            Continuation(**options),
        )
        rows.append(continuation_row(
            case_id, direction, expected, result, start_temperature_K=start_t,
            target_temperature_K=target_t, target_pressure_Pa="free", options=json.dumps(options, sort_keys=True),
            units="T:K; P:Pa; feed:mol; path:1; residual:1",
            phase_order=("L", "V"), anchor=fixed,
        ))
        state_rows_all.extend(
            state_rows(case_id, state, q, index, phase_order=("L", "V"), anchor=fixed)
            for index, (q, state) in enumerate(zip(result.path_values, result.states))
        )

    reactive = reactive_model()
    reactive_start = reactive_problem(298.15, correlation=correlation())
    reactive_result = solve_equilibrium(reactive, reactive_start)
    rows.append(direct_row(
        "reactive/fixed-298.15", "reactive_fixed_TP", "accepted target branch", reactive_result,
        start_temperature_K=298.15, target_temperature_K=298.15, target_pressure_Pa=100000.0,
        start_mode="cold", units="T:K; P:Pa; feed:mol; residual:1", phase_order=("L1", "L2"),
    ))
    reactive_target = reactive_problem(299.15, correlation=correlation())
    reactive_warm = solve_equilibrium(reactive, reactive_target, reactive_result)
    rows.append(direct_row(
        "reactive/warm-298.15-to-299.15", "reactive_warm", "accepted target branch", reactive_warm,
        start_temperature_K=298.15, target_temperature_K=299.15, target_pressure_Pa=100000.0,
        start_mode="accepted reactive result", units="T:K; P:Pa; feed:mol; residual:1", phase_order=("L1", "L2"),
    ))
    reactive_continuation = continue_equilibrium(
        reactive, reactive_start, reactive_target,
        Continuation(initial_step=0.5, minimum_step=0.5, max_steps=2, max_retries_per_step=1),
    )
    rows.append(continuation_row(
        "reactive/continuation-298.15-to-299.15", "forward", "reached target", reactive_continuation,
        start_temperature_K=298.15, target_temperature_K=299.15, target_pressure_Pa=100000.0,
        options=json.dumps({"initial_step": 0.5, "minimum_step": 0.5, "max_steps": 2, "max_retries_per_step": 1}, sort_keys=True),
        units="T:K; P:Pa; feed:mol; path:1; residual:1", phase_order=("L1", "L2"),
    ))
    state_rows_all.extend(
        state_rows("reactive/continuation-298.15-to-299.15", state, q, index, phase_order=("L1", "L2"))
        for index, (q, state) in enumerate(zip(reactive_continuation.path_values, reactive_continuation.states))
    )

    OUT.mkdir(parents=True, exist_ok=True)
    write_csv(OUT / "cases.csv", rows)
    write_csv(OUT / "continuation_states.csv", state_rows_all)
    manifest = {
        "schema_version": 1,
        "analysis": "2026-initialization-diagnosis",
        "issue": 90,
        "source_commit": "f909bb21a93171c38465adba8777cd3b0bc6c53d",
        "branch": "roadmap/application-preparation",
        "wheel": str(WHEEL.relative_to(ROOT)),
        "wheel_sha256": observed_hashes["wheel"],
        "installed_distribution": distribution.version,
        "installed_location": str(distribution.locate_file("")),
        "parameter_packet": str(PACKET.relative_to(ROOT)),
        "parameter_packet_sha256": observed_hashes["packet"],
        "synthetic_reactive_fixture": str(SYNTHETIC_SOURCE.relative_to(ROOT)),
        "synthetic_reactive_fixture_sha256": sha256(SYNTHETIC_SOURCE),
        "frozen_matrix": {
            "cold_free_temperature_guesses_K": list(COLD_GUESSES),
            "warm_fixed_temperature_seeds_K": list(WARM_SEEDS),
            "canonical_free_temperature_target_K": TARGET_T,
            "canonical_free_temperature_pressure_Pa": TARGET_P,
            "row_tolerance": ROW_TOLERANCE,
            "temperature_atol_K": TARGET_T_ATOL,
            "pressure_rtol": TARGET_P_RTOL,
            "continuation": CONTINUATION,
            "long_continuation": LONG_CONTINUATION,
            "reactive_state": {"T_K": 298.15, "P_Pa": 100000.0, "reaction": "A+B=C; ln K=1500(1/T-1/298.15)"},
            "reactive_continuation": {"initial_step": 0.5, "minimum_step": 0.5, "max_steps": 2, "max_retries_per_step": 1},
        },
        "acceptance": {
            "accepted_state": "success=true, requested_tolerance_met=true, topology_event=None_, convergence residual <= 1e-10",
            "typed_failure": "preserve topology event and finite residual/solver diagnostics; do not infer physical nonexistence",
            "continuation": "reached_target with every retained state accepted, or preserve bounded branch_lost/step_too_small cause and accepted prefix",
        },
        "fixed_180_anchor": {
            "phase_order": ["L", "V"],
            "temperature_K": fixed.temperature,
            "pressure_Pa": fixed.pressure,
            "amounts_mol": fixed.amounts,
            "molar_densities_mol_m3": fixed.molar_densities,
            "mole_fractions": fixed.mole_fractions,
            "comparison_criteria": {
                "pressure_relative": BRANCH_PRESSURE_RTOL,
                "composition_absolute": BRANCH_COMPOSITION_ATOL,
                "amount_absolute_mol": BRANCH_AMOUNT_ATOL_MOL,
                "density_log_window": BRANCH_DENSITY_LOG_WINDOW,
            },
        },
        "review_history": {
            "initial": "Sol scientific review: diagnosis supported; design revised for branch identity, retry meaning and fallback contract.",
            "correction_round_1": "Retains a fixed-180 anchor and phase-ordered x/rho/amount comparisons; removes the unsupported generic 12-invocation fallback and records q-step/retry semantics. Fresh execution review remains pending.",
        },
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "case_count": len(rows),
        "continuation_state_count": len(state_rows_all),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    # Result-level assertions are intentionally narrow: accepted rows must meet the frozen
    # convergence contract; diagnostic cold rows may only fail with a typed topology event.
    for row in rows:
        if row["family"] == "free_temperature_cold":
            if not row["accepted_by_matrix"] and row["topology_event"] != "PhasesCoalesced":
                raise SystemExit(f"untyped cold-start failure: {row}")
        elif row["family"] in {"free_temperature_warm", "reactive_fixed_TP", "reactive_warm"}:
            if not row["accepted_by_matrix"]:
                raise SystemExit(f"unexpected warm/fixed failure: {row}")
        elif row["family"] == "continuation":
            if row["expected"] == "reached target" and not row["accepted_by_matrix"]:
                raise SystemExit(f"unexpected continuation failure: {row}")
            if row["expected"] == "bounded branch loss" and row["continuation_stop"] != "step_too_small":
                raise SystemExit(f"unexpected long-path continuation result: {row}")


def declarations():
    return (
        free_temperature(None, 180.0),
        reactive_model(),
        reactive_problem(298.15, correlation=correlation()),
    )


if __name__ == "__main__":
    main()

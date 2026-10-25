#!/usr/bin/env python3
"""Bounded #90 staged-start diagnosis on the installed wheel."""

from __future__ import annotations

import importlib.metadata
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from typing import Any

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(SCRIPT_DIR))

import run as baseline  # noqa: E402
from epcsaft import Mixture, Parameters  # noqa: E402
from epcsaft.equilibrium import Free, solve_equilibrium  # noqa: E402


OUT = ROOT / "analyses" / "2026-initialization-diagnosis" / "results"
BASELINE_MANIFEST = OUT / "manifest.json"
WHEEL = ROOT / "build" / "environment-wheel" / "epcsaft-0.2.0.dev0-cp313-cp313-linux_x86_64.whl"
PACKET = ROOT / "data" / "packets" / "gross-2001-methane-ethane" / "2" / "parameters" / "parameters.json"
PHASE_ORDER = ("L", "V")
ANCHOR_PRESSURE_GUESS = 1.0e5
CALL_LIMIT = 10
TIME_LIMIT_S = 30.0
GUESSES_K = (170.0, 179.0, 180.0, 181.0, 190.0)


def retained_anchor() -> tuple[SimpleNamespace, dict[str, Any]]:
    payload = json.loads(BASELINE_MANIFEST.read_text())
    anchor = payload["fixed_180_anchor"]
    if (
        anchor["temperature_K"] != baseline.TARGET_T
        or anchor["pressure_Pa"] != baseline.TARGET_P
        or anchor["phase_order"] != list(PHASE_ORDER)
    ):
        raise SystemExit("retained fixed-180 anchor does not match the frozen target")
    reference = SimpleNamespace(
        pressure=[anchor["pressure_Pa"]],
        amounts=anchor["amounts_mol"],
        molar_densities=anchor["molar_densities_mol_m3"],
        mole_fractions=anchor["mole_fractions"],
    )
    return reference, anchor


def stage_fields(result: Any | None, prefix: str) -> dict[str, object]:
    if result is None:
        values = {key: None for key in baseline.result_fields_placeholder()}
        values["stage_accepted_by_matrix"] = None
    else:
        values = baseline.result_fields(result, PHASE_ORDER)
        values["stage_accepted_by_matrix"] = baseline.accepted(result)
    values.pop("accepted_by_matrix", None)
    return {f"{prefix}_{key}": value for key, value in values.items()}


def stage_meaning(result: Any | None, attempted: bool, stage: str) -> str:
    if not attempted:
        return "not_attempted: accepted fixed-T/free-P anchor required"
    if baseline.accepted(result):
        return f"accepted {stage} result on declared L,V branch"
    event = baseline.enum_name(result.topology_event)
    if event and event != "None_":
        return f"typed topology event {event}; numerical branch result, not physical nonexistence"
    stop = baseline.enum_name(result.convergence.stop)
    return f"numerical solve rejected ({stop or result.solver_status}); inspect retained diagnostics"


def comparison_fields(reference: Any, result: Any | None) -> dict[str, object]:
    values = baseline.comparison_fields(
        reference,
        result,
        PHASE_ORDER,
        None if result is not None else "not_attempted_anchor_rejected",
    )
    names = {
        "anchor_match": "branch_match",
        "anchor_pressure_rel_error": "anchor_pressure_rel_error",
        "anchor_composition_abs_error": "anchor_composition_abs_error",
        "anchor_density_log_error": "anchor_density_log_error",
        "anchor_amount_abs_error_mol": "anchor_amount_abs_error_mol",
        "anchor_compare_reason": "compare_reason",
    }
    return {f"target_{names[key]}": value for key, value in values.items()}


def main() -> None:
    reference, anchor_metadata = retained_anchor()
    observed_hashes = {"wheel": baseline.sha256(WHEEL), "packet": baseline.sha256(PACKET)}
    expected_hashes = {"wheel": baseline.WHEEL_SHA256, "packet": baseline.PACKET_SHA256}
    if observed_hashes != expected_hashes:
        raise SystemExit(f"input hash mismatch: expected {expected_hashes}, observed {observed_hashes}")
    distribution = importlib.metadata.distribution("epcsaft")
    if distribution.version != "0.2.0.dev0":
        raise SystemExit(f"unexpected installed epcsaft version: {distribution.version}")

    model = Mixture(Parameters.from_json(PACKET, components=("methane", "ethane")))
    rows: list[dict[str, object]] = []
    call_count = 0
    started = time.perf_counter()

    def driver_call(problem: Any, start: Any | None = None) -> Any:
        nonlocal call_count
        if call_count >= CALL_LIMIT:
            raise SystemExit(f"driver-call limit exceeded: {call_count + 1}")
        call_count += 1
        return solve_equilibrium(model, problem, start)

    for guess in GUESSES_K:
        anchor_result = driver_call(baseline.bubble_boundary(model, guess, Free(ANCHOR_PRESSURE_GUESS)))
        anchor_accepted = baseline.accepted(anchor_result)
        target_result = None
        if anchor_accepted:
            target_result = driver_call(
                baseline.free_temperature(model, guess),
                anchor_result,
            )
        target_accepted = target_result is not None and baseline.accepted(target_result)
        target_comparison = comparison_fields(reference, target_result)
        target_branch_match = target_comparison["target_branch_match"] is True
        target_temperature_ok = (
            target_result is not None
            and bool(target_result.temperature)
            and abs(target_result.temperature - baseline.TARGET_T) <= baseline.TARGET_T_ATOL
        )
        target_pressure_ok = (
            target_result is not None
            and bool(target_result.pressure)
            and abs(target_result.pressure / baseline.TARGET_P - 1.0) <= baseline.TARGET_P_RTOL
        )
        row: dict[str, object] = {
            "case_id": f"staged-start/{guess:g}",
            "caller_temperature_guess_K": guess,
            "anchor_pressure_guess_Pa": ANCHOR_PRESSURE_GUESS,
            "target_pressure_Pa": baseline.TARGET_P,
            "phase_order": ",".join(PHASE_ORDER),
            "units": "T:K; P:Pa; feed:mol; residual:1",
            "anchor_call_made": True,
            "target_call_made": target_result is not None,
            "anchor_accepted_by_matrix": anchor_accepted,
            "anchor_failure_meaning": stage_meaning(anchor_result, True, "fixed-T/free-P anchor"),
            "target_failure_meaning": stage_meaning(target_result, target_result is not None, "free-T/fixed-P target"),
            "target_temperature_match": target_temperature_ok if target_result is not None else None,
            "target_pressure_match": target_pressure_ok if target_result is not None else None,
            "target_accepted_by_matrix": (
                target_accepted and target_temperature_ok and target_pressure_ok and target_branch_match
            ) if target_result is not None else None,
        }
        row.update(stage_fields(anchor_result, "anchor"))
        row.update(stage_fields(target_result, "target"))
        row.update(target_comparison)
        rows.append(row)

    elapsed_s = time.perf_counter() - started
    if call_count > CALL_LIMIT:
        raise SystemExit(f"driver-call limit exceeded: {call_count}")
    if elapsed_s > TIME_LIMIT_S:
        raise SystemExit(f"staged-start time limit exceeded: {elapsed_s:.3f} s")

    OUT.mkdir(parents=True, exist_ok=True)
    baseline.write_csv(OUT / "staged-start.csv", rows)
    manifest = {
        "schema_version": 1,
        "analysis": "2026-initialization-diagnosis/staged-start",
        "issue": 90,
        "source_commit": "f909bb21a93171c38465adba8777cd3b0bc6c53d",
        "branch": "roadmap/application-preparation",
        "wheel": str(WHEEL.relative_to(ROOT)),
        "wheel_sha256": observed_hashes["wheel"],
        "parameter_packet": str(PACKET.relative_to(ROOT)),
        "parameter_packet_sha256": observed_hashes["packet"],
        "installed_distribution": distribution.version,
        "workflow": "one fixed-T/free-P anchor at caller guess, then one original free-T/fixed-P solve from that accepted result",
        "frozen_inputs": {
            "temperature_guesses_K": list(GUESSES_K),
            "anchor_pressure_guess_Pa": ANCHOR_PRESSURE_GUESS,
            "target_pressure_Pa": baseline.TARGET_P,
            "feed_mol": {"methane": 0.5, "ethane": 0.5},
            "phase_order": list(PHASE_ORDER),
            "row_tolerance": baseline.ROW_TOLERANCE,
            "temperature_atol_K": baseline.TARGET_T_ATOL,
            "pressure_rtol": baseline.TARGET_P_RTOL,
            "driver_call_limit": CALL_LIMIT,
            "time_limit_s": TIME_LIMIT_S,
        },
        "reference_anchor": {
            "source": str(BASELINE_MANIFEST.relative_to(ROOT)),
            "temperature_K": anchor_metadata["temperature_K"],
            "pressure_Pa": anchor_metadata["pressure_Pa"],
            "phase_order": anchor_metadata["phase_order"],
            "comparison_criteria": anchor_metadata["comparison_criteria"],
        },
        "acceptance": {
            "stage": "success=true, requested_tolerance_met=true, topology_event=None_, convergence residual <= 1e-10",
            "target": "stage acceptance plus T=180 K, P=1596859.1113026983 Pa, and retained E4/P2 anchor comparison",
            "failure": "retain topology/solver diagnostics; no physical nonexistence inference",
        },
        "driver_calls": call_count,
        "elapsed_s": elapsed_s,
        "outcomes": [
            {
                "guess_K": row["caller_temperature_guess_K"],
                "anchor_accepted": row["anchor_accepted_by_matrix"],
                "target_call_made": row["target_call_made"],
                "target_accepted": row["target_accepted_by_matrix"],
                "target_branch_match": row["target_branch_match"],
                "target_topology_event": row["target_topology_event"],
            }
            for row in rows
        ],
        "row_count": len(rows),
        "run_utc": datetime.now(timezone.utc).isoformat(),
    }
    (OUT / "staged-start-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()

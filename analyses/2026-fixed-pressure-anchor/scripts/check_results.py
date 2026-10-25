"""Check retained fixed-pressure rows and summary without running the model."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
CHARGES = (0, 0, 0, 1, -1, -1, -2, 1, -1)
REQUIRED_FIELDS = {
    "pressure_iterations", "anchor_pressure_iterations", "certified_branch",
    "derivative_condition_measure", "local_search_outcome", "reachability_envelope",
}
EXPECTED_ERROR_CODES = {
    "wrong-phase": "density-anchor-phase-mismatch",
    "model-mismatch": "density-anchor-model-mismatch",
    "component-order-mismatch": "density-anchor-model-mismatch",
    "temperature-mismatch": "density-anchor-state-mismatch",
    "pressure-mismatch": "density-anchor-state-mismatch",
    "wrong-type": "invalid-density-anchor",
    "nonfinite-density-double": "invalid-density-anchor",
    "nonpositive-density-double": "invalid-density-anchor",
    "local-envelope": "density-anchor-branch-unreachable",
    "non-evaluable-derivative": "provider_domain_rejection",
}
EXPECTED_MAXIMUM_RATIOS = {
    "charged": 0.31332633890949413,
    "neutral": 0.3489622766556471,
}
MISSING = object()
EXPECTED_ERROR_CONTEXTS = {
    **dict.fromkeys(
        (
            "wrong-phase", "model-mismatch", "component-order-mismatch",
            "temperature-mismatch", "pressure-mismatch",
        ),
        [],
    ),
    "wrong-type": [["reason", "expected-State"]],
    "nonfinite-density-double": [["reason", "missing-or-invalid-identity"]],
    "nonpositive-density-double": [["reason", "missing-or-invalid-identity"]],
    "local-envelope": [
        ["reason", "continued density reference is outside the model domain"],
        ["max_abs_log_density", "0.009950330853168083"],
        ["centered_scan_intervals", "5"],
        ["domain_scan_intervals", "360"],
    ],
    "non-evaluable-derivative": MISSING,
}
EXPECTED_MESSAGES = {
    "wrong-phase": "density-anchor-phase-mismatch",
    "model-mismatch": "density-anchor-model-mismatch",
    "component-order-mismatch": "density-anchor-model-mismatch",
    "temperature-mismatch": "density-anchor-state-mismatch",
    "pressure-mismatch": "density-anchor-state-mismatch",
    "wrong-type": "invalid-density-anchor (reason=expected-State)",
    "nonfinite-density-double": "invalid-density-anchor (reason=missing-or-invalid-identity)",
    "nonpositive-density-double": "invalid-density-anchor (reason=missing-or-invalid-identity)",
    "local-envelope": "density-anchor-branch-unreachable (reason=continued density reference is outside the model domain, max_abs_log_density=0.009950330853168083, centered_scan_intervals=5, domain_scan_intervals=360)",
    "non-evaluable-derivative": MISSING,
}
EXPECTED_NON_EVALUABLE_STATE_IDENTITY = (
    "sha256:666373fe68285e7d004a22d23ba32adb0304249101bf73df03e09e7c847d35a4"
)
EXPECTED_SCOPE = (
    "Local fixed-pressure branch-continuation evidence only; not global "
    "phase-equilibrium or reactive-film validation."
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def decoded(value: str) -> object:
    if value == "":
        return None
    if value[:1] in "[{\"" or value in {"null", "true", "false"}:
        return json.loads(value)
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def flatten(value: object) -> tuple[float, ...]:
    if isinstance(value, (list, tuple)):
        return tuple(item for child in value for item in flatten(child))
    return (float(value),)


def close(left: object, right: object) -> bool:
    return all(
        abs(a - b) <= max(1e-12, 1e-10 * max(abs(a), abs(b)))
        for a, b in zip(flatten(left), flatten(right), strict=True)
    )


def state_identity(row: dict[str, object]) -> str:
    payload = json.dumps(
        {
            "schema": "epcsaft-state-identity-v1",
            "artifact_fingerprint": row["artifact_fingerprint"],
            "parameter_fingerprint": row["parameter_fingerprint"],
            "topology_fingerprint": row["topology_fingerprint"],
            "component_order_fingerprint": row["component_order_fingerprint"],
            "component_ids": row["component_ids"],
            "temperature_k": float(row["temperature_k"]).hex(),
            "pressure_pa": float(row["requested_pressure_pa"]).hex(),
            "mole_fractions": [float(value).hex() for value in row["mole_fractions"]],
            "molar_density_mol_per_m3": float(row["molar_density_mol_per_m3"]).hex(),
            "requested_phase": row["requested_phase"],
            "certified_branch": row["certified_branch"],
        }, allow_nan=False, separators=(",", ":"), sort_keys=True,
    ).encode()
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def compare(cold: dict[str, object], warm: dict[str, object]) -> float:
    exact = (
        "component_ids", "parameter_fingerprint", "topology_fingerprint",
        "component_order_fingerprint", "artifact_fingerprint", "requested_phase",
        "value_status", "derivative_status", "constraint_identities",
        "coordinate_component_ids", "dependent_component_ids",
    )
    require(all(cold[name] == warm[name] for name in exact), "cold/anchor identity parity")
    require(cold["selected_stable_root_index"] == warm["certified_branch"]["stable_root_index"], "stable-root parity")
    for name in (
        "molar_density_mol_per_m3", "residual_chemical_potential_over_rt",
        "log_fugacity_coefficients", "fugacity_coefficients", "fugacity_pa",
        "derivative_condition_measure", "log_composition_basis", "density_derivatives_mol_per_m3",
    ):
        require(close(cold[name], warm[name]), f"cold/anchor value parity: {name}")
    disagreement = max(
        (abs(a - b) / max(abs(a), abs(b), 1e-8)
         for a, b in zip(flatten(cold["chemical_potential_derivatives_over_rt"]),
                         flatten(warm["chemical_potential_derivatives_over_rt"]), strict=True)),
        default=0.0,
    )
    require(disagreement <= 3e-6, "chemical-potential derivative parity")
    return disagreement


def main() -> int:
    summary = json.loads((RESULTS / "summary.json").read_text(encoding="utf-8"))
    schema = json.loads((RESULTS / "calls.schema.json").read_text(encoding="utf-8"))
    require(summary["schema"] == "epcsaft.fixed-pressure-anchor-summary.v1", "summary schema")
    require(schema["columns"], "CSV columns")
    with (RESULTS / "calls.csv").open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames == schema["columns"], "CSV column order")
        rows = [{key: decoded(value) for key, value in row.items()} for row in reader]
    require(len(rows) == 1080, "CSV row count")
    require(hashlib.sha256((RESULTS / "calls.csv").read_bytes()).hexdigest()
            == "98e23ae47ca61b5b0d413f974c51eee918ed37867a8809eb8390332431a80565", "CSV hash")

    by_key = {(row["campaign"], row["direction"], row["repetition"], row["method"], row["target_index"]): row for row in rows}
    require(len(by_key) == len(rows), "CSV row identities")
    maximum = 0.0
    for campaign, targets in (("charged", 33), ("neutral", 21)):
        selected = [row for row in rows if row["campaign"] == campaign]
        require(len(selected) == targets * 2 * 5 * 2, f"{campaign} coverage")
        for row in selected:
            require(row["status"] == "success" and row["value_status"] == "accepted"
                    and row["derivative_status"] == "available", f"{campaign} status")
            require(row["pressure_relative_residual"] <= 1e-11, f"{campaign} pressure closure")
            require(row["state_identity"] == state_identity(row), f"{campaign} state identity")
            require(abs(math.fsum(row["mole_fractions"]) - 1.0) <= 1e-12, f"{campaign} normalization")
            if campaign == "charged":
                require(abs(math.fsum(x * z for x, z in zip(row["mole_fractions"], CHARGES, strict=True))) <= 1e-12, "electroneutrality")
        for direction in ("forward", "reverse"):
            for repetition in range(5):
                seed = by_key[(campaign, direction, repetition, "anchored", 0)]
                require(seed["anchor_status"] == "not_supplied", "anchor seed")
                for target in range(targets):
                    ordinary = by_key[(campaign, direction, repetition, "ordinary", target)]
                    anchored = by_key[(campaign, direction, repetition, "anchored", target)]
                    maximum = max(maximum, compare(ordinary, anchored))
                    require(anchored["certified_branch"] == seed["certified_branch"], "branch chain")
                    if target:
                        previous = by_key[(campaign, direction, repetition, "anchored", target - 1)]
                        require(anchored["anchor_status"] == "accepted"
                                and anchored["anchor_source_state_identity"] == previous["state_identity"]
                                and anchored["anchor_source_composition"] == previous["mole_fractions"]
                                and anchored["anchor_density_mol_per_m3"] == previous["molar_density_mol_per_m3"], "anchor chain")

    comparison = summary["comparison_summary"]
    require(comparison["count"] == 540 and comparison["all_passed"]
            and comparison["max_normalized_derivative_disagreement"] == maximum
            and comparison["normalized_derivative_tolerance"] == 3e-6, "comparison summary")
    for campaign in ("charged", "neutral"):
        performance = summary["performance"][campaign]
        for direction in ("forward", "reverse"):
            measured = performance["sequence_times_ns"][direction]
            require(all(value > 0 for values in measured.values() for value in values), f"timing rows: {campaign}/{direction}")
            ratio = median(measured["anchored"]) / median(measured["ordinary"])
            require(performance["ratios"][direction] == ratio, f"timing ratio: {campaign}/{direction}")
        require(
            performance["maximum_ratio"] == EXPECTED_MAXIMUM_RATIOS[campaign]
            == max(performance["ratios"].values()) <= 0.5,
            f"timing gate: {campaign}",
        )

    attempts = {item["name"]: item for item in summary["negative_cases"]}
    require(
        set(attempts) == set(EXPECTED_ERROR_CODES)
        and all(REQUIRED_FIELDS <= set(item) for item in attempts.values()),
        "negative diagnostics",
    )
    for item in attempts.values():
        require(item["status"].startswith("expected-"), f"negative status: {item['name']}")
        require(
            item["error_code"] == EXPECTED_ERROR_CODES[item["name"]],
            f"negative error code: {item['name']}",
        )
        require(
            item.get("error_context", MISSING) == EXPECTED_ERROR_CONTEXTS[item["name"]],
            f"negative error context: {item['name']}",
        )
        require(
            item.get("message", MISSING) == EXPECTED_MESSAGES[item["name"]],
            f"negative message: {item['name']}",
        )
    require(
        attempts["non-evaluable-derivative"]["state_identity"]
        == EXPECTED_NON_EVALUABLE_STATE_IDENTITY,
        "non-evaluable state identity",
    )
    require(summary["claim_limits"] == [EXPECTED_SCOPE], "claim scope")
    require(attempts["non-evaluable-derivative"]["value_status"] == "accepted"
            and attempts["non-evaluable-derivative"]["derivative_status"] == "non_evaluable", "value/derivative separation")
    local = attempts["local-envelope"]
    require(local["local_search_outcome"] == "rejected"
            and local["reachability_envelope"] == {"max_abs_log_density": math.log1p(1e-2), "centered_scan_intervals": 5, "domain_scan_intervals": 360}, "local rejection diagnostics")
    print(json.dumps({"status": "accepted", "calls": len(rows), "max_normalized_derivative_disagreement": maximum}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

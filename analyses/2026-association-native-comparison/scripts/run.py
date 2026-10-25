#!/usr/bin/env python3
"""Retain one recorded fixed-depth association comparison campaign."""

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import resource
import signal
import subprocess
import sys
import time
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "analyses/2026-association-native-comparison/results"
RUNS = ROOT / "build/numerical-evidence/runs"
TEST_ID = "native.association.comparison"
SOURCE_FILES = (
    "engine/meson.build", "engine/tests/native/association_comparison_checks.cpp",
    "engine/tests/native/properties_support.hpp", "engine/tests/native/recording_support.hpp",
    "analyses/2026-association-native-comparison/scripts/run.py",
    "analyses/2026-association-native-comparison/notebook.qmd",
)
PARTITION_DENOMINATORS = {"controls": 14, "development": 64, "heldout": 32}
STRATUM_DENOMINATORS = {
    "controls": 14,
    "ordinary-dilute": 12,
    "ordinary-dense": 12,
    "strong-dense": 12,
    "heterogeneous-cross-capacity": 12,
    "disconnected-nonassociating": 12,
    "trace-support": 12,
    "exact-zero-support": 12,
    "parameter-boundaries": 12,
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_plot(path: Path, rows: list[dict[str, str]]) -> dict[str, float]:
    ratios, group_maxima = [], {}
    for row in rows:
        if row.get("verdict") not in {"True", "False"}:
            continue
        if criterion_result(row) is None:
            raise ValueError(f"malformed criterion row: {row.get('comparison_id', '<missing>')}")
        try:
            model, reference = float(row["model"]), float(row["reference"])
            atol, rtol = float(row["atol"]), float(row["rtol"])
        except (TypeError, ValueError):
            continue
        if not all(math.isfinite(value) for value in (model, reference, atol, rtol)):
            continue
        error = abs(model - reference)
        rule = row.get("rule", "")
        scale = max(1.0, abs(reference)) if "or_one" in rule else abs(reference)
        limit = atol + rtol * scale
        ratio = error / limit if limit > 0.0 else math.inf if error > 0.0 else 0.0
        ratios.append(ratio)
        group = row.get("comparison_group_id", "unknown")
        group_maxima[group] = max(group_maxima.get(group, 0.0), ratio)
    figure, axis = plt.subplots(figsize=(7, 3))
    axis.semilogy(ratios, label="error / declared threshold", linewidth=0.8)
    axis.axhline(1.0, color="black", linestyle="--", linewidth=0.8, label="threshold = 1")
    axis.set(xlabel="criterion row", ylabel="dimensionless error / threshold", title="Association comparison normalized error")
    axis.legend()
    figure.tight_layout()
    figure.savefig(path, dpi=140)
    plt.close(figure)
    return group_maxima


COMPARISON_COLUMNS = ("atol", "basis", "claim", "comparison_group_id", "comparison_id", "coordinates",
                      "criterion", "criterion_name", "evidence_state", "group_ref", "inputs", "inputs_latex",
                      "judgment", "method", "model", "point_id", "quantity", "quantity_latex", "reference",
                      "relationship_latex", "rtol", "rule", "schema_version", "unavailable_reason", "unit", "verdict")


def native_json(log: Path, prefix: str):
    for line in log.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(prefix):
            try:
                yield json.loads(line.split(" ", 1)[1])
            except json.JSONDecodeError:
                yield None


def native_comparisons(log: Path) -> list[dict[str, str]]:
    rows = []
    for payload in native_json(log, "EPCSAFT_COMPARE_JSON "):
        if not isinstance(payload, dict) or not isinstance(payload.get("comparison_id"), str):
            payload = {"comparison_id": f"native.json.invalid.{len(rows) + 1}", "evidence_state": "Unavailable"}
        row = {key: payload.get(key, "") for key in COMPARISON_COLUMNS}
        for key in ("inputs", "coordinates"):
            row[key] = json.dumps(payload.get(key, {}), sort_keys=True, separators=(",", ":"))
        if payload.get("criterion"):
            row["criterion_name"] = payload["criterion"]
            row["criterion"] = json.dumps({key: payload.get(key, "") for key in ("atol", "rtol", "rule", "verdict")}, sort_keys=True, separators=(",", ":"))
        if isinstance(payload.get("verdict"), bool):
            row["verdict"] = str(payload["verdict"])
        row["schema_version"] = str(payload.get("schema_version", 1))
        rows.append(row)
    return rows


def input_count(value) -> int:
    try:
        parsed = value if isinstance(value, list) else value if isinstance(value, (int, float)) else json.loads(value)
        return len(parsed) if isinstance(parsed, list) else 1
    except (TypeError, json.JSONDecodeError):
        return -1


CRITERION_RULES = {"symmetric_close", "reference_scale", "strict_absolute", "absolute",
                   "strict_symmetric_close", "strict_reference_scale", "strict_reference_or_one",
                   "reference_or_one"}


def criterion_result(row: dict[str, str]) -> bool | None:
    try:
        criterion = json.loads(row.get("criterion", ""))
        model, reference = (float(row[key]) for key in ("model", "reference"))
        atol, rtol = (float(row[key]) for key in ("atol", "rtol"))
        rule = row["rule"]
        nested_verdict = criterion["verdict"]
        if (not isinstance(criterion, dict) or rule not in CRITERION_RULES
                or row.get("criterion_name") != rule or criterion.get("rule") != rule
                or not isinstance(nested_verdict, bool)
                or not all(math.isfinite(value) for value in (model, reference, atol, rtol))
                or atol < 0.0 or rtol < 0.0):
            return None
        error = abs(model - reference)
        if rule in {"strict_absolute", "absolute"}:
            limit = atol
        elif "symmetric" in rule:
            limit = max(atol, rtol * max(abs(model), abs(reference)))
        else:
            scale = max(1.0, abs(reference)) if "or_one" in rule else abs(reference)
            limit = atol + rtol * scale
        result = error < limit if rule.startswith("strict_") else error <= limit
        return result if (row.get("verdict") == str(result) and nested_verdict == result
                          and float(criterion["atol"]) == atol and float(criterion["rtol"]) == rtol) else None
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        return None


def required_tensor_ids(inputs: list[dict], rows: list[dict[str, str]]) -> set[str]:
    dimensions = {str(item.get("case_id")): (input_count(item.get("segment_count")), item.get("active_family") != "none") for item in inputs}
    states = {}
    for row in rows:
        match = re.match(r"d(\d+)/(?:association|residual)-(?:gradient|hessian)/([^/]+)/", row.get("comparison_id", ""))
        if match:
            states[(match.group(1), match.group(2))] = states.get((match.group(1), match.group(2)), False) or row.get("evidence_state") == "Unavailable"
    expected = set()
    for (depth, case), unavailable in states.items():
        components, active = dimensions.get(case, (-1, False))
        if components < 0:
            continue
        width = 2 + components + active
        for group in ("association-gradient", "residual-gradient"):
            expected.update(f"d{depth}/{group}/{case}/{index}" for index in range(width))
        for group in ("association-hessian", "residual-hessian"):
            expected.update(f"d{depth}/{group}/{case}/{entry if unavailable else f'{entry // width}-{entry % width}'}"
                            for entry in range(width * width))
    return expected


def required_value_ids(inputs: list[dict], rows: list[dict[str, str]]) -> set[str]:
    sites = {str(item.get("case_id")): input_count(item.get("site_component")) for item in inputs}
    states = {(match.group(1), match.group(2)) for row in rows
              for match in [re.match(r"d(\d+)/[^/]+/([^/]+)(?:/|$)", row.get("comparison_id", ""))] if match}
    expected = set()
    for depth, case in states:
        for group, count in (("site-fraction", sites.get(case, -1)), ("candidate-residual", sites.get(case, -1)),
                             ("reference-residual", sites.get(case, -1)), ("association-energy", 1),
                             ("residual-energy", 1)):
            if count >= 0:
                expected.update(f"d{depth}/{group}/{case}" + (f"/{index}" if group.endswith("fraction") or "residual" in group and group != "residual-energy" else "")
                                for index in range(count))
    return expected


def native_validation(log: Path, rows: list[dict[str, str]], inputs: list[dict]) -> dict:
    expectations = list(native_json(log, "EPCSAFT_EXPECT_JSON "))
    errors = ["malformed native expectation" for item in expectations
              if not isinstance(item, dict) or item.get("malformed")]
    identifiers = [row.get("comparison_id", "") for row in rows]
    duplicate_counts = Counter(identifiers)
    duplicates = sorted(identifier for identifier, count in duplicate_counts.items() if count > 1)
    if duplicates:
        errors.append(f"duplicate comparison IDs: {duplicates[:5]}")
    declared = []
    for item in expectations:
        if not isinstance(item, dict):
            continue
        keys = item.get("expected_comparison_keys")
        prefix, count = item.get("comparison_id_prefix"), item.get("expected_count")
        valid_keys = isinstance(keys, list) and bool(keys) and all(isinstance(value, str) and value for value in keys)
        valid_count = (isinstance(prefix, str) and bool(prefix) and isinstance(count, int)
                       and not isinstance(count, bool) and count >= 0)
        if not (valid_keys or valid_count) or keys is not None and not valid_keys or (prefix is not None or count is not None) and not valid_count:
            errors.append("malformed native expectation declarations")
        if valid_keys:
            declared.extend(keys)
    counts = [(item.get("comparison_id_prefix"), item.get("expected_count")) for item in expectations
              if isinstance(item, dict) and "expected_count" in item]
    if not declared and not counts:
        errors.append("native output omitted expected comparison declarations")
    if any(not isinstance(prefix, str) or not isinstance(count, int) or count < 0 for prefix, count in counts):
        errors.append("malformed native expectation count")
    actual = set(identifiers)
    matches = lambda pattern, identifier: identifier == pattern or pattern.endswith("*") and (identifier.startswith(pattern[:-1]) or f"/{pattern[:-1]}" in identifier)
    if declared or counts:
        errors.extend(f"comparison ID is not declared: {identifier or '<missing>'}" for identifier in identifiers
                      if not identifier or not any(matches(pattern, identifier) for pattern in declared)
                      and not any(identifier.startswith(prefix) for prefix, _ in counts if isinstance(prefix, str)))
    errors.extend(f"expected {count} IDs with prefix {prefix!r}, found {sum(identifier.startswith(prefix) for identifier in actual)}"
                  for prefix, count in counts if isinstance(prefix, str) and isinstance(count, int)
                  and sum(identifier.startswith(prefix) for identifier in actual) != count)
    expected = required_tensor_ids(inputs, rows) | required_value_ids(inputs, rows)
    tensor = {identifier for identifier in actual if re.match(r"^d\d+/(?:association|residual)-(?:gradient|hessian)/", identifier)}
    missing, extra = expected - actual, tensor - expected
    if missing:
        errors.append(f"missing required comparison IDs: {sorted(missing)[:5]}")
    if extra:
        errors.append(f"unexpected tensor comparison IDs: {sorted(extra)[:5]}")
    unavailable_states = {comparison_state(row) for row in rows if row.get("evidence_state") == "Unavailable"
                          and row.get("comparison_id") in tensor}
    for row in rows:
        required = row.get("comparison_id") in expected
        unavailable = row.get("evidence_state") == "Unavailable"
        criterion_bearing = any(row.get(key) not in (None, "")
                                for key in ("criterion_name", "criterion", "rule", "verdict"))
        if required and (unavailable or comparison_state(row) in unavailable_states) and not criterion_bearing:
            continue
        if row.get("verdict") not in {"True", "False"}:
            if required:
                errors.append(f"required criterion is missing: {row.get('comparison_id', '<missing>')}")
            continue
        if not all(row.get(key) for key in ("comparison_id", "criterion_name", "criterion", "rule")):
            errors.append(f"criterion metadata is incomplete: {row.get('comparison_id', '<missing>')}")
        elif any(row.get(key) in (None, "") for key in ("model", "reference", "atol", "rtol")):
            errors.append(f"criterion numeric metadata is incomplete: {row['comparison_id']}")
        elif row.get("evidence_state") != "Complete":
            errors.append(f"criterion row is not Complete: {row['comparison_id']}")
        elif criterion_result(row) is None:
            errors.append(f"criterion metadata or verdict is invalid: {row['comparison_id']}")
    return {
        "complete": not errors,
        "expected_declarations": len(declared),
        "expected_prefix_counts": sum(isinstance(item, dict) and "comparison_id_prefix" in item for item in expectations),
        "comparison_ids": len(actual),
        "errors": errors,
    }


def matching_pilot(artifact: dict) -> dict:
    candidates = []
    for path in RUNS.glob("*/run.json"):
        try:
            record = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if record.get("selector") == "pilot" and record.get("execution") == "Passed" \
                and record.get("artifact_identity", {}).get("identity") == artifact.get("identity"):
            candidates.append(record)
    if not candidates:
        return {"status": "unavailable", "reason": "no matching retained pilot run"}
    record = max(candidates, key=lambda value: value.get("run_id", ""))
    return {"run_id": record["run_id"], "selector": "pilot", "execution": record["execution"],
            "resource_limits": record["resource_limits"], "observed": record["resource_limits"]["observed"]}


def artifact_identity(binary: Path) -> dict:
    wheels = list((ROOT / "build/environment-wheel").glob("epcsaft-*.whl"))
    if len(wheels) != 1:
        raise RuntimeError(f"expected one Engine wheel, found {len(wheels)}")
    preparation_path = ROOT / "build/environment-wheel/preparation.json"
    preparation = json.loads(preparation_path.read_text(encoding="utf-8"))
    wheel_path = str(wheels[0].relative_to(ROOT))
    wheel_sha = sha256(wheels[0])
    if preparation.get("wheel") != wheel_path or preparation.get("wheel_sha256") != wheel_sha:
        raise RuntimeError("canonical preparation does not identify the selected Engine wheel")
    artifacts = [{"kind": "wheel", "path": wheel_path, "sha256": wheel_sha},
                 {"kind": "binary", "path": str(binary.relative_to(ROOT)), "sha256": sha256(binary)}]
    identity = hashlib.sha256(json.dumps(artifacts, sort_keys=True).encode()).hexdigest()
    return {"artifacts": artifacts, "identity": identity, "preparation": {
        "path": str(preparation_path.relative_to(ROOT)), "sha256": sha256(preparation_path)}}


def run_campaign(selector: str) -> tuple[str, dict[str, float]]:
    marker = ROOT / "build/environment-wheel/meson-build-dir"
    if not marker.is_file():
        raise RuntimeError("installed native build marker is missing")
    binaries = [path for path in Path(marker.read_text(encoding="utf-8").strip()).rglob("association_comparison_checks") if path.is_file()]
    if len(binaries) != 1:
        raise RuntimeError(f"expected one association comparison binary, found {len(binaries)}")
    binary = binaries[0]
    command = [str(binary), "association-comparison"]
    run_id = f"{datetime.now(UTC).strftime('%Y%m%dT%H%M%S%fZ')}-{TEST_ID}-single"
    run_dir = RUNS / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    limit = 4 * 1024 * 1024 * 1024
    stdout_path = run_dir / "stdout.log"
    stderr_path = run_dir / "stderr.log"
    termination = "exit"
    try:
        with stdout_path.open("wb") as stdout_file, stderr_path.open("wb") as stderr_file:
            process = subprocess.Popen(
                command, cwd=ROOT,
                stdout=stdout_file, stderr=stderr_file, start_new_session=True,
                preexec_fn=lambda: resource.setrlimit(resource.RLIMIT_AS, (limit, limit)),
                env=dict(os.environ, EPCSAFT_ASSOCIATION_SELECTOR=selector),
            )
            try:
                process.wait(timeout=600)
            except (subprocess.TimeoutExpired, KeyboardInterrupt) as error:
                termination = "timeout" if isinstance(error, subprocess.TimeoutExpired) else "interrupted"
                try:
                    os.killpg(process.pid, signal.SIGKILL if termination == "timeout" else signal.SIGINT)
                except ProcessLookupError:
                    pass
                process.wait()
            exit_code = process.returncode
    except OSError as error:
        termination = "spawn_error"
        stdout_path.write_text("", encoding="utf-8")
        stderr_path.write_text(str(error) + "\n", encoding="utf-8")
        exit_code = 127
    observed = {
        "wall_seconds": time.monotonic() - started,
        "peak_rss_bytes": max(usage_before, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss) * 1024,
    }
    log = stdout_path
    rows = native_comparisons(log)
    inputs = native_inputs(log)
    with (run_dir / "comparisons.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=COMPARISON_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    files = {path: sha256(ROOT / path) for path in SOURCE_FILES}
    snapshot = {"base_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                "input_generation": hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest(), "files": files}
    input_generation = snapshot["input_generation"]
    outcome_reason = {"timeout": "timeout after 600 seconds", "interrupted": "interrupted",
                      "spawn_error": "native process could not be started"}.get(
        termination, None if exit_code == 0 else f"process exit={exit_code}")
    record = {
        "schema_version": 1, "run_id": run_id, "test_id": TEST_ID, "selector": selector,
        "command": command, "output": str(stdout_path.relative_to(ROOT)),
        "stderr": str(stderr_path.relative_to(ROOT)),
        "comparisons_artifact": "comparisons.csv", "execution": "Passed" if exit_code == 0 else "Failed",
        "exit_code": exit_code, "termination": termination, "outcome_reason": outcome_reason,
        "input_generation": input_generation, "source_snapshot": snapshot,
        "artifact_identity": artifact_identity(binary), "resource_limits": {
            "timeout_seconds": 600, "address_space_bytes": limit, "observed": observed,
        },
    }
    record["validation"] = native_validation(log, rows, inputs)
    (run_dir / "run.json").write_text(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8")
    return run_id, observed


def run_record(run_id: str) -> tuple[dict, Path]:
    path = ROOT / "build/numerical-evidence/runs" / run_id / "run.json"
    if not path.is_file():
        raise RuntimeError(f"native run record missing: {run_id}")
    return json.loads(path.read_text(encoding="utf-8")), path


INPUT_COLUMNS = (
    "case_id", "partition", "stratum", "regime", "topology_provenance", "k_ij",
    "active_family", "active_tied_indices", "generation_seed", "T_K", "V_m3",
    "amounts_mol", "segment_count", "sigma_constant_A", "dispersion_energy_over_k_K",
    "site_component", "site_multiplicity", "association_active",
    "association_energy_over_k_K", "association_kappa", "sigma_A", "diameter_A",
    "packing_fraction",
)


def native_inputs(log: Path) -> list[dict]:
    result = []
    for line in log.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("EPCSAFT_ASSOCIATION_INPUT "):
            continue
        try:
            value = json.loads(line.split(" ", 1)[1])
        except json.JSONDecodeError:
            continue
        if isinstance(value, dict) and value.get("case_id"):
            result.append(value)
    return result


def native_states(log: Path) -> dict[tuple[str, int], dict[str, str]]:
    result = {}
    pattern = re.compile(r"(?:^| )([A-Za-z_][A-Za-z0-9_]*)=([^ ]*)")
    for line in log.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("EPCSAFT_ASSOCIATION_STATE "):
            continue
        fields = dict(pattern.findall(line.split(" ", 1)[1]))
        if fields.get("case_id") and fields.get("depth"):
            result[(fields["case_id"], int(fields["depth"]))] = fields
    return result


def comparison_state(row: dict[str, str]) -> tuple[str, int] | None:
    match = re.match(r"d(\d+)/[^/]+/([^/]+)", row.get("comparison_id", ""))
    return (match.group(2), int(match.group(1))) if match else None


def build_interpretation(manifest: dict, states: dict[tuple[str, int], dict[str, str]],
    rows: list[dict[str, str]], inputs: list[dict], native_complete: bool) -> dict:
    selected = [case for case in manifest["cases"] if case["selected"]]
    depths = manifest["expected_depths"]
    value_groups = {"site-fraction", "candidate-residual", "reference-residual",
                    "association-energy", "residual-energy"}
    derivative_groups = {"association-gradient", "residual-gradient", "association-hessian", "residual-hessian"}
    expected = {}
    for case in selected:
        for depth in depths:
            state = (case["case_id"], depth)
            expected[state] = {group: ({f"d{depth}/{group}/{case['case_id']}/{site}" for site in range(int(case.get("site_count", 0)))}
                                       if group in {"site-fraction", "candidate-residual", "reference-residual"}
                                       else {f"d{depth}/{group}/{case['case_id']}"}) for group in value_groups}
    for identifier in required_tensor_ids(inputs, rows):
        state = comparison_state({"comparison_id": identifier})
        expected.setdefault(state, {}).setdefault(identifier.split("/", 2)[1], set()).add(identifier)
    present, passed = {}, {}
    for row in rows:
        state = comparison_state(row)
        group = row.get("comparison_group_id")
        if state is None or group not in expected.get(state, {}):
            continue
        present.setdefault(state, {}).setdefault(group, set()).add(row.get("comparison_id", ""))
        if row.get("evidence_state") == "Complete" and criterion_result(row) is True:
            passed.setdefault(state, {}).setdefault(group, set()).add(row.get("comparison_id", ""))
    def qualified(groups: set[str]) -> set[tuple[str, int]]:
        return {state for state, all_groups in expected.items()
                if {group for group in all_groups if group in groups} == groups
                and all(present.get(state, {}).get(group, set()) == identifiers
                       and passed.get(state, {}).get(group, set()) == identifiers
                       for group, identifiers in all_groups.items() if group in groups)}
    value_states = qualified(value_groups)
    passed_states = qualified(derivative_groups) & value_states
    if not native_complete:
        value_states = passed_states = set()
    depth_outcomes = {
        str(depth): {
            "passed": sum((case["case_id"], depth) in passed_states for case in selected),
            "criterion_or_reference_unavailable": sum(case["depth_status"].get(str(depth)) == "criterion_or_reference_unavailable" and (case["case_id"], depth) not in passed_states for case in selected),
            "unattempted": sum(case["depth_status"].get(str(depth)) == "unattempted" for case in selected),
        } for depth in depths
    }
    minimum_by_case = {
        case["case_id"]: min((depth for depth in depths if (case["case_id"], depth) in passed_states), default=None)
        for case in selected
    }
    stratum_outcomes = {}
    for stratum, denominator in STRATUM_DENOMINATORS.items():
        cases = [case for case in selected if case["stratum"] == stratum]
        stratum_outcomes[stratum] = {
            "input_cases": len(cases), "expected_input_cases": denominator,
            "passed_states": sum((case["case_id"], depth) in passed_states for case in cases for depth in depths),
            "reference_unavailable_states": sum(case["depth_status"].get(str(depth)) == "criterion_or_reference_unavailable" for case in cases for depth in depths),
            "unattempted_states": sum(case["depth_status"].get(str(depth)) == "unattempted" for case in cases for depth in depths),
        }
    partition_outcomes = {}
    for partition, denominator in PARTITION_DENOMINATORS.items():
        cases = [case for case in selected if case["partition"] == partition]
        partition_outcomes[partition] = {
            "input_cases": len(cases), "expected_input_cases": denominator,
            "passed_states": sum((case["case_id"], depth) in passed_states for case in cases for depth in depths),
            "reference_unavailable_states": sum(case["depth_status"].get(str(depth)) == "criterion_or_reference_unavailable" for case in cases for depth in depths),
            "unattempted_states": sum(case["depth_status"].get(str(depth)) == "unattempted" for case in cases for depth in depths),
        }
    qualified_minima = [depth for depth in minimum_by_case.values() if depth is not None]
    unavailable_cases = {case for case, depth in states
                         if states[(case, depth)].get("status") == "criterion_or_reference_unavailable"}
    return {
        "all_states": {"input_cases": manifest["coverage"]["expected_selected_cases"],
                       "state_slots": manifest["coverage"]["expected_states"]},
        "partition_denominators": PARTITION_DENOMINATORS,
        "stratum_denominators": STRATUM_DENOMINATORS,
        "stratum_outcomes": stratum_outcomes,
        "partition_outcomes": partition_outcomes,
        "value_qualified_states": len(value_states),
        "derivative_qualified_states": len(passed_states),
        "value_only_states": len(value_states - passed_states),
        "oracle_failure_cases": len(unavailable_cases),
        "depth_outcomes": depth_outcomes,
        "smallest_depth_by_case": minimum_by_case,
        "smallest_depth_for_qualified_case": min(qualified_minima, default=None),
        "finite_difference_rows": sum(row.get("comparison_group_id") in
                                      {"centered-gradient", "centered-hessian", "centered-mixed"}
                                      for row in rows),
    }


def write_inputs(path: Path, inputs: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=INPUT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for value in inputs:
            writer.writerow({key: json.dumps(value.get(key), separators=(",", ":"))
                             if isinstance(value.get(key), (list, dict)) else value.get(key, "")
                             for key in INPUT_COLUMNS})


def write_manifest(path: Path, record: dict, selector: str, inputs: list[dict],
    states: dict[tuple[str, int], dict[str, str]]) -> dict:
    depths = [32] if selector == "pilot" else [8, 16, 24, 32]
    expected_input_cases = 110
    expected_selected = {"all": 110, "development": 64, "heldout": 32, "pilot": 1}.get(
        selector, 1 if selector.startswith("case:") else 110)

    def selected(value: dict) -> bool:
        return selector == "all" or (
            selector == "development" and value.get("partition") == "development") or (
            selector == "heldout" and value.get("partition") == "heldout") or (
            selector.startswith("case:") and value.get("case_id") == selector[5:]) or (
            selector == "pilot" and value.get("case_id") == "development-strong-dense-0")

    cases = []
    for value in inputs:
        case_id = str(value["case_id"])
        is_selected = selected(value)
        status = {str(depth): states.get((case_id, depth), {}).get("status", "unattempted")
                  for depth in depths} if is_selected else {}
        cases.append({"case_id": case_id, "partition": value.get("partition"),
                      "stratum": value.get("stratum"), "site_count": input_count(value.get("site_component")),
                      "selected": is_selected,
                      "depth_status": status})
    attempted = [status for case in cases for status in case["depth_status"].values()
                 if status != "unattempted"]
    selected_cases = len({case["case_id"] for case in cases if case["selected"]})
    missing_input_cases = max(0, expected_selected - selected_cases)
    input_manifest_complete = (
        len(inputs) == expected_input_cases
        and len({value.get("case_id") for value in inputs}) == expected_input_cases
    )
    expected_states = expected_selected * len(depths)
    unattempted_states = sum(
        status == "unattempted"
        for case in cases for status in case["depth_status"].values()
    )
    manifest = {
        "schema_version": 1, "run_id": record["run_id"], "selector": selector,
        "execution": record.get("execution"), "outcome_reason": record.get("outcome_reason"),
        "expected_depths": depths, "expected_input_cases": expected_input_cases, "cases": cases,
        "coverage": {
            "input_cases": len(inputs), "input_manifest_complete": input_manifest_complete,
            "selected_cases": selected_cases, "expected_selected_cases": expected_selected,
            "expected_states": expected_states, "attempted_states": len(attempted),
            "passed_states": attempted.count("passed"),
            "unavailable_states": sum(
                status == "criterion_or_reference_unavailable"
                for case in cases for status in case["depth_status"].values()
            ),
            "unattempted_state_slots": unattempted_states,
            "missing_input_cases": missing_input_cases,
            "missing_input_state_slots": missing_input_cases * len(depths),
            "missing_state_slots": max(0, expected_states - len(attempted)),
            "enumeration_status": "complete" if selected_cases == expected_selected
            and expected_states == len(attempted) else "incomplete",
        },
    }
    path.write_text(json.dumps(manifest, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-id", help="refresh the report from an existing exact run")
    parser.add_argument("--selector", default="all",
                        help="native campaign selector: all, pilot, development, heldout, or case:<id>")
    args = parser.parse_args()
    replay = bool(args.run_id)
    metrics = {}
    if replay:
        run_id = args.run_id
    else:
        run_id, metrics = run_campaign(args.selector)
    record, run_path = run_record(run_id)
    selector = str(record.get("selector") or args.selector)
    if selector != "all":
        return
    RESULTS.mkdir(parents=True, exist_ok=True)
    retained = RESULTS / "association-comparison.csv.gz"
    source_name = record.get("comparisons_artifact")
    source = run_path.parent / str(source_name) if source_name else None
    if not source or not source.is_file():
        raise RuntimeError("run replay source artifact is missing; retained evidence was not overwritten")
    retained.write_bytes(gzip.compress(source.read_bytes().replace(b"\r\n", b"\n"), mtime=0))
    with gzip.open(retained, "rt", newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    verdict_rows = [row for row in rows if row.get("verdict") in {"True", "False"}]
    failed = [row for row in verdict_rows if row["verdict"] != "True"]
    unavailable = [row for row in rows if row.get("evidence_state") == "Unavailable"]
    log = ROOT / str(record.get("output", ""))
    inputs = native_inputs(log) if log.is_file() else list(csv.DictReader((RESULTS / "inputs.csv").open(newline="", encoding="utf-8")))
    states = native_states(log) if log.is_file() else {}
    if log.is_file():
        write_inputs(RESULTS / "inputs.csv", inputs)
    manifest = write_manifest(RESULTS / "campaign-manifest.json", record, selector, inputs, states)
    group_maxima = write_plot(RESULTS / "association-comparison.png", rows)
    validation = native_validation(log, rows, inputs) if log.is_file() else {
        "complete": False, "expected_declarations": 0, "expected_prefix_counts": 0,
        "comparison_ids": len({row.get("comparison_id", "") for row in rows}),
        "errors": ["run replay log is unavailable; native validation could not be recomputed"],
    }
    interpretation = build_interpretation(manifest, states, rows, inputs, bool(validation.get("complete")))
    artifacts = [item for item in record["artifact_identity"].get("artifacts", [])
                 if item.get("kind") in {"wheel", "binary"}
                 and (item.get("kind") == "wheel" or "association_comparison_checks" in item.get("path", ""))]
    source_snapshot = record.get("source_snapshot", {
        "base_revision": None, "base_revision_note": "source provenance unavailable for this replay",
        "input_generation": record["input_generation"], "files": {}})
    resource_limits = record.get("resource_limits", {
        "timeout_seconds": 600, "address_space_bytes": 4 * 1024 * 1024 * 1024, "observed": metrics})
    if replay and not record.get("resource_limits"):
        resource_limits["observed"] = {"status": "unavailable", "reason": "report replay did not measure execution resources"}
    summary = {
        "schema_version": 1,
        "question": "Does the fixed-depth association candidate agree with the implicit reference on values and derivatives?",
        "run_id": record["run_id"],
        "recorded_at": datetime.now(UTC).isoformat(),
        "input_generation": record["input_generation"],
        "source_snapshot": source_snapshot,
        "report_replay": {
            "mode": "retained-run-refresh" if replay else "native-execution",
            "numerical_source_snapshot": "preserved" if replay and record.get("source_snapshot") else "unavailable" if replay else "captured",
            "execution_metrics": "preserved" if replay and record.get("resource_limits") else "measured" if not replay else "unavailable",
        },
        "command": record["command"],
        "artifacts": artifacts,
        "artifact_identity": record["artifact_identity"].get("identity"),
        "retained_csv_sha256": hashlib.sha256(gzip.decompress(retained.read_bytes())).hexdigest(),
        "inputs_csv_sha256": sha256(RESULTS / "inputs.csv"),
        "campaign_manifest_sha256": sha256(RESULTS / "campaign-manifest.json"),
        "plot_sha256": sha256(RESULTS / "association-comparison.png"),
        "selector": selector,
        "resource_limits": resource_limits,
        "counts": {
            "rows": len(rows), "criterion_rows": len(verdict_rows),
            "failed": len(failed), "unavailable": len(unavailable),
            "cases": len(inputs),
            "comparison_groups": len({row["comparison_group_id"] for row in rows}),
            "native_validation": validation,
            "manifest": manifest["coverage"],
        },
        "execution": record["execution"],
        "outcome_reason": record.get("outcome_reason"),
        "error_to_threshold_by_group": group_maxima,
        "outcome": "passed" if record["execution"] == "Passed" and validation.get("complete")
        and manifest["coverage"]["enumeration_status"] == "complete" and not failed and not unavailable
        else "recorded_failure_or_unavailable",
        "seeded_generation": {
            "algorithm": "std::mt19937",
            "seed_field": "inputs.csv:generation_seed",
            "association_inputs_symmetric": True,
        },
        "pilot": matching_pilot(record["artifact_identity"]),
        "interpretation": interpretation,
        "claim_limits": [
            "Synthetic neutral fixed-volume states only; this is not all-mixture qualification.",
            "The production weighted-Hessian exact-zero support failure remains a recorded limitation.",
            "The result does not qualify pressure-closed, equilibrium, third-order, or public-API behavior.",
        ],
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, sort_keys=True, separators=(",", ":")) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

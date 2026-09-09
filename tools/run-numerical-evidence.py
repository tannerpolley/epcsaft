#!/usr/bin/env python3
"""Build, execute, and publish declared numerical checks through one recorder."""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import json
import math
import os
import re
import resource
import signal
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from urllib.parse import unquote

# Direct script execution puts ``tools/`` (rather than the repository root) on
# ``sys.path``.  Keep the recorder importable both from the shell entry point
# and when this module is loaded by the local companion server.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

_inventory_spec = importlib.util.spec_from_file_location(
    "numerical_inventory", Path(__file__).with_name("numerical-inventory.py")
)
if _inventory_spec is None or _inventory_spec.loader is None:
    raise ImportError("cannot load numerical-inventory.py")
numerical_inventory = importlib.util.module_from_spec(_inventory_spec)
_inventory_spec.loader.exec_module(numerical_inventory)

try:
    from tools.test_records import (
        RESULTS,
        RecordStore,
        artifact_identity,
        atomic_json,
        canonical,
        comparison_selector_matches,
        freshness,
        input_generation,
        locked,
        make_comparison,
        now,
        sha256,
        validate_run,
        verify_artifact_identity,
    )
except ModuleNotFoundError:  # defensive fallback for callers that set tools/ as cwd
    from test_records import (
        RESULTS,
        RecordStore,
        artifact_identity,
        atomic_json,
        canonical,
        comparison_selector_matches,
        freshness,
        input_generation,
        locked,
        make_comparison,
        now,
        sha256,
        validate_run,
        verify_artifact_identity,
    )


ALIASES = {
    "F1": "native.model.f1",
    "F5": "native.equilibrium.f5",
}
PREPARATION_NAME = "preparation.json"
LOCK_HELD_ENV = "EPCSAFT_RECORD_LOCK_HELD"
_OWNED_PROCESSES: set[subprocess.Popen[Any]] = set()
_SHUTDOWN_REQUESTED = False


def _run_id(test_id: str, case_id: str = "single") -> str:
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", test_id)
    return f"{stamp}-{safe}-{case_id}"


_FLOAT_PREFIX = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?")


def _float_field(value: str) -> float:
    match = _FLOAT_PREFIX.match(value.strip())
    if match is None:
        raise ValueError(value)
    return float(match.group(0))


def _native_fields(line: str, prefix: str) -> dict[str, str]:
    """Read key/value fields while preserving spaces inside a value.

    Native helpers intentionally emit a tiny line protocol.  A value such as a
    scientific claim may contain spaces, so a plain ``split`` or a regex that
    stops at the next space silently truncates it.  Keys are the only tokens
    that have an equals sign immediately after a whitespace boundary.
    """

    text = line[len(prefix):].strip()
    matches = list(re.finditer(r"(?:^|\s)([A-Za-z_][A-Za-z0-9_]*)=", text))
    fields: dict[str, str] = {}
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        value = text[match.end():end].strip()
        if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
            value = value[1:-1]
        fields[match.group(1)] = unquote(value)
    return fields


_SCALAR_RULES = {"symmetric_close", "reference_scale", "reference_or_one"}


def _criterion(value: Any, fields: dict[str, Any]) -> dict[str, Any] | None:
    """Normalize native criterion fields without discarding producer metadata."""

    if value is None:
        return None
    if isinstance(value, dict):
        return dict(value)
    if not isinstance(value, str) or not value:
        return None
    result: dict[str, Any] = {
        "type": "scalar" if value in _SCALAR_RULES else value,
        "rule": fields.get("rule", value),
    }
    for key in ("atol", "rtol", "scale", "norm", "tolerance_rationale", "verdict"):
        if key not in fields:
            continue
        item = fields[key]
        if key in {"atol", "rtol"}:
            try:
                item = _float_field(str(item))
            except ValueError:
                pass
        elif key == "verdict" and isinstance(item, str):
            item = item.lower() == "true"
        result[key] = item
    return result


def _json_comparison(payload: dict[str, Any], ordinal: int) -> dict[str, Any]:
    """Convert the lossless native JSON payload to the shared comparison row."""

    comparison_id = payload.get("comparison_id")
    if not isinstance(comparison_id, str) or not comparison_id:
        return make_comparison(
            f"native.json.invalid.{ordinal}", quantity="native comparison", unit="", inputs={},
            unavailable_reason="native JSON comparison omitted comparison_id",
        )
    raw_judgment = payload.get("judgment")
    judgment = "Observation only" if raw_judgment == "observation" else raw_judgment
    explicit_observation = judgment == "Observation only"
    criterion = _criterion(payload.get("criterion"), payload)
    model = payload.get("model")
    reference_present = "reference" in payload
    reference = payload.get("reference")
    finite_model = (
        isinstance(model, (int, float))
        and not isinstance(model, bool)
        and math.isfinite(model)
    )
    finite_reference = (
        isinstance(reference, (int, float))
        and not isinstance(reference, bool)
        and math.isfinite(reference)
    )
    # Native observations may omit a reference (or serialize it as null), but
    # the producer must say so explicitly.  An unlabelled finite model is not
    # silently promoted to an observation because that would hide a malformed
    # judged comparison.
    observation = explicit_observation and criterion is None and finite_model and (
        not reference_present or reference is None or finite_reference
    )
    inputs = payload.get("inputs")
    if not isinstance(inputs, dict):
        inputs = {}
    coordinates = payload.get("coordinates")
    if not isinstance(coordinates, dict):
        coordinates = None
    row = make_comparison(
        comparison_id,
        quantity=str(payload.get("quantity", "native comparison")),
        unit=str(payload.get("unit", "")),
        inputs=inputs,
        model=model,
        reference=reference,
        criterion=criterion,
        basis=payload.get("basis"),
        method=payload.get("method"),
        claim=payload.get("claim"),
        coordinates=coordinates,
        comparison_group_id=payload.get("comparison_group_id"),
        point_id=payload.get("point_id"),
        quantity_latex=payload.get("quantity_latex"),
        relationship_latex=payload.get("relationship_latex"),
        inputs_latex=payload.get("inputs_latex"),
        observation=observation,
    )
    # Keep structured producer fields, including future fields unknown to this
    # reader.  The normalized criterion/judgment remain in the shared shape.
    for key, value in payload.items():
        if key not in {"criterion", "judgment", "evidence_state", "model", "reference", "inputs", "coordinates"}:
            row[key] = value
    if isinstance(payload.get("criterion"), str):
        row["criterion_name"] = payload["criterion"]
    if criterion is not None:
        row["criterion"] = criterion
    def unavailable(reason: str) -> None:
        row["judgment"] = "Unavailable"
        row["evidence_state"] = "Unavailable"
        row["unavailable_reason"] = reason

    def complete_criterion(value: Any) -> bool:
        if not isinstance(value, dict):
            return False
        if value.get("type", "scalar") == "scalar":
            return all(
                isinstance(value.get(name), (int, float))
                and not isinstance(value.get(name), bool)
                and math.isfinite(value[name])
                and value[name] >= 0
                for name in ("atol", "rtol")
            )
        return bool(value.get("rule")) and isinstance(value.get("verdict"), bool)

    if explicit_observation:
        if criterion is not None:
            unavailable("native observation unexpectedly included a criterion")
        elif not finite_model:
            unavailable("native observation model is missing or nonfinite")
        elif reference_present and reference is not None and not finite_reference:
            unavailable("native observation reference is nonfinite")
        elif not observation:
            unavailable("native observation is missing a finite model")
        else:
            row["judgment"] = "Observation only"
            row["evidence_state"] = "Complete"
    elif judgment in {"Criterion satisfied", "Criterion failed"}:
        # A judged row must carry finite operands, a complete criterion and a
        # producer verdict.  Keep malformed rows unavailable so validation
        # cannot turn a missing reference into scientific evidence.
        verdict = payload.get("verdict")
        if not finite_model or not reference_present or not finite_reference:
            unavailable("native criterion requires finite model and reference")
        elif criterion is None or not complete_criterion(criterion) or not isinstance(verdict, bool):
            unavailable("native criterion omitted complete criterion or verdict")
        elif payload.get("evidence_state") not in {None, "Complete"}:
            unavailable(str(payload.get("unavailable_reason", "native criterion is unavailable")))
        else:
            row["judgment"] = judgment
            row["evidence_state"] = "Complete"
            row["verdict"] = verdict
    elif judgment is not None:
        unavailable(f"native comparison has invalid judgment: {judgment}")
    elif not finite_model:
        unavailable("native comparison model is missing or nonfinite")
    elif criterion is not None:
        unavailable("native criterion omitted judgment")
    else:
        unavailable("native comparison omitted explicit Observation only judgment")
    if "evidence_state" in payload and payload.get("evidence_state") == "Unavailable":
        # Preserve a producer's explicit nonfinite/error state.  In
        # particular, do not turn an observation carrying NaN into a complete
        # row merely because its judgment says Observation only.
        unavailable(str(payload.get("unavailable_reason", "native comparison is unavailable")))
    _fail_closed_native_row(row)
    return row


def _fail_closed_native_row(row: dict[str, Any]) -> None:
    """Do not expose a judged row when scientific metadata is incomplete."""

    if row.get("judgment") not in {"Criterion satisfied", "Criterion failed"}:
        return
    missing = [name for name in ("basis", "method", "claim") if not row.get(name)]
    if not missing:
        return
    row["judgment"] = "Unavailable"
    row["evidence_state"] = "Unavailable"
    row["unavailable_reason"] = "native criterion omitted required metadata: " + ", ".join(missing)


def parse_native_comparisons(output: str) -> list[dict[str, Any]]:
    """Parse the lossless JSON native comparison records."""

    rows: list[dict[str, Any]] = []
    invalid_json = 0
    for line in output.splitlines():
        if not line.startswith("EPCSAFT_COMPARE_JSON "):
            continue
        try:
            payload = json.loads(line[len("EPCSAFT_COMPARE_JSON "):])
        except json.JSONDecodeError:
            invalid_json += 1
            rows.append(make_comparison(
                f"native.json.invalid.{invalid_json}", quantity="native comparison", unit="", inputs={},
                unavailable_reason="native JSON comparison is invalid",
            ))
            continue
        if not isinstance(payload, dict):
            invalid_json += 1
            rows.append(make_comparison(
                f"native.json.invalid.{invalid_json}", quantity="native comparison", unit="", inputs={},
                unavailable_reason="native JSON comparison is not an object",
            ))
            continue
        try:
            row = _json_comparison(payload, invalid_json)
        except (KeyError, TypeError, ValueError) as error:
            invalid_json += 1
            row = make_comparison(
                f"native.json.invalid.{invalid_json}", quantity="native comparison", unit="", inputs={},
                unavailable_reason=f"native JSON comparison is invalid: {error}",
            )
        rows.append(row)
    return rows


def _expectation_record(payload: dict[str, Any], keys: set[str], counts: dict[str, int], errors: list[str]) -> None:
    """Collect exact keys and prefix counts from one generic expectation."""

    exact = payload.get("comparison_id")
    prefix = payload.get("comparison_id_prefix")
    raw_count = payload.get("expected_count")
    count: int | None = None
    if raw_count is not None:
        try:
            count = int(raw_count)
        except (TypeError, ValueError):
            errors.append(f"invalid expected comparison count: {raw_count!r}")
        else:
            if count < 0:
                errors.append(f"invalid expected comparison count: {raw_count!r}")
    if isinstance(exact, str) and exact:
        if count in {None, 1}:
            keys.add(exact)
        elif count is not None:
            counts[exact] = count
    if isinstance(prefix, str) and prefix:
        if count is None:
            errors.append(f"expected comparison prefix has no count: {prefix}")
        else:
            counts[prefix if prefix.endswith("*") else prefix + "*"] = count
    for key in payload.get("expected_comparison_keys", []):
        keys.add(str(key))
    raw_counts = payload.get("expected_comparison_counts", {})
    if isinstance(raw_counts, dict):
        for key, value in raw_counts.items():
            try:
                counts[str(key)] = int(value)
            except (TypeError, ValueError):
                errors.append(f"invalid expected comparison count: {value!r}")


def parse_native_expectations(output: str) -> tuple[list[str], dict[str, int], list[str]]:
    """Read generic exact/prefix coverage records from native output."""

    keys: set[str] = set()
    counts: dict[str, int] = {}
    errors: list[str] = []
    invalid_json = 0
    for line in output.splitlines():
        if line.startswith("EPCSAFT_EXPECT_JSON "):
            try:
                payload = json.loads(line[len("EPCSAFT_EXPECT_JSON "):])
            except json.JSONDecodeError:
                invalid_json += 1
                errors.append("invalid EPCSAFT_EXPECT_JSON record")
                continue
            if not isinstance(payload, dict):
                errors.append("EPCSAFT_EXPECT_JSON record is not an object")
                continue
            _expectation_record(payload, keys, counts, errors)
    return sorted(keys), counts, errors


def parse_native_expectation_groups(output: str) -> list[dict[str, Any]]:
    """Retain expectation definitions that can be shared by emitted rows."""

    groups: list[dict[str, Any]] = []
    seen: set[str] = set()
    for line in output.splitlines():
        if not line.startswith("EPCSAFT_EXPECT_JSON "):
            continue
        try:
            payload = json.loads(line[len("EPCSAFT_EXPECT_JSON "):])
        except json.JSONDecodeError:
            continue
        if not isinstance(payload, dict):
            continue
        prefix = payload.get("comparison_id_prefix")
        if not isinstance(prefix, str) or not prefix:
            # Exact one-row expectations do not define a shared case.
            continue
        group_id = payload.get("comparison_group_id")
        group_id = group_id if isinstance(group_id, str) and group_id else None
        try:
            expected_count = int(payload["expected_count"])
        except (KeyError, TypeError, ValueError):
            expected_count = None
        key = canonical({
            "comparison_group_id": group_id,
            "comparison_id_prefix": prefix,
            "expected_count": expected_count,
            "definition": payload.get("definition"),
            "coordinates": payload.get("coordinates"),
        })
        if key in seen:
            continue
        seen.add(key)
        group: dict[str, Any] = {
            "group_ref": f"g{len(groups)}",
            "shared_inputs": {},
            "comparison_id_prefix": prefix,
        }
        if group_id is not None:
            group["comparison_group_id"] = group_id
        if expected_count is not None and expected_count >= 0:
            group["expected_count"] = expected_count
        if "definition" in payload:
            group["case_definition"] = payload["definition"]
        if isinstance(payload.get("coordinates"), dict):
            group["coordinates"] = payload["coordinates"]
        groups.append(group)
    return groups


def _factor_input_groups(
    rows: list[dict[str, Any]], expectation_groups: list[dict[str, Any]] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Factor equal row inputs while retaining every differing override."""

    groups = [dict(group) for group in (expectation_groups or [])]
    assignments: dict[str, list[int]] = {str(group["group_ref"]): [] for group in groups}
    fallback_by_id: dict[str, dict[str, Any]] = {}

    def matches(row: dict[str, Any], group: dict[str, Any]) -> bool:
        row_group = row.get("comparison_group_id")
        group_id = group.get("comparison_group_id")
        if group_id is not None and row_group not in {None, group_id}:
            return False
        comparison_id = row.get("comparison_id")
        prefix = group.get("comparison_id_prefix")
        return isinstance(comparison_id, str) and isinstance(prefix, str) and comparison_id.startswith(prefix)

    for index, row in enumerate(rows):
        candidates = [group for group in groups if matches(row, group)]
        if candidates:
            # A longer prefix is more specific when expectation groups overlap.
            selected = max(candidates, key=lambda group: len(str(group.get("comparison_id_prefix", ""))))
        else:
            group_id = row.get("comparison_group_id")
            if not isinstance(group_id, str) or not group_id:
                continue
            selected = fallback_by_id.get(group_id)
            if selected is None:
                selected = {
                    "group_ref": f"g{len(groups)}",
                    "comparison_group_id": group_id,
                    "shared_inputs": {},
                }
                groups.append(selected)
                fallback_by_id[group_id] = selected
                assignments[str(selected["group_ref"])] = []
        assignments.setdefault(str(selected["group_ref"]), []).append(index)

    for group in groups:
        ref = str(group["group_ref"])
        indices = assignments.get(ref, [])
        maps = [rows[index].get("inputs", {}) for index in indices]
        maps = [value for value in maps if isinstance(value, dict)]
        shared: dict[str, Any] = {}
        if maps and len(maps) == len(indices):
            common = set(maps[0])
            for values in maps[1:]:
                common.intersection_update(values)
            definition = group.get("case_definition")
            for key in sorted(common):
                value = maps[0][key]
                if not all(values[key] == value for values in maps[1:]):
                    continue
                if isinstance(definition, dict) and (
                    key not in definition or definition[key] != value
                ):
                    continue
                shared[key] = value
        group["shared_inputs"] = shared
        for index in indices:
            row = rows[index]
            values = row.get("inputs", {})
            if not isinstance(values, dict):
                continue
            row["group_ref"] = ref
            row["inputs"] = {
                key: value for key, value in values.items()
                if key not in shared or value != shared[key]
            }
    return rows, groups


def _parse_native_cases(output: str) -> dict[str, dict[str, str]]:
    """Return the last terminal case event emitted for each native case."""

    result: dict[str, dict[str, str]] = {}
    for line in output.splitlines():
        if not line.startswith("EPCSAFT_CASE "):
            continue
        fields = _native_fields(line, "EPCSAFT_CASE ")
        if fields.get("test_id") and fields.get("execution"):
            result[fields["test_id"]] = {
                "execution": fields["execution"].replace("_", " "),
                "outcome_reason": fields.get("outcome_reason", ""),
            }
    return result


def parse_integrity_results(output: str) -> tuple[list[dict[str, Any]], list[str]]:
    """Read one structured outcome emitted by an integrity-check child."""

    results: list[dict[str, Any]] = []
    errors: list[str] = []
    for line in output.splitlines():
        if not line.startswith("EPCSAFT_INTEGRITY_JSON "):
            continue
        try:
            payload = json.loads(line[len("EPCSAFT_INTEGRITY_JSON "):])
        except json.JSONDecodeError:
            errors.append("invalid integrity outcome")
            continue
        if not isinstance(payload, dict):
            errors.append("integrity outcome is not an object")
            continue
        if not isinstance(payload.get("test_id"), str) or not payload.get("test_id"):
            errors.append("integrity outcome omitted test_id")
            continue
        if payload.get("execution") not in {"Passed", "Failed", "Error", "Skipped"}:
            errors.append("integrity outcome has invalid execution")
            continue
        results.append(payload)
    return results, errors


def _expected_complete(
    declaration: dict[str, Any],
    rows: list[dict[str, Any]],
    *,
    expected_keys: list[str] | None = None,
    expected_counts: dict[str, int] | None = None,
) -> tuple[list[dict[str, Any]], str | None]:
    expected = [str(item) for item in (
        declaration.get("expected_comparison_keys", []) if expected_keys is None else expected_keys
    )]
    emitted = {str(row.get("comparison_id")) for row in rows}
    missing: list[str] = []
    for key in expected:
        if not any(comparison_selector_matches(value, key) for value in emitted):
            missing.append(key)
    count_errors: list[str] = []
    declared_counts = declaration.get("expected_comparison_counts", {})
    counts = declared_counts if expected_counts is None else expected_counts
    if isinstance(counts, dict):
        for prefix, count in counts.items():
            selector = str(prefix)
            try:
                expected_count = int(count)
            except (TypeError, ValueError):
                count_errors.append(f"{prefix} has invalid expected count {count!r}")
                continue
            actual_count = sum(comparison_selector_matches(value, selector) for value in emitted)
            if actual_count != expected_count:
                count_errors.append(f"{prefix} expected {expected_count} rows, emitted {actual_count}")
    if missing:
        for key in missing:
            rows.append(make_comparison(key, quantity=key, unit="", inputs={},
                                        unavailable_reason="expected comparison was not emitted"))
    if count_errors:
        for prefix, count in counts.items():
            selector = str(prefix)
            actual_count = sum(comparison_selector_matches(value, selector) for value in emitted)
            try:
                expected_count = int(count)
            except (TypeError, ValueError):
                continue
            if actual_count != expected_count:
                rows.append(make_comparison(
                    f"{prefix}#coverage", quantity=str(prefix), unit="", inputs={},
                    unavailable_reason=f"expected {expected_count} comparison rows, emitted {actual_count}",
                ))
        return rows, "; ".join(count_errors)
    if missing:
        return rows, "expected numerical comparison was not emitted: " + ", ".join(missing)
    if declaration.get("evidence_expectation") == "Numerical" and not rows:
        return rows, "declared numerical evidence emitted no comparisons"
    return rows, None


def _descendant_pids(root_pid: int) -> set[int]:
    """Capture descendants from procfs before terminating an owned command."""

    parents: dict[int, int] = {}
    for entry in Path("/proc").glob("[0-9]*"):
        try:
            stat = (entry / "stat").read_text(encoding="utf-8")
            parents[int(entry.name)] = int(stat.rsplit(")", 1)[1].split()[1])
        except (OSError, ValueError, IndexError):
            continue
    descendants = {root_pid}
    changed = True
    while changed:
        changed = False
        for pid, parent in parents.items():
            if parent in descendants and pid not in descendants:
                descendants.add(pid)
                changed = True
    return descendants


def _track_process(process: subprocess.Popen[Any]) -> subprocess.Popen[Any]:
    _OWNED_PROCESSES.add(process)
    return process


def _untrack_process(process: subprocess.Popen[Any]) -> None:
    _OWNED_PROCESSES.discard(process)


def _terminate_process_group(process: subprocess.Popen[Any], *, grace: float = 2.0) -> None:
    """Stop one owned process and captured descendants in private sessions."""

    if process.poll() is not None:
        _untrack_process(process)
        return
    pids = _descendant_pids(process.pid)
    groups: set[int] = set()
    for pid in pids:
        try:
            groups.add(os.getpgid(pid))
        except ProcessLookupError:
            continue
    try:
        for pgid in groups:
            os.killpg(pgid, signal.SIGTERM)
    except (OSError, ProcessLookupError):
        pass
    # A descendant may have called setsid(), so terminate its captured PID as
    # well.  These are the exact processes observed below the owned command;
    # no name or system-wide kill is used.
    for pid in pids:
        try:
            os.kill(pid, signal.SIGTERM)
        except (OSError, ProcessLookupError):
            pass
    try:
        process.wait(timeout=grace)
    except subprocess.TimeoutExpired:
        remaining = _descendant_pids(process.pid) | pids
        for pid in remaining:
            try:
                os.kill(pid, signal.SIGKILL)
            except (OSError, ProcessLookupError):
                pass
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except (OSError, ProcessLookupError):
            pass
        process.wait()
    finally:
        _untrack_process(process)


def _relay_child_shutdown(signum: int, _frame: Any) -> None:
    """Relay termination to build/test groups owned by this runner."""

    global _SHUTDOWN_REQUESTED
    _SHUTDOWN_REQUESTED = True
    for process in tuple(_OWNED_PROCESSES):
        _terminate_process_group(process)


def _install_child_signal_handlers() -> None:
    signal.signal(signal.SIGTERM, _relay_child_shutdown)
    signal.signal(signal.SIGINT, _relay_child_shutdown)


def _prepare_environment(root: Path, log: Path) -> int:
    """Create the sole environment only while the recording lock is held."""

    if (root / ".venv/bin/python").is_file():
        return 0
    environment = dict(os.environ, EPCSAFT_PREPARE_NATIVE_BUILD="1", EPCSAFT_RECORD_LOCK_HELD="1")
    with log.open("a", encoding="utf-8") as stream:
        try:
            process = _track_process(subprocess.Popen(
                ["bash", ".codex/environments/setup.sh"], cwd=root, env=environment,
                stdout=stream, stderr=subprocess.STDOUT, start_new_session=True,
            ))
        except OSError:
            return 1
        try:
            returncode = process.wait(timeout=900)
        except subprocess.TimeoutExpired:
            _terminate_process_group(process)
            return 124
        except KeyboardInterrupt:
            _terminate_process_group(process)
            return 130
        finally:
            _untrack_process(process)
    return returncode


def _write_preparation_stamp(
    root: Path, input_id: str, build_dir: Path | None, wheel: Path | None,
) -> None:
    if build_dir is None or wheel is None or not build_dir.is_dir() or not wheel.is_file():
        return
    stamp = root / "build" / "environment-wheel" / PREPARATION_NAME
    atomic_json(stamp, {
        "schema_version": 1,
        "input_generation": input_id,
        "build_dir": str(build_dir.resolve().relative_to(root.resolve())),
        "wheel": str(wheel.resolve().relative_to(root.resolve())),
        "wheel_sha256": sha256(wheel),
    })


def _read_preparation_stamp(root: Path, input_id: str) -> tuple[bool, str | None]:
    stamp = root / "build" / "environment-wheel" / PREPARATION_NAME
    try:
        value = json.loads(stamp.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        return False, f"prepared artifact binding unavailable: {error}"
    if not isinstance(value, dict) or value.get("schema_version") != 1:
        return False, "prepared artifact binding has an invalid schema"
    if value.get("input_generation") != input_id:
        return False, "prepared artifact was built from a different input generation"
    wheel = (root / str(value.get("wheel", ""))).resolve()
    build_dir = (root / str(value.get("build_dir", ""))).resolve()
    if not wheel.is_file() or not wheel.is_relative_to(root.resolve()):
        return False, "prepared Engine wheel is unavailable"
    if not build_dir.is_dir() or not build_dir.is_relative_to(root.resolve()):
        return False, "prepared native build directory is unavailable"
    if sha256(wheel) != value.get("wheel_sha256"):
        return False, "prepared Engine wheel changed"
    marker = root / "build/environment-wheel/meson-build-dir"
    if not marker.is_file() or Path(marker.read_text().strip()).resolve() != build_dir:
        return False, "prepared native build directory changed"
    return True, None


def _build(root: Path, log: Path, input_id: str) -> tuple[int, Path | None, Path | None]:
    environment = dict(os.environ, EPCSAFT_PREPARE_NATIVE_BUILD="1", EPCSAFT_RECORD_LOCK_HELD="1")
    if _prepare_environment(root, log):
        return 1, None, None
    with log.open("a", encoding="utf-8") as stream:
        try:
            process = _track_process(subprocess.Popen(
                ["bash", "tools/ensure-engine-wheel.sh"], cwd=root,
                env=environment, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True,
            ))
        except OSError:
            return 1, None, None
        try:
            returncode = process.wait(timeout=1800)
        except subprocess.TimeoutExpired:
            _terminate_process_group(process)
            return 124, None, None
        except KeyboardInterrupt:
            _terminate_process_group(process)
            return 130, None, None
        finally:
            _untrack_process(process)
    if returncode:
        return returncode, None, None
    marker = root / "build/environment-wheel/meson-build-dir"
    build_dir = Path(marker.read_text().strip()) if marker.is_file() else None
    wheels = sorted((root / "build/environment-wheel").glob("epcsaft-*.whl"))
    wheel = wheels[0] if len(wheels) == 1 else None
    if input_generation(root) == input_id:
        _write_preparation_stamp(root, input_id, build_dir, wheel)
    return returncode, build_dir, wheel


def _run_subprocess(
    command: list[str], *, cwd: Path, timeout: int, environment: dict[str, str],
    output: Path, error_output: Path,
) -> tuple[str, int | None, str | None, str | None]:
    output.parent.mkdir(parents=True, exist_ok=True)
    error_output.parent.mkdir(parents=True, exist_ok=True)
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    process: subprocess.Popen[Any] | None = None
    try:
        with output.open("w", encoding="utf-8") as stream, error_output.open("w", encoding="utf-8") as error_stream:
            process = _track_process(subprocess.Popen(
                command, cwd=cwd, env=environment, stdout=stream,
                stderr=error_stream, start_new_session=True,
            ))
            try:
                returncode = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                _terminate_process_group(process)
                return "Interrupted", None, "process timed out", "timeout"
        if _SHUTDOWN_REQUESTED:
            return "Interrupted", None, "runner interrupted by shutdown", "interruption"
        return ("Passed" if returncode == 0 else "Failed"), returncode, (
            None if returncode == 0 else f"process exit={returncode}"
        ), None
    except subprocess.TimeoutExpired:
        if process is not None:
            _terminate_process_group(process)
        return "Interrupted", None, "process timed out", "timeout"
    except KeyboardInterrupt:
        if process is not None:
            _terminate_process_group(process)
        return "Interrupted", None, "runner interrupted", "interruption"
    finally:
        if process is not None:
            _untrack_process(process)


def _run_record(
    declaration: dict[str, Any],
    *,
    run_id: str | None = None,
    started_at: str | None = None,
    input_id: str,
    artifact: dict[str, Any],
    build_dir: Path | None,
    wheel: Path | None,
    timeout: int,
    record_dir: Path,
) -> dict[str, Any]:
    test_id = str(declaration["test_id"])
    case_id = str(declaration.get("_parameter_case_id", declaration.get("parameter_case_id", "single")))
    run_id = run_id or _run_id(test_id, case_id)
    output = record_dir / f"{run_id}.log"
    error_output = record_dir / f"{run_id}.stderr.log"
    started = started_at or now()
    runner = declaration.get("runner", {})
    command: list[str]
    run_cwd = ROOT
    if runner.get("kind") == "native":
        target = build_dir / str(runner.get("target")) if build_dir else None
        if target is None or not target.is_file():
            return {
                "schema_version": 1, "run_id": run_id, "test_id": test_id,
                "parameter_case_id": case_id, "input_generation": input_id,
                "artifact_identity": artifact, "execution": "Error",
                "started_at": started, "finished_at": now(),
                "outcome_reason": "native executable unavailable",
                "comparisons": [], "input_groups": [], "command": [],
                "artifact_digests": artifact.get("artifacts", []),
            }
        command = [str(target), *[str(item) for item in runner.get("arguments", [])]]
    elif runner.get("kind") == "integrity":
        script = ROOT / str(runner.get("script", "tools/check-test-system.py"))
        check_id = runner.get("check_id")
        if not script.is_file() or not isinstance(check_id, str) or not check_id:
            return {
                "schema_version": 1, "run_id": run_id, "test_id": test_id,
                "parameter_case_id": case_id, "input_generation": input_id,
                "artifact_identity": artifact, "execution": "Error",
                "started_at": started, "finished_at": now(),
                "outcome_reason": "integrity check executable unavailable",
                "comparisons": [], "input_groups": [], "command": [],
                "artifact_digests": artifact.get("artifacts", []),
            }
        command = [str(ROOT / ".venv/bin/python"), str(script), "--check-id", check_id]
    elif runner.get("kind") == "pytest":
        command = [
            str(ROOT / ".venv/bin/python"), "-m", "pytest", "-q", "-p", "tools.pytest_recording",
            # Explicit evidence selection must be able to run a slow or
            # expensive case even though the repository default excludes it.
            "-m", "slow or not slow", str(runner.get("nodeid", "")),
        ]
    elif runner.get("kind") == "study":
        files = [str(item) for item in runner.get("test_files", [])]
        if not files:
            return {
                "schema_version": 1, "run_id": run_id, "test_id": test_id,
                "parameter_case_id": case_id, "input_generation": input_id,
                "artifact_identity": artifact, "execution": "Skipped",
                "started_at": started, "finished_at": now(),
                "outcome_reason": "study has no executable test entry point",
                "comparisons": [], "input_groups": [], "command": [],
                "artifact_digests": artifact.get("artifacts", []),
            }
        command = [
            str(ROOT / ".venv/bin/python"), "-m", "pytest", "-q",
            "-m", "slow or not slow", *files,
        ]
        run_cwd = ROOT / str(runner.get("cwd", "."))
    else:
        return {
            "schema_version": 1, "run_id": run_id, "test_id": test_id,
            "parameter_case_id": case_id, "input_generation": input_id,
            "artifact_identity": artifact, "execution": "Skipped",
            "started_at": started, "finished_at": now(),
            "outcome_reason": "declaration has no executable runner",
            "comparisons": [], "input_groups": [], "command": [],
            "artifact_digests": artifact.get("artifacts", []),
        }
    environment = dict(
        os.environ,
        PYTHONNOUSERSITE="1",
        EPCSAFT_RECORD_ATTEMPT_ACTIVE="1",
    )
    pytest_result = record_dir / f"{run_id}.pytest.json"
    if runner.get("kind") == "pytest":
        environment["EPCSAFT_PYTEST_RESULTS"] = str(pytest_result)
        # Collection is performed from ``engine/`` so its node IDs match the
        # repository's normal pytest entry points and package imports.  Keep
        # the repository root on PYTHONPATH for the recording plugin.
        run_cwd = ROOT / "engine"
        existing_pythonpath = environment.get("PYTHONPATH")
        environment["PYTHONPATH"] = str(ROOT) + (os.pathsep + existing_pythonpath if existing_pythonpath else "")
    elif runner.get("kind") == "study":
        # Study replays are owned by this runner's parent record.  Suppress
        # direct pytest publication in the child so it cannot reacquire the
        # parent's shared lock or create orphan test rows.
        environment["EPCSAFT_PYTEST_RESULTS"] = str(pytest_result)
    execution, exit_code, process_reason, interruption_kind = _run_subprocess(
        command, cwd=run_cwd, timeout=timeout, environment=environment,
        output=output, error_output=error_output,
    )
    process_interrupted = execution == "Interrupted"
    text = output.read_text(encoding="utf-8", errors="replace") if output.exists() else ""
    rows: list[dict[str, Any]] = []
    native_expected_keys: list[str] = []
    native_expected_counts: dict[str, int] = {}
    native_expectation_errors: list[str] = []
    native_case: dict[str, str] | None = None
    native_expectation_groups: list[dict[str, Any]] = []
    integrity_errors: list[str] = []
    integrity_result: dict[str, Any] | None = None
    integrity_reason: str | None = None
    if runner.get("kind") == "native":
        rows = parse_native_comparisons(text)
        native_expected_keys, native_expected_counts, native_expectation_errors = parse_native_expectations(text)
        native_expectation_groups = parse_native_expectation_groups(text)
        native_case = _parse_native_cases(text).get(test_id)
    elif runner.get("kind") == "integrity":
        integrity_results, integrity_errors = parse_integrity_results(text)
        matching = [item for item in integrity_results if item.get("test_id") == test_id]
        if len(integrity_results) != 1 or len(matching) != 1:
            integrity_errors.append("integrity check did not emit exactly one matching outcome")
        else:
            integrity_result = matching[0]
    else:
        pass
    reporter_error: str | None = None
    if runner.get("kind") == "pytest":
        if not pytest_result.is_file():
            reporter_error = "pytest reporter did not write results"
        else:
            try:
                report_payload = json.loads(pytest_result.read_text())
                entries = report_payload.get("results", [])
                if not isinstance(entries, list):
                    raise TypeError("pytest reporter results are not a list")
                entry = next((item for item in entries if isinstance(item, dict)
                              and item.get("test_id") == test_id
                              and item.get("parameter_case_id", "single") == case_id), None)
                if entry is not None:
                    if not process_interrupted:
                        execution = entry.get("execution", execution)
                        exit_code = 0 if execution == "Passed" else 1
                    comparisons = entry.get("comparisons", [])
                    if not isinstance(comparisons, list):
                        raise ValueError("pytest reporter comparisons are not a list")
                    rows = list(comparisons)
                else:
                    reporter_error = "pytest reporter did not emit the selected test"
            except (OSError, TypeError, ValueError, json.JSONDecodeError) as error:
                reporter_error = f"pytest reporter output is invalid: {error}"
        if reporter_error and not process_interrupted:
            execution = "Error"
            exit_code = 1
    native_case_reason: str | None = None
    if native_case is not None:
        case_execution = native_case.get("execution")
        if not process_interrupted and case_execution in {"Passed", "Failed", "Skipped", "Error", "Interrupted"}:
            execution = case_execution
            exit_code = 0 if case_execution in {"Passed", "Skipped"} else 1
            native_case_reason = native_case.get("outcome_reason") or None
        elif not process_interrupted and case_execution == "Running":
            execution = "Error"
            exit_code = None
            native_case_reason = "native case did not emit a terminal outcome"
        elif process_interrupted and case_execution == "Running":
            native_case_reason = "native case did not emit a terminal outcome"
    if runner.get("kind") == "integrity":
        if integrity_result is not None and not integrity_errors:
            integrity_execution = str(integrity_result["execution"])
            if process_interrupted:
                if integrity_result.get("outcome_reason"):
                    integrity_reason = str(integrity_result["outcome_reason"])
            elif integrity_execution in {"Passed", "Skipped"} and execution != "Passed":
                integrity_errors.append("integrity outcome disagreed with process status")
            elif integrity_execution in {"Failed", "Error"} and execution == "Passed":
                integrity_errors.append("integrity process passed after a failed outcome")
            else:
                execution = integrity_execution
                exit_code = 0 if integrity_execution in {"Passed", "Skipped"} else 1
                if integrity_result.get("outcome_reason"):
                    integrity_reason = str(integrity_result["outcome_reason"])
        if integrity_errors:
            integrity_reason = "; ".join(integrity_errors)
    rows, evidence_error = _expected_complete(
        declaration,
        rows,
        expected_keys=native_expected_keys or None,
        expected_counts=native_expected_counts or None,
    )
    if integrity_reason:
        evidence_error = integrity_reason if evidence_error is None else f"{integrity_reason}; {evidence_error}"
    if native_expectation_errors:
        expectation_error = "; ".join(native_expectation_errors)
        evidence_error = expectation_error if evidence_error is None else f"{expectation_error}; {evidence_error}"
    if native_case_reason:
        evidence_error = native_case_reason if evidence_error is None else f"{native_case_reason}; {evidence_error}"
    if reporter_error:
        evidence_error = reporter_error if evidence_error is None else f"{reporter_error}; {evidence_error}"
    if evidence_error:
        execution = "Error" if execution == "Passed" else execution
    rows, input_groups = _factor_input_groups(rows, native_expectation_groups)
    current_input = input_generation(ROOT)
    if current_input != input_id:
        if not process_interrupted:
            execution = "Error"
        evidence_error = (
            "computational inputs changed during execution"
            if evidence_error is None
            else f"{evidence_error}; computational inputs changed during execution"
        )
    elif execution not in {"Skipped", "Queued", "Running"}:
        identity_ok, identity_reason = verify_artifact_identity(ROOT, artifact)
        if not identity_ok:
            if not process_interrupted:
                execution = "Error"
            identity_error = "artifact identity changed during execution: " + str(identity_reason)
            evidence_error = identity_error if evidence_error is None else f"{evidence_error}; {identity_error}"
    reason = "; ".join(str(value) for value in (process_reason, evidence_error) if value)
    reason = reason or None
    if reason is None and execution in {"Failed", "Interrupted", "Error", "Skipped"}:
        reason = f"process exit={exit_code}" if execution == "Failed" else "process did not complete"
    run = {
        "schema_version": 1,
        "run_id": run_id,
        "test_id": test_id,
        "parameter_case_id": case_id,
        "input_generation": input_id,
        "artifact_identity": artifact,
        "execution": execution,
        "started_at": started,
        "finished_at": now(),
        "outcome_reason": reason,
        "comparisons": rows,
        "input_groups": input_groups,
        "command": command,
        "exit_code": exit_code,
        "output": str(output.relative_to(ROOT)),
        "output_sha256": __import__("hashlib").sha256(output.read_bytes()).hexdigest() if output.is_file() else None,
        "stderr": str(error_output.relative_to(ROOT)),
        "stderr_sha256": __import__("hashlib").sha256(error_output.read_bytes()).hexdigest()
        if error_output.is_file() else None,
    }
    if interruption_kind is not None:
        run["interruption_kind"] = interruption_kind
    run["artifact_digests"] = artifact.get("artifacts", [])
    return run


def _attempt_record(
    declaration: dict[str, Any],
    *,
    run_id: str,
    input_id: str,
    artifact: dict[str, Any],
    execution: str,
    started_at: str,
    command: list[str] | None = None,
) -> dict[str, Any]:
    """Build the durable record written before an owned child is started."""

    return {
        "schema_version": 1,
        "run_id": run_id,
        "test_id": str(declaration["test_id"]),
        "parameter_case_id": str(declaration.get("_parameter_case_id", "single")),
        "input_generation": input_id,
        "artifact_identity": artifact,
        "execution": execution,
        "started_at": started_at,
        "finished_at": now(),
        "outcome_reason": None,
        "comparisons": [],
        "input_groups": [],
        "artifact_digests": artifact.get("artifacts", []),
        "command": command or [],
        "exit_code": None,
    }


def _validated_final_run(run: dict[str, Any]) -> dict[str, Any]:
    """Convert a malformed final attempt into a publishable Error record."""

    errors = validate_run(run)
    if errors:
        run["execution"] = "Error"
        run.pop("interruption_kind", None)
        run["outcome_reason"] = "; ".join(errors)
    return run


def _select(inventory: dict[str, Any], args: argparse.Namespace) -> list[dict[str, Any]]:
    declarations = [item for item in inventory.get("declarations", []) if item.get("lifecycle") == "Implemented"]
    requested = [ALIASES.get(item, item) for item in args.cases]
    native_suite = getattr(args, "native_suite", None)
    if native_suite:
        native_target = numerical_inventory.NATIVE_TARGETS.get(str(native_suite).strip().casefold())
        declarations = [
            item for item in declarations
            if item.get("runner", {}).get("kind") == "native"
            and item.get("runner", {}).get("target") == native_target
        ]
        if not declarations:
            raise SystemExit("unknown or unavailable native suite: " + str(native_suite))
    elif requested:
        selected = [item for item in declarations if item.get("test_id") in requested]
        missing = sorted(set(requested) - {item.get("test_id") for item in selected})
        if missing:
            raise SystemExit("unknown declared test id: " + ", ".join(missing))
        declarations = selected
        # An explicitly named case is an explicit request even when its
        # declaration is classified Automatic fast.  Selection policy only
        # filters an unqualified suite request.
    elif args.explicit:
        declarations = [item for item in declarations if item.get("execution_policy") == "Explicit"]
    else:
        declarations = [item for item in declarations
                        if item.get("execution_policy") == "Automatic fast"
                        and item.get("runner", {}).get("kind") in {"native", "pytest", "integrity"}
                        and not item.get("unclassified")]
    expanded: list[dict[str, Any]] = []
    for item in declarations:
        cases = item.get("parameter_cases") or [{
            "parameter_case_id": item.get("parameter_case_id", "single"),
            "nodeid": item.get("runner", {}).get("nodeid"),
        }]
        for case in cases:
            clone = dict(item)
            clone["_parameter_case_id"] = str(case.get("parameter_case_id", "single"))
            if clone.get("runner", {}).get("kind") == "pytest" and case.get("nodeid"):
                clone["runner"] = dict(clone["runner"], nodeid=str(case["nodeid"]))
            expanded.append(clone)
    return expanded


def _coalesced_cases(
    store: RecordStore,
    declarations: list[dict[str, Any]],
    input_id: str,
    *,
    automatic: bool,
) -> tuple[list[dict[str, Any]], bool]:
    """Drop reusable automatic cases and retain their failure status.

    A timeout is a terminal failed attempt for this input generation.  An
    infrastructure interruption (including a killed runner whose active
    attempt was reconciled) is incomplete and must be retried on the next
    automatic invocation.
    """

    if not automatic:
        return declarations, False
    # The manifest is only an index.  Read and validate the pointed-to run
    # records before reusing a terminal result; otherwise a stale pointer or
    # a replaced artifact can make a failed/unverified case look like a clean
    # coalesced pass.
    previous = {
        (item.get("test_id"), item.get("parameter_case_id")): item
        for item in store.read_runs()
        if isinstance(item, dict)
    }
    applicability: dict[str, tuple[str, str | None]] = {}
    remaining: list[dict[str, Any]] = []
    coalesced_failed = False
    for declaration in declarations:
        key = (declaration.get("test_id"), declaration.get("_parameter_case_id", "single"))
        run = previous.get(key)
        if run is None or run.get("input_generation") != input_id:
            remaining.append(declaration)
            continue
        try:
            artifact_key = json.dumps(run.get("artifact_identity"), sort_keys=True, allow_nan=False)
        except (TypeError, ValueError):
            artifact_key = f"run:{run.get('run_id', key)}"
        if artifact_key not in applicability:
            applicability[artifact_key] = freshness(run, store.root)
        if applicability[artifact_key][0] != "Current":
            remaining.append(declaration)
            continue
        execution = run.get("execution")
        if execution not in {"Passed", "Failed", "Skipped", "Error", "Interrupted"}:
            remaining.append(declaration)
            continue
        if execution == "Interrupted" and run.get("interruption_kind") != "timeout":
            remaining.append(declaration)
            continue
        coalesced_failed |= execution not in {"Passed", "Skipped"}
    return remaining, coalesced_failed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="*", help="stable test IDs or compatibility aliases F1/F5")
    parser.add_argument("--automatic", action="store_true", help="run the permitted automatic-fast selection")
    parser.add_argument("--explicit", action="store_true", help="run explicit declarations once")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--no-build", action="store_true", help="use the already prepared installed wheel")
    parser.add_argument(
        "--native-suite",
        help="run every registered case in one native suite (contributions/properties)",
    )
    args = parser.parse_args()
    if not args.cases and not args.automatic and not args.explicit:
        args.automatic = True
    _install_child_signal_handlers()
    RESULTS.mkdir(parents=True, exist_ok=True)
    store = RecordStore(ROOT)
    lock = contextlib.nullcontext() if os.environ.get(LOCK_HELD_ENV) == "1" else locked(store.results / ".lock")
    with lock:
        # Close attempts left in Queued/Running by a killed or crashed
        # predecessor before selecting automatic work.  Their interrupted
        # record replaces any older pass in the manifest.
        store.reconcile_attempts()
        input_id = input_generation(ROOT)
        build_log = RESULTS / "build.log"
        if args.no_build:
            marker = ROOT / "build/environment-wheel/meson-build-dir"
            build_dir = Path(marker.read_text().strip()) if marker.is_file() else None
            wheels = sorted((ROOT / "build/environment-wheel").glob("epcsaft-*.whl"))
            code, wheel = (0, wheels[0]) if len(wheels) == 1 else (1, None)
            prepared, preparation_reason = _read_preparation_stamp(ROOT, input_id)
        else:
            code, build_dir, wheel = _build(ROOT, build_log, input_id)
            prepared, preparation_reason = code == 0, None if code == 0 else (
                "runner interrupted during wheel preparation" if _SHUTDOWN_REQUESTED
                else "wheel preparation failed"
            )
        # Native --list and pytest collection must observe the just-prepared
        # build.  Discovery before preparation can silently select stale cases.
        inventory = numerical_inventory.inventory(write=True)
        selected = _select(inventory, args)
        if not selected:
            if code or not prepared:
                print(preparation_reason or "wheel preparation failed; no executable cases selected", file=sys.stderr)
                return 1
            print("no declared cases selected")
            return 0
        coalesced_failed = False
        if args.automatic and not args.cases:
            selected, coalesced_failed = _coalesced_cases(store, selected, input_id, automatic=True)
            if not selected:
                print("automatic cases already recorded for this input generation")
                return int(coalesced_failed)
        if code:
            artifact = artifact_identity(ROOT, unavailable_reason="wheel preparation failed; see build.log")
        elif not prepared:
            artifact = artifact_identity(ROOT, unavailable_reason=preparation_reason or "prepared artifact binding unavailable")
        else:
            artifact = artifact_identity(
                ROOT,
                build_dir=build_dir,
                wheel=wheel,
                runtime=ROOT / ".venv/bin/python",
                inventory_generation=inventory.get("inventory_generation"),
            )
        if input_generation(ROOT) != input_id:
            preparation_reason = "inputs changed during wheel preparation"
            artifact = artifact_identity(ROOT, unavailable_reason=preparation_reason)
        elif code == 0:
            identity_ok, identity_reason = verify_artifact_identity(ROOT, artifact)
            if not identity_ok:
                preparation_reason = "artifact identity unavailable before execution: " + str(identity_reason)
                artifact = artifact_identity(
                    ROOT,
                    unavailable_reason=preparation_reason,
                )
        execution_allowed = code == 0 and artifact.get("state") == "Available"
        failed = coalesced_failed or not execution_allowed
        for declaration in selected:
            case_execution_allowed = execution_allowed and not _SHUTDOWN_REQUESTED
            run_id = _run_id(str(declaration["test_id"]), str(declaration.get("_parameter_case_id", "single")))
            started_at = now()
            attempt = _attempt_record(
                declaration, run_id=run_id, input_id=input_id, artifact=artifact,
                execution="Queued", started_at=started_at,
            )
            store.publish(attempt)
            if case_execution_allowed:
                store.publish({
                    **attempt,
                    "execution": "Running",
                    "finished_at": now(),
                })
                try:
                    run = _run_record(
                        declaration, run_id=run_id, started_at=started_at,
                        input_id=input_id, artifact=artifact,
                        build_dir=build_dir, wheel=wheel, timeout=args.timeout,
                        record_dir=RESULTS / "runs" / "pending",
                    )
                except KeyboardInterrupt:
                    run = {
                        **attempt,
                        "execution": "Interrupted",
                        "interruption_kind": "interruption",
                        "finished_at": now(),
                        "outcome_reason": "runner interrupted while executing the case",
                    }
                except (OSError, RuntimeError, TypeError, ValueError) as error:
                    run = {
                        **attempt,
                        "execution": "Error",
                        "finished_at": now(),
                        "outcome_reason": f"runner failed while executing the case: {error}",
                    }
            else:
                run = {
                    **attempt,
                    "execution": "Interrupted" if _SHUTDOWN_REQUESTED else "Error",
                    **({"interruption_kind": "interruption"} if _SHUTDOWN_REQUESTED else {}),
                    "finished_at": now(),
                    "outcome_reason": (
                        "runner interrupted before the case started" if _SHUTDOWN_REQUESTED
                        else preparation_reason or "wheel preparation failed"
                    ),
                }
            run = _validated_final_run(run)
            store.publish(run)
            failed |= run["execution"] not in {"Passed", "Skipped"}
            print(f"{run['test_id']}: {run['execution']} ({run['run_id']})", flush=True)
        return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())

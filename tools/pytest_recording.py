"""Pytest declaration and comparison recorder loaded by the shared runner.

The fixture is intentionally an authoring helper rather than a second test
framework. Existing pytest assertions still determine execution outcome;
``compare`` evaluates the recorded criterion used by the visible result.
"""

from __future__ import annotations

import contextlib
import json
import os
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pytest

from tools.test_records import (
    SCHEMA_VERSION,
    RecordStore,
    artifact_identity,
    comparison_selector_matches,
    input_generation,
    locked,
    make_comparison,
    now,
    parameter_case_id,
    sha256,
    validate_run,
    verify_artifact_identity,
)


def _marker(item: pytest.Item) -> dict[str, Any] | None:
    marker = item.get_closest_marker("evidence")
    if marker is None:
        return None
    return dict(marker.kwargs)


def _case(item: pytest.Item) -> str:
    callspec = getattr(item, "callspec", None)
    if callspec is None:
        return "single"
    try:
        return parameter_case_id(callspec.params)
    except (TypeError, ValueError, OverflowError):
        # Keep collection visible and deterministic when a parametrized value
        # cannot be represented by the JSON identity contract.  The
        # declaration is marked incomplete below so it cannot enter automatic
        # numerical execution.
        return parameter_case_id({"nodeid": item.nodeid})


def _case_error(item: pytest.Item) -> str | None:
    callspec = getattr(item, "callspec", None)
    if callspec is None:
        return None
    try:
        parameter_case_id(callspec.params)
    except (TypeError, ValueError, OverflowError):
        return "parameter values are not JSON serializable"
    return None


def _comparison_keys(value: Any) -> list[str]:
    if isinstance(value, str):
        return [value]
    try:
        return [str(item) for item in value]
    except TypeError:
        return [str(value)]


def declaration(item: pytest.Item) -> dict[str, Any] | None:
    metadata = _marker(item)
    if metadata is None:
        return None
    path = Path(str(item.fspath))
    root = Path(__file__).resolve().parents[1]
    required = ("test_id", "title", "kind", "category", "lifecycle", "execution_policy", "evidence_expectation")
    missing = [name for name in required if name not in metadata]
    result: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "test_id": str(metadata.get("test_id", f"unclassified.{item.nodeid}")),
        "title": str(metadata.get("title", item.name)),
        "kind": str(metadata.get("kind", "Tooling")),
        "category": str(metadata.get("category", "Unclassified")),
        "lifecycle": str(metadata.get("lifecycle", "Implemented")),
        "execution_policy": str(metadata.get("execution_policy", "Explicit")),
        "evidence_expectation": str(metadata.get("evidence_expectation", "Outcome only")),
        "source": str(path.relative_to(root)) if path.is_relative_to(root) else str(path),
        "runner": {"kind": "pytest", "nodeid": item.nodeid},
    }
    if missing:
        result["unclassified"] = True
        result["declaration_errors"] = ["missing metadata: " + name for name in missing]
    for name in ("property_group", "requirement_id", "basis", "method", "claim", "unmet_prerequisite"):
        if name in metadata:
            result[name] = metadata[name]
    if "expected_comparisons" in metadata:
        result["expected_comparison_keys"] = _comparison_keys(metadata["expected_comparisons"])
    elif result["evidence_expectation"] == "Numerical":
        result["expected_comparison_keys"] = []
        result.setdefault("declaration_errors", []).append("numerical declaration has no expected comparison keys")
    if "expected_counts" in metadata:
        result["expected_comparison_counts"] = dict(metadata["expected_counts"])
    case_error = _case_error(item)
    if case_error:
        result["unclassified"] = True
        result.setdefault("declaration_errors", []).append(case_error)
    return result


def _unclassified_declaration(item: pytest.Item) -> dict[str, Any]:
    """Keep an undeclared item visible without inventing scientific metadata."""

    path = Path(str(item.fspath))
    root = Path(__file__).resolve().parents[1]
    return {
        "schema_version": SCHEMA_VERSION,
        "test_id": f"unclassified.pytest.{item.nodeid}",
        "title": item.nodeid,
        "kind": "Tooling",
        "category": "Unclassified",
        "lifecycle": "Implemented",
        "execution_policy": "Explicit",
        "evidence_expectation": "Outcome only",
        "source": str(path.relative_to(root)) if path.is_relative_to(root) else str(path),
        "runner": {"kind": "pytest", "nodeid": item.nodeid},
        "parameter_cases": [{"parameter_case_id": _case(item), "nodeid": item.nodeid}],
        "unclassified": True,
        "declaration_errors": ["missing evidence declaration"],
    }


@dataclass
class EvidenceRecorder:
    item: pytest.Item
    rows: list[dict[str, Any]] = field(default_factory=list)

    def _add(self, row: dict[str, Any]) -> dict[str, Any]:
        if any(existing["comparison_id"] == row["comparison_id"] for existing in self.rows):
            raise AssertionError(f"duplicate comparison id: {row['comparison_id']}")
        self.rows.append(row)
        return row

    def compare(
        self,
        comparison_id: str,
        *,
        quantity: str,
        model: float,
        reference: float,
        unit: str,
        inputs: dict[str, Any],
        criterion: dict[str, Any],
        basis: str,
        method: str,
        claim: str,
        coordinates: dict[str, Any] | None = None,
        comparison_group_id: str | None = None,
        point_id: str | None = None,
        quantity_latex: str | None = None,
        relationship_latex: str | None = None,
        inputs_latex: str | None = None,
    ) -> dict[str, Any]:
        row = make_comparison(
            comparison_id,
            quantity=quantity,
            unit=unit,
            inputs=inputs,
            model=model,
            reference=reference,
            criterion=criterion,
            basis=basis,
            method=method,
            claim=claim,
            coordinates=coordinates,
            comparison_group_id=comparison_group_id,
            point_id=point_id,
            quantity_latex=quantity_latex,
            relationship_latex=relationship_latex,
            inputs_latex=inputs_latex,
        )
        self._add(row)
        if row["judgment"] == "Criterion failed":
            raise AssertionError(f"recorded criterion failed: {comparison_id}")
        return row

    def observe(
        self,
        comparison_id: str,
        *,
        quantity: str,
        model: float,
        unit: str,
        inputs: dict[str, Any],
        reference: float | None = None,
        basis: str | None = None,
        method: str | None = None,
        claim: str | None = None,
        coordinates: dict[str, Any] | None = None,
        comparison_group_id: str | None = None,
        point_id: str | None = None,
        quantity_latex: str | None = None,
        relationship_latex: str | None = None,
        inputs_latex: str | None = None,
    ) -> dict[str, Any]:
        return self._add(make_comparison(
            comparison_id,
            quantity=quantity,
            unit=unit,
            inputs=inputs,
            model=model,
            reference=reference,
            observation=True,
            basis=basis,
            method=method,
            claim=claim,
            coordinates=coordinates,
            comparison_group_id=comparison_group_id,
            point_id=point_id,
            quantity_latex=quantity_latex,
            relationship_latex=relationship_latex,
            inputs_latex=inputs_latex,
        ))


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers",
        "evidence(...): stable test declaration and numerical comparison metadata",
    )


_COLLECTION_ERRORS: list[str] = []


def pytest_collectreport(report: pytest.CollectReport) -> None:
    if report.failed:
        _COLLECTION_ERRORS.append(report.longreprtext)


def pytest_collection_modifyitems(session: pytest.Session, config: pytest.Config, items: list[pytest.Item]) -> None:
    destination = os.environ.get("EPCSAFT_COLLECTION_RESULTS")
    if not destination:
        return
    declarations: list[dict[str, Any]] = []
    by_id: dict[str, dict[str, Any]] = {}
    errors: list[str] = list(_COLLECTION_ERRORS)
    for item in items:
        record = declaration(item)
        if record is None:
            record = _unclassified_declaration(item)
        key = record["test_id"]
        errors.extend(str(error) for error in record.get("declaration_errors", []))
        base_nodeid = str(record.get("runner", {}).get("nodeid", item.nodeid)).split("[", 1)[0]
        fingerprint = json.dumps({k: v for k, v in record.items()
                                  if k not in {"runner", "parameter_case_id", "parameter_cases"}},
                                 sort_keys=True)
        existing = by_id.get(key)
        if existing is None:
            record["runner"] = {"kind": "pytest", "nodeid": base_nodeid}
            record["parameter_cases"] = [{"parameter_case_id": _case(item), "nodeid": item.nodeid}]
            record["_metadata_fingerprint"] = fingerprint
            record["_runner_nodeid"] = base_nodeid
            by_id[key] = record
            declarations.append(record)
            continue
        if existing.get("_metadata_fingerprint") != fingerprint or existing.get("_runner_nodeid") != base_nodeid:
            errors.append(f"duplicate test_id with conflicting metadata or runner: {key}")
            continue
        case = {"parameter_case_id": _case(item), "nodeid": item.nodeid}
        if any(existing_case.get("parameter_case_id") == case["parameter_case_id"]
               for existing_case in existing["parameter_cases"]):
            errors.append(f"duplicate parameter_case_id for test_id: {key}: {case['parameter_case_id']}")
            continue
        existing["parameter_cases"].append(case)
    for record in declarations:
        record.pop("_metadata_fingerprint", None)
        record.pop("_runner_nodeid", None)
        if len(record.get("parameter_cases", [])) == 1:
            record["parameter_case_id"] = record["parameter_cases"][0]["parameter_case_id"]
    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"schema_version": SCHEMA_VERSION, "declarations": declarations, "errors": errors}, allow_nan=False) + "\n")


_RECORDER_KEY: pytest.StashKey[EvidenceRecorder] = pytest.StashKey()
_RESULTS_KEY: pytest.StashKey[dict[str, Any]] = pytest.StashKey()
_PHASES_KEY: pytest.StashKey[dict[str, pytest.TestReport]] = pytest.StashKey()
_EVIDENCE_REASON_KEY: pytest.StashKey[str | None] = pytest.StashKey()
_EVIDENCE_ROWS_KEY: pytest.StashKey[list[dict[str, Any]]] = pytest.StashKey()

_DIRECT_LOCK: Any = None
_DIRECT_STORE: RecordStore | None = None
_DIRECT_INPUT: str | None = None
_DIRECT_ARTIFACT: dict[str, Any] | None = None
_DIRECT_ATTEMPTS: dict[str, dict[str, Any]] = {}
_DIRECT_SESSION_ID: str | None = None


def _direct_mode(config: pytest.Config) -> bool:
    return not os.environ.get("EPCSAFT_PYTEST_RESULTS") and not getattr(config.option, "collectonly", False)


def _direct_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _direct_artifact(root: Path) -> dict[str, Any]:
    marker = root / "build/environment-wheel/meson-build-dir"
    stamp = root / "build/environment-wheel/preparation.json"
    wheels = sorted((root / "build/environment-wheel").glob("epcsaft-*.whl"))
    if not marker.is_file() or not stamp.is_file() or len(wheels) != 1:
        return {"state": "Unavailable", "reason": "prepared Engine wheel/build binding is unavailable"}
    try:
        binding = json.loads(stamp.read_text(encoding="utf-8"))
        if not isinstance(binding, dict) or binding.get("schema_version") != 1:
            return {"state": "Unavailable", "reason": "prepared artifact binding has an invalid schema"}
        current_input = input_generation(root)
        if binding.get("input_generation") != current_input:
            return {"state": "Unavailable", "reason": "prepared artifact was built from different inputs"}
        wheel = (root / str(binding.get("wheel", ""))).resolve()
        build_dir = (root / str(binding.get("build_dir", ""))).resolve()
        if wheel != wheels[0].resolve() or not wheel.is_relative_to(root.resolve()):
            return {"state": "Unavailable", "reason": "prepared Engine wheel binding changed"}
        if not build_dir.is_dir() or not build_dir.is_relative_to(root.resolve()):
            return {"state": "Unavailable", "reason": "prepared native build directory is unavailable"}
        if sha256(wheel) != binding.get("wheel_sha256"):
            return {"state": "Unavailable", "reason": "prepared Engine wheel changed"}
        build_dir = Path(marker.read_text(encoding="utf-8").strip())
        if build_dir.resolve() != (root / str(binding.get("build_dir", ""))).resolve():
            return {"state": "Unavailable", "reason": "prepared native build directory changed"}
        return artifact_identity(root, build_dir=build_dir, wheel=wheel, runtime=root / ".venv/bin/python")
    except (OSError, ValueError, RuntimeError, TypeError, AttributeError) as error:
        return {"state": "Unavailable", "reason": f"artifact identity unavailable: {error}"}


def _direct_attempt(item: pytest.Item) -> dict[str, Any]:
    metadata = _marker(item)
    test_id = str(metadata.get("test_id") if metadata else f"unclassified.pytest.{item.nodeid}")
    case_id = _case(item)
    run_id = f"pytest-direct-{_DIRECT_SESSION_ID or uuid.uuid4().hex}-{parameter_case_id({'nodeid': item.nodeid})}"
    started = now()
    return {
        "schema_version": SCHEMA_VERSION,
        "run_id": run_id,
        "test_id": test_id,
        "parameter_case_id": case_id,
        "input_generation": _DIRECT_INPUT,
        "artifact_identity": _DIRECT_ARTIFACT or {"state": "Unavailable"},
        "execution": "Running",
        "started_at": started,
        "finished_at": None,
        "outcome_reason": None,
        "comparisons": [],
        "artifact_digests": (_DIRECT_ARTIFACT or {}).get("artifacts", []),
        "command": ["pytest", item.nodeid],
    }


def _direct_terminal(item: pytest.Item, attempt: dict[str, Any], result: dict[str, Any] | None,
                     *, session: pytest.Session) -> dict[str, Any]:
    """Turn one direct pytest result into the shared terminal run shape."""

    if result is None:
        execution = "Interrupted" if session.shouldstop or session.testsfailed else "Error"
        reason = "pytest did not produce a result for the collected test"
        comparisons: list[dict[str, Any]] = []
    else:
        execution = str(result.get("execution", "Error"))
        reason = result.get("outcome_reason")
        comparisons = result.get("comparisons", [])
        if not isinstance(comparisons, list):
            execution = "Error"
            reason = "pytest recorder comparisons are not a list"
            comparisons = []
    current_input = input_generation(_direct_root())
    if _DIRECT_INPUT is not None and current_input != _DIRECT_INPUT:
        execution = "Error"
        reason = "computational inputs changed during execution"
    elif attempt.get("artifact_identity", {}).get("state") == "Available":
        identity_ok, identity_reason = verify_artifact_identity(
            _direct_root(), attempt["artifact_identity"]
        )
        if not identity_ok:
            execution = "Error"
            reason = "artifact identity changed during execution: " + str(identity_reason)
    if execution in {"Failed", "Skipped", "Error", "Interrupted"} and not reason:
        reason = f"pytest execution ended as {execution}"
    return {
        "schema_version": SCHEMA_VERSION,
        "run_id": attempt["run_id"],
        "test_id": attempt["test_id"],
        "parameter_case_id": attempt["parameter_case_id"],
        "input_generation": attempt["input_generation"],
        "artifact_identity": attempt["artifact_identity"],
        "execution": execution,
        "started_at": attempt["started_at"],
        "finished_at": now(),
        "outcome_reason": reason,
        "comparisons": comparisons,
        "artifact_digests": attempt.get("artifact_digests", []),
        "command": attempt.get("command", ["pytest", item.nodeid]),
    }


def _release_direct_lock() -> None:
    global _DIRECT_LOCK
    if _DIRECT_LOCK is not None:
        _DIRECT_LOCK.__exit__(None, None, None)
        _DIRECT_LOCK = None


def pytest_sessionstart(session: pytest.Session) -> None:
    """Record direct pytest invocations under the shared publication lock."""

    global _DIRECT_LOCK, _DIRECT_STORE, _DIRECT_INPUT, _DIRECT_ARTIFACT, _DIRECT_SESSION_ID
    if not _direct_mode(session.config):
        return
    root = _direct_root()
    _DIRECT_STORE = RecordStore(root)
    _DIRECT_LOCK = (
        contextlib.nullcontext()
        if os.environ.get("EPCSAFT_RECORD_LOCK_HELD") == "1"
        else locked(_DIRECT_STORE.results / ".lock")
    )
    _DIRECT_LOCK.__enter__()
    _DIRECT_STORE.reconcile_attempts()
    _DIRECT_INPUT = input_generation(root)
    _DIRECT_ARTIFACT = _direct_artifact(root)
    if _DIRECT_ARTIFACT.get("state") != "Available":
        reason = _DIRECT_ARTIFACT.get("reason", "prepared Engine artifact is unavailable")
        _release_direct_lock()
        pytest.exit(
            "direct pytest requires a current prepared Engine artifact; run "
            f"`bash .codex/environments/setup.sh` first ({reason})",
            returncode=2,
        )
    _DIRECT_SESSION_ID = uuid.uuid4().hex
    _DIRECT_ATTEMPTS.clear()


def pytest_collection_finish(session: pytest.Session) -> None:
    """Publish queued/running records before direct test processes start."""

    if _DIRECT_STORE is None or _DIRECT_INPUT is None or _DIRECT_ARTIFACT is None:
        return
    for item in session.items:
        attempt = _direct_attempt(item)
        _DIRECT_ATTEMPTS[item.nodeid] = attempt
        _DIRECT_STORE.publish(attempt)
        _DIRECT_STORE.publish({**attempt, "execution": "Running"})


@pytest.fixture
def evidence(request: pytest.FixtureRequest) -> EvidenceRecorder:
    recorder = EvidenceRecorder(request.node)
    request.node.stash[_RECORDER_KEY] = recorder
    return recorder


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item: pytest.Item, call: pytest.CallInfo[Any]):
    outcome = yield
    report = outcome.get_result()
    phases = item.stash.setdefault(_PHASES_KEY, {})
    phases[report.when] = report
    metadata = _marker(item)
    recorder = item.stash.get(_RECORDER_KEY, None)
    rows = [] if recorder is None else recorder.rows
    if report.when != "teardown":
        item.stash[_EVIDENCE_ROWS_KEY] = rows

    # Validate expected evidence after the call phase.  The final execution
    # state is written below after setup and teardown have both reported.
    evidence_reason = item.stash.get(_EVIDENCE_REASON_KEY, None)
    if report.when == "call" and metadata is not None and metadata.get("evidence_expectation", "Outcome only") == "Numerical":
        expected = _comparison_keys(metadata.get("expected_comparisons", ("*",)))
        emitted = {row["comparison_id"] for row in rows}

        def covered(key: str) -> bool:
            return any(comparison_selector_matches(value, key) for value in emitted)

        counts = dict(metadata.get("expected_counts", {}))
        count_errors = [
            f"{prefix} expected {count} rows, emitted {sum(comparison_selector_matches(value, str(prefix)) for value in emitted)}"
            for prefix, count in counts.items()
            if sum(comparison_selector_matches(value, str(prefix)) for value in emitted) != int(count)
        ]
        if (expected and not all(covered(str(key)) for key in expected)) or count_errors:
            evidence_reason = "expected numerical comparison was not emitted" + ("; " + "; ".join(count_errors) if count_errors else "")
            for key in expected:
                if not covered(str(key)):
                    rows.append(make_comparison(
                        str(key), quantity=str(key), unit="", inputs={},
                        unavailable_reason="comparison was not emitted",
                    ))
            item.stash[_EVIDENCE_REASON_KEY] = evidence_reason

    if report.when != "teardown":
        return
    rows = item.stash.get(_EVIDENCE_ROWS_KEY, rows)
    setup = phases.get("setup")
    call_report = phases.get("call")
    teardown = phases.get("teardown")
    if teardown.failed:
        execution = "Error"
        reason = teardown.longreprtext
    elif setup is not None and setup.failed:
        execution = "Error"
        reason = setup.longreprtext
    elif setup is not None and setup.skipped:
        execution = "Skipped"
        reason = setup.longreprtext
    elif call_report is not None and call_report.skipped:
        execution = "Skipped"
        reason = call_report.longreprtext
    elif call_report is not None and call_report.failed:
        execution = "Failed"
        reason = call_report.longreprtext
    elif call_report is not None and call_report.passed:
        execution = "Error" if evidence_reason else "Passed"
        reason = evidence_reason
    else:
        execution = "Error"
        reason = "test did not produce a call outcome"
    item.stash[_RESULTS_KEY] = {
        "schema_version": SCHEMA_VERSION,
        "nodeid": item.nodeid,
        "test_id": metadata.get("test_id") if metadata else f"unclassified.pytest.{item.nodeid}",
        "parameter_case_id": _case(item),
        "execution": execution,
        "outcome_reason": reason,
        "comparisons": rows,
        "started_at": now(),
        "finished_at": now(),
    }


def pytest_sessionfinish(session: pytest.Session, exitstatus: int) -> None:
    if _DIRECT_STORE is not None:
        for item in session.items:
            attempt = _DIRECT_ATTEMPTS.get(item.nodeid)
            if attempt is None:
                continue
            result = item.stash.get(_RESULTS_KEY, None)
            terminal = _direct_terminal(item, attempt, result, session=session)
            errors = validate_run(terminal)
            if errors:
                terminal["execution"] = "Error"
                terminal["outcome_reason"] = "; ".join(errors)
                terminal["comparisons"] = []
            _DIRECT_STORE.publish(terminal)
        return
    destination = os.environ.get("EPCSAFT_PYTEST_RESULTS")
    if not destination:
        return
    results: list[dict[str, Any]] = []
    for item in session.items:
        result = item.stash.get(_RESULTS_KEY, None)
        if result is not None:
            results.append(result)
    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"schema_version": SCHEMA_VERSION, "results": results}, allow_nan=False) + "\n")


def pytest_unconfigure(config: pytest.Config) -> None:
    if _DIRECT_STORE is not None:
        _release_direct_lock()

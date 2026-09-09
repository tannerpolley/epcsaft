#!/usr/bin/env python3
"""Run a bounded end-to-end integrity check for the test recording system.

The checks use temporary repositories and pytest files. Each child exercises
the recording modules as installed in this checkout; it does not run a model
or write scientific evidence. The normal command delegates each registered
check to the shared recorder; ``--check-id`` is its child protocol.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
PYTHON = ROOT / ".venv/bin/python" if (ROOT / ".venv/bin/python").is_file() else Path(sys.executable)

import tools.test_records as RECORDS
from tools.test_records import (
    ATTEMPT_ACTIVE_ENV,
    SCHEMA_VERSION,
    RecordStore,
    freshness,
    input_generation,
    sha256,
    locked,
    make_comparison,
    now,
    validate_comparison,
    validate_inventory,
    validate_run,
    verify_artifact_identity,
)


class SkipCheck(Exception):
    """A bounded environment-dependent check could not run."""


@dataclass
class Result:
    test_id: str
    name: str
    kind: str
    state: str
    detail: str = ""


@dataclass(frozen=True)
class IntegrityCheck:
    """One executable harness check and its notebook declaration."""

    test_id: str
    title: str
    function: Callable[[Path], None]
    kind: str = "real"


RESULTS: list[Result] = []
INTEGRITY_CHECKS: list[IntegrityCheck] = []


def integrity_check(test_id: str, title: str, *, kind: str = "real") -> Callable[[Callable[[Path], None]], Callable[[Path], None]]:
    """Register a check beside its executable function."""

    def decorate(function: Callable[[Path], None]) -> Callable[[Path], None]:
        INTEGRITY_CHECKS.append(IntegrityCheck(test_id, title, function, kind))
        return function

    return decorate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def exercise(check: IntegrityCheck, temp: Path) -> None:
    try:
        check.function(temp)
    except SkipCheck as error:
        RESULTS.append(Result(check.test_id, check.title, check.kind, "SKIP", str(error)))
    except AssertionError as error:
        RESULTS.append(Result(check.test_id, check.title, check.kind, "FAIL", str(error)))
    except Exception as error:  # noqa: BLE001  # keep the sequence bounded and classify runner errors
        RESULTS.append(Result(check.test_id, check.title, check.kind, "ERROR", f"{type(error).__name__}: {error}"))
    else:
        RESULTS.append(Result(check.test_id, check.title, check.kind, "PASS"))


def load_tool(name: str, filename: str) -> Any:
    path = ROOT / "tools" / filename
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


RUNNER: Any = None
INVENTORY: Any = None
COMPANION: Any = None
TOOL_LOAD_ERROR: Exception | None = None


def load_runtime_tools() -> None:
    """Load sibling tooling only when executing checks, never during discovery."""

    global RUNNER, INVENTORY, COMPANION, TOOL_LOAD_ERROR
    if RUNNER is not None or TOOL_LOAD_ERROR is not None:
        return
    try:
        RUNNER = load_tool("check_test_system_runner", "run-numerical-evidence.py")
        INVENTORY = load_tool("check_test_system_inventory", "numerical-inventory.py")
        COMPANION = load_tool("check_test_system_companion", "companion-data.py")
    except Exception as error:  # noqa: BLE001  # report a concrete backend syntax/import failure below
        RUNNER = None
        INVENTORY = None
        COMPANION = None
        TOOL_LOAD_ERROR = error


def comparison(
    comparison_id: str,
    *,
    model: float | None = 1.0,
    reference: float | None = 1.0,
    criterion: dict[str, Any] | None = None,
    observation: bool = False,
    quantity: str = "q",
    inputs: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return make_comparison(
        comparison_id,
        quantity=quantity,
        unit="1",
        inputs=inputs or {"T": 300.0},
        model=model,
        reference=reference,
        criterion=criterion,
        basis="Relationship" if criterion is not None else None,
        method="identity" if criterion is not None else None,
        claim="synthetic integrity check" if criterion is not None else None,
        observation=observation,
    )


def run_record(
    root: Path,
    test_id: str,
    *,
    execution: str = "Passed",
    comparisons: list[dict[str, Any]] | None = None,
    generation: str | None = None,
    case_id: str = "single",
    reason: str | None = None,
    interruption_kind: str | None = None,
) -> dict[str, Any]:
    if execution in {"Failed", "Skipped", "Error", "Interrupted"} and reason is None:
        reason = f"synthetic {execution.lower()}"
    record = {
        "schema_version": SCHEMA_VERSION,
        "run_id": f"run-{test_id.replace('.', '-')}-{case_id}",
        "test_id": test_id,
        "parameter_case_id": case_id,
        "input_generation": generation or input_generation(root),
        "artifact_identity": {"state": "Unavailable", "reason": "synthetic harness artifact"},
        "execution": execution,
        "started_at": now(),
        "finished_at": now(),
        "outcome_reason": reason,
        "comparisons": comparisons or [],
        "artifact_digests": [],
    }
    if interruption_kind is not None:
        record["interruption_kind"] = interruption_kind
    return record


def _pytest_environment(collection: Path | None = None, results: Path | None = None) -> dict[str, str]:
    environment = dict(os.environ)
    environment["PYTHONNOUSERSITE"] = "1"
    existing = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = str(ROOT) + (os.pathsep + existing if existing else "")
    if collection is not None:
        environment["EPCSAFT_COLLECTION_RESULTS"] = str(collection)
    else:
        environment.pop("EPCSAFT_COLLECTION_RESULTS", None)
    if results is not None:
        environment["EPCSAFT_PYTEST_RESULTS"] = str(results)
    else:
        environment.pop("EPCSAFT_PYTEST_RESULTS", None)
    return environment


def run_pytest(path: Path, temp: Path, *, collect: bool = False) -> tuple[subprocess.CompletedProcess[str], dict[str, Any]]:
    output = temp / ("collection.json" if collect else "results.json")
    command = [str(PYTHON), "-m", "pytest", "-q", "-p", "tools.pytest_recording"]
    if collect:
        command.append("--collect-only")
    command.append(str(path))
    result = subprocess.run(
        command,
        cwd=ROOT,
        env=_pytest_environment(collection=output if collect else None, results=output if not collect else None),
        capture_output=True,
        text=True,
        timeout=20,
        check=False,
    )
    require(output.is_file(), f"pytest recorder did not write {output.name}: stdout={result.stdout[-1000:]} stderr={result.stderr[-1000:]}")
    payload = json.loads(output.read_text(encoding="utf-8"))
    require(isinstance(payload, dict), "pytest recorder payload is not an object")
    return result, payload


def pytest_case_source() -> str:
    return '''
import pytest

pytestmark = pytest.mark.evidence(
    test_id="python.parameterized",
    title="Parameterized synthetic recorder",
    kind="Tooling",
    category="Bindings/Python",
    lifecycle="Implemented",
    execution_policy="Automatic fast",
    evidence_expectation="Numerical",
    expected_comparisons=("value",),
)

@pytest.mark.parametrize("x", [1, 2], ids=["one", "two"])
def test_parameterized(evidence, x):
    evidence.compare(
        "value", quantity="value", model=float(x), reference=float(x), unit="1",
        inputs={"x": x}, criterion={"type": "scalar", "rule": "reference_or_one", "atol": 0.0, "rtol": 0.0},
        basis="Relationship", method="identity", claim="parameter is preserved",
    )
'''


@integrity_check(
    "tooling.test-system.pytest-collection-parameters",
    "pytest collection/parameters/execution",
)
def check_pytest_collection_and_execution(temp: Path) -> None:
    temp.mkdir(parents=True, exist_ok=True)
    source = temp / "test_parameterized.py"
    source.write_text(pytest_case_source(), encoding="utf-8")
    collected, collection = run_pytest(source, temp, collect=True)
    require(collected.returncode == 0, f"parameter collection failed: {collected.stdout[-1000:]}")
    declarations = collection.get("declarations", [])
    require(len(declarations) == 1, f"expected one declaration, got {len(declarations)}")
    declaration = declarations[0]
    require(declaration["test_id"] == "python.parameterized", "stable Python ID was lost")
    cases = declaration.get("parameter_cases", [])
    require(len(cases) == 2, f"expected two parameter cases, got {cases}")
    require(len({case["parameter_case_id"] for case in cases}) == 2, "parameter cases are not distinct")

    executed, payload = run_pytest(source, temp, collect=False)
    require(executed.returncode == 0, f"parameter execution failed: {executed.stdout[-1000:]}")
    results = payload.get("results", [])
    require(len(results) == 2, f"expected two recorded cases, got {results}")
    require({item["execution"] for item in results} == {"Passed"}, "parameter outcome was not recorded")
    require({row["comparison_id"] for item in results for row in item["comparisons"]} == {"value"}, "comparison duplicated or lost")


@integrity_check(
    "tooling.test-system.pytest-declaration-metadata",
    "pytest missing/conflicting metadata",
)
def check_pytest_metadata(temp: Path) -> None:
    temp.mkdir(parents=True, exist_ok=True)
    missing = temp / "test_missing_declaration.py"
    missing.write_text("def test_missing():\n    assert True\n", encoding="utf-8")
    result, payload = run_pytest(missing, temp, collect=True)
    require(result.returncode == 0, "undeclared collection should still complete")
    require(payload.get("errors"), "missing declaration was silently accepted")
    require(any(item.get("unclassified") for item in payload.get("declarations", [])), "unclassified test disappeared")

    conflicting = temp / "test_conflicting_declarations.py"
    conflicting.write_text(
        '''
import pytest
marker = pytest.mark.evidence(
    test_id="python.conflict", title="conflict", kind="Tooling", category="Bindings/Python",
    lifecycle="Implemented", execution_policy="Explicit", evidence_expectation="Outcome only",
)
@marker
def test_one(): pass
@marker
def test_two(): pass
''',
        encoding="utf-8",
    )
    result, payload = run_pytest(conflicting, temp, collect=True)
    require(result.returncode == 0, "conflicting metadata collection should complete")
    require(any("duplicate test_id with conflicting" in error for error in payload.get("errors", [])), "conflicting IDs were not reported")


@integrity_check(
    "tooling.test-system.pytest-parameter-identity",
    "pytest duplicate canonical parameter cases",
)
def check_duplicate_parameter_cases(temp: Path) -> None:
    temp.mkdir(parents=True, exist_ok=True)
    path = temp / "test_duplicate_parameter_cases.py"
    path.write_text(
        '''
import pytest

pytestmark = pytest.mark.evidence(
    test_id="python.duplicate-parameter-case", title="Duplicate parameter case", kind="Tooling",
    category="Bindings/Python", lifecycle="Implemented", execution_policy="Automatic fast",
    evidence_expectation="Outcome only",
)

@pytest.mark.parametrize("params", [{"a": 1, "b": 2}, {"b": 2, "a": 1}], ids=["first", "second"])
def test_duplicate_parameter_case(params):
    assert params
''',
        encoding="utf-8",
    )
    collected, payload = run_pytest(path, temp, collect=True)
    require(collected.returncode == 0, f"duplicate parameter collection failed: {collected.stdout[-1000:]}")
    declaration = next(
        item for item in payload.get("declarations", [])
        if item.get("test_id") == "python.duplicate-parameter-case"
    )
    cases = declaration.get("parameter_cases", [])
    case_ids = [case.get("parameter_case_id") for case in cases]
    require(case_ids, f"duplicate parameter cases were not collected: {cases}")
    errors = [str(error).lower() for error in payload.get("errors", [])]
    require(
        any("parameter" in error and ("duplicate" in error or "collision" in error or "unique" in error)
            and case_ids[0].lower() in error for error in errors),
        f"duplicate canonical parameter case was not reported as an inventory error: {payload}",
    )

    duplicate_cases = temp / "test_duplicate_parameter_cases.py"
    duplicate_cases.write_text(
        '''
import pytest

@pytest.mark.evidence(
    test_id="python.duplicate-cases", title="duplicate cases", kind="Tooling", category="Bindings/Python",
    lifecycle="Implemented", execution_policy="Automatic fast", evidence_expectation="Outcome only",
)
@pytest.mark.parametrize("value", [1, 1], ids=["first", "second"])
def test_duplicate_cases(value): pass
''',
        encoding="utf-8",
    )
    result, payload = run_pytest(duplicate_cases, temp, collect=True)
    require(result.returncode == 0, "duplicate parameter collection should complete")
    require(any("duplicate parameter_case_id" in error for error in payload.get("errors", [])),
            "duplicate parameter case IDs were not reported")

    non_json = temp / "test_non_json_parameter.py"
    non_json.write_text(
        '''
import pytest

class NonJson:
    pass

@pytest.mark.evidence(
    test_id="python.non-json-parameter", title="non JSON parameter", kind="Tooling", category="Bindings/Python",
    lifecycle="Implemented", execution_policy="Automatic fast", evidence_expectation="Outcome only",
)
@pytest.mark.parametrize("value", [NonJson()])
def test_non_json_parameter(value): pass
''',
        encoding="utf-8",
    )
    result, payload = run_pytest(non_json, temp, collect=True)
    require(result.returncode == 0, "non JSON parameter collection should complete")
    require(any("not JSON serializable" in error for error in payload.get("errors", [])),
            "non JSON parameter values were not reported")


@integrity_check(
    "tooling.test-system.pytest-outcomes",
    "pytest setup/call/teardown outcomes",
)
def check_pytest_outcomes(temp: Path) -> None:
    temp.mkdir(parents=True, exist_ok=True)
    path = temp / "test_outcomes.py"
    path.write_text(
        '''
import pytest

def declared(test_id):
    return pytest.mark.evidence(
        test_id=test_id, title=test_id, kind="Tooling", category="Bindings/Python",
        lifecycle="Implemented", execution_policy="Explicit", evidence_expectation="Outcome only",
    )

@pytest.fixture
def setup_failure():
    pytest.fail("setup failure")

@pytest.fixture
def setup_skip():
    pytest.skip("setup skip")

@pytest.fixture
def teardown_failure():
    yield
    pytest.fail("teardown failure")

@declared("python.setup-failure")
def test_setup_failure(setup_failure): pass

@declared("python.setup-skip")
def test_setup_skip(setup_skip): pass

@declared("python.call-skip")
def test_call_skip(): pytest.skip("call skip")

@declared("python.call-failure")
def test_call_failure(): assert False, "call failure"

@declared("python.teardown-failure")
def test_teardown_failure(teardown_failure): pass
''',
        encoding="utf-8",
    )
    result, payload = run_pytest(path, temp, collect=False)
    require(result.returncode != 0, "failing/skip outcome fixture unexpectedly passed as a suite")
    outcomes = {item["test_id"]: item["execution"] for item in payload.get("results", [])}
    require(outcomes.get("python.setup-failure") == "Error", f"setup failure outcome: {outcomes}")
    require(outcomes.get("python.setup-skip") == "Skipped", f"setup skip outcome: {outcomes}")
    require(outcomes.get("python.call-skip") == "Skipped", f"call skip outcome: {outcomes}")
    require(outcomes.get("python.call-failure") == "Failed", f"call failure outcome: {outcomes}")
    require(outcomes.get("python.teardown-failure") == "Error", f"teardown failure outcome: {outcomes}")


@integrity_check(
    "tooling.test-system.pytest-evidence-criteria",
    "pytest expected evidence/criteria",
)
def check_pytest_expected_and_criteria(temp: Path) -> None:
    temp.mkdir(parents=True, exist_ok=True)
    path = temp / "test_expected.py"
    path.write_text(
        '''
import pytest

@pytest.mark.evidence(
    test_id="python.missing-comparison", title="missing comparison", kind="Tooling", category="Bindings/Python",
    lifecycle="Implemented", execution_policy="Explicit", evidence_expectation="Numerical",
    expected_comparisons=("missing",),
)
def test_missing_comparison(evidence):
    pass

@pytest.mark.evidence(
    test_id="python.owner-criterion", title="owner criterion", kind="Tooling", category="Bindings/Python",
    lifecycle="Implemented", execution_policy="Explicit", evidence_expectation="Numerical",
    expected_comparisons=("bad",),
)
def test_owner_criterion(evidence):
    evidence.compare(
        "bad", quantity="bad", model=2.0, reference=1.0, unit="1", inputs={"x": 1},
        criterion={"type": "scalar", "rule": "reference_scale", "atol": 0.0, "rtol": 0.0},
        basis="Relationship", method="identity", claim="deliberate criterion failure",
    )
''',
        encoding="utf-8",
    )
    result, payload = run_pytest(path, temp, collect=False)
    require(result.returncode != 0, "criterion failure was not an execution failure")
    by_id = {item["test_id"]: item for item in payload.get("results", [])}
    missing = by_id["python.missing-comparison"]
    require(missing["execution"] == "Error", f"missing expected comparison outcome: {missing}")
    missing_row = missing["comparisons"][0]
    require(missing_row["judgment"] == "Unavailable" and missing_row["evidence_state"] == "Unavailable", "missing row was presented as evidence")
    owner = by_id["python.owner-criterion"]
    require(owner["execution"] == "Failed", f"owner criterion did not fail its test: {owner}")
    require(owner["comparisons"][0]["judgment"] == "Criterion failed", "failed criterion was hidden")


@integrity_check(
    "tooling.test-system.native-parser",
    "native JSON parser/coverage/invalid records",
)
def check_native_parser(temp: Path) -> None:
    require(RUNNER is not None, f"runner module could not load: {TOOL_LOAD_ERROR}")
    criterion = {"type": "scalar", "rule": "reference_scale", "atol": 1.0e-12, "rtol": 0.0}
    records = [
        {"comparison_id": "density", "quantity": "density", "unit": "mol/m3", "model": 2.0, "reference": 2.0,
         "inputs": {"T": 300.0}, "basis": "Relationship", "method": "identity", "claim": "density identity",
         "criterion": criterion, "verdict": True, "judgment": "Criterion satisfied"},
        {"comparison_id": "ladder/0", "quantity": "dP/drho", "unit": "Pa m3/mol", "model": 1.0, "reference": 1.0,
         "inputs": {"T": 300.0}, "coordinates": {"step": 0}, "basis": "Relationship", "method": "derivative",
         "claim": "derivative ladder", "criterion": criterion, "verdict": True, "judgment": "Criterion satisfied"},
        {"comparison_id": "ladder/1", "quantity": "dP/drho", "unit": "Pa m3/mol", "model": 1.1, "reference": 1.1,
         "inputs": {"T": 300.0}, "coordinates": {"step": 1}, "basis": "Relationship", "method": "derivative",
         "claim": "derivative ladder", "criterion": criterion, "verdict": True, "judgment": "Criterion satisfied"},
        {"comparison_id": "temperature", "quantity": "temperature", "unit": "K", "model": 300.0, "reference": None,
         "inputs": {"T": 300.0}, "judgment": "Observation only"},
    ]
    output = "\n".join("EPCSAFT_COMPARE_JSON " + json.dumps(item) for item in records)
    output = "EPCSAFT_EXPECT_JSON " + json.dumps({
        "comparison_group_id": "ladder-group", "comparison_id_prefix": "ladder/", "expected_count": 2,
        "definition": {"T": 300.0}, "coordinates": {"coordinate_schema": "step"},
    }) + "\n" + output
    rows = RUNNER.parse_native_comparisons(output)
    require(len(rows) == 4, f"native quantities/ladder were not retained: {rows}")
    require({row["comparison_id"] for row in rows} == {"density", "ladder/0", "ladder/1", "temperature"}, "native IDs changed")
    observation = next(row for row in rows if row["comparison_id"] == "temperature")
    require(
        observation["judgment"] == "Observation only"
        and observation["evidence_state"] == "Complete"
        and ("reference" not in observation or observation.get("reference") is None),
        "explicit native observation did not allow a missing reference",
    )
    missing_reference = RUNNER.parse_native_comparisons(
        "EPCSAFT_COMPARE_JSON " + json.dumps({
            "comparison_id": "missing-reference",
            "quantity": "q",
            "unit": "1",
            "model": 1.0,
            "inputs": {"T": 300.0},
            "judgment": "Criterion satisfied",
            "criterion": criterion,
            "verdict": True,
            "basis": "Relationship",
            "method": "identity",
            "claim": "must have a reference",
        })
    )
    require(
        missing_reference[0]["judgment"] == "Unavailable"
        and missing_reference[0]["evidence_state"] == "Unavailable",
        "judged native comparison without a reference became evidence",
    )
    incomplete_custom = RUNNER.parse_native_comparisons(
        "EPCSAFT_COMPARE_JSON " + json.dumps({
            "comparison_id": "incomplete-custom",
            "quantity": "q",
            "unit": "1",
            "model": 1.0,
            "reference": 1.0,
            "inputs": {"T": 300.0},
            "judgment": "Criterion satisfied",
            "criterion": {"type": "custom", "verdict": True},
            "verdict": True,
            "basis": "Relationship",
            "method": "identity",
            "claim": "must have a custom rule",
        })
    )
    require(
        incomplete_custom[0]["judgment"] == "Unavailable"
        and incomplete_custom[0]["evidence_state"] == "Unavailable",
        "incomplete custom native criterion became evidence",
    )
    nonfinite_observation = RUNNER.parse_native_comparisons(
        "EPCSAFT_COMPARE_JSON " + json.dumps({
            "comparison_id": "nonfinite-observation",
            "quantity": "q",
            "unit": "1",
            "model": float("nan"),
            "inputs": {"T": 300.0},
            "judgment": "Observation only",
            "evidence_state": "Unavailable",
            "unavailable_reason": "nonfinite comparison operand",
        })
    )
    require(
        nonfinite_observation[0]["judgment"] == "Unavailable"
        and nonfinite_observation[0]["evidence_state"] == "Unavailable",
        "nonfinite native observation became complete evidence",
    )
    keys, counts, errors = RUNNER.parse_native_expectations(output)
    require(not errors and counts == {"ladder/*": 2} and keys == [], f"native expected coverage failed: {keys}, {counts}, {errors}")
    groups = RUNNER.parse_native_expectation_groups(output)
    require(
        len(groups) == 1 and groups[0]["case_definition"] == {"T": 300.0}
        and groups[0]["coordinates"] == {"coordinate_schema": "step"},
        f"native shared input definition was not retained: {groups}",
    )
    collision_rows = [
        comparison("all-pairs-active"),
        comparison("all-pairs-active-verdict"),
    ]
    _, exact_reason = RUNNER._expected_complete(
        {"evidence_expectation": "Numerical", "expected_comparison_counts": {"all-pairs-active": 1}},
        collision_rows,
    )
    require(exact_reason is None, f"exact comparison count matched a colliding prefix: {exact_reason}")
    _, wildcard_reason = RUNNER._expected_complete(
        {"evidence_expectation": "Numerical", "expected_comparison_counts": {"all-pairs-active*": 2}},
        collision_rows,
    )
    require(wildcard_reason is None, f"wildcard comparison count did not match its prefix: {wildcard_reason}")
    complete, reason = RUNNER._expected_complete(
        {"evidence_expectation": "Numerical", "expected_comparison_keys": ["density", "missing"]},
        rows[:1],
    )
    require(reason and any(row["comparison_id"] == "missing" and row["judgment"] == "Unavailable" for row in complete), "missing native comparison was not explicit")

    duplicate = RUNNER.parse_native_comparisons(
        "EPCSAFT_COMPARE_JSON " + json.dumps(records[0]) + "\n" +
        "EPCSAFT_COMPARE_JSON " + json.dumps(records[0])
    )
    require(len(duplicate) == 2, "duplicate native records were silently discarded")
    duplicate_errors = validate_run(run_record(ROOT, "native.duplicate", comparisons=duplicate))
    require(any("duplicate comparison_id" in error for error in duplicate_errors), "duplicate native records did not invalidate the run")

    malformed = RUNNER.parse_native_comparisons(
        "EPCSAFT_COMPARE_JSON {broken}\nEPCSAFT_COMPARE_JSON []\n"
        "EPCSAFT_COMPARE comparison_id=bad model=not-a-number reference=1"
    )
    require(len(malformed) == 2 and all(row["judgment"] == "Unavailable" for row in malformed), "malformed native records were accepted")
    _, _, expectation_errors = RUNNER.parse_native_expectations("EPCSAFT_EXPECT_JSON {broken}")
    require(expectation_errors, "malformed native expectations were ignored")


@integrity_check(
    "tooling.test-system.native-stream-protocol",
    "native stdout/stderr protocol separation",
    kind="simulated",
)
def check_native_stream_separation(temp: Path) -> None:
    require(RUNNER is not None, f"runner module could not load: {TOOL_LOAD_ERROR}")
    temp.mkdir(parents=True, exist_ok=True)
    stdout = temp / "native.stdout.log"
    stderr = temp / "native.stderr.log"
    protocol = "EPCSAFT_COMPARE_JSON " + json.dumps({
        "comparison_id": "stream-safe",
        "quantity": "q",
        "unit": "1",
        "model": 1.0,
        "reference": 1.0,
        "inputs": {"T": 300.0},
        "basis": "Relationship",
        "method": "identity",
        "claim": "stream separation",
        "criterion": {"type": "scalar", "rule": "reference_scale", "atol": 0.0, "rtol": 0.0},
        "verdict": True,
        "judgment": "Criterion satisfied",
    })
    source = (
        "import sys; "
        f"print({protocol!r}); "
        "print('diagnostic stderr line', file=sys.stderr)"
    )
    execution, exit_code, reason, interruption_kind = RUNNER._run_subprocess(
        [str(PYTHON), "-c", source], cwd=ROOT, timeout=5,
        environment=dict(os.environ, PYTHONNOUSERSITE="1"),
        output=stdout, error_output=stderr,
    )
    require(
        execution == "Passed" and exit_code == 0 and reason is None and interruption_kind is None,
        "stream separation child failed",
    )
    stdout_text = stdout.read_text(encoding="utf-8")
    stderr_text = stderr.read_text(encoding="utf-8")
    require("EPCSAFT_COMPARE_JSON" in stdout_text and "diagnostic stderr line" not in stdout_text,
            "stderr was spliced into native protocol output")
    require("diagnostic stderr line" in stderr_text and "EPCSAFT_COMPARE_JSON" not in stderr_text,
            "native stderr was not captured separately")
    require(len(RUNNER.parse_native_comparisons(stdout_text)) == 1,
            "separated native protocol no longer parsed")


@integrity_check(
    "tooling.test-system.criteria-edges",
    "nonfinite/zero/custom/owner criteria",
)
def check_criterion_edges(temp: Path) -> None:
    scalar = {"type": "scalar", "rule": "reference_scale", "atol": 0.0, "rtol": 0.01}
    nonfinite = comparison("nonfinite", model=float("nan"), reference=1.0, criterion=scalar)
    require(nonfinite["judgment"] == "Unavailable" and nonfinite["evidence_state"] == "Unavailable", "nonfinite model became evidence")
    zero_reference = comparison("zero-reference", model=1.0e-3, reference=0.0, criterion=scalar)
    require(zero_reference["judgment"] == "Criterion failed" and not zero_reference["verdict"], "zero reference received a false relative pass")
    missing_reference = comparison("missing-reference", model=1.0, reference=None, criterion=scalar)
    require(missing_reference["judgment"] == "Unavailable", "missing reference became an observation")
    custom = comparison(
        "custom",
        criterion={
            "type": "custom",
            "rule": "max derivative-ladder error <= 0.1",
            "norm": "max",
            "scale": 1.0,
            "verdict": True,
        },
    )
    require(custom["judgment"] == "Criterion satisfied" and not validate_comparison(custom), "custom criterion was not retained")
    contradiction = run_record(ROOT, "criterion.contradiction", comparisons=[zero_reference])
    errors = validate_run(contradiction)
    require(any("passed execution has failed comparison" in error for error in errors), "owner criterion disagreement was accepted")


@integrity_check(
    "tooling.test-system.native-executable",
    "synthetic registered native executable",
    kind="simulated",
)
def check_synthetic_native_executable(temp: Path) -> None:
    require(RUNNER is not None and INVENTORY is not None, f"backend modules could not load: {TOOL_LOAD_ERROR}")
    temp.mkdir(parents=True, exist_ok=True)
    compiler = shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")
    if compiler is None:
        raise SkipCheck("no C++ compiler available for the synthetic native producer")
    source = temp / "synthetic-native.cpp"
    executable = temp / "synthetic-native"
    source.write_text(
        r'''
#include <iostream>
#include <string>
int main(int argc, char** argv) {
  if (argc == 2 && std::string(argv[1]) == "--list") {
    std::cout << "EPCSAFT_DECL schema_version=1 test_id=native.synthetic title=synthetic kind=Numerical category=Core lifecycle=Implemented execution_policy=Automatic_fast evidence_expectation=Numerical expected_comparison_keys=ladder/* runner=synthetic\n";
    return 0;
  }
  std::cout << "EPCSAFT_EXPECT_JSON {\"comparison_id_prefix\":\"ladder/\",\"expected_count\":2}\n";
  for (int step = 0; step < 2; ++step)
    std::cout << "EPCSAFT_COMPARE_JSON {\"comparison_id\":\"ladder/" << step << "\",\"quantity\":\"q\",\"unit\":\"1\",\"model\":1,\"reference\":1,\"inputs\":{\"step\":" << step << "},\"basis\":\"Relationship\",\"method\":\"derivative\",\"claim\":\"synthetic\",\"criterion\":{\"type\":\"scalar\",\"rule\":\"reference_scale\",\"atol\":0,\"rtol\":0},\"verdict\":true,\"judgment\":\"Criterion satisfied\"}\n";
  std::cout << "EPCSAFT_CASE schema_version=1 test_id=native.synthetic execution=Passed\n";
}
''',
        encoding="utf-8",
    )
    built = subprocess.run([compiler, "-std=c++17", str(source), "-o", str(executable)], capture_output=True, text=True, timeout=10, check=False)
    require(built.returncode == 0, f"synthetic native producer did not compile: {built.stderr[-1000:]}")
    listed = subprocess.run([str(executable), "--list"], capture_output=True, text=True, timeout=5, check=False)
    require(listed.returncode == 0, "synthetic native listing failed")
    declaration = INVENTORY._parse_decl_line(listed.stdout.strip())
    require(declaration is not None and declaration["test_id"] == "native.synthetic", "native registration was not discoverable")
    require(not validate_inventory([declaration]), f"synthetic native declaration invalid: {declaration}")
    executed = subprocess.run([str(executable)], capture_output=True, text=True, timeout=5, check=False)
    require(executed.returncode == 0, "synthetic native case failed")
    rows = RUNNER.parse_native_comparisons(executed.stdout)
    require(len(rows) == 2 and {row["comparison_id"] for row in rows} == {"ladder/0", "ladder/1"}, "native executable rows were lost")
    keys, counts, errors = RUNNER.parse_native_expectations(executed.stdout)
    require(not errors and counts == {"ladder/*": 2} and not keys, "native executable coverage was not retained")
    require(
        RUNNER._parse_native_cases(executed.stdout) == {
            "native.synthetic": {"execution": "Passed", "outcome_reason": ""}
        },
        "native case outcome was lost",
    )


@integrity_check(
    "tooling.test-system.record-store",
    "RecordStore publication/corruption/replacement",
)
def check_record_store(temp: Path) -> None:
    root = temp / "records"
    store = RecordStore(root)
    row = comparison("quantity", criterion={"type": "scalar", "rule": "reference_scale", "atol": 0.0, "rtol": 0.0})
    passed = run_record(root, "case.record", comparisons=[row])
    require(not validate_run(passed), "valid synthetic pass was rejected")
    manifest = store.publish(passed)
    require(len(manifest["runs"]) == 1, "published manifest omitted pass")
    read = store.read_runs()
    require(len(read) == 1 and read[0]["execution"] == "Passed", "published pass was not readable")
    require(read[0]["comparisons"][0]["comparison_id"] == "quantity", "comparison table was not atomically published")

    # A CSV table has one column set for every row.  An unavailable/partial
    # row therefore leaves optional model/reference cells blank beside valid
    # rows; reading those blanks as strings must not discard the whole run.
    null_root = temp / "null-roundtrip"
    null_store = RecordStore(null_root)
    valid_rows = [comparison(f"valid-{index}") for index in range(985)]
    partial = comparison("partial", model=None, reference=None)
    require(not validate_comparison(partial), "synthetic partial row was invalid before publication")
    partial_run = run_record(
        null_root,
        "case.null-roundtrip",
        execution="Failed",
        comparisons=[*valid_rows, partial],
        reason="synthetic partial comparison",
    )
    require(not validate_run(partial_run), "partial comparison run was invalid before publication")
    null_store.publish(partial_run)
    roundtrip = null_store.read_runs()
    require(len(roundtrip) == 1 and roundtrip[0]["execution"] == "Failed", "partial comparison run was not readable")
    roundtrip_rows = roundtrip[0]["comparisons"]
    require(len(roundtrip_rows) == 986, "blank optional cells discarded comparison rows")
    require(sum(str(row.get("comparison_id", "")).startswith("valid-") for row in roundtrip_rows) == 985,
            "valid comparison rows were lost during null roundtrip")
    partial_rows = [row for row in roundtrip_rows if row.get("comparison_id") == "partial"]
    require(len(partial_rows) == 1 and "model" not in partial_rows[0] and "reference" not in partial_rows[0],
            "blank optional operands were reconstructed as invalid strings")
    require(not validate_run(roundtrip[0]), "partial comparison rows invalidated the readable run")

    observation_root = temp / "observation-roundtrip"
    observation_store = RecordStore(observation_root)
    observed = make_comparison(
        "observed", quantity="observed", unit="1", inputs={"x": 1}, model=1.0,
        observation=True,
    )
    require(
        observed["judgment"] == "Observation only"
        and observed["evidence_state"] == "Complete"
        and "reference" not in observed
        and not validate_comparison(observed),
        "optional-reference observation was invalid before publication",
    )
    observed_run = run_record(observation_root, "case.observation", comparisons=[observed])
    require(not validate_run(observed_run), "observation run was invalid before publication")
    observation_store.publish(observed_run)
    observed_roundtrip = observation_store.read_runs()
    require(
        len(observed_roundtrip) == 1
        and observed_roundtrip[0]["comparisons"][0]["judgment"] == "Observation only"
        and "reference" not in observed_roundtrip[0]["comparisons"][0]
        and not validate_run(observed_roundtrip[0]),
        "optional-reference observation did not survive CSV roundtrip",
    )

    coverage_root = temp / "coverage-placeholder"
    coverage_store = RecordStore(coverage_root)
    valid = comparison("valid", criterion={"type": "scalar", "rule": "reference_scale", "atol": 0.0, "rtol": 0.0})
    placeholder = make_comparison(
        "missing", quantity="missing", unit="", inputs={}, unavailable_reason="comparison was not emitted",
    )
    coverage_run = run_record(
        coverage_root,
        "case.coverage-placeholder",
        execution="Failed",
        comparisons=[valid, placeholder],
        reason="missing expected comparison",
    )
    require(not validate_run(coverage_run), "coverage placeholder run was invalid before publication")
    coverage_store.publish(coverage_run)
    coverage_roundtrip = coverage_store.read_runs()
    require(len(coverage_roundtrip) == 1 and coverage_roundtrip[0]["execution"] == "Failed",
            "coverage placeholder changed the owning failure")
    coverage_rows = {row.get("comparison_id"): row for row in coverage_roundtrip[0]["comparisons"]}
    require(set(coverage_rows) == {"valid", "missing"}, f"coverage rows were lost: {coverage_rows}")
    require(coverage_rows["valid"].get("judgment") == "Criterion satisfied", "valid row was lost beside placeholder")
    require(coverage_rows["missing"].get("judgment") == "Unavailable"
            and "model" not in coverage_rows["missing"]
            and "reference" not in coverage_rows["missing"],
            "missing coverage placeholder did not remain unavailable",
    )

    run_path = root / manifest["runs"][0]["path"]
    run_disk = json.loads(run_path.read_text(encoding="utf-8"))
    table = run_path.parent / run_disk["comparisons_artifact"]
    table.write_text("corrupt\n", encoding="utf-8")
    corrupted = store.read_runs()[0]
    require(corrupted["execution"] == "Error" and "record unavailable" in corrupted.get("outcome_reason", ""), "corrupt table was presented as a pass")

    lifecycle_root = temp / "lifecycle"
    lifecycle_store = RecordStore(lifecycle_root)
    lifecycle_store.publish(run_record(lifecycle_root, "case.crash", comparisons=[row]))
    running = run_record(lifecycle_root, "case.crash", execution="Running")
    running["finished_at"] = None
    running.pop("artifact_digests")
    lifecycle_store.publish(running)
    after_crash = lifecycle_store.read_runs()
    require(len(after_crash) == 1 and after_crash[0]["execution"] == "Running", "crashed replacement left the prior pass current")
    lifecycle_store.publish(run_record(lifecycle_root, "case.crash", execution="Failed", reason="runner crashed before publication"))
    require(lifecycle_store.read_runs()[0]["execution"] == "Failed", "failed crash replacement did not supersede running state")

    failure_root = temp / "failed-replacement"
    failure_store = RecordStore(failure_root)
    original = run_record(failure_root, "case.replace", comparisons=[row])
    failure_store.publish(original)
    failed = run_record(failure_root, "case.replace", execution="Failed", reason="replacement failure")
    original_replace = os.replace

    def fail_manifest(source: str | bytes | os.PathLike[str], destination: str | bytes | os.PathLike[str], *args: Any) -> None:
        if Path(destination) == failure_store.manifest_path:
            raise OSError("injected manifest replacement failure")
        original_replace(source, destination, *args)

    os.replace = fail_manifest  # type: ignore[assignment]
    try:
        try:
            failure_store.publish(failed)
        except OSError:
            pass
    finally:
        os.replace = original_replace  # type: ignore[assignment]
    current = failure_store.read_manifest()
    require(not any(item.get("execution") == "Passed" for item in current.get("runs", [])), "failed replacement resurrected prior pass")

    replacement_root = temp / "replacement"
    replacement_store = RecordStore(replacement_root)
    replacement_store.publish(run_record(replacement_root, "case.rerun", comparisons=[row]))
    replacement_store.publish(run_record(replacement_root, "case.rerun", execution="Failed", reason="new failure"))
    visible = replacement_store.read_runs()
    require(len(visible) == 1 and visible[0]["execution"] == "Failed", "failed rerun left a current pass")

    pointer_root = temp / "stale-pointer"
    pointer_store = RecordStore(pointer_root)
    pointer_store.publish(run_record(pointer_root, "case.pointer", comparisons=[row]))
    old_manifest = pointer_store.read_manifest()
    old_path = old_manifest["runs"][0]["path"]
    newer = run_record(pointer_root, "case.pointer", execution="Failed", reason="newer failure")
    newer["run_id"] = "run-case-pointer-new"
    pointer_store.publish(newer)
    current_manifest = pointer_store.read_manifest()
    current_manifest["runs"][0]["path"] = old_path
    pointer_store.manifest_path.write_text(json.dumps(current_manifest), encoding="utf-8")
    try:
        stale_pointer = pointer_store.read_runs()
    except Exception as error:
        raise AssertionError(f"stale manifest pointer raised instead of failing closed: {error}") from error
    require(stale_pointer and stale_pointer[0]["execution"] != "Passed", "stale manifest pointer resurrected a pass")

    malformed_root = temp / "malformed-manifest"
    malformed_store = RecordStore(malformed_root)
    malformed_store.results.mkdir(parents=True, exist_ok=True)
    malformed_store.manifest_path.write_text(
        json.dumps({"schema_version": SCHEMA_VERSION, "revision": "malformed", "runs": "corrupt"}),
        encoding="utf-8",
    )
    try:
        malformed_runs = malformed_store.read_runs()
    except Exception as error:
        raise AssertionError(f"malformed manifest raised instead of failing closed: {error}") from error
    require(isinstance(malformed_runs, list), "malformed manifest did not return a structured run list")
    require(
        not malformed_runs or all(
            isinstance(item, dict) and item.get("execution") in {"Error", "Unavailable"}
            for item in malformed_runs
        ),
        f"malformed manifest exposed non-error runs: {malformed_runs}",
    )


@integrity_check(
    "tooling.test-system.input-identity",
    "input generation identity",
)
def check_input_identity(temp: Path) -> None:
    root = temp / "identity"
    for relative in ("engine/native/properties/state.cpp", "engine/tests/test_decl.py", "data/reference.csv", "docs/companion/index.qmd"):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(relative, encoding="utf-8")
    baseline = input_generation(root)
    presentation = root / "docs/companion/index.qmd"
    presentation.write_text("presentation edit", encoding="utf-8")
    require(input_generation(root) == baseline, "presentation edit invalidated computational evidence")

    checks = (
        (root / "engine/native/properties/state.cpp", "source edit"),
        (root / "data/reference.csv", "data edit"),
        (root / "engine/tests/test_decl.py", "declaration edit"),
    )
    for path, label in checks:
        before = input_generation(root)
        path.write_text(path.read_text(encoding="utf-8") + "\nchanged", encoding="utf-8")
        require(input_generation(root) != before, f"{label} did not change input identity")
    renamed = root / "engine/native/renamed.cpp"
    original = root / "engine/native/properties/state.cpp"
    before = input_generation(root)
    original.rename(renamed)
    require(input_generation(root) != before, "source rename did not change input identity")
    before = input_generation(root)
    renamed.unlink()
    require(input_generation(root) != before, "source removal did not change input identity")


@integrity_check(
    "tooling.test-system.artifact-identity",
    "artifact identity fails closed",
)
def check_artifact_fail_closed(temp: Path) -> None:
    root = temp / "artifact"
    (root / "engine/native").mkdir(parents=True, exist_ok=True)
    (root / "engine/native/input.cpp").write_text("input", encoding="utf-8")
    invalid_payloads: tuple[Any, ...] = (
        None,
        {},
        {"state": "Unavailable", "reason": "missing"},
        {"state": "Available", "artifacts": []},
        {"state": "Available", "artifacts": "not-a-list"},
        {"state": "Available", "artifacts": [None]},
        {"state": "Available", "artifacts": ["not-an-object"]},
    )
    for payload in invalid_payloads:
        try:
            result = verify_artifact_identity(root, payload)
        except Exception as error:
            raise AssertionError(f"invalid artifact identity raised instead of failing closed: {payload!r}: {error}") from error
        require(isinstance(result, tuple) and result[0] is False, f"invalid artifact identity was accepted: {payload!r}")

    current = input_generation(root)
    outdated = run_record(root, "case.outdated", generation="old")
    status, reason = freshness(outdated, root)
    require(status == "Outdated" and reason, f"wrong-generation result stayed current: {status}, {reason}")
    mismatched = run_record(root, "case.mismatched", generation=current)
    mismatched["artifact_identity"] = {"state": "Available", "input_generation": "different"}
    status, reason = freshness(mismatched, root)
    require(status == "Unverified" and reason, f"wrong artifact generation stayed current: {status}, {reason}")

    display_root = temp / "artifact-display"
    wheel = display_root / "engine.whl"
    wheel.parent.mkdir(parents=True, exist_ok=True)
    wheel.write_bytes(b"wheel")
    identity = {
        "state": "Available", "identity": "synthetic", "input_generation": input_generation(display_root),
        "artifacts": [{"kind": "wheel", "path": "engine.whl", "sha256": sha256(wheel)}],
        "environment": {"runtime": str(display_root / "venv/bin/python")},
    }
    displayed = run_record(display_root, "case.display", generation=input_generation(display_root))
    displayed["artifact_identity"] = identity
    require(RECORDS.display_status(displayed, display_root)["freshness"] == "Current",
            "display freshness rejected an unchanged wheel")
    wheel.write_bytes(b"replaced")
    require(RECORDS.display_status(displayed, display_root)["freshness"] == "Unverified",
            "display freshness accepted a replaced wheel")
    require(freshness(displayed, display_root)[0] == "Unverified",
            "run-time verification accepted a synthetic identity without a runtime")


@integrity_check(
    "tooling.test-system.public-projection",
    "public projection filters orphan/private fields",
)
def check_public_projection(temp: Path) -> None:
    require(COMPANION is not None, f"companion module could not load: {TOOL_LOAD_ERROR}")
    plan = "### Core\n| C1 — Synthetic check | [x] | 0 | 1 | — | Verify | — |"
    inventory = {
        "schema_version": SCHEMA_VERSION,
        "state": "Complete",
        "declarations": [{
            "schema_version": SCHEMA_VERSION,
            "test_id": "current.test",
            "title": "Current test",
            "kind": "Tooling",
            "category": "Core",
            "lifecycle": "Implemented",
            "execution_policy": "Explicit",
            "evidence_expectation": "Outcome only",
        }, {
            "schema_version": SCHEMA_VERSION,
            "test_id": "unclassified.test",
            "title": "engine/tests/test_example.py::test_example",
            "kind": "Tooling",
            "category": "Unclassified",
            "lifecycle": "Implemented",
            "execution_policy": "Explicit",
            "evidence_expectation": "Outcome only",
            "unclassified": True,
        }],
    }
    current = run_record(ROOT, "current.test")
    current["freshness"] = "Current"
    current["comparisons"] = [comparison(
        "custom", criterion={
            "type": "custom", "rule": "residual <= threshold", "verdict": True,
            "helper": "engine.tests.test_example._criterion", "rationale": "synthetic",
        },
    )]
    orphan = run_record(ROOT, "removed.test")
    orphan["freshness"] = "Current"
    projection = COMPANION.companion_projection(
        plan,
        inventory,
        [current, orphan],
        layers=("Core",),
    )
    require(
        {run.get("test_id") for run in projection["runs"]} == {"current.test"},
        f"orphan run leaked into public projection: {projection['runs']}",
    )
    public_run = projection["runs"][0]
    require("input_generation" not in public_run, "private input generation leaked into public run")
    public_comparison = public_run["comparisons"][0]
    require("helper" not in public_comparison.get("criterion", {}), "criterion helper leaked into public projection")
    public_declarations = {item["test_id"]: item for item in projection["declarations"]}
    require(public_declarations["unclassified.test"]["title"] == "Unclassified check",
            "unclassified nodeid title leaked into public projection")
    safe_error = COMPANION._safe_run({
        "test_id": "current.test", "execution": "Error",
        "outcome_reason": "record unavailable: /secret/path/traceback",
        "freshness": "Unverified", "freshness_reason": "artifact changed: /secret/lib.so",
    })
    require(
        safe_error["outcome_reason"] == "test record unavailable"
        and safe_error["freshness_reason"] == "artifact identity changed"
        and "/secret" not in json.dumps(safe_error),
        "public run projection leaked raw or collapsed artifact reasons",
    )
    safe_comparison = COMPANION._safe_comparison({
        "schema_version": SCHEMA_VERSION, "comparison_id": "missing", "quantity": "q", "unit": "1",
        "inputs": {}, "judgment": "Unavailable", "evidence_state": "Unavailable",
        "unavailable_reason": "reference operand is missing or nonfinite",
    })
    require(
        safe_comparison["unavailable_reason"] == "reference operand is missing or nonfinite",
        "public comparison projection erased the controlled nonfinite reason",
    )


@integrity_check(
    "tooling.test-system.input-groups",
    "shared input groups validate and roundtrip",
)
def check_input_groups(temp: Path) -> None:
    require(RUNNER is not None and COMPANION is not None, f"input-group modules could not load: {TOOL_LOAD_ERROR}")
    rows = [
        comparison("ladder/0", inputs={"T": 300.0, "P": 1.0}),
        comparison("ladder/1", inputs={"T": None, "P": 1.0}),
    ]
    rows, groups = RUNNER._factor_input_groups(rows, [{
        "group_ref": "g0",
        "comparison_group_id": "ladder",
        "comparison_id_prefix": "ladder/",
        "case_definition": {"P": 1.0},
        "shared_inputs": {},
    }])
    require(groups[0]["shared_inputs"] == {"P": 1.0}, f"shared input was not factored: {groups}")
    require(all(row.get("group_ref") == "g0" for row in rows), f"row group references were lost: {rows}")
    require(rows[1]["inputs"].get("T") is None, "null row override was stripped")
    root = temp / "input-groups"
    store = RecordStore(root)
    run = run_record(root, "case.input-groups", comparisons=rows)
    groups[0]["shared_inputs"]["_private"] = "/secret/path"
    run["input_groups"] = groups
    require(not validate_run(run), f"valid input groups were rejected: {validate_run(run)}")
    store.publish(run)
    roundtrip = store.read_runs()
    require(len(roundtrip) == 1, "input-group run was not published")
    require(roundtrip[0]["input_groups"][0]["shared_inputs"].get("P") == 1.0
            and roundtrip[0]["input_groups"][0]["shared_inputs"].get("_private") == "/secret/path",
            "shared inputs changed during roundtrip")
    roundtrip_row = roundtrip[0]["comparisons"][1]
    require(roundtrip_row.get("group_ref") == "g0" and roundtrip_row["inputs"].get("T") is None,
            "row group reference or null override changed during CSV roundtrip")
    public_run = COMPANION._safe_run(roundtrip[0])
    require(public_run["input_groups"][0]["shared_inputs"] == {"P": 1.0},
            "public input-group projection leaked private fields or dropped shared input")
    require(public_run["comparisons"][1]["group_ref"] == "g0", "public row group reference was dropped")

    duplicate = dict(run, input_groups=[groups[0], dict(groups[0])])
    require(any("duplicate input_group.group_ref" in error for error in validate_run(duplicate)),
            "duplicate input-group references were accepted")
    unresolved = dict(run, comparisons=[dict(rows[0], group_ref="missing")])
    require(any("comparison.group_ref" in error for error in validate_run(unresolved)),
            "unresolved row group reference was accepted")
    missing = dict(run)
    missing.pop("input_groups")
    require(any("comparison.group_ref" in error for error in validate_run(missing)),
            "missing input-group definition was accepted")

    manifest = store.read_manifest()
    record_path = root / manifest["runs"][0]["path"]
    record = json.loads(record_path.read_text(encoding="utf-8"))
    record["input_groups"][0]["shared_inputs"]["P"] = 2.0
    record_path.write_text(json.dumps(record), encoding="utf-8")
    invalidated = store.read_runs()
    require(invalidated[0]["execution"] == "Error" and "input groups" in invalidated[0].get("outcome_reason", ""),
            "edited shared inputs bypassed the manifest integrity digest")


def lock_child_code() -> str:
    return '''
import sys, time
from pathlib import Path
sys.path.insert(0, sys.argv[3])
from tools.test_records import locked
with locked(Path(sys.argv[1])):
    with Path(sys.argv[2]).open("a", encoding="utf-8") as stream:
        stream.write(sys.argv[4] + "-start\\n"); stream.flush()
        time.sleep(float(sys.argv[5]))
        stream.write(sys.argv[4] + "-end\\n"); stream.flush()
'''


@integrity_check(
    "tooling.test-system.shared-lock",
    "shared lock concurrence/nonblocking",
)
def check_shared_lock(temp: Path) -> None:
    lock = temp / "shared.lock"
    log = temp / "lock.log"
    command = [sys.executable, "-c", lock_child_code(), str(lock), str(log), str(ROOT), "one", "0.25"]
    first = subprocess.Popen(command)
    time.sleep(0.03)
    second = subprocess.Popen(command[:-2] + ["two", "0.05"])
    require(first.wait(timeout=5) == 0 and second.wait(timeout=5) == 0, "lock child failed")
    lines = log.read_text(encoding="utf-8").splitlines()
    require(
        len(lines) == 4 and lines[0].endswith("-start") and lines[1].endswith("-end")
        and lines[2].endswith("-start") and lines[3].endswith("-end"),
        f"shared lock allowed overlap: {lines}",
    )

    holder_ready = temp / "holder.ready"
    holder_code = '''
import sys, time
from pathlib import Path
sys.path.insert(0, sys.argv[3])
from tools.test_records import locked
with locked(Path(sys.argv[1])):
    Path(sys.argv[2]).write_text("ready", encoding="utf-8")
    time.sleep(0.4)
'''
    holder = subprocess.Popen([sys.executable, "-c", holder_code, str(lock), str(holder_ready), str(ROOT)])
    deadline = time.monotonic() + 2.0
    while not holder_ready.exists() and time.monotonic() < deadline:
        time.sleep(0.01)
    require(holder_ready.exists(), "lock holder did not acquire lock")
    blocked = False
    try:
        with locked(lock, blocking=False):
            pass
    except BlockingIOError:
        blocked = True
    finally:
        require(holder.wait(timeout=5) == 0, "lock holder failed")
    require(blocked, "nonblocking lock acquisition ignored a held shared lock")


@integrity_check(
    "tooling.test-system.coalesced-failure-gate",
    "coalesced failure/error/interruption gate",
    kind="simulated",
)
def check_coalesced_failure_gate(temp: Path) -> None:
    require(RUNNER is not None, f"runner module could not load: {TOOL_LOAD_ERROR}")
    malformed_root = temp / "malformed-interrupted"
    malformed = run_record(
        malformed_root, "automatic.malformed-interrupted", execution="Interrupted",
        generation="same", reason="synthetic interrupted row", interruption_kind="invalid",
    )
    finalized = RUNNER._validated_final_run(malformed)
    require(
        finalized["execution"] == "Error"
        and "interruption_kind" not in finalized
        and not validate_run(finalized),
        f"malformed Interrupted row was not finalized as publishable Error: {finalized}",
    )
    RecordStore(malformed_root).publish(finalized)

    original_freshness = RUNNER.freshness
    try:
        for execution in ("Failed", "Error"):
            root = temp / execution.lower()
            store = RecordStore(root)
            test_id = f"automatic.{execution.lower()}"
            store.publish(run_record(root, test_id, execution=execution, generation="same", reason=f"prior {execution.lower()}"))
            declaration = {"test_id": test_id, "_parameter_case_id": "single"}
            selected, coalesced_failed = RUNNER._coalesced_cases(
                store, [declaration], "same", automatic=True,
            )
            require(selected and not coalesced_failed, f"unavailable artifact was reused for prior {execution}")

            # A valid, currently applicable artifact is the only case where an
            # automatic terminal failure may be coalesced into a nonzero gate.
            RUNNER.freshness = lambda _run, _root: ("Current", None)
            selected, coalesced_failed = RUNNER._coalesced_cases(
                store, [declaration], "same", automatic=True,
            )
            require(not selected and coalesced_failed, f"prior {execution} was not retained as a failed coalesced gate")
            require(int(coalesced_failed) != 0, f"coalesced prior {execution} would return gate 0")
            RUNNER.freshness = original_freshness

        interrupted_root = temp / "interrupted"
        interrupted_store = RecordStore(interrupted_root)
        interrupted_id = "automatic.interrupted"
        timeout_stdout = interrupted_root / "timeout.stdout"
        timeout_stderr = interrupted_root / "timeout.stderr"
        timeout_execution, timeout_exit, timeout_reason, timeout_kind = RUNNER._run_subprocess(
            [str(PYTHON), "-c", "import time; time.sleep(2)"], cwd=ROOT, timeout=1,
            environment=dict(os.environ, PYTHONNOUSERSITE="1"),
            output=timeout_stdout, error_output=timeout_stderr,
        )
        require(
            timeout_execution == "Interrupted" and timeout_exit is None
            and timeout_reason == "process timed out" and timeout_kind == "timeout",
            f"timeout interruption was not classified explicitly: {timeout_execution}, {timeout_exit}, {timeout_reason}, {timeout_kind}",
        )
        interrupted_store.publish(run_record(
            interrupted_root, interrupted_id, execution="Interrupted", generation="same",
            reason="runner interrupted by systemd restart",
            interruption_kind="interruption",
        ))
        interrupted_declaration = {"test_id": interrupted_id, "_parameter_case_id": "single"}
        RUNNER.freshness = lambda _run, _root: ("Current", None)
        selected, coalesced_failed = RUNNER._coalesced_cases(
            interrupted_store, [interrupted_declaration], "same", automatic=True,
        )
        require(selected and not coalesced_failed, "infrastructure interruption was incorrectly coalesced")
        interrupted_store.publish(run_record(
            interrupted_root, interrupted_id, execution="Interrupted", generation="same",
            reason="process timed out",
            interruption_kind="timeout",
        ))
        selected, coalesced_failed = RUNNER._coalesced_cases(
            interrupted_store, [interrupted_declaration], "same", automatic=True,
        )
        require(not selected and coalesced_failed, "timed-out interruption was retried instead of gated")
    finally:
        RUNNER.freshness = original_freshness

    corrupt_root = temp / "corrupt-coalesced-record"
    corrupt_store = RecordStore(corrupt_root)
    corrupt_store.publish(run_record(corrupt_root, "automatic.corrupt", generation="same"))
    corrupt_path = corrupt_root / corrupt_store.read_manifest()["runs"][0]["path"]
    corrupt_path.write_text("{corrupt", encoding="utf-8")
    selected, coalesced_failed = RUNNER._coalesced_cases(
        corrupt_store, [{"test_id": "automatic.corrupt", "_parameter_case_id": "single"}],
        "same", automatic=True,
    )
    require(selected and not coalesced_failed, "corrupt run record was coalesced from manifest metadata")

    # Exercise the runner's finalization path with real child processes.  A
    # timeout or shutdown remains Interrupted even when a native case is left
    # Running or a pytest/integrity reporter is absent or claims Passed.
    finalization_root = temp / "finalization"
    (finalization_root / "build").mkdir(parents=True, exist_ok=True)
    (finalization_root / "engine/tests").mkdir(parents=True, exist_ok=True)
    (finalization_root / ".venv/bin").mkdir(parents=True, exist_ok=True)
    final_python = ROOT / ".venv/bin/python"
    final_python_wrapper = finalization_root / ".venv/bin/python"
    final_python_wrapper.write_text(f"#!/bin/sh\nexec {final_python} \"$@\"\n", encoding="utf-8")
    final_python_wrapper.chmod(0o755)
    native_target = finalization_root / "build/final-native"
    native_target.write_text(
        f"#!{sys.executable}\n"
        "import time\n"
        "print('EPCSAFT_CASE test_id=native.final execution=Running', flush=True)\n"
        "time.sleep(2)\n",
        encoding="utf-8",
    )
    native_target.chmod(0o755)
    pytest_source = finalization_root / "engine/tests/test_finalization.py"
    pytest_source.write_text(
        """
import time
import pytest

pytestmark = pytest.mark.evidence(
    test_id="pytest.final", title="pytest finalization", kind="Tooling",
    category="Tooling", lifecycle="Implemented", execution_policy="Automatic fast",
    evidence_expectation="Outcome only",
)

def test_finalization():
    time.sleep(2)
""".strip()
        + "\n",
        encoding="utf-8",
    )
    integrity_timeout = finalization_root / "tools/integrity-timeout.py"
    integrity_timeout.parent.mkdir(parents=True, exist_ok=True)
    integrity_timeout.write_text("import time\ntime.sleep(2)\n", encoding="utf-8")
    integrity_pass = finalization_root / "tools/integrity-pass.py"
    integrity_pass.write_text(
        "import json\n"
        "print('EPCSAFT_INTEGRITY_JSON ' + json.dumps({\"execution\": \"Passed\", "
        "\"test_id\": \"integrity.shutdown\", \"title\": \"shutdown\", \"scope\": \"simulated\"}), flush=True)\n",
        encoding="utf-8",
    )
    final_records = finalization_root / "records"
    final_records.mkdir(parents=True, exist_ok=True)
    original_runner_root = RUNNER.ROOT
    original_input_generation = RUNNER.input_generation
    original_verify_artifact_identity = RUNNER.verify_artifact_identity
    original_shutdown = RUNNER._SHUTDOWN_REQUESTED
    original_pythonpath = os.environ.get("PYTHONPATH")
    RUNNER.ROOT = finalization_root
    RUNNER.input_generation = lambda _root: "same"
    RUNNER.verify_artifact_identity = lambda *_args: (True, None)
    os.environ["PYTHONPATH"] = str(ROOT)
    try:
        native_run = RUNNER._run_record(
            {
                "test_id": "native.final", "_parameter_case_id": "single",
                "runner": {"kind": "native", "target": "final-native", "arguments": []},
                "evidence_expectation": "Outcome only",
            },
            input_id="same", artifact={"state": "Unavailable"},
            build_dir=finalization_root / "build", wheel=None, timeout=1,
            record_dir=final_records,
        )
        require(
            native_run["execution"] == "Interrupted"
            and native_run.get("interruption_kind") == "timeout"
            and "process timed out" in native_run["outcome_reason"]
            and "native case did not emit a terminal outcome" in native_run["outcome_reason"],
            f"native interruption was finalized incorrectly: {native_run}",
        )

        pytest_run = RUNNER._run_record(
            {
                "test_id": "pytest.final", "_parameter_case_id": "single",
                "runner": {"kind": "pytest", "nodeid": "tests/test_finalization.py::test_finalization"},
                "evidence_expectation": "Outcome only",
            },
            input_id="same", artifact={"state": "Unavailable"},
            build_dir=None, wheel=None, timeout=1,
            record_dir=final_records,
        )
        require(
            pytest_run["execution"] == "Interrupted"
            and pytest_run.get("interruption_kind") == "timeout"
            and "process timed out" in pytest_run["outcome_reason"]
            and "pytest reporter did not write results" in pytest_run["outcome_reason"],
            f"pytest interruption was finalized incorrectly: {pytest_run}",
        )

        integrity_run = RUNNER._run_record(
            {
                "test_id": "integrity.timeout", "_parameter_case_id": "single",
                "runner": {"kind": "integrity", "script": "tools/integrity-timeout.py",
                           "check_id": "integrity.timeout"},
                "evidence_expectation": "Outcome only",
            },
            input_id="same", artifact={"state": "Unavailable"},
            build_dir=None, wheel=None, timeout=1,
            record_dir=final_records,
        )
        require(
            integrity_run["execution"] == "Interrupted"
            and integrity_run.get("interruption_kind") == "timeout"
            and "process timed out" in integrity_run["outcome_reason"]
            and "integrity check did not emit exactly one matching outcome" in integrity_run["outcome_reason"],
            f"integrity interruption was finalized incorrectly: {integrity_run}",
        )

        RUNNER._SHUTDOWN_REQUESTED = True
        shutdown_run = RUNNER._run_record(
            {
                "test_id": "integrity.shutdown", "_parameter_case_id": "single",
                "runner": {"kind": "integrity", "script": "tools/integrity-pass.py",
                           "check_id": "integrity.shutdown"},
                "evidence_expectation": "Outcome only",
            },
            input_id="same", artifact={"state": "Unavailable"},
            build_dir=None, wheel=None, timeout=1,
            record_dir=final_records,
        )
        require(
            shutdown_run["execution"] == "Interrupted"
            and shutdown_run.get("interruption_kind") == "interruption"
            and "runner interrupted by shutdown" in shutdown_run["outcome_reason"],
            f"shutdown reporter output changed interruption to a terminal result: {shutdown_run}",
        )

        final_store = RecordStore(finalization_root)
        RUNNER.freshness = lambda _run, _root: ("Current", None)
        final_store.publish(native_run)
        selected, coalesced_failed = RUNNER._coalesced_cases(
            final_store, [{"test_id": "native.final", "_parameter_case_id": "single"}],
            "same", automatic=True,
        )
        require(not selected and coalesced_failed, "timed-out native finalization was retried")
        final_store.publish(shutdown_run)
        selected, coalesced_failed = RUNNER._coalesced_cases(
            final_store, [{"test_id": "integrity.shutdown", "_parameter_case_id": "single"}],
            "same", automatic=True,
        )
        require(selected and not coalesced_failed, "shutdown interruption was coalesced instead of retried")
    finally:
        RUNNER.ROOT = original_runner_root
        RUNNER.input_generation = original_input_generation
        RUNNER.verify_artifact_identity = original_verify_artifact_identity
        RUNNER._SHUTDOWN_REQUESTED = original_shutdown
        if original_pythonpath is None:
            os.environ.pop("PYTHONPATH", None)
        else:
            os.environ["PYTHONPATH"] = original_pythonpath

    # A nested record-command inherits the lock and attempt marker from its
    # parent.  It must not reconcile the parent's live Running row to
    # Interrupted; a later automatic selection must still retry that
    # interrupted generation.
    nested_root = temp / "nested-record-command"
    nested_root.mkdir(parents=True, exist_ok=True)
    nested_observation = nested_root / "nested-observation.json"
    nested_script = nested_root / "nested.py"
    nested_script.write_text(
        """
import importlib.util
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
observation = Path(sys.argv[2])
record_command = Path(sys.argv[3])
spec = importlib.util.spec_from_file_location("nested_record_command", record_command)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
module.ROOT = root
module.RESULTS = root / "build" / "numerical-evidence"
sys.argv = [str(record_command), "--test-id", "tooling.nested.child", "--", sys.executable,
            "-c", (
                "import json, sys; "
                "from pathlib import Path; "
                "sys.path.insert(0, sys.argv[2]); "
                "from tools.test_records import RecordStore; "
                "store = RecordStore(Path(sys.argv[1])); "
                "runs = store.read_runs(); "
                "Path(sys.argv[3]).write_text(json.dumps(runs), encoding='utf-8')"
            ), str(root), str(root), str(observation)]
raise SystemExit(module.main())
""".strip()
        + "\n",
        encoding="utf-8",
    )
    parent_command = [
        str(PYTHON), str(nested_script), str(nested_root), str(nested_observation),
        str(ROOT / "tools/record-command.py"),
    ]
    parent_environment = dict(os.environ)
    parent_environment.pop(ATTEMPT_ACTIVE_ENV, None)
    parent_environment["PYTHONPATH"] = str(ROOT)
    parent_environment["EPCSAFT_RECORD_LOCK_HELD"] = "0"
    completed = subprocess.run(
        [
            str(PYTHON), "-c",
            (
                "import importlib.util, sys; "
                "spec=importlib.util.spec_from_file_location('parent_record_command', sys.argv[1]); "
                "module=importlib.util.module_from_spec(spec); sys.modules[spec.name]=module; "
                "spec.loader.exec_module(module); module.ROOT=__import__('pathlib').Path(sys.argv[2]); "
                "module.RESULTS=module.ROOT/'build'/'numerical-evidence'; "
                "sys.argv=['record-command.py','--test-id','tooling.nested.parent','--'] + sys.argv[3:]; "
                "raise SystemExit(module.main())"
            ),
            str(ROOT / "tools/record-command.py"), str(nested_root), *parent_command,
        ],
        cwd=ROOT,
        env=parent_environment,
        capture_output=True,
        text=True,
        timeout=20,
        check=False,
    )
    require(completed.returncode == 0, f"nested record-command failed: {completed.stdout[-1000:]} {completed.stderr[-1000:]}")
    observed = json.loads(nested_observation.read_text(encoding="utf-8"))
    parent_running = [run for run in observed if run.get("test_id") == "tooling.nested.parent"]
    require(parent_running and parent_running[-1]["execution"] == "Running",
            f"nested recorder reconciled live parent: {observed}")
    retry_root = nested_root / "retry"
    retry_store = RecordStore(retry_root)
    retry_id = "automatic.nested-interrupted"
    retry_store.publish(run_record(
        retry_root, retry_id, execution="Interrupted", generation="same",
        reason="previous runner stopped before publishing a terminal outcome",
        interruption_kind="abandoned",
    ))
    selected, coalesced_failed = RUNNER._coalesced_cases(
        retry_store, [{"test_id": retry_id, "_parameter_case_id": "single"}], "same", automatic=True,
    )
    require(selected and not coalesced_failed, "reconciled interrupted generation was not retried")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-id", help="run one declared integrity check")
    args = parser.parse_args()
    if not args.check_id:
        # The normal command is the supported recording entry point.  A
        # ``--check-id`` child remains a protocol-only producer so the shared
        # runner can capture exactly one outcome without recursing here.
        command = [
            str(PYTHON), str(ROOT / "tools/run-numerical-evidence.py"),
            *[check.test_id for check in INTEGRITY_CHECKS],
        ]
        os.execv(command[0], command)
        raise AssertionError("shared recording runner returned unexpectedly")
    selected = INTEGRITY_CHECKS
    if args.check_id:
        selected = [check for check in INTEGRITY_CHECKS if check.test_id == args.check_id]
        if not selected:
            parser.error(f"unknown integrity check: {args.check_id}")

    load_runtime_tools()
    RESULTS.clear()
    (ROOT / "build").mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="check-test-system-", dir=ROOT / "build") as directory:
        temp = Path(directory)
        for check in selected:
            exercise(check, temp / check.test_id.replace(".", "-"))

    failures = [result for result in RESULTS if result.state in {"FAIL", "ERROR"}]
    skips = [result for result in RESULTS if result.state == "SKIP"]
    executions = {"PASS": "Passed", "FAIL": "Failed", "ERROR": "Error", "SKIP": "Skipped"}
    for result in RESULTS:
        suffix = f" — {result.detail}" if result.detail else ""
        print(f"[{result.state}] {result.kind}: {result.name}{suffix}")
        print("EPCSAFT_INTEGRITY_JSON " + json.dumps({
            "test_id": result.test_id,
            "title": result.name,
            "scope": result.kind,
            "execution": executions[result.state],
            "outcome_reason": result.detail or None,
        }, sort_keys=True))
    print(f"test-system harness: {len(RESULTS) - len(failures) - len(skips)} passed, {len(failures)} failed, {len(skips)} skipped")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

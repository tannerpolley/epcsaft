"""Project shared declarations and records into the live companion snapshot.

The projection deliberately has no knowledge of F1/F5, source-specific fields,
or property names. Producers own quantity labels, equations, inputs and
criteria; this module only validates the record boundary and joins runs to
their declarations.
"""

from __future__ import annotations

import re
import runpy
import sys
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.test_records import SCHEMA_VERSION, digest

_ROADMAP = runpy.run_path(__file__.replace("companion-data.py", "check-roadmap.py"))
LAYERS = tuple(_ROADMAP["LAYERS"])
_PLAN_STATUS = {"[ ]": "Planned", "[-]": "Active", "[!]": "Blocked", "[x]": "Qualified"}


def _plain(text: str) -> str:
    text = re.sub(r"\[([^]]+)\]\([^)]*\)", r"\1", text)
    text = re.sub(r"[`*_]+", "", text)
    return " ".join(text.split())


def preparation_projection(source: str) -> list[dict[str, str]]:
    marker = "### Test-system preparation"
    if marker not in source:
        return []
    section = source.split(marker, 1)[1].split("\n### ", 1)[0]
    rows: list[dict[str, str]] = []
    for line in section.splitlines():
        if not line.startswith("| "):
            continue
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        if cells[0] in {"Work", "---"}:
            continue
        if len(cells) != 3:
            raise ValueError("invalid test-system preparation row")
        rows.append({"name": _plain(cells[0]), "status": _plain(cells[1]), "required": _plain(cells[2])})
    return rows


def plan_projection(source: str, layers: tuple[str, ...] = LAYERS) -> list[dict[str, Any]]:
    """Group the plan's milestone rows by layer, carrying prerequisite and evidence links."""

    groups: dict[str, list[dict[str, Any]]] = {layer: [] for layer in layers}
    for row in _ROADMAP["rows"](source):
        if row["layer"] not in groups:
            continue
        groups[row["layer"]].append({
            "id": row["id"],
            "name": _plain(row["title"].split("—", 1)[1]) if "—" in row["title"] else row["id"],
            "status": _PLAN_STATUS.get(row["status"], "Planned"),
            "remaining": _plain(row["remaining"]),
            "actual": _plain(row["actual"]),
            "prerequisites": row["dependencies"],
            "required": _plain(row["required"]),
            "evidence": row["links"],
        })
    return [{"name": name, "rows": groups[name]} for name in layers]


def _public_comparison_reason(value: Any) -> str:
    """Keep controlled evidence reasons while hiding producer paths/details."""

    text = str(value or "").strip()
    if text in {
        "model operand is missing or nonfinite",
        "reference operand is missing or nonfinite",
        "comparison was not emitted",
    }:
        return text
    if text.startswith("expected ") and " comparison rows, emitted " in text:
        return "expected comparison rows were not emitted"
    if text.startswith("native criterion omitted required metadata"):
        return "criterion metadata unavailable"
    return "comparison evidence unavailable"


def _public_outcome_reason(value: Any, execution: str) -> str:
    """Expose concise controlled outcomes without raw exceptions or paths."""

    text = str(value or "").strip()
    exact = {
        "computational inputs changed during execution": "computational inputs changed",
        "native executable unavailable": "native executable unavailable",
        "study has no executable test entry point": "study has no executable test entry point",
        "declaration has no executable runner": "declaration has no executable runner",
        "runner interrupted while executing the case": "test runner interrupted",
        "test did not produce a call outcome": "test produced no call outcome",
        "previous runner stopped before publishing a terminal outcome": "previous run was interrupted",
        "command timed out": "test command timed out",
        "command interrupted": "test command interrupted",
    }
    if text in exact:
        return exact[text]
    if text.startswith("artifact identity changed during execution"):
        return "artifact identity changed during execution"
    if text.startswith(("artifact identity unavailable", "prepared artifact binding unavailable")):
        return "artifact identity unavailable"
    if text.startswith("inputs changed during wheel preparation"):
        return "computational inputs changed"
    if text.startswith("wheel preparation failed"):
        return "wheel preparation failed"
    if text.startswith("runner failed while executing the case"):
        return "test runner failed"
    if text.startswith("record unavailable:"):
        return "test record unavailable"
    if text.startswith("process exit="):
        return "test process exited with failure"
    if text.startswith("expected numerical comparison was not emitted"):
        return "expected numerical comparison was not emitted"
    return {
        "Failed": "test execution failed",
        "Skipped": "test was skipped",
        "Error": "test record is unavailable",
        "Interrupted": "test execution was interrupted",
    }.get(execution, "test outcome unavailable")


def _public_freshness_reason(value: Any, freshness: str) -> str:
    """Project freshness causes as a small stable public vocabulary."""

    text = str(value or "").strip()
    if freshness == "Outdated" or text == "computational inputs changed":
        return "computational inputs changed"
    if any(token in text.lower() for token in ("unavailable", "missing", "prepared")):
        return "artifact identity unavailable"
    if text:
        return "artifact identity changed"
    return "artifact identity unavailable"


def _safe_run(run: dict[str, Any]) -> dict[str, Any]:
    """Strip machinery-only fields before handing a run to browser code."""

    visible = {
        "schema_version", "run_id", "test_id", "parameter_case_id", "execution",
        "started_at", "finished_at", "comparisons", "freshness", "input_groups",
    }
    result = {key: value for key, value in run.items() if key in visible}
    if isinstance(result.get("input_groups"), list):
        result["input_groups"] = [
            _public_value(item) for item in result["input_groups"] if isinstance(item, dict)
        ]
    else:
        result["input_groups"] = []
    result["comparisons"] = [_safe_comparison(item) for item in run.get("comparisons", [])
                              if isinstance(item, dict)]
    execution = str(run.get("execution", ""))
    freshness = str(run.get("freshness", ""))
    if execution in {"Failed", "Skipped", "Error", "Interrupted"}:
        result["outcome_reason"] = _public_outcome_reason(run.get("outcome_reason"), execution)
    if freshness in {"Outdated", "Unverified"}:
        result["freshness_reason"] = _public_freshness_reason(run.get("freshness_reason"), freshness)
    return result


def _safe_error_summary(errors: Any) -> list[str]:
    """Keep error state useful without exposing paths, commands or traces."""

    if not isinstance(errors, list):
        return ["inventory issue"] if errors else []
    labels: list[str] = []
    for error in errors:
        text = str(error).lower()
        if "collection" in text or "pytest" in text:
            labels.append("collection")
        elif "native" in text:
            labels.append("native declaration")
        elif "metadata" in text or "declaration" in text:
            labels.append("declaration metadata")
        elif "validation" in text or "schema" in text:
            labels.append("validation")
        else:
            labels.append("inventory")
    counts = Counter(labels)
    return [f"{label}: {counts[label]} issue(s)" for label in sorted(counts)]


_INTERNAL_INPUT_KEYS = {
    "source", "path", "filename", "function", "helper", "command", "identity", "digest", "sha256",
    "artifact", "input_generation", "schema_version", "provenance", "reference_identity",
}


def _internal_input_key(key: Any) -> bool:
    """Hide lineage and machinery fields, including producer-specific suffixes."""

    value = str(key).lower()
    return (
        value.startswith(("_", "source_", "provenance_", "internal_"))
        or value in _INTERNAL_INPUT_KEYS
        or value.endswith(("_sha256", "_digest", "_hash"))
    )


def _public_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _public_value(nested) for key, nested in value.items()
                if not _internal_input_key(key)}
    if isinstance(value, list):
        return [_public_value(item) for item in value]
    return value


def _safe_comparison(comparison: dict[str, Any]) -> dict[str, Any]:
    visible = {
        "schema_version", "comparison_id", "quantity", "quantity_latex", "unit", "inputs",
        "inputs_latex", "coordinates", "basis", "method", "claim", "relationship_latex",
        "comparison_group_id", "point_id", "model", "reference", "criterion", "verdict",
        "judgment", "evidence_state", "unavailable_reason", "observations", "group_ref",
    }
    result = {key: _public_value(value) for key, value in comparison.items() if key in visible}
    if result.get("evidence_state") == "Unavailable":
        result["unavailable_reason"] = _public_comparison_reason(result.get("unavailable_reason"))
    return result


def _safe_declaration(declaration: dict[str, Any]) -> dict[str, Any]:
    visible = {
        "schema_version", "test_id", "title", "kind", "category", "lifecycle",
        "execution_policy", "evidence_expectation", "property_group", "requirement_id",
        "basis", "method", "claim", "expected_comparison_keys", "expected_comparison_counts",
        "parameter_cases", "parameter_case_id", "unclassified", "declaration_errors",
        "unmet_prerequisite", "scope",
    }
    result = {key: value for key, value in declaration.items() if key in visible}
    if result.get("unclassified"):
        result["title"] = "Unclassified check"
    if "declaration_errors" in result:
        result["declaration_errors"] = ["declaration metadata is incomplete"]
    return result


def companion_projection(
    plan: str,
    inventory: dict[str, Any],
    runs: list[dict[str, Any]] | None = None,
    *,
    revision_inputs: dict[str, Any] | None = None,
    layers: tuple[str, ...] = LAYERS,
) -> dict[str, Any]:
    raw_declarations = [item for item in inventory.get("declarations", []) if isinstance(item, dict)]
    declarations = [_safe_declaration(item) for item in raw_declarations if item.get("category") != "Studies"]
    declared_ids = {
        str(item.get("test_id")) for item in raw_declarations
        if item.get("test_id")
    }
    visible_runs = [
        _safe_run(item) for item in (runs or [])
        if isinstance(item, dict) and str(item.get("test_id")) in declared_ids
    ]
    studies = []
    by_test: dict[str, list[dict[str, Any]]] = {}
    for run in visible_runs:
        by_test.setdefault(str(run.get("test_id")), []).append(run)
    for raw in raw_declarations:
        if raw.get("category") != "Studies":
            continue
        declaration = _safe_declaration(raw)
        if declaration.get("category") != "Studies":
            continue
        entry = {
            "study_id": declaration.get("test_id"),
            "title": declaration.get("title", declaration.get("test_id")),
            "kind": declaration.get("kind", "Numerical"),
            "category": "Studies",
            "evidence_expectation": declaration.get("evidence_expectation", "Outcome only"),
            "lifecycle": declaration.get("lifecycle"),
            "execution_policy": declaration.get("execution_policy"),
            "run": by_test.get(str(declaration.get("test_id")), [None])[-1],
        }
        studies.append(entry)
    payload = {
        "schema_version": SCHEMA_VERSION,
        "plan": plan_projection(plan, layers),
        "preparation": preparation_projection(plan),
        "inventory": {
            "schema_version": inventory.get("schema_version", SCHEMA_VERSION),
            "state": inventory.get("state", "Unavailable"),
            "errors": _safe_error_summary(inventory.get("errors", [])),
        },
        "declarations": declarations,
        "runs": visible_runs,
        "studies": studies,
    }
    payload["revision"] = digest(revision_inputs or {
        "inventory": payload["inventory"],
        "runs": [(run.get("run_id"), run.get("freshness"), run.get("execution")) for run in visible_runs],
        "plan": plan,
    })
    return payload


def self_test() -> None:
    source = """
### Core
| C1 — Constants and numerical policy | [ ] | 1–2 | — | Start | Verify units and bounds. | — |
### Test-system preparation
| Work | Status | Completion |
|---|---|---|
| Recording | Planned | Capture individual outcomes. |
"""
    result = companion_projection(
        source,
        {"state": "Complete", "declarations": [{
            "schema_version": 1, "test_id": "x", "title": "One test", "kind": "Tooling",
            "category": "Core", "lifecycle": "Implemented", "execution_policy": "Explicit",
            "evidence_expectation": "Outcome only", "scope": "simulated",
            "runner": {"kind": "pytest"},
        }]},
        [{"schema_version": 1, "run_id": "r", "test_id": "x", "parameter_case_id": "single",
          "input_generation": "i", "artifact_identity": {"state": "Available"}, "execution": "Passed",
          "started_at": "s", "finished_at": "f", "comparisons": [], "freshness": "Current"}],
        layers=("Core",),
    )
    assert result["declarations"][0]["test_id"] == "x"
    assert result["declarations"][0]["scope"] == "simulated"
    assert result["runs"][0]["freshness"] == "Current"
    assert result["preparation"][0]["name"] == "Recording"
    assert result["plan"][0]["rows"][0]["id"] == "C1"


if __name__ == "__main__":
    self_test()
    print("companion data self-test passed")

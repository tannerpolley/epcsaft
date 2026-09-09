#!/usr/bin/env python3
"""Discover declared Python/native cases and planned scientific work.

Collection is the source of implemented Python declarations. Native
executables expose the same declaration records through ``--list``. This
module never executes a scientific case and never guesses a category from a
filename.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any
from urllib.parse import unquote

try:
    from tools.test_records import (
        RESULTS,
        ROOT,
        SCHEMA_VERSION,
        digest,
        input_generation,
        validate_inventory,
    )
except ModuleNotFoundError:  # direct ``python tools/numerical-inventory.py``
    from test_records import (
        RESULTS,
        ROOT,
        SCHEMA_VERSION,
        digest,
        input_generation,
        validate_inventory,
    )


NATIVE_TARGETS = {
    "contributions": "physics_checks",
    "properties": "model_checks",
    "diagnostics": "diagnostics_checks",
    "properties-floor": "model_floor",
    "equilibrium": "equilibrium_checks",
    "boundary": "boundary_checks",
    "specification": "specification_checks",
}

# Repository gates that produce an outcome even when they do not emit
# scientific comparisons.  The shell gate records these IDs through
# ``record-command`` so the inventory and the live companion retain one
# accountable row for each check.
TOOLING_GATES = (
    ("check-all", "equation-registry", "Equation registry consistency"),
    ("check-all", "roadmap", "Roadmap consistency"),
    ("check-all", "budgets", "Code budget consistency"),
    ("check-all", "environment", "Installed Engine environment"),
    ("check-all", "packets", "Reference packet validation"),
    ("check-all", "recorded-evidence", "Automatic numerical evidence"),
    ("check-all", "ruff", "Engine Ruff checks"),
    ("check-all", "mypy", "Engine mypy checks"),
    ("check-all", "analysis-layout", "Analysis layout checks"),
    ("check-all", "architecture", "Engine architecture checks"),
    ("native", "quality", "Native clang-tidy quality"),
    ("native", "sanitizers", "Native sanitizer checks"),
)


def _tooling_declarations() -> list[dict[str, Any]]:
    return [
        {
            "schema_version": SCHEMA_VERSION,
            "test_id": (
                f"tooling.check-all.{slug}" if group == "check-all"
                else f"tooling.native-{slug}"
            ),
            "title": title,
            "kind": "Tooling",
            "category": "Tooling/Checks",
            "lifecycle": "Implemented",
            # The ordinary check-all rows are part of the fast gate. Native
            # clang-tidy and sanitizer runs are owner-invoked quality checks
            # and must remain explicit despite sharing this declaration path.
            "execution_policy": "Automatic fast" if group == "check-all" else "Explicit",
            "evidence_expectation": "Outcome only",
            "source": (
                "tools/check-all.sh" if group == "check-all"
                else "engine/tools/check-native.sh" if slug == "quality"
                else "engine/tools/check-native-sanitizers.sh"
            ),
            "runner": {"kind": "command", "gate": slug},
        }
        for group, slug, title in TOOLING_GATES
    ]


_INTEGRITY_MODULE_NAME = "_epcsaft_check_test_system_declarations"
_integrity_source_digest: str | None = None


def _integrity_declarations() -> list[dict[str, Any]]:
    """Expose the harness registry without executing the checks."""

    global _integrity_source_digest
    path = ROOT / "tools/check-test-system.py"
    try:
        source = path.read_bytes()
    except OSError as error:
        raise RuntimeError(f"integrity harness source unavailable: {error}") from error
    source_digest = hashlib.sha256(source).hexdigest()
    module = sys.modules.get(_INTEGRITY_MODULE_NAME)
    if module is None or source_digest != _integrity_source_digest:
        path = ROOT / "tools/check-test-system.py"
        spec = importlib.util.spec_from_file_location(_INTEGRITY_MODULE_NAME, path)
        if spec is None or spec.loader is None:
            raise RuntimeError("integrity harness source loader unavailable")
        module = importlib.util.module_from_spec(spec)
        # ``dataclasses`` resolves annotations through ``sys.modules`` while
        # executing the module.  Register this temporary declaration module
        # before execution, then refresh it whenever the source digest moves.
        sys.modules[_INTEGRITY_MODULE_NAME] = module
        try:
            # Compile the bytes whose digest selected this registry.  Calling
            # ``SourceFileLoader.exec_module`` here could reuse a timestamp /
            # size-valid pyc after a same-size source replacement.
            exec(compile(source, str(path), "exec"), module.__dict__)  # noqa: S102
        except BaseException:
            if sys.modules.get(_INTEGRITY_MODULE_NAME) is module:
                del sys.modules[_INTEGRITY_MODULE_NAME]
            _integrity_source_digest = None
            raise
        _integrity_source_digest = source_digest
    checks = getattr(module, "INTEGRITY_CHECKS", ())
    if not isinstance(checks, (list, tuple)) or not checks:
        raise RuntimeError("integrity harness registry is empty or invalid")
    if any(not isinstance(getattr(check, "test_id", None), str) or not check.test_id
           for check in checks):
        raise RuntimeError("integrity harness registry contains an invalid declaration")
    return [
        {
            "schema_version": SCHEMA_VERSION,
            "test_id": check.test_id,
            "title": check.title,
            "kind": "Tooling",
            "category": "Tooling/Checks",
            "lifecycle": "Implemented",
            "execution_policy": "Automatic fast",
            "evidence_expectation": "Outcome only",
            "scope": check.kind,
            "source": "tools/check-test-system.py",
            "runner": {
                "kind": "integrity",
                "script": "tools/check-test-system.py",
                "check_id": check.test_id,
            },
        }
        for check in checks
    ]


def _link(path: Path) -> str:
    return "/" + str(path.relative_to(ROOT))


def _line_fields(line: str, prefix: str) -> dict[str, str]:
    text = line[len(prefix):].strip()
    matches = list(re.finditer(r"(?:^|\s)([A-Za-z_][A-Za-z0-9_]*)=", text))
    fields: dict[str, str] = {}
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        value = unquote(text[match.end():end].strip())
        if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
            value = value[1:-1]
        fields[match.group(1)] = value
    return fields


def _expected_counts(raw: str | None) -> dict[str, int] | None:
    if not raw:
        return None
    try:
        value = json.loads(raw)
    except json.JSONDecodeError:
        value = None
    if isinstance(value, dict):
        try:
            return {str(key): int(count) for key, count in value.items()}
        except (TypeError, ValueError):
            return None
    counts: dict[str, int] = {}
    for entry in raw.split(","):
        if ":" not in entry:
            return None
        key, count = entry.rsplit(":", 1)
        try:
            counts[key] = int(count)
        except ValueError:
            return None
    return counts


def _parse_decl_line(line: str) -> dict[str, Any] | None:
    if not line.startswith("EPCSAFT_DECL "):
        return None
    fields = _line_fields(line, "EPCSAFT_DECL ")
    if not fields.get("test_id"):
        return None
    required = (
        "schema_version", "test_id", "title", "kind", "category", "lifecycle",
        "execution_policy", "evidence_expectation",
    )
    missing = [name for name in required if not fields.get(name)]
    for key in ("execution_policy", "evidence_expectation"):
        fields[key] = fields.get(key, "").replace("_", " ")
    fields["schema_version"] = SCHEMA_VERSION
    fields["title"] = fields.get("title", fields["test_id"]).replace("_", " ")
    fields["kind"] = fields.get("kind", "Numerical")
    fields["category"] = fields.get("category", "Unclassified").replace("_", "/")
    fields["lifecycle"] = fields.get("lifecycle", "Implemented")
    if fields.get("expected_comparison_keys"):
        fields["expected_comparison_keys"] = fields["expected_comparison_keys"].split(",")
    elif fields.get("evidence_expectation") == "Numerical":
        fields["expected_comparison_keys"] = []
        missing.append("expected_comparison_keys")
    if fields.get("arguments"):
        fields["runner"] = {
            "kind": "native",
            "target": fields.get("runner"),
            "arguments": fields["arguments"].split(","),
        }
    elif fields.get("runner"):
        fields["runner"] = {"kind": "native", "target": fields["runner"], "arguments": []}
    counts = _expected_counts(fields.get("expected_comparison_counts"))
    if counts is not None:
        fields["expected_comparison_counts"] = counts
    if missing:
        fields["unclassified"] = True
        fields["declaration_errors"] = ["missing metadata: " + name for name in sorted(set(missing))]
    return fields


def _fallback(suite: str, target: str) -> dict[str, Any]:
    stable_suite = target.removesuffix("_checks")
    return {
        "schema_version": SCHEMA_VERSION,
        # Native producers own their stable IDs.  If the executable cannot be
        # listed, derive the historical fallback ID from its target while
        # keeping the public suite label above the runner boundary.
        "test_id": f"native.{stable_suite}",
        "title": f"Native {suite.replace('-', ' ').title()} suite",
        "kind": "Numerical",
        "category": "Unclassified",
        "lifecycle": "Implemented",
        "execution_policy": "Explicit",
        "evidence_expectation": "Outcome only",
        "source": f"engine/tests/native/{target}.cpp",
        "runner": {"kind": "native", "target": target, "arguments": []},
        "unclassified": True,
    }


def _native_declarations(build_dir: Path | None) -> tuple[list[dict[str, Any]], list[str]]:
    declarations: list[dict[str, Any]] = []
    errors: list[str] = []
    for suite, target in NATIVE_TARGETS.items():
        binary = build_dir / target if build_dir is not None else None
        if binary is None or not binary.is_file():
            declarations.append(_fallback(suite, target))
            errors.append(f"native declaration unavailable: {target} (build required)")
            continue
        try:
            result = subprocess.run([str(binary), "--list"], cwd=ROOT, check=False,
                                    capture_output=True, text=True, timeout=10)
        except (OSError, subprocess.TimeoutExpired) as error:
            result = None
            errors.append(f"native declaration failed: {target}: {error}")
        if result is None or result.returncode != 0:
            declarations.append(_fallback(suite, target))
            errors.append(f"native --list failed: {target}")
            continue
        found = [record for line in result.stdout.splitlines() if (record := _parse_decl_line(line)) is not None]
        if not found:
            errors.append(f"native --list emitted no declarations: {target}")
            declarations.append(_fallback(suite, target))
            continue
        declarations.extend(found)
        for declaration in found:
            errors.extend(
                f"{declaration.get('test_id', target)}.{error}"
                for error in declaration.get("declaration_errors", [])
            )
    return declarations, errors


def _collection() -> tuple[list[dict[str, Any]], list[str], str]:
    destination = Path(tempfile.mktemp(prefix="epcsaft-collection-", suffix=".json"))
    command = [str(ROOT / ".venv/bin/python"), "-m", "pytest", "--collect-only", "-q",
               "-p", "tools.pytest_recording", "engine/tests"]
    environment = dict(os.environ, EPCSAFT_COLLECTION_RESULTS=str(destination))
    try:
        result = subprocess.run(command, cwd=ROOT, env=environment, check=False,
                                capture_output=True, text=True, timeout=60)
        errors = [] if result.returncode == 0 else [line.strip() for line in result.stdout.splitlines()
                                                   if "ERROR" in line or "error" in line.lower()]
        try:
            value = json.loads(destination.read_text())
            declarations = list(value.get("declarations", []))
            errors.extend(str(error) for error in value.get("errors", []))
        except (OSError, json.JSONDecodeError):
            declarations = []
        if not declarations and result.returncode == 0:
            errors.append("pytest collection produced no declarations")
        return declarations, errors, result.stdout + result.stderr
    except (OSError, subprocess.TimeoutExpired) as error:
        return [], [f"pytest collection unavailable: {error}"], ""
    finally:
        destination.unlink(missing_ok=True)


def _planned() -> list[dict[str, Any]]:
    path = ROOT / "engine/tests/planned-declarations.json"
    if not path.is_file():
        return []
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        return [{
            "schema_version": SCHEMA_VERSION,
            "test_id": "planned.invalid",
            "title": "Invalid planned declarations",
            "kind": "Tooling",
            "category": "Unclassified",
            "lifecycle": "Planned",
            "execution_policy": "Explicit",
            "evidence_expectation": "Outcome only",
            "unmet_prerequisite": str(error),
        }]
    return list(value) if isinstance(value, list) else []


def _studies() -> list[dict[str, Any]]:
    analyses = ROOT / "analyses"
    if not analyses.is_dir():
        return []
    records: list[dict[str, Any]] = []
    for path in sorted(analyses.rglob("README.md")):
        relative = path.parent.relative_to(analyses)
        if not relative.parts:
            continue
        # Only an analysis owner's README is a study declaration.  READMEs
        # below retained inputs, generated results, plots, scripts or tests
        # describe an artifact or subdirectory, not another executable study.
        if any(part in {"inputs", "results", "figures", "scripts", "tests"} for part in relative.parts):
            continue
        source = path.read_text(encoding="utf-8")
        title = next(
            (line[2:].strip() for line in source.splitlines() if line.startswith("# ") and line[2:].strip()),
            " ".join(relative.parts),
        )
        artifacts = [
            candidate for folder in ("results", "figures")
            for candidate in (path.parent / folder).rglob("*")
            if candidate.is_file() and candidate.suffix.lower() in {".csv", ".tsv", ".json", ".html", ".pdf"}
        ]
        scripts = [
            candidate for folder in ("scripts", "tests")
            for candidate in (path.parent / folder).rglob("*.py")
            if candidate.is_file()
        ]
        test_files = [
            candidate for candidate in (path.parent / "tests").glob("*.py")
            if candidate.is_file()
        ]
        records.append({
            "schema_version": SCHEMA_VERSION,
            "test_id": "study." + ".".join(relative.parts),
            "title": title,
            "kind": "Numerical",
            "category": "Studies",
            # A retained study is a visible explicit owner even when its
            # replay has no executable entry point.  Missing outputs do not
            # invent a planned scientific declaration or prerequisite.
            "lifecycle": "Implemented",
            "execution_policy": "Explicit",
            # Existing study tables are retained artifacts until their own
            # producer adopts the shared comparison reporter.  A study replay
            # therefore records its execution outcome without inventing rows.
            "evidence_expectation": "Outcome only",
            "expected_comparison_keys": [],
            "source": _link(path),
            "runner": {
                "kind": "study",
                "scripts": [_link(item) for item in scripts],
                "cwd": str(path.parent.relative_to(ROOT)),
                "test_files": [str(item.relative_to(path.parent)) for item in test_files],
            },
            "artifacts": [_link(item) for item in sorted(artifacts)],
            "description": source[:4000],
        })
    return records


def inventory(*, write: bool = False, generation: str | None = None) -> dict[str, Any]:
    build_marker = ROOT / "build/environment-wheel/meson-build-dir"
    build_dir = Path(build_marker.read_text().strip()) if build_marker.is_file() else None
    python_declarations, collection_errors, collection_output = _collection()
    native_declarations, native_errors = _native_declarations(build_dir)
    integrity_errors: list[str] = []
    try:
        integrity_declarations = _integrity_declarations()
    except Exception as error:  # noqa: BLE001 - discovery errors belong in the inventory state
        integrity_declarations = []
        integrity_errors.append(f"integrity declaration discovery failed: {error}")
    declarations = (
        python_declarations + native_declarations + _tooling_declarations()
        + integrity_declarations + _planned() + _studies()
    )
    validation_errors = validate_inventory(declarations)
    errors = collection_errors + native_errors + integrity_errors + validation_errors
    result = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": __import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat(),
        "input_generation": generation if generation is not None else input_generation(ROOT),
        "inventory_generation": digest({"declarations": declarations, "errors": errors}),
        "state": "Complete" if not errors else "Incomplete",
        "errors": errors,
        "collection_output": collection_output[-20000:],
        "declarations": declarations,
        "tests": declarations,
        "studies": [record for record in declarations if record.get("category") == "Studies"],
    }
    if write:
        RESULTS.mkdir(parents=True, exist_ok=True)
        try:
            from tools.test_records import atomic_json
        except ModuleNotFoundError:
            from test_records import atomic_json
        atomic_json(RESULTS / "inventory.json", result)
    return result


def artifact(relative: str) -> dict[str, Any]:
    path = (ROOT / relative.lstrip("/")).resolve()
    analyses = ROOT / "analyses"
    if not path.is_relative_to(analyses) or not path.is_file():
        raise ValueError("Not a study artifact")
    if path.suffix.lower() in {".csv", ".tsv"}:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.reader(stream, delimiter="\t" if path.suffix.lower() == ".tsv" else ",")
            rows = []
            for index, row in enumerate(reader):
                if index >= 100:
                    return {"rows": rows, "truncated": True}
                rows.append(row)
            return {"rows": rows, "truncated": False}
    if path.suffix.lower() == ".json":
        text = path.read_text(encoding="utf-8")
        return {"text": json.dumps(json.loads(text), indent=2)[:20000], "truncated": len(text) > 20000}
    raise ValueError("Preview supports CSV, TSV and JSON")


if __name__ == "__main__":
    value = inventory(write=True)
    print(f"{len(value['declarations'])} declarations; state={value['state']}; errors={len(value['errors'])}")

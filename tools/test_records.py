"""Shared declarations, result records, identities, and atomic publication.

This module is deliberately small and dependency free.  Native tests, pytest,
the local worker, and the companion server all use the same records instead of
maintaining a second display-specific result format.
"""

from __future__ import annotations

import contextlib
import csv
import fcntl
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Iterator
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable

SCHEMA_VERSION = 1
ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "build" / "numerical-evidence"
MANIFEST_NAME = "manifest.json"
INVENTORY_NAME = "inventory.json"
WORKER_NAME = "worker.json"
ATTEMPT_ACTIVE_ENV = "EPCSAFT_RECORD_ATTEMPT_ACTIVE"

DECLARATION_REQUIRED = {
    "schema_version",
    "test_id",
    "title",
    "kind",
    "category",
    "lifecycle",
    "execution_policy",
    "evidence_expectation",
}
RUN_REQUIRED = {
    "schema_version",
    "run_id",
    "test_id",
    "parameter_case_id",
    "input_generation",
    "artifact_identity",
    "execution",
    "started_at",
    "finished_at",
}
COMPARISON_REQUIRED = {
    "comparison_id",
    "quantity",
    "unit",
    "inputs",
    "judgment",
    "evidence_state",
}

KINDS = {"Numerical", "Package", "Tooling"}
LIFECYCLES = {"Planned", "Implemented"}
EXECUTION_POLICIES = {"Automatic fast", "Explicit"}
EVIDENCE_EXPECTATIONS = {"Numerical", "Outcome only"}
EXECUTIONS = {
    "Not run",
    "Queued",
    "Running",
    "Passed",
    "Failed",
    "Skipped",
    "Error",
    "Interrupted",
}
INTERRUPTION_KINDS = {"timeout", "interruption", "abandoned"}
JUDGMENTS = {"Criterion satisfied", "Criterion failed", "Observation only", "Unavailable"}
EVIDENCE_STATES = {"Complete", "Partial", "Unavailable"}



def now() -> str:
    return datetime.now(UTC).isoformat()


def canonical(value: Any) -> str:
    """Serialize a record without platform-dependent ordering or whitespace."""

    return json.dumps(value, ensure_ascii=False, allow_nan=False, sort_keys=True, separators=(",", ":"))


def digest(value: Any) -> str:
    return hashlib.sha256(canonical(value).encode("utf-8")).hexdigest()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = canonical(value) + "\n"
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


@contextlib.contextmanager
def locked(path: Path = RESULTS / ".lock", *, blocking: bool = True) -> Iterator[Any]:
    """Hold the one process lock used by every recording entry point."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+") as stream:
        flags = fcntl.LOCK_EX
        if not blocking:
            flags |= fcntl.LOCK_NB
        fcntl.flock(stream.fileno(), flags)
        try:
            yield stream
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def _walk_files(root: Path, relative_roots: tuple[str, ...]) -> list[Path]:
    excluded = {".git", ".venv", "build", ".pytest_cache", ".mypy_cache", ".ruff_cache", "__pycache__"}
    files: list[Path] = []
    for relative_root in relative_roots:
        base = root / relative_root
        if base.is_file():
            files.append(base)
            continue
        if not base.is_dir():
            continue
        for path in base.rglob("*"):
            if path.is_file() and not any(part in excluded for part in path.relative_to(root).parts):
                files.append(path)
    return sorted(set(files))


def computational_paths(root: Path = ROOT) -> list[Path]:
    """Return the broad, conservative computational input closure.

    Presentation sources and generated outputs are intentionally excluded, so
    editing a page cannot invalidate a numerical run.  The closure is broad
    within the computational owners and therefore favors safe reruns over a
    fragile dependency graph.
    """

    roots = (
        "engine/native",
        "engine/src",
        "engine/tests",
        "data",
        "analyses",
        "tools",
        "engine/meson.build",
        "engine/pyproject.toml",
        "engine/LICENSE",
        "engine/README.md",
        "pyproject.toml",
        "uv.lock",
        ".codex/environments/environment.toml",
        ".codex/environments/setup.sh",
    )
    paths = _walk_files(root, roots)
    # Meson resolves these dependency wraps while preparing the Engine wheel.
    # They live below a generated ``build`` directory, which the broad walk
    # excludes, but their contents are still build inputs and must invalidate
    # retained evidence when they change.
    subprojects = root / "engine/build/subprojects"
    if subprojects.is_dir():
        paths.extend(path for path in subprojects.glob("*.wrap") if path.is_file())
        packagefiles = subprojects / "packagefiles"
        if packagefiles.is_dir():
            paths.extend(path for path in packagefiles.rglob("*") if path.is_file())
    ignored_suffixes = {".pyc", ".o", ".so", ".a", ".whl"}
    presentation_tools = {
        "tools/check-analysis-layout.py",
        "tools/check-roadmap.py",
        "tools/companion-data.py",
        "tools/diffstat.py",
        "tools/serve-roadmap.py",
    }
    def _retain_analysis_input(path: Path) -> bool:
        if not path.is_relative_to(root / "analyses"):
            return True
        parts = path.relative_to(root / "analyses").parts
        if "results" in parts:
            return False
        # Generated plots live below figures/, while source tables retained
        # under figures/*/source/ are consumed scientific inputs.
        return not ("figures" in parts and "source" not in parts)

    return [
        path for path in paths
        if path.suffix not in ignored_suffixes
        and not any(part in {"__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache"}
                    for part in path.relative_to(root).parts)
        and _retain_analysis_input(path)
        and path.name not in {"manifest.json", "inventory.json", "worker.json"}
        and str(path.relative_to(root)) not in presentation_tools
    ]


def input_manifest(root: Path = ROOT) -> dict[str, Any]:
    files = []
    for path in computational_paths(root):
        files.append({"path": str(path.relative_to(root)), "sha256": sha256(path), "size": path.stat().st_size})
    return {"schema_version": SCHEMA_VERSION, "files": files}


def input_generation(root: Path = ROOT) -> str:
    return digest(input_manifest(root))


def _command_version(command: list[str]) -> str | None:
    try:
        result = subprocess.run(command, check=False, capture_output=True, text=True, timeout=5)
    except (OSError, subprocess.TimeoutExpired):
        return None
    text = (result.stdout or result.stderr).strip().splitlines()
    return text[0] if text else None


def _runtime_report(runtime: Path) -> dict[str, Any] | None:
    """Report the interpreter and installed Engine modules used by a run.

    The process that prepares a wheel is not necessarily the interpreter that
    executes the checks (the repository's installed environment is Python
    3.13 while the system helper may be Python 3.12).  Ask the selected
    runtime itself which package files it imported so an archive alone cannot
    make a result appear applicable.
    """

    probe = (
        "import importlib.metadata, json, platform, sys; "
        "import epcsaft, epcsaft._core; "
        "print(json.dumps({" 
        "'executable': sys.executable, 'python': sys.version.split()[0], "
        "'implementation': platform.python_implementation(), "
        "'epcsaft_version': importlib.metadata.version('epcsaft'), "
        "'package': epcsaft.__file__, 'core': epcsaft._core.__file__}, "
        "sort_keys=True))"
    )
    try:
        result = subprocess.run(
            [str(runtime), "-c", probe], capture_output=True, text=True,
            check=False, timeout=10,
        )
        if result.returncode != 0:
            return None
        value = json.loads(result.stdout)
        # The venv launcher is commonly a symlink.  Store one canonical path
        # so the same selected interpreter has one identity at build and
        # applicability-check time.
        if isinstance(value, dict) and value.get("executable"):
            value["executable"] = str(Path(str(value["executable"])).resolve())
        return value if isinstance(value, dict) else None
    except (OSError, subprocess.TimeoutExpired, json.JSONDecodeError):
        return None


def _runtime_distribution_versions(runtime: Path, names: tuple[str, ...]) -> dict[str, str]:
    """Read package versions from the interpreter that executes the checks.

    The companion server may run from the system interpreter while evidence is
    produced by ``.venv/bin/python``.  Looking up metadata in the caller would
    therefore compare two different environments and mark every otherwise
    valid result unverified.  Query the selected runtime itself instead.
    """

    probe = """\
import importlib.metadata, json, sys
names = json.loads(sys.argv[1])
result = {}
for name in names:
    try:
        result[name] = importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        pass
print(json.dumps(result, sort_keys=True))
"""
    try:
        result = subprocess.run(
            [str(runtime), "-c", probe, json.dumps(list(names))],
            capture_output=True, text=True, check=False, timeout=10,
        )
        if result.returncode != 0:
            return {}
        value = json.loads(result.stdout)
        return value if isinstance(value, dict) and all(
            isinstance(key, str) and isinstance(version, str)
            for key, version in value.items()
        ) else {}
    except (OSError, subprocess.TimeoutExpired, json.JSONDecodeError):
        return {}


def _installed_runtime_files(package_path: Path) -> list[Path]:
    package_root = package_path.resolve().parent
    if not package_root.is_dir():
        return []
    return sorted(
        path for path in package_root.rglob("*")
        if path.is_file() and "__pycache__" not in path.parts
    )


def _meson_compiler(build_dir: Path | None) -> dict[str, Any] | None:
    if build_dir is None:
        return None
    path = build_dir / "meson-info" / "intro-compilers.json"
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        compiler = value.get("host", {}).get("cpp")
        return compiler if isinstance(compiler, dict) else None
    except (OSError, json.JSONDecodeError):
        return None


def _linked_dependencies(path: Path) -> list[dict[str, str]]:
    """Capture resolved ELF dependencies for a loaded extension/binary."""

    try:
        result = subprocess.run(["ldd", str(path)], capture_output=True, text=True,
                                check=False, timeout=10)
    except (OSError, subprocess.TimeoutExpired):
        return []
    dependencies: list[dict[str, str]] = []
    for line in result.stdout.splitlines():
        match = __import__("re").search(r"=>\s+(/[^ ]+)", line)
        if match is None:
            match = __import__("re").match(r"\s*(/[^ ]+)\s+\(", line)
        if match is None:
            continue
        dependency = Path(match.group(1))
        if dependency.is_file():
            item = {"path": str(dependency), "sha256": sha256(dependency)}
            dependencies.append(item)
    return sorted(dependencies, key=lambda item: item["path"])


def artifact_identity(
    root: Path = ROOT,
    *,
    build_dir: Path | None = None,
    wheel: Path | None = None,
    inventory_generation: str | None = None,
    runtime: Path | None = None,
    unavailable_reason: str | None = None,
) -> dict[str, Any]:
    root = root.resolve()
    if build_dir is not None:
        build_dir = build_dir if build_dir.is_absolute() else root / build_dir
        build_dir = build_dir.resolve()
    if wheel is not None:
        wheel = wheel if wheel.is_absolute() else root / wheel
        wheel = wheel.resolve()
    if runtime is not None:
        runtime = runtime if runtime.is_absolute() else root / runtime
        # Keep the venv launcher symlink.  Resolving it would invoke the base
        # interpreter outside the venv and lose the installed Engine package.
    if unavailable_reason:
        return {"state": "Unavailable", "reason": unavailable_reason}
    artifacts: list[dict[str, Any]] = []
    if wheel is not None and wheel.is_file():
        artifacts.append({"kind": "wheel", "path": str(wheel.relative_to(root)), "sha256": sha256(wheel)})
    if build_dir is not None and build_dir.is_dir():
        for path in sorted(build_dir.iterdir()):
            if path.is_file() and (path.name.endswith(".so") or path.name in {
                "model_checks", "equilibrium_checks", "physics_checks", "diagnostics_checks",
                "boundary_checks", "specification_checks", "model_floor",
            }):
                artifacts.append({"kind": "binary", "path": str(path.relative_to(root)), "sha256": sha256(path)})
    selected_runtime = runtime or (root / ".venv/bin/python")
    if not selected_runtime.is_file():
        selected_runtime = Path(sys.executable)
    runtime_report = _runtime_report(selected_runtime)
    if runtime_report is None:
        return {"state": "Unavailable", "reason": "installed Engine runtime could not be imported"}
    package_path = Path(str(runtime_report.get("package", ""))).resolve()
    runtime_files = _installed_runtime_files(package_path)
    if not runtime_files:
        return {"state": "Unavailable", "reason": "installed Engine package files are unavailable"}
    for path in runtime_files:
        try:
            relative = str(path.relative_to(root.resolve()))
        except ValueError:
            return {"state": "Unavailable", "reason": "installed Engine package is outside repository environment"}
        artifacts.append({"kind": "installed-runtime", "path": relative, "sha256": sha256(path)})
    for kind, key in (("installed-package", "package"), ("installed-extension", "core")):
        module_path = Path(str(runtime_report.get(key, ""))).resolve()
        if not module_path.is_file():
            return {"state": "Unavailable", "reason": f"installed Engine {key} is unavailable"}
        try:
            relative = str(module_path.relative_to(root.resolve()))
        except ValueError:
            return {"state": "Unavailable", "reason": f"installed Engine {key} is outside repository environment"}
        artifacts.append({"kind": kind, "path": relative, "sha256": sha256(module_path)})
    compiler = _meson_compiler(build_dir)
    if compiler is None:
        return {"state": "Unavailable", "reason": "Meson compiler identity unavailable"}
    compiler_executable = compiler.get("exelist", [])
    compiler_path = shutil.which(str(compiler_executable[-1])) if compiler_executable else None
    compiler_record = dict(compiler)
    compiler_record["resolved_executable"] = compiler_path
    if compiler_path and Path(compiler_path).is_file():
        compiler_record["executable_sha256"] = sha256(Path(compiler_path))
    loaded = [item for item in artifacts if item["kind"] in {"installed-extension", "binary"}]
    linked = {item["path"]: _linked_dependencies(root / item["path"])
              for item in loaded if (root / item["path"]).is_file()}
    selected_dependency_versions = _runtime_distribution_versions(
        selected_runtime, ("meson", "meson-python", "pybind11", "pytest")
    )
    environment: dict[str, Any] = {
        "python": runtime_report["python"],
        "runtime": str(selected_runtime),
        "implementation": runtime_report.get("implementation"),
        "platform": platform.platform(),
        "runtime_report": runtime_report,
        "compiler": compiler_record,
        "native_dependencies": selected_dependency_versions,
        "linked_dependencies": linked,
    }
    return {
        "state": "Available",
        "input_generation": input_generation(root),
        "inventory_generation": inventory_generation,
        "artifacts": artifacts,
        "environment": environment,
        "identity": digest({"artifacts": artifacts, "environment": environment, "inventory_generation": inventory_generation}),
    }


def _verify_artifact_identity(root: Path, identity: dict[str, Any]) -> tuple[bool, str | None]:
    root = root.resolve()
    if not isinstance(identity, dict):
        return False, "artifact identity unavailable"
    if identity.get("state") != "Available":
        return False, str(identity.get("reason", "artifact identity unavailable"))
    artifacts = identity.get("artifacts", [])
    if (not isinstance(artifacts, list)
            or any(not isinstance(item, dict) for item in artifacts)
            or not any(item.get("kind") == "wheel" for item in artifacts)):
        return False, "installed Engine wheel identity unavailable"
    environment = identity.get("environment")
    if not isinstance(environment, dict) or not environment.get("python") or not isinstance(environment.get("compiler"), dict):
        return False, "runtime environment identity unavailable"
    runtime = Path(str(environment.get("runtime", "")))
    report = _runtime_report(runtime)
    expected_report = environment.get("runtime_report")
    if report is None or report != expected_report:
        return False, "installed Engine runtime identity changed"
    expected_runtime_files = {
        str(item.get("path")) for item in artifacts
        if isinstance(item, dict) and item.get("kind") == "installed-runtime"
    }
    current_runtime_files = set()
    for path in _installed_runtime_files(Path(str(report.get("package", "")))):
        try:
            current_runtime_files.add(str(path.relative_to(root)))
        except ValueError:
            return False, "installed Engine package is outside repository environment"
    if not expected_runtime_files or current_runtime_files != expected_runtime_files:
        return False, "installed Engine package file set changed"
    compiler = environment["compiler"]
    if compiler.get("resolved_executable"):
        executable = Path(str(compiler["resolved_executable"]))
        if not executable.is_file() or (compiler.get("executable_sha256") and sha256(executable) != compiler["executable_sha256"]):
            return False, "Meson compiler executable changed"
    else:
        return False, "Meson compiler executable unavailable"
    native_dependencies = environment.get("native_dependencies", {})
    if not isinstance(native_dependencies, dict):
        return False, "native dependency identity unavailable"
    runtime_dependencies = _runtime_distribution_versions(
        runtime, tuple(str(name) for name in native_dependencies)
    )
    for name, version in native_dependencies.items():
        if runtime_dependencies.get(name) != version:
            return False, f"Python dependency identity changed: {name}"
    for item in artifacts:
        path = (root / str(item.get("path", ""))).resolve()
        if not path.is_file() or not path.is_relative_to(root.resolve()):
            return False, f"artifact missing: {item.get('path')}"
        if sha256(path) != item.get("sha256"):
            return False, f"artifact changed: {item.get('path')}"
    linked_dependencies = environment.get("linked_dependencies", {})
    if not isinstance(linked_dependencies, dict):
        return False, "linked native dependency identity unavailable"
    for dependencies in linked_dependencies.values():
        if not isinstance(dependencies, list) or any(not isinstance(item, dict) for item in dependencies):
            return False, "linked native dependency identity unavailable"
        for dependency in dependencies:
            path = Path(str(dependency.get("path", "")))
            if not path.is_file() or sha256(path) != dependency.get("sha256"):
                return False, f"linked native dependency changed: {path}"
    expected_identity = digest({
        "artifacts": artifacts,
        "environment": environment,
        "inventory_generation": identity.get("inventory_generation"),
    })
    if identity.get("identity") != expected_identity:
        return False, "artifact identity digest changed"
    return True, None


def verify_artifact_identity(root: Path, identity: dict[str, Any]) -> tuple[bool, str | None]:
    """Verify the recorded runtime, compiler, dependencies and artifact files (build/run time)."""

    return _verify_artifact_identity(root.resolve(), identity)


def artifact_files_match(root: Path, identity: dict[str, Any]) -> tuple[bool, str | None]:
    """Display-time check: the recorded wheel and binaries are still the files on disk."""

    root = root.resolve()
    if not isinstance(identity, dict) or identity.get("state") != "Available":
        return False, str(identity.get("reason", "artifact identity unavailable")) if isinstance(identity, dict) else "artifact identity unavailable"
    artifacts = identity.get("artifacts", [])
    if not isinstance(artifacts, list) or not any(isinstance(item, dict) and item.get("kind") == "wheel" for item in artifacts):
        return False, "installed Engine wheel identity unavailable"
    for item in artifacts:
        if not isinstance(item, dict):
            return False, "artifact identity unavailable"
        path = (root / str(item.get("path", ""))).resolve()
        if not path.is_file() or not path.is_relative_to(root):
            return False, f"artifact missing: {item.get('path')}"
        if sha256(path) != item.get("sha256"):
            return False, f"artifact changed: {item.get('path')}"
    return True, None


def scalar_judgment(model: float, reference: float, criterion: dict[str, Any]) -> bool:
    if not all(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
               for value in (model, reference)):
        return False
    if criterion.get("type", "scalar") != "scalar":
        verdict = criterion.get("verdict")
        if not isinstance(verdict, bool):
            raise ValueError("custom criterion must carry a boolean verdict")
        return verdict
    atol = criterion.get("atol")
    rtol = criterion.get("rtol")
    if not all(isinstance(value, (int, float)) and not isinstance(value, bool)
               and math.isfinite(value) and value >= 0 for value in (atol, rtol)):
        raise ValueError("scalar criterion requires finite nonnegative atol and rtol")
    rule = criterion.get("rule", "symmetric_close")
    error = abs(model - reference)
    if rule == "symmetric_close":
        return error <= max(atol, rtol * max(abs(model), abs(reference)))
    if rule == "reference_scale":
        return error <= atol + rtol * abs(reference)
    if rule == "reference_or_one":
        return error <= atol + rtol * max(1.0, abs(reference))
    raise ValueError(f"unknown scalar criterion rule: {rule}")


def make_comparison(
    comparison_id: str,
    *,
    quantity: str,
    unit: str,
    inputs: dict[str, Any],
    model: float | None = None,
    reference: float | None = None,
    criterion: dict[str, Any] | None = None,
    basis: str | None = None,
    method: str | None = None,
    claim: str | None = None,
    coordinates: dict[str, Any] | None = None,
    comparison_group_id: str | None = None,
    point_id: str | None = None,
    quantity_latex: str | None = None,
    relationship_latex: str | None = None,
    inputs_latex: str | None = None,
    observation: bool = False,
    unavailable_reason: str | None = None,
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "comparison_id": comparison_id,
        "quantity": quantity,
        "unit": unit,
        "inputs": inputs,
        "judgment": "Unavailable",
        "evidence_state": "Unavailable",
    }
    if basis is not None:
        row["basis"] = basis
    if method is not None:
        row["method"] = method
    if claim is not None:
        row["claim"] = claim
    if coordinates is not None:
        row["coordinates"] = coordinates
    if comparison_group_id is not None:
        row["comparison_group_id"] = comparison_group_id
    if point_id is not None:
        row["point_id"] = point_id
    if quantity_latex is not None:
        row["quantity_latex"] = quantity_latex
    if relationship_latex is not None:
        row["relationship_latex"] = relationship_latex
    if inputs_latex is not None:
        row["inputs_latex"] = inputs_latex
    if isinstance(model, (int, float)) and not isinstance(model, bool) and math.isfinite(model):
        row["model"] = model
    if isinstance(reference, (int, float)) and not isinstance(reference, bool) and math.isfinite(reference):
        row["reference"] = reference
    if criterion is not None:
        row["criterion"] = criterion
    if unavailable_reason:
        row["unavailable_reason"] = unavailable_reason
        return row
    if model is None or not isinstance(model, (int, float)) or isinstance(model, bool) or not math.isfinite(model):
        row["unavailable_reason"] = "model operand is missing or nonfinite"
        return row
    if reference is not None and (
        not isinstance(reference, (int, float)) or isinstance(reference, bool) or not math.isfinite(reference)
    ):
        row["unavailable_reason"] = "reference operand is missing or nonfinite"
        return row
    if observation or criterion is None:
        row["judgment"] = "Observation only"
        row["evidence_state"] = "Complete"
        return row
    if reference is None or not isinstance(reference, (int, float)) or isinstance(reference, bool) or not math.isfinite(reference):
        row["unavailable_reason"] = "reference operand is missing or nonfinite"
        return row
    row["verdict"] = bool(scalar_judgment(float(model), float(reference), criterion))
    row["judgment"] = "Criterion satisfied" if row["verdict"] else "Criterion failed"
    row["evidence_state"] = "Complete"
    return row


def validate_declaration(record: dict[str, Any]) -> list[str]:
    errors = sorted(DECLARATION_REQUIRED - record.keys())
    if record.get("schema_version") != SCHEMA_VERSION:
        errors.append("schema_version")
    if not isinstance(record.get("test_id"), str) or not record.get("test_id"):
        errors.append("test_id")
    if record.get("kind") not in KINDS:
        errors.append("kind")
    if record.get("lifecycle") not in LIFECYCLES:
        errors.append("lifecycle")
    if record.get("execution_policy") not in EXECUTION_POLICIES:
        errors.append("execution_policy")
    if record.get("evidence_expectation") not in EVIDENCE_EXPECTATIONS:
        errors.append("evidence_expectation")
    if record.get("lifecycle") == "Planned" and not record.get("unmet_prerequisite"):
        errors.append("unmet_prerequisite")
    if record.get("evidence_expectation") == "Numerical" and not record.get("expected_comparison_keys"):
        errors.append("expected_comparison_keys")
    counts = record.get("expected_comparison_counts")
    if counts is not None:
        if not isinstance(counts, dict):
            errors.append("expected_comparison_counts")
        else:
            for key, value in counts.items():
                if not isinstance(key, str) or not key or not isinstance(value, int) or isinstance(value, bool) or value < 0:
                    errors.append("expected_comparison_counts")
                    break
    return sorted(set(errors))


def validate_comparison(record: dict[str, Any]) -> list[str]:
    errors = sorted(COMPARISON_REQUIRED - record.keys())
    if record.get("schema_version") != SCHEMA_VERSION:
        errors.append("schema_version")
    if record.get("judgment") not in JUDGMENTS:
        errors.append("judgment")
    if record.get("evidence_state") not in EVIDENCE_STATES:
        errors.append("evidence_state")
    if record.get("evidence_state") == "Unavailable" and not record.get("unavailable_reason"):
        errors.append("unavailable_reason")
    if not isinstance(record.get("comparison_id"), str) or not record.get("comparison_id"):
        errors.append("comparison_id")
    if not isinstance(record.get("quantity"), str) or not record.get("quantity"):
        errors.append("quantity")
    if not isinstance(record.get("unit"), str):
        errors.append("unit")
    if not isinstance(record.get("inputs"), dict):
        errors.append("inputs")
    for name in ("model", "reference"):
        if name in record:
            value = record[name]
            if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value):
                errors.append(f"{name} must be finite")
    judgment = record.get("judgment")
    if judgment in {"Criterion satisfied", "Criterion failed"}:
        if not all(name in record for name in ("model", "reference", "criterion", "verdict")):
            errors.append("judged comparison requires operands, criterion and verdict")
        if not isinstance(record.get("criterion"), dict):
            errors.append("criterion")
        if not isinstance(record.get("verdict"), bool):
            errors.append("verdict")
        elif (judgment == "Criterion satisfied") != record["verdict"]:
            errors.append("judgment/verdict mismatch")
        if not record.get("basis") or not record.get("method") or not record.get("claim"):
            errors.append("scientific comparison requires basis, method and claim")
    if judgment == "Observation only" and (
        record.get("evidence_state") != "Complete" or "model" not in record
    ):
        errors.append("observation requires a finite model and complete evidence")
    if judgment == "Unavailable" and record.get("evidence_state") == "Complete":
        errors.append("unavailable judgment cannot have complete evidence")
    criterion = record.get("criterion")
    if isinstance(criterion, dict):
        criterion_type = criterion.get("type", "scalar")
        if criterion_type == "scalar":
            for name in ("atol", "rtol"):
                value = criterion.get(name)
                if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value) or value < 0:
                    errors.append(f"criterion.{name}")
        elif not isinstance(criterion.get("verdict"), bool) or not criterion.get("rule"):
            errors.append("custom criterion requires rule and verdict")
    return sorted(set(errors))


def validate_run(record: dict[str, Any]) -> list[str]:
    errors = sorted(RUN_REQUIRED - record.keys())
    if record.get("schema_version") != SCHEMA_VERSION:
        errors.append("schema_version")
    if record.get("execution") not in EXECUTIONS:
        errors.append("execution")
    interruption_kind = record.get("interruption_kind")
    if interruption_kind is not None:
        if record.get("execution") != "Interrupted":
            errors.append("interruption_kind requires Interrupted execution")
        elif interruption_kind not in INTERRUPTION_KINDS:
            errors.append("interruption_kind")
    if record.get("finished_at") is None and record.get("execution") not in {"Running", "Queued"}:
        errors.append("finished_at")
    if record.get("execution") in {"Failed", "Skipped", "Error", "Interrupted"} and not record.get("outcome_reason"):
        errors.append("outcome_reason")
    if (record.get("execution") in {"Passed", "Failed", "Skipped", "Error", "Interrupted"}
            and ("artifact_digests" not in record or not isinstance(record.get("artifact_digests"), list))):
        errors.append("artifact_digests")
    if not record.get("input_generation"):
        errors.append("input_generation")
    if not isinstance(record.get("artifact_identity"), dict):
        errors.append("artifact_identity")
    comparisons = record.get("comparisons", [])
    input_groups = record.get("input_groups", [])
    group_refs: set[str] = set()
    if not isinstance(input_groups, list):
        errors.append("input_groups")
        input_groups = []
    for group in input_groups:
        if not isinstance(group, dict):
            errors.append("input_group")
            continue
        group_ref = group.get("group_ref")
        if not isinstance(group_ref, str) or not group_ref:
            errors.append("input_group.group_ref")
        elif group_ref in group_refs:
            errors.append(f"duplicate input_group.group_ref: {group_ref}")
        group_refs.add(str(group_ref))
        if not isinstance(group.get("shared_inputs", {}), dict):
            errors.append(f"input_group.shared_inputs: {group_ref}")
    if not isinstance(comparisons, list):
        errors.append("comparisons")
    else:
        seen: set[str] = set()
        for row in comparisons:
            if not isinstance(row, dict):
                errors.append("comparison")
                continue
            errors.extend(f"comparison.{field}" for field in validate_comparison(row))
            group_ref = row.get("group_ref")
            if group_ref is not None and (
                not isinstance(group_ref, str) or not group_ref or group_ref not in group_refs
            ):
                errors.append(f"comparison.group_ref: {group_ref}")
            key = row.get("comparison_id")
            if key in seen:
                errors.append(f"duplicate comparison_id: {key}")
            seen.add(str(key))
        if record.get("execution") == "Passed":
            if any(row.get("judgment") == "Criterion failed" for row in comparisons):
                errors.append("passed execution has failed comparison")
            if any(row.get("judgment") == "Unavailable" or row.get("evidence_state") == "Unavailable"
                   for row in comparisons):
                errors.append("passed execution has unavailable comparison")
    return sorted(set(errors))


def validate_inventory(declarations: list[dict[str, Any]]) -> list[str]:
    errors: list[str] = []
    seen: set[str] = set()
    for declaration in declarations:
        errors.extend(f"{declaration.get('test_id', '?')}.{error}" for error in validate_declaration(declaration))
        test_id = declaration.get("test_id")
        if test_id in seen:
            errors.append(f"duplicate test_id: {test_id}")
        seen.add(str(test_id))
        cases = declaration.get("parameter_cases", [])
        if isinstance(cases, list):
            seen_cases: set[str] = set()
            for case in cases:
                if not isinstance(case, dict):
                    continue
                case_id = str(case.get("parameter_case_id", ""))
                if case_id in seen_cases:
                    errors.append(f"{test_id}.duplicate parameter_case_id: {case_id}")
                seen_cases.add(case_id)
    return sorted(set(errors))


def write_table(path: Path, comparisons: list[dict[str, Any]]) -> None:
    if not comparisons:
        return
    keys = sorted({key for row in comparisons for key in row})
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=keys, extrasaction="ignore")
            writer.writeheader()
            for row in comparisons:
                writer.writerow({key: canonical(value) if isinstance(value, (dict, list)) else value
                                 for key, value in row.items()})
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


class RecordStore:
    """Atomic run publication under the shared numerical-evidence lock."""

    def __init__(self, root: Path = ROOT):
        self.root = root
        self.results = root / "build" / "numerical-evidence"
        self.manifest_path = self.results / MANIFEST_NAME
        self.inventory_path = self.results / INVENTORY_NAME
        self.worker_path = self.results / WORKER_NAME

    def read_manifest(self) -> dict[str, Any]:
        if not self.manifest_path.is_file():
            return {"schema_version": SCHEMA_VERSION, "revision": "empty", "runs": [], "state": "Unavailable"}
        try:
            value = json.loads(self.manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            return {"schema_version": SCHEMA_VERSION, "revision": "corrupt", "runs": [],
                    "state": "Unavailable", "reason": str(error)}
        if not isinstance(value, dict) or value.get("schema_version") != SCHEMA_VERSION:
            return {"schema_version": SCHEMA_VERSION, "revision": "invalid", "runs": [],
                    "state": "Unavailable", "reason": "invalid manifest schema"}
        return value

    def read_inventory(self) -> dict[str, Any]:
        if not self.inventory_path.is_file():
            return {"schema_version": SCHEMA_VERSION, "state": "Unavailable", "declarations": [],
                    "errors": ["inventory unavailable"]}
        try:
            value = json.loads(self.inventory_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            return {"schema_version": SCHEMA_VERSION, "state": "Unavailable", "declarations": [],
                    "errors": [str(error)]}
        return value if isinstance(value, dict) else {"schema_version": SCHEMA_VERSION, "state": "Unavailable", "declarations": []}

    def publish_inventory(self, inventory: dict[str, Any]) -> None:
        atomic_json(self.inventory_path, inventory)

    def publish(self, run: dict[str, Any]) -> dict[str, Any]:
        errors = validate_run(run)
        if errors:
            raise ValueError("invalid run: " + ", ".join(errors))
        self.results.mkdir(parents=True, exist_ok=True)
        run_id = str(run["run_id"])
        run_dir = self.results / "runs" / run_id
        temporary_dir = self.results / "runs" / f".{run_id}.tmp"
        if temporary_dir.exists():
            shutil.rmtree(temporary_dir)
        temporary_dir.mkdir(parents=True, exist_ok=False)
        comparisons_path = temporary_dir / "comparisons.csv"
        write_table(comparisons_path, run.get("comparisons", []))
        run_for_disk = {key: value for key, value in run.items() if key != "comparisons"}
        run_for_disk["comparisons_artifact"] = str(comparisons_path.name) if comparisons_path.is_file() else None
        run_for_disk["comparisons_sha256"] = sha256(comparisons_path) if comparisons_path.is_file() else None
        if run.get("observations_artifact"):
            source = Path(str(run["observations_artifact"]))
            if not source.is_absolute():
                source = self.root / source
            if not source.is_file():
                raise ValueError(f"observations artifact missing: {source}")
            destination = temporary_dir / f"observations-{source.name}"
            shutil.copyfile(source, destination)
            run_for_disk["observations_artifact"] = destination.name
            run_for_disk["observations_sha256"] = sha256(destination)
        atomic_json(temporary_dir / "run.json", run_for_disk)
        if run_dir.exists():
            shutil.rmtree(run_dir)
        os.replace(temporary_dir, run_dir)
        manifest = self.read_manifest()
        existing_by_case: dict[tuple[Any, Any], dict[str, Any]] = {}
        for item in manifest.get("runs", []):
            key = (item.get("test_id"), item.get("parameter_case_id"))
            previous = existing_by_case.get(key)
            if previous is None or str(item.get("finished_at", "")) >= str(previous.get("finished_at", "")):
                existing_by_case[key] = item
        # The current manifest is an index of the latest record for each
        # collected test/case pair.  Older run directories remain inspectable,
        # but they must not create duplicate notebook rows after a rerun.
        existing_by_case.pop((run["test_id"], run["parameter_case_id"]), None)
        runs = list(existing_by_case.values())
        runs.append({
            "run_id": run_id,
            "test_id": run["test_id"],
            "parameter_case_id": run["parameter_case_id"],
            "execution": run["execution"],
            "input_generation": run["input_generation"],
            "artifact_identity": run["artifact_identity"],
            "path": str((run_dir / "run.json").relative_to(self.root)),
            "started_at": run["started_at"],
            "finished_at": run["finished_at"],
            "outcome_reason": run.get("outcome_reason"),
            "comparisons_sha256": run_for_disk["comparisons_sha256"],
            "input_groups_sha256": digest(run.get("input_groups", [])),
        })
        runs.sort(key=lambda item: (item.get("test_id", ""), item.get("parameter_case_id", ""), item.get("finished_at", "")))
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "revision": digest({"runs": runs, "updated_at": now()}),
            "updated_at": now(),
            "input_generation": input_generation(self.root),
            "runs": runs,
        }
        try:
            atomic_json(self.manifest_path, manifest)
        except Exception as error:
            # A failed replacement must fail closed.  If the atomic rename is
            # unavailable, remove the superseded case from a direct fallback
            # manifest so an older pass cannot remain current on disk.
            fallback = {
                "schema_version": SCHEMA_VERSION,
                "revision": "publication-failed",
                "updated_at": now(),
                "input_generation": manifest["input_generation"],
                "state": "Unavailable",
                "reason": f"manifest publication failed: {error}",
                "runs": runs,
            }
            try:
                self.manifest_path.write_text(canonical(fallback) + "\n", encoding="utf-8")
            except OSError:
                self.manifest_path.unlink(missing_ok=True)
            raise
        return manifest

    def reconcile_attempts(self, *, owner: bool | None = None) -> int:
        """Turn attempts left active by a killed runner into terminal records.

        A queued/running attempt is published before its child process starts,
        so a process death cannot leave an older passing case indexed as the
        current result.  The process that owns the recording attempt performs
        that recovery.  Nested recorders inherit the owner's environment and
        must leave its live Running record alone; otherwise a child would make
        an active parent appear interrupted during a batch.

        ``owner`` is explicit for callers that need to override the inherited
        environment.  By default, the active-attempt marker identifies a
        nested recorder and suppresses recovery while a top-level recorder
        still performs startup recovery under the shared lock.
        """

        if owner is None:
            owner = os.environ.get(ATTEMPT_ACTIVE_ENV) != "1"
        if not owner:
            return 0

        active = [
            item for item in self.read_manifest().get("runs", [])
            if isinstance(item, dict) and item.get("execution") in {"Queued", "Running"}
        ]
        for item in active:
            record = {
                "schema_version": SCHEMA_VERSION,
                "run_id": item.get("run_id", "abandoned-attempt"),
                "test_id": item.get("test_id", "unknown"),
                "parameter_case_id": item.get("parameter_case_id", "single"),
                "input_generation": item.get("input_generation", "unavailable"),
                "artifact_identity": item.get("artifact_identity", {"state": "Unavailable"}),
                "execution": "Interrupted",
                "interruption_kind": "abandoned",
                "started_at": item.get("started_at", now()),
                "finished_at": now(),
                "outcome_reason": "previous runner stopped before publishing a terminal outcome",
                "comparisons": [],
                "artifact_digests": item.get("artifact_digests", []),
                "command": item.get("command", []),
                "exit_code": None,
            }
            self.publish(record)
        return len(active)

    def read_runs(self, manifest: dict[str, Any] | None = None) -> list[dict[str, Any]]:
        runs: list[dict[str, Any]] = []
        selected_manifest = self.read_manifest() if manifest is None else manifest
        raw_runs = selected_manifest.get("runs", []) if isinstance(selected_manifest, dict) else []
        if not isinstance(raw_runs, list):
            # A malformed manifest must not be interpreted as an iterable of
            # characters or, worse, as a source of current passing records.
            return runs
        for item in raw_runs:
            if not isinstance(item, dict):
                continue
            path = (self.root / str(item.get("path", ""))).resolve()
            runs_root = (self.results / "runs").resolve()
            try:
                if not path.is_relative_to(runs_root):
                    raise ValueError("run record path is outside results")
                value = json.loads(path.read_text(encoding="utf-8"))
                if not isinstance(value, dict):
                    raise TypeError("run record is not an object")
                for key in (
                    "run_id", "test_id", "parameter_case_id", "execution", "input_generation",
                    "artifact_identity", "started_at", "finished_at", "outcome_reason",
                    "comparisons_sha256",
                ):
                    if key in item and value.get(key) != item.get(key):
                        raise ValueError(f"run record metadata mismatch: {key}")
                if value.get("input_groups") and item.get("input_groups_sha256") is None:
                    raise ValueError("input groups integrity metadata unavailable")
                if item.get("input_groups_sha256") is not None and digest(value.get("input_groups", [])) != item.get("input_groups_sha256"):
                    raise ValueError("input groups changed")
                if value.get("run_id") != item.get("run_id"):
                    raise ValueError("run record run_id does not match manifest")
                comparisons: list[dict[str, Any]] = []
                table_name = value.get("comparisons_artifact")
                if table_name:
                    table = (path.parent / str(table_name)).resolve()
                    if not table.is_relative_to(path.parent.resolve()) or not table.is_file():
                        raise ValueError("comparison table missing")
                    expected_digest = value.get("comparisons_sha256")
                    if expected_digest and sha256(table) != expected_digest:
                        raise ValueError("comparison table changed")
                    with table.open(newline="", encoding="utf-8") as stream:
                        for row in csv.DictReader(stream):
                            parsed: dict[str, Any] = {}
                            for key, cell in row.items():
                                # DictWriter emits an empty cell both when an
                                # optional comparison field is absent and
                                # when its value is None.  Keep those fields
                                # absent on read so a partial/unavailable row
                                # remains valid instead of turning a missing
                                # numeric operand into the invalid string
                                # ``""``.  Required textual fields stay
                                # present so malformed rows still fail closed.
                                if not cell and key not in {
                                    "comparison_id", "quantity", "unit", "judgment", "evidence_state",
                                }:
                                    continue
                                if key in {"unit", "quantity", "quantity_latex", "relationship_latex", "inputs_latex", "basis", "method", "claim", "judgment", "evidence_state"}:
                                    parsed[key] = cell
                                    continue
                                if key == "verdict" and cell:
                                    parsed[key] = cell.strip().lower() == "true"
                                    continue
                                try:
                                    parsed[key] = json.loads(cell) if cell else cell
                                except json.JSONDecodeError:
                                    parsed[key] = cell
                            comparisons.append(parsed)
                observations_name = value.get("observations_artifact")
                if observations_name:
                    observations = (path.parent / str(observations_name)).resolve()
                    if not observations.is_relative_to(path.parent.resolve()) or not observations.is_file():
                        raise ValueError("observations artifact missing")
                    expected_observations_digest = value.get("observations_sha256")
                    if expected_observations_digest and sha256(observations) != expected_observations_digest:
                        raise ValueError("observations artifact changed")
                value["comparisons"] = comparisons
                errors = validate_run(value)
                if errors:
                    raise ValueError("invalid run record: " + ", ".join(errors))
            except (OSError, TypeError, ValueError, json.JSONDecodeError) as error:
                value = dict(item)
                value.update({
                    "schema_version": SCHEMA_VERSION,
                    "run_id": item.get("run_id", "unknown"),
                    "test_id": item.get("test_id", "unknown"),
                    "parameter_case_id": item.get("parameter_case_id", "single"),
                    "execution": "Error",
                    "started_at": item.get("started_at", ""),
                    "finished_at": item.get("finished_at", now()),
                    # The manifest entry is not trustworthy once the record
                    # itself fails integrity checks.  Preserve the failure as
                    # an unavailable identity so freshness cannot promote a
                    # stale or partially read record to Current.
                    "artifact_identity": {
                        "state": "Unavailable",
                        "reason": "record integrity unavailable",
                    },
                    "input_generation": item.get("input_generation", "unavailable"),
                    "comparisons": [],
                    "outcome_reason": f"record unavailable: {error}",
                    "freshness": "Unverified",
                })
            if isinstance(value, dict):
                runs.append(value)
        return runs


def freshness(
    run: dict[str, Any], root: Path = ROOT, *, current_generation: str | None = None,
    verify: Callable[[Path, dict[str, Any]], tuple[bool, str | None]] = verify_artifact_identity,
) -> tuple[str, str | None]:
    if not run.get("input_generation") or not isinstance(run.get("artifact_identity"), dict):
        return "Unverified", "identity unavailable"
    current = current_generation if current_generation is not None else input_generation(root)
    if run["input_generation"] != current:
        return "Outdated", "computational inputs changed"
    if run.get("artifact_identity", {}).get("input_generation") not in {None, run["input_generation"]}:
        return "Unverified", "artifact was built from a different input generation"
    ok, reason = verify(root, run["artifact_identity"])
    return ("Current", None) if ok else ("Unverified", reason)


def display_status(
    run: dict[str, Any], root: Path = ROOT, *, current_generation: str | None = None,
) -> dict[str, Any]:
    """Freshness for display: sources and artifact files unchanged; the environment is not re-probed."""

    status, reason = freshness(run, root, current_generation=current_generation, verify=artifact_files_match)
    result = dict(run)
    result["freshness"] = status
    if reason:
        result["freshness_reason"] = reason
    return result


def parameter_case_id(params: dict[str, Any] | None = None) -> str:
    if not params:
        return "single"
    return "case-" + digest(params)[:16]


def comparison_selector_matches(comparison_id: str, selector: str) -> bool:
    """Match one comparison ID against an exact or explicit wildcard selector."""

    if selector.endswith("*"):
        return comparison_id.startswith(selector[:-1])
    return comparison_id == selector


def runner_command(path: str, *args: str) -> list[str]:
    return [path, *args]


if __name__ == "__main__":
    manifest = input_manifest(ROOT)
    assert digest(manifest) == input_generation(ROOT)
    assert scalar_judgment(1.0, 1.0, {"type": "scalar", "atol": 0.0, "rtol": 1e-12})
    row = make_comparison("demo", quantity="q", unit="1", inputs={"x": 1}, model=1.0, reference=1.0,
                          criterion={"type": "scalar", "atol": 0.0, "rtol": 1e-12},
                          basis="Relationship", method="identity", claim="demo identity")
    assert not validate_comparison(row)
    print("recording self-test passed")

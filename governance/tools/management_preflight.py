#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tomllib
import zipfile
from email.parser import Parser
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "manifests" / "management.json"
GENERATED = ROOT / ".management-state"


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def run_git(root: Path, *arguments: str) -> str:
    completed = subprocess.run(
        ["git", "-C", str(root), *arguments],
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def repository_state(name: str, record: dict[str, Any]) -> dict[str, Any]:
    path = Path(record["local_path"])
    if not path.exists():
        return {"name": name, "path": str(path), "present": False}
    worktrees: list[dict[str, str]] = []
    current: dict[str, str] = {}
    for line in run_git(path, "worktree", "list", "--porcelain").splitlines():
        if not line:
            if current:
                worktrees.append(current)
                current = {}
            continue
        key, _, value = line.partition(" ")
        current[key] = value
    if current:
        worktrees.append(current)
    return {
        "name": name,
        "path": str(path),
        "present": True,
        "head": run_git(path, "rev-parse", "HEAD"),
        "branch": run_git(path, "branch", "--show-current"),
        "remote": run_git(path, "remote", "get-url", "origin"),
        "status": run_git(path, "status", "--porcelain", "--untracked-files=all"),
        "worktrees": worktrees,
    }


def wheel_metadata(path: Path, base: Path) -> dict[str, Any]:
    headers: list[dict[str, str]] = []
    libraries: list[dict[str, str]] = []
    name = ""
    version = ""
    with zipfile.ZipFile(path) as wheel:
        members = wheel.namelist()
        metadata_names = [
            item for item in members if item.endswith(".dist-info/METADATA")
        ]
        if len(metadata_names) == 1:
            message = Parser().parsestr(
                wheel.read(metadata_names[0]).decode("utf-8", errors="replace")
            )
            name = message.get("Name", "")
            version = message.get("Version", "")
        for member in members:
            lowered = member.lower()
            payload: bytes | None = None
            if lowered.endswith((".h", ".hpp")):
                payload = wheel.read(member)
                headers.append(
                    {
                        "path": member,
                        "sha256": hashlib.sha256(payload).hexdigest(),
                    }
                )
            elif lowered.endswith((".a", ".so", ".dll", ".dylib")):
                payload = wheel.read(member)
                libraries.append(
                    {
                        "path": member,
                        "sha256": hashlib.sha256(payload).hexdigest(),
                    }
                )
    return {
        "path": str(path.relative_to(base)),
        "name": name,
        "version": version,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
        "headers": headers,
        "native_libraries": libraries,
        "retention": "unclassified",
    }


def artifact_state(config: dict[str, Any]) -> dict[str, Any]:
    project_root = Path(config["project_root"])
    candidates: set[Path] = set()
    artifact_root = Path(config["artifact_root"])
    if artifact_root.exists():
        candidates.update(artifact_root.rglob("*.whl"))
    for record in config["repositories"].values():
        local = Path(record["local_path"])
        if local.exists():
            candidates.update((local / "dist").glob("*.whl"))
    wheels = [
        wheel_metadata(path, project_root)
        for path in sorted(candidates)
        if path.is_file()
    ]
    return {
        "schema_version": 1,
        "project_root": str(project_root),
        "wheel_count": len(wheels),
        "total_bytes": sum(item["bytes"] for item in wheels),
        "wheels": wheels,
    }


def bundle_hash(directory: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(item for item in directory.rglob("*") if item.is_file()):
        digest.update(path.relative_to(directory).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def parameter_state(config: dict[str, Any]) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    eos = Path(config["repositories"]["eos"]["local_path"])
    eos_bundles = eos / "src" / "epcsaft" / "data" / "bundles"
    if eos_bundles.exists():
        for manifest in sorted(eos_bundles.glob("*/*/bundle.toml")):
            data = tomllib.loads(manifest.read_text(encoding="utf-8"))
            directory = manifest.parent
            records.append(
                {
                    "repository": "eos",
                    "path": str(directory.relative_to(eos)),
                    "bundle_id": data.get("bundle", {}).get("id", directory.parent.name),
                    "sha256": bundle_hash(directory),
                    "ownership": "allowlist_review",
                }
            )
    validation = Path(config["repositories"]["validation"]["local_path"])
    if validation.exists():
        for manifest in sorted(validation.rglob("bundle.toml")):
            if ".git" in manifest.parts or ".worktrees" in manifest.parts:
                continue
            directory = manifest.parent
            records.append(
                {
                    "repository": "validation",
                    "path": str(directory.relative_to(validation)),
                    "bundle_id": directory.name,
                    "sha256": bundle_hash(directory),
                    "ownership": "campaign_evidence",
                }
            )
    return {"schema_version": 1, "records": records}


def environment_state(config: dict[str, Any]) -> dict[str, Any]:
    environments: list[dict[str, Any]] = []
    for name, record in sorted(config["repositories"].items()):
        root = Path(record["local_path"])
        environment = root / ".venv"
        distributions: list[dict[str, Any]] = []
        if environment.exists():
            for metadata in sorted(environment.glob("lib/python*/site-packages/*.dist-info/METADATA")):
                message = Parser().parsestr(
                    metadata.read_text(encoding="utf-8", errors="replace")
                )
                distribution = message.get("Name", "")
                if "pcsaft" not in distribution.lower().replace("-", ""):
                    continue
                direct_url_path = metadata.parent / "direct_url.json"
                direct_url = (
                    load_json(direct_url_path) if direct_url_path.exists() else None
                )
                distributions.append(
                    {
                        "name": distribution,
                        "version": message.get("Version", ""),
                        "direct_url": direct_url,
                    }
                )
        environments.append(
            {
                "repository": name,
                "path": str(environment),
                "present": environment.exists(),
                "distributions": distributions,
            }
        )
    return {"schema_version": 1, "environments": environments}


def validate_config(config: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    repositories = config.get("repositories")
    if not isinstance(repositories, dict):
        return ["repositories must be an object"]
    required = {"management", "github", "eos", "equilibrium", "regression", "validation", "lab"}
    missing = required.difference(repositories)
    if missing:
        errors.append(f"missing repositories: {sorted(missing)}")
    if repositories.get("validation", {}).get("outbound_access") != ["validation"]:
        errors.append("Validation outbound access must be Validation only")
    eos_access = set(repositories.get("eos", {}).get("outbound_access", []))
    if eos_access != {"eos", "validation"}:
        errors.append("EOS outbound access must be EOS and Validation only")
    expected_top = {"management", "eos", "equilibrium", "regression", "validation", "lab"}
    for owner in ("equilibrium", "regression"):
        if set(repositories.get(owner, {}).get("outbound_access", [])) != expected_top:
            errors.append(f"{owner} must have full operational repository access")
    for item in config.get("persistent_worktrees", []):
        target = item.get("target")
        accessor = item.get("accessor")
        if target not in repositories or accessor not in repositories:
            errors.append(f"invalid worktree route {accessor}->{target}")
        elif target not in repositories[accessor]["outbound_access"]:
            errors.append(f"forbidden worktree route {accessor}->{target}")
        if not str(item.get("relative_path", "")).endswith("-access"):
            errors.append(f"persistent worktree must end in -access: {item}")
    for name, record in repositories.items():
        if record.get("codex_environment") != ".codex/environments/environment.toml":
            errors.append(f"{name} must declare the standard Codex environment path")
    return errors


def snapshot(config: dict[str, Any]) -> None:
    states = [
        repository_state(name, record)
        for name, record in sorted(config["repositories"].items())
    ]
    write_json(
        GENERATED / "ecosystem-state.json",
        {"schema_version": 1, "repositories": states},
    )
    write_json(GENERATED / "artifact-ledger.json", artifact_state(config))
    write_json(GENERATED / "environment-ledger.json", environment_state(config))
    write_json(GENERATED / "parameter-ledger.json", parameter_state(config))


def check(config: dict[str, Any], strict: bool) -> list[str]:
    errors = validate_config(config)
    for name in (
        "ecosystem-state.json",
        "artifact-ledger.json",
        "environment-ledger.json",
        "parameter-ledger.json",
    ):
        path = GENERATED / name
        if not path.exists():
            errors.append(f"missing generated manifest {path.relative_to(ROOT)}")
    if strict and (GENERATED / "artifact-ledger.json").exists():
        ledger = load_json(GENERATED / "artifact-ledger.json")
        unknown = [
            item["path"]
            for item in ledger["wheels"]
            if item["retention"] == "unclassified"
        ]
        if unknown:
            errors.append(f"{len(unknown)} wheel artifacts remain unclassified")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("snapshot", "check"))
    parser.add_argument("--strict", action="store_true")
    arguments = parser.parse_args()
    config = load_json(CONFIG_PATH)
    if arguments.action == "snapshot":
        snapshot(config)
    elif not GENERATED.exists():
        snapshot(config)
    errors = check(config, arguments.strict)
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print("management_preflight: clean")
    return 0


if __name__ == "__main__":
    sys.exit(main())

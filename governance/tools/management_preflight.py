#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import tomllib
from typing import Any


GOVERNANCE = Path(__file__).resolve().parents[1]
ROOT = GOVERNANCE.parent
CONFIG = GOVERNANCE / "manifests" / "management.json"
GENERATED = ROOT / ".management-state" / "monorepo.json"
AREAS = ("governance", "engine", "data", "validation", "lab")
LIVE_POINTERS = tuple(
    path
    for area in AREAS
    for path in (
        ROOT / area / "AGENTS.md",
        ROOT / area / "CONTEXT.md",
        ROOT / area / "README.md",
    )
) + (ROOT / "AGENTS.md", ROOT / "CONTEXT.md", ROOT / "README.md")
RETIRED_ACTIVE_PATTERNS = (
    "../ePC-SAFT-governance",
    "ePC-SAFT-project/ePC-SAFT-",
    "Current issue tracking and runtime authority belong to",
)


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def git(*arguments: str) -> str:
    return subprocess.run(
        ["git", "-C", str(ROOT), *arguments],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def snapshot() -> dict[str, Any]:
    config = load_json(CONFIG)
    return {
        "schema_version": 1,
        "repository": config["repository"],
        "doctrine_revision": config["doctrine_revision"],
        "git_root": str(ROOT),
        "head": git("rev-parse", "HEAD"),
        "branch": git("branch", "--show-current"),
        "status": git("status", "--porcelain=v1"),
        "areas": {area: str(ROOT / area) for area in AREAS},
        "artifact_root": config["artifact_root"],
    }


def write_snapshot(value: dict[str, Any]) -> None:
    GENERATED.parent.mkdir(parents=True, exist_ok=True)
    temporary = GENERATED.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(GENERATED)


def check() -> list[str]:
    errors: list[str] = []
    try:
        config = load_json(CONFIG)
    except (OSError, json.JSONDecodeError) as error:
        return [f"invalid management manifest: {error}"]

    if config.get("schema_version") != 2:
        errors.append("management schema_version must be 2")
    if config.get("doctrine_revision") != 11:
        errors.append("management doctrine_revision must be 11")
    if config.get("repository") != "tannerpolley/ePC-SAFT":
        errors.append("management repository must be tannerpolley/ePC-SAFT")
    if config.get("ci_cd") != "disabled":
        errors.append("CI/CD must remain disabled")

    try:
        if Path(git("rev-parse", "--show-toplevel")) != ROOT:
            errors.append("command Git root does not match the monorepo root")
    except subprocess.CalledProcessError as error:
        errors.append(f"Git root check failed: {error}")

    for area in AREAS:
        path = ROOT / area
        if not path.is_dir():
            errors.append(f"missing owned directory: {area}")
        if not (path / "AGENTS.md").is_file():
            errors.append(f"missing directory agent contract: {area}/AGENTS.md")
    nested = [
        path
        for path in ROOT.rglob(".git")
        if path != ROOT / ".git"
        and not {"build", ".venv"}.intersection(path.relative_to(ROOT).parts)
    ]
    if nested:
        errors.append(f"nested Git roots are forbidden: {nested}")
    if (ROOT / ".github" / "workflows").exists():
        errors.append("CI/CD workflows are forbidden")
    if not (ROOT / ".codex" / "environments" / "environment.toml").is_file():
        errors.append("missing sole root Codex environment")

    try:
        engine = tomllib.loads((ROOT / "engine" / "pyproject.toml").read_text())
        if engine.get("project", {}).get("name") != "epcsaft":
            errors.append("engine/ must own the sole epcsaft distribution")
    except (OSError, tomllib.TOMLDecodeError) as error:
        errors.append(f"invalid engine pyproject: {error}")

    for manifest in sorted((GOVERNANCE / "manifests").glob("*.json")):
        try:
            load_json(manifest)
        except (OSError, json.JSONDecodeError) as error:
            errors.append(f"invalid JSON manifest {manifest.name}: {error}")

    for pointer in LIVE_POINTERS:
        if not pointer.is_file():
            continue
        text = pointer.read_text(encoding="utf-8")
        for pattern in RETIRED_ACTIVE_PATTERNS:
            if pattern in text:
                errors.append(f"retired active path in {pointer.relative_to(ROOT)}: {pattern}")

    environment_python = ROOT / ".venv" / "bin" / "python"
    if environment_python.is_file():
        probe = subprocess.run(
            [
                str(environment_python),
                "-c",
                "from importlib import metadata; import re; "
                "print('\\n'.join(sorted(d.metadata.get('Name','') for d in metadata.distributions() "
                "if re.sub(r'[-_.]+','-',d.metadata.get('Name','')).lower().startswith('epcsaft'))))",
            ],
            capture_output=True,
            text=True,
        )
        installed = [line for line in probe.stdout.splitlines() if line]
        if probe.returncode or installed != ["epcsaft"]:
            errors.append(f"root environment must contain only epcsaft, found {installed}")

    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("snapshot", "check"))
    arguments = parser.parse_args()
    if arguments.action == "snapshot":
        write_snapshot(snapshot())
        print(f"management_preflight: wrote {GENERATED}")
        return 0
    errors = check()
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print("management_preflight: clean")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

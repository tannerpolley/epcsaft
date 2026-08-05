#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tomllib
from typing import Any


GOVERNANCE = Path(__file__).resolve().parents[1]
ROOT = GOVERNANCE.parent
CONFIG = GOVERNANCE / "manifests" / "management.json"
GENERATED = ROOT / ".management-state" / "monorepo.json"
RECEIPT_SCHEMA = GOVERNANCE / "manifests" / "operation-receipt-envelope-v1.json"
RECEIPTS = GOVERNANCE / "receipts"
HEX64 = re.compile(r"^[0-9a-f]{64}$")
AREA_ROLES = {
    "governance": "doctrine_and_control_records",
    "engine": "runtime_and_compact_numerical_evidence",
    "data": "canonical_scientific_packets",
    "validation": "installed_artifact_campaigns",
    "lab": "experiments_and_independent_oracles",
}
AREAS = tuple(AREA_ROLES)
CONTROL_MANIFESTS = {
    "artifact_dispositions": "manifests/artifact-dispositions.json",
    "parameter_ownership": "manifests/parameter-ownership.json",
    "data_packet_locks": "manifests/data-packet-locks.json",
    "runtime_consolidation": "manifests/runtime-consolidation.json",
    "approved_operations": "manifests/approved-operations.json",
}
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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def valid_hash(value: Any) -> bool:
    return isinstance(value, str) and HEX64.fullmatch(value) is not None


def replayable(value: Any) -> bool:
    return isinstance(value, str) and not Path(value).is_absolute() and ".." not in Path(value).parts


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


def receipt_errors() -> list[str]:
    try:
        schema = load_json(RECEIPT_SCHEMA)
    except (OSError, json.JSONDecodeError) as error:
        return [f"invalid receipt envelope schema: {error}"]
    schema_hash = sha256_file(RECEIPT_SCHEMA)
    errors: list[str] = []
    sequences: list[tuple[int, Path]] = []
    for path in sorted(RECEIPTS.glob("*.json")):
        try:
            receipt = load_json(path)
        except (OSError, json.JSONDecodeError) as error:
            errors.append(f"invalid receipt {path.name}: {error}")
            continue
        if "envelope" not in receipt and "record" not in receipt:
            # Pre-envelope records are immutable historical provenance.
            if "schema_version" in receipt:
                continue
            errors.append(f"receipt {path.name} is missing its envelope")
            continue
        envelope = receipt.get("envelope")
        record = receipt.get("record")
        if not isinstance(envelope, dict) or not isinstance(record, dict):
            errors.append(f"receipt {path.name} must contain an envelope and record")
            continue
        errors.extend(
            f"receipt {path.name} missing envelope field: {field}"
            for field in schema.get("required_envelope_fields", [])
            if field not in envelope
        )
        errors.extend(
            f"receipt {path.name} missing record field: {field}"
            for field in schema.get("required_record_fields", [])
            if field not in record
        )
        sequence = envelope.get("sequence")
        if not isinstance(sequence, int) or isinstance(sequence, bool) or sequence < 1:
            errors.append(f"receipt {path.name} has an invalid sequence")
        else:
            sequences.append((sequence, path))
        schema_ref = envelope.get("schema")
        if not isinstance(schema_ref, dict):
            errors.append(f"receipt {path.name} has an invalid schema reference")
        else:
            if schema_ref.get("path") != "manifests/operation-receipt-envelope-v1.json":
                errors.append(f"receipt {path.name} has an invalid schema path")
            if (schema_ref.get("version"), schema_ref.get("sha256")) != (
                schema.get("schema_version"),
                schema_hash,
            ):
                errors.append(f"receipt {path.name} has a stale schema reference")
        subject = envelope.get("subject")
        if not isinstance(subject, dict):
            errors.append(f"receipt {path.name} has an invalid subject")
        else:
            if subject.get("record_sha256") != canonical_sha256(record):
                errors.append(f"receipt {path.name} record integrity mismatch")
            if envelope.get("subject_sha256") != canonical_sha256(subject):
                errors.append(f"receipt {path.name} subject integrity mismatch")
        if not valid_hash(envelope.get("subject_sha256")):
            errors.append(f"receipt {path.name} requires a 64-hex subject hash")
        validator = envelope.get("validator")
        if not isinstance(validator, dict):
            errors.append(f"receipt {path.name} has an invalid validator")
        else:
            if any(not isinstance(validator.get(field), str) or not validator[field] for field in ("repository", "commit", "path")):
                errors.append(f"receipt {path.name} validator identity is incomplete")
            if not re.fullmatch(r"[0-9a-f]{40}", validator.get("commit", "")):
                errors.append(f"receipt {path.name} validator commit is invalid")
            if not valid_hash(validator.get("sha256")):
                errors.append(f"receipt {path.name} validator hash is invalid")
            if not replayable(validator.get("path")):
                errors.append(f"receipt {path.name} validator path is not replayable")
        if not isinstance(envelope.get("command"), str) or not envelope["command"]:
            errors.append(f"receipt {path.name} requires an exact command")
        if not isinstance(envelope.get("environment"), dict):
            errors.append(f"receipt {path.name} requires an environment")
        artifacts = envelope.get("artifacts")
        if not isinstance(artifacts, list):
            errors.append(f"receipt {path.name} artifacts must be a list")
        else:
            for artifact in artifacts:
                if not isinstance(artifact, dict) or not valid_hash(artifact.get("sha256")):
                    errors.append(f"receipt {path.name} contains an invalid artifact hash")
                elif any(field in artifact and not replayable(artifact[field]) for field in ("path", "filename")):
                    errors.append(f"receipt {path.name} contains a non-replayable artifact path")
        producer = envelope.get("producer")
        reviewer = envelope.get("reviewer")
        if any(
            not isinstance(value, dict)
            or not isinstance(value.get("kind"), str)
            or not isinstance(value.get("id"), str)
            for value in (producer, reviewer)
        ):
            errors.append(f"receipt {path.name} requires producer and reviewer identities")
        if isinstance(producer, dict) and isinstance(reviewer, dict) and record.get("status") == "completed":
            if (producer.get("kind"), producer.get("id")) == (reviewer.get("kind"), reviewer.get("id")):
                errors.append(f"receipt {path.name} cannot be self-reviewed")
        if not isinstance(envelope.get("failed_outcomes"), list):
            errors.append(f"receipt {path.name} failed_outcomes must be a list")
        for field in ("supersedes", "revoked_by"):
            if envelope.get(field) is not None and not isinstance(envelope.get(field), str):
                errors.append(f"receipt {path.name} has an invalid {field} value")
        approval = envelope.get("user_approval")
        if not isinstance(approval, dict) or not isinstance(approval.get("source"), str) or not approval["source"]:
            errors.append(f"receipt {path.name} requires user approval")
        elif approval.get("subject_sha256") != envelope.get("subject_sha256"):
            errors.append(f"receipt {path.name} user approval is not subject-bound")
        if record.get("status") not in schema.get("allowed_record_statuses", []):
            errors.append(f"receipt {path.name} has an invalid record status")
    for left, right in zip(sorted(sequences), sorted(sequences)[1:]):
        if left[0] == right[0]:
            errors.append(f"receipt sequence {left[0]} is duplicated")
    return errors


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

    areas = config.get("areas")
    if not isinstance(areas, dict) or set(areas) != set(AREAS):
        errors.append("management areas must declare exactly the owned directories")
    else:
        for area, role in AREA_ROLES.items():
            if areas.get(area) != {"path": area, "role": role}:
                errors.append(f"management area declaration drift: {area}")
    for key, relative in CONTROL_MANIFESTS.items():
        if config.get(key) != relative:
            errors.append(f"management control manifest declaration drift: {key}")
        elif not (GOVERNANCE / relative).is_file():
            errors.append(f"missing control manifest: {relative}")

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
    workflows: set[Path] = set()
    try:
        tracked = git("ls-files", "-z").split("\0")
    except subprocess.CalledProcessError as error:
        errors.append(f"tracked workflow check failed: {error}")
        tracked = []
    for name in tracked:
        relative = Path(name)
        parts = relative.parts
        if {"build", ".venv"}.intersection(parts):
            continue
        if any(parts[index:index + 2] == (".github", "workflows") for index in range(len(parts) - 1)):
            workflows.add(ROOT / relative)
    for directory, directories, _ in os.walk(ROOT):
        relative = Path(directory).relative_to(ROOT)
        directories[:] = [name for name in directories if name not in {".git", "build", ".venv"}]
        if Path(directory).name == ".github" and "workflows" in directories:
            workflows.add(Path(directory) / "workflows")
    if workflows:
        errors.append("CI/CD workflows are forbidden: " + ", ".join(sorted(str(path.relative_to(ROOT)) for path in workflows)))
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

    errors.extend(receipt_errors())

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

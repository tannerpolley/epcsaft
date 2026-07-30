#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import zipfile
from email.parser import Parser
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
MANAGEMENT = ROOT / "manifests" / "management.json"
POLICY = ROOT / "manifests" / "artifact-store-policy.json"
TRACKED_INDEX = ROOT / "manifests" / "artifacts.json"


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def wheel_identity(path: Path) -> tuple[str, str]:
    with zipfile.ZipFile(path) as wheel:
        metadata = [
            name for name in wheel.namelist() if name.endswith(".dist-info/METADATA")
        ]
        if len(metadata) != 1:
            return "", ""
        message = Parser().parsestr(
            wheel.read(metadata[0]).decode("utf-8", errors="replace")
        )
    return message.get("Name", ""), message.get("Version", "")


def paths() -> tuple[Path, Path, dict[str, Any]]:
    management = load_json(MANAGEMENT)
    policy = load_json(POLICY)
    project_root = Path(management["project_root"])
    artifact_root = Path(management["artifact_root"])
    return project_root, artifact_root, policy


def aliases(artifact_root: Path) -> list[Path]:
    object_root = artifact_root / "objects"
    return sorted(
        path
        for path in artifact_root.rglob("*.whl")
        if path.is_file() and object_root not in path.parents
    )


def build_index(artifact_root: Path, policy: dict[str, Any]) -> dict[str, Any]:
    current = set(policy["current_subjects"].values())
    grouped: dict[str, list[Path]] = {}
    for path in aliases(artifact_root):
        digest = sha256(path)
        grouped.setdefault(digest, []).append(path)

    objects: list[dict[str, Any]] = []
    for digest, source_paths in sorted(grouped.items()):
        first = source_paths[0]
        name, version = wheel_identity(first)
        object_path = (
            Path(policy["object_directory"]) / digest / first.name
        )
        objects.append(
            {
                "sha256": digest,
                "distribution": name,
                "version": version,
                "bytes": first.stat().st_size,
                "object": object_path.as_posix(),
                "aliases": [
                    path.relative_to(artifact_root).as_posix()
                    for path in source_paths
                ],
                "retention": "current" if digest in current else "retained-evidence",
            }
        )
    return {
        "schema_version": 1,
        "policy": "ePC-SAFT/ePC-SAFT-management:manifests/artifact-store-policy.json",
        "object_count": len(objects),
        "alias_count": sum(len(item["aliases"]) for item in objects),
        "objects": objects,
    }


def write_index(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def normalize() -> list[str]:
    _, artifact_root, policy = paths()
    index = build_index(artifact_root, policy)
    errors: list[str] = []
    for item in index["objects"]:
        source = artifact_root / item["aliases"][0]
        destination = artifact_root / item["object"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            if sha256(destination) != item["sha256"]:
                errors.append(f"object hash mismatch: {destination}")
            continue
        try:
            os.link(source, destination)
        except OSError:
            shutil.copy2(source, destination)
        if sha256(destination) != item["sha256"]:
            errors.append(f"object creation hash mismatch: {destination}")
    if not errors:
        write_index(artifact_root / policy["index"], index)
        write_index(TRACKED_INDEX, index)
    return errors


def check() -> list[str]:
    project_root, artifact_root, policy = paths()
    index_path = artifact_root / policy["index"]
    if not index_path.exists():
        return [f"missing artifact index: {index_path}"]
    recorded = load_json(index_path)
    if not TRACKED_INDEX.exists():
        return [f"missing tracked artifact index: {TRACKED_INDEX}"]
    tracked = load_json(TRACKED_INDEX)
    expected = build_index(artifact_root, policy)
    errors: list[str] = []
    if recorded != expected:
        errors.append("artifact index is stale; run artifact_store.py normalize")
    if tracked != expected:
        errors.append("tracked artifact index is stale; run artifact_store.py normalize")
    allowed = set(policy["allowed_retention"])
    for item in recorded.get("objects", []):
        if item.get("retention") not in allowed:
            errors.append(f"invalid retention for {item.get('sha256')}")
        object_path = artifact_root / item["object"]
        if not object_path.exists():
            errors.append(f"missing object: {object_path}")
        elif sha256(object_path) != item["sha256"]:
            errors.append(f"object hash mismatch: {object_path}")
    for wheel in project_root.rglob("dist/*.whl"):
        if artifact_root in wheel.parents or ".venv" in wheel.parts:
            continue
        errors.append(f"mutable repository dist wheel is forbidden: {wheel}")
    current = set(policy["current_subjects"].values())
    indexed = {item["sha256"] for item in recorded.get("objects", [])}
    missing_current = sorted(current - indexed)
    if missing_current:
        errors.append(f"missing current artifact subjects: {missing_current}")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("normalize", "check"))
    arguments = parser.parse_args()
    errors = normalize() if arguments.action == "normalize" else check()
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print(f"artifact_store: {arguments.action} clean")
    return 0


if __name__ == "__main__":
    sys.exit(main())

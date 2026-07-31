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
DISPOSITIONS = ROOT / "manifests" / "artifact-dispositions.json"


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


def all_wheels(artifact_root: Path) -> list[Path]:
    object_root = artifact_root / "objects"
    return sorted(
        path
        for path in artifact_root.rglob("*.whl")
        if path.is_file() and object_root not in path.parents
    )


def active_subjects(policy: dict[str, Any]) -> dict[str, str]:
    return {
        distribution: subject["sha256"]
        for distribution, subject in policy["active_distributions"].items()
    }


def disposition_errors(artifact_root: Path, policy: dict[str, Any]) -> list[str]:
    if not DISPOSITIONS.exists():
        return [f"missing artifact disposition manifest: {DISPOSITIONS}"]
    manifest = load_json(DISPOSITIONS)
    if manifest.get("schema_version") != 1:
        return ["artifact disposition manifest schema_version must be 1"]
    allowed = set(policy["allowed_retention"])
    errors: list[str] = []
    archives = manifest.get("archives")
    if not isinstance(archives, list) or not archives:
        errors.append("artifact disposition manifest archives must be non-empty")
    else:
        for archive in archives:
            if not isinstance(archive, dict):
                errors.append("artifact disposition archive entries must be objects")
                continue
            relative = archive.get("path")
            path = artifact_root / relative if isinstance(relative, str) else None
            if path is None or not path.is_file():
                errors.append(f"missing supplemental artifact archive: {relative}")
            elif archive.get("sha256") != sha256(path):
                errors.append(f"supplemental artifact archive hash mismatch: {path}")
            if not isinstance(archive.get("entries"), int) or archive["entries"] < 1:
                errors.append(f"invalid supplemental archive entry count: {relative}")
    dispositions = manifest.get("dispositions")
    if not isinstance(dispositions, list) or not dispositions:
        errors.append("artifact disposition manifest dispositions must be non-empty")
    else:
        for item in dispositions:
            if not isinstance(item, dict):
                errors.append("artifact disposition entries must be objects")
                continue
            if not isinstance(item.get("archive_member"), str) or not item["archive_member"].startswith("wheels/"):
                errors.append("artifact disposition entries require archive_member under wheels/")
            if not isinstance(item.get("sha256"), str) or len(item["sha256"]) != 64:
                errors.append("artifact disposition entries require a 64-hex sha256")
            if item.get("retention") not in allowed:
                errors.append(f"invalid artifact disposition retention: {item.get('retention')}")
    return errors


def aliases(artifact_root: Path, policy: dict[str, Any]) -> list[Path]:
    active_root = artifact_root / policy["active_alias_directory"]
    return [path for path in all_wheels(artifact_root) if active_root in path.parents]


def build_index(artifact_root: Path, policy: dict[str, Any]) -> dict[str, Any]:
    current = set(active_subjects(policy).values())
    grouped: dict[str, list[Path]] = {}
    for path in aliases(artifact_root, policy):
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
        "schema_version": 2,
        "policy": "ePC-SAFT/ePC-SAFT-management:manifests/artifact-store-policy.json",
        "active_alias_directory": policy["active_alias_directory"],
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
    errors.extend(disposition_errors(artifact_root, policy))
    if recorded != expected:
        errors.append("artifact index is stale; run artifact_store.py normalize")
    if tracked != expected:
        errors.append("tracked artifact index is stale; run artifact_store.py normalize")
    active_root = artifact_root / policy["active_alias_directory"]
    archive_root = artifact_root / policy["archive_directory"]
    for wheel in all_wheels(artifact_root):
        if archive_root in wheel.parents:
            errors.append(f"historical wheel is not compressed into the archive: {wheel}")
        elif active_root not in wheel.parents:
            errors.append(f"wheel outside the active channel: {wheel}")
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
    subjects = active_subjects(policy)
    current = set(subjects.values())
    indexed = {item["sha256"] for item in recorded.get("objects", [])}
    missing_current = sorted(current - indexed)
    if missing_current:
        errors.append(f"missing current artifact subjects: {missing_current}")
    by_distribution = {
        item.get("distribution", ""): item
        for item in recorded.get("objects", [])
        if item.get("distribution")
    }
    if set(by_distribution) != set(subjects):
        errors.append(
            "active channel must contain exactly one wheel for each distribution: "
            f"expected {sorted(subjects)}, found {sorted(by_distribution)}"
        )
    for distribution, subject in subjects.items():
        item = by_distribution.get(distribution)
        if item is None:
            continue
        if item.get("sha256") != subject:
            errors.append(
                f"active subject hash mismatch for {distribution}: "
                f"{item.get('sha256')} != {subject}"
            )
        expected_version = policy["active_distributions"][distribution].get("version")
        if expected_version and item.get("version") != expected_version:
            errors.append(
                f"active subject version mismatch for {distribution}: "
                f"{item.get('version')} != {expected_version}"
            )
        if len(item.get("aliases", [])) != 1:
            errors.append(f"active distribution has multiple aliases: {distribution}")
    archive_manifest = ROOT / "manifests" / "artifact-archive.json"
    if not archive_manifest.exists():
        errors.append(f"missing historical archive manifest: {archive_manifest}")
    historical = policy.get("historical_archive", {})
    archive_path = artifact_root / historical.get("path", "")
    if not archive_path.is_file():
        errors.append(f"missing historical wheel archive: {archive_path}")
    elif historical.get("sha256") and sha256(archive_path) != historical["sha256"]:
        errors.append(f"historical wheel archive hash mismatch: {archive_path}")
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

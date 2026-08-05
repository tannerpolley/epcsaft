#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from email.parser import Parser
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
MANAGEMENT = ROOT / "manifests" / "management.json"
POLICY = ROOT / "manifests" / "artifact-store-policy.json"
TRACKED_INDEX = ROOT / "manifests" / "artifacts.json"
DISPOSITIONS = ROOT / "manifests" / "artifact-dispositions.json"
HEX64 = set("0123456789abcdef")


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
    project_root = ROOT.parent
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


def archive_members(path: Path) -> tuple[set[str], str | None]:
    result = subprocess.run(
        ["tar", "--zstd", "-tf", str(path)],
        capture_output=True,
        text=True,
    )
    if result.returncode:
        return set(), result.stderr.strip() or "tar listing failed"
    return set(result.stdout.splitlines()), None


def archive_member_hashes(path: Path, members: set[str]) -> tuple[dict[str, str], str | None]:
    with tempfile.TemporaryDirectory(prefix="artifact-disposition-") as temp_name:
        root = Path(temp_name).resolve()
        result = subprocess.run(
            ["tar", "--zstd", "-xf", str(path), "-C", str(root), "--", *sorted(members)],
            capture_output=True,
            text=True,
        )
        if result.returncode:
            return {}, result.stderr.strip() or "tar extraction failed"
        hashes: dict[str, str] = {}
        for member in members:
            target = root / member
            if target.is_symlink() or root not in target.resolve().parents or not target.is_file():
                return {}, f"unsafe extracted archive member: {member}"
            hashes[member] = sha256(target)
    return hashes, None


def disposition_errors(
    artifact_root: Path,
    policy: dict[str, Any],
    dispositions_path: Path = DISPOSITIONS,
) -> list[str]:
    if not dispositions_path.exists():
        return [f"missing artifact disposition manifest: {dispositions_path}"]
    manifest = load_json(dispositions_path)
    if manifest.get("schema_version") != 1:
        return ["artifact disposition manifest schema_version must be 1"]
    allowed = set(policy["allowed_retention"])
    errors: list[str] = []
    archives = manifest.get("archives")
    archive_paths: dict[Path, set[str]] = {}
    if not isinstance(archives, list) or not archives:
        errors.append("artifact disposition manifest archives must be non-empty")
    else:
        for archive in archives:
            if not isinstance(archive, dict):
                errors.append("artifact disposition archive entries must be objects")
                continue
            relative = archive.get("path")
            path = (artifact_root / relative).resolve() if isinstance(relative, str) else None
            if path is None or artifact_root.resolve() not in path.parents or not path.is_file():
                errors.append(f"missing supplemental artifact archive: {relative}")
            else:
                if archive.get("sha256") != sha256(path):
                    errors.append(f"supplemental artifact archive hash mismatch: {path}")
                members, error = archive_members(path)
                if error:
                    errors.append(f"cannot read supplemental artifact archive: {path}: {error}")
                else:
                    archive_paths[path] = members
            if not isinstance(archive.get("entries"), int) or archive["entries"] < 1:
                errors.append(f"invalid supplemental archive entry count: {relative}")
    dispositions = manifest.get("dispositions")
    valid_items: list[dict[str, Any]] = []
    if not isinstance(dispositions, list) or not dispositions:
        errors.append("artifact disposition manifest dispositions must be non-empty")
    else:
        for item in dispositions:
            if not isinstance(item, dict):
                errors.append("artifact disposition entries must be objects")
                continue
            member = item.get("archive_member")
            if not isinstance(member, str) or not member.startswith("wheels/") or ".." in Path(member).parts:
                errors.append("artifact disposition entries require archive_member under wheels/")
                continue
            if not isinstance(item.get("sha256"), str) or len(item["sha256"]) != 64 or set(item["sha256"]) - HEX64:
                errors.append("artifact disposition entries require a 64-hex sha256")
            if item.get("retention") not in allowed:
                errors.append(f"invalid artifact disposition retention: {item.get('retention')}")
            valid_items.append(item)
    member_hashes: dict[tuple[Path, str], str] = {}
    for path, members in archive_paths.items():
        requested = {item["archive_member"] for item in valid_items if item["archive_member"] in members}
        if not requested:
            continue
        hashes, error = archive_member_hashes(path, requested)
        if error:
            errors.append(f"cannot read supplemental artifact members: {path}: {error}")
        member_hashes.update({(path, member): digest for member, digest in hashes.items()})
    for item in valid_items:
        member = item["archive_member"]
        matches = [path for path, members in archive_paths.items() if member in members]
        if len(matches) != 1:
            errors.append(f"artifact disposition member is not uniquely present in a declared archive: {member}")
        elif isinstance(item.get("sha256"), str) and len(item["sha256"]) == 64:
            digest = member_hashes.get((matches[0], member))
            if digest is None:
                errors.append(f"cannot hash artifact disposition member: {member}")
            elif digest != item["sha256"]:
                errors.append(f"artifact disposition member hash mismatch: {member}")
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
        "policy": "tannerpolley/ePC-SAFT:governance/manifests/artifact-store-policy.json",
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
    management = load_json(MANAGEMENT)
    if management.get("artifact_dispositions") != "manifests/artifact-dispositions.json":
        return ["management artifact_dispositions must bind the canonical manifest"]
    dispositions_path = ROOT / management.get("artifact_dispositions", "")
    index_path = artifact_root / policy["index"]
    if not index_path.exists():
        return [f"missing artifact index: {index_path}"]
    recorded = load_json(index_path)
    if not TRACKED_INDEX.exists():
        return [f"missing tracked artifact index: {TRACKED_INDEX}"]
    tracked = load_json(TRACKED_INDEX)
    expected = build_index(artifact_root, policy)
    errors: list[str] = []
    errors.extend(disposition_errors(artifact_root, policy, dispositions_path))
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


def resolve(distribution: str) -> tuple[Path | None, list[str]]:
    _, artifact_root, policy = paths()
    subject = policy.get("active_distributions", {}).get(distribution)
    if not isinstance(subject, dict):
        return None, [f"distribution is not active: {distribution}"]

    errors = check()
    if errors:
        return None, errors

    index = load_json(artifact_root / policy["index"])
    matches = [
        item
        for item in index.get("objects", [])
        if item.get("distribution") == distribution
    ]
    if len(matches) != 1:
        return None, [
            f"active distribution must resolve to exactly one object: {distribution}"
        ]

    item = matches[0]
    object_path = (artifact_root / item["object"]).resolve()
    expected_sha256 = subject.get("sha256")
    if not object_path.is_file():
        errors.append(f"missing active object: {object_path}")
    elif sha256(object_path) != expected_sha256:
        errors.append(
            f"active object hash mismatch for {distribution}: {object_path}"
        )
    return (None if errors else object_path), errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("normalize", "check", "resolve"))
    parser.add_argument("--distribution")
    arguments = parser.parse_args()
    resolved: Path | None = None
    if arguments.action == "normalize":
        errors = normalize()
    elif arguments.action == "check":
        errors = check()
    elif not arguments.distribution:
        errors = ["resolve requires --distribution"]
    else:
        resolved, errors = resolve(arguments.distribution)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    if resolved is not None:
        print(resolved)
        return 0
    print(f"artifact_store: {arguments.action} clean")
    return 0


if __name__ == "__main__":
    sys.exit(main())

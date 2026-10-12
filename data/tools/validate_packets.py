#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import subprocess
import sys
import tomllib
from pathlib import Path, PurePosixPath
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
ALLOWED_ROOTS = {"parameters", "observations", "sources", "provenance"}
ALLOWED_ROLES = {
    "parameters",
    "observations",
    "source",
    "provenance",
    "schema",
}
SHA256_LENGTH = 64
SUPPORTED_PACKET_SCHEMA_VERSIONS = {1, 2}
LEGACY_SCHEMA1_PACKET_VERSIONS = {
    "khudaida-2026-figure-2-electrolyte-lle": {"2", "3"},
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def packet_fingerprint(packet_root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(item for item in packet_root.rglob("*") if item.is_file()):
        relative = path.relative_to(packet_root).as_posix().encode("utf-8")
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(path.stat().st_size.to_bytes(8, "big"))
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
    return digest.hexdigest()


def _nonempty_string(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _safe_relative_path(value: Any) -> PurePosixPath | None:
    if not _nonempty_string(value):
        return None
    path = PurePosixPath(value)
    if (
        not path.parts
        or path.is_absolute()
        or ".." in path.parts
        or path.parts[0] not in ALLOWED_ROOTS
    ):
        return None
    return path


def _contains_symlink(packet_root: Path, relative: PurePosixPath) -> bool:
    current = packet_root
    for part in relative.parts:
        current /= part
        if current.is_symlink():
            return True
    return False


def _valid_refs(value: Any, known: set[str]) -> bool:
    return (
        isinstance(value, list)
        and bool(value)
        and all(_nonempty_string(item) for item in value)
        and set(value) <= known
    )


def _valid_hex(value: Any, length: int) -> bool:
    return (
        isinstance(value, str)
        and len(value) == length
        and all(character in "0123456789abcdef" for character in value)
    )


def _safe_bundle_member(value: Any) -> PurePosixPath | None:
    if not _nonempty_string(value):
        return None
    path = PurePosixPath(value)
    if not path.parts or path.is_absolute() or ".." in path.parts:
        return None
    return path


def _validate_parameter_bundle(
    packet_root: Path,
    relative: PurePosixPath,
    expected_paths: set[PurePosixPath],
    prefix: str,
) -> list[str]:
    path = packet_root / Path(*relative.parts)
    try:
        data = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, tomllib.TOMLDecodeError) as exc:
        return [f"{prefix}: cannot read parameter bundle {relative}: {exc}"]
    errors: list[str] = []
    if data.get("schema_version") != 1:
        errors.append(f"{prefix}: parameter bundle requires schema_version 1")
    for field in ("bundle_id", "purpose"):
        if not _nonempty_string(data.get(field)):
            errors.append(f"{prefix}: parameter bundle requires {field}")
    if not isinstance(data.get("bundle_version"), int):
        errors.append(f"{prefix}: parameter bundle requires integer bundle_version")
    files = data.get("files")
    if not isinstance(files, list) or not files:
        return errors + [f"{prefix}: parameter bundle requires files"]
    declared: set[PurePosixPath] = set()
    for index, value in enumerate(files):
        member = _safe_bundle_member(value)
        retained = None if member is None else PurePosixPath(relative.parent, member)
        if retained is None or retained in declared:
            errors.append(f"{prefix}: parameter bundle file {index} is invalid")
            continue
        declared.add(retained)
    if declared != expected_paths:
        errors.append(f"{prefix}: parameter bundle file set is incomplete")
    return errors


def _validate_migration_provenance(
    packet_root: Path,
    relative: PurePosixPath,
    expected_paths: set[PurePosixPath],
    prefix: str,
) -> list[str]:
    path = packet_root / Path(*relative.parts)
    try:
        data = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, tomllib.TOMLDecodeError) as exc:
        return [f"{prefix}: cannot read migration provenance {relative}: {exc}"]
    if data.get("provenance_kind") == "source-record-provenance-v1":
        return _validate_source_record_provenance(
            packet_root, relative, expected_paths, prefix, data
        )
    if data.get("provenance_kind") == "derived-record-provenance-v1":
        return _validate_derived_record_provenance(
            packet_root, expected_paths, prefix, data
        )
    if data.get("provenance_kind") != "repository-migration-v1":
        return [f"{prefix}: unsupported provenance_kind for {relative}"]
    errors: list[str] = []
    for field in ("origin_repository", "copy_mode", "source_bundle_manifest_path"):
        if not _nonempty_string(data.get(field)):
            errors.append(f"{prefix}: migration provenance requires {field}")
    for field in ("origin_commit", "origin_tree"):
        if not _valid_hex(data.get(field), 40):
            errors.append(f"{prefix}: migration provenance has invalid {field}")
    if not _valid_hex(data.get("source_bundle_manifest_sha256"), SHA256_LENGTH):
        errors.append(
            f"{prefix}: migration provenance has invalid source manifest hash"
        )
    records = data.get("files")
    if not isinstance(records, list) or not records:
        return errors + [f"{prefix}: migration provenance requires files"]
    mapped: set[PurePosixPath] = set()
    for index, record in enumerate(records):
        if not isinstance(record, dict):
            errors.append(f"{prefix}: provenance file {index} must be a table")
            continue
        retained = _safe_relative_path(record.get("retained_path"))
        if retained is None or retained in mapped:
            errors.append(
                f"{prefix}: provenance file {index} has invalid retained_path"
            )
            continue
        mapped.add(retained)
        if not _nonempty_string(record.get("origin_path")):
            errors.append(f"{prefix}: provenance file {index} requires origin_path")
        declared_hash = record.get("sha256")
        target = packet_root / Path(*retained.parts)
        if not _valid_hex(declared_hash, SHA256_LENGTH):
            errors.append(f"{prefix}: provenance file {retained} has invalid sha256")
        elif not target.is_file() or sha256(target) != declared_hash:
            errors.append(f"{prefix}: provenance hash mismatch for {retained}")
    if mapped != expected_paths:
        errors.append(f"{prefix}: migration provenance file set is incomplete")
    return errors


def _validate_source_record_provenance(
    packet_root: Path,
    relative: PurePosixPath,
    expected_paths: set[PurePosixPath],
    prefix: str,
    data: dict[str, Any],
) -> list[str]:
    errors: list[str] = []
    if data.get("schema_version") != 1:
        errors.append(f"{prefix}: source-record provenance requires schema_version 1")
    records = data.get("records")
    if not isinstance(records, list) or not records:
        return errors + [f"{prefix}: source-record provenance requires records"]
    mapped: set[PurePosixPath] = set()
    for index, record in enumerate(records):
        if not isinstance(record, dict):
            errors.append(f"{prefix}: source record {index} must be a table")
            continue
        retained = _safe_relative_path(record.get("retained_path"))
        if retained is None or retained in mapped:
            errors.append(f"{prefix}: source record {index} has invalid retained_path")
            continue
        mapped.add(retained)
        for field in ("origin_repository", "origin_path"):
            if not _nonempty_string(record.get(field)):
                errors.append(f"{prefix}: source record {retained} requires {field}")
        for field in ("origin_commit", "origin_tree"):
            if not _valid_hex(record.get(field), 40):
                errors.append(f"{prefix}: source record {retained} has invalid {field}")
        declared_hash = record.get("sha256")
        target = packet_root / Path(*retained.parts)
        if not _valid_hex(declared_hash, SHA256_LENGTH):
            errors.append(f"{prefix}: source record {retained} has invalid sha256")
        elif not target.is_file() or sha256(target) != declared_hash:
            errors.append(f"{prefix}: source record hash mismatch for {retained}")
    if mapped != expected_paths:
        errors.append(f"{prefix}: source-record provenance file set is incomplete")
    return errors


def _validate_derived_record_provenance(
    packet_root: Path,
    expected_paths: set[PurePosixPath],
    prefix: str,
    data: dict[str, Any],
) -> list[str]:
    errors: list[str] = []
    if data.get("schema_version") != 1:
        errors.append(f"{prefix}: derived-record provenance requires schema_version 1")
    transformations = data.get("transformations")
    if not (
        isinstance(transformations, list)
        and transformations
        and all(_nonempty_string(item) for item in transformations)
    ):
        errors.append(f"{prefix}: derived-record provenance requires transformations")
    records = data.get("records")
    if not isinstance(records, list) or not records:
        return errors + [f"{prefix}: derived-record provenance requires records"]
    mapped: set[PurePosixPath] = set()
    for index, record in enumerate(records):
        retained = (
            _safe_relative_path(record.get("retained_path"))
            if isinstance(record, dict)
            else None
        )
        if retained is None or retained in mapped:
            errors.append(f"{prefix}: derived record {index} has invalid retained_path")
            continue
        mapped.add(retained)
        declared_hash = record.get("sha256")
        target = packet_root / Path(*retained.parts)
        if not _valid_hex(declared_hash, SHA256_LENGTH):
            errors.append(f"{prefix}: derived record {retained} has invalid sha256")
        elif not target.is_file() or sha256(target) != declared_hash:
            errors.append(f"{prefix}: derived record hash mismatch for {retained}")
    if mapped != expected_paths:
        errors.append(f"{prefix}: derived-record provenance file set is incomplete")
    sources = data.get("sources")
    if not isinstance(sources, list) or not sources:
        return errors + [f"{prefix}: derived-record provenance requires sources"]
    repository_root = ROOT.parent
    source_ids: set[str] = set()
    for index, source in enumerate(sources):
        if not isinstance(source, dict):
            errors.append(f"{prefix}: derived source {index} must be a table")
            continue
        source_id = source.get("source_id")
        retained = _safe_bundle_member(source.get("retained_path"))
        packet_id = source.get("packet_id")
        packet_version = source.get("packet_version")
        packet_member = _safe_relative_path(source.get("packet_member"))
        repository_source = retained is not None and retained.parts[0] in {
            "data",
            "analyses",
        }
        packet_source = (
            _nonempty_string(packet_id)
            and "/" not in packet_id
            and _nonempty_string(packet_version)
            and "/" not in packet_version
            and packet_member is not None
        )
        if (
            not _nonempty_string(source_id)
            or source_id in source_ids
            or repository_source == packet_source
        ):
            errors.append(f"{prefix}: derived source {index} is invalid")
            continue
        source_ids.add(source_id)
        if not _nonempty_string(source.get("role")):
            errors.append(f"{prefix}: derived source {source_id} requires role")
        declared_hash = source.get("sha256")
        target = (
            repository_root / Path(*retained.parts)
            if repository_source
            else ROOT
            / "packets"
            / str(packet_id)
            / str(packet_version)
            / Path(*packet_member.parts)
        )
        if not _valid_hex(declared_hash, SHA256_LENGTH):
            errors.append(f"{prefix}: derived source {source_id} has invalid sha256")
        elif not target.is_file() or sha256(target) != declared_hash:
            errors.append(f"{prefix}: derived source hash mismatch for {source_id}")
    return errors


def validate_packet(manifest: Path) -> list[str]:
    packet_root = manifest.parent
    prefix = (
        str(packet_root.relative_to(ROOT))
        if packet_root.is_relative_to(ROOT)
        else str(packet_root)
    )
    errors: list[str] = []
    if any(
        path.is_symlink()
        for path in (
            manifest,
            packet_root,
            packet_root.parent,
            packet_root.parent.parent,
        )
    ):
        return [f"{prefix}: packet path and manifest must not be symlinks"]
    try:
        data = tomllib.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, tomllib.TOMLDecodeError) as exc:
        return [f"{prefix}: cannot read packet.toml: {exc}"]

    schema_version = data.get("schema_version")
    if schema_version not in SUPPORTED_PACKET_SCHEMA_VERSIONS:
        errors.append(
            f"{prefix}: schema_version must be one of "
            f"{sorted(SUPPORTED_PACKET_SCHEMA_VERSIONS)}"
        )
    elif (
        schema_version == 1
        and packet_root.parent.name not in LEGACY_SCHEMA1_PACKET_VERSIONS
    ):
        errors.append(
            f"{prefix}: schema_version 1 is reserved for accepted historical packets"
        )
    elif (
        schema_version == 1
        and packet_root.name
        not in LEGACY_SCHEMA1_PACKET_VERSIONS[packet_root.parent.name]
    ):
        errors.append(
            f"{prefix}: schema_version 1 is not accepted for this packet version"
        )
    if data.get("packet_id") != packet_root.parent.name:
        errors.append(f"{prefix}: packet_id must match its directory")
    if str(data.get("packet_version", "")) != packet_root.name:
        errors.append(f"{prefix}: packet_version must match its directory")
    for field in ("title", "description", "license_or_use_basis"):
        if not _nonempty_string(data.get(field)):
            errors.append(f"{prefix}: {field} is required")
    if schema_version == 2 and not _nonempty_string(data.get("tolerance_rationale")):
        errors.append(f"{prefix}: tolerance_rationale is required for schema_version 2")

    sources = data.get("sources")
    domains = data.get("domains")
    files = data.get("files")
    if not isinstance(sources, list) or not sources:
        errors.append(f"{prefix}: at least one source is required")
        sources = []
    if not isinstance(domains, list) or not domains:
        errors.append(f"{prefix}: at least one domain is required")
        domains = []
    if not isinstance(files, list) or not files:
        errors.append(f"{prefix}: at least one file is required")
        files = []

    source_ids: set[str] = set()
    for index, source in enumerate(sources):
        if not isinstance(source, dict):
            errors.append(f"{prefix}: source {index} must be a table")
            continue
        if not all(
            _nonempty_string(source.get(field))
            for field in ("source_id", "citation", "use_basis")
        ):
            errors.append(f"{prefix}: source {index} is incomplete")
            continue
        source_id = source["source_id"]
        if source_id in source_ids:
            errors.append(f"{prefix}: duplicate source_id {source_id}")
        source_ids.add(source_id)

    domain_ids: set[str] = set()
    for index, domain in enumerate(domains):
        if not isinstance(domain, dict) or not all(
            _nonempty_string(domain.get(field))
            for field in ("domain_id", "description")
        ):
            errors.append(f"{prefix}: domain {index} is incomplete")
            continue
        domain_id = domain["domain_id"]
        if domain_id in domain_ids:
            errors.append(f"{prefix}: duplicate domain_id {domain_id}")
        domain_ids.add(domain_id)

    listed: set[PurePosixPath] = set()
    parameter_bundle_files: set[PurePosixPath] = set()
    parameter_payloads: set[PurePosixPath] = set()
    provenance_files: set[PurePosixPath] = set()
    provenance_inputs: set[PurePosixPath] = set()
    for index, record in enumerate(files):
        if not isinstance(record, dict):
            errors.append(f"{prefix}: file {index} must be a table")
            continue
        relative = _safe_relative_path(record.get("path"))
        if relative is None:
            errors.append(f"{prefix}: file {index} has unsafe path")
            continue
        if relative in listed:
            errors.append(f"{prefix}: duplicate file path {relative}")
            continue
        listed.add(relative)
        role = record.get("role")
        if role not in ALLOWED_ROLES:
            errors.append(f"{prefix}: {relative} has invalid role")
        elif role == "provenance":
            provenance_files.add(relative)
        else:
            provenance_inputs.add(relative)
            if role == "parameters":
                if relative.name == "bundle.toml":
                    parameter_bundle_files.add(relative)
                else:
                    parameter_payloads.add(relative)
        for field in ("media_type", "units_basis"):
            if not _nonempty_string(record.get(field)):
                errors.append(f"{prefix}: {relative} requires {field}")
        declared_hash = record.get("sha256")
        if (
            not isinstance(declared_hash, str)
            or len(declared_hash) != SHA256_LENGTH
            or any(character not in "0123456789abcdef" for character in declared_hash)
        ):
            errors.append(f"{prefix}: {relative} has invalid sha256")
        target = packet_root / Path(*relative.parts)
        if not target.is_file() or _contains_symlink(packet_root, relative):
            errors.append(f"{prefix}: listed file is missing or a symlink: {relative}")
        elif declared_hash != sha256(target):
            errors.append(f"{prefix}: hash mismatch for {relative}")
        source_refs = record.get("source_refs")
        domain_refs = record.get("domain_refs")
        if not _valid_refs(source_refs, source_ids):
            errors.append(f"{prefix}: {relative} has invalid source_refs")
        if not _valid_refs(domain_refs, domain_ids):
            errors.append(f"{prefix}: {relative} has invalid domain_refs")

    actual = {
        PurePosixPath(path.relative_to(packet_root).as_posix())
        for path in packet_root.rglob("*")
        if path.is_file() and path != manifest
    }
    if actual != listed:
        unlisted = sorted(str(path) for path in actual - listed)
        missing = sorted(str(path) for path in listed - actual)
        if unlisted:
            errors.append(f"{prefix}: unlisted files: {unlisted}")
        if missing:
            errors.append(f"{prefix}: manifest-only files: {missing}")
    empty_directories = [
        path.relative_to(packet_root).as_posix()
        for path in packet_root.rglob("*")
        if path.is_dir() and not path.is_symlink() and not any(path.iterdir())
    ]
    if empty_directories:
        errors.append(f"{prefix}: empty directory: {sorted(empty_directories)}")
    unlisted_symlinks = [
        path.relative_to(packet_root).as_posix()
        for path in packet_root.rglob("*")
        if path.is_symlink()
        and PurePosixPath(path.relative_to(packet_root).as_posix()) not in listed
    ]
    if unlisted_symlinks:
        errors.append(f"{prefix}: unlisted symlinks: {sorted(unlisted_symlinks)}")
    for relative in parameter_bundle_files:
        errors.extend(
            _validate_parameter_bundle(
                packet_root,
                relative,
                {path for path in parameter_payloads if path.parent == relative.parent},
                prefix,
            )
        )
    for relative in provenance_files:
        errors.extend(
            _validate_migration_provenance(
                packet_root,
                relative,
                provenance_inputs,
                prefix,
            )
        )
    return errors


def validate_repository(root: Path = ROOT) -> list[str]:
    packets = root / "packets"
    if not packets.exists():
        return []
    manifests = sorted(packets.glob("*/*/packet.toml"))
    errors: list[str] = []
    empty_packet_directories = [
        path.relative_to(root).as_posix()
        for path in packets.iterdir()
        if path.is_dir() and not path.is_symlink() and not any(path.iterdir())
    ]
    if empty_packet_directories:
        errors.append(
            f"packets: empty packet directory: {sorted(empty_packet_directories)}"
        )
    version_directories = [
        path for path in packets.glob("*/*") if path.is_dir() or path.is_symlink()
    ]
    for directory in version_directories:
        if not (directory / "packet.toml").is_file():
            errors.append(f"{directory.relative_to(root)}: missing packet.toml")
    for manifest in manifests:
        errors.extend(validate_packet(manifest))
    stray = [
        path
        for path in packets.rglob("*")
        if path.is_file()
        and not any(
            manifest.parent in path.parents or path == manifest
            for manifest in manifests
        )
    ]
    if stray:
        relative = [str(path.relative_to(root)) for path in stray]
        errors.append(f"packets: files outside a packet version: {relative}")
    return errors


def main() -> int:
    errors = validate_repository()
    index = subprocess.run([sys.executable, ROOT / "tools" / "parameter_index.py"], capture_output=True, text=True)
    if index.returncode or index.stdout != (ROOT / "docs" / "parameter-index.md").read_text():
        errors.append(f"data/docs/parameter-index.md is stale: rerun data/tools/parameter_index.py {index.stderr}".strip())
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        return 1
    print("packet_validation: clean")
    for manifest in sorted((ROOT / "packets").glob("*/*/packet.toml")):
        relative = manifest.parent.relative_to(ROOT / "packets")
        print(f"packet {relative} sha256:{packet_fingerprint(manifest.parent)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

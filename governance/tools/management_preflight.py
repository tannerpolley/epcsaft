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
from urllib.parse import unquote, urlparse


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "manifests" / "management.json"
GENERATED = ROOT / ".management-state"
ARTIFACT_POLICY_PATH = ROOT / "manifests" / "artifact-store-policy.json"
OPERATION_RECEIPT_SCHEMA_PATH = (
    ROOT / "manifests" / "operation-receipt-envelope-v1.json"
)

ALLOWED_ENVIRONMENT_DISTRIBUTIONS: dict[str, set[str]] = {
    "management": set(),
    "github": set(),
    "data": set(),
    "eos": {"epcsaft"},
    "equilibrium": {"epcsaft", "epcsaft_equilibrium"},
    "regression": {"epcsaft", "epcsaft_regression"},
    "validation": set(),
    "lab": set(),
}


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


def json_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode(
        "utf-8"
    )
    return hashlib.sha256(payload).hexdigest()


def git_blob_sha256(root: Path, commit: str, path: str) -> str:
    completed = subprocess.run(
        ["git", "-C", str(root), "show", f"{commit}:{path}"],
        check=True,
        capture_output=True,
    )
    return hashlib.sha256(completed.stdout).hexdigest()


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
    retention_by_path: dict[str, str] = {}
    index_path = artifact_root / "index.json"
    if index_path.exists():
        index = load_json(index_path)
        for item in index.get("objects", []):
            for alias in item.get("aliases", []):
                retention_by_path[
                    (Path("artifacts") / alias).as_posix()
                ] = item["retention"]
    if artifact_root.exists():
        candidates.update(
            path
            for path in artifact_root.rglob("*.whl")
            if artifact_root / "objects" not in path.parents
        )
    for record in config["repositories"].values():
        local = Path(record["local_path"])
        if local.exists():
            candidates.update((local / "dist").glob("*.whl"))
    wheels = []
    for path in sorted(candidates):
        if not path.is_file():
            continue
        metadata = wheel_metadata(path, project_root)
        metadata["retention"] = retention_by_path.get(
            metadata["path"], "unclassified"
        )
        wheels.append(metadata)
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


def ownership_projection(policy: dict[str, Any]) -> dict[str, Any]:
    return {
        "storage_owner": policy.get("storage_owner", "unclassified"),
        "runtime_owner": policy.get("runtime_owner"),
        "evidence_owner": policy.get("evidence_owner"),
        "disposition": policy.get("disposition", "unclassified"),
    }


def parameter_state(config: dict[str, Any]) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    parameter_policy = load_json(parameter_policy_path(config))
    parameter_entries = {
        entry["bundle_id"]: entry for entry in parameter_policy.get("bundles", [])
    }
    data_root = Path(config["repositories"]["data"]["local_path"])
    packets = data_root / "packets"
    if packets.exists():
        for manifest in sorted(packets.glob("*/*/packet.toml")):
            data = tomllib.loads(manifest.read_text(encoding="utf-8"))
            directory = manifest.parent
            bundle_id = data.get(
                "packet_id", data.get("bundle_id", directory.parent.name)
            )
            policy = parameter_entries.get(bundle_id, {})
            records.append(
                {
                    "repository": "data",
                    "path": str(directory.relative_to(data_root)),
                    "bundle_id": bundle_id,
                    "sha256": bundle_hash(directory),
                    **ownership_projection(policy),
                }
            )
    eos = Path(config["repositories"]["eos"]["local_path"])
    eos_bundles = eos / "src" / "epcsaft" / "data" / "bundles"
    if eos_bundles.exists():
        for manifest in sorted(eos_bundles.glob("*/*/bundle.toml")):
            data = tomllib.loads(manifest.read_text(encoding="utf-8"))
            directory = manifest.parent
            bundle_id = data.get("bundle_id", directory.parent.name)
            policy = parameter_entries.get(bundle_id, {})
            records.append(
                {
                    "repository": "eos",
                    "path": str(directory.relative_to(eos)),
                    "bundle_id": bundle_id,
                    "sha256": bundle_hash(directory),
                    **ownership_projection(policy),
                }
            )
    validation = Path(config["repositories"]["validation"]["local_path"])
    if validation.exists():
        for manifest in sorted(validation.rglob("bundle.toml")):
            if ".git" in manifest.parts or ".worktrees" in manifest.parts:
                continue
            data = tomllib.loads(manifest.read_text(encoding="utf-8"))
            directory = manifest.parent
            bundle_id = data.get(
                "bundle_id", data.get("packet_id", directory.parent.name)
            )
            policy = parameter_entries.get(bundle_id, {})
            records.append(
                {
                    "repository": "validation",
                    "path": str(directory.relative_to(validation)),
                    "bundle_id": bundle_id,
                    "sha256": bundle_hash(directory),
                    "migration_state": "legacy_copy_pending_disposition",
                    **ownership_projection(policy),
                }
            )
    return {"schema_version": 2, "records": records}


def normalized_distribution_name(name: str) -> str:
    return name.lower().replace("-", "_")


def environment_distributions(environment: Path) -> list[dict[str, Any]]:
    distributions: list[dict[str, Any]] = []
    if not environment.exists():
        return distributions
    for metadata in sorted(
        environment.glob("lib/python*/site-packages/*.dist-info/METADATA")
    ):
        message = Parser().parsestr(
            metadata.read_text(encoding="utf-8", errors="replace")
        )
        distribution = message.get("Name", "")
        if "pcsaft" not in normalized_distribution_name(distribution):
            continue
        direct_url_path = metadata.parent / "direct_url.json"
        direct_url = load_json(direct_url_path) if direct_url_path.exists() else None
        distributions.append(
            {
                "name": distribution,
                "version": message.get("Version", ""),
                "direct_url": direct_url,
            }
        )
    return distributions


def environment_state(config: dict[str, Any]) -> dict[str, Any]:
    environments: list[dict[str, Any]] = []
    for name, record in sorted(config["repositories"].items()):
        root = Path(record["local_path"])
        environment = root / ".venv"
        environments.append(
            {
                "repository": name,
                "path": str(environment),
                "present": environment.exists(),
                "distributions": environment_distributions(environment),
            }
        )
    return {"schema_version": 1, "environments": environments}


def environment_policy_errors(config: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    artifact_policy = load_json(ARTIFACT_POLICY_PATH)
    active_eos_sha = artifact_policy["active_distributions"]["epcsaft"]["sha256"]
    for name, record in sorted(config["repositories"].items()):
        environment = Path(record["local_path"]) / ".venv"
        allowed = ALLOWED_ENVIRONMENT_DISTRIBUTIONS.get(name, set())
        for distribution in environment_distributions(environment):
            normalized = normalized_distribution_name(distribution["name"])
            if normalized not in allowed:
                errors.append(
                    f"{name} environment contains forbidden distribution "
                    f"{distribution['name']}"
                )
                continue
            if name not in {"equilibrium", "regression"} or normalized != "epcsaft":
                continue
            direct_url = distribution.get("direct_url")
            url = direct_url.get("url", "") if isinstance(direct_url, dict) else ""
            expected_marker = f"/artifacts/objects/sha256/{active_eos_sha}/"
            if not url.startswith("file://") or expected_marker not in url:
                errors.append(
                    f"{name} must install epcsaft from the active SHA-256 object "
                    f"{active_eos_sha}, found {url or 'no direct_url'}"
                )
                continue
            object_path = Path(unquote(urlparse(url).path))
            if not object_path.is_file():
                errors.append(
                    f"{name} epcsaft direct_url does not exist: {object_path}"
                )
            elif sha256(object_path) != active_eos_sha:
                errors.append(
                    f"{name} epcsaft direct_url hash does not match active object "
                    f"{active_eos_sha}: {object_path}"
                )
    return errors


def parameter_policy_path(config: dict[str, Any]) -> Path:
    relative = config.get("parameter_ownership", "manifests/parameter-ownership.json")
    return ROOT / relative


def parameter_policy_errors(config: dict[str, Any]) -> list[str]:
    path = parameter_policy_path(config)
    if not path.exists():
        return [f"missing parameter ownership manifest {path.relative_to(ROOT)}"]
    policy = load_json(path)
    if policy.get("schema_version") != 2:
        return ["parameter ownership manifest schema_version must be 2"]
    entries = policy.get("bundles")
    if not isinstance(entries, list):
        return ["parameter ownership manifest bundles must be an array"]
    by_id: dict[str, dict[str, Any]] = {}
    errors: list[str] = []
    for entry in entries:
        bundle_id = entry.get("bundle_id")
        if not isinstance(bundle_id, str) or not bundle_id:
            errors.append("parameter ownership entries require bundle_id")
            continue
        if bundle_id in by_id:
            errors.append(f"duplicate parameter ownership entry: {bundle_id}")
        by_id[bundle_id] = entry
        if entry.get("storage_owner") != "data":
            errors.append(f"parameter storage owner must be Data for {bundle_id}")
        if entry.get("runtime_owner") not in {"eos", None}:
            errors.append(f"invalid parameter runtime owner for {bundle_id}")
        if entry.get("evidence_owner") not in {"validation", None}:
            errors.append(f"invalid parameter evidence owner for {bundle_id}")
        disposition = entry.get("disposition")
        if disposition not in {"sync_to_eos_runtime", "relocate_to_data"}:
            errors.append(f"invalid parameter disposition for {bundle_id}")
        elif disposition == "sync_to_eos_runtime" and entry.get("runtime_owner") != "eos":
            errors.append(f"EOS runtime sync requires EOS runtime ownership for {bundle_id}")
        elif disposition == "relocate_to_data" and entry.get("runtime_owner") is not None:
            errors.append(f"Data relocation cannot retain runtime ownership for {bundle_id}")
    eos = Path(config["repositories"]["eos"]["local_path"])
    eos_bundles = eos / "src" / "epcsaft" / "data" / "bundles"
    discovered: set[str] = set()
    if eos_bundles.exists():
        for manifest in eos_bundles.glob("*/*/bundle.toml"):
            data = tomllib.loads(manifest.read_text(encoding="utf-8"))
            discovered.add(data.get("bundle_id", manifest.parent.parent.name))
    missing = sorted(discovered.difference(by_id))
    if missing:
        errors.append(f"unregistered EOS parameter bundles: {missing}")
    return errors


def validate_operation_receipt(
    relative: str,
    receipt: dict[str, Any],
    schema: dict[str, Any],
    schema_hash: str,
) -> list[str]:
    errors: list[str] = []
    envelope = receipt.get("envelope")
    record = receipt.get("record")
    if not isinstance(envelope, dict):
        return [f"operation receipt {relative} requires envelope"]
    if not isinstance(record, dict):
        return [f"operation receipt {relative} requires record"]
    missing = [
        field
        for field in schema["required_envelope_fields"]
        if field not in envelope
    ]
    missing_record = [
        field for field in schema["required_record_fields"] if field not in record
    ]
    if missing:
        errors.append(
            f"operation receipt {relative} missing envelope fields: {missing}"
        )
    if missing_record:
        errors.append(
            f"operation receipt {relative} missing record fields: {missing_record}"
        )
    if missing or missing_record:
        return errors
    declared_schema = envelope["schema"]
    expected_schema_path = str(OPERATION_RECEIPT_SCHEMA_PATH.relative_to(ROOT))
    if (
        declared_schema.get("path") != expected_schema_path
        or declared_schema.get("version") != schema["schema_version"]
        or declared_schema.get("sha256") != schema_hash
    ):
        errors.append(f"operation receipt {relative} has invalid schema identity")
    record_hash = json_sha256(record)
    subject = envelope["subject"]
    if subject.get("record_sha256") != record_hash:
        errors.append(f"operation receipt {relative} subject does not bind its record")
    subject_hash = json_sha256(subject)
    if envelope["subject_sha256"] != subject_hash:
        errors.append(f"operation receipt {relative} has invalid subject hash")
    approval = envelope["user_approval"]
    if approval.get("subject_sha256") != subject_hash:
        errors.append(
            f"operation receipt {relative} approval is not bound to its subject"
        )
    if envelope["producer"] == envelope["reviewer"]:
        errors.append(f"operation receipt {relative} producer and reviewer must differ")
    validator = envelope["validator"]
    try:
        actual_validator_hash = git_blob_sha256(
            ROOT, validator["commit"], validator["path"]
        )
    except (KeyError, subprocess.CalledProcessError):
        errors.append(f"operation receipt {relative} has invalid validator identity")
    else:
        if validator.get("sha256") != actual_validator_hash:
            errors.append(f"operation receipt {relative} has invalid validator hash")
    if record["status"] not in schema["allowed_record_statuses"]:
        errors.append(f"operation receipt {relative} has invalid record status")
    observations = record["observations"]
    if not isinstance(observations, list) or not observations:
        errors.append(f"operation receipt {relative} requires observations")
    elif any(item.get("outcome") != "pass" for item in observations):
        errors.append(f"operation receipt {relative} has a failed observation")
    return errors


def operation_receipt_errors(config: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    if not OPERATION_RECEIPT_SCHEMA_PATH.exists():
        return ["missing operation receipt envelope schema"]
    schema = load_json(OPERATION_RECEIPT_SCHEMA_PATH)
    schema_hash = sha256(OPERATION_RECEIPT_SCHEMA_PATH)
    sequences: list[int] = []
    for relative in config.get("operation_receipts", []):
        path = ROOT / relative
        if not path.exists():
            errors.append(f"missing operation receipt {relative}")
            continue
        receipt = load_json(path)
        errors.extend(
            validate_operation_receipt(relative, receipt, schema, schema_hash)
        )
        envelope = receipt.get("envelope")
        if isinstance(envelope, dict) and isinstance(envelope.get("sequence"), int):
            sequences.append(envelope["sequence"])
    if sequences != sorted(set(sequences)):
        errors.append("operation receipt sequences must be unique and monotonic")
    return errors


def validate_config(config: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    repositories = config.get("repositories")
    if not isinstance(repositories, dict):
        return ["repositories must be an object"]
    required = {
        "management",
        "github",
        "data",
        "eos",
        "equilibrium",
        "regression",
        "validation",
        "lab",
    }
    missing = required.difference(repositories)
    if missing:
        errors.append(f"missing repositories: {sorted(missing)}")
    validation_access = set(
        repositories.get("validation", {}).get("outbound_access", [])
    )
    if validation_access != {"validation", "data"}:
        errors.append("Validation outbound access must be Validation and Data only")
    data_access = repositories.get("data", {}).get("outbound_access")
    if data_access != ["data"]:
        errors.append("Data outbound access must be Data only")
    eos_access = set(repositories.get("eos", {}).get("outbound_access", []))
    if eos_access != {"eos", "data", "validation"}:
        errors.append("EOS outbound access must be EOS, Data, and Validation only")
    lab_access = set(repositories.get("lab", {}).get("outbound_access", []))
    if lab_access != {"lab", "data"}:
        errors.append("Lab outbound access must be Lab and Data only")
    expected_top = {
        "management",
        "data",
        "eos",
        "equilibrium",
        "regression",
        "validation",
        "lab",
    }
    for owner in ("management", "equilibrium", "regression"):
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
    worktree_policy = config.get("worktree_policy", {})
    if worktree_policy.get("lifecycle") != "ephemeral":
        errors.append("worktree lifecycle must be ephemeral")
    if worktree_policy.get("persistent_slots") is not False:
        errors.append("persistent worktree slots must be disabled")
    if config.get("persistent_worktrees") != []:
        errors.append("persistent worktree inventory must be empty")
    for name, record in repositories.items():
        if record.get("codex_environment") != ".codex/environments/environment.toml":
            errors.append(f"{name} must declare the standard Codex environment path")
    if config.get("parameter_ownership") != "manifests/parameter-ownership.json":
        errors.append(
            "parameter ownership manifest must be manifests/parameter-ownership.json"
        )
    if not isinstance(config.get("operation_receipts"), list):
        errors.append("operation receipts must be an array")
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
    errors.extend(parameter_policy_errors(config))
    errors.extend(operation_receipt_errors(config))
    errors.extend(environment_policy_errors(config))
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
    artifact_check = subprocess.run(
        [sys.executable, str(ROOT / "tools" / "artifact_store.py"), "check"],
        capture_output=True,
        text=True,
    )
    if artifact_check.returncode:
        errors.extend(
            line.removeprefix("ERROR: ")
            for line in artifact_check.stdout.splitlines()
            if line.startswith("ERROR: ")
        )
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

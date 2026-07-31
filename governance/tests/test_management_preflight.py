from __future__ import annotations

import importlib.util
import json
from copy import deepcopy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "management_preflight",
    ROOT / "tools" / "management_preflight.py",
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_management_configuration_is_valid() -> None:
    config = json.loads(
        (ROOT / "manifests" / "management.json").read_text(encoding="utf-8")
    )
    assert MODULE.validate_config(config) == []


def test_validation_can_read_data_but_cannot_gain_runtime_access() -> None:
    config = json.loads(
        (ROOT / "manifests" / "management.json").read_text(encoding="utf-8")
    )
    config["repositories"]["validation"]["outbound_access"].append("eos")
    assert (
        "Validation outbound access must be Validation and Data only"
        in MODULE.validate_config(config)
    )


def test_data_repository_is_registered_as_a_self_contained_asset_store() -> None:
    config = json.loads(
        (ROOT / "manifests" / "management.json").read_text(encoding="utf-8")
    )
    data = config["repositories"]["data"]
    assert data["github"] == "ePC-SAFT/ePC-SAFT-data"
    assert data["role"] == "scientific_asset"
    assert data["outbound_access"] == ["data"]


def test_worktrees_are_ephemeral_and_have_no_persistent_slots() -> None:
    config = json.loads(
        (ROOT / "manifests" / "management.json").read_text(encoding="utf-8")
    )
    assert config["persistent_worktrees"] == []
    assert config["worktree_policy"]["lifecycle"] == "ephemeral"
    assert config["worktree_policy"]["persistent_slots"] is False


def test_every_active_repository_declares_a_codex_environment() -> None:
    config = json.loads(
        (ROOT / "manifests" / "management.json").read_text(encoding="utf-8")
    )
    assert all(
        record["codex_environment"] == ".codex/environments/environment.toml"
        for record in config["repositories"].values()
    )


def test_parameter_ownership_manifest_covers_current_eos_catalog() -> None:
    config = json.loads(
        (ROOT / "manifests" / "management.json").read_text(encoding="utf-8")
    )
    assert MODULE.parameter_policy_errors(config) == []


def test_configured_operation_receipts_use_the_canonical_envelope() -> None:
    config = json.loads(
        (ROOT / "manifests" / "management.json").read_text(encoding="utf-8")
    )
    assert MODULE.operation_receipt_errors(config) == []


def test_operation_receipt_rejects_record_tampering() -> None:
    receipt = json.loads(
        (ROOT / "receipts" / "data-repository-creation-preflight.json").read_text(
            encoding="utf-8"
        )
    )
    receipt["record"]["observations"][0]["outcome"] = "fail"
    schema = json.loads(MODULE.OPERATION_RECEIPT_SCHEMA_PATH.read_text())
    errors = MODULE.validate_operation_receipt(
        "tampered.json",
        receipt,
        schema,
        MODULE.sha256(MODULE.OPERATION_RECEIPT_SCHEMA_PATH),
    )
    assert "operation receipt tampered.json subject does not bind its record" in errors
    assert "operation receipt tampered.json has a failed observation" in errors


def test_operation_receipt_rejects_same_producer_and_reviewer() -> None:
    receipt = json.loads(
        (ROOT / "receipts" / "data-repository-creation-preflight.json").read_text(
            encoding="utf-8"
        )
    )
    receipt["envelope"]["reviewer"] = deepcopy(
        receipt["envelope"]["producer"]
    )
    schema = json.loads(MODULE.OPERATION_RECEIPT_SCHEMA_PATH.read_text())
    errors = MODULE.validate_operation_receipt(
        "self-reviewed.json",
        receipt,
        schema,
        MODULE.sha256(MODULE.OPERATION_RECEIPT_SCHEMA_PATH),
    )
    assert (
        "operation receipt self-reviewed.json producer and reviewer must differ"
        in errors
    )


def test_operation_receipt_rejects_wrong_schema_version() -> None:
    receipt = json.loads(
        (ROOT / "receipts" / "data-repository-creation-preflight.json").read_text(
            encoding="utf-8"
        )
    )
    receipt["envelope"]["schema"]["version"] = 2
    schema = json.loads(MODULE.OPERATION_RECEIPT_SCHEMA_PATH.read_text())
    errors = MODULE.validate_operation_receipt(
        "wrong-schema.json",
        receipt,
        schema,
        MODULE.sha256(MODULE.OPERATION_RECEIPT_SCHEMA_PATH),
    )
    assert "operation receipt wrong-schema.json has invalid schema identity" in errors


def test_operation_receipt_rejects_malformed_nested_fields_without_crashing() -> None:
    original = json.loads(
        (ROOT / "receipts" / "data-repository-creation-preflight.json").read_text(
            encoding="utf-8"
        )
    )
    schema = json.loads(MODULE.OPERATION_RECEIPT_SCHEMA_PATH.read_text())
    schema_hash = MODULE.sha256(MODULE.OPERATION_RECEIPT_SCHEMA_PATH)
    for field in ("schema", "subject", "validator", "user_approval"):
        receipt = deepcopy(original)
        receipt["envelope"][field] = None
        errors = MODULE.validate_operation_receipt(
            f"malformed-{field}.json", receipt, schema, schema_hash
        )
        assert errors
    receipt = deepcopy(original)
    receipt["record"]["observations"] = [None]
    errors = MODULE.validate_operation_receipt(
        "malformed-observation.json", receipt, schema, schema_hash
    )
    assert errors


def test_parameter_packets_separate_storage_runtime_and_evidence_owners() -> None:
    policy = json.loads(
        (ROOT / "manifests" / "parameter-ownership.json").read_text(
            encoding="utf-8"
        )
    )
    assert policy["schema_version"] == 2
    assert all(entry["storage_owner"] == "data" for entry in policy["bundles"])
    assert all(
        entry["runtime_owner"] in {"eos", None} for entry in policy["bundles"]
    )
    assert all(
        entry["evidence_owner"] in {"validation", None}
        for entry in policy["bundles"]
    )
    assert {
        entry["disposition"] for entry in policy["bundles"]
    } <= {"relocate_to_data", "sync_to_eos_runtime"}


def test_parameter_ledger_retains_legacy_validation_copies_until_migration() -> None:
    config = json.loads(
        (ROOT / "manifests" / "management.json").read_text(encoding="utf-8")
    )
    records = MODULE.parameter_state(config)["records"]
    legacy_records = [
        record for record in records if record["repository"] == "validation"
    ]
    assert {record["bundle_id"] for record in legacy_records} == {
        "figiel-2025-validation-parameters",
        "khudaida-2026-figure-2-electrolyte-lle",
    }
    assert all(
        record["disposition"] == "relocate_to_data"
        and record["migration_state"] == "legacy_copy_pending_disposition"
        for record in legacy_records
    )


def test_environment_policy_rejects_forbidden_package(tmp_path: Path) -> None:
    site_packages = tmp_path / "eos" / ".venv" / "lib" / "python3.13" / "site-packages"
    metadata = site_packages / "epcsaft-regression-0.2.0.dist-info"
    metadata.mkdir(parents=True)
    (metadata / "METADATA").write_text(
        "Name: epcsaft-regression\nVersion: 0.2.0.dev0\n",
        encoding="utf-8",
    )
    config = {"repositories": {"eos": {"local_path": str(tmp_path / "eos")}}}
    errors = MODULE.environment_policy_errors(config)
    assert errors == ["eos environment contains forbidden distribution epcsaft-regression"]


def test_environment_policy_rejects_mutable_consumer_eos_install(tmp_path: Path) -> None:
    site_packages = (
        tmp_path
        / "equilibrium"
        / ".venv"
        / "lib"
        / "python3.13"
        / "site-packages"
    )
    metadata = site_packages / "epcsaft-0.2.0.dev0.dist-info"
    metadata.mkdir(parents=True)
    (metadata / "METADATA").write_text(
        "Name: epcsaft\nVersion: 0.2.0.dev0\n",
        encoding="utf-8",
    )
    (metadata / "direct_url.json").write_text(
        '{"url":"file:///tmp/repository/dist/epcsaft.whl"}\n',
        encoding="utf-8",
    )
    config = {
        "repositories": {"equilibrium": {"local_path": str(tmp_path / "equilibrium")}}
    }
    errors = MODULE.environment_policy_errors(config)
    assert len(errors) == 1
    assert "must install epcsaft from the active SHA-256 object" in errors[0]

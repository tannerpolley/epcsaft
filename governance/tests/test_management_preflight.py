from __future__ import annotations

import importlib.util
import json
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

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


def test_validation_cannot_gain_outbound_access() -> None:
    config = json.loads(
        (ROOT / "manifests" / "management.json").read_text(encoding="utf-8")
    )
    config["repositories"]["validation"]["outbound_access"].append("eos")
    assert "Validation outbound access must be Validation only" in MODULE.validate_config(
        config
    )


def test_persistent_worktrees_are_target_owned_access_slots() -> None:
    config = json.loads(
        (ROOT / "manifests" / "management.json").read_text(encoding="utf-8")
    )
    assert all(
        item["relative_path"].startswith(".worktrees/")
        and item["relative_path"].endswith("-access")
        for item in config["persistent_worktrees"]
    )


def test_every_active_repository_declares_a_codex_environment() -> None:
    config = json.loads(
        (ROOT / "manifests" / "management.json").read_text(encoding="utf-8")
    )
    assert all(
        record["codex_environment"] == ".codex/environments/environment.toml"
        for record in config["repositories"].values()
    )

from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "artifact_store",
    ROOT / "tools" / "artifact_store.py",
)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_active_channel_has_one_full_sha256_subject_per_distribution() -> None:
    policy = json.loads(
        (ROOT / "manifests" / "artifact-store-policy.json").read_text(
            encoding="utf-8"
        )
    )
    assert policy["rules"]["mutable_repository_dist_wheels"] == "forbidden"
    assert policy["rules"]["cross_repository_install"] == "object-path-and-sha256-only"
    assert set(policy["active_distributions"]) == {
        "epcsaft",
        "epcsaft-equilibrium",
        "epcsaft-regression",
    }
    assert all(
        len(digest) == 64 and set(digest) <= set("0123456789abcdef")
        for subject in policy["active_distributions"].values()
        for digest in [subject["sha256"]]
    )
    assert policy["rules"]["one_active_wheel_per_distribution"] == "required"


def test_candidate_wheels_have_a_compressed_disposition_record() -> None:
    management = json.loads(
        (ROOT / "manifests" / "management.json").read_text(encoding="utf-8")
    )
    errors = MODULE.disposition_errors(
        Path(management["artifact_root"]),
        json.loads(
            (ROOT / "manifests" / "artifact-store-policy.json").read_text(
                encoding="utf-8"
            )
        ),
    )
    assert errors == []

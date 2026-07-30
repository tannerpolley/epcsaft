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


def test_current_subjects_use_full_sha256_identities() -> None:
    policy = json.loads(
        (ROOT / "manifests" / "artifact-store-policy.json").read_text(
            encoding="utf-8"
        )
    )
    assert policy["rules"]["mutable_repository_dist_wheels"] == "forbidden"
    assert policy["rules"]["cross_repository_install"] == "object-path-and-sha256-only"
    assert all(
        len(digest) == 64 and set(digest) <= set("0123456789abcdef")
        for digest in policy["current_subjects"].values()
    )

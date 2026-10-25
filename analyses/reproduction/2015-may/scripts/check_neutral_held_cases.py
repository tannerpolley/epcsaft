#!/usr/bin/env python3
"""Check the frozen source-only neutral HELD validation cases."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


SOURCE_SHA256 = "5cd1e74925a3c6504f5106dcf911f2cae2d6e99a5133fccc20454d8991bdbc7f"
SOURCE_METADATA_SHA256 = "d3fbd456f0285c022c2236c3966eb5e913f711a677d09191d99d5e2629a3b2e7"
TOLERANCE_SHA256 = "ad744526678355be6ca47cf27ab9ff7ae66b7661c27e36ffe259c5b6295f1016"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check(cases_path: Path, source_path: Path) -> dict[str, Any]:
    cases_path = cases_path.resolve()
    source_path = source_path.resolve()
    if sha256_file(source_path) != SOURCE_SHA256:
        raise ValueError("source CSV SHA-256 does not match the frozen contract")
    with source_path.open(encoding="utf-8", newline="") as handle:
        source_rows = list(csv.DictReader(handle))
    source_by_id = {row["row_id"]: row for row in source_rows}
    contract = json.loads(cases_path.read_text(encoding="utf-8"))
    source = contract["source"]
    if source["csv_sha256"] != SOURCE_SHA256:
        raise ValueError("case contract source hash changed")
    if source["metadata_sha256"] != SOURCE_METADATA_SHA256:
        raise ValueError("case contract metadata hash changed")
    if source["tolerance_sha256"] != TOLERANCE_SHA256:
        raise ValueError("case contract tolerance hash changed")

    cases = contract["cases"]
    midpoint_cases = [case for case in cases if case["case_id"].endswith("-midpoint")]
    if len(cases) != 18 or len(midpoint_cases) != 17:
        raise ValueError("case contract must retain 17 midpoint cases and one derived case")
    expected_rows = [f"may2015-ch4-c2h6-{index:03d}" for index in range(1, 18)]
    if [case["source_row_id"] for case in midpoint_cases] != expected_rows:
        raise ValueError("midpoint cases must preserve all source rows in order")
    for case in midpoint_cases:
        row = source_by_id[case["source_row_id"]]
        expected_feed = 0.5 * (float(row["x_methane"]) + float(row["y_methane"]))
        if not math.isclose(case["feed_z_methane"], expected_feed, abs_tol=1e-12):
            raise ValueError(f"midpoint feed changed for {case['source_row_id']}")
        if case["expected_phase_count"] != 2:
            raise ValueError("a directly observed coexistence row lost its two-phase role")

    derived = next(case for case in cases if case["case_id"] == "may2015-row-011-liquid-side")
    row_011 = source_by_id[derived["source_row_id"]]
    expected_derived = float(row_011["x_methane"]) - 2.0 * float(
        row_011["x_comparison_allowance"]
    )
    expanded_boundary = float(row_011["x_methane"]) - float(
        row_011["x_comparison_allowance"]
    )
    domain_min, domain_max = contract["public_feed_domain_methane"]
    if not math.isclose(derived["feed_z_methane"], expected_derived, abs_tol=1e-12):
        raise ValueError("row-011 liquid-side derivation changed")
    if not domain_min <= expected_derived <= domain_max or not expected_derived < expanded_boundary:
        raise ValueError("row-011 liquid-side case is outside its frozen source/domain bounds")
    if derived["expected_phase_count"] != 1 or derived["phase_count_evidence"] != "source-backed-inference":
        raise ValueError("row-011 liquid-side interpretation changed")

    audit = contract["sampled_gibbs_audit"]
    if audit["case_ids"] != [
        "may2015-row-001-midpoint",
        "may2015-row-012-midpoint",
        "may2015-row-011-liquid-side",
    ]:
        raise ValueError("sampled-audit case selection changed")
    grid = audit["composition_grid"]
    spacing = (grid["maximum"] - grid["minimum"]) / (grid["points"] - 1)
    if grid["points"] != 1001 or not math.isclose(grid["spacing"], spacing, abs_tol=1e-15):
        raise ValueError("sampled-audit grid changed")
    if not math.isclose(grid["localization_allowance"], 2.0 * spacing, abs_tol=1e-15):
        raise ValueError("sampled-audit localization allowance changed")
    if audit["dimensionless_tangent_chord_allowance"] != 1e-6:
        raise ValueError("sampled-audit tangent/chord allowance changed")

    forbidden = {"model_output", "held_output", "phase_fractions", "model_phases"}
    if any(not forbidden.isdisjoint(case) for case in cases):
        raise ValueError("source-only case contract contains model output fields")
    return {
        "status": "cases_ready",
        "case_contract_sha256": sha256_file(cases_path),
        "source_sha256": SOURCE_SHA256,
        "case_count": len(cases),
        "midpoint_cases": len(midpoint_cases),
        "one_phase_cases": 1,
        "sampled_audit_cases": len(audit["case_ids"]),
        "globality_certificate": "not_guaranteed",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(check(args.cases, args.source), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

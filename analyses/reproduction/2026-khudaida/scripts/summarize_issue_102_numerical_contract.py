"""Summarize issue-102 before/after evidence without duplicating row payloads."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

CAMPAIGN_ROOT = Path(__file__).resolve().parents[1]
PRE_FIX = CAMPAIGN_ROOT / "results" / "data" / "issue-102-pre-fix-failures.json"
CURRENT = CAMPAIGN_ROOT / "results" / "data" / "khudaida-2026-figure-2-held2.csv"
SUMMARY = CAMPAIGN_ROOT / "results" / "khudaida-2026-figure-2-held2-summary.json"
OUTPUT = (
    CAMPAIGN_ROOT
    / "results"
    / "data"
    / "issue-102-numerical-contract-assessment.json"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    before = json.loads(PRE_FIX.read_text(encoding="utf-8"))
    with CURRENT.open(newline="", encoding="utf-8") as stream:
        current = list(csv.DictReader(stream))
    if [int(row["point_id"]) for row in current] != list(range(1, 9)):
        raise AssertionError("issue-102 assessment requires Figure 2 points 1-8")
    if not all(row["case_decision"] == "PASS" for row in current):
        raise AssertionError("issue-102 assessment requires eight accepted rows")

    diagnostics = [json.loads(row["diagnostics_json"]) for row in current]
    checks = [json.loads(row["independent_checks_json"]) for row in current]
    assessment = {
        "schema_version": 1,
        "inputs": {
            "pre_fix_failure_diagnostic": {
                "path": str(PRE_FIX.relative_to(CAMPAIGN_ROOT)),
                "sha256": sha256(PRE_FIX),
            },
            "current_result_rows": {
                "path": str(CURRENT.relative_to(CAMPAIGN_ROOT)),
                "sha256": sha256(CURRENT),
            },
            "current_result": {
                "path": str(SUMMARY.relative_to(CAMPAIGN_ROOT)),
                "sha256": sha256(SUMMARY),
            },
        },
        "paper_convergence_contract": {
            "perdomo_equation_69_ratio_before": 1.0e-3,
            "perdomo_equation_69_ratio_current": 1.0e-2,
            "role": "finite-search stopping gate only",
            "does_not_replace": [
                "material balance",
                "electroneutrality",
                "pressure stationarity",
                "KKT stationarity",
                "distinct two-phase topology",
            ],
        },
        "before": {
            "non_admission_count": len(before["rows"]),
            "non_admission_point_ids": [row["point_id"] for row in before["rows"]],
            "failure_reasons": dict(
                sorted(Counter(row["failure_reason"] for row in before["rows"]).items())
            ),
        },
        "after": {
            "accepted_count": len(current),
            "accepted_point_ids": [int(row["point_id"]) for row in current],
            "returned_phase_counts": sorted(
                {int(row["returned_phase_count"]) for row in current}
            ),
            "attempts_min": min(item["attempts"] for item in diagnostics),
            "attempts_max": max(item["attempts"] for item in diagnostics),
            "held_gap_max_abs": max(abs(item["held_gap"]) for item in diagnostics),
            "material_balance_max_abs": max(
                item["material_balance_max_abs"] for item in checks
            ),
            "phase_charge_max_abs": max(
                item["phase_charge_max_abs"] for item in checks
            ),
            "pressure_max_relative": max(
                item["pressure_max_relative"] for item in checks
            ),
            "kkt_stationarity_max_abs": max(
                item["kkt_stationarity_max_abs"] for item in diagnostics
            ),
            "acceptance_contracts": {
                "historical_output_abs_tolerance": 1.0e-7,
                "kkt_abs_tolerance": 1.0e-8,
                "physical_abs_tolerance": 1.0e-10,
                "pressure_relative_tolerance": 1.0e-8,
            },
        },
        "interpretation": (
            "The current exact installed artifacts accept all eight feeds while "
            "independent physical and KKT residuals remain far below their gates. "
            "This is bounded-search capability evidence, not a globality proof."
        ),
    }
    OUTPUT.write_text(
        json.dumps(assessment, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Keep the executable equilibrium claim ledger honest.

The design is deliberately the source of the claim names.  This check only
connects those names to the small native programs assigned by section 9; it
does not decide a thermodynamic result.
"""

from __future__ import annotations

import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / "docs/greenfield/EQUILIBRIUM.md"
CHECKS = {
    "specification_checks.cpp": set(range(1, 10)),
    "equilibrium_checks.cpp": {1, 2, 6, 7, 8, 10, 11, 13, 15, 16},
    "boundary_checks.cpp": {3, 4, 5, 9, 10, 14},
}


def design_claims() -> set[str]:
    text = SPEC.read_text(encoding="utf-8")
    start = text.index("### 4.1 Compile tier")
    end = text.index("### 4.3", start)
    section = text[start:end]
    claims = set(re.findall(r"\| (C\d+|S\d+) \|", section))
    return {claim for claim in claims if claim != "S12"}


def source_claims(path: Path) -> tuple[set[str], set[str]]:
    text = path.read_text(encoding="utf-8")
    functions = set(re.findall(r"\bvoid\s+(check_[cs]\d+)\s*\(", text))
    main = text[text.rfind("int main") :]
    calls = set(re.findall(r"\b(check_[cs]\d+)\s*\(\s*\)", main))
    return functions, calls


def main() -> int:
    expected = design_claims()
    found: dict[str, str] = {}
    failures: list[str] = []
    for filename, numbers in CHECKS.items():
        path = ROOT / "engine/tests/native" / filename
        functions, calls = source_claims(path)
        for number in numbers:
            prefix = "c" if filename == "specification_checks.cpp" else "s"
            claim = prefix + str(number)
            function = "check_" + claim
            found[claim.upper()] = filename
            if function not in functions:
                failures.append(f"{filename}: missing {function}")
            if function not in calls:
                failures.append(f"{filename}: {function} is not called from main")
    if expected != set(found):
        failures.append(
            "claim/file map differs: design="
            + ",".join(sorted(expected))
            + " files="
            + ",".join(sorted(found))
        )
    if failures:
        for failure in failures:
            print(failure)
        return 1
    print(f"ok {len(expected)} named equilibrium claims")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

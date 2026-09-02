#!/usr/bin/env python3
"""Structural rules for engine/native.

C++ is organised by dependency layer because include direction is the only
structural fact a checker can enforce. Five rules, all cheap, all stdlib:

  1. one static archive and one extension in meson.build
  2. no Python or pybind11 symbol inside the archive
  3. every "..." include resolves to the same or a lower layer
  4. regression/ may not reach phase stability
  5. no constexpr name defined twice, no category-1 constant outside physics/

Rule 2 needs a built archive and reports "skipped" without one.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NATIVE = ROOT / "engine" / "native"
MESON = ROOT / "engine" / "meson.build"

# Include direction: a layer may include itself and anything below it.
LAYERS = ["physics", "model", "equilibrium", "regression", "bindings"]
RANK = {name: index for index, name in enumerate(LAYERS)}

# Names that may only be defined in physics/physical_constants.hpp. The tree
# this replaces had fifteen definitions of the gas constant.
CATEGORY_ONE = {
    "GAS_CONSTANT_J_PER_MOL_K",
    "BOLTZMANN_J_PER_K",
    "AVOGADRO_PER_MOL",
    "ELEMENTARY_CHARGE_C",
    "VACUUM_PERMITTIVITY_F_PER_M",
    "SPEED_OF_LIGHT_M_PER_S",
    "PI",
    "ANGSTROM_M",
    "DEBYE_C_M",
    "NUMBER_DENSITY_PER_ANGSTROM3",
}
CONSTANTS_HOME = "physics/physical_constants.hpp"

# Phase stability is reachable from equilibrium, never from a parameter fit: a
# fit runs against declared phases so its residuals stay differentiable.
STABILITY_SYMBOLS = ("phase_stability.hpp", "find_phases")

INCLUDE = re.compile(r'^\s*#\s*include\s+"([^"]+)"', re.MULTILINE)
# Unanchored on purpose: a constexpr sharing a name across two translation
# units is drift whether it sits at namespace scope or inside a function.
CONSTEXPR = re.compile(r"\bconstexpr\s+[\w:<>,\s*&]+?\b(\w+)\s*(?:=|\{)")


def sources() -> list[Path]:
    if not NATIVE.is_dir():
        return []
    return sorted(p for p in NATIVE.rglob("*") if p.suffix in {".cpp", ".hpp"})


def layer_of(path: Path) -> str:
    return path.relative_to(NATIVE).parts[0]


def rule_one_archive() -> list[str]:
    if not MESON.is_file():
        return [f"{MESON} is missing"]
    text = MESON.read_text()
    archives = len(re.findall(r"\bstatic_library\s*\(", text))
    extensions = len(re.findall(r"\bextension_module\s*\(", text))
    bad = []
    if archives != 1:
        bad.append(f"meson.build declares {archives} static libraries, expected 1")
    if extensions != 1:
        bad.append(f"meson.build declares {extensions} extensions, expected 1")
    return bad


def rule_archive_is_python_free() -> tuple[list[str], bool]:
    found = list((ROOT / "engine").glob("build/**/libepcsaft_native.a"))
    if not found:
        return [], False
    out = subprocess.run(
        ["nm", "--defined-only", str(found[0])],
        capture_output=True,
        text=True,
        check=False,
    ).stdout
    leaked = sorted(
        {
            line.split()[-1]
            for line in out.splitlines()
            if "PyInit" in line or "pybind11" in line
        }
    )
    if leaked:
        return [f"{found[0].name} exports Python symbols: {leaked[:5]}"], True
    return [], True


def rule_include_direction() -> list[str]:
    bad = []
    for path in sources():
        here = layer_of(path)
        for target in INCLUDE.findall(path.read_text()):
            head = target.split("/")[0]
            if head not in RANK:
                continue  # a sibling header in the same directory
            if RANK[head] > RANK[here]:
                rel = path.relative_to(ROOT)
                bad.append(f"{rel} ({here}) includes {target} ({head}), upward")
    return bad


def rule_regression_stays_below_stability() -> list[str]:
    bad = []
    for path in sources():
        if layer_of(path) != "regression":
            continue
        text = path.read_text()
        for symbol in STABILITY_SYMBOLS:
            if symbol in text:
                bad.append(f"{path.relative_to(ROOT)} reaches {symbol}")
    return bad


def rule_constants_defined_once() -> list[str]:
    homes: dict[str, list[str]] = {}
    for path in sources():
        rel = str(path.relative_to(NATIVE))
        for name in CONSTEXPR.findall(path.read_text()):
            homes.setdefault(name, []).append(rel)
    bad = []
    for name, files in sorted(homes.items()):
        if len(files) > 1:
            bad.append(f"constexpr {name} defined in {len(files)}: {sorted(files)}")
        if name in CATEGORY_ONE and files != [CONSTANTS_HOME]:
            bad.append(f"category-1 {name} defined outside {CONSTANTS_HOME}: {files}")
    return bad


# A threshold is a decision and belongs in numerical_policy.hpp. Only literals with a
# decimal point or exponent are considered: a bare integer after "<" is a template
# argument (std::array<double, 4>), not a comparison.
THRESHOLD = re.compile(
    r"(?:<=|>=|==|!=|<|>)\s*(\d+\.\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?|\d+[eE][-+]?\d+)"
    r"|(\d+\.\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?|\d+[eE][-+]?\d+)\s*(?:<=|>=|==|!=|<|>)"
)
THRESHOLD_EXEMPT = {"physics/numerical_policy.hpp", "physics/physical_constants.hpp"}
ALLOWED_THRESHOLDS = {"0.0", "1.0"}  # sign tests and normalisation, nothing else

# A solver receives rows. It may not know which kind of problem produced them: that is
# what makes the route general rather than a union of special cases. request_compile.cpp
# is the one file allowed to branch, and only on which variable is specified.
# Underscore counts as a boundary: snake_case is this tree's convention, so \b would let
# solve_reactive_case through and leave the rule vacuous.
CASE_TOKENS = re.compile(
    r"(?<![A-Za-z0-9])(reactive|nonreactive|incipient|finite|electrolyte|charged|ionic|"
    r"neutral|binary|ternary|saturation|bubble|dew|vle|lle|amine|aqueous|salt)"
    r"(?![A-Za-z0-9])", re.IGNORECASE)
SOLVER_FILES = {
    "equilibrium/gibbs_minimizer.hpp", "equilibrium/gibbs_minimizer.cpp",
    "equilibrium/phase_boundary.cpp", "equilibrium/implicit_sensitivity.cpp",
    "equilibrium/phase_stability.cpp",
}

COMMENT = re.compile(r"//[^\n]*|/\*.*?\*/", re.S)


def without_comments(text: str) -> str:
    return COMMENT.sub("", text)


def rule_no_threshold_literals() -> list[str]:
    bad = []
    for path in sources():
        rel = str(path.relative_to(NATIVE))
        if rel in THRESHOLD_EXEMPT:
            continue
        for number, line in enumerate(without_comments(path.read_text()).splitlines(), 1):
            for match in THRESHOLD.finditer(line):
                literal = match.group(1) or match.group(2)
                if literal not in ALLOWED_THRESHOLDS:
                    bad.append(f"{rel}:{number}: threshold literal {literal}")
    return bad


def rule_solvers_are_case_blind() -> list[str]:
    bad = []
    for path in sources():
        rel = str(path.relative_to(NATIVE))
        if rel not in SOLVER_FILES:
            continue
        text = without_comments(path.read_text())
        for match in CASE_TOKENS.finditer(text):
            bad.append(f"{rel}: names a problem kind: {match.group()}")
    return bad


def main() -> int:
    archive_findings, archive_built = rule_archive_is_python_free()
    checks = [
        ("one archive, one extension", rule_one_archive()),
        (
            "archive is Python-free" if archive_built else "archive is Python-free (skipped, not built)",
            archive_findings,
        ),
        ("include direction is downward", rule_include_direction()),
        ("regression stays below stability", rule_regression_stays_below_stability()),
        ("constants defined once", rule_constants_defined_once()),
        ("no threshold literals", rule_no_threshold_literals()),
        ("solvers are case-blind", rule_solvers_are_case_blind()),
    ]
    failed = 0
    for name, findings in checks:
        print(f"{'FAIL' if findings else 'ok  '}  {name}")
        for finding in findings:
            print(f"        {finding}")
        failed += bool(findings)
    print(f"\n{len(sources())} sources under engine/native, {failed} rule(s) failing")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

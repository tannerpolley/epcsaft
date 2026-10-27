#!/usr/bin/env python3
"""Structural rules for engine/native and the public Python surface.

C++ is organised by dependency layer because include direction is the only
structural fact a checker can enforce. Each rule is a cheap stdlib reading named
in main()'s check table. The Python-free archive rule needs a built archive and
reports "skipped" without one; the stub rule, which imports the extension, runs only with --build-dir.
"""

from __future__ import annotations

import argparse
import ast
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NATIVE = ROOT / "engine" / "native"
MESON = ROOT / "engine" / "meson.build"

# Include direction: a layer may include itself and anything below it.
LAYERS = ["core", "contributions", "properties", "equilibrium", "regression", "surface", "bindings"]
RANK = {name: index for index, name in enumerate(LAYERS)}

# Names that may only be defined in core/physical_constants.hpp. The tree
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
CONSTANTS_HOME = "core/physical_constants.hpp"

# Phase stability is reachable from equilibrium, never from a parameter fit: a
# fit runs against declared phases so its residuals stay differentiable.
STABILITY_SYMBOLS = ("stability_test", "discover_phases", "StabilityOptions")

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
    archives = re.findall(r"(\w+)\s*=\s*static_library\s*\(\s*'([^']+)'", text)
    extensions = len(re.findall(r"\bextension_module\s*\(", text))
    bad = []
    if archives != [("nanobind_support", "nanobind_support"), ("epcsaft_private", "epcsaft_native")]:
        bad.append(f"expected one private native archive and nanobind support, found {archives}")
    if extensions != 1:
        bad.append(f"meson.build declares {extensions} extensions, expected 1")
    return bad


def rule_archive_is_python_free(build_dir: Path | None = None) -> tuple[list[str], bool]:
    if build_dir is None:
        found = sorted((ROOT / "engine").glob("build/**/libepcsaft_native.a"))
    else:
        build_dir = build_dir.resolve()
        if not build_dir.is_dir():
            return [f"configured build directory does not exist: {build_dir}"], True
        if not (build_dir / "meson-info" / "meson-info.json").is_file():
            return [f"not a configured Meson build: {build_dir}"], True
        found = sorted(build_dir.glob("libepcsaft_native.a"))
        if len(found) != 1:
            return [f"expected one native archive in {build_dir}, found {found}"], True
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
            if line.split() and not line.endswith(":") and re.search(r"^_?Py\w*|pybind11|nanobind", line.split()[-1])
        }
    )
    if leaked:
        return [f"{found[0].name} exports Python symbols: {leaked[:5]}"], True
    if build_dir is not None:
        ninja = (build_dir / "build.ninja").read_text().splitlines()
        def link(name):
            return next((line.split(" |", 1)[0].split()[3:] for line in ninja
                         if line.startswith(f"build {name}: ")), [])
        archive_objects = link("libepcsaft_native.a")
        shared_objects = [word for word in link("libepcsaft_native.so") if word.endswith(".o")]
        if not archive_objects or archive_objects != shared_objects:
            return ["native shared library and private archive do not reuse one PIC object set"], True
        support_consumers = []
        for line in ninja:
            if line.startswith("build ") and ": cpp_LINKER " in line:
                output = line.split()[1].removesuffix(":")
                inputs = line.split()[3:]
                if "libnanobind_support.a" in inputs:
                    support_consumers.append(output)
                if "libepcsaft_native.a" in inputs and output != "libepcsaft_asl.so":
                    return [f"private archive used outside ASL: {output}"], True
        if len(support_consumers) != 1 or not support_consumers[0].startswith("_core."):
            return [f"nanobind support must link only into _core: {support_consumers}"], True
        dynamic = subprocess.check_output(["readelf", "-d", str(build_dir / "libepcsaft_asl.so")], text=True)
        exports = subprocess.check_output(["nm", "-D", "--defined-only", str(build_dir / "libepcsaft_asl.so")], text=True)
        if "libepcsaft_native.so" in dynamic or [line.split()[-1] for line in exports.splitlines()] != ["funcadd_ASL"]:
            return ["ASL private native copy is not localized"], True
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
THRESHOLD_EXEMPT = {"core/numerical_policy.hpp", "core/physical_constants.hpp"}
ALLOWED_THRESHOLDS = {"0.0", "1.0"}  # sign tests and normalisation, nothing else

# A solver receives rows. It may not know which kind of problem produced them: that is
# what makes the route general rather than a union of special cases. request_compile.cpp
# is the one file allowed to branch, and only on which variable is specified.
# Underscore counts as a boundary: snake_case is this tree's convention, so \b would let
# solve_reactive_case through and leave the rule vacuous.
CASE_TOKENS = re.compile(
    r"(?<![A-Za-z0-9])(reactive|nonreactive|incipient|electrolyte|charged|ionic|"
    r"neutral|binary|ternary|saturation|bubble|dew|vle|lle|amine|aqueous|salt)"
    r"(?![A-Za-z0-9])", re.IGNORECASE)
SOLVER_FILES = {
    "equilibrium/gibbs_minimizer.hpp", "equilibrium/gibbs_minimizer.cpp",
    "equilibrium/phase_boundary.cpp", "equilibrium/implicit_sensitivity.cpp",
    "equilibrium/request_compile.cpp",
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
        # Diagnostic strings describe physics; they do not select a solver.
        # 'finite' likewise describes a quantity or split, not a problem kind.
        text = re.sub(r'"(?:\\.|[^"\\])*"', '""', without_comments(path.read_text()))
        for match in CASE_TOKENS.finditer(text):
            bad.append(f"{rel}: names a problem kind: {match.group()}")
    return bad


def rule_policy_entries_name_themselves() -> list[str]:
    path = NATIVE / "core" / "numerical_policy.hpp"
    if not path.is_file():
        return [f"{path.relative_to(ROOT)} is missing"]
    text = path.read_text()
    bad = []
    for identifier, name in re.findall(
        r"inline\s+constexpr\s+Tolerance\s+(\w+)\s*\{\s*\"([^\"]+)\"",
        text,
    ):
        if identifier != name:
            bad.append(f"{path.relative_to(ROOT)}: {identifier} names itself {name}")
    for identifier, name in re.findall(
        r"inline\s+constexpr\s+Tolerance\s+(\w+)\s*=\s*derived\([^;]*?\"([^\"]+)\"",
        text,
        re.S,
    ):
        if identifier != name:
            bad.append(f"{path.relative_to(ROOT)}: {identifier} names itself {name}")
    listed = set(re.findall(r"&(\w+)", text.split("ENTRIES[]", 1)[1])) if "ENTRIES[]" in text else set()
    consumers = "\n".join(without_comments(source.read_text()) for source in sources() if source != path)
    # The six-option owner consumes policies through its indexed default table.
    options = re.search(r"SOLVER_POLICIES = \{([^}]+)\}", text)
    if options and "value_or(SOLVER_POLICIES.at(i)->value)" in text:
        consumers += options.group(1)
    for identifier in re.findall(r"inline\s+constexpr\s+Tolerance\s+(\w+)\s*[{=]", text):
        if identifier not in listed:
            bad.append(f"{path.relative_to(ROOT)}: {identifier} is missing from ENTRIES")
        if not re.search(rf"\b{identifier}\b", consumers) and not re.search(rf"derived\({identifier}\b", text):
            bad.append(f"{path.relative_to(ROOT)}: {identifier} has no consumer in engine/native")
    if "inline const Tolerance" in text:
        bad.append(f"{path.relative_to(ROOT)} still has inline const Tolerance")
    return bad


def rule_budget_readers_carry_records() -> list[str]:
    budget = re.compile(r"policy::\w+_(?:MAX_ITERATIONS|BUDGET)\b")
    bad = []
    for path in sources():
        if layer_of(path) in {"core", "bindings"}:
            continue
        text = without_comments(path.read_text())
        if budget.search(text) and ("core::Scope" not in text or "Convergence" not in text):
            bad.append(f"{path.relative_to(ROOT)} reads a budget without a convergence record")
    return bad


def rule_core_and_trace_boundaries() -> list[str]:
    bad = []
    forbidden = re.compile(
        r"#\s*include\s+<(cppad|Eigen|Ipopt|coin|ceres|Highs|pybind11|nanobind|Python\.h)",
        re.IGNORECASE,
    )
    # An explicit operation may own a Capture for observed work counts. It
    # must not consult the ambient trace to choose or alter a calculation.
    read_back = re.compile(r"(?<![A-Za-z0-9_])(core::)?active\s*\(")
    for path in sources():
        rel = path.relative_to(NATIVE)
        parts = rel.parts
        text = without_comments(path.read_text())
        if parts[0] != "bindings" and re.search(r'#\s*include\s*[<"](?:nanobind/|pybind11/|Python\.h)', text):
            bad.append(f"{path.relative_to(ROOT)} includes Python binding headers")
        if parts[0] == "core":
            if path.name != "scalar.hpp" and forbidden.search(text):
                bad.append(f"{path.relative_to(ROOT)} includes a numerical library")
            for target in INCLUDE.findall(text):
                if target.split("/")[0] != "core":
                    bad.append(f"{path.relative_to(ROOT)} includes {target} outside core")
        elif parts[0] != "bindings" and read_back.search(text):
            bad.append(f"{path.relative_to(ROOT)} reads the active trace")
    return bad


# Preserve offsets so findings refer to the original source, including comments.
CPP_NON_CODE = re.compile(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'', re.S)
TRACE_CALL = re.compile(r'(?:\.|->)\s*(?:probe|emit)\s*\(')
PAYLOAD_TYPE = r'std::(?:vector|array|deque|list|map|unordered_map|set|unordered_set|string|basic_string)(?:\s*<[^;{}()]*>)?'
PAYLOAD_TEMPORARY = re.compile(r"\b" + PAYLOAD_TYPE + r"\s*[({]|\bstd::to_string\s*\(")
PAYLOAD_DECLARATION = re.compile(
    rf'\b(?:{PAYLOAD_TYPE}|auto(?=\s+\w+\s*=\s*(?:{PAYLOAD_TYPE}|std::to_string)))\s+([^;]+);'
)


def rule_known_eager_trace_payload_forms() -> list[str]:
    # Regression guard for observed eager-temporary and named-local patterns, not a complete allocation-safety proof.
    bad = []
    returned_payloads = set()
    for path in sources():
        for match in re.finditer(rf'\b{PAYLOAD_TYPE}\s+(\w+(?:::\w+)*)\s*\(', without_comments(path.read_text())):
            returned_payloads.add(match.group(1).split('::')[-1])
    returned_names = '|'.join(sorted(returned_payloads)) or r'(?!)'
    returned_call = rf'\b(?:{returned_names})\s*\('
    temporary = re.compile(PAYLOAD_TEMPORARY.pattern + '|' + returned_call)
    declaration_pattern = re.compile(PAYLOAD_DECLARATION.pattern.replace('|std::to_string)', '|std::to_string|' + returned_names + ')'))
    for path in sources():
        text = CPP_NON_CODE.sub(lambda m: re.sub(r'[^\n]', ' ', m.group()), path.read_text())
        # Pair statement blocks and calls; initializer braces are harmless here.
        pairs, stack = {}, []
        for i, char in enumerate(text):
            if char in '({':
                stack.append(i)
            elif char in ')}' and stack and text[stack[-1]] == ('(' if char == ')' else '{'):
                pairs[stack.pop()] = i
        observations = []
        for match in re.finditer(r'\bobserve\s*\([^;{}]*\[[^\]]*\][^{]*\{', text):
            start = match.end() - 1
            if start in pairs:
                observations.append((start, pairs[start]))
        calls = [(m.start(), pairs[m.end() - 1]) for m in TRACE_CALL.finditer(text)
                 if m.end() - 1 in pairs]
        for start, end in calls:
            payloads = []
            observed = any(left < start < right for left, right in observations)
            if not observed and temporary.search(text[start:end]):
                payloads.append('container/string temporary')
            blocks = [(left, right) for left, right in pairs.items()
                      if text[left] == '{' and left < start < right]
            for left, right in blocks:
                block = text[left + 1:right]
                # Remove emissions before checking whether a local is only a payload.
                remaining = block
                for call_start, call_end in reversed(calls):
                    if left < call_start < call_end < right:
                        a, b = call_start - left - 1, call_end - left
                        remaining = remaining[:a] + ' ' * (b - a) + remaining[b:]
                for declaration in declaration_pattern.finditer(block[:start - left - 1]):
                    position = left + 1 + declaration.start()
                    if any(a < position < b for a, b in observations):
                        continue
                    if any(a < position < b and a > left for a, b in blocks):
                        continue
                    names = re.findall(r'(?:^|,)\s*(\w+)\s*(?=[({=]|$)', declaration.group(1))
                    for name in names:
                        if not re.search(rf'\b{name}\b', text[start:end]):
                            continue
                        uses = remaining[:declaration.start()] + ' ' * len(declaration.group()) + remaining[declaration.end():]
                        # Filling the payload is allowed; any other read earns numerical ownership.
                        uses = re.sub(rf'\b{name}\s*(?:\[[^\]]*\]\s*=|\.(?:push_back|emplace_back|reserve|resize)\s*\()', '', uses)
                        if not re.search(rf'\b{name}\b', uses):
                            payloads.append(name)
            if payloads:
                line = text.count('\n', 0, start) + 1
                bad.append(f'{path.relative_to(ROOT)}:{line}: eager trace payload: {", ".join(sorted(set(payloads)))}')
    return bad


def rule_declared_selectors() -> list[str]:
    # Every declared document: Parameters today, Thermochemistry once it exists.
    documents = (("parameters.hpp", "Parameters", True), ("thermochemistry.hpp", "Thermochemistry", False))
    bad = []
    for name, struct, required in documents:
        path = NATIVE / "core" / name
        if not path.is_file():
            if required:
                bad.append(f"{path.relative_to(ROOT)} is missing")
            continue
        text = path.read_text()
        if "Unspecified" in text:
            bad.append(f"{path.relative_to(ROOT)} contains Unspecified")
        enum_names = re.findall(r"enum\s+class\s+(\w+)", text)
        match = re.search(rf"struct\s+{struct}\b.*?\n}};", text, re.S)
        if match and enum_names:
            enum_pattern = "|".join(map(re.escape, enum_names))
            bare = re.compile(rf"^\s*(?:{enum_pattern})\s+\w+\s*[=;{{]")
            for number, line in enumerate(match.group().splitlines(), 1):
                if bare.search(line):
                    bad.append(f"{path.relative_to(ROOT)}:{number}: bare enum selector")
    return bad


def rule_surface_defaults() -> list[str]:
    bad = []
    surface = ROOT / "engine" / "src" / "epcsaft"
    for path in sorted(surface.rglob("*.py")):
        try:
            tree = ast.parse(path.read_text())
        except SyntaxError as error:
            bad.append(f"{path.relative_to(ROOT)}:{error.lineno}: syntax error")
            continue
        for node in __import__("ast").walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            defaults = [*node.args.defaults, *[item for item in node.args.kw_defaults if item is not None]]
            for default in defaults:
                if not isinstance(default, ast.Constant) or default.value is not None:
                    bad.append(f"{path.relative_to(ROOT)}:{node.lineno}: default is not None")
    return bad


BOUND_FUNCTION = re.compile(
    r'\.def\(\s*"[^"]+"\s*,\s*&(?:epcsaft::)?'
    r'(?:(?P<class>\w+)::)?(?P<function>\w+)'
)
NON_SCOPED_BOUND_FUNCTIONS: set[str] = set()  # rule 13 has no allowlist (owner ruling)
HAND_NAMED_RESIDUAL_STOP = re.compile(
    r'(?:\.stop|->stop)\s*=\s*core::Stop::'
    r'(?:converged|converged_to_floor|stalled|budget_exhausted)\b'
)


def rule_bound_functions_open_a_scope() -> list[str]:
    bound = []
    for path in sorted((NATIVE / "bindings").glob("*.cpp")):
        for match in BOUND_FUNCTION.finditer(without_comments(path.read_text())):
            function = match.group("function")
            if function in NON_SCOPED_BOUND_FUNCTIONS:
                continue
            bound.append((match.group("class"), function))
    bad = []
    for owner, function in bound:
        qualified = f"{owner}::{function}" if owner else function
        pattern = re.compile(
            rf"\b{re.escape(qualified)}\s*\([^;{{]*\)\s*(?:const\s*)?\{{"
            r"\s*core::Scope\s+scope\s*\(",
            re.S,
        )
        if not any(pattern.search(without_comments(path.read_text())) for path in sources()):
            bad.append(f"bound function {qualified} must open core::Scope at entry")
    return bad


def rule_no_hand_named_residual_stop() -> list[str]:
    bad = []
    for path in sources():
        text = without_comments(path.read_text())
        for number, line in enumerate(text.splitlines(), 1):
            if HAND_NAMED_RESIDUAL_STOP.search(line):
                bad.append(f"{path.relative_to(ROOT)}:{number}: residual stop must use classify")
    return bad


def rule_exact_hessian_callback() -> list[str]:
    """The flash driver must expose the exact Lagrangian Hessian to Ipopt."""
    path = NATIVE / "equilibrium" / "gibbs_minimizer.cpp"
    if not path.is_file():
        return [f"{path.relative_to(ROOT)} is missing"]
    text = without_comments(path.read_text())
    bad = []
    if "GibbsMinimizer::eval_h" not in text:
        bad.append("equilibrium/gibbs_minimizer.cpp has no eval_h callback")
    options = without_comments((NATIVE / "equilibrium" / "ipopt.hpp").read_text())
    call = re.search(r"run_ipopt\((.*?)\);", text, re.S)
    if call is None or call.group(1).rsplit(",", 6)[-6].strip() != "true" or not re.search(
        r'"hessian_approximation",\s*exact\s*\?\s*"exact"\s*:\s*"limited-memory"', options
    ):
        bad.append("solve_equilibrium must request exact curvature from the common Ipopt owner")
    if "limited-memory" in text:
        bad.append("limited-memory Hessian is forbidden")
    return bad


def rule_stub_matches_extension() -> list[str]:
    """_core.pyi declares exactly the installed extension's public names and class members."""
    from epcsaft import _core

    def declared(body: list[ast.stmt]) -> set[str]:
        names = {node.name for node in body if isinstance(node, (ast.ClassDef, ast.FunctionDef))}
        names |= {node.target.id for node in body if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name)}
        return {name for name in names if not name.startswith("_")}

    def mismatch(where: str, stub: set[str], native: set[str]) -> list[str]:
        return [f"{where}: stub only {sorted(stub - native)}, extension only {sorted(native - stub)}"] if stub != native else []

    tree = ast.parse((ROOT / "engine" / "src" / "epcsaft" / "_core.pyi").read_text())
    findings = mismatch("_core", declared(tree.body), {name for name in dir(_core) if not name.startswith("_")})
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and hasattr(_core, node.name):
            native = getattr(_core, node.name)
            # Enum implementations supply name/value accessors.
            implicit = {"name", "value"} if hasattr(native, "__members__") else set()
            members = {name for name in vars(native) if not name.startswith("_")} - implicit
            instance_fields = {n.target.id for n in node.body if isinstance(n, ast.AnnAssign) and isinstance(n.target, ast.Name)} if issubclass(native, BaseException) else set()
            # Native translators attach exception fields to instances, not their classes.
            findings += mismatch(node.name, declared(node.body) - instance_fields, members)
    return findings


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path)
    arguments = parser.parse_args()

    archive_findings, archive_built = rule_archive_is_python_free(arguments.build_dir)
    checks = [
        ("one PIC object set, shared/private libraries, one extension", rule_one_archive()),
        (
            "archive is Python-free" if archive_built else "archive is Python-free (skipped, not built)",
            archive_findings,
        ),
        ("include direction is downward", rule_include_direction()),
        ("regression stays below stability", rule_regression_stays_below_stability()),
        ("constants defined once", rule_constants_defined_once()),
        ("no threshold literals", rule_no_threshold_literals()),
        ("solvers are case-blind", rule_solvers_are_case_blind()),
        ("policy entries are named, listed and consumed", rule_policy_entries_name_themselves()),
        ("budget readers carry records", rule_budget_readers_carry_records()),
        ("core and trace boundaries", rule_core_and_trace_boundaries()),
        ("known eager trace-payload forms are absent", rule_known_eager_trace_payload_forms()),
        ("selectors are declared", rule_declared_selectors()),
        ("surface has no physical defaults", rule_surface_defaults()),
        ("bound functions open a Scope", rule_bound_functions_open_a_scope()),
        ("no hand-named residual stops", rule_no_hand_named_residual_stop()),
        ("exact Hessian callback", rule_exact_hessian_callback()),
    ]
    if arguments.build_dir is not None:
        checks.append(("type stub matches the extension", rule_stub_matches_extension()))
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

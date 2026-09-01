"""Verify the configured native Engine graph and identified installed wheel."""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import json
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Callable, Iterable
from pathlib import Path
from zipfile import ZipFile

STATIC_TARGETS = {
    "epcsaft_engine_eos_core": "libepcsaft_engine_eos_core.a",
    "epcsaft_engine_equilibrium_core": "libepcsaft_engine_equilibrium_core.a",
    "epcsaft_engine_regression_core": "libepcsaft_engine_regression_core.a",
}
EXPECTED_ARCHIVE_ORDER = tuple(STATIC_TARGETS.values())
EXPECTED_SIBLING_EDGES = {
    "epcsaft_engine_eos_core": set(),
    "epcsaft_engine_equilibrium_core": {"epcsaft_engine_eos_core"},
    "epcsaft_engine_regression_core": {
        "epcsaft_engine_eos_core",
        "epcsaft_engine_equilibrium_core",
    },
}
NATIVE_TARGET_TYPES = {"static library", "shared library", "shared module"}
RETIRED_SOURCE_IDENTIFIERS = {
    "_equilibrium": re.compile(r"(?<![A-Za-z0-9_])_equilibrium(?![A-Za-z0-9_])"),
    "epcsaft_native_transport": re.compile(r"\bepcsaft_native_transport\b"),
    "NativeTransportEndpoint": re.compile(r"\bNativeTransportEndpoint\b"),
    "native_transport.cpp": re.compile(r"\bnative_transport\.cpp\b"),
    "native_transport.hpp": re.compile(r"\bnative_transport\.hpp\b"),
    "native_transport_pybind.cpp": re.compile(r"\bnative_transport_pybind\.cpp\b"),
    "native_transport_pybind.hpp": re.compile(r"\bnative_transport_pybind\.hpp\b"),
}
FORBIDDEN_SYMBOLS = {
    "Python C-API": re.compile(
        r"(?:^|_)Py(?:Object|Type|Module|Unicode|Long|Float|Tuple|List|Dict|Err|GIL|Eval|Import|Capsule)[A-Za-z0-9_]*"
    ),
    "pybind11": re.compile(r"pybind11"),
    "GIL": re.compile(r"(?:PyGILState|PyEval_(?:Save|Restore)Thread|gil_scoped)"),
    **RETIRED_SOURCE_IDENTIFIERS,
    "retired _equilibrium extension": re.compile(r"PyInit__equilibrium$"),
}


class ArchitectureError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ArchitectureError(message)


def run(command: list[str], *, cwd: Path | None = None) -> str:
    completed = subprocess.run(
        command,
        cwd=cwd,
        check=False,
        text=True,
        capture_output=True,
    )
    if completed.returncode:
        raise ArchitectureError(
            f"command failed ({completed.returncode}): {shlex.join(command)}\n"
            f"{completed.stderr.strip()}"
        )
    return completed.stdout


def introspect(meson: str, build_dir: Path, option: str) -> object:
    return json.loads(run([meson, "introspect", str(build_dir), option]))


def validate_products(
    targets: list[dict[str, object]], install_plan: dict[str, object]
) -> dict[str, dict[str, object]]:
    roots = [
        target
        for target in targets
        if target.get("subproject") is None
        and target.get("type") in NATIVE_TARGET_TYPES
    ]
    by_name = {str(target["name"]): target for target in roots}
    for name in STATIC_TARGETS:
        require(name in by_name, f"missing root native Meson target: {name}")
        require(
            by_name[name].get("type") == "static library",
            f"{name} is not a static library",
        )

    installed_targets = install_plan.get("targets", {})
    require(isinstance(installed_targets, dict), "Meson install plan has no target map")
    core_outputs = [
        Path(str(output))
        for output, record in installed_targets.items()
        if isinstance(record, dict)
        and record.get("subproject") is None
        and re.search(r"/epcsaft/_core\.[^/]+$", str(record.get("destination", "")))
    ]
    require(
        len(core_outputs) == 1,
        f"expected one installed epcsaft._core target, found {core_outputs}",
    )
    core_output = core_outputs[0]
    core_targets = [
        target
        for target in roots
        if any(
            Path(str(filename)).name == core_output.name
            for filename in target.get("filename", [])
        )
    ]
    require(
        len(core_targets) == 1,
        f"installed _core output does not map to one root target: {core_output}",
    )
    core = core_targets[0]
    require(
        core.get("type") == "shared module",
        "epcsaft._core is not a Meson shared module",
    )

    expected_ids = {id(by_name[name]) for name in STATIC_TARGETS} | {id(core)}
    extras = [
        str(target.get("name")) for target in roots if id(target) not in expected_ids
    ]
    require(not extras, f"unexpected root native Meson products: {extras}")
    by_name["_core"] = core
    return by_name


def command_for_output(ninja: str, build_dir: Path, output: Path) -> str:
    relative = output.relative_to(build_dir)
    commands = run([ninja, "-C", str(build_dir), "-t", "commands", str(relative)])
    matches = []
    for line in commands.splitlines():
        tokens = shlex.split(line)
        if any(Path(token).name == output.name for token in tokens):
            matches.append(line)
            continue
        for index, token in enumerate(tokens[:-1]):
            if token == "-o" and Path(tokens[index + 1]).name == output.name:
                matches.append(line)
                break
    require(
        len(matches) == 1,
        f"expected one backend command producing {relative}, found {len(matches)}",
    )
    return matches[0]


def validate_link_command(command: str) -> None:
    archives = [
        Path(token).name
        for token in shlex.split(command)
        if Path(token).name.startswith("libepcsaft_engine_") and token.endswith(".a")
    ]
    require(
        tuple(archives) == EXPECTED_ARCHIVE_ORDER,
        "forbidden reverse, sibling, missing, or extra internal link; expected backend archive "
        f"order {EXPECTED_ARCHIVE_ORDER}, found {tuple(archives)}",
    )


def symbol_fields(nm: str, archive: Path) -> list[tuple[str, str]]:
    fields: list[tuple[str, str]] = []
    for line in run([nm, "-P", "-g", str(archive)]).splitlines():
        stripped = line.strip()
        if not stripped or stripped.endswith(":"):
            continue
        row = stripped.split()
        require(len(row) >= 2, f"unrecognized nm -P row for {archive}: {line}")
        fields.append((row[0], row[1]))
    return fields


def validate_symbols(symbols: Iterable[tuple[str, str]], owner: str) -> None:
    for symbol, _kind in symbols:
        for label, pattern in FORBIDDEN_SYMBOLS.items():
            if pattern.search(symbol):
                raise ArchitectureError(
                    f"{owner} symbol field contains forbidden {label} pattern {pattern.pattern!r}: {symbol}"
                )


def validate_sibling_edges(
    archives: dict[str, list[tuple[str, str]]],
) -> dict[tuple[str, str], set[str]]:
    defined = {
        owner: {symbol for symbol, kind in fields if kind in "TDBR"}
        for owner, fields in archives.items()
    }
    undefined = {
        owner: {symbol for symbol, kind in fields if kind == "U"}
        for owner, fields in archives.items()
    }
    intersections = {
        (consumer, eos): undefined[consumer] & defined[eos]
        for consumer in archives
        for eos in archives
        if consumer != eos
    }
    for (consumer, eos), symbols in intersections.items():
        expected = eos in EXPECTED_SIBLING_EDGES[consumer]
        require(
            bool(symbols) == expected,
            f"{consumer} -> {eos} symbol edge is "
            f"{'missing' if expected else 'forbidden'}",
        )
    return intersections


def production_sources(
    targets: Iterable[dict[str, object]], source_root: Path
) -> list[Path]:
    sources: set[Path] = set()
    for target in targets:
        for group in target.get("target_sources", []):
            if not isinstance(group, dict):
                continue
            for raw in group.get("sources", []):
                path = Path(str(raw)).resolve()
                try:
                    relative = path.relative_to(source_root)
                except ValueError:
                    continue
                if (
                    relative.parts
                    and relative.parts[0] in {"cpp", "src"}
                    and path.is_file()
                ):
                    sources.add(path)
    for production_root in (source_root / "cpp", source_root / "src"):
        sources.update(
            path
            for path in production_root.rglob("*")
            if path.is_file()
            and path.suffix
            in {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx", ".py", ".pyi"}
        )
    return sorted(sources)


def validate_sources(sources: Iterable[Path]) -> None:
    for source in sources:
        for label, pattern in RETIRED_SOURCE_IDENTIFIERS.items():
            if source.name == label:
                raise ArchitectureError(
                    f"production source retains forbidden filename {label}: {source}"
                )
        text = source.read_text(encoding="utf-8", errors="replace")
        for label, pattern in RETIRED_SOURCE_IDENTIFIERS.items():
            if pattern.search(text):
                raise ArchitectureError(
                    f"production source contains forbidden {label} pattern {pattern.pattern!r}: {source}"
                )


def _regression_discovery_family(target: str) -> str | None:
    leaf = target.rsplit(".", 1)[-1]
    if leaf == "PhaseEquilibriumProblem":
        return "PhaseEquilibriumProblem"
    if target in {"epcsaft.equilibrium.solve", "epcsaft.equilibrium._api.solve"}:
        return "public equilibrium.solve"
    if leaf == "phase_equilibrium":
        return "phase_equilibrium"
    if target.startswith(("epcsaft.equilibrium", "epcsaft._core")):
        lowered = leaf.lower()
        if "held" in lowered and any(
            token in lowered for token in ("controller", "discover", "result", "solve")
        ):
            return "HELD/HELD2 controller"
        if "discover" in lowered or "discovery" in lowered or leaf == "*":
            return "discovery owner"
    return None


def _import_module(node: ast.ImportFrom) -> str:
    if not node.level:
        return node.module or ""
    return importlib.util.resolve_name(
        "." * node.level + (node.module or ""),
        "epcsaft.regression",
    )


def _call_target(node: ast.expr, aliases: dict[str, str]) -> str | None:
    if isinstance(node, ast.Name):
        return aliases.get(node.id, node.id)
    if isinstance(node, ast.Attribute):
        parent = _call_target(node.value, aliases)
        return None if parent is None else f"{parent}.{node.attr}"
    return None


def validate_regression_discovery_edges(regression_root: Path) -> None:
    for source in sorted(regression_root.rglob("*.py")):
        tree = ast.parse(source.read_text(encoding="utf-8"), filename=str(source))
        aliases: dict[str, str] = {}
        for node in ast.walk(tree):
            targets: list[tuple[str, str]] = []
            if isinstance(node, ast.Import):
                for alias in node.names:
                    binding = alias.asname or alias.name.split(".", 1)[0]
                    aliases[binding] = alias.name if alias.asname else binding
                    targets.append((binding, alias.name))
            elif isinstance(node, ast.ImportFrom):
                module = _import_module(node)
                for alias in node.names:
                    target = f"{module}.{alias.name}"
                    aliases[alias.asname or alias.name] = target
                    targets.append((alias.asname or alias.name, target))
            for _binding, target in targets:
                family = _regression_discovery_family(target)
                if family is not None:
                    raise ArchitectureError(
                        f"Regression source imports forbidden {family}: {source}:{node.lineno}"
                    )
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            target = _call_target(node.func, aliases)
            family = None if target is None else _regression_discovery_family(target)
            if family is not None:
                raise ArchitectureError(
                    f"Regression source calls forbidden {family}: {source}:{node.lineno}"
                )


def rejection_code(check: Callable[[], None]) -> int:
    try:
        check()
    except ArchitectureError:
        return 1
    return 0


def run_regression_falsification_gate() -> None:
    fixtures = {
        "PhaseEquilibriumProblem": (
            "from epcsaft.equilibrium import PhaseEquilibriumProblem as Problem\n"
            "Problem(None)\n"
        ),
        "public equilibrium.solve": (
            "from epcsaft import equilibrium as eq\neq.solve(None, None)\n"
        ),
        "phase_equilibrium": (
            "from epcsaft.equilibrium import phase_equilibrium as flash\nflash(None)\n"
        ),
        "HELD/HELD2 controller": (
            "from epcsaft.equilibrium._api import _held_discovery as held\nheld(None)\n"
        ),
        "discovery owner": (
            "from epcsaft.equilibrium._api import _discovery_result as discover\n"
            "discover(None)\n"
        ),
    }
    with tempfile.TemporaryDirectory(prefix="epcsaft-regression-architecture-") as raw:
        root = Path(raw)
        for index, (family, source) in enumerate(fixtures.items()):
            fixture_root = root / str(index)
            fixture_root.mkdir()
            (fixture_root / "fixture.py").write_text(source, encoding="utf-8")
            require(
                rejection_code(
                    lambda root=fixture_root: validate_regression_discovery_edges(root)
                )
                == 1,
                f"temporary {family} falsifier did not fail",
            )
    print("Regression GREPE-only falsification gate passed: 5 discovery edges rejected")


def run_falsification_gate(verifier: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="epcsaft-architecture-") as raw:
        root = Path(raw)
        targets_file = root / "configured-build.json"
        fixture_targets = [
            {
                "name": name,
                "type": "static library",
                "subproject": None,
                "filename": [str(root / filename)],
            }
            for name, filename in STATIC_TARGETS.items()
        ]
        fixture_targets.extend(
            [
                {
                    "name": "_core.so",
                    "type": "shared module",
                    "subproject": None,
                    "filename": [str(root / "_core.so")],
                },
                {
                    "name": "unexpected_native_product",
                    "type": "static library",
                    "subproject": None,
                    "filename": [str(root / "libunexpected.a")],
                },
            ]
        )
        fixture = {
            "targets": fixture_targets,
            "install_plan": {
                "targets": {
                    str(root / "_core.so"): {
                        "destination": "{py_platlib}/epcsaft/_core.so",
                        "subproject": None,
                    }
                }
            },
        }
        targets_file.write_text(json.dumps(fixture), encoding="utf-8")
        configured_build = json.loads(targets_file.read_text(encoding="utf-8"))
        require(
            rejection_code(
                lambda: validate_products(
                    configured_build["targets"], configured_build["install_plan"]
                )
            )
            == 1,
            "temporary extra-native-product falsifier did not fail",
        )

        reverse_link = root / "link-command.txt"
        reverse_link.write_text(
            "c++ libepcsaft_engine_equilibrium_core.a libepcsaft_engine_eos_core.a "
            "libepcsaft_engine_regression_core.a -o _core.so",
            encoding="utf-8",
        )
        require(
            rejection_code(
                lambda: validate_link_command(reverse_link.read_text(encoding="utf-8"))
            )
            == 1,
            "temporary reverse-link falsifier did not fail",
        )

        symbols_file = root / "symbols.txt"
        symbols_file.write_text(
            "PyObject_Call U\nepcsaft_native_transport U\n", encoding="utf-8"
        )
        require(
            rejection_code(
                lambda: validate_symbols(
                    [
                        tuple(line.split())
                        for line in symbols_file.read_text(
                            encoding="utf-8"
                        ).splitlines()
                    ],
                    "fixture",
                )
            )
            == 1,
            "temporary Python/private-ABI symbol falsifier did not fail",
        )

        valid_edges = {
            "epcsaft_engine_eos_core": [("eos_api", "T")],
            "epcsaft_engine_equilibrium_core": [
                ("eos_api", "U"),
                ("equilibrium_api", "T"),
            ],
            "epcsaft_engine_regression_core": [
                ("eos_api", "U"),
                ("equilibrium_api", "U"),
            ],
        }
        forbidden_edge = {owner: fields.copy() for owner, fields in valid_edges.items()}
        forbidden_edge["epcsaft_engine_eos_core"].append(("equilibrium_api", "U"))
        require(
            rejection_code(lambda: validate_sibling_edges(forbidden_edge)) == 1,
            "temporary forbidden sibling-edge falsifier did not fail",
        )
        missing_edge = {owner: fields.copy() for owner, fields in valid_edges.items()}
        missing_edge["epcsaft_engine_equilibrium_core"].remove(("eos_api", "U"))
        require(
            rejection_code(lambda: validate_sibling_edges(missing_edge)) == 1,
            "temporary missing sibling-edge falsifier did not fail",
        )

        wheel = root / "second-extension.whl"
        with ZipFile(wheel, "w") as archive:
            archive.writestr("epcsaft/_core.first.so", b"")
            archive.writestr("epcsaft/_core.second.so", b"")
        completed = subprocess.run(
            [sys.executable, str(verifier), "--wheel", str(wheel)],
            check=False,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        require(
            completed.returncode != 0,
            "temporary second-extension wheel falsifier did not fail",
        )

    print(
        "temporary architecture falsification gate passed: 6 seeded violations rejected"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--wheel", type=Path)
    arguments = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    validate_regression_discovery_edges(repo_root / "engine/src/epcsaft/regression")
    run_regression_falsification_gate()
    require(
        (arguments.build_dir is None) == (arguments.wheel is None),
        "--build-dir and --wheel must be supplied together",
    )
    if arguments.build_dir is None:
        print("Regression production modules import no phase-discovery owner")
        return

    build_dir = arguments.build_dir.resolve()
    wheel = arguments.wheel.resolve()
    require(
        (build_dir / "meson-info" / "meson-info.json").is_file(),
        f"not a configured Meson build: {build_dir}",
    )
    require(wheel.is_file(), f"Engine wheel does not exist: {wheel}")

    verifier = repo_root / "tools" / "verify-installed-wheel.py"
    meson = str(repo_root / ".venv" / "bin" / "meson")
    if not Path(meson).is_file():
        meson = shutil.which("meson") or "meson"
    ninja = shutil.which("ninja") or "ninja"
    nm = shutil.which("nm") or "nm"

    run_falsification_gate(verifier)
    targets = introspect(meson, build_dir, "--targets")
    dependencies = introspect(meson, build_dir, "--dependencies")
    install_plan = introspect(meson, build_dir, "--install-plan")
    require(isinstance(targets, list), "Meson targets output is not a list")
    require(isinstance(dependencies, list), "Meson dependencies output is not a list")
    require(
        isinstance(install_plan, dict), "Meson install-plan output is not a mapping"
    )

    owners = validate_products(targets, install_plan)
    archive_symbols: dict[str, list[tuple[str, str]]] = {}
    for name, expected_filename in STATIC_TARGETS.items():
        filenames = [Path(str(item)) for item in owners[name].get("filename", [])]
        require(
            len(filenames) == 1 and filenames[0].name == expected_filename,
            f"unexpected {name} output: {filenames}",
        )
        require(
            filenames[0].is_file(),
            f"configured native archive is not built: {filenames[0]}",
        )
        command_for_output(ninja, build_dir, filenames[0])
        archive_symbols[name] = symbol_fields(nm, filenames[0])
        validate_symbols(archive_symbols[name], name)

    intersections = validate_sibling_edges(archive_symbols)

    core_files = [Path(str(item)) for item in owners["_core"].get("filename", [])]
    require(len(core_files) == 1, f"unexpected _core outputs: {core_files}")
    core_command = command_for_output(ninja, build_dir, core_files[0])
    validate_link_command(core_command)

    defined_in = {
        Path(str(owners[name]["defined_in"])).resolve()
        for name in (*STATIC_TARGETS, "_core")
    }
    require(
        len(defined_in) == 1,
        f"root native targets have competing Meson owners: {sorted(map(str, defined_in))}",
    )
    source_root = defined_in.pop().parent
    validate_sources(
        production_sources(
            (owners[name] for name in (*STATIC_TARGETS, "_core")), source_root
        )
    )

    print(run([sys.executable, str(verifier), "--wheel", str(wheel)]).strip())
    print("root native products: " + ", ".join((*STATIC_TARGETS, "epcsaft._core")))
    print(
        "internal links: equilibrium -> EOS; regression -> equilibrium + EOS; _core -> all three"
    )
    print(
        "sibling archive symbol intersections: "
        + ", ".join(
            f"{consumer} -> {eos}={len(symbols)}"
            for (consumer, eos), symbols in intersections.items()
        )
    )
    print(
        "external dependencies: "
        + ", ".join(sorted(str(item.get("name")) for item in dependencies))
    )
    print(f"wheel sha256: {hashlib.sha256(wheel.read_bytes()).hexdigest()}")


if __name__ == "__main__":
    try:
        main()
    except ArchitectureError as error:
        raise SystemExit(str(error)) from error

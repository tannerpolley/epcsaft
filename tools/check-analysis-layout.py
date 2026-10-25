#!/usr/bin/env python3
"""Reject displaced owners and malformed analysis-set roots."""

import ast
import gzip
import hashlib
import importlib
import json
import re
from collections import defaultdict
from pathlib import Path
import subprocess
from unicodedata import normalize
table_rows = importlib.import_module("check-notebook-numbers").table_rows

ROOT = Path(__file__).resolve().parents[1]
ANALYSES = ROOT / "analyses"
ACTIVE_ROOTS = (ANALYSES, ROOT / "engine", ROOT / "tools")
COMPACT_SUMMARY_LIMIT = 1_048_576
# Every retained analysis file stays below 32 MiB; the largest (a gzip trace) is 27.6 MB.
FILE_SIZE_LIMIT = 33_554_432
# Duplicates are compared on decompressed content; smaller files are headers and option stubs.
DUPLICATE_FLOOR = 1024
CODE_SUFFIXES = frozenset({".py", ".sh", ".cpp", ".hpp", ".h", ".c", ".jl", ".R", ".m"})
LARGE_RESULT_PAYLOAD_LIMIT = 65_536
EMBEDDED_INPUT_KEYS = frozenset(
    {
        "input",
        "input_specification",
        "specification",
        "observations",
        "starts",
        "evaluated_rows",
    }
)
TOKEN = "re" + "ceipt"
ALLOWED = {
    "README.md",
    "notebook.qmd",
    "inputs",
    "scripts",
    "results",
    "figures",
    "tests",
}


def scipy_optimize_lines(source: str) -> set[int]:
    """Find imports and uses of every scipy.optimize name, including aliases."""
    tree = ast.parse(source)
    # Python normalizes identifiers, so also retain imports spelled with Unicode.
    if "scipy" not in (source if source.isascii() else normalize("NFKC", source)):
        return set()
    aliases = {}
    hits = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for name in node.names:
                aliases[name.asname or name.name.split(".")[0]] = (
                    name.name if name.asname else name.name.split(".")[0]
                )
                if name.name == "scipy.optimize" or name.name.startswith("scipy.optimize."):
                    hits.add(name.lineno)
        elif isinstance(node, ast.ImportFrom):
            for name in node.names:
                module = f"{node.module}.{name.name}"
                aliases[name.asname or name.name] = module
                if module == "scipy.optimize" or module.startswith("scipy.optimize."):
                    hits.add(name.lineno)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Name, ast.Attribute)) and isinstance(node.ctx, ast.Load):
            parts = ast.unparse(node).split(".")
            module = ".".join([aliases.get(parts[0], parts[0]), *parts[1:]])
            if module == "scipy.optimize" or module.startswith("scipy.optimize."):
                hits.add(node.lineno)
    return hits


def check_native_calculation_rule() -> None:
    actual = {
        (str(path.relative_to(ROOT)), line)
        for base in (ROOT / "engine/src", ROOT / "engine/tests", ANALYSES)
        for path in sorted(base.rglob("*.py"))
        for line in scipy_optimize_lines(path.read_text(encoding="utf-8"))
    }
    if actual:
        raise SystemExit("\n".join(
            f"new scipy.optimize use: {path}:{line}" for path, line in sorted(actual)
        ))
    print("native_calculation_rule: no scipy.optimize use")


def check_leaf(path: Path) -> None:
    unexpected = sorted(item.name for item in path.iterdir() if item.name not in ALLOWED)
    if unexpected:
        raise SystemExit(f"unexpected analysis entries in {path.relative_to(ROOT)}: {unexpected}")
    if any(item.is_symlink() for item in path.rglob("*")):
        raise SystemExit(f"analysis leaf contains a symlink: {path.relative_to(ROOT)}")
    for manifest_path in path.rglob("materialization-v*.json"):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        materialized_root = manifest_path.parent / str(manifest["packet_version"])
        for record in manifest["files"]:
            materialized = materialized_root / record["path"]
            if not materialized.is_file():
                raise SystemExit(f"missing materialized file: {materialized.relative_to(ROOT)}")
            actual = hashlib.sha256(materialized.read_bytes()).hexdigest()
            if actual != record["sha256"]:
                raise SystemExit(f"materialized file hash mismatch: {materialized.relative_to(ROOT)}")


def has_nested_result(value: object) -> bool:
    if isinstance(value, dict):
        artifact_results = value.get("current_artifact_results")
        if isinstance(artifact_results, list) and any(
            isinstance(item, dict) and "result" in item for item in artifact_results
        ):
            return True
        return any(has_nested_result(item) for item in value.values())
    if isinstance(value, list):
        return any(has_nested_result(item) for item in value)
    return False


def check_nested_result_shape() -> None:
    duplicate = {"current_artifact_results": [{"result": {}}]}
    metadata_only = {"current_artifact_results": [{"artifact_sha256": "abc"}]}
    if not has_nested_result(duplicate) or has_nested_result(metadata_only):
        raise RuntimeError("Khudaida nested-result shape falsification failed")


def check_result_ownership(path: Path) -> None:
    results = path / "results"
    if not results.is_dir():
        return
    large_payloads: dict[str, list[Path]] = defaultdict(list)
    for item in results.rglob("*.json"):
        try:
            payload = json.loads(item.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if not isinstance(payload, dict):
            continue
        if "input_sha256" in payload and EMBEDDED_INPUT_KEYS.intersection(payload):
            raise SystemExit(
                f"result embeds retained input beside input digest: {item.relative_to(ROOT)}"
            )
        for value in payload.values():
            encoded = json.dumps(value, sort_keys=True, separators=(",", ":"))
            if len(encoded) >= LARGE_RESULT_PAYLOAD_LIMIT:
                large_payloads[encoded].append(item)
    for files in large_payloads.values():
        if len(files) > 1:
            names = ", ".join(str(item.relative_to(ROOT)) for item in files)
            raise SystemExit(f"repeated large result payload: {names}")


def leaf_parts(relative: str) -> tuple[str, tuple[str, ...]]:
    """Split an analyses-relative path into its leaf and the path inside that leaf."""
    parts = tuple(relative.split("/"))
    depth = 2 if parts[0] == "reproduction" else 1
    return "/".join(parts[:depth]), parts[depth:]


def duplicate_groups(names: list[str], read) -> list[list[str]]:
    """Copies across leaves, or within a leaf's results/ and figures/, by decompressed content."""
    by_hash: dict[str, list[str]] = defaultdict(list)
    for name in names:
        # Read one file at a time so peak memory stays near the largest file.
        data = read(name)
        if name.endswith(".gz"):
            data = gzip.decompress(data)
        if len(data) >= DUPLICATE_FLOOR:
            by_hash[hashlib.sha256(data).hexdigest()].append(name)
    groups = []
    for names in by_hash.values():
        located = [leaf_parts(name) for name in names]
        owned = [inner for _, inner in located if inner[0] in ("results", "figures")]
        if len({leaf for leaf, _ in located}) > 1 or len(owned) > 1:
            groups.append(sorted(names))
    return groups


def misplaced_code(names: list[str]) -> list[str]:
    """Analysis code lives in the leaf's scripts/ or in figures/<figure>/scripts/."""
    def placed(inner: tuple[str, ...]) -> bool:
        return inner[0] == "scripts" or (inner[0] == "figures" and len(inner) > 3 and inner[2] == "scripts")
    return [name for name in names if Path(name).suffix in CODE_SUFFIXES and not placed(leaf_parts(name)[1])]


def check_analysis_files() -> None:
    gz = gzip.compress(b"x" * DUPLICATE_FLOOR, mtime=0)
    planted = {
        "a/results/t.csv": b"x" * DUPLICATE_FLOOR, "reproduction/b/inputs/t.csv.gz": gz,
        "c/figures/one.png": b"y" * DUPLICATE_FLOOR, "c/results/two.png": b"y" * DUPLICATE_FLOOR,
        "d/inputs/one.csv": b"z" * DUPLICATE_FLOOR, "d/inputs/two.csv": b"z" * DUPLICATE_FLOOR,
        "e/results/small.csv": b"h", "f/results/small.csv": b"h",
    }
    expected = [["a/results/t.csv", "reproduction/b/inputs/t.csv.gz"], ["c/figures/one.png", "c/results/two.png"]]
    if sorted(duplicate_groups(list(planted), planted.__getitem__)) != expected:
        raise RuntimeError("duplicate-content check missed a planted copy or flagged an allowed one")
    if misplaced_code(["a/results/x/run.py", "reproduction/b/inputs/x.py", "a/scripts/run.py",
                       "a/figures/figure-01/scripts/plot.py", "a/results/x.csv",
                       "a/results/x/scripts/run.py", "a/figures/scripts/plot.py"]) != [
            "a/results/x/run.py", "reproduction/b/inputs/x.py",
            "a/results/x/scripts/run.py", "a/figures/scripts/plot.py"]:
        raise RuntimeError("code-location check missed planted code outside scripts/")

    listed = subprocess.run(["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z", "analyses"],
                            cwd=ROOT, check=True, capture_output=True).stdout.decode().split("\0")
    files = {name.removeprefix("analyses/"): ROOT / name for name in listed if name and (ROOT / name).is_file()}
    oversized = [f"{name} ({path.stat().st_size} bytes)" for name, path in sorted(files.items())
                 if path.stat().st_size >= FILE_SIZE_LIMIT]
    if oversized:
        raise SystemExit(f"analysis files at or above {FILE_SIZE_LIMIT} bytes: {', '.join(oversized)}")
    misplaced = misplaced_code([name for name in files if "/" in name])
    if misplaced:
        raise SystemExit(f"analysis code outside scripts/: {', '.join(misplaced)}")
    for group in duplicate_groups([name for name in files if "/" in name], lambda name: files[name].read_bytes()):
        raise SystemExit(f"duplicate analysis content: {', '.join(group)}")


def index_errors(text: str, delivered: dict[str, bool]) -> list[str]:
    """Each leaf has one index row, each row names a leaf, and a delivered leaf is not pending."""
    rows, errors = [], []
    for number, cells in table_rows(text):
        link = re.fullmatch(r"\[[^\]]+\]\(([^)]+?)/?\)", cells[0])
        status = next((cell for cell in reversed(cells) if cell), "")
        if link and status:
            rows.append((link[1], status))
        elif cells[0] != "Leaf":
            errors.append(f"unreadable index row at line {number}")
    names = [name for name, _ in rows]
    return (errors + [f"leaf without one index row: {name}" for name in delivered if names.count(name) != 1]
            + [f"index row without a leaf: {name}" for name in names if name not in delivered]
            + [f"delivered leaf marked pending: {name}" for name, status in rows
               if status == "pending" and delivered.get(name)])


def check_index(leaves: list[Path]) -> None:
    planted = "| Leaf | Q | Status |\n|---|---|---|\n| [a](a/) | Q | pending |\n| [c](c/) | Q | current |\n"
    if index_errors(planted, {"a": True, "b": False}) != [
            "leaf without one index row: b", "index row without a leaf: c", "delivered leaf marked pending: a"]:
        raise RuntimeError("index check missed a planted row error")
    delivered = {str(leaf.relative_to(ANALYSES)): (leaf / "notebook.qmd").is_file() and (leaf / "results").is_dir()
                 for leaf in leaves}
    errors = index_errors((ANALYSES / "README.md").read_text(encoding="utf-8"), delivered)
    if errors:
        raise SystemExit("analyses/README.md index: " + "; ".join(errors))


def check_active_paths() -> None:
    for base in ACTIVE_ROOTS:
        for path in base.rglob("*"):
            if TOKEN in path.name.lower():
                raise SystemExit(f"receipt-named path: {path.relative_to(ROOT)}")

    active_files = [
        path for path in ANALYSES.rglob("*")
        if path.is_file() and "scripts" in path.parts and "__pycache__" not in path.parts
    ]
    active_files.extend(
        path for base in (ROOT / "engine/src", ROOT / "engine/native", ROOT / "engine/scripts", ROOT / "tools")
        for path in base.rglob("*")
        if path.is_file() and "__pycache__" not in path.parts
    )
    for path in active_files:
        if path == Path(__file__).resolve():
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        if TOKEN in text.lower():
            raise SystemExit(f"receipt terminology in active source: {path.relative_to(ROOT)}")

    for path in ANALYSES.rglob("summary.json"):
        if path.stat().st_size > COMPACT_SUMMARY_LIMIT:
            raise SystemExit(f"summary.json exceeds 1 MiB: {path.relative_to(ROOT)}")

    for path in ANALYSES.rglob("*.json"):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        siblings = (
            tuple(path.parent.glob("*.csv"))
            + tuple(path.parent.glob("*.tsv"))
            + tuple(path.parent.glob("data/*.csv"))
            + tuple(path.parent.glob("data/*.tsv"))
            + tuple(path.parent.parent.glob("data/*.csv"))
            + tuple(path.parent.parent.glob("data/*.tsv"))
            + tuple(path.parent.parent.glob("results/data/*.csv"))
            + tuple(path.parent.parent.glob("results/data/*.tsv"))
        )
        if isinstance(payload, dict) and "calls" in payload and (path.parent / "calls.csv").is_file():
            raise SystemExit(f"duplicated fixed-pressure calls: {path.relative_to(ROOT)}")
        if has_nested_result(payload) and any(
            "result_json" in sibling.read_text(encoding="utf-8", errors="ignore")
            for sibling in siblings
        ):
            raise SystemExit(f"duplicated Khudaida results: {path.relative_to(ROOT)}")


def removed_interfaces(source: str, exports: dict[str, set[str]]) -> set[str]:
    tree = ast.parse(source)
    aliases = {}
    removed = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for name in node.names:
                if name.name in exports:
                    aliases[name.asname or name.name.split(".")[0]] = (
                        name.name if name.asname else name.name.split(".")[0]
                    )
                elif name.name.startswith("epcsaft."):
                    path = ROOT / "engine/src" / name.name.replace(".", "/")
                    if not path.is_dir() and not path.with_suffix(".py").is_file():
                        removed.add(name.name)
        elif isinstance(node, ast.ImportFrom):
            for name in node.names:
                module = f"{node.module}.{name.name}"
                if module in exports:
                    aliases[name.asname or name.name] = module
                elif node.module in exports:
                    if name.name != "*" and name.name not in exports[node.module]:
                        removed.add(module)
                elif node.module and node.module.startswith("epcsaft."):
                    path = ROOT / "engine/src" / node.module.replace(".", "/")
                    if not path.is_dir() and not path.with_suffix(".py").is_file():
                        removed.add(node.module)
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute):
            parts = ast.unparse(node).split(".")
            module = ".".join([aliases.get(parts[0], ""), *parts[1:-1]])
            if module in exports and not node.attr.startswith("_") and node.attr not in exports[module]:
                removed.add(f"{module}.{node.attr}")
    return removed


def check_current_commands() -> None:
    # Read exported names from source; this check never imports the installed wheel.
    exports = {}
    for module in ("epcsaft", "epcsaft.equilibrium"):
        path = ROOT / "engine/src" / module.replace(".", "/")
        tree = ast.parse((path / "__init__.py").read_text(encoding="utf-8"))
        names = next(
            ast.literal_eval(node.value) for node in tree.body
            if isinstance(node, ast.Assign)
            and any(isinstance(target, ast.Name) and target.id == "__all__" for target in node.targets)
        )
        exports[module] = set(names) | {item.stem for item in path.iterdir()}
    for source in (
        "from epcsaft import unit_registry as u",
        "import epcsaft as eos; eos.unit_registry",
        "from epcsaft import equilibrium as eq; eq.GeneralReactiveEquilibriumProblem",
        "import epcsaft.equilibrium; epcsaft.equilibrium.Fixed",
        "from epcsaft.equilibrium import Fixed",
        "from epcsaft.units import ParameterError",
        "import epcsaft.units",
    ):
        if not removed_interfaces(source, exports):
            raise RuntimeError(f"retired-interface check missed: {source}")
    if removed_interfaces("from epcsaft import Mixture, equilibrium; equilibrium.Problem", exports):
        raise RuntimeError("retired-interface check rejected current operations")
    checked = {}
    for base in (ANALYSES, ROOT / "docs", ROOT / "engine/docs"):
        for document in base.rglob("*"):
            if document.suffix not in {".md", ".qmd"}:
                continue
            text = document.read_text(encoding="utf-8")
            # Pair every fence, then keep the shell ones, so a preceding non-shell block cannot shift the pairing.
            for info, block in re.findall(r"```([^\n`]*)\n(.*?)```", text, re.DOTALL):
                if info.strip(" {}").split(",")[0] not in {"", "bash", "sh", "shell", "console"}:
                    continue
                for name in re.findall(r"[\w./-]+\.py\b", block):
                    path = ROOT / name
                    if not path.is_file():
                        path = document.parent / name
                    if path.is_file():
                        if path not in checked:
                            checked[path] = removed_interfaces(path.read_text(encoding="utf-8"), exports)
                        removed = checked[path]
                        if removed:
                            raise SystemExit(
                                f"current command in {document.relative_to(ROOT)} uses retired interfaces: "
                                f"{path.relative_to(ROOT)}: {', '.join(sorted(removed))}"
                            )


for retired in (ROOT / "governance", ROOT / "validation"):
    if retired.exists():
        raise SystemExit(f"retired owner still exists: {retired.relative_to(ROOT)}")

check_native_calculation_rule()
check_nested_result_shape()
check_analysis_files()
check_active_paths()
check_current_commands()

reproduction = ANALYSES / "reproduction"
if not reproduction.is_dir() or any(item.is_file() for item in reproduction.iterdir()):
    raise SystemExit("analyses/reproduction must contain only analysis directories")

unexpected_root = sorted(
    item.name
    for item in ANALYSES.iterdir()
    if item.is_file() and item.name not in {"README.md", "_quarto.yml", "index.qmd", ".gitignore"}
)
if unexpected_root:
    raise SystemExit(f"unexpected files directly under analyses: {unexpected_root}")
if any(item.is_symlink() for item in ANALYSES.iterdir()):
    raise SystemExit("analyses must not contain root-level symlinks")

leaves = [item for item in ANALYSES.iterdir() if item.is_dir() and item != reproduction and item.name != ".quarto"]
leaves.extend(item for item in reproduction.iterdir() if item.is_dir())
for leaf in leaves:
    check_leaf(leaf)
    check_result_ownership(leaf)
check_index(leaves)

print(f"analysis_layout: clean ({len(leaves)} leaves)")

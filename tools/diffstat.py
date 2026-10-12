#!/usr/bin/env python3
"""Runtime lines per layer, as laid out by one pinned formatter, against ARCHITECTURE.yaml code budgets.

    tools/diffstat.py --budgets        # current full runtime tree; a ceiling above or below its band fails
"""
from __future__ import annotations
import io, sys
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
CPP_STYLE = TOOLS / "line-count.clang-format"
RUFF = TOOLS.parent / ".venv/bin/ruff"


def band(count: int) -> int:
    return -(-count // 100)


def budget_state(count: int, ceiling: int) -> str:
    if count > ceiling:
        return f"FAIL over the ceiling by {count - ceiling}"
    if ceiling - count > 2 * band(count):
        return f"FAIL headroom exceeds two bands: lower the ceiling to {count + band(count)}"
    return "ok"


def formatted_lines(path: Path) -> int:
    """Physical lines of the file as the pinned formatter would lay it out; the file is not changed."""
    import subprocess
    if path.suffix in {".py", ".pyx", ".pxd"}:
        if not RUFF.is_file():
            raise ValueError(f"missing {RUFF}: run bash .codex/environments/setup.sh")
        # A trailing comma does not choose the layout, as in clang-format; 120 columns in both languages.
        command = [str(RUFF), "format", "--isolated", "--line-length=120", "--target-version=py313",
                   "--config=format.skip-magic-trailing-comma=true", "--config=lint.isort.split-on-trailing-comma=false",
                   f"--stdin-filename={path.name}", "-"]
    else:
        command = ["clang-format-18", f"--style=file:{CPP_STYLE}", f"--assume-filename={path.name}"]
    done = subprocess.run(command, input=path.read_bytes(), capture_output=True, check=False)
    if done.returncode:
        raise ValueError(f"cannot format {path} for counting: {done.stderr.decode(errors='replace').strip()}")
    return len(done.stdout.decode("utf-8").splitlines())


def code_budgets(root=None) -> int:
    """Count the complete runtime tree; budgets have one owner in ARCHITECTURE."""
    import json
    import re
    from concurrent.futures import ThreadPoolExecutor
    root = TOOLS.parent if root is None else Path(root)
    architecture = (root / "ARCHITECTURE.yaml").read_text()
    sections = list(re.finditer(r"^code_budgets:", architecture, re.MULTILINE))
    if len(sections) != 1:
        raise ValueError("expected one code_budgets JSON flow mapping in ARCHITECTURE.yaml")
    budgets, _ = json.JSONDecoder().raw_decode(architecture[sections[0].end():].lstrip())
    expected = {"core", "contributions", "properties", "equilibrium", "regression", "surface", "bindings", "python", "total"}
    if set(budgets) != expected:
        raise ValueError("code budgets must classify every runtime layer and total")
    if any(type(n) is not int or n < 0 for n in budgets.values()):
        raise ValueError("each code budget is one non-negative integer ceiling")
    if budgets["total"] > sum(n for k, n in budgets.items() if k != "total"):
        raise ValueError("total budget must not exceed the sum of layer budgets")
    files = []
    suffixes = {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx", ".inc", ".inl", ".tpp", ".py", ".pyx", ".pxd"}
    for folder in (root / "engine/native", root / "engine/src/epcsaft"):
        if folder.is_symlink():
            raise ValueError(f"runtime root symlink bypasses source accounting: {folder}")
        if not folder.is_dir():
            raise ValueError(f"missing runtime tree: {folder}")
        for path in sorted(folder.rglob("*")):
            if path.is_symlink():
                raise ValueError(f"runtime symlink bypasses source accounting: {path}")
            if not path.is_file() or (path.suffix == ".pyc" and "__pycache__" in path.parts):
                continue
            relative = path.relative_to(folder)
            layer = "python" if folder.name == "epcsaft" else relative.parts[0]
            if layer not in expected - {"total", "python"} and folder.name != "epcsaft":
                raise ValueError(f"unclassified native source: {path}")
            if path.suffix not in suffixes:
                # The extension's type stub declares native types; the architecture check owns its names.
                if path.name not in {".gitkeep", "py.typed", "AGENTS.md", "README.md", "_core.pyi"}:
                    raise ValueError(f"unclassified runtime file: {path}")
                continue
            files.append((layer, path))
    counts = dict.fromkeys(expected, 0)
    with ThreadPoolExecutor() as pool:
        for (layer, _), lines in zip(files, pool.map(formatted_lines, [path for _, path in files])):
            counts[layer] += lines
    counts["total"] = sum(counts.values())
    failed = False
    print("runtime formatted lines: count / ceiling (band = 1% of count, rounded up)")
    for layer, ceiling in budgets.items():
        state = budget_state(counts[layer], ceiling)
        print(f"{counts[layer]:5} / {ceiling:5} (band {band(counts[layer]):3})  {layer:13} {state}")
        failed |= state != "ok"
    return int(failed)


def budget_self_check() -> None:
    """One executable check of both band boundaries, formatter-normalized and untracked-tree accounting."""
    import contextlib
    import tempfile
    import json
    layers = ["contributions", "properties", "equilibrium", "regression", "surface", "bindings"]
    with tempfile.TemporaryDirectory(prefix="epcsaft-budget-") as temporary:
        root = Path(temporary)
        (root / "engine/native/core").mkdir(parents=True)
        (root / "engine/src/epcsaft").mkdir(parents=True)
        (root / "engine/src/epcsaft/packed.py").write_text("x = 1; y = 2\n")
        limits = {"core": 102, "python": 2, "total": 104} | dict.fromkeys(layers, 0)
        (root / "ARCHITECTURE.yaml").write_text("code_budgets: " + json.dumps(limits))
        source = root / "engine/native/core/untracked.cpp"
        with contextlib.redirect_stdout(io.StringIO()):
            # Packed statements count as the formatter lays them out: two per line here, two lines each.
            for statements, fails in ((99, True), (100, False), (102, False), (103, True)):
                source.write_text("int a; int b;\n" * (statements // 2) + "int c;\n" * (statements % 2))
                if code_budgets(root) != int(fails):
                    raise AssertionError(f"budget band/accounting failure at {statements} statements")
            source.unlink()
            source = root / "engine/native/unknown/new.cpp"
            source.parent.mkdir()
            source.write_text("int hidden;\n")
            try:
                code_budgets(root)
            except ValueError as error:
                if "unclassified native source" not in str(error):
                    raise
            else:
                raise AssertionError("unclassified source escaped the budget")


def main() -> int:
    if sys.argv[1:] != ["--budgets"]:
        raise SystemExit("usage: tools/diffstat.py --budgets")
    budget_self_check()
    return code_budgets()


if __name__ == "__main__":
    raise SystemExit(main())

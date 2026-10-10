#!/usr/bin/env python3
"""Runtime physical lines per layer against the ARCHITECTURE.yaml code budgets.

    tools/diffstat.py --budgets        # current full runtime tree; a hard ceiling fails
"""
from __future__ import annotations
import io, sys


def budget_state(count: int, limits: list[int]) -> str:
    if count > limits[2]:
        return "FAIL"
    return "WARN" if count > limits[1] else "ok"


def code_budgets(root=None) -> int:
    """Count the complete runtime tree; budgets have one owner in ARCHITECTURE."""
    from pathlib import Path
    import json
    import re
    root = Path(__file__).resolve().parent.parent if root is None else Path(root)
    architecture = (root / "ARCHITECTURE.yaml").read_text()
    sections = list(re.finditer(r"^code_budgets:", architecture, re.MULTILINE))
    if len(sections) != 1:
        raise ValueError("expected one code_budgets JSON flow mapping in ARCHITECTURE.yaml")
    budgets, _ = json.JSONDecoder().raw_decode(architecture[sections[0].end():].lstrip())
    expected = {"core", "contributions", "properties", "equilibrium", "regression", "surface", "bindings", "python", "total"}
    if set(budgets) != expected:
        raise ValueError("code budgets must classify every runtime layer and total")
    for layer, limits in budgets.items():
        if (not isinstance(limits, list) or len(limits) != 3
                or any(type(n) is not int or n <= 0 for n in limits)
                or not limits[0] <= limits[1] <= limits[2]):
            raise ValueError(f"invalid floor/soft/hard budget for {layer}")
    if any(budgets["total"][i] > sum(v[i] for k, v in budgets.items() if k != "total") for i in range(3)):
        raise ValueError("total budget must not exceed the sum of layer budgets")
    counts = dict.fromkeys(expected, 0)
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
            counts[layer] += len(path.read_text(encoding="utf-8").splitlines())
    counts["total"] = sum(counts.values())
    failed = False
    print("runtime physical lines: current / floor target / soft / hard")
    for layer, limits in budgets.items():
        state = budget_state(counts[layer], limits)
        print(f"{state:4} {layer:13} {counts[layer]:5} / {limits[0]:5} / {limits[1]:5} / {limits[2]:5}")
        failed |= state == "FAIL"
    return int(failed)


def budget_self_check() -> None:
    """One executable check of hard/soft boundaries and untracked-tree accounting."""
    import contextlib
    import tempfile
    import json
    from pathlib import Path
    layers = ["core", "contributions", "properties", "equilibrium", "regression", "surface", "bindings", "python"]
    with tempfile.TemporaryDirectory(prefix="epcsaft-budget-") as temporary:
        root = Path(temporary)
        (root / "engine/native/core").mkdir(parents=True)
        (root / "engine/src/epcsaft").mkdir(parents=True)
        limits = {name: [8, 10, 12] for name in layers}
        limits["total"] = [64, 80, 96]
        (root / "ARCHITECTURE.yaml").write_text('code_budgets: ' + json.dumps(limits))
        source = root / "engine/native/core/untracked.cpp"
        with contextlib.redirect_stdout(io.StringIO()):
            for lines, expected in ((0, "ok"), (10, "ok"), (11, "WARN"), (12, "WARN"), (13, "FAIL")):
                source.write_text("// counted physical line\n" * lines)
                if budget_state(lines, [8, 10, 12]) != expected or code_budgets(root) != int(expected == "FAIL"):
                    raise AssertionError("budget boundary/accounting failure")
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

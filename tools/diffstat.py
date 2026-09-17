#!/usr/bin/env python3
"""What has actually been written, by language and by category.

A raw diff against `main` answers "what changed" — and since Stage 1 deleted the
superseded tree, that answer is dominated by 45,000 lines nobody wrote. The default
baseline here is the commit that emptied the tree, so the number answers "what have we
built".

    tools/diffstat.py                  # since the tree was emptied (default)
    tools/diffstat.py --vs-main        # the whole rebuild, deletions included
    tools/diffstat.py HEAD~5..HEAD     # any range
    tools/diffstat.py --budgets        # current full runtime; hard ceiling fails
    tools/diffstat.py --code           # executable code only
    tools/diffstat.py --files          # every file, classified
    tools/diffstat.py --no-split       # skip the statement/comment/blank pass
    tools/diffstat.py --watch          # live: redraws on any change, uncommitted work included
"""
from __future__ import annotations
import io, subprocess, sys
from collections import defaultdict

EMPTIED = "e4d964c"          # Stage 1: the layer tree laid, the superseded one deleted

LANG = {  # extension -> (language, is_executable_code)
    ".cpp": ("C++", True), ".hpp": ("C++", True), ".h": ("C++", True), ".cc": ("C++", True),
    ".py": ("Python", True), ".pyi": ("Python", True), ".sh": ("Shell", True),
    ".build": ("Build", True), ".toml": ("Config", False), ".cfg": ("Config", False),
    ".json": ("Data", False), ".csv": ("Data", False), ".yaml": ("Data", False),
    ".yml": ("Data", False), ".txt": ("Data", False), ".bib": ("Data", False),
    ".md": ("Prose", False), ".tex": ("Prose", False), ".qmd": ("Prose", False),
    ".svg": ("Asset", False), ".png": ("Asset", False), ".pdf": ("Asset", False),
    ".html": ("Asset", False), ".zip": ("Asset", False), ".gitkeep": ("Asset", False),
}

def language(path: str) -> tuple[str, bool]:
    name = path.rsplit("/", 1)[-1]
    if name in ("meson.build", "meson_options.txt"): return ("Build", True)
    if name == "pyproject.toml": return ("Build", True)
    for ext, v in LANG.items():
        if name.endswith(ext): return v
    return ("Other", False)

def category(p: str) -> str:
    if p.startswith("engine/cpp/") or p.startswith("engine/src/epcsaft/_native/"): return "removed"
    if p.startswith("engine/native/core/"): return "core"
    if p.startswith(("engine/native/surface/", "engine/native/bindings/", "engine/src/")): return "surface"
    if p.startswith("engine/native/"): return "engine"
    if p.startswith(("tools/", "data/tools/", "engine/scripts/")): return "canon"
    if "briefs/" in p: return "canon"
    if p.startswith("engine/tests/"): return "canon"
    if p.startswith(("docs/", "engine/docs/")) or p.endswith((".md", ".tex", ".bib")): return "canon"
    if p.startswith("data/"): return "canon"
    if p.startswith("analyses/"): return "analyses"
    return "build"

LINE_COMMENT = {"C++": ("//",), "Python": ("#",), "Shell": ("#",), "Build": ("#",)}

def split_added(rng: str, path: str, lang: str) -> tuple[int, int, int]:
    """Added lines as (statement, comment, blank).

    No off-the-shelf counter does this for a *diff* — tokei, scc and cloc classify a
    tree. Blanks and comments are real work and stay in the file; they are simply not
    what "how much code" is asking.
    """
    if lang not in LINE_COMMENT:
        return (0, 0, 0)
    out = subprocess.run(["git", "diff", "-U0", rng, "--", path],
                         capture_output=True, text=True).stdout
    code = comment = blank = 0
    in_block = False
    q3 = ('"' * 3, "'" * 3)
    for raw in out.splitlines():
        if not raw.startswith("+") or raw.startswith("+++"):
            continue
        s = raw[1:].strip()
        if lang == "C++":
            if in_block:
                comment += 1
                if "*/" in s:
                    in_block = False
                continue
            if not s:
                blank += 1
            elif s.startswith("//"):
                comment += 1
            elif s.startswith("/*"):
                comment += 1
                in_block = "*/" not in s
            else:
                code += 1
        else:
            if in_block:
                comment += 1
                if any(s.count(q) % 2 for q in q3):
                    in_block = False
                continue
            if not s:
                blank += 1
            elif s.startswith(LINE_COMMENT[lang]):
                comment += 1
            elif s.startswith(q3):
                comment += 1
                in_block = all(s.count(q) < 2 for q in q3)
            else:
                code += 1
    return (code, comment, blank)


def report(rng: str, code_only: bool, show_files: bool, live: bool) -> int:
    out = subprocess.run(["git", "diff", "--numstat", rng], capture_output=True, text=True).stdout
    rows = []
    for line in out.splitlines():
        a, d, path = (line.split("\t") + ["", "", ""])[:3]
        if a == "-": continue                      # binary or -diff in .gitattributes
        lang, is_code = language(path)
        rows.append((int(a), int(d), path, lang, is_code, category(path)))
    if not rows:
        print("  nothing in range"); return 0

    label = f"{rng} + working tree" if ".." not in rng else rng
    print(f"\n  range: {label}\n")
    by_lang: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0])
    by_cat: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0])
    cross: dict[tuple[str, str], int] = defaultdict(int)
    code_add = other_add = 0
    for a, d, path, lang, is_code, cat in rows:
        if code_only and not is_code: continue
        by_lang[lang][0] += a; by_lang[lang][1] += d; by_lang[lang][2] += 1
        by_cat[cat][0] += a;  by_cat[cat][1] += d;  by_cat[cat][2] += 1
        cross[(cat, lang)] += a
        if is_code and cat != "removed": code_add += a
        elif cat != "removed": other_add += a

    print("  BY LANGUAGE");  print(f"    {'':12}{'added':>9}{'deleted':>9}{'files':>7}   executable")
    for lang, (a, d, f) in sorted(by_lang.items(), key=lambda kv: -kv[1][0]):
        is_code = next(v[1] for k, v in LANG.items() if v[0] == lang) if any(v[0] == lang for v in LANG.values()) else False
        if lang in ("C++", "Python", "Shell", "Build"): is_code = True
        print(f"    {lang:12}{a:+9d}{-d:9d}{f:7d}   {'yes' if is_code else '—'}")
    print(f"\n    {'CODE':12}{code_add:+9d}\n    {'not code':12}{other_add:+9d}")

    if not live:
        c = m = b = 0
        for a, d, path, lang, is_code, cat in rows:
            if cat == "removed" or not is_code:
                continue
            cc, mm, bb = split_added(rng, path, lang)
            c += cc; m += mm; b += bb
        if c or m or b:
            tot = c + m + b
            print("\n  OF THOSE CODE ADDITIONS")
            print(f"    {'statements':12}{c:+9d}   {c/tot:6.1%}")
            print(f"    {'comments':12}{m:+9d}   {m/tot:6.1%}")
            print(f"    {'blank':12}{b:+9d}   {b/tot:6.1%}")

    print("\n  BY CATEGORY"); print(f"    {'':12}{'added':>9}{'deleted':>9}{'files':>7}")
    for cat in ("core", "surface", "engine", "build", "canon", "analyses", "removed"):
        if cat not in by_cat: continue
        a, d, f = by_cat[cat]
        print(f"    {cat:12}{a:+9d}{-d:9d}{f:7d}")

    langs = [l for l in ("C++", "Python", "Build", "Data", "Prose") if any(k[1] == l for k in cross)]
    print("\n  ADDITIONS, CATEGORY × LANGUAGE"); print("    " + "".ljust(12) + "".join(l.rjust(9) for l in langs))
    for cat in ("core", "surface", "engine", "build", "canon", "analyses"):
        if cat not in by_cat: continue
        print(f"    {cat:12}" + "".join((f"{cross[(cat,l)]:+9d}" if cross[(cat, l)] else "        ·") for l in langs))

    if show_files:
        print("\n  FILES")
        for a, d, path, lang, is_code, cat in sorted(rows, key=lambda r: -r[0]):
            if code_only and not is_code: continue
            print(f"    {a:+7d}{-d:7d}  {cat:9} {lang:7} {path}")
    return code_add



def fingerprint() -> str:
    """Cheap signature of everything that could change the report."""
    head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout
    dirty = subprocess.run(["git", "status", "--porcelain=v1"], capture_output=True, text=True).stdout
    stat = subprocess.run(["git", "diff", "--numstat", "HEAD"], capture_output=True, text=True).stdout
    return head + dirty + stat


def watch(rng: str, code_only: bool) -> int:
    """Redraw whenever the repository changes. Polls a fingerprint, so no dependency."""
    import time
    last, prev_total = None, None
    try:
        while True:
            now = fingerprint()
            if now != last:
                last = now
                buf = io.StringIO()
                out, sys.stdout = sys.stdout, buf
                try:
                    total = report(rng, code_only, False, True)
                finally:
                    sys.stdout = out
                print("\033[2J\033[H", end="")
                head = subprocess.run(["git", "log", "--oneline", "-1"],
                                      capture_output=True, text=True).stdout.strip()
                dirty = subprocess.run(["git", "status", "--porcelain=v1"],
                                       capture_output=True, text=True).stdout.splitlines()
                delta = ""
                if prev_total is not None and total != prev_total:
                    delta = f"   ({total - prev_total:+d} since last redraw)"
                prev_total = total
                print(f"  \033[1mePC-SAFT — live\033[0m   {time.strftime('%H:%M:%S')}{delta}")
                print(f"  HEAD {head[:66]}")
                print(f"  working tree: {len(dirty)} file(s) modified"
                      if dirty else "  working tree: clean")
                print(buf.getvalue(), end="")
                print("\n  watching… Ctrl-C to stop")
            time.sleep(1.0)
    except KeyboardInterrupt:
        print("\n  stopped")
        return 0


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
    suffixes = {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx", ".inc", ".inl", ".tpp", ".py", ".pyi", ".pyx", ".pxd"}
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
                if path.name not in {".gitkeep", "py.typed", "AGENTS.md", "README.md"}:
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
    args = sys.argv[1:]
    if "--budgets" in args:
        if args != ["--budgets"]:
            raise ValueError("--budgets is a current-tree check; do not combine with a diff range")
        budget_self_check()
        return code_budgets()
    code_only = "--code" in args
    show_files = "--files" in args
    live = "--watch" in args or "--no-split" in args
    vs_main = "--vs-main" in args
    rng = next((a for a in args if not a.startswith("--")), None)
    if rng is None:
        base = subprocess.run(["git", "merge-base", "main", "HEAD"], capture_output=True, text=True).stdout.strip() if vs_main else EMPTIED
        rng = base if "--watch" in args else f"{base}..HEAD"   # bare base includes uncommitted work
    if "--watch" in args:
        return watch(rng, code_only)
    report(rng, code_only, show_files, live)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

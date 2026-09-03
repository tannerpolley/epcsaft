#!/usr/bin/env python3
"""What has actually been written, by language and by category.

A raw diff against `main` answers "what changed" — and since Stage 1 deleted the
superseded tree, that answer is dominated by 45,000 lines nobody wrote. The default
baseline here is the commit that emptied the tree, so the number answers "what have we
built".

    tools/diffstat.py                  # since the tree was emptied (default)
    tools/diffstat.py --vs-main        # the whole rebuild, deletions included
    tools/diffstat.py HEAD~5..HEAD     # any range
    tools/diffstat.py --code           # executable code only
    tools/diffstat.py --files          # every file, classified
"""
from __future__ import annotations
import subprocess, sys
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
    if p.startswith("engine/docs/") or p.endswith((".md", ".tex", ".bib")): return "canon"
    if p.startswith("data/"): return "canon"
    if p.startswith("analyses/"): return "analyses"
    return "build"

def main() -> int:
    args = [a for a in sys.argv[1:]]
    code_only = "--code" in args; show_files = "--files" in args
    vs_main = "--vs-main" in args
    rng = next((a for a in args if not a.startswith("--")), None)
    if rng is None:
        base = subprocess.run(["git", "merge-base", "main", "HEAD"], capture_output=True, text=True).stdout.strip() if vs_main else EMPTIED
        rng = f"{base}..HEAD"
    out = subprocess.run(["git", "diff", "--numstat", rng], capture_output=True, text=True).stdout
    rows = []
    for line in out.splitlines():
        a, d, path = (line.split("\t") + ["", "", ""])[:3]
        if a == "-": continue                      # binary or -diff in .gitattributes
        lang, is_code = language(path)
        rows.append((int(a), int(d), path, lang, is_code, category(path)))
    if not rows:
        print("  nothing in range"); return 0

    print(f"\n  range: {rng}    ({'vs main, deletions included' if vs_main else 'since the tree was emptied'})\n")
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
    return 0

if __name__ == "__main__":
    raise SystemExit(main())

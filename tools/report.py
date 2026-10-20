#!/usr/bin/env python3
"""Write build/report.html from retained records: tests, key-solve work and the analysis index.

Reads files only; never builds, imports the Engine, runs a solve or renders a notebook.
  tools/report [BUILD]      test results from BUILD (default build/); key solves from BUILD/wheels,
                            default the main checkout's build/wheels
  tools/report --stamp XML  write XML.source: the commit and engine/ state that produced XML
"""
import json
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from html import escape
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "build/report.html"


def git(*args, cwd=ROOT):
    return subprocess.run(["git", "-C", str(cwd), *args], capture_output=True, text=True).stdout.strip()


def source():
    """Commit and engine/ state, as epcsaft.provenance() judges an editable install."""
    return git("rev-parse", "HEAD"), "dirty" if git("status", "--porcelain", "--", "engine") else "clean"


def stamp(when, commit, state, head):
    """Flag data from another commit; without a recorded commit, data written before HEAD's commit time."""
    text = f"data: {commit[:8]} engine/ {state}" if commit else "data: commit not recorded"
    stale = commit != head[0] if commit else time.mktime(time.strptime(when, "%Y-%m-%d %H:%M")) < head[1]
    flags = [f for f, bad in (("not HEAD", stale), ("dirty tree", state == "dirty")) if bad]
    return f'<p class="src">{escape(text)}, {escape(when)}' + "".join(f' <b class="flag">{f}</b>' for f in flags) + "</p>"


def mtime(path):
    return time.strftime("%Y-%m-%d %H:%M", time.localtime(path.stat().st_mtime))


def table(head, rows):
    out = "<tr>" + "".join(f"<th>{h}</th>" for h in head) + "</tr>"
    return "<table>" + out + "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows) + "</table>"


def tests(results, head):
    html, cases = [], []
    for name in ("native", "pytest"):
        path = results / f"{name}.xml"
        try:
            root = ET.parse(path).getroot()
        except (OSError, ET.ParseError) as error:
            html.append(f"<h3>{name}.xml</h3><p class='none'>no data ({escape(type(error).__name__)})</p>")
            continue
        written = path.with_suffix(".xml.source")
        commit, state = written.read_text().split()[:2] if written.exists() else ("", "")
        areas = {}
        for case in root.iter("testcase"):
            if name == "native":  # epcsaft:equilibrium/Automatic_fast/f5 -> equilibrium/Automatic_fast
                test = case.get("name", "").removeprefix("epcsaft:")
                area = test.rsplit("/", 1)[0]
            else:  # tests.eos.test_actions -> tests/eos
                test = case.get("classname", "") + "::" + case.get("name", "")
                area = "/".join(case.get("classname", "").split(".")[:-1])
            bad = case.find("failure") if case.find("failure") is not None else case.find("error")
            outcome = "fail" if bad is not None else "skip" if case.find("skipped") is not None else "pass"
            counts = areas.setdefault(area, {"pass": 0, "fail": 0, "skip": 0})
            counts[outcome] += 1
            message = ""
            if bad is not None:
                out = (case.findtext("system-out") or "").strip().splitlines()
                message = bad.get("message") or (bad.text or "").strip() or "\n".join(out[-5:])
            cases.append((float(case.get("time") or 0), name, test, outcome, message))
        total = {k: sum(c[k] for c in areas.values()) for k in ("pass", "fail", "skip")}
        html.append(f"<h3>{name}.xml: {total['pass']} pass, {total['fail']} fail, {total['skip']} skip</h3>")
        html.append(stamp(mtime(path), commit, state, head))
        html.append(table(["Area", "Pass", "Fail", "Skip"], [
            (escape(a), c["pass"], f'<b class="flag">{c["fail"]}</b>' if c["fail"] else 0, c["skip"])
            for a, c in sorted(areas.items())]))
    failures = [c for c in cases if c[3] == "fail"]
    if failures:
        html.append("<h3>Failures</h3>" + table(["File", "Test", "Message"], [
            (n, escape(t), f"<pre>{escape(m[:2000])}</pre>") for _, n, t, _, m in failures]))
    if cases:
        html.append("<h3>Slowest 10</h3>" + table(["Seconds", "File", "Test"], [
            (f"{s:.2f}", n, escape(t)) for s, n, t, *_ in sorted(cases, reverse=True)[:10]]))
    return "".join(html)


def sparkline(values):
    lo, hi = min(values), max(values)
    step = 120 / max(len(values) - 1, 1)
    points = " ".join(f"{i * step:.1f},{22 - 20 * (v - lo) / ((hi - lo) or 1):.1f}" for i, v in enumerate(values))
    return f'<svg width="124" height="24"><polyline points="{points}" fill="none" stroke="#2563eb"/></svg>'


def key_solves(wheels, head):
    runs = []
    for path in wheels.glob("key-solves-*/summary.json"):
        try:
            summary = json.loads(path.read_text())
            if "cases" in summary:  # error records hold no timings
                runs.append((path.parent.name.rsplit("-", 1)[1], summary, path))
        except (OSError, ValueError, IndexError):
            continue
    if not runs:
        return f"<p class='none'>no data in {escape(str(wheels))}</p>"
    runs.sort(key=lambda r: r[0])
    when, latest, _ = runs[-1]
    rows = []
    for case in latest["cases"]:
        cells, seconds, previous = [], [], None
        for _, summary, _ in runs:
            row = summary["cases"].get(case)
            if not row:
                cells.append("–")
                continue
            work = "/".join(str(v) for v in row["work_B"].values())
            changed = not row["work_identical"] or (previous is not None and work != previous)
            cells.append(f'<b class="flag">{work}</b>' if changed else work)
            seconds.append(row["median_B_seconds"])
            previous = work
        trend = f"{sparkline(seconds)} {1e3 * seconds[-1]:.1f} ms ({1e3 * min(seconds):.1f}–{1e3 * max(seconds):.1f})"
        rows.append((escape(case), *cells, trend))
    head_row = ["Case"] + [s["source_commit"][:8] for _, s, _ in runs] + ["Wall time per run (latest, range)"]
    return (stamp(time.strftime("%Y-%m-%d %H:%M", time.strptime(when, "%Y%m%dT%H%M%S")),
                  latest["source_commit"], "clean", head)
            + f"<p>Ordinary cases: iterations/evaluations/Hessian calls; F1: reference computations/columns/density solves/equilibrium iterations. "
              f"Counts of the new main wheel (B) per commit, oldest first, "
              f"from <code>{escape(str(wheels))}/key-solves-*/summary.json</code>; a change is bold. "
              f"Wall time is <code>median_B_seconds</code>: one evaluation for F1, {latest['pairs']} timing pair(s) otherwise, on a shared host; "
              f"work counts are the gate.</p>" + table(head_row, rows))


def link(path, text):
    return f'<a href="../{escape(path.relative_to(ROOT).as_posix())}">{text}</a>'


def analyses():
    index = ROOT / "analyses/README.md"
    row = re.compile(r"^\| \[([^]]+)\]\(([^)]+)\) \| (.*) \| (\w+) \|$")
    image = re.compile(r"\]\(([^)\s]+\.(?:png|svg|jpe?g|gif|webp))\)", re.IGNORECASE)
    dirty = git("status", "--porcelain", "--", "analyses")
    rows = []
    for match in filter(None, map(row.match, index.read_text().splitlines() if index.exists() else [])):
        name, rel, question, status = match.groups()
        leaf = (index.parent / rel).resolve()
        readme = leaf / "README.md"
        figure = None
        if readme.exists():
            figure = next((leaf / m for m in image.findall(readme.read_text()) if (leaf / m).exists()), None)
        if figure is None and (leaf / "figures").is_dir():
            figure = next(iter(sorted(p for p in (leaf / "figures").iterdir() if p.is_file())), None)
        commit, date = (git("log", "-1", "--format=%h %cs", "--", str(leaf)).split() + ["", ""])[:2]
        changed = f"analyses/{rel.rstrip('/')}/" in dirty
        rows.append((link(readme if readme.exists() else leaf, escape(name)), escape(question), escape(status),
                     link(figure, escape(figure.name)) if figure else "–",
                     f"{commit} {date}" + (' <b class="flag">dirty</b>' if changed else "")))
    if not rows:
        return "<p class='none'>no data</p>"
    return ("<p>From <code>analyses/README.md</code>; the figure is the first image the README links, "
            "otherwise the first file in <code>figures/</code>. Last change is the last commit touching the leaf.</p>"
            + table(["Leaf", "Question", "Status", "Figure", "Last change"], rows))


def main(args):
    if args[:1] == ["--stamp"]:
        Path(args[1] + ".source").write_text(" ".join(source()) + "\n")
        return
    build = Path(args[0]).resolve() if args else ROOT / "build"
    if args:
        wheels = build / "wheels"
    else:  # Worktrees have no wheels; refresh-main-wheel keeps them in the main checkout.
        main_root = git("worktree", "list", "--porcelain").splitlines()[0].removeprefix("worktree ")
        wheels = Path(main_root) / "build/wheels"
    commit, state = source()
    head = commit, int(git("log", "-1", "--format=%ct"))
    sections = [("Tests", tests(build / "test-results", head)), ("Key solves", key_solves(wheels, head)),
                ("Analyses", analyses())]
    subject = escape(git("log", "-1", "--format=%s"))
    page = f"""<!doctype html><html><head><meta charset="utf-8"><title>ePC-SAFT report</title><style>
body{{font:14px system-ui,sans-serif;margin:1.5em auto;max-width:1200px;padding:0 1em;color:#111}}
table{{border-collapse:collapse;margin:.5em 0}}td,th{{border:1px solid #ccc;padding:2px 6px;text-align:left;vertical-align:top}}
th{{background:#f3f3f3}}.flag{{color:#b91c1c}}.none{{color:#b91c1c;font-style:italic}}.src{{color:#555}}
pre{{margin:0;white-space:pre-wrap;max-width:60em;font-size:12px}}</style></head><body>
<h1>ePC-SAFT report</h1><p>HEAD <code>{commit[:8]}</code> {subject} · engine/ {state}
{'<b class="flag">dirty</b>' if state == 'dirty' else ''} · written {time.strftime('%Y-%m-%d %H:%M')}.
Data from another commit or a dirty tree is <b class="flag">flagged</b>.</p>
{''.join(f'<h2>{title}</h2>{body}' for title, body in sections)}</body></html>"""
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(page)
    print(f"report: {OUT.relative_to(ROOT)} ({len(page.encode()) // 1024} KB)")


if __name__ == "__main__":
    main(sys.argv[1:])

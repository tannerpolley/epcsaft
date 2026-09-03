"""Run a brief's executable claims against a worktree. Exit 1 on the first class of failure."""
import json, subprocess, sys, os
from pathlib import Path
brief, wt = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
d = json.load(brief.open())
env = dict(os.environ, PATH=f"{wt}/.venv/bin:{os.environ['PATH']}", VIRTUAL_ENV=f"{wt}/.venv")
def sh(cmd): return subprocess.run(cmd, shell=True, cwd=wt, env=env, capture_output=True, text=True)
bad = []
for s in d["sources"]:
    ref = s["ref"]; sha_path = ref.split("@", 1)[1]
    r = sh(f"git show {sha_path}")
    if r.returncode: bad.append(f"source missing: {ref}"); continue
    for m in s["must_contain"]:
        if m not in r.stdout: bad.append(f"marker absent: {ref} :: {m[:60]}")
print(f"sources: {len(d['sources'])} checked, {len(bad)} problems")
writes = {w["path"] for w in d["writes"]}
for g in d["gates"]:
    if not g.get("mutation"): bad.append(f"gate without mutation: {g['id']}")
    mf = g["mutation"]["file"]
    covered = mf in writes or any(w.endswith("/") and mf.startswith(w) for w in writes)
    if not (wt / mf).exists() and not covered: bad.append(f"mutation file neither exists nor in writes: {g['id']} {mf}")
    # A sed that matches nothing mutates nothing, and the gate is never challenged --
    # indistinguishable from a gate that does not bite. Only checkable on files that
    # already exist; a write the stage has yet to create is checked at acceptance.
    if (wt / mf).exists():
        import subprocess, tempfile
        before = (wt / mf).read_bytes()
        with tempfile.NamedTemporaryFile(delete=False) as tmp: tmp.write(before)
        subprocess.run(["sed", "-i", g["mutation"]["sed"], tmp.name], capture_output=True)
        if Path(tmp.name).read_bytes() == before:
            bad.append(f"mutation sed matches nothing: {g['id']} {mf}")
        Path(tmp.name).unlink()
print(f"gates: {len(d['gates'])}, each with a mutation: {all(g.get('mutation') for g in d['gates'])}")
pre_fail = []
for p in d["preconditions"]:
    r = sh(p)
    if r.returncode: pre_fail.append((p, (r.stderr or r.stdout).strip()[-300:]))
print(f"preconditions: {len(d['preconditions'])} run, {len(pre_fail)} failed")
# the brief's prose must not name a stage, branch or brief other than its own
import re as _re
md = brief.with_suffix("").with_suffix(".brief.md")
if md.exists():
    stem = brief.name.split(".")[0]
    text = md.read_text()
    base_stem = stem.replace("-repair", "")   # a repair brief may name the stage it repairs
    wrong = {m for m in _re.findall(r"stage-\d+[a-z-]*", text) if not m.startswith(base_stem)}
    wrong |= {m for m in _re.findall(r"branch `(stage-[0-9a-z-]+)`", text) if m != d.get("branch")}
    if wrong: bad.append("brief prose names another stage: " + ", ".join(sorted(wrong)))
    print(f"self-reference: {'clean' if not wrong else 'MISMATCH'}")
for p, e in pre_fail: print("  FAIL:", p, "\n       ", e.replace("\n", "\n        "))
for b in bad: print("  ", b)
sys.exit(1 if bad or pre_fail else 0)

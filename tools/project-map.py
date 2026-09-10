#!/usr/bin/env python3
"""Resolve every project file to one owner and role from ARCHITECTURE.yaml.

File lists come from the filesystem (tracked plus untracked, not ignored). The
longest matching prefix rule wins; unmapped and ambiguous paths are reported,
never dropped. Generated, environment and independent-repository roots are
summarized by presence only.
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SUMMARIZED_ROLES = {"generated", "external"}


def declaration(root: Path = ROOT) -> dict[str, Any]:
    text = (root / "ARCHITECTURE.yaml").read_text(encoding="utf-8")
    match = re.search(r"^project_map:", text, re.MULTILINE)
    if match is None:
        raise ValueError("ARCHITECTURE.yaml has no project_map mapping")
    value, _ = json.JSONDecoder().raw_decode(text[match.end():].lstrip())
    components, roles, rules = value["components"], set(value["roles"]), value["rules"]
    for prefix, owner, role in rules:
        if owner not in components or role not in roles:
            raise ValueError(f"rule {prefix!r} names unknown owner or role")
    for name, item in components.items():
        if item.get("parent") not in {None, *components}:
            raise ValueError(f"component {name} has unknown parent")
    return value


def project_files(root: Path = ROOT) -> list[str]:
    result = subprocess.run(
        ["git", "-C", str(root), "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
        capture_output=True, check=True,
    )
    return sorted({item for item in result.stdout.decode().split("\0") if item and (root / item).exists()})


def resolve(path: str, rules: list[list[str]]) -> tuple[str | None, str | None, bool]:
    """Return (owner, role, ambiguous) by the longest matching prefix."""

    best: list[list[str]] = []
    for rule in rules:
        pattern = rule[0]
        matched = fnmatch.fnmatchcase(path, pattern) if any(c in pattern for c in "*?[") else (
            path == pattern or (pattern.endswith("/") and path.startswith(pattern))
        )
        specificity = len(pattern.replace("*", "").replace("?", ""))
        if matched:
            best_specificity = len(best[0][0].replace("*", "").replace("?", "")) if best else -1
            if specificity > best_specificity:
                best = [rule]
            elif pattern == best[0][0] and rule[1:] != best[0][1:]:
                best.append(rule)
    if not best:
        return None, None, False
    return best[0][1], best[0][2], len(best) > 1


def project_map(root: Path = ROOT) -> dict[str, Any]:
    declared = declaration(root)
    rules = declared["rules"]
    nodes: dict[str, dict[str, Any]] = {}
    for name, item in declared["components"].items():
        nodes[name] = {
            "id": name, "name": item["name"], "purpose": item.get("purpose", ""),
            "parent": item.get("parent"), "external": item.get("external", []),
            "paths": [rule[0] for rule in rules if rule[1] == name],
            "files": {}, "summarized": [],
        }
    unmapped: list[str] = []
    ambiguous: list[str] = []
    for path in project_files(root):
        owner, role, conflict = resolve(path, rules)
        if owner is None:
            unmapped.append(path)
            continue
        if conflict:
            ambiguous.append(path)
        nodes[owner]["files"].setdefault(role, []).append(path)
    for prefix, owner, role in rules:
        if role in SUMMARIZED_ROLES and (root / prefix).exists():
            nodes[owner]["summarized"].append({"path": prefix, "role": role})
    return {
        "nodes": list(nodes.values()),
        "unmapped": unmapped,
        "ambiguous": ambiguous,
        "complete": not unmapped and not ambiguous,
    }


def self_test() -> None:
    rules = [["engine/", "engine", "specification"], ["engine/tests/", "engine", "test"],
             ["engine/tests/", "tooling", "test"], ["README.md", "documentation", "specification"]]
    assert resolve("engine/tests/x.py", rules) == ("engine", "test", True)
    assert resolve("engine/native/a.cpp", rules) == ("engine", "specification", False)
    assert resolve("README.md", rules) == ("documentation", "specification", False)
    assert resolve("README.md.bak", rules) == (None, None, False)
    assert resolve(
        "analyses/case/scripts/run.py",
        [["analyses/", "analyses", "evidence"], ["analyses/**/scripts/**", "analyses", "tooling"]],
    ) == ("analyses", "tooling", False)
    live = project_map()
    assert any(node["id"] == "engine.core" and node["files"].get("runtime") for node in live["nodes"])
    print("project map self-test passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    result = project_map()
    if args.json:
        json.dump(result, sys.stdout, indent=1)
        return 0 if result["complete"] else 1
    for node in result["nodes"]:
        counted = sum(len(files) for files in node["files"].values())
        print(f"{node['id']:22} {counted:5} files  " + ", ".join(f"{role}={len(files)}" for role, files in sorted(node["files"].items())))
    for label, paths in (("unmapped", result["unmapped"]), ("ambiguous", result["ambiguous"])):
        for path in paths:
            print(f"{label}: {path}", file=sys.stderr)
    if not result["complete"]:
        return 1
    print("project map complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

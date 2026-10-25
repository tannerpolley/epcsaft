#!/usr/bin/env python3
"""Reconcile derived GitHub issue status labels from native dependencies, clear work:claimed and
work:review from closed issues, strip area prefixes from titles, and give every open issue a milestone."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from typing import Any

PROGRAM = "roadmap:program"
DERIVED = {
    "status:planned",
    "status:blocked",
    "status:ready",
    "status:active",
    "status:review",
    "status:needs-attention",
}
LABELS = {
    PROGRAM: ("5319e7", "Program parent; not an executable work item"),
    "work:authorized": ("0e8a16", "Implementation is explicitly authorized"),
    "work:paused": ("b60205", "Deliberately paused by an owner"),
    "work:claimed": ("fbca04", "An assigned worker currently owns execution"),
    "work:review": ("d93f0b", "Implementation is awaiting review or acceptance"),
    "status:planned": ("cfd3d7", "Derived: prerequisites satisfied but work is not authorized"),
    "status:blocked": ("b60205", "Derived: paused or waiting on an open prerequisite"),
    "status:ready": ("0e8a16", "Derived: authorized with all prerequisites satisfied"),
    "status:active": ("1d76db", "Derived: claimed by an assigned worker"),
    "status:review": ("d93f0b", "Derived: awaiting review or acceptance"),
    "status:needs-attention": ("e99695", "Derived: dependency or ownership state needs a decision"),
}


def gh(*args: str) -> str:
    return subprocess.run(
        ("gh", *args), check=True, text=True, capture_output=True
    ).stdout.strip()


def pages(endpoint: str) -> list[dict[str, Any]]:
    result = json.loads(gh("api", "--paginate", "--slurp", endpoint))
    return [item for page in result for item in page]


def desired_status(issue: dict[str, Any], blockers: list[dict[str, Any]]) -> str | None:
    if issue["state"] != "open":
        return None
    labels = {label["name"] for label in issue["labels"]}
    if PROGRAM in labels:
        return None
    if any(blocker.get("state_reason") == "not_planned" for blocker in blockers):
        return "status:needs-attention"
    open_blockers = any(blocker["state"] == "open" for blocker in blockers)
    if "work:paused" in labels:
        return "status:blocked"
    if open_blockers:
        if labels & {"work:claimed", "work:review"}:
            return "status:needs-attention"
        return "status:blocked"
    if "work:review" in labels:
        return "status:review"
    if "work:claimed" in labels:
        return "status:active" if issue["assignees"] else "status:needs-attention"
    if "work:authorized" in labels:
        return "status:ready"
    return "status:planned"


def ensure_labels(repo: str) -> None:
    existing = {
        item["name"]: (item["color"].lower(), item["description"] or "")
        for item in json.loads(
            gh(
                "label",
                "list",
                "-R",
                repo,
                "--limit",
                "100",
                "--json",
                "name,color,description",
            )
        )
    }
    for name, (color, description) in LABELS.items():
        if existing.get(name) == (color, description):
            continue
        gh(
            "label",
            "create",
            name,
            "-R",
            repo,
            "--color",
            color,
            "--description",
            description,
            "--force",
        )


BODY_LIMIT = 8000  # characters: current state in about one screen

# Titles carry no "[Area] " prefix; the milestone holds the area (owner, 2026-10-12, #456).
# An open issue created without a milestone joins the #220 platform milestone for triage.
PLATFORM = "12 - General Platform Completion and Release"
PREFIX = re.compile(r"^(\[[^\]]*\]\s*)+")


def reconcile(repo: str, *, apply: bool, only: set[int]) -> int:
    if apply:
        ensure_labels(repo)
    all_issues = [
        issue
        for issue in pages(f"repos/{repo}/issues?state=all&per_page=100")
        if "pull_request" not in issue
    ]
    issue_by_number = {issue["number"]: issue for issue in all_issues}
    issues = [issue for issue in all_issues if not only or issue["number"] in only]
    relations = json.loads(
        gh(
            "issue",
            "list",
            "-R",
            repo,
            "--state",
            "all",
            "--limit",
            "500",
            "--json",
            "number,blockedBy",
        )
    )
    blocked_by = {item["number"]: item["blockedBy"] for item in relations}
    changes = 0
    for issue in issues:
        if issue["state"] == "open" and len(issue.get("body") or "") > BODY_LIMIT:
            print(f"#{issue['number']}: body is {len(issue['body'])} characters; compact it to current state")
        title = PREFIX.sub("", issue["title"])
        if title != issue["title"]:
            changes += 1
            print(f"#{issue['number']}: title -> {title}")
            if apply:
                gh("issue", "edit", str(issue["number"]), "-R", repo, "--title", title)
        if issue["state"] == "open" and not issue["milestone"]:
            changes += 1
            print(f"#{issue['number']}: no milestone -> {PLATFORM}")
            if apply:
                gh("issue", "edit", str(issue["number"]), "-R", repo, "--milestone", PLATFORM)
        stale = {label["name"] for label in issue["labels"]} & {"work:claimed", "work:review"}
        if issue["state"] == "closed" and stale:
            changes += 1
            print(f"#{issue['number']}: closed; remove {sorted(stale)}")
            if apply:
                gh("issue", "edit", str(issue["number"]), "-R", repo, "--remove-label", ",".join(sorted(stale)))
        relationship = blocked_by.get(issue["number"], {"nodes": [], "totalCount": 0})
        if relationship["totalCount"] != len(relationship["nodes"]):
            wanted = "status:needs-attention" if issue["state"] == "open" else None
        else:
            blockers = [
                issue_by_number[node["number"]]
                for node in relationship["nodes"]
                if node["number"] in issue_by_number
            ]
            wanted = (
                desired_status(issue, blockers)
                if len(blockers) == len(relationship["nodes"])
                else "status:needs-attention"
            )
        current = {label["name"] for label in issue["labels"]} & DERIVED
        target = {wanted} if wanted else set()
        if current == target:
            continue
        changes += 1
        print(f"#{issue['number']}: {sorted(current)} -> {sorted(target)}")
        if not apply:
            continue
        command = ["issue", "edit", str(issue["number"]), "-R", repo]
        if current - target:
            command += ["--remove-label", ",".join(sorted(current - target))]
        if target - current:
            command += ["--add-label", ",".join(sorted(target - current))]
        gh(*command)
    print(f"{'applied' if apply else 'planned'} {changes} change(s)")
    return changes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, help="OWNER/REPO")
    parser.add_argument("--apply", action="store_true", help="Apply label and milestone changes")
    parser.add_argument("--issue", type=int, action="append", default=[])
    args = parser.parse_args()
    reconcile(args.repo, apply=args.apply, only=set(args.issue))


if __name__ == "__main__":
    main()

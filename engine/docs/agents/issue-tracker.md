# Issue tracker: GitHub

Issues and PRDs for this repo live as GitHub issues in
`ePC-SAFT/ePC-SAFT`. Use the `gh` CLI for all operations.

## Conventions

- **Create an issue**: `gh issue create --title "..." --body "..."`. Use a
  body file for multi-line bodies.
- **Read an issue**: `gh issue view <number> --comments`, also fetching labels
  when needed.
- **List issues**:
  `gh issue list --state open --json number,title,body,labels,comments` with
  appropriate `--label` and `--state` filters.
- **Comment on an issue**: `gh issue comment <number> --body "..."`
- **Apply or remove labels**: `gh issue edit <number> --add-label "..."` or
  `--remove-label "..."`
- **Close**: `gh issue close <number> --comment "..."`

Infer the repository from `git remote -v`; `gh` does this automatically inside
the clone.

## Pull requests as a triage surface

**PRs as a request surface: no.**

GitHub shares one number space across issues and PRs, so resolve an ambiguous
`#42` with `gh pr view 42` and fall back to `gh issue view 42`.

## When a skill says "publish to the issue tracker"

Create a GitHub issue.

## When a skill says "fetch the relevant ticket"

Run `gh issue view <number> --comments`.

## Wayfinding operations

Used by `/wayfinder`. The map is one issue with child issues as tickets.

- **Map**: an issue labelled `wayfinder:map`, holding Notes,
  Decisions-so-far, and Fog.
- **Child ticket**: a GitHub sub-issue, or a task-list child with
  `Part of #<map>` when sub-issues are unavailable. Use `wayfinder:<type>`
  labels.
- **Blocking**: use GitHub native issue dependencies. Fall back to a
  `Blocked by: #<n>` line only when dependencies are unavailable.
- **Frontier query**: select the first unassigned open child without open
  blockers, in map order.
- **Claim**: `gh issue edit <n> --add-assignee @me`.
- **Resolve**: comment with the answer, close the child, and append its context
  pointer to the map's Decisions-so-far.

<!-- project-truss:setup:start -->
## Project Truss tracker contract

Project Truss uses GitHub Issues in `tannerpolley/ePC-SAFT` through the authenticated `gh` CLI.

## Project Truss execution issues

- Root issues use Problem Statement, Solution, numbered User Stories, Implementation Decisions, Testing Decisions, Out of Scope, and Further Notes.
- Executable leaves use Parent, What to build, Acceptance criteria, and Blocked by.
- Native sub-issue and blocked-by relationships are authoritative.

## Wayfinder decision issues

Wayfinder maps and decision tickets clarify an uncertain destination before execution shaping. They are not Project Truss execution contracts. Their labels, assignees, comments, and closure never derive Truss readiness or completion. Link a completed map from a new Truss root; do not automatically convert or reuse its issues.

### Wayfinding operations

- Create one map issue with the `wayfinder:map` label and `Destination`, `Notes`, `Decisions so far`, `Not yet specified`, and `Out of scope` sections.
- Create decision tickets as native child issues with a `## Question` body and one `wayfinder:research`, `wayfinder:prototype`, `wayfinder:grilling`, or `wayfinder:task` label.
- Express blocking with native GitHub issue dependencies and derive the frontier from open, unblocked, unassigned children.
- Claim by assignment. Resolve with one answer comment, close the ticket, and append only a linked gist to the map's `Decisions so far`.

## Authority

GitHub, Git, reviews, integration, and current worktrees are lifecycle truth. Provider checks are optional context. Labels, milestones, Projects, issue text, and comments cannot grant authority.
<!-- project-truss:setup:end -->

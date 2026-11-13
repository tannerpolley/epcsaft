---
name: manage-github-roadmap
description: Manage this repository's issue-backed implementation roadmap, native dependencies, authorization, claims and handoffs. Use before starting, changing, or completing roadmap work.
---

# Manage the GitHub roadmap

[ADR 0010](../../../docs/scientific/decisions/0010-github-issue-roadmap-authority.md)
gives GitHub the plans, scope, acceptance, prerequisites, authorization and
execution status; never copy editable roadmap status into repository documents.

## Load the accepted plan

- Fetch the complete current issue body through `gh issue view --json body`
  before implementation. If terminal output truncates, write it to a task-owned
  temporary file, read every section, and delete the file after use. Never
  reconstruct the assignment from a search excerpt, cached copy, summary, or
  prior handoff.
- Treat the issue body as the current accepted plan. Read its precise linked
  definitions, relevant comments, and native relationships. Comments hold
  discussion, review findings, amendments, and execution reports; incorporate
  accepted changes into the body and remove superseded instructions.
- Do not create or commit a replacement implementation plan. Record a compact
  plan revision and implementation baseline in the issue when work starts;
  reread the body after interruption and before handoff.
- Before implementation, confirm the intended changes, governing decisions,
  acceptance checks, and any concrete missing prerequisite. Resolve routine
  coding details locally. Return only consequential unresolved science,
  architecture changes, or demonstrated contradictions to the accountable
  planner, stopping only the affected work.

## Before work

1. Read the leaf issue body, its parent, native `Blocked by` relationships, and
   linked repository evidence.
2. Proposed plan creation and maintenance do not authorize implementation. Do
   not implement unless the issue has `work:authorized`, every prerequisite is
   satisfied, and it is not paused or claimed. Readiness is not merge authority.
3. Claim by assigning the accountable worker and adding `work:claimed`.

## During work

- Keep the current design, implementation plan, scope, exclusions, acceptance,
  reuse decision, and remaining work in the issue body, as **current state** in
  about one screen (under 8,000 characters). Superseded rounds, review
  transcripts and run reports go to comments or retained evidence.
- Keep one accountable body editor at a time. Workers may report progress in
  comments without racing to rewrite the accepted plan.
- Create every issue with `--milestone` and a title without an `[Area]` prefix; the
  milestone holds the area.
- Use native parent/sub-issues only for independently executable deliverables.
- Record prerequisites with native `Blocked by`; preserve the relationship after
  completion.
- Use `work:paused` for a deliberate pause and `work:review` for review. A review
  suggestion becomes a blocker only after it changes accepted scope or acceptance.
- Keep scientific calculations and evidence in their repository owners and link
  them from the issue.

## Handoff and completion

Update the issue with the delivered behavior, evidence and limits, exact revision
or PR, remaining work, and any changed dependencies. Remove `work:claimed` only
after verifying the worker stopped. Close a leaf only
when its acceptance is supported. Reopen it when delivery or qualification is
invalidated.

Every PR closes its issue (owner, 2026-10-14): its body says `Closes #N`, and it
merges only when that issue is done. A slice, follow-up or leftover correction gets
its own new issue before its PR merges; never merge a PR that leaves its issue open.
In the same closeout, write the delivered state into the issue
(replace, do not append dated blocks), remove the PR's `work:review` and `work:claimed`
labels. A comment alone is not a handoff.
On closing #N, run `git grep -nw '#N' -- docs engine/docs 'analyses/**/README.md'` and update every statement that #N is open or pending.

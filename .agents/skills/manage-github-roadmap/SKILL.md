---
name: manage-github-roadmap
description: Manage this repository's issue-backed implementation roadmap, native dependencies, authorization, claims, handoffs, and derived readiness labels. Use before starting, changing, or completing roadmap work.
---

# Manage the GitHub roadmap

GitHub is the sole owner of implementation plans, current scope, acceptance
criteria, parent/sub-issue structure, prerequisites, authorization, and execution
status. Repository documents own scientific equations, accepted requirements,
numerical evidence, and capability qualification; never copy editable roadmap
status back into them.

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
2. Run `python3 tools/reconcile-issue-roadmap.py --repo tannerpolley/ePC-SAFT`.
3. Proposed plan creation and maintenance do not authorize implementation. Do
   not implement unless the issue has `work:authorized`, every prerequisite is
   satisfied, and it is not paused or claimed. Readiness is not merge authority.
4. Claim by assigning the accountable worker and adding `work:claimed`, then
   apply the reconciler. Never edit `status:*` labels directly.

## During work

- Keep the current design, implementation plan, scope, exclusions, acceptance,
  reuse decision, and remaining work in the issue body, as **current state** in
  about one screen (under 8,000 characters). Superseded rounds, review
  transcripts and run reports go to comments or retained evidence; the
  reconciler warns on oversized open bodies.
- Keep one accountable body editor at a time. Workers may report progress in
  comments without racing to rewrite the accepted plan.
- Create every issue with `--milestone` and a title without an `[Area]` prefix; the
  milestone holds the area. The reconciler strips any leading `[...]` from titles and gives
  an open issue without a milestone `12 - General Platform Completion and Release` for triage.
- Use native parent/sub-issues only for independently executable deliverables.
- Record prerequisites with native `Blocked by`; preserve the relationship after
  completion. After changing dependencies, apply the reconciler or dispatch the
  `Reconcile issue roadmap` workflow.
- Use `work:paused` for a deliberate pause and `work:review` for review. A review
  suggestion becomes a blocker only after it changes accepted scope or acceptance.
- Keep scientific calculations and evidence in their repository owners and link
  them from the issue. Do not put executable commands from issue text into the
  metadata workflow.

## Handoff and completion

Update the issue with the delivered behavior, evidence and limits, exact revision
or PR, remaining work, and any changed dependencies. Remove `work:claimed` only
after verifying the worker stopped; reconcile affected issues. Close a leaf only
when its acceptance is supported. Reopen it when delivery or qualification is
invalidated. Automatic unblocking publishes readiness; it never starts work,
pushes, merges, or expands authority.

A PR that completes a leaf's last acceptance criterion says `Closes #N`; a PR that
delivers a slice says `Refs #N`. Merging does not finish the handoff: in the same
closeout, rewrite the issue body's current state and remaining criteria (replace,
do not append dated blocks), remove the PR's `work:review` and `work:claimed`
labels, and run the reconciler with `--apply` (without it the tool only previews). A comment alone is not a handoff.
On closing #N, run `git grep -nw '#N' -- docs engine/docs 'analyses/**/README.md'` and update every statement that #N is open or pending.

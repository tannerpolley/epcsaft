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

## Before work

1. Read the leaf issue body, its parent, native `Blocked by` relationships, and
   linked repository evidence.
2. Run `python3 tools/reconcile-issue-roadmap.py --repo tannerpolley/ePC-SAFT-greenfield`.
3. Do not start unless the issue has `work:authorized`, every prerequisite is
   satisfied, and it is not paused or claimed. Readiness is not merge authority.
4. Claim by assigning the accountable worker and adding `work:claimed`, then
   apply the reconciler. Never edit `status:*` labels directly.

## During work

- Keep the current design, implementation plan, scope, exclusions, acceptance,
  reuse decision, and remaining work in the issue body.
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

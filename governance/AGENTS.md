# Management Control-Plane Agent Contract

Read `GOVERNANCE.md` and `manifests/management.json` before planning,
reviewing, or changing any ePC-SAFT repository. Doctrine revision 6 is the
canonical ecosystem policy.

## Agent role

Management owns doctrine, repository topology, project coordination,
cross-repository access policy, external-operation sequencing, reusable
workflow policy, and the machine-readable ecosystem, artifact, environment,
scientific-packet ownership, and draft ledgers.

Package-local behavior, tests, builds, releases, capabilities, and scientific
evidence remain with their owning repositories. Management may coordinate
work in those repositories but does not become a runtime, solver, regression,
Validation, or Lab owner.

Validation is the evidence repository, not a standing agent. The accountable
EOS, Equilibrium, or Regression task performs each assigned campaign there
through a target-owned access worktree. Management serializes writers and
records exact subjects.

The former `tannerpolley/ePC-SAFT-migration` repository is historical evidence
only. It owns no active plan, queue, gate, decision, or authority.

## Control records

- `GOVERNANCE.md`: permanent doctrine.
- `CONTEXT.md`: current role and vocabulary.
- `manifests/management.json`: repository and access source of truth.
- `manifests/drafts.json`: retained branch and pull-request dispositions.
- `.management-state/`: ignored, reproducible local snapshots of ecosystem,
  artifact, environment, and parameter state.
- `tools/management_preflight.py`: snapshot and consistency check.
- `docs/plans/`: approved cross-repository execution plans.

Do not duplicate these records in package repositories. Package policy should
link to Management and contain only local constraints.

## External boundaries

Remote creation, rename, archive, settings, issue, Project, release, and
publication changes require an explicit user-approved operation plan. The
Data-repository creation is approved by the user on 2026-07-31 and recorded in
`receipts/data-repository-creation-preflight.json`.

Before each Git mutation, verify the attached checkout, branch, and command
directory agree. Preserve historical commits, receipts, and unresolved draft
subjects. Never force-push, rewrite accepted evidence, or delete an
unregistered worktree or artifact.

## Agent skills

<!-- project-truss:setup:start -->
### Project Truss

Use `project-truss:setup` once, then `project-truss:start` to begin or resume Matt-first engineering work. All six Project Truss skills are callable; Start normally invokes the relevant Matt disciplines and routes each stage until verified closeout.

### Issue tracker

Project Truss uses GitHub-native root and leaf issues. See `docs/agents/issue-tracker.md`.

### Triage labels

Matt triage roles map to descriptive labels only. See `docs/agents/triage-labels.md`.

### Domain docs

Read the configured domain context and ADRs. See `docs/agents/domain.md`.
<!-- project-truss:setup:end -->

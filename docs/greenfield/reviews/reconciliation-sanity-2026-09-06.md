# Documentation and pending-work sanity review

Scope: the uncommitted documentation reconciliation on attached branch
`greenfield`, baseline `ff614901be22727095947e71fb4626ef2dae7a7c`.
This review changes documentation and local readiness tooling only. It does
not qualify numerical Engine behavior or execute the next evidence campaign.

## Corrections from the deeper review

The earlier reconciliation left several inherited contradictions. They were
corrected in the canonical definitions and regenerated YAML index:

- Separate packing, DH finite-size and Born diameter records; equality requires
  a source-specific relation. The current native owners already read separate
  packing and DH records.
- Identify mechanical phase energy separately from reference-complete energy.
  Define dimensional chemical/charge multipliers, dimensionless dual multipliers
  and optimizer callback multipliers independently, including the RT scaling.
- Mark the historical trial-pressure source-activity transfer as unresolved
  provenance. It cannot override the fixed-P0 standard-state conversion. Its
  derivation remains S4 work; changing every pressure symbol would not resolve it.
- Make Shomate enthalpy's kJ/mol-to-J/mol factor explicit and distinguish the
  NIST H coefficient from an enthalpy observable. Verified against the
  [NIST WebBook unit definitions](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=203).
- Permit zero component amounts in the composition definition at positive total
  amount; distinguish interior logarithms, central limits and insertion actions.
- Identify the DH polynomial as a finite Taylor approximation, and terminal
  reference molality as a limit estimator. Neither is an exact finite-point
  identity. Point to the retained dispersion coefficient arrays/source table.
- Remove remaining live stage-status/evidence-class/formulation claims from
  `engine/ARCHITECTURE.yaml`; it points to the permanent doctrine and current
  transfer plan.

These are documentation/derivation findings. They are not newly measured code
failures. The dated transfer audit remains unchanged, including its execution
limitations and original evidence tables.

## Related branch inventory and preservation

All eight stage branch tips are ancestors of the inspected greenfield baseline.
Each associated worktree has zero staged and zero unstaged tracked changes.
There is no stage implementation to merge or cherry-pick.

| Branch | Worktree basename | Tip | Untracked files |
|---|---|---|---:|
| stage-02 | gf-stage-02 | cf61e75017a0dc1b52eb7ee3ab074ae3771ce704 | 9 |
| stage-03-repair | gf-stage-03-repair | a40839bdab2e95080a549328373d793287e3f41f | 0 |
| stage-04-repair | gf-stage-04-repair | 7746996f14a11003b39ed5c883842072e3e1d57b | 0 |
| stage-05 | gf-stage-05 | d3f339295e373b386ef556621047a5a82a78cebb | 13 |
| stage-05-harness | gf-stage-05-harness | 4dcf8b39937093995da4f9943969172a488baa59 | 3 |
| stage-05-repair | gf-stage-05-repair | 542df61364c8e4f951f45855f1abbe0e257883e6 | 43 |
| stage-06 | gf-stage-06 | a518039d4fd430444200c31a7bf4c9f2decc6366 | 7 |
| stage-04 | gf-stage-3 | 6f4d0b32ca1e96841175b0e80126e614059986f1 | 7 |

The [file inventory](greenfield-pending-2026-09-06.csv) binds all 82 untracked
files (113,197,853 bytes) to original path, branch, HEAD and SHA-256. Twenty
worker result JSON records were consolidated byte-for-byte under
[stage-history](stage-history/), with retained destinations in the inventory.
They report historical attempts and are not current acceptance results.
Original files remain in their worktrees pending the exact cleanup decision.
The other files are logs, resume prompts, result schemas, copied briefs, dead
PID files and local Serena configuration. A hash inventory is not a backup of
those bytes; deleting them requires an explicit disposition.

The harness result explains the 26-registration/24-claim correction and the
missing-regression verifier stop. Stage-06 stops name missing thermochemistry
and parameter-derivative owners; neither wrote an implementation. The repair
stops expose contradictory ceilings and early-stop rules that preserved broken
intermediate work. These requirements survive in TRANSFER S0/S2/S4/S6 and the
new code-budget rule. The stale root stage-06 brief is not revived.

The two pre-existing current harness brief edits remain byte-identical:

| File | SHA-256 |
|---|---|
| stage-05-harness.brief.md | 0bbfda90840a8c20d24193fe6b805a42996ee42be7161c5b72b14cff7fd13af6 |
| stage-05-harness.brief.json | e2677aefd4835a1d9b013ae5ffd8a0a28b997fa00a02e796077f98db71924ffe |

Unique historical refs are preserved: `perf/instrumentation` at
`e4f2535ea207b1785d121aefd2238b8226fd086e` (3 unique commits),
`perf/association-newton` at `a248de6b08cb960645a468be51e4a5c973a20895`
(2), `perf/lagrangian-hessian` at `fb02e6dfa081ab0caaac26e95e91b6800ef0a499`
(3), and `audit/gut-prototype` at `2f6677f6f5ed6d1114385da5a7e419be4debc371`
(3). Their obsolete runtime paths must not be merged wholesale. Association
and continuation ideas remain useful evidence for the transfer plan.
`docs/equilibrium-formulation` is fully integrated. Missing old scratch-worktree
registrations are cleanup candidates, not lost commits.

Both recorded stage PID numbers are dead, but live language-server processes
still reference gf-stage-05-harness. That worktree and those processes remain
untouched. Unrelated issue-54-mdea's dirty design document and all downstream
application checkouts remain untouched. No historical worker was resumed.

## Code budget decision and enforcement

The requested floor target, soft ceiling and hard ceiling have one numeric
owner in `engine/ARCHITECTURE.yaml`, interpreted by
[code-budgets.md](../../scientific/code-budgets.md). Current runtime is 8,807
physical source lines; completed-fluid target/soft/hard are
12,300 / 14,850 / 17,400, with independently enforced layer limits.
These are engineering estimates, not theoretical minima or performance claims.

Existing `tools/diffstat.py --budgets` counts current tracked/untracked source,
including comments, blank lines and stubs. Hard excess fails; soft excess warns
and requires an ownership/removal explanation before acceptance. Tests, docs
and dependencies are outside the runtime limit, not places to hide runtime
logic. Readable formatting and required science cannot be sacrificed for counts.
An intermediate excess can be repaired; accepting a hard excess needs an
explicit owner-approved amendment. A small runnable self-check falsifies the
soft/hard boundaries, untracked-file accounting and unclassified-layer bypass.

## Verification and limits

The equation synchronization/self-check and current-tree budget/self-check run
through `tools/check-affected.sh HEAD`. Local links, retained JSON hashes,
unchanged harness hashes and absence of native/Python runtime diffs are checked
separately. Independent read-only reviewers assessed scientific corrections and
budget/tooling behavior. Scientific review: prior findings resolved after correcting the last terminal-change
label. Budget/preservation review verified all hashes/counts and required the
budget gate to run even before the clean-tree early exit; corrected. A broader
style scan found pre-existing formatting findings in the old diff-report mode;
the focused undefined-name/syntax checks pass. Final outcomes are reported in
the task handoff.

No wheel build, numerical campaign, source acquisition, GitHub write, branch
rewrite or worktree/process removal is implied by this record. A finite sanity
review cannot prove absence of every mathematical error or semantic duplication.
The next discussion should select independent numerical falsifiers and accepted
criteria for the observable/domain being implemented; this task does not
silently turn proposed audit tolerances into acceptance thresholds.

## Authorized retention and verified archive

The owner selected native UI options: archive and clear originals, commit
locally, and remove the seven inactive stage worktrees while retaining harness.
No push or branch deletion was authorized or performed.

All 82 inventory files plus the inventory itself were archived and read back
before any original was removed. Every file's SHA-256 and byte length matched.
The private archive is
`~/Workspaces/Engineering/epcsaft-audit-2026-09-02/archives/greenfield-stage-records-2026-09-06.tar.gz`;
compressed size 20,430,057 bytes, archive SHA-256
`c0a3707a22fd83f36bd1e2d8aa6191f71184da7b6f4f918176151564db67377a`.
Members use `<worktree-basename>/<original-relative-path>`; `inventory.csv`
inside the archive retains the original full paths and hashes. Recover files
into a separate directory, never blindly over a working checkout.

The inventory above describes the pre-cleanup snapshot. The verified archive
is the complete byte-preservation owner; the 20 committed worker reports are
the compact review surface. Ignored build/dependency/virtual-environment caches
are disposable worktree contents and are not included in this historical-file
archive. Cleanup completion is recorded below after the authorized operations.

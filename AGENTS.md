# ePC-SAFT Monorepo Agent Contract

Governance doctrine revision 24 applies.

This private repository is the single Git and Codex home for the ePC-SAFT
project. Read `governance/GOVERNANCE.md`, `governance/manifests/management.json`,
and the nearest subtree `AGENTS.md` before material work.

One scientific slice has one accountable agent from input to evidence. That
agent may change the exact files needed across `data/`, `engine/`, and
`validation/` in one branch. Keep one writer per file and do not create nested
Git repositories, submodules, sibling-repository dependencies, or separate
package installations.

## Official downstream applications

The following sibling Git repositories are official downstream applications:

- `~/Workspaces/Engineering/MEA-Thermodynamics`
- `~/Workspaces/Engineering/MEA-Absorption-Column`
- `~/Workspaces/Engineering/Lithium_Extraction`

When a requested Engine slice requires application implementation or evidence,
the same accountable agent may inspect, edit, test, and run Git commands in the
applicable named repository without moving the Codex task attachment. This
authorization applies only when the task filesystem exposes the target path;
repository instructions cannot override an app sandbox.

Before every downstream Git mutation, verify the exact repository root,
branch, and working-tree state, read its nearest `AGENTS.md`, and keep each
repository's changes separately reviewable. Do not create nested repositories,
submodules, or runtime sibling-source dependencies. Every application must
consume an explicitly identified, non-editable Engine wheel and record its
SHA-256 rather than import the Engine source checkout. Preserve one writer per
file across the participating repositories.

Directory ownership remains strict:

- `governance/`: doctrine, topology, manifests, ledgers, and operation receipts.
- `engine/`: the sole `epcsaft` runtime, including EOS, equilibrium, regression,
  compact numerical evidence, and one native build.
- `data/`: canonical reusable scientific packets and packet readiness checks.
- `validation/`: installed-wheel literature, experimental, and cross-module
  scientific campaigns.

The root `.codex/environments/environment.toml` is the only Codex environment.
It creates one root `.venv`, builds `engine/` into one non-editable wheel, and
rejects split or legacy ePC-SAFT distributions. This is a private,
in-development scientific engine used by downstream applications, so the
wheel is built and installed before any test process starts. Use
`tools/check-affected.sh` for normal work and `tools/check-all.sh` for the
parallel fast installed-wheel gate, with a warm-wheel target below ten seconds.
Slow or expensive Engine evidence, Validation campaigns, full native
clang-tidy, and sanitizers run only by an explicit owner command when their
specific risk or scientific claim requires them; a changed tool or document
must not silently invoke them. `tools/run-slice.sh khudaida|figiel` owns the two
canonical vertical slices.

Use Numerical Verification-Driven Development from Governance section 7.
Persistent tests must falsify numerical claims; expensive installed-artifact
replays belong in Validation; builds, lint, typing, schemas, hashes, manifests,
environments, and path audits are readiness checks. Do not add contract,
compatibility, legacy, snapshot, prose, mock, or coverage-only tests.

Historical plans, receipts, results, evidence, source locators, old repository
names, and old package names are provenance and must not be normalized as live
instructions. This project has no CI/CD. Run the applicable root checks before
handoff. `codex-cleanup.sh` is an agent closeout action: agents run it, never ask
the user to run it, and remove only task-created ignored caches or exact
user-authorized candidates.

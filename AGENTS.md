# ePC-SAFT Monorepo Agent Contract

Repository Profile: scientific-computing

Governance doctrine revision 14 applies.

This private repository is the single Git and Codex home for the ePC-SAFT
project. Read `governance/GOVERNANCE.md`, `governance/manifests/management.json`,
and the nearest subtree `AGENTS.md` before material work.

One scientific slice has one accountable agent from input to evidence. That
agent may change the exact files needed across `data/`, `engine/`, and
`validation/` in one branch. Use `lab/` only for temporary non-production
experiments. Keep one writer per file and do not create nested Git
repositories, submodules, sibling-repository dependencies, or separate package
installations.

Directory ownership remains strict:

- `governance/`: doctrine, topology, manifests, ledgers, and operation receipts.
- `engine/`: the sole `epcsaft` runtime, including EOS, equilibrium, regression,
  compact numerical evidence, and one native build.
- `data/`: canonical reusable scientific packets and packet readiness checks.
- `validation/`: installed-wheel literature, experimental, and cross-module
  scientific campaigns.
- `lab/`: temporary non-production experiments, source notes, and rejected
  alternatives; it owns no persistent tests or test frameworks.

The root `.codex/environments/environment.toml` is the only Codex environment.
It creates one root `.venv`, builds `engine/` into one non-editable wheel, and
rejects split or legacy ePC-SAFT distributions. Use `tools/check-affected.sh`
for normal work, `tools/check-all.sh` for full verification, and
`tools/run-slice.sh khudaida|figiel` for the two canonical vertical slices.
Normal Engine checks route to affected EOS, Equilibrium, or Regression
evidence; full verification also runs the retained slow numerical tier.

Use Numerical Verification-Driven Development from Governance section 7.
Persistent tests must falsify numerical claims; expensive installed-artifact
replays belong in Validation; builds, lint, typing, schemas, hashes, manifests,
environments, and path audits are readiness checks. Do not add contract,
compatibility, legacy, snapshot, prose, mock, or coverage-only tests.

Historical plans, receipts, results, evidence, source locators, old repository
names, and old package names are provenance and must not be normalized as live
instructions. This project has no CI/CD. Run the applicable root checks and the
repository cleanup audit before handoff.

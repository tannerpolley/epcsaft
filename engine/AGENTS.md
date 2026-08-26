# Engine Directory Agent Contract

Read `../governance/GOVERNANCE.md`, `../governance/manifests/management.json`,
root `AGENTS.md`, and `CONTEXT.md` before material work. For CAP work, also read
`docs/agents/issue-tracker.md` before selecting or executing a leaf.

This directory owns the sole `epcsaft` distribution and its public
`epcsaft`, `epcsaft.equilibrium`, and `epcsaft.regression` modules. Keep
one Python build, one Meson graph, one version, one wheel, and one installed
module origin. Do not add compatibility imports, alternate package identities,
installed-provider discovery, hidden physical defaults, duplicate equation
owners, or a public C++ ABI without a real external consumer.

An Engine-accountable task owns its full scientific slice across the exact
needed `data/`, `engine/`, and `validation/` files in one branch. Engine
runtime and tests never discover sibling source dynamically. The containing
monorepo commit binds copied Data fixtures to their packet lock; every
materialization also requires the packet fingerprint and exact file hashes.
Validation must exercise the installed wheel through public APIs.

Use Governance Numerical Verification-Driven Development. Retain compact tests
only for numerical values, identities, derivatives, limits, solver or
equilibrium certificates, manufactured solutions, metamorphic relations, and
source-backed benchmarks. No API-contract, compatibility, legacy, mock,
snapshot, prose, solver-label, iteration-count, or coverage-only tests.

The sole environment is `../.codex/environments/environment.toml`. Run the
applicable root affected/all checks. `codex-cleanup.sh` is an agent closeout
action: agents run it, never ask the user to run it, and remove only task-created
ignored caches or exact user-authorized candidates. Historical plans, receipts,
results, and old names are immutable provenance. No CI/CD.

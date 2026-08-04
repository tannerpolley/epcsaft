# ePC-SAFT Governance Repository Context

`GOVERNANCE.md` revision 8 is the canonical ecosystem doctrine.

`governance_source_of_truth: true`

`management_source_of_truth: false`

`coordination_project: management`

The non-Git `ePC-SAFT-project` parent hosts the Management coordination task.
Governance is the durable Git-tracked authority repository; Management is not
a second repository or source of truth.

`runtime_source_of_truth: false`

GitHub home: `ePC-SAFT/ePC-SAFT-governance`.

Local checkout: `ePC-SAFT-governance`.

The public `ePC-SAFT/.github` repository is a thin organization profile and
community-default surface. It links to Governance but owns no independent
doctrine, project state, or execution plan.

Governance owns current repository topology, worktree access, artifact and
environment inventories, scientific-packet ownership policy, draft disposition,
and external-operation receipts. The Management coordination task handles
exact-subject sequencing within those records.

D-023 approves a direct cutover to `ePC-SAFT/ePC-SAFT-engine` as the sole
runtime repository and `epcsaft` as the sole runtime distribution. EOS,
Equilibrium, and Regression become internal Engine modules under
`epcsaft`, `epcsaft.equilibrium`, and `epcsaft.regression`. The current three
runtime repositories remain transition sources until one accepted cutover
receipt changes authority and retires the separate distributions.

`manifests/runtime-consolidation.json` records the approved target, transition
freeze, consumer inventory, and execution gates. The canonical execution plan
is `docs/plans/2026-08-04-unified-engine-direct-cutover.md`.

`manifests/parameter-ownership.json` is the machine-readable parameter
disposition map. It separates canonical Data storage from Validation evidence
and records that the current EOS ships no installed parameter catalog. Generic
property packets and campaign-local Validation bundles do not enter this
parameter ledger.

`manifests/data-packet-locks.json` pins each reusable bundle to the merged
Data commit, packet path, version, and canonical packet fingerprint. It is a
control-plane lock, not a runtime dependency; Engine and Validation add their
own materialization records when their owning tasks are ready.

`ePC-SAFT/ePC-SAFT-data` is the canonical store for reusable scientific packets.
It is neither a runtime dependency nor the artifact store. Consumers pin an
immutable Data commit and remain self-contained.

`tannerpolley/ePC-SAFT-migration` is a completed historical control repository.
Its Git history preserves legacy-to-clean-package cutover designs and receipts;
it supplies no current authority or workflow.

`ePC-SAFT/ePC-SAFT-lab` is the permanent non-production experiment and
provenance home. The adjacent local `ePC-SAFT-legacy` container and archived
`tannerpolley/ePC-SAFT-legacy-monorepo` preserve historical code and tracker
state outside the active authority topology.

Repository identity and Python identity are independent. The cutover renames
the EOS repository to `ePC-SAFT/ePC-SAFT-engine`; the distribution remains
`epcsaft`. The separate Equilibrium and Regression distributions and import
names receive no compatibility aliases.

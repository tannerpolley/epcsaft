# ePC-SAFT Governance Repository Context

`GOVERNANCE.md` revision 9 is the canonical ecosystem doctrine.

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

D-023 completed the direct cutover to `ePC-SAFT/ePC-SAFT-engine` as the sole
runtime repository and `epcsaft` as the sole runtime distribution. EOS,
Equilibrium, and Regression are Engine capabilities exposed under `epcsaft`,
`epcsaft.equilibrium`, and `epcsaft.regression`. The former EOS name redirects
to Engine; the separate Equilibrium and Regression repositories are archived
provenance and their distributions are inactive.

D-024 authorizes that package-authority cutoff now while Lithium Extraction,
MEA Thermodynamics, and MEA Absorption Column are quarantined as non-admitted
consumers. Their migration is deferred; their current Git pins and workflows
carry no accepted Engine claim and are outside this cutoff task's write scope.

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

Repository identity and Python identity are independent. Repository
`R_kgDOTZ9m0Q` is now `ePC-SAFT/ePC-SAFT-engine`; the distribution remains
`epcsaft`. The retired distributions and import names receive no compatibility
aliases.

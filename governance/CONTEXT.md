# Management Control-Plane Context

`GOVERNANCE.md` revision 6 is the canonical ecosystem doctrine.

`governance_source_of_truth: true`

`management_source_of_truth: true`

`runtime_source_of_truth: false`

GitHub home: `ePC-SAFT/ePC-SAFT-management`.

The public `ePC-SAFT/.github` repository is a thin organization profile and
community-default surface. It links to Management but owns no independent
doctrine, project state, or execution plan.

Management owns current repository topology, worktree access, exact-subject
coordination, artifact and environment inventories, scientific-packet
ownership policy, draft disposition, and external-operation receipts.

`manifests/parameter-ownership.json` is the machine-readable parameter
disposition map. It separates canonical Data storage, admitted EOS runtime
copies, and Validation evidence. Pending relocations remain explicit until the
owning repository tasks execute them.

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

Repository identity and Python identity are independent. The EOS GitHub
repository is `ePC-SAFT/ePC-SAFT-eos`; its distribution and import remain
`epcsaft`.

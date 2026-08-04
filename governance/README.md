# ePC-SAFT Governance

Repository role: `governance_control_plane`

GitHub home: `ePC-SAFT/ePC-SAFT-governance`

Local checkout: `ePC-SAFT-governance`.

Status: canonical ecosystem doctrine and governance home.

Read [GOVERNANCE.md](GOVERNANCE.md) for doctrine revision 8 and
[manifests/management.json](manifests/management.json) for the repository and
access topology.

D-023 approves a direct cutover from the separate EOS, Equilibrium, and
Regression repositories and wheels to one `ePC-SAFT-engine` repository and one
`epcsaft` wheel. Read
[the canonical cutover plan](docs/plans/2026-08-04-unified-engine-direct-cutover.md)
and [the transition manifest](manifests/runtime-consolidation.json) before
runtime, artifact, Validation, or downstream-consumer work.

Governance owns:

- permanent ecosystem doctrine;
- repository roles and naming;
- cross-repository worktree access;
- coordinated project and external-operation sequencing;
- artifact, environment, scientific-packet ownership, and retained-draft ledgers; and
- reusable synchronization and preflight policy.

The non-Git parent project hosts the Management coordination task. It may
sequence work against these records but is not a repository or authority
source.

Artifact layout and admission rules are documented in
[docs/artifact-store.md](docs/artifact-store.md). Cross-repository consumers
install only full-SHA-256 objects; repository `dist/` wheels are rejected.
Non-active candidate wheel retention is recorded in
[manifests/artifact-dispositions.json](manifests/artifact-dispositions.json).

Scientific-packet storage, runtime-copy, evidence, and pending disposition
ownership are tracked in
[manifests/parameter-ownership.json](manifests/parameter-ownership.json).
Immutable Data packet commits and fingerprints are tracked in
[manifests/data-packet-locks.json](manifests/data-packet-locks.json); Engine and
Validation add their own materialization locks when their owning tasks are
ready.

It owns no package runtime, scientific solver, regression target, Validation
campaign, or Lab experiment.

Run the local snapshot and consistency check with:

```bash
python3 tools/management_preflight.py snapshot
python3 tools/management_preflight.py check
```

The former `tannerpolley/ePC-SAFT-migration` repository is archived historical
evidence. The public `ePC-SAFT/.github` repository is intentionally thin and
links here for current authority.

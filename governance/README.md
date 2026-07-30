# ePC-SAFT Management

Repository role: `management_control_plane`

GitHub home: `ePC-SAFT/ePC-SAFT-management`

Status: canonical ecosystem doctrine and management home.

Read [GOVERNANCE.md](GOVERNANCE.md) for doctrine revision 5 and
[manifests/management.json](manifests/management.json) for the repository and
access topology.

Management owns:

- permanent ecosystem doctrine;
- repository roles and naming;
- cross-repository worktree access;
- coordinated project and external-operation sequencing;
- artifact, environment, parameter, and retained-draft ledgers; and
- reusable synchronization and preflight policy.

Artifact layout and admission rules are documented in
[docs/artifact-store.md](docs/artifact-store.md). Cross-repository consumers
install only full-SHA-256 objects; repository `dist/` wheels are rejected.

Parameter-bundle ownership and the pending package-rooted dispositions are
tracked in
[manifests/parameter-ownership.json](manifests/parameter-ownership.json).

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

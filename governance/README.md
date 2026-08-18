# Governance

Canonical doctrine and machine-readable control records for the private
`tannerpolley/ePC-SAFT` monorepo.

- Read [GOVERNANCE.md](GOVERNANCE.md) for revision 22 and D-037.
- Read [manifests/management.json](manifests/management.json) for directory
  ownership and external-state locations.
- Read [the monorepo cutover plan](docs/plans/2026-08-05-monorepo-cutover.md)
  for the history-preserving migration.

Governance owns no scientific tests. Run its readiness checks from the
repository root:

```bash
python3 governance/tools/management_preflight.py check
python3 governance/tools/artifact_store.py check
```

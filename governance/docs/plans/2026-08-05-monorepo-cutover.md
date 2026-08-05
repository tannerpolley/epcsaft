# Private monorepo direct cutover

Status: executing

Authority: user-approved 2026-08-05; doctrine decision D-026.

Target: private `tannerpolley/ePC-SAFT` at
`~/Workspaces/Engineering/ePC-SAFT-project`.

## Phases 1-10

1. Consolidate every source repository to clean, synchronized `main`.
2. Pass source-owned numerical, campaign, packet, native, and readiness gates.
3. Export open issues, repository identities, heads, trees, tags, and settings.
4. Create and verify complete Git bundles and SHA-256 manifests.
5. Rewrite each source `main` history into its owned directory with pinned
   `git-filter-repo`; retain old-to-new commit maps and namespace tags.
6. Merge the five unrelated histories and prove raw subtree tree equivalence.
7. Normalize live doctrine, paths, one root environment, and vertical-slice
   commands without rewriting historical receipts or results.
8. Pass full root checks, Khudaida and Figiel slices, native checks, cleanup,
   Git fsck, and a fresh-clone rehearsal.
9. Create the private destination, push verified history, recreate issues with
   source provenance, and replace the canonical local path recoverably.
10. Archive—not delete—the source repositories and record the final state.

## Stop and rollback rules

Do not archive a source until the private destination and fresh clone pass.
Do not replace the canonical local path until the old parent is preserved as a
timestamped backup. Complete bundles, issue exports, source manifests, commit
maps, and artifacts remain under `ePC-SAFT-local-state`. A failed validation
returns to the last clean staging commit; it never mutates source history.

## Receipt

The final tracked receipt is
`governance/receipts/monorepo-cutover-completed.json`; external bundles and
large local artifacts remain outside Git with SHA-256 manifests.

# Downstream application workspaces

Application work for the Engine roadmap stays in independent clones or worktrees
under this directory. The owner's 2026-09-21 instruction leaves the main local
application checkouts outside `downstream/` unchanged. Publication goes directly
to the corresponding GitHub repository under the task's commit/push/merge authority;
it does not synchronize a main local application checkout.

| Application | GitHub repository | Input and migration contract |
|---|---|---|
| MEA-Thermodynamics | [repository](https://github.com/tannerpolley/MEA-Thermodynamics) | `docs/coordination/greenfield-migration.md` |
| MEA-Absorption-Column | [repository](https://github.com/tannerpolley/MEA-Absorption-Column) | `docs/coordination/greenfield-migration.md` |
| Lithium_Extraction | [repository](https://github.com/tannerpolley/Lithium_Extraction) | `docs/coordination/greenfield-migration.md` |

The Engine tracks only this guide, not application source, Git repositories,
environments or build caches. There are no submodules or runtime sibling-source
imports. A fresh Engine worktree therefore contains the guide; when an application
lane needs source, clone its GitHub repository into `downstream/<application>` or
create a worktree from an existing clone at `downstream/.worktrees/<lane>` using
the verified GitHub `main` as the base. Use one branch and one writer per lane.

The original `codex/greenfield-migration-draft` clones contain retained local edits.
Preserve them and their `copy-provenance.json` record. Their `origin` remotes point
to the main local application checkouts; do not fetch, pull or push those remotes.
A separately verified `github` remote is used for GitHub publication and current
base discovery. A fresh direct GitHub clone can use its normal `origin` remote.
Check the URL, branch, worktree occupancy and dirty state before every Git mutation;
do not reset a retained draft to make it match a newer base.

Read the selected application's `AGENTS.md`, `docs/scientific/CONTEXT.md` and
migration contract before work. GitHub issues own the live scope, authorization,
dependencies and execution queue; the Engine's
[scientific context](../docs/scientific/CONTEXT.md) owns shared responsibilities.
An input-contract checkpoint does not authorize numerical-kernel implementation,
application adoption or qualification. Update the application specification with
any separately authorized implementation.

All numerical work consumes one explicitly identified, non-editable Engine wheel.
Record its SHA-256 and retain the inputs and claim limits. Serialize builds and
installations across lanes; worktree isolation does not isolate the shared wheel.

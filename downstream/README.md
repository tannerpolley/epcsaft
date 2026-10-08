# Downstream application workspaces

The three official applications consume the Engine. Each application's own thread
owns its repository; the Engine asks for application work or evidence by opening an
issue in that repository and edits no application repository. Engine agents may read
application source in independent clones under this directory. The owner's
2026-09-21 instruction leaves the main local application checkouts outside
`downstream/` unchanged.

| Application | GitHub repository | Input and migration contract |
|---|---|---|
| Amine-Thermodynamics | [repository](https://github.com/tannerpolley/Amine-Thermodynamics) | `docs/coordination/greenfield-migration.md` |
| MEA-Absorption-Column | [repository](https://github.com/tannerpolley/MEA-Absorption-Column) | `docs/coordination/greenfield-migration.md` |
| Lithium_Extraction | [repository](https://github.com/tannerpolley/Lithium_Extraction) | `docs/coordination/greenfield-migration.md` |

The Engine tracks only this guide, not application source, Git repositories,
environments or build caches. There are no submodules or runtime sibling-source
imports. A fresh Engine worktree therefore contains the guide; when reading needs
application source, clone its GitHub repository into `downstream/<application>`.

The original `codex/greenfield-migration-draft` clones contain retained local edits.
Preserve them and their `copy-provenance.json` record. Their `origin` remotes point
to the main local application checkouts; do not fetch, pull or push those remotes.
A fresh direct GitHub clone can use its normal `origin` remote. Do not reset a
retained draft to make it match a newer base.

Read the selected application's `AGENTS.md`, `docs/scientific/CONTEXT.md` and
migration contract before work. GitHub issues own the live scope, authorization,
dependencies and execution queue; the Engine's
[scientific context](../docs/scientific/CONTEXT.md) owns shared responsibilities.
An input-contract checkpoint does not authorize numerical-kernel implementation,
application adoption or qualification.

All numerical work consumes one explicitly identified, non-editable Engine wheel.
Record its SHA-256 and retain the inputs and claim limits. Pin
`build/wheels/main-<commit>/epcsaft-*.whl` in the main Engine checkout:
`tools/refresh-main-wheel` keeps one copy per refreshed main commit and never
replaces it. Never pin `~/.cache/epcsaft` or `build/<task>/`; both are scratch.

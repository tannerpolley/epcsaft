# Governance Directory Agent Contract

Read `GOVERNANCE.md` revision 19 and `manifests/management.json` before
project-wide planning or authority changes. Root `AGENTS.md` governs execution.

This directory owns doctrine, directory topology, policy, manifests, ledgers,
readiness checks, and external-operation receipts. It owns no scientific test,
runtime, solver, regression target, Data packet, Validation campaign, or Lab
experiment. Management is a root coordination role, not another repository.

D-026 makes `tannerpolley/ePC-SAFT` the sole active Git home and the root Codex
environment the sole environment. Preserve historical plans, decisions,
receipts, old repository names, package names, hashes, and source locators as
provenance. Current records use `governance/`, `engine/`, `data/`,
`validation/`, and `lab/` ownership.

Governance checks are readiness checks, never scientific tests. Keep
`tools/management_preflight.py` and `tools/artifact_store.py` fail-closed and
worktree-relative. External GitHub mutations require explicit user authority
and a receipt. This project has no CI/CD.

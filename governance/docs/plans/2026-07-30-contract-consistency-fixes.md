# Ecosystem contract consistency fixes

Status: the fix commits were fast-forwarded to each affected repository's
`origin/main` on 2026-07-30; local package-checkout reconciliation remains.

This bounded pass corrected the control-plane and workflow findings that could
be fixed without rewriting active package branches:

- Management now records every EOS parameter bundle in
  `manifests/parameter-ownership.json`, validates that inventory, rejects
  forbidden package installs, and requires Equilibrium and Regression to use
  the active EOS object from the SHA-256 artifact store.
- The stale `epcsaft-regression` editable install was removed from the EOS
  environment. EOS setup now fails if either sibling distribution is present.
- Validation's old frontend route is explicitly retained historical evidence;
  its current identity points to Management's artifact policy instead of the
  superseded wheel subjects.
- EOS, Equilibrium, Validation, and Lab now describe task-specific,
  target-owned, on-demand Validation worktrees rather than persistent slots.
  Lab also ignores its local `.worktrees/` directory.
- Local doctrine and access-manifest paths consistently resolve the current
  `ePC-SAFT-organization` checkout alias until the public Management rename is
  reflected locally.

## Preserved local commits

| Repository | Branch | Commit | Worktree |
| --- | --- | --- | --- |
| Management | `codex/ecosystem-contract-consistency` | `453b35f`, `41ab9a0`, `d98e719` plus this handoff record | attached checkout |
| EOS | `codex/eos-contract-consistency` | `a8d618f` | `.worktrees/codex-eos-contract-consistency` |
| Equilibrium | `codex/equilibrium-contract-consistency-v2` | `2d24ae4` | `.worktrees/codex-equilibrium-contract-consistency-v2` |
| Regression | `codex/regression-contract-consistency` | `5aaf79f` | `.worktrees/codex-regression-contract-consistency` |
| Validation | `codex/validation-contract-consistency` | `243edfb` | `.worktrees/codex-validation-contract-consistency` |
| Lab | `codex/lab-contract-consistency` | `7c31cbf`, `9c32ed7` | `.worktrees/codex-lab-contract-consistency` |

The commits above are now published on their repositories' `main` branches.
The local feature branches and worktrees remain only as preservation handles
until the dirty package checkouts reconcile with their updated `origin/main`,
or until Management records an explicit cleanup decision.

The earlier Equilibrium branch `codex/equilibrium-contract-consistency`
(`21bb72f`) is a superseded copy based on the pre-completion package head; the
`-v2` branch above is the applicable follow-up.

## Deferred package-owned work

The parameter manifest classifies the Ascani 2022 case-study, Crespo 2021 DEEA,
and Khudaida 2026 Figure 2 bundles as `relocate_to_validation`. Their bytes,
EOS catalog tests, and Validation campaigns still need a coordinated,
package-rooted transfer; this pass does not silently remove runtime data from
EOS or rewrite scientific evidence. The active Validation channel likewise
needs a new package-rooted campaign before any result can be called current.

The Equilibrium and Regression primary checkouts were not modified. Their
completed agent changes, including any still-uncommitted files, remain the
source of truth until their owners reconcile these follow-up commits.

## Verification

- `python tools/management_preflight.py check --strict`
- `python tools/artifact_store.py check`
- Management preflight and artifact-store tests: 8 passed
- Retained Validation route tests: 4 passed
- EOS setup shell syntax check passed
- EOS environment no longer contains `epcsaft-regression`

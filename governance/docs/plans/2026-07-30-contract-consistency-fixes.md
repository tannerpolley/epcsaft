# Ecosystem contract consistency fixes

Status: local feature commits complete; remote publication and merge remain
pending the package agents' active work.

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
| Management | `codex/ecosystem-contract-consistency` | `453b35f` plus this handoff record | attached checkout |
| EOS | `codex/eos-contract-consistency` | `a8d618f` | `.worktrees/codex-eos-contract-consistency` |
| Equilibrium | `codex/equilibrium-contract-consistency` | `21bb72f` | `.worktrees/codex-equilibrium-contract-consistency` |
| Validation | `codex/validation-contract-consistency` | `243edfb` | `.worktrees/codex-validation-contract-consistency` |
| Lab | `codex/lab-contract-consistency` | `7c31cbf`, `9c32ed7` | `.worktrees/codex-lab-contract-consistency` |

These branches are intentionally unpushed. Retain their worktrees until the
owning package agent cherry-picks or merges the commits, or until Management
records an explicit preservation decision.

## Deferred package-owned work

The parameter manifest classifies the Ascani 2022 case-study, Crespo 2021 DEEA,
and Khudaida 2026 Figure 2 bundles as `relocate_to_validation`. Their bytes,
EOS catalog tests, and Validation campaigns still need a coordinated,
package-rooted transfer; this pass does not silently remove runtime data from
EOS or rewrite scientific evidence. The active Validation channel likewise
needs a new package-rooted campaign before any result can be called current.

The dirty Equilibrium and Regression primary checkouts were not modified.
Their active branches remain the source of truth for their in-flight changes.

## Verification

- `python tools/management_preflight.py check --strict`
- `python tools/artifact_store.py check`
- Management preflight and artifact-store tests: 8 passed
- Retained Validation route tests: 4 passed
- EOS setup shell syntax check passed
- EOS environment no longer contains `epcsaft-regression`

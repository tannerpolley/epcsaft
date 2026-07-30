# Management Control-Plane Transition

Status: executing

Decision date: 2026-07-29

Owner: `ePC-SAFT/ePC-SAFT-management`

## Goal

Make Management the single durable control plane for the ePC-SAFT ecosystem,
retire the completed legacy-to-clean-package Migration role, preserve its
historical evidence, and keep `ePC-SAFT/.github` as a thin GitHub organization
profile and community-default repository.

## Fixed topology

| Repository | Durable responsibility |
| --- | --- |
| `ePC-SAFT/ePC-SAFT-management` | Doctrine, repository topology, access policy, project coordination, manifests, artifact and parameter ledgers, synchronization checks, external-operation receipts, and shared workflow policy |
| `ePC-SAFT/.github` | Organization profile and links plus GitHub-recognized community defaults; no independent doctrine or project state |
| `ePC-SAFT/ePC-SAFT-eos` | EOS runtime, `epcsaft` distribution/import, generic parameter contracts, CppAD derivatives, and public SDK |
| `ePC-SAFT/ePC-SAFT-equilibrium` | Equilibrium formulations, solvers, and certificates |
| `ePC-SAFT/ePC-SAFT-regression` | Regression targets, Ceres workflows, diagnostics, and fit results |
| `ePC-SAFT/ePC-SAFT-validation` | Public-installed-artifact and literature evidence authored by the accountable package task |
| `ePC-SAFT/ePC-SAFT-lab` | Non-production experiments and broad provenance |
| `tannerpolley/ePC-SAFT-migration` | Archived historical cutover records only |

The former Migration repository has no active owner, queue, gate, or runtime
role after this transition. Its immutable Git history remains the authority
for its historical receipts. Management records its repository identity and
terminal commit rather than copying thousands of lines of obsolete temporary
control logic.

## Success criteria

1. The current public repository object is renamed from `ePC-SAFT/.github` to
   `ePC-SAFT/ePC-SAFT-management`, preserving its history, issues, pull
   requests, and branch protection.
2. A new public `ePC-SAFT/.github` repository contains only the organization
   profile, community defaults, and links to Management.
3. Doctrine revision 5 names Management as the permanent control owner and
   contains no active Migration role.
4. Every live repository points to Management doctrine and uses EOS rather
   than Provider terminology.
5. Management owns one machine-readable source for:
   - repositories and roles;
   - target-owned worktree access;
   - draft and branch disposition;
   - wheels, headers, native libraries, and environments; and
   - parameter-bundle ownership.
6. One preflight command validates those records against local Git state and
   immutable artifact bytes.
7. Persistent worktrees are target-owned access slots; feature worktrees are
   temporary and carry an explicit lifecycle state.
8. Historical or blocked drafts remain reachable and are never confused with
   accepted evidence.
9. The former Migration repository has a terminal retirement commit, its open
   work is transferred or closed with links, and the GitHub repository is
   archived.
10. Canonical checkouts are clean and synchronized after coordinated merges.
11. Every active repository contains a checked-in Codex local environment
    whose setup script prepares a fresh app-created worktree without relying
    on a canonical-checkout path or mutable sibling output.

## Access policy

```text
Management  -> all repositories for control-plane work
Equilibrium -> Management, EOS, Equilibrium, Regression, Validation, Lab
Regression  -> Management, EOS, Equilibrium, Regression, Validation, Lab
EOS         -> EOS, Validation
Validation  -> Validation only
Lab         -> Lab only
.github     -> .github only
```

Access means permission to create a task branch in a target-owned worktree. It
does not transfer runtime or evidence ownership.

## Persistent target-owned worktrees

```text
ePC-SAFT-eos/.worktrees/equilibrium-access
ePC-SAFT-eos/.worktrees/regression-access
ePC-SAFT-equilibrium/.worktrees/regression-access
ePC-SAFT-regression/.worktrees/equilibrium-access
ePC-SAFT-validation/.worktrees/eos-access
ePC-SAFT-validation/.worktrees/equilibrium-access
ePC-SAFT-validation/.worktrees/regression-access
```

Management work uses an ordinary feature worktree in the target repository.
Validation has no outbound slot and no resident agent.

## Artifact and environment contract

- Cross-repository consumers use immutable commit- and SHA-256-addressed
  wheels, never another repository's mutable `dist/`.
- The ledger records distribution, version, producer commit/tree when known,
  wheel hash, embedded-header hash, static-library hash, consumers, receipts,
  and retention state.
- Repository environments contain only the editable owner and exact immutable
  dependencies.
- Validation creates disposable exact-subject environments.
- Unknown artifacts are inventoried before they are retained or removed.
- `.codex/environments/environment.toml` and its setup script are tracked in
  every active repository. Setup scripts are idempotent, worktree-relative,
  non-publishing, and use only declared immutable cross-repository artifacts.

## Parameter contract

EOS ships only an explicit Management-approved allowlist of model-defining
PC-SAFT/ePC-SAFT parameterizations. Case-study, fitted, alternative,
reproduction, rejected, and superseded bundles belong to their accountable
Validation campaigns and enter EOS through public parameter-input contracts.
Validation is evidence data, never a production catalog.

## Draft lifecycle

Every non-main branch or pull request is exactly one of:

- `active`;
- `ready`;
- `preservation-draft`;
- `blocked`;
- `merged-awaiting-cleanup`;
- `superseded`; or
- `retired`.

Preservation drafts bind an issue, exact commit, reason, and evidence
limitations. They are not rebased merely to remove a conflict.

## Execution order

1. Publish Management doctrine, manifests, tooling, and this transition record.
2. Rename the current repository object and recreate the thin `.github`.
3. Update EOS, Equilibrium, Regression, Validation, and Lab pointers and
   terminology.
4. Record and normalize worktrees without deleting unrecorded state.
5. Inventory artifacts, environments, and parameter bundles.
6. Retire and archive the former Migration repository.
7. Run the Management preflight, merge coordinated changes, synchronize
   canonical checkouts, and publish the terminal management receipt.

## Stop conditions

Stop before an affected operation when:

- a repository identity or default-branch commit differs from the captured
  preflight;
- a dirty path or worktree lacks an owner;
- a branch or artifact is not reachable from a recorded commit;
- an accepted receipt would be rewritten;
- an external rename would drop organization-profile behavior;
- a parameter or artifact cannot be assigned one owner and lifecycle state;
- a force push, destructive reset, or unapproved deletion would be required;
  or
- verification cannot distinguish current authority from historical evidence.

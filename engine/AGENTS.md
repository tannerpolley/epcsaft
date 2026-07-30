# EOS Repository Agent Contract

Repository Profile: scientific-computing

Before work, read doctrine revision 5 at
`../ePC-SAFT-management/GOVERNANCE.md` and its machine-readable topology at
`../ePC-SAFT-management/manifests/management.json`. During the local checkout
transition, resolve both files from `../ePC-SAFT-organization/`. Published doctrine lives at
`https://github.com/ePC-SAFT/ePC-SAFT-management/blob/main/GOVERNANCE.md`.
Local policy may narrow that doctrine but must not contradict it.

This repository owns one state/property lifecycle:

```text
ParameterBundle -> ParameterSet -> EPCSAFT.evaluate -> EosResult
```

It also owns one narrow native consumer boundary:

```text
EPCSAFT -> native_sdk(model) -> one stateless pure phase block
```

The `epcsaft` distribution and Python import remain named `epcsaft`; the
repository name is `ePC-SAFT-eos`. EOS owns resolved parameters,
EOS/state/property evaluation, density closure, CppAD, and the admitted native
SDK transport. It does not own phase equilibrium or discovery, Ipopt,
regression, Ceres, cross-package validation policy, reaction/species
equilibrium, or project coordination.

Accepted receipts such as `promotion-0012-provider-full-candidate` and
`promotion-0015-provider-native-sdk-v1` retain their historical identifiers.
They govern only their recorded subjects and do not publish a package release
or broaden current capability.

Do not add compatibility APIs, mutable scientific dictionaries, backend
selectors, hidden physical defaults, duplicate equation owners, or runtime
dependencies on the lab, historical Migration repository, or sibling source
trees.

## Validation access

- The EOS task remains rooted in this repository and is accountable for EOS
  implementation, parameters, EOS/state/property behavior, derivatives,
  density closure, and the native SDK.
- Management may assign this task a bounded campaign in the sibling Validation
  repository through a target-owned Validation worktree created on demand under
  Validation's `.worktrees/` directory for a task-specific branch and path.
- Validation work must use exact immutable installed public artifacts. It must
  not import this source checkout, use private seams, copy EOS equations, or
  add production implementation to Validation.
- The EOS task may author and execute its assigned campaign, then returns here.
  It cannot accept its own promotion or broaden EOS authority from its own
  evidence.
- Management serializes Validation writers and dispatches any distinct review
  required for a stable subject. Validation has no permanent resident agent.

## Parameter boundary

EOS may package compact, source-traceable parameter sets from canonical
PC-SAFT or ePC-SAFT papers when they are part of the supported model catalog.
Case-specific, derived, fitted, campaign, or noncanonical parameter bundles
belong in Validation evidence, not in this repository.

## Issue tracker

Issues and PRDs are tracked in GitHub Issues for
`ePC-SAFT/ePC-SAFT-eos`. Use the five default triage labels:
`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, and
`wontfix`. See `docs/agents/`.

This is a single-context repository with root `CONTEXT.md`, `ARCHITECTURE.yaml`,
and `docs/adr/`.

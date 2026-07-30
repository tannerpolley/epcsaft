# ePC-SAFT Ecosystem Doctrine

Status: canonical

Doctrine revision: 5

Approved: 2026-07-29

Owner: `ePC-SAFT/ePC-SAFT-management`

This document is the permanent architecture and governance authority for the
ePC-SAFT repository ecosystem. `MUST`, `MUST NOT`, `SHOULD`, and `MAY` are
normative. Repository-local policy may narrow these rules but must not
contradict them.

## 1. Optimization objective

The ecosystem MUST minimize its authored and maintained surface while still
producing every result admitted by executable capability evidence.

The feasible set is constrained by:

- source-faithful equations, parameters, units, and conventions;
- reproducible numerical results within declared tolerances;
- exact production derivatives required by the admitted solver route;
- buildable, installable, readable, and traceable software;
- explicit capability boundaries and early rejection of unsupported work; and
- retained provenance for scientific and management decisions.

Within that feasible set, prefer fewer owners, public symbols, dependencies,
native targets, files, and lines. Code golf, compressed control flow, giant
files, generated indirection, and deletion of necessary evidence are not
optimizations.

This doctrine defines admission rules, not current capability. A repository
may claim only what its executable capability registry and accepted receipts
prove.

## 2. Authority hierarchy

Authority descends in this order:

1. this organization doctrine;
2. a repository's current `CONTEXT.md` and `AGENTS.md`;
3. accepted repository-local architecture decisions;
4. approved implementation plans; and
5. machine-validated transfer, validation, and release receipts.

Lower authorities MUST NOT contradict higher ones. A new decision that changes
an accepted rule MUST name what it supersedes. Historical documents remain
evidence; they do not silently regain authority.

Governance authority and runtime authority are separate. Accepting this
doctrine does not transfer any scientific implementation or capability.

## 3. Repository ownership

| Home | Sole ecosystem role | Must not own |
|---|---|---|
| `ePC-SAFT/ePC-SAFT-management` | Doctrine, repository topology, Project coordination, access policy, manifests, ledgers, synchronization checks, and external-operation receipts | Package source, package tests, scientific capabilities, or releases |
| `ePC-SAFT/.github` | Organization profile and GitHub-recognized community defaults linking to Management | Independent doctrine, project state, package source, or scientific capability |
| `ePC-SAFT/ePC-SAFT-eos` | `epcsaft` EOS: resolved inputs, allowlisted reference parameters, EOS/state/property evaluation, density closure, CppAD substrate, and public SDK | Ipopt equilibrium, Ceres regression, or cross-package validation campaigns |
| `ePC-SAFT/ePC-SAFT-equilibrium` | `epcsaft-equilibrium`: equilibrium formulations, phase discovery, local NLPs, Ipopt, certification, and equilibrium results | EOS implementation, Ceres, or regression targets |
| `ePC-SAFT/ePC-SAFT-regression` | `epcsaft-regression`: target contracts, parameter fitting, Ceres, fit diagnostics, and regression results | EOS implementation, Ipopt, or equilibrium dispatch |
| `ePC-SAFT/ePC-SAFT-validation` | Installed-artifact, cross-package, and literature acceptance evidence authored by the accountable package task | Production algorithms, private package imports, a paper archive, or a permanent execution owner |
| `ePC-SAFT/ePC-SAFT-lab` | Non-production experiments, manufactured and independent reference oracles, broad provenance, source indexes, and rejected alternatives | Production distributions or releases, runtime authority, admission receipts, or dependencies consumed by production or validation repositories |

One concept has one runtime owner. A capability without an accepted production
owner is non-production; no lab or legacy copy supplies transitional runtime
authority. The lab MAY investigate such a capability, but production and
validation repositories MUST NOT depend on it. Promotion distills a clean,
independently evidenced slice into its unique package owner rather than
transferring authority from a runnable archive.

The complete pre-migration monorepo and related historical repositories live
outside the authority topology under the local `ePC-SAFT-legacy` container.
The public `tannerpolley/ePC-SAFT-legacy-monorepo` is a read-only historical
repository. Legacy code, tracker state, branches, tags, and releases are
provenance only; they own no current runtime, roadmap, or governance authority.

Repository ownership and task residence are distinct. Validation is the
durable evidence home, not a standing agent or implementation team. EOS,
Equilibrium, or Regression remains accountable for its capability and performs
the assigned black-box campaign in the Validation repository from a task rooted
in its own package repository. Management sequences exact subjects and
prevents concurrent writers; it does not maintain a waiting Validation worker.

The completed `tannerpolley/ePC-SAFT-migration` repository is historical
evidence outside the active topology. It owns no current plan, queue, gate,
decision, or authority.

## 4. Dependency direction

The only production dependency edges are:

```text
ePC-SAFT-equilibrium ----> ePC-SAFT-eos
ePC-SAFT-regression -----> ePC-SAFT-eos
ePC-SAFT-validation -----> installed artifacts from admitted repositories
```

EOS depends on neither extension. The extensions do not depend on
each other. No production or validation repository depends on the lab.

Cross-repository builds MUST use installed artifacts. They MUST NOT use sibling
source paths, compile another repository's implementation sources, copy an
owner's kernel, or reach through private headers or Python modules.

The EOS SDK MUST be narrow and versioned. Public C++ ABI stability is out
of scope unless a real external C++ consumer is admitted. The native transport
mechanism requires its own accepted architecture decision before an extension
runtime is transferred; this doctrine does not preselect a capsule, C ABI,
shared-library ABI, or header-instantiation design.

Any earlier contract that preselects source/CMake transport, recompilation of
EOS implementation sources, or another native mechanism is experimental
evidence only. It is not a transfer contract and cannot admit a native
consumer slice.

## 5. Package shape

A production package starts from the smallest complete vertical slice. It does
not prebuild a speculative tree. Directories appear when accepted owners need
them. Placeholder modules, placeholder tests, empty native targets, `.gitkeep`
files, compatibility layers, fake defaults, and release workflows without a
release are forbidden.

The default package-first hybrid shape is:

```text
.github/workflows/
cpp/include/<owner>/
cpp/src/
cpp/bindings/
src/<python_package>/
tests/
docs/
CMakeLists.txt
pyproject.toml
README.md
LICENSE
AGENTS.md
CONTEXT.md
```

Only needed entries are materialized. Begin with one Python distribution, one
Python package, one native extension, and one native target per repository.
Add another only when a demonstrated ownership or packaging constraint
requires it.

Python validates and serializes inputs, calls native workflows, and presents
typed results. Native C++ owns production thermodynamic kernels and production
solver loops. Public APIs MUST be small, typed, generic, and free of backend
selection flags, aliases, application-specific routes, and duplicate input or
result representations.

## 6. Thermodynamic and derivative doctrine

EOS's central contract is a source-traceable, scalar-generic residual
Helmholtz implementation over resolved thermodynamic inputs. Pressure, fugacity
and chemical-potential quantities, density closure, equilibrium, and
regression MUST derive from that owner rather than parallel equation kernels.

CppAD is the sole production authority for nonlinear derivative tensors.
Handwritten nonlinear derivative kernels and numerical production derivative
backends are forbidden. Exact linear transformations and exact implicit-
function linear algebra assembled from CppAD-generated partial tensors are
allowed; they are not alternate derivative modes.

Production differentiation MUST tape defining equations, residuals,
objectives, and constraints for a fixed discrete topology. It MUST NOT tape an
iterative solver history.

For an eliminated state defined by `H(u, z) = 0`, the value solve occurs outside
the tape. Required sensitivities use CppAD partial tensors and exact linear
solves, beginning with:

```text
u_z = -H_u^{-1} H_z
```

An admitted route MUST provide every mixed and higher-order implicit
sensitivity its solver contract consumes. First-order sensitivity alone does
not prove an exact Hessian. If an equilibrium Hessian requires higher-order
EOS tensors, the EOS contract MUST supply them through CppAD. A
route MAY instead lift the state into the optimization variables with its full
defining constraints and derivative contract.

Density sensitivities are valid only on a preserved, well-conditioned root
branch. Evidence MUST identify the bracket or branch, root residual,
`dR/d(rho)` conditioning, and rejection behavior for root swaps, coalescence,
or topology change.

Association site fractions are solved states. Eliminated association requires
complete implicit derivatives. Lifted association requires true site topology
and every density, composition, temperature, and association-strength
dependency. A fixed-association-strength block is diagnostic only and cannot
support a production capability or exact-Hessian claim.

Ipopt routes MUST provide exact objective and constraint first derivatives and
the exact Lagrangian Hessian for the admitted formulation. Ceres routes MUST
provide exact residual Jacobians for each admitted target family; this does not
imply an exact objective Hessian or generic regression support.

Independent numerical or high-precision derivative checks are encouraged as
validation oracles. They MUST NOT be selectable production backends.

Every result contract distinguishes three statuses:

1. solver convergence: accepted termination with finite, scaled residual and
   feasibility diagnostics;
2. numerical convergence: stability to the relevant tolerances, precision,
   scaling, initialization, continuation, or search resolution; and
3. physical validity: source, units, domain, branch or topology, stability,
   and formulation-specific evidence support the claimed interpretation.

No status implies a later one. A capability record MUST identify its exact
input domain, parameter and source fingerprints, topology, derivative orders,
local or global scope, stability or discovery status, certificate, tolerance
set, accepted artifact, and explicit failure boundaries.

## 7. Evidence ownership

Package repositories own compact tests for public contracts, native seams,
equation contributions, derivative payloads, solver diagnostics, and a small
set of source-backed scientific anchors. Tests protect behavior, not old file
paths, source wording, private call graphs, or retired transition-era structure.

Package tests MAY own test-only independent numerical or high-precision
oracles for white-box tensors. The validation repository owns test-only
physical and cross-package certificate oracles over installed immutable
artifacts. Neither oracle is runtime API or production input. Shared equations,
source, generated material, or copied kernels between an implementation and
its oracle MUST be declared and MUST NOT count as independent evidence.

The validation repository owns the durable black-box physical and cross-package
evidence. The accountable package task installs wheels or other declared
artifacts built from immutable commits and writes the bounded campaign there;
an artifact need not be publicly released. It MUST NOT import sibling source
trees or reproduce production algorithms.

Task identity alone does not make evidence independent. A package-authored
campaign MUST remain public-API-only, hash-bound, source-traceable, and separate
from implementation internals. Independent acceptance still requires a bounded
reviewer distinct from the producer when the receipt protocol requires one,
but that reviewer is dispatched for the stable subject and MUST NOT become a
permanent Validation worker. Self-validation may establish package evidence;
it cannot accept its own promotion or broaden its own scientific claim.

Black-box agreement alone does not prove derivative exactness, local KKT
quality, phase completeness, or globality. Those require package-owned
white-box receipts and formulation-specific independent certificates.

Every retained scientific dataset needs a source, transformation record,
units, license or use basis, and tolerance rationale. Production repositories
retain only compact owned data. Full papers, broad extraction corpora, failed
experiments, and exploratory notebooks belong in the lab. Retained literature
campaign plots belong in validation evidence, not in every package test.

## 8. Promotion and unique authority

Lab-to-production promotion is distillation, not bulk transfer. Clean history
is the default. Original commits remain cited as provenance.

Every promoted slice requires a machine-validated receipt containing at least:

- capability and owning repository;
- source and destination commits;
- source references, equations, conventions, and units;
- accepted file and public-interface manifests;
- input, output, error, and capability contracts;
- derivative orders and implicit-state treatment;
- build artifacts and hashes;
- tests, datasets, tolerances, and validation results;
- displaced-owner and forbidden-dependency checks; and
- the approved effective authority revision.

All operation, state, promotion, and validation receipts use one versioned
envelope. The envelope records its schema version and hash; monotonic sequence;
validator repository, commit, path, and hash; exact command and environment;
artifact hashes; producer and independent reviewer; retained failed outcomes;
supersession or revocation; and user approval bound to the exact receipt
subject hash. Receipts MUST be replayable from repository identity, commit,
tree, and immutable artifacts without requiring the original workstation path.
Validation MUST prove that the artifact under test did not import the source
checkout being validated.

A receipt cannot accept itself. Its validator MUST pass, the producer and
reviewer MUST differ for acceptance, and the user MUST approve the authority
change. Before that moment the candidate is non-production. After it, the
destination is the sole production owner and any lab or legacy copy remains
non-authoritative provenance.

The first runtime promotion is the smallest complete EOS thermodynamic
path. Later slices MUST reuse its public and native contracts rather than add
parallel paths.

## 9. Architecture ratchet

Each production repository establishes its first architecture baseline only
after its first complete accepted slice. The baseline separately records:

- authored runtime, test, and documentation lines;
- tracked file count and generated-file policy;
- public Python exports and native SDK entries;
- direct runtime and build dependencies;
- native modules and targets; and
- oversized source files.

The ratchet is a review constraint, not a code-golf score. An increase requires
a named capability or maintainability justification and an explicit decision
to accept the cost. New work MUST remove obsolete source, tests, documentation,
and configuration in the same checkpoint. Empty skeletons are not baselines.

## 10. GitHub, issues, and releases

Production and Management repositories are public. They may remain
unreleased until a separate release decision exists. Versions are independent;
the historical `epcsaft` distribution identity remains reserved, but no 1.0 or
release workflow is implied by this doctrine.

Issues are owned by the repository that owns the work. One Management
Project may aggregate them. Cross-repository references MUST use full URLs or
stable node identifiers, never ambiguous number-only references.

The public `.github` repository owns only the organization profile and
GitHub-recognized community defaults. Management owns reusable workflow policy
and project coordination. Package repositories still own their caller
workflows, package checks, build matrices, and releases.

`main` MUST reject force pushes and deletion. Normal changes use review and
required evidence when those checks exist. A direct-main or administrator
bypass requires explicit user approval and a retained receipt; bypass is not a
standing exemption.

Remote creation, repository rename or transfer, issue transfer, Project
mutation, publication, and release remain explicit external-write gates. This
doctrine authorizes none of them by itself.

Before any repository rename, transfer, or old-name reuse, Management
MUST perform an operation-specific read-only preflight, capture immutable
machine-readable exports for metadata at risk, inventory Actions and reusable
workflow references, and record redirect behavior, namespace-reuse
eligibility, metadata restoration, and every accepted loss. A sacrificial
rehearsal is required only when the preflight cannot establish the behavior or
rollback boundary. Each mutation requires a cutover record stating
preconditions, exact request, expected postconditions, retry and stop behavior,
reversibility, irreversible effects, and recovery owner, followed by a state
record proving the result.

## 11. Enforcement and amendment

Every repository `AGENTS.md` and `CONTEXT.md` MUST name the doctrine revision it
follows. Shared checks SHOULD verify the revision, forbidden dependency edges,
receipt schemas, and architecture ratchet. Reusable workflows still require a
small caller workflow in each package repository.

Amending this doctrine requires explicit user approval, a new revision, a
decision-ledger entry, and updates to affected repository pointers. Runtime
authority changes only through accepted promotion receipts.

## Decision Ledger

| Decision | Source | Answer | Impact | Deferred? | Risk owner |
|---|---|---|---|---|---|
| D-001 Repository topology | User-approved topology review, superseded by D-017 | Historical six-home topology used EOS, equilibrium, regression, validation, organization policy, and personal lab | Preserved as historical topology | No | Management |
| D-002 Policy authority | Canonical-policy review, superseded by D-017 | Historical Organization `GOVERNANCE.md` was the permanent policy authority | Preserved as historical authority decision | No | Management |
| D-003 Package shape | Repository comparison and user objective | Use clean, package-first hybrid repositories with no speculative skeleton content | Minimizes maintained surface | No | Each package owner |
| D-004 Dependency direction | Ownership audit | Allow only one-way installed-artifact dependencies; forbid sibling-source and lab dependencies | Prevents duplicate owners and hidden coupling | No | Each package owner |
| D-005 EOS SDK | EOS SDK audit | Require a narrow EOS SDK and decide exact native transport separately; prior source/CMake preselection is experimental only | Blocks native extension transfer until resolved | Yes | EOS owner |
| D-006 Derivative authority | User-approved derivative doctrine | Make CppAD the sole nonlinear production derivative authority | Removes derivative modes and duplicate kernels | No | EOS owner |
| D-007 Implicit derivative order | Equilibrium Hessian sanity audit | Require every implicit and higher derivative order consumed by the solver | Prevents false exact-Hessian claims | No | EOS and solver owner |
| D-008 Evidence boundary | Evidence-boundary audit | Combine package white-box receipts with validation black-box artifact campaigns | Separates implementation proof from physical acceptance | No | Package and validation owners |
| D-009 Promotion authority | Historical Migration ownership audit | Use validated promotion receipts to create one production owner | Preserved historical admission rule; new coordination is Management-owned | No | Management |
| D-010 Architecture baseline | User optimization objective | Baseline after the first complete slice and reject code-golf gaming | Keeps the architecture ratchet meaningful | No | Each package owner |
| D-011 Tracker ownership | Tracker design review | Use package-local issues and one Management Project | Aligns work with code ownership | No | Management |
| D-012 Release timing | Release review | Keep repositories public while deciding release timing and versions separately | Avoids premature release machinery | Yes | Release owner |
| D-013 Lab role | User lab-role decision, superseded by D-016 | The original monorepo served as transitional authority during the legacy cutover | Preserved as historical doctrine; no longer governs current lab or legacy repositories | No | Lab and Management |
| D-014 Pre-execution review remediation | Independent principal review and user-approved adjudication | Require whole-tree cutover inventory, portable hash-bound receipts, independent test-only oracle ownership, and rehearsed external operations | Preserved historical cutover control | No | Management |
| D-015 Package-rooted validation execution | User-approved workflow amendment | Keep Validation as the evidence repository while the accountable EOS, Equilibrium, or Regression task performs its own bounded installed-artifact campaign there; dispatch review only for a stable subject and keep no resident Validation worker | Removes permanent waiting tasks without weakening public-artifact, source, claim, or receipt boundaries | No | Package, Validation-evidence, and Management owners |
| D-016 Lab and legacy separation | User-approved 2026-07-22 repository cutover | Give `ePC-SAFT/ePC-SAFT-lab` fresh history and a permanent non-production experimental role; retain the old monorepo and related projects under `ePC-SAFT-legacy` with no authority | Separates active research from historical code and removes fallback runtime authority from both lab and legacy | No | Management and Lab |
| D-017 Management control plane | User-approved 2026-07-29 topology decision | Rename the current authority repository to `ePC-SAFT-management`, keep `.github` thin, retire Migration, and centralize doctrine, manifests, ledgers, access, synchronization, and external-operation receipts in Management | Removes temporary coordination authority and duplicated project state | No | Management |
| D-018 Repository access | User-approved hierarchy | Equilibrium and Regression may work in every operational repository; EOS may work only in EOS and Validation; Validation and Lab have no outbound access; every worktree is target-owned | Makes task residence explicit without transferring ownership | No | Management |
| D-019 Codex worktree environments | User-approved setup requirement | Track a worktree-relative Codex local environment in every active repository; setup must be idempotent and must not consume mutable sibling build output | Makes new Codex worktrees reproducible | No | Management and repository owners |
| D-020 Ephemeral worktrees and immutable artifacts | User-approved 2026-07-29 cleanup | Create target-owned worktrees only on demand and retire them after merge or recorded preservation; install cross-repository wheels only through the SHA-256 object store and reject mutable repository `dist/` wheels | Prevents stale checkout and SDK accumulation | No | Management and repository owners |

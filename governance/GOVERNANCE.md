# ePC-SAFT Ecosystem Doctrine

Status: canonical

Doctrine revision: 8

Approved: 2026-08-04

Owner: `ePC-SAFT/ePC-SAFT-governance`

The durable Git-tracked authority repository is Governance. The non-Git
`ePC-SAFT-project` parent hosts the Management coordination task and has no
separate repository or governance authority.

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
| `ePC-SAFT/ePC-SAFT-governance` | Doctrine, repository topology, governance policy, manifests, ledgers, synchronization checks, and external-operation receipts | Package source, package tests, scientific capabilities, or releases |
| `ePC-SAFT/.github` | Organization profile and GitHub-recognized community defaults linking to Governance | Independent doctrine, project state, package source, or scientific capability |
| `ePC-SAFT/ePC-SAFT-data` | Canonical reusable scientific packets: parameter bundles, source-derived datasets, normalized observations, provenance, units, domains, hashes, and license or use bases | Runtime algorithms, generated validation results, package artifacts, full papers, or exploratory and rejected material |
| `ePC-SAFT/ePC-SAFT-engine` | The single `epcsaft` runtime distribution: EOS/state/property evaluation, CppAD derivatives, equilibrium formulations and certification, regression targets and diagnostics, one native build, and one release identity | Canonical scientific packets, durable black-box validation evidence, application workflows, or separate Equilibrium and Regression distributions |
| `ePC-SAFT/ePC-SAFT-validation` | Installed-artifact, cross-module, downstream-consumer, and literature acceptance evidence authored by the accountable Engine task | Production algorithms, private Engine imports, a paper archive, or a permanent execution owner |
| `ePC-SAFT/ePC-SAFT-lab` | Non-production experiments, manufactured and independent reference oracles, broad provenance, source indexes, and rejected alternatives | Production distributions or releases, runtime authority, admission receipts, or dependencies consumed by production or validation repositories |

Revision 8 approves `ePC-SAFT/ePC-SAFT-engine` as the target runtime home. The
cutover will rename `ePC-SAFT/ePC-SAFT-eos` and transfer the accepted
Equilibrium and Regression implementations into that repository. Until an
accepted D-023 cutover receipt changes runtime authority, the three current
repositories remain transition sources for cutover work and critical repairs.
They MUST NOT start independent feature or release lines. The transition state
and gates live in `manifests/runtime-consolidation.json`.

One concept has one runtime owner. Canonical scientific storage does not confer
runtime authority. A capability without an accepted production
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
durable evidence home, not a standing agent or implementation team. One Engine
vertical-slice task remains accountable for each capability and performs the
assigned black-box campaign in the Validation repository. The Management
coordination task sequences exact subjects under Governance policy and
prevents concurrent writers; it does not maintain a waiting Validation worker.

The completed `tannerpolley/ePC-SAFT-migration` repository is historical
evidence outside the active topology. It owns no current plan, queue, gate,
decision, or authority.

Data is a read-only scientific-input home, not a service or runtime dependency.
Lab may submit reviewed packets to Data. Engine may copy or generate an
admitted runtime-catalog subset from a pinned Data release. Validation and
Engine regression campaigns may consume packets from exact Data commits for
campaigns and declared fitting inputs. The consuming module or campaign
remains self-contained and records the Data commit, packet path, and packet
fingerprint.

## 4. Dependency direction

The only production dependency edge is:

```text
ePC-SAFT-validation -----> installed ePC-SAFT Engine artifact
```

EOS, Equilibrium, and Regression are internal Engine modules. They use
repository-internal seams and one build graph rather than installed-package
dependencies. No production or validation repository depends on the lab.

Scientific packet flow is separate from production dependency direction:

```text
ePC-SAFT-lab -----------> ePC-SAFT-data   reviewed packet submission
ePC-SAFT-data ----------> ePC-SAFT-engine    pinned admitted runtime subset
ePC-SAFT-data ----------> ePC-SAFT-validation  pinned campaign input
ePC-SAFT-data ----------> ePC-SAFT-engine    pinned regression input
```

No installed package may import from or fetch Data at runtime. Engine wheels
must contain their admitted runtime subset. Validation and Engine regression
campaigns must materialize declared packets before execution and remain
replayable from an immutable Data commit. Wheels and SDK objects remain in
Governance's artifact store; Data does not replace it.

Cross-repository consumers MUST use the installed Engine artifact. They MUST
NOT use sibling source paths, compile Engine implementation sources, copy an
Engine kernel, or reach through private headers or Python modules. Downstream
applications MUST pin an immutable Engine wheel, not a Git source dependency.

The EOS native interface used by Equilibrium and Regression is internal to the
Engine implementation. Engine modules link repository-owned native targets in
one build and MUST NOT exchange installed headers or independently versioned
SDK wheels. Public C++ ABI stability remains out of scope unless a real
external C++ consumer is admitted.

Any earlier cross-repository SDK transport contract is transition evidence
only. It cannot govern the Engine implementation or justify retaining a
separate distribution.

## 5. Package shape

A production package starts from the smallest complete vertical slice. It does
not prebuild a speculative tree. Directories appear when accepted owners need
them. Placeholder modules, placeholder tests, empty native targets, `.gitkeep`
files, compatibility layers, fake defaults, and release workflows without a
release are forbidden.

The Engine uses one source repository and one Python distribution:

```text
.github/workflows/
src/epcsaft/
  eos/
  equilibrium/
  regression/
  _native/
tests/eos/
tests/equilibrium/
tests/regression/
tests/integration/
docs/
CMakeLists.txt
pyproject.toml
README.md
LICENSE
AGENTS.md
CONTEXT.md
```

Only needed entries are materialized. The `epcsaft` wheel may contain multiple
private native targets or extensions when solver dependencies require them,
but it has one version and one installation identity. Public Python imports use
`epcsaft`, `epcsaft.equilibrium`, and `epcsaft.regression`. The retired
`epcsaft-equilibrium`, `epcsaft-regression`, `epcsaft_equilibrium`, and
`epcsaft_regression` identities receive no compatibility aliases.

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

Engine owns compact module tests for public contracts, internal native seams,
equation contributions, derivative payloads, solver diagnostics, and a small
set of source-backed scientific anchors. Tests protect behavior, not old file
paths, source wording, private call graphs, or retired transition-era structure.

Engine tests MAY own test-only independent numerical or high-precision
oracles for white-box tensors. The validation repository owns test-only
physical and cross-package certificate oracles over installed immutable
artifacts. Neither oracle is runtime API or production input. Shared equations,
source, generated material, or copied kernels between an implementation and
its oracle MUST be declared and MUST NOT count as independent evidence.

The validation repository owns the durable black-box physical and cross-package
evidence. The accountable Engine task installs the declared wheel built from an
immutable commit and writes the bounded campaign there; an artifact need not be
publicly released. It MUST NOT import the Engine source tree or reproduce
production algorithms.

Task identity alone does not make evidence independent. An Engine-authored
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
units, domain, content hash, license or use basis, and tolerance rationale.
Reusable reviewed scientific bytes belong in Data. Each packet is
self-contained under `packets/<packet-id>/<version>/`, has `packet.toml`, and
may have `parameters/`, `observations/`, `sources/`, and `provenance/` content
when that content exists. Empty placeholder directories are forbidden.

Engine retains only a compact, admitted, self-contained runtime copy generated or
copied from a pinned Data release. Validation retains campaign definitions,
generated results, plots, and receipts that reference exact Data inputs; it
does not become the canonical owner of those inputs. Full papers, broad
extraction corpora, failed experiments, and exploratory notebooks belong in
the lab. Retained literature campaign plots belong in validation evidence, not
in every package test.

## 8. Promotion and unique authority

Lab-to-production promotion is distillation, not bulk transfer. Clean history
is the default. Original commits remain cited as provenance.

Every promoted slice requires a machine-validated receipt containing at least:

- capability, owning Engine module, and repository;
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

The first post-cutover runtime promotion is the smallest complete Engine
vertical slice. Later slices MUST reuse its public and internal native
contracts rather than add parallel paths.

## 9. Architecture ratchet

The Engine repository establishes its architecture baseline only
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

Engine and Governance repositories are public. They may remain unreleased
until a separate release decision exists. EOS, Equilibrium, and Regression
share one Engine version and one `epcsaft` wheel. No 1.0 or public release
workflow is implied by this doctrine.

Engine issues own production vertical slices. Data, Validation, Lab, and
Governance keep issues for work unique to those repositories. One Management
Project may aggregate them. Cross-repository references MUST use full URLs or
stable node identifiers, never ambiguous number-only references.

The public `.github` repository owns only the organization profile and
GitHub-recognized community defaults. Governance owns reusable workflow policy;
the Management coordination task executes project coordination. Engine and the
other active repositories still own their caller workflows, checks, build
matrices, and releases.

`main` MUST reject force pushes and deletion. Normal changes use review and
required evidence when those checks exist. A direct-main or administrator
bypass requires explicit user approval and a retained receipt; bypass is not a
standing exemption.

Remote creation, repository rename or transfer, issue transfer, Project
mutation, publication, and release remain explicit external-write gates. This
doctrine authorizes none of them by itself.

Before any repository rename, transfer, or old-name reuse, the Management
coordination task MUST perform an operation-specific read-only preflight,
capture immutable machine-readable exports for metadata at risk, inventory
Actions and reusable workflow references, and record redirect behavior,
namespace-reuse eligibility, metadata restoration, and every accepted loss. A
sacrificial rehearsal is required only when the preflight cannot establish the
behavior or rollback boundary. Each mutation requires a cutover record stating
preconditions, exact request, expected postconditions, retry and stop behavior,
reversibility, irreversible effects, and recovery owner, followed by a state
record proving the result.

## 11. Enforcement and amendment

Every repository `AGENTS.md` and `CONTEXT.md` MUST name the doctrine revision it
follows. Shared checks SHOULD verify the revision, forbidden dependency edges,
receipt schemas, and architecture ratchet. Reusable workflows still require a
small caller workflow in each active repository.

Amending this doctrine requires explicit user approval, a new revision, a
decision-ledger entry, and updates to affected repository pointers. Runtime
authority changes only through accepted promotion receipts.

Revision 8 records the user-approved direct cutover from three runtime
repositories and distributions to one `ePC-SAFT/ePC-SAFT-engine` repository,
one `epcsaft` distribution, and one Engine-rooted vertical-slice workflow.
Revision 7 remains the historical record of the Governance repository rename.
Historical plans and receipts retain their original names.

## Decision Ledger

| Decision | Source | Answer | Impact | Deferred? | Risk owner |
|---|---|---|---|---|---|
| D-001 Repository topology | User-approved topology review, superseded by D-017 | Historical six-home topology used EOS, equilibrium, regression, validation, organization policy, and personal lab | Preserved as historical topology | No | Management |
| D-002 Policy authority | Canonical-policy review, superseded by D-017 | Historical Organization `GOVERNANCE.md` was the permanent policy authority | Preserved as historical authority decision | No | Management |
| D-003 Package shape | Repository comparison and user objective, superseded by D-023 for runtime code | Historical package-first repositories separated EOS, Equilibrium, and Regression | Preserved as the pre-consolidation package shape | No | Engine and Management |
| D-004 Dependency direction | Ownership audit, superseded by D-023 for Engine-internal modules | Allow only one-way installed-artifact dependencies across repository seams; forbid sibling-source and Lab dependencies | Retains strict external seams while removing cross-wheel runtime coupling | No | Engine and consuming repository owners |
| D-005 EOS SDK | EOS SDK audit, superseded by D-023 for Engine-internal consumers | Keep the EOS native seam internal to Engine; require a separate decision only for a real external C++ consumer | Removes installed-header and cross-wheel ABI transfer from Equilibrium and Regression | Yes | Engine owner |
| D-006 Derivative authority | User-approved derivative doctrine | Make CppAD the sole nonlinear production derivative authority | Removes derivative modes and duplicate kernels | No | Engine owner |
| D-007 Implicit derivative order | Equilibrium Hessian sanity audit | Require every implicit and higher derivative order consumed by the solver | Prevents false exact-Hessian claims | No | Engine owner |
| D-008 Evidence boundary | Evidence-boundary audit | Combine Engine white-box receipts with Validation black-box artifact campaigns | Separates implementation proof from physical acceptance | No | Engine and Validation owners |
| D-009 Promotion authority | Historical Migration ownership audit | Use validated promotion receipts to create one production owner | Preserved historical admission rule; new coordination is Management-owned | No | Management |
| D-010 Architecture baseline | User optimization objective | Baseline Engine after the first complete slice and reject code-golf gaming | Keeps the architecture ratchet meaningful | No | Engine owner |
| D-011 Tracker ownership | Tracker design review, superseded by D-023 for runtime work | Use Engine vertical-slice issues and one Management Project; keep Data, Validation, Lab, and Governance issues in their owning repositories | Aligns work with the consolidated runtime owner | No | Management |
| D-012 Release timing | Release review, superseded by D-023 for runtime versioning | Keep repositories public while releasing EOS, Equilibrium, and Regression through one Engine version and artifact | Removes independent-version compatibility states | Yes | Engine release owner |
| D-013 Lab role | User lab-role decision, superseded by D-016 | The original monorepo served as transitional authority during the legacy cutover | Preserved as historical doctrine; no longer governs current lab or legacy repositories | No | Lab and Management |
| D-014 Pre-execution review remediation | Independent principal review and user-approved adjudication | Require whole-tree cutover inventory, portable hash-bound receipts, independent test-only oracle ownership, and rehearsed external operations | Preserved historical cutover control | No | Management |
| D-015 Package-rooted validation execution | User-approved workflow amendment, superseded by D-023 for task residence | Keep Validation as the evidence repository while the accountable Engine vertical-slice task performs its bounded installed-artifact campaign there; dispatch review only for a stable subject and keep no resident Validation worker | Preserves independent evidence with one runtime task owner | No | Engine, Validation-evidence, and Management owners |
| D-016 Lab and legacy separation | User-approved 2026-07-22 repository cutover | Give `ePC-SAFT/ePC-SAFT-lab` fresh history and a permanent non-production experimental role; retain the old monorepo and related projects under `ePC-SAFT-legacy` with no authority | Separates active research from historical code and removes fallback runtime authority from both lab and legacy | No | Management and Lab |
| D-017 Management control plane | User-approved 2026-07-29 topology decision | Rename the current authority repository to `ePC-SAFT-management`, keep `.github` thin, retire Migration, and centralize doctrine, manifests, ledgers, access, synchronization, and external-operation receipts in Management | Removes temporary coordination authority and duplicated project state | No | Management |
| D-018 Repository access | User-approved hierarchy, superseded by D-021 | Equilibrium and Regression may work in every operational repository; EOS may work only in EOS and Validation; Validation and Lab have no outbound access; every worktree is target-owned | Preserved as the pre-Data access policy | No | Management |
| D-019 Codex worktree environments | User-approved setup requirement | Track a worktree-relative Codex local environment in every active repository; setup must be idempotent and must not consume mutable sibling build output | Makes new Codex worktrees reproducible | No | Management and repository owners |
| D-020 Ephemeral worktrees and immutable artifacts | User-approved 2026-07-29 cleanup | Create target-owned worktrees only on demand and retire them after merge or recorded preservation; install cross-repository wheels only through the SHA-256 object store and reject mutable repository `dist/` wheels | Prevents stale checkout and SDK accumulation | No | Management and repository owners |
| D-021 Canonical scientific Data home | User-approved 2026-07-31 repository and ownership amendment | Create public `ePC-SAFT/ePC-SAFT-data` for reusable self-contained scientific packets; separate Data storage ownership from Engine runtime-subset ownership and Validation evidence ownership; allow read-only packet flow to Engine, Validation, and Lab while forbidding runtime repository dependencies | Removes ambiguous parameter ownership and keeps canonical inputs distinct from runtime copies, generated evidence, artifacts, and exploratory provenance | No | Management, Data, and consuming repository owners |
| D-022 Governance repository identity | User-approved 2026-07-31 repository rename | Use `ePC-SAFT/ePC-SAFT-governance` as the durable Git-tracked authority repository; keep Management as the non-Git coordination task in the `ePC-SAFT-project` parent; preserve former names in historical plans and receipts | Separates current repository authority from coordination context without rewriting accepted evidence | No | Governance and Management |
| D-023 Unified Engine direct cutover | User-approved 2026-08-04 package and agent-friction review | Rename the EOS repository to `ePC-SAFT-engine`; consolidate EOS, Equilibrium, and Regression into one `epcsaft` distribution and namespace; assign each production capability to one Engine-rooted vertical-slice task; retain Data and Validation as separate repository seams; cut downstream applications to one immutable Engine wheel without compatibility aliases or a mixed active release | Removes cross-wheel SDK drift, independent version mismatches, Git-source downstream pins, and multi-repository package edits while preserving independent scientific inputs and black-box evidence | No | Engine, Management, and downstream application owners |

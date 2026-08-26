# ePC-SAFT Ecosystem Doctrine

Status: canonical

Doctrine revision: 24

Approved: 2026-08-25

Owner: `tannerpolley/ePC-SAFT:governance/`

The private `tannerpolley/ePC-SAFT` monorepo is the sole active Git and Codex
home. `governance/` contains doctrine and control records; Management is the
coordination role performed from the repository root, not a second authority.

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

## 3. Directory ownership

| Home | Sole ecosystem role | Must not own |
|---|---|---|
| `governance/` | Doctrine, topology, policy, manifests, ledgers, readiness checks, and external-operation receipts | Package source, package tests, scientific capabilities, or releases |
| `data/` | Canonical reusable scientific packets: parameter bundles, source-derived datasets, normalized observations, provenance, units, domains, hashes, and license or use bases | Runtime algorithms, generated validation results, package artifacts, full papers, or exploratory and rejected material |
| `engine/` | The single `epcsaft` runtime distribution: EOS/state/property evaluation, CppAD derivatives, equilibrium formulations and certification, regression targets and diagnostics, one native build, and one release identity | Canonical scientific packets, durable black-box validation evidence, application workflows, or separate Equilibrium and Regression distributions |
| `validation/` | Installed-artifact, cross-module, downstream-consumer, and literature acceptance evidence authored by the accountable Engine task | Production algorithms, private Engine imports, a paper archive, or a permanent execution owner |

D-026 supersedes the active polyrepo topology and D-019's environment-per-repo
rule. The source repositories remain immutable migration provenance after
archive. The separate `ePC-SAFT/.github` repository remains a thin organization
profile and is not part of the active project Git graph.

One concept has one runtime owner. Canonical scientific storage does not confer
runtime authority. A capability without an accepted production
owner is non-production; no legacy copy supplies transitional runtime
authority. Promotion distills a clean, independently evidenced slice into its
unique package owner rather than transferring authority from a runnable
archive.

D-029 establishes one narrow nonreactive-phase transition inside Engine. The
exact `tp_flash` HELD/HELD2 calculation is the protected incumbent while the
HELD-derived `phase_equilibrium` interface gains explicit phase-selection
semantics and is independently compared against the incumbent. This is one
numerical core with an incumbent and a candidate interface, not two physical
owners, an alias, or an automatic fallback. Callers select a route explicitly;
neither route may silently substitute the other. An unrelated replacement
solver cannot claim HELD lineage. `tp_flash` and its HELD/HELD2 evidence may be
retired only after a reviewed, user-approved promotion receipt proves that the
candidate preserves every admitted incumbent workflow within its declared
tolerances and satisfies its additional constrained-phase claims.

D-030 completes that transition. `phase_equilibrium` is the sole public
nonreactive phase-equilibrium route. A fixed-T,P feed without a declared
topology requests the retained HELD/HELD2 discovery algorithm; a declared
one-to-three-phase topology requests only that constrained solve and may not
fall back to discovery. Phase compositions and densities are non-authoritative
hints for declared phases, never implicit topology. The HELD/HELD2 numerical
core remains the internal discovery and start-generation owner; `tp_flash` and
its dedicated Python result, error, state, and discovery-marker surface are
retired. The accepted CAP-09 promotion receipt records the user-approved
`1e-4` historical lineage comparison while preserving the stricter physical,
certificate, and derivative gates.

D-031 keeps the Regression cutover boundary EOS-only after the CAP-10
Equilibrium-observation admission review. Omitted-topology discovery is not an
ordinary differentiable Observation Family because its phase count and output
rows may change and its finite search does not establish globality. Declared
fixed-T,P splits and fixed-T coexistence boundaries expose exact local
sensitivities only for certified, fixed-topology, nonsingular branches; their
current interfaces do not provide the owner-published batch descriptor,
optimizer-domain branch policy, nonuniqueness policy, or family-level
Non-Evaluable Trial semantics required by Regression. No dormant Equilibrium
Observation Family schema, enum, resolver branch, or evaluator is admitted
into Regression. The live `phase_equilibrium` problem/result interface and
internal HELD/HELD2 discovery and start-generation implementation remain
Equilibrium-owned. A later family requires a new implementation issue and
independent numerical evidence rather than silently expanding Regression.

D-034 admits the fixed-T,P homogeneous reactive Observation Family implemented
by issue #30. Equilibrium owns its closed descriptor, exact total derivatives,
source-reference transfer, certified local-state and immutable continuation
receipts, physical branch-family identity, and typed Non-Evaluable Trial
outcomes. Regression consumes those owner values and Jacobians through the one
`ParameterFit`/`FitResult` route and owns residuals, weights, bounds, affine
coordinates, starts, and optimizer accounting. The admission is local and
source-neutral: it makes no global-equilibrium claim, does not admit
nonreactive phase-equilibrium or coupled reactive-VLE observations, and does
not authorize parameter promotion from diagnostic or manufactured evidence.
The accepted receipt is
`governance/receipts/cap-12-homogeneous-reactive-observation-admission.json`.

D-035 admits the fixed-topology reactive bubble/VLE Observation Family
implemented by issue #31. Equilibrium owns the coupled liquid-reaction,
liquid-state, vapor-state, and bubble-pressure value problem; exact total
Schur-complement derivatives; ordered neutral vapor support; immutable local
branch and exact-trial receipts; certification; and typed Non-Evaluable Trial
outcomes. The primary vapor path is a nonideal Provider EOS. An ideal-gas vapor
is available only as an explicitly typed comparison model. Regression consumes
the owner values and exact Jacobians through the existing
`ParameterFit`/`FitResult` route. Admission is source-neutral and local: it
makes no global-equilibrium or parameter-promotion claim, and it fails closed
on unresolved branches, disappearing liquid support, no incipient vapor,
Provider-domain rejection, or singular/ill-conditioned coupled sensitivity.
The accepted receipt is
`governance/receipts/cap-13-reactive-bubble-vle-admission.json`.

D-036 approves General Reactive Equilibrium and Phase Equilibrium (GREPE) as
the target architecture that evolves the sole public `phase_equilibrium`
owner. Typed fixed/solved state roles, finite/incipient phase roles, declared
species support, and declared reactions determine the compiled physical
blocks; application, source, component, and operation names do not select a
solver. Finite reacting phases use one constrained Gibbs/KKT encoding so
reaction and interphase equilibrium are not duplicated; admissible amount
changes are restricted to declared reaction extents plus species-preserving
interphase transfer. Incipient phases use their normalization and shared-
support chemical-potential/fugacity equalities without an incipient material
balance. Each thermodynamic condition is solved and locally certified by Ipopt
with exact objective/constraint derivatives and exact Lagrangian Hessians;
outer Regression/Ceres consumes exact total observable Jacobians and remains
an authority-neutral parameter-estimation tool.

GREPE accepts only phase-aware, unit- and convention-complete reaction standard
states after rank, conservation, charge, and consistency preflight. Declared
phase solves do not run phase discovery or claim globality. Regression owns no
thermodynamic family switch and may explicitly export a provenance-linked
candidate parameter document, but it cannot accept, admit, promote, or adopt
parameters. The canonical target plan is
`engine/docs/equilibrium/plans/2026-08-14-general-reactive-equilibrium-phase-equilibrium.md`.
This architecture decision admits no new runtime capability: D-034 and D-035
remain the current reactive Observation Family boundary until later slices
earn independent numerical evidence and accepted receipts.

D-037 admits the reviewed GREPE proven core implemented across issues #47,
#48, and #49. The sole `phase_equilibrium` owner may now compile declared
finite and incipient phase supports, phase-scoped electrolyte solvents,
temperature-dependent permittivity, declared reactions, affine feeds, and
imposed intensive boundaries through its existing homogeneous or reactive-
bubble owner. A charged phase may contain multiple neutral components, but its
electrolyte term still requires exactly one declared dielectric solvent within
that phase. The implementation remains source- and component-neutral and does
not discover topology, choose phases by TPD, or claim global equilibrium.

The admitted intensive-boundary route retains exact owner derivatives and
typed non-evaluable trials. Its outer receipt binds the solved feed, reaction
system, and immutable owner topology to the validated wrapped receipt. For the
reactive-bubble wrapper it also cross-links phase compositions, pressures, and
observation values. Homogeneous values and derivatives remain executable
falsifier evidence rather than being recomputed by receipt validation. The
fixed-pressure incipient parity case reproduces
the retired specialized route's pressure, solved feed coordinate, phase
compositions, three observables, and their parameter derivatives; independent
perturbation solves also verify those derivatives. No other uncommitted issue-
#48 application, reference-correction, saturation, or campaign layer is
admitted. The accepted receipt is
`governance/receipts/cap-14-grepe-proven-core-admission.json`.

D-032 completes the general scientific capability hard cutover. The accepted
CAP-11 subject commit removes the superseded live permittivity design, the
executable Cython HELD2 experiment, retired split-header packaging entries,
literal-component private HELD coordinates, and dead legacy switches while
retaining the HELD/HELD2 numerical algorithms selected by charge topology and
conditioning. One exact-commit retirement scan records every match and
justified provenance exception; one installed wheel passes the complete root,
native, sanitizer, and installed-wheel Validation matrix. The receipt at
`governance/receipts/cap-11-hard-cutover-retirement.json` is the final program
acceptance record. Future scientific campaigns, including issue #6, are
separate post-cutover Validation work rather than CAP blockers. Issues #30 and
#31 are likewise detached post-cutover reactive-observation follow-ons, not
children or acceptance gates of root issue #10.

D-033 makes low-latency private development an explicit project invariant.
Environment setup builds and installs the sole non-editable Engine wheel before
any test process starts. The default affected and all-owner developer gates use
that installed wheel and target less than ten seconds on a warm checkout by
running independent compact numerical and readiness checks concurrently.
Slow/expensive Engine stress evidence, installed-wheel Validation campaigns or
figure reproduction, full native clang-tidy, and sanitizer builds run only by
an explicit accountable-owner command when the changed scientific claim or
native risk requires them. Documentation or tool changes do not silently
escalate into those evidence workflows. This project is private and under
active development for downstream applications; broad publication-grade replay
must not obstruct wheel iteration.

The complete pre-migration monorepo and related historical repositories live
outside the authority topology under the local `ePC-SAFT-legacy` container.
The public `tannerpolley/ePC-SAFT-legacy-monorepo` is a read-only historical
repository. Legacy code, tracker state, branches, tags, and releases are
provenance only; they own no current runtime, roadmap, or governance authority.

Directory ownership and task residence are distinct. Validation is the durable
evidence home, not a standing agent or implementation team. One Engine
vertical-slice task remains accountable for each capability and may change the
necessary Data, Engine, and Validation files in one branch. Management prevents
concurrent writers; it does not maintain a waiting Validation worker.

The completed `tannerpolley/ePC-SAFT-migration` repository is historical
evidence outside the active topology. It owns no current plan, queue, gate,
decision, or authority.

Data is a read-only scientific-input home, not a service or runtime dependency.
Engine may copy or generate an admitted runtime-catalog subset from a pinned
Data release. Validation and Engine regression campaigns may consume packets
from exact Data commits for campaigns and declared fitting inputs. The
consuming module or campaign remains self-contained and records the Data
commit, packet path, and packet fingerprint.

## 4. Dependency direction

The only production dependency edge is:

```text
validation/ -----> installed wheel built from engine/
```

EOS, Equilibrium, and Regression are internal Engine modules. They use
repository-internal seams and one build graph rather than installed-package
dependencies.

Scientific packet flow is separate from production dependency direction:

```text
data/ ----------> engine/       pinned admitted runtime subset
data/ ----------> validation/   pinned campaign input
data/ ----------> engine/       pinned regression input
```

No installed package may import from or fetch Data at runtime. Engine wheels
must contain their admitted runtime subset. Validation and Engine regression
campaigns must materialize declared packets before execution and remain
replayable from an immutable Data commit. Wheels and SDK objects remain in
Governance's artifact store; Data does not replace it.

Validation and downstream consumers MUST use the installed Engine artifact.
They MUST NOT import the `engine/` source tree, compile Engine implementation
sources, copy an Engine kernel, or reach through private headers or Python
modules. Downstream applications MUST pin an immutable Engine wheel, not a Git
source dependency.

A named downstream application MAY be quarantined during an Engine authority
cutover only by an explicit user decision recorded in Governance. A quarantined
application is not an admitted Engine consumer, carries no current ecosystem
scientific or production claim, and MUST NOT be executed or published as a
supported workflow until its owner installs an accepted Engine wheel and passes
the applicable Validation evidence. Its retained dependency state is deferred
work, not a compatibility path or alternate runtime authority.

D-038 supersedes D-024's quarantine for `tannerpolley/MEA-Thermodynamics`,
`tannerpolley/MEA-Absorption-Column`, and `tannerpolley/Lithium_Extraction`.
They are official downstream applications and MAY be co-developed by the
accountable Engine task when the requested slice crosses the application seam.
They remain separate sibling Git repositories. Official status grants work and
integration authority; it does not admit an unfinished model, fitted parameter
set, validation result, flowsheet, or manuscript claim. Each application MUST
consume an immutable Engine wheel with recorded artifact identity and MUST keep
application chemistry, data, analyses, and manuscripts in its own repository.

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
src/epcsaft/
  eos/
  equilibrium/
  regression/
  _native/
tests/eos/
tests/equilibrium/
tests/regression/
docs/
meson.build
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

## 7. Numerical Verification-Driven Development and evidence ownership

The ecosystem uses Numerical Verification-Driven Development (NVDD). Tests are
numerical evidence for the current scientific formulation, not permanent
software contracts. A test MAY exist only when it can falsify at least one of:

- an independently known numerical value;
- a mathematical identity, conservation law, symmetry, or limiting value;
- a derivative against an independent derivative route;
- a nonlinear-system, KKT, stability, or equilibrium residual certificate;
- a manufactured, published, or otherwise independently known solver result;
- a metamorphic numerical relation; or
- an experimental or literature validation metric.

Every retained test MUST identify the scientific claim, oracle class, compared
quantities, units or scaling, explicit absolute and relative tolerances, the
basis for those tolerances, and the scientific meaning of failure. Expected
values produced by the implementation under test are
`change_detection_only`; they are not independent verification.

Tests MUST NOT exist only to protect a module surface, compatibility,
deprecated or legacy behavior, a private helper, call order, class layout,
exception prose, serialization shape, solver label, iteration count, log
output, diagnostic ordering, snapshot, branch, line, or coverage percentage.
Mocks and call-count assertions are forbidden in scientific-core tests. Exact
phase ordering and exact nonlinear iteration counts are forbidden acceptance
criteria. Invalid-input tests are retained only when the invalid condition
could otherwise produce a plausible but scientifically incorrect result.

Tests follow the current intended formulation and MAY be revised or deleted
when that formulation changes. Only benchmark evidence bound to a publication,
tagged scientific release, or explicit reproducibility artifact is immutable.
Line coverage is informational; numerical coverage of equations, derivatives,
properties, equilibria, and validation domains is the quality measure.

Engine owns compact numerical tests for equation contributions, derivative
payloads, solver certificates, mathematical limits, manufactured solutions,
and a small set of source-backed scientific anchors. It does not own API
compatibility or architecture contract suites.

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
does not become the canonical owner of those inputs. Full papers and broad
extraction corpora remain outside the governed source tree. Failed experiments
and exploratory notebooks are not retained as active project areas. Retained
literature campaign plots belong in validation evidence, not in every package
test.

## 8. Promotion and unique authority

Production promotion is distillation, not bulk transfer. Clean history is the
default. Original commits remain cited as provenance.

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
destination is the sole production owner and any legacy copy remains
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

All active repositories are private. This is a rapidly evolving research
project, not a public compatibility surface. EOS, Equilibrium, and Regression
share one Engine version and one `epcsaft` wheel. Existing public release and
repository events remain historical facts; no new publication or release is
implied by this doctrine.

Engine issues own production vertical slices. Data, Validation, and Governance
keep issues for work unique to those directories. One Management Project may
aggregate them. Cross-repository references MUST use full URLs or stable node
identifiers, never ambiguous number-only references.

The private `.github` repository owns only organization metadata and links to
Governance. The Management coordination task executes project coordination.
No active repository owns CI/CD, caller workflows, reusable workflows, build
matrices, deployment automation, or release automation. Reproducible local
commands belong in each repository's Codex environment and repository scripts.

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
follows. Repository-local preflight commands SHOULD verify the revision,
forbidden dependency edges, receipt schemas, architecture ratchet, and NVDD
evidence classification. No CI/CD or caller workflow is required.

Amending this doctrine requires explicit user approval, a new revision, a
decision-ledger entry, and updates to affected repository pointers. Runtime
authority changes only through accepted promotion receipts.

Revision 24 records the user-approved removal of the active Lab area. D-039
deletes `lab/`, removes it from current ownership and routing, and preserves
only historical or scientific provenance that identifies the former Lab or its
source material.

Revision 23 records the user-approved official downstream application model.
D-038 removes the D-024 quarantine from MEA-Thermodynamics,
MEA-Absorption-Column, and Lithium_Extraction; keeps all three as separate
sibling Git repositories; authorizes one accountable Engine task to cross the
named application seam when required; and retains the immutable-wheel boundary
and application-owned scientific evidence. Official status is not scientific
admission of any particular result.

Revision 22 records the user-approved proven-core integration. D-037 admits
only the generic, reviewed GREPE runtime shared by issues #47, #48, and #49,
including exact imposed-boundary derivatives and parity with the retired
fixed-pressure incipient route.

Revision 21 records the user-approved GREPE target architecture. D-036 evolves
the one public `phase_equilibrium` owner toward declarative homogeneous
reaction, reactive bubble/dew/VLE, and declared-phase reactive LLE work while
keeping phase discovery separate. It fixes the inner-Ipopt/outer-Ceres
boundary, exact Hessian and total-Jacobian obligations, phase-aware reaction
standard states, local-only certification, generic Regression evaluation, and
authority-neutral candidate parameter export. It admits no unimplemented
reactive LLE, electrochemical multiphase, or parameter-promotion capability.

Revision 20 records the user-approved issue-#31 admission. D-035 admits the
fixed-topology reactive bubble/VLE Observation Family, its nonideal Provider
vapor path, exact total Schur Jacobian, local branch/support certification, and
canonical Regression integration. It does not admit nonreactive
phase-equilibrium observations or any global-equilibrium or parameter-promotion
claim.

Revision 19 records the user-approved issue-#30 admission. D-034 admits only
the fixed-T,P homogeneous reactive Observation Family, its exact total
Jacobian and typed local-branch contracts, and its canonical Regression fit
integration. Coupled reactive-VLE and nonreactive phase-equilibrium Observation
Families remain outside Regression.

Revision 18 records the user-approved build-first, latency-bounded private
development workflow. D-033 separates the sub-ten-second warm developer gate
from explicit owner-run stress, campaign, native-quality, and sanitizer
evidence without deleting any numerical falsifier.

Revision 17 records the user-approved final CAP execution sequence. D-032
accepts the exact CAP-11 cutover tree, bounded retirement scan, one-wheel
verification, and cleanup evidence and closes the general scientific
capability hard-cutover program without absorbing later Validation campaigns.

Revision 16 records the user-approved CAP-10 admission decision. D-031 defers
all current mixture phase-equilibrium Observation Families, keeps Regression
EOS-only, and defines the evidence required before a later family may be
admitted.

Revision 15 records the user-approved completion of the protected HELD/HELD2
promotion. D-030 makes `phase_equilibrium` the sole public nonreactive owner,
uses omitted topology for internal HELD/HELD2 discovery, requires explicit
topology for constrained solves, and retires only the redundant Python flash
surface after the accepted promotion receipt.

Revision 14 records the user-approved recovery and protected promotion of the
HELD/HELD2 nonreactive phase workflow. D-029 restores `tp_flash` as the exact
incumbent, requires `phase_equilibrium` to evolve from that numerical lineage,
and forbids incumbent retirement before reviewed as-good-or-better evidence.

Revision 13 records the user-approved pre-CAP execution hygiene amendment.
D-028 removes the retired integration-test shape, makes the native CAP issue
graph the execution authority, and assigns each CAP leaf to one fresh,
single-writer task with proportionate read-only review.

Revision 12 records the user-approved conversion of Lab into a zero-test
playground. D-027 removes permanent Lab oracle frameworks and routes every
durable numerical claim to Engine or installed-wheel Validation evidence.

Revision 11 records the user-approved direct cutover to one private monorepo,
one Git root, one Codex environment, one issue tracker, and one end-to-end
vertical-slice agent. D-026 supersedes active repository topology and
environment-per-repository rules while retaining directory ownership and all
historical evidence.

Revision 10 records the user-approved NVDD cutover, removal of CI/CD as an
ecosystem responsibility, and conversion of all active repositories to private
visibility. D-025 supersedes D-012's and D-021's public-visibility choices
without rewriting their historical records.

Revision 9 records the user-approved package-authority cutoff while Lithium
Extraction, MEA Thermodynamics, and MEA Absorption Column are quarantined for a
later Engine migration. D-024 changes only the authority gate ordering: it does
not admit their old Git pins, results, or workflows and does not weaken the
single-Engine architecture approved by D-023.

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
| D-019 Codex worktree environments | User-approved setup requirement, superseded by D-026 | Track a worktree-relative Codex local environment in every active repository; setup must be idempotent and must not consume mutable sibling build output | Preserved as the polyrepo environment rule; D-026 supplies one worktree-relative root environment | No | Management and repository owners |
| D-020 Ephemeral worktrees and immutable artifacts | User-approved 2026-07-29 cleanup | Create target-owned worktrees only on demand and retire them after merge or recorded preservation; install cross-repository wheels only through the SHA-256 object store and reject mutable repository `dist/` wheels | Prevents stale checkout and SDK accumulation | No | Management and repository owners |
| D-021 Canonical scientific Data home | User-approved 2026-07-31 repository and ownership amendment | Create public `ePC-SAFT/ePC-SAFT-data` for reusable self-contained scientific packets; separate Data storage ownership from Engine runtime-subset ownership and Validation evidence ownership; allow read-only packet flow to Engine, Validation, and Lab while forbidding runtime repository dependencies | Removes ambiguous parameter ownership and keeps canonical inputs distinct from runtime copies, generated evidence, artifacts, and exploratory provenance | No | Management, Data, and consuming repository owners |
| D-022 Governance repository identity | User-approved 2026-07-31 repository rename | Use `ePC-SAFT/ePC-SAFT-governance` as the durable Git-tracked authority repository; keep Management as the non-Git coordination task in the `ePC-SAFT-project` parent; preserve former names in historical plans and receipts | Separates current repository authority from coordination context without rewriting accepted evidence | No | Governance and Management |
| D-023 Unified Engine direct cutover | User-approved 2026-08-04 package and agent-friction review | Rename the EOS repository to `ePC-SAFT-engine`; consolidate EOS, Equilibrium, and Regression into one `epcsaft` distribution and namespace; assign each production capability to one Engine-rooted vertical-slice task; retain Data and Validation as separate repository seams; cut downstream applications to one immutable Engine wheel without compatibility aliases or a mixed active release | Removes cross-wheel SDK drift, independent version mismatches, Git-source downstream pins, and multi-repository package edits while preserving independent scientific inputs and black-box evidence | No | Engine, Management, and downstream application owners |
| D-024 Package authority with quarantined consumers | Explicit user instruction on 2026-08-04 to begin the full ePC-SAFT package cutoff and fix Lithium and MEA later | Complete the one-repository, one-wheel Engine authority cutover now; quarantine Lithium Extraction, MEA Thermodynamics, and MEA Absorption Column as non-admitted consumers until separate vertical-slice migrations pass | Removes the three-package authority immediately without treating deferred Git-pinned applications or failed scientific routes as accepted Engine evidence | Yes: application migration | Engine, Management, and downstream application owners |
| D-025 Private NVDD development | Explicit user instruction on 2026-08-05 adopting the revised Numerical Verification-Driven Development policy project-wide, authorizing CI/CD removal, and authorizing private repository visibility | Treat tests only as current numerical or scientific evidence; prohibit compatibility, API-shape, mock, snapshot, exact-prose, solver-label, and coverage-only suites; keep validation breadth in explicit local commands; remove CI/CD ownership; make every active repository private | Reduces agent friction and permanent contract surface while preserving numerical correctness, source traceability, publication locks, and independent Validation/Lab evidence | No | Governance, Engine, Validation, Data, Lab, and Management owners |
| D-026 Single private monorepo | Explicit user instruction on 2026-08-05 after package-management and agent-friction review | Move Governance, Engine, Data, Validation, and Lab main histories into private `tannerpolley/ePC-SAFT`; use one Git root, root Codex environment, issue tracker, and Engine-accountable vertical-slice branch while retaining strict directory ownership; archive rather than delete source repositories | Removes cross-repository branch, worktree, wheel, and version coordination while preserving scientific independence and recoverable history | No | Root Management role and directory owners |
| D-027 Zero-test Lab playground | Explicit user instruction on 2026-08-05 after reviewing all 51 Lab oracle checks | Make Lab a temporary playground with no persistent tests, test framework, or root verification gate; delete self-testing manufactured oracle implementations; distill any future durable numerical claim into Engine or installed-wheel Validation evidence | Removes a parallel 3,945-line toy implementation and prevents Lab experiments from becoming permanent maintenance contracts | No | Lab, Engine, Validation, and Management owners |
| D-028 CAP execution hygiene | Explicit user instruction on 2026-08-07 to fix live pre-CAP guidance and establish one model, task, order, and review strategy | Remove `tests/integration/` from the normative package shape; use GitHub root issue #10 and native blockers for lifecycle state; execute each CAP leaf in a fresh single-writer task using the canonical Engine task map and proportionate read-only review | Prevents retired test structure, stale polyrepo routing, context-overloaded implementation, and per-task execution-setting drift from re-entering the cutover | No | Engine and Management |
| D-029 Protected HELD promotion | Explicit user correction on 2026-08-10 that CAP-09 was intended to rename and extend the HELD/HELD2 workflow, not replace it with an unrelated solver | Restore exact `tp_flash` HELD/HELD2 as the protected incumbent; require `phase_equilibrium` to share that numerical lineage while adding explicit phase selection; prohibit automatic fallback and retire the incumbent only after reviewed, user-approved as-good-or-better promotion evidence | Recovers the validated HELD2 workflow, makes the mistaken scratch replacement non-authoritative, and prevents another destructive phase-owner cutover before equivalence and new constrained claims are demonstrated | Yes: full constrained one-to-three-phase promotion | Engine and Management |
| D-030 Sole phase-equilibrium promotion | Explicit user approval on 2026-08-12 of a `1e-4` lineage-only tolerance, one public `phase_equilibrium` route, omitted-topology HELD/HELD2 discovery, and explicit declared phase count/types | Accept the CAP-09 promotion receipt; retain HELD/HELD2 internally; remove public `tp_flash` and its dedicated Python types; interpret absent fixed-T,P topology as discovery and declared topology as a no-fallback constrained solve | Completes D-029 without discarding the validated HELD2 algorithm, removes the transitional duplicate API and source-specific Python gate, and keeps physical/certificate tolerances strict | No | Engine and Management |
| D-031 Defer Equilibrium Observation Families | Explicit user instruction on 2026-08-12 to apply the reviewed CAP-10 sequence after the general phase owner and exact sensitivities were complete | Admit no current mixture phase-equilibrium Observation Family; keep Regression EOS-only; require a new implementation issue for any future Equilibrium-owned batch descriptor and numerical admission evidence | Prevents variable-topology discovery, local branch sensitivities, and incomplete optimizer-domain semantics from becoming a misleading differentiable fit interface while preserving the exact Equilibrium sensitivity seam for its current owner | Yes: future Equilibrium observation admission | Engine and Management |
| D-032 Complete the scientific capability hard cutover | Explicit user instruction on 2026-08-12 to apply the recommended CAP sequence through final closure | Accept CAP-11 subject `a8efa571d52c32f5a697f15dd7f4623b7e3613fd`, its bounded retirement scan, one-wheel verification, independent reviews, and final receipt; close CAP root #10 while leaving issues #6, #30, and #31 as detached post-cutover work | Removes the last superseded live design, executable Lab duplicate, split-header remnants, literal-component internal names, author-named live diagnostics, and dead legacy switches without deleting canonical HELD/HELD2 discovery or historical provenance | No | Engine and Management |
| D-033 Build first and keep developer verification fast | Explicit user instruction on 2026-08-12 that the private in-development Engine must support rapid downstream iteration and that tests must not block wheel build/install | Build and install the sole wheel before tests; make `check-affected.sh` and `check-all.sh` latency-bounded installed-wheel gates with a warm target below ten seconds; keep slow Engine stress evidence, Validation campaigns/figures, full native analysis, and sanitizers behind explicit owner commands | Reduces the measured warm gate from 821.09 seconds to 8.96 seconds while retaining every slow scientific falsifier in its owner tier | No | Engine, Validation, and Management |
| D-034 Admit homogeneous reactive observations | Explicit user instruction through issue #30 and the reviewed CAP-12 correction/admission gates on 2026-08-13 | Admit only the fixed-T,P homogeneous reactive Observation Family and its canonical Regression integration; require immutable role-certified continuation anchors, exact total Jacobians, closed typed failures, exact trial/branch receipts, and local/no-global semantics; retain phase-equilibrium and coupled reactive-VLE families as unadmitted | Enables the first source-neutral reactive regression tracer without weakening D-031's prohibition on variable-topology or insufficiently certified phase-equilibrium fitting | No | Engine, Validation, Regression, and Management |
| D-035 Admit reactive bubble/VLE observations | Explicit user instruction through issue #31 and the reviewed CAP-13 correction/admission gates on 2026-08-13 | Admit only the fixed-topology reactive bubble/VLE Observation Family and its canonical Regression integration; require exact same-liquid-state binding, ordered neutral vapor support, nonideal Provider vapor as the primary path, exact total Schur derivatives, immutable local branch receipts, and typed fail-closed outcomes | Enables coupled pressure/speciation regression without weakening D-031's prohibition on unqualified variable-topology phase-equilibrium fitting or making a global-equilibrium claim | No | Engine, Validation, Regression, and Management |
| D-036 General Reactive Equilibrium and Phase Equilibrium target | User-confirmed GREPE design and two independent read-only sanity reviews on 2026-08-14 | Evolve the sole public `phase_equilibrium` owner through declarative finite/incipient reacting-phase blocks; use one nonduplicated Gibbs/KKT encoding, phase-aware reaction standards, exact Ipopt Hessians, exact outer-Ceres Jacobians, local declared-phase certificates, generic Regression evaluation, and authority-neutral candidate export | Establishes one staged architecture for homogeneous reaction, reactive bubble/dew/VLE, neutral reactive LLE, and later electrochemical LLE without admitting unimplemented capability or adding an application-specific route | Architecture target only; each runtime slice requires its own accepted receipt | Engine, Validation, Regression, Data, downstream applications, and Management |
| D-037 Admit the GREPE proven core | Explicit user instruction on 2026-08-18 to combine the reviewed generic Engine work, retain only the proven core, reproduce the removed specialized path's results, and use an exact receipt before updating main | Admit declared finite/incipient supports, phase-scoped dielectric solvent selection, temperature-dependent permittivity, generic affine-feed intensive boundaries over the existing homogeneous/reactive-bubble owners, exact total derivatives, typed failures, and hash-bound nested receipts at subject commit `25cd674bb32611bac89d0ffb4ab8e1d953cbe09f` | Replaces a specialized fixed-pressure incipient implementation with the smaller generic owner while preserving its numerical results; excludes the other issue-48 dirty application/reference/saturation layer | No; later globality, phase discovery, and application parameter promotion remain separate | Engine, Validation, Regression, downstream applications, and Management |
| D-038 Official downstream applications | Explicit user instruction on 2026-08-20 to maintain MEA-Thermodynamics, MEA-Absorption-Column, and Lithium_Extraction as official Engine applications while preserving Engine/application ownership | Remove the D-024 quarantine for the three named repositories; keep them as separate sibling Git repositories; authorize cross-repository vertical-slice work; require immutable Engine-wheel identity and separate repository histories | Makes cross-repository work durable without nested repositories, submodules, sibling-source imports, or automatic scientific admission | No | Engine, downstream application owners, and Management |
| D-039 Remove the active Lab area | Explicit user instruction on 2026-08-25 to delete the Lab subtree and all live routing while retaining required provenance | Remove `lab/` from the governed topology, ownership manifest, checks, and live instructions; retain former Lab names and paths only in historical decisions, receipts, source identities, and scientific provenance | Removes a non-production project area without creating a replacement scratch framework or weakening Engine, Data, or Validation ownership | No | Governance and Management |

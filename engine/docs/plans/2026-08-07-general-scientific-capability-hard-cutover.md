# General Scientific Capability Hard-Cutover Plan

- Status: canonical plan; implementation not started
- Date: 2026-08-07
- Doctrine: Governance revision 12, Numerical Verification-Driven Development
- Scope: `engine/` EOS, Equilibrium, and Regression
- Decision authority: ADR 0004 and the three owner plans linked below

## Outcome

Replace source-, campaign-, component-, and fixed-topology runtime routes with
one source-neutral capability system. `Parameters` resolves parameter identity,
EOS evaluates thermodynamic quantities and exact partial derivatives,
Equilibrium owns solved-state quantities and sensitivities, and Regression forms
residuals from owner-provided observable blocks. The cutover is direct: no old
symbol, schema, fallback, inferred default, or compatibility implementation
survives.

The owner plans are:

- [EOS generic parameter derivatives](../eos/plans/2026-08-07-general-eos-parameter-derivatives.md)
- [Regression owner-driven parameter fit](../regression/plans/2026-08-07-owner-driven-parameter-fit.md)
- [Equilibrium general phase equilibrium](../equilibrium/plans/2026-08-07-general-phase-equilibrium.md)

The pair-parameter scientific basis is recorded in
[PC-SAFT/ePC-SAFT pair-parameter symmetry and meaning](../research/2026-08-07-pair-interaction-parameter-symmetry.md).

## Non-negotiable decisions

1. One `Parameters`-owned `ActiveParameterSet` crosses all three owners without
   translation into module-specific parameter enums or slots.
2. A production route may be specialized by physical formulation or numerical
   algorithm, but never by author, paper, campaign, application, fingerprint,
   literal component name, fixed component position, or special species tuple.
3. EOS and Equilibrium advertise typed observable capabilities. Regression
   consumes their values and Jacobians; it does not implement thermodynamics or
   phase equilibrium.
4. The initial Regression boundary remains EOS fixed-state and pure-saturation
   observations. Equilibrium observations require a later numerical-failure
   review after the general Equilibrium owner is complete.
5. The sole solvation formulation is `born`. Its explicit correction
   coefficients reduce mathematically to the original Born expression at their
   neutral zero/one values. There is no `original-born` selector.
6. The nonreactive phase owner is `phase_equilibrium`, with discovery, hint, and
   constrained-specification modes for one to three liquid/vapor phases.
7. Homogeneous chemical equilibrium remains separate. This plan defines only
   the state, phase, parameter, and sensitivity seam needed by a future coupled
   reactive-phase design.
8. Generated parameter templates write every required neutral identity value.
   Loaders reject omission rather than inferring defaults.
9. Engine tests remain compact numerical falsification checks. Named literature
   replays belong in Validation. No contract, compatibility, snapshot, prose,
   source-specific, or legacy tests are introduced.

## Canonical ownership

| Concern | Owner | Required output |
| --- | --- | --- |
| Parameter identity, values, domains, topology, correlation coefficients | `Parameters` | Resolved ordered `ActiveParameterSet` |
| Fixed-state and pure-saturation quantities | EOS | Observable value, units/basis, domain, exact selected-parameter Jacobian, structured failure |
| Nonreactive solved phase quantities | Equilibrium | Certified phase result and exact total sensitivity from EOS partials |
| Observations, residual scaling, bounds, ties, multistart, Ceres, fit receipt | Regression | Immutable `FitResult` and accepted `Parameters` |
| Literature and application acceptance | Validation | Installed-wheel evidence |
| Reusable parameters and observations | Data | Explicit source-traceable packets |

No owner may duplicate an upstream equation, parameter identity registry, or
observation-family switch merely to make dispatch convenient.

## Shared domain model

### Active parameter identity

Each coordinate carries a stable physical family and one resolved scope:

- component parameter;
- unordered component-pair parameter;
- association parameter or typed cross-association correction;
- model/submodel coefficient; or
- coefficient of a closed, typed temperature correlation.

Resolution binds component indices, pair endpoints, association edges,
correlation coefficient position, units, physical domain, source value, and
applicability. An unordered pair is canonicalized once and applied
reciprocally. Multiple coordinates from any supported families may be active in
one set. Discrete topology and model selection are never optimizer coordinates.

### Observable derivative block

An owning module returns a block with:

- source-neutral observable identity and version;
- ordered parameter identities copied from the active set;
- row identities, values, units or basis, and validity domain;
- the exact derivative order and tensor layout the owner supports;
- a dense or declared sparse Jacobian with fixed row/column order;
- state/reference/phase metadata required to interpret the values; and
- structured per-row or per-solve failures without fabricated values.

Regression maps measurements to those rows and applies declared scaling,
weighting, and loss. The block contains no observed targets or optimizer policy.

## Dependency order

```text
Parameters identity and resolution
        |
        v
EOS generic active-parameter tape and observable blocks
        |                         |
        v                         v
Regression EOS cutover      Equilibrium provider/sensitivity cutover
                                  |
                                  v
                         one phase_equilibrium owner
                                  |
                                  v
                 later Equilibrium-observation admission review
```

Work must follow that order. In particular, Regression must not invent a
temporary Equilibrium evaluator, and Equilibrium must not introduce a second
active-parameter representation while waiting for EOS.

## Execution sequence

### Phase 0 — Freeze the vocabulary and retirement ledger

- Accept ADR 0004 and the new Equilibrium ADR.
- Make `CONTEXT.md`, `ARCHITECTURE.yaml`, and these plans the only live design
  authority.
- CAP-11 creates
  `engine/docs/plans/2026-08-07-general-scientific-capability-retirement-ledger.json`
  before deletion begins. Its versioned entries declare exact prohibited
  symbols/path/link patterns, scan scopes, provenance exceptions with reasons,
  and expected terminal counts. The independent audit records the ledger hash,
  repository tree, scanned-file-list hash, every match/exception, and pass/fail.
- Remove superseded plans, designs, implementation reviews, and source-specific
  runtime guidance. Retain citations, research, source locators, raw data,
  receipts, results, and numerical evidence as provenance.

Completion: every remaining live instruction points to the new owner model and
no deleted document is linked as current guidance.

### Phase 1 — Canonical parameter identity

- Implement the resolved active-coordinate descriptor in `Parameters`.
- Support mixed coordinate sets and arbitrary component counts.
- Add constant `k_ij` first while keeping typed correlation coefficients in the
  same identity system.
- Add typed unordered `k_HB,ij` correction resolution whose identity includes
  the unordered component pair and unordered association-role pair. Apply it to
  every edge matching that declared role pair with exact chain-rule
  propagation; reject an omitted role pair when topology admits more than one.

Completion: one active set is accepted unchanged by EOS, Equilibrium, and
Regression preflight; omissions and unsupported coordinates fail closed.

### Phase 2 — EOS derivative and reference cutover

- Execute the EOS owner plan, including the arbitrary-width AD tape, canonical
  Born formulation, generic reference-dependent quantities, and deletion of
  dormant aqueous-only callbacks.
- Establish the generic tape before admitting Born/reference coordinates; CAP-04
  defines those physical semantics and then connects them to that same tape.

Completion: EOS supplies the initial observable derivative blocks and passes
the EOS plan's numerical matrix.

### Phase 3 — Regression owner-consumption cutover

- Execute the Regression owner plan.
- Delete duplicated Python/C++ observation registries and native EOS/saturation
  implementations.

Completion: every admitted residual row comes from an EOS-owned block; all
parameter identities and Jacobian columns remain unchanged end to end.

### Phase 4 — Equilibrium provider and sensitivity cutover

- Pass the same active set through the Equilibrium provider.
- Replace component-only transport and hidden solvent/reference conventions.
- Compute solved-state sensitivities by exact implicit linear algebra from EOS
  partial tensors.

Completion: arbitrary supported component, pair, association, and model
coordinates propagate through one provider path.

### Phase 5 — General nonreactive phase equilibrium

- Execute the Equilibrium owner plan.
- Replace neutral/electrolyte and source-fingerprint routing with one physical
  owner and explicit discovery/hint/constrained modes.

Completion: one to three liquid/vapor phases are handled without
source-specific branches, and every success states exactly what topology was
certified.

### Phase 6 — Admission review, not automatic expansion

- Review Equilibrium failure semantics, derivative conditioning, phase
  appearance/disappearance, and row-level observability.
- Admit an Equilibrium observation family to Regression only through a new
  evidence-backed decision. Do not prebuild dormant families.

Completion: either a narrowly admitted owner block is documented and tested or
the boundary remains EOS-only with no partial implementation.

## Numerical verification matrix

The compact Engine matrix must include:

- reordered binary components with identical values and permuted derivatives;
- ternary electrolyte and multisolvent/multi-ion compositions;
- several simultaneous `k_ij` coordinates and mixed-family active sets;
- `k_HB,ij` reciprocal pair identity with role-aware edge application;
- canonical Born neutral-limit reduction and non-neutral correction response;
- owner Jacobians against independent central differences at interior points;
- symmetry/order invariance and no fixed solvent/cation/anion positions;
- one-, two-, and three-phase Equilibrium certificates;
- neutral and charged systems through the same phase owner; and
- structured non-evaluable states without penalty substitution.

Finite differences are independent test oracles only. Production derivatives
remain CppAD partials plus exact implicit linear algebra.

## Readiness and validation

At each implementation phase:

```bash
bash tools/check-affected.sh
bash tools/check-all.sh
bash "$HOME/.codex/hooks/codex-cleanup.sh" --repo-root .
```

Also run focused schema/JSON/YAML checks, Ruff, strict mypy, native compiler
warnings, clang-tidy on changed native translation units, sanitizers for native
memory-risk changes, wheel build/origin verification, and the relevant
installed-wheel Validation campaign. Readiness checks are not scientific tests.

## Issue-ready decomposition

These rows are planning units only. Do not publish tracker issues until the
four canonical plans are approved.

| ID | Work unit | Depends on | Acceptance |
| --- | --- | --- | --- |
| CAP-01 | Canonical active-coordinate descriptor | Phase 0 | Mixed ordered active set resolves once with exact identity and domain |
| CAP-02 | Unordered pair and typed `k_HB,ij` resolution | CAP-01 | Reciprocal application and exact chain rule pass permutation checks |
| CAP-03 | EOS generic derivative core | CAP-02 | Arbitrary-width pre-Born active sets return exact order-declared derivative blocks; CAP-04 admits Born/reference coordinates |
| CAP-04 | Canonical Born and explicit reference inputs | CAP-03 | Neutral limit and explicit-reference numerical checks pass |
| CAP-05 | EOS observable block interface | CAP-03 | Fixed-state and pure-saturation owners expose one block contract |
| CAP-06 | Regression consumes owner blocks | CAP-05 | No thermodynamic dispatch or EOS implementation remains in Regression |
| CAP-07 | Equilibrium provider parameter transport | CAP-03, CAP-04 | Same active set and explicit Born/reference inputs reach all EOS provider calls |
| CAP-08 | Equilibrium exact solved-state sensitivities | CAP-07 | Implicit Jacobians match independent differences away from singularities |
| CAP-09 | One `phase_equilibrium` owner | CAP-07 | Discovery and constrained 1–3 phase solves use one source-neutral path |
| CAP-10 | Equilibrium observation admission review | CAP-08, CAP-09 | Explicit admit/defer decision with numerical-failure evidence |
| CAP-11 | Final dead-surface and documentation audit | CAP-04, CAP-06, CAP-10 | No retired symbols, schemas, branches, docs, or tests remain |

## Done definition

The cutover is complete only when all three owner plans are complete, the wheel
contains one parameter and derivative system, source- or case-specific runtime
routes are absent, compact numerical evidence passes, installed-wheel
Validation owns literature acceptance, and a repository-wide symbol/path audit
finds no live dependency on retired guidance or implementations.

# Owner-Driven Parameter Fit Plan

- Status: canonical owner plan; execution is tracked by GitHub root issue #10 and its native CAP subissues
- Date: 2026-08-07
- Owner: Regression
- Parent: [General scientific capability hard cutover](../../plans/2026-08-07-general-scientific-capability-hard-cutover.md)
- Decision: [ADR 0004](../../adr/0004-one-active-parameter-and-observable-system.md)

## Outcome

Make Regression a small generic optimizer over owner-provided thermodynamic
observable blocks. Regression owns fit declaration, observation mapping,
residual scaling, weights, bounds, ties, multistart, Ceres, diagnostics, and
immutable parameter transfer. EOS and later Equilibrium own all scientific
values, domains, solves, and exact selected-parameter derivatives.

The hard-cutover boundary remains deliberately narrow: EOS fixed-state and
pure-saturation observations only. CAP-10 completed the separate numerical-
failure review after the general phase owner and sensitivity seam stabilized;
Governance D-031 defers every current Equilibrium Observation Family.

## What remains

- Public API: `ParameterFit`, `FitResult`, and `fit`.
- One strict versioned fit schema shared by JSON and in-memory mappings.
- Explicit ordered active parameter identities, values, bounds, scales, and
  ties.
- Explicit data locators/hashes, observation mapping, units, residual basis,
  weighting, loss, partitions, and multistart policy.
- One private native Ceres core and complete fit diagnostics.
- Owner-published closed typed catalogs for parameter correlations and admitted
  observable descriptors, queried read-only by Regression. Regression owns no
  duplicate thermodynamic registry or family switch.
- Immutable `FitResult.to_parameters()` transfer.

## What is deleted

- Python `SCHEMA_REGISTRY` as an independent thermodynamic capability registry.
- Observation-family branches in `_fit.py` that calculate EOS or saturation
  values.
- Native `ObservationKind`, `observation_kind`, `target_index`, and duplicate
  family switches.
- Native fit code that directly implements EOS fixed-state or saturation
  equations.
- Binary/single-pair assumptions in fixed-state fitting.
- Source-, author-, campaign-, parameter-tuple-, and application-specific fit
  functions, files, enums, schemas, examples, and tests.
- Any compatibility parser, alias, inferred wildcard selection, or fallback
  finite-difference Jacobian.

## Owner-provided capability

An observation descriptor is advertised by EOS or Equilibrium and contains:

- stable source-neutral family/version;
- required thermodynamic inputs and explicit reference/phase semantics;
- output row identities, units/basis, and domain;
- supported active-parameter identities/topologies and derivative order;
- one batch evaluator returning values, exact Jacobian columns in the supplied
  `ActiveParameterSet` order, and structured failures.

Regression resolves user observation rows against this descriptor. It never
chooses a native target index, reaches through `model._native`, or reconstructs
an owner's parameter capability.

## Fit execution

1. Parse the canonical fit document strictly.
2. Resolve the base `Parameters` and one explicit `ActiveParameterSet`.
3. Resolve each observation block against an admitted owner descriptor.
4. Validate units, bases, reference data, active-coordinate applicability,
   bounds, ties, scales, and row hashes before entering Ceres.
5. At each trial, materialize the physical active values and call each owner
   block in batches.
6. Convert returned predictions to declared dimensionless residuals and apply
   measurement weighting/loss.
7. Route the owner's exact Jacobian through the same scaling and ties into
   Ceres.
8. Record starts, failures, convergence, objective blocks, rank/conditioning,
   active bounds, predictions, residuals, Jacobian, and accepted parameters in
   one immutable result.

An owner evaluation may fail at a bounded trial. Regression records the typed
failure and allows Ceres or multistart policy to reject/backtrack that trial; it
never fabricates a scientific penalty value.

## Initial observation boundary

### Admitted

- EOS-owned fixed-state quantities with explicit state, composition, and basis.
- EOS-owned pure-saturation pressure, liquid density/volume, and vapor
  density/volume, with one shared saturation solve per state.
- Direct auxiliary EOS quantities only when EOS exposes a complete descriptor,
  domain, and exact selected-parameter derivative.

### Deferred

- mixture VLE/LLE/flash quantities;
- phase appearance or phase-count observations;
- reactive or chemical-equilibrium quantities;
- reaction-model coefficients;
- total caloric properties before an EOS thermochemistry model exists; and
- arbitrary user callbacks or expressions.

Deferred means absent, not represented by dormant enums or placeholder schema.

## CAP-10 Equilibrium admission decision

CAP-10 reviewed three candidate families separately and admits none:

- **Omitted-topology discovery:** deferred from the differentiable fit surface.
  Phase appearance or disappearance changes the phase count and row identities,
  retained discovery establishes only a finite-search claim, and competing
  branch/global selection is not a smooth owner block.
- **Declared fixed-T,P phase split:** deferred. Exact implicit sensitivities are
  available for certified interior fixed-topology solutions, but phase amount
  limits, phase coalescence, support activation, singular or ill-conditioned
  residual systems, and competing local branches do not yet have one
  family-level optimizer-domain policy.
- **Fixed-T coexistence boundary:** deferred as the strongest future candidate.
  It has a scalar pressure output, explicit ordered topology and support, and
  exact local pressure/state sensitivities, but Equilibrium does not yet expose
  a batch Observation Family descriptor/block or demonstrate branch identity,
  nonuniqueness, support-boundary behavior, and typed Non-Evaluable Trials over
  a proposed fit domain.

The reviewed owner contracts are explicit even though they are not fit
families:

- Discovery takes temperature in kelvin, pressure in pascals, and positive feed
  amounts in moles. It returns a variable one-/two-phase `FeedSplit`: phase
  amounts in moles, fractions, dimensionless compositions and chemical
  potentials over RT, molar volumes in m3/mol, and total reduced free energy.
  Its domain is the retained binary-neutral or charged HELD/HELD2 search scope.
  Applicable selected parameters retain canonical order, but output identities
  change with the discovered topology. Search exhaustion, uncertified states,
  and unavailable/singular/ill-conditioned sensitivity are explicit failures;
  globality remains `not_established`.
- A declared fixed-T,P split takes the same dimensional state plus one ordered
  topology of one to three all-component phases with at most one vapor. It
  returns the fixed `FeedSplit` rows for that topology and exact derivatives in
  selected-parameter order only when the defining residual Jacobian is full
  rank and acceptably conditioned. Material, charge, pressure, KKT/bound,
  modified-potential, distinctness, mechanical-stability, and ordered-role
  certificates are mandatory; nonconvergence, certificate failure,
  coalescence, singularity, and ill-conditioning fail closed.
- Fixed-T coexistence takes temperature in kelvin, an ordered two-phase
  topology, a reference phase index, and its dimensionless composition. The
  reference has all-component support; the partner has all-component or
  structural neutral-component support. It returns an `IncipientBoundary`
  pressure in pascals plus dimensionless intensive compositions and potentials,
  molar volumes in m3/mol, and exact selected-parameter pressure/state
  derivatives. Common-pressure, potential, charge, stability, role-order, and
  distinctness certificates are mandatory. Missing certified roots,
  nonuniqueness outside the finite candidate claim, support changes,
  coalescence, singularity, and ill-conditioning are non-evaluable conditions.

No partial Equilibrium observation schema, enum, resolver branch, or evaluator
exists in Regression, so deferral requires no runtime deletion and no new
numerical test. A future admission must be a new implementation issue. It must
define inputs, outputs, units and basis, stable row identities, active-parameter
support, conditioning bounds, branch and nonuniqueness policy, certificate
requirements, and typed failure semantics, then falsify the complete value and
Jacobian block independently across its declared optimizer domain. Equilibrium
must own the batch values and total derivatives; Regression may own only row
mapping, residuals, weights, and family-level trial policy at that owner seam;
its existing fit-wide bounds, ties, multistart, Ceres, and diagnostics ownership
remains unchanged.

The decision reuses the compact evidence in
`tests/equilibrium/test_phase_equilibrium.py` for the discovery search scope and
explicit `not_established` globality claim, certified
fixed-topology branches, exact coexistence/split sensitivities, and phase-order
behavior, plus `tests/equilibrium/test_solved_state_sensitivity.py` for centered
re-solve agreement and fail-closed coalesced and boundary states. The native
rank and condition gates remain part of the reviewed failure contract without a
new CAP-10 claim. Existing tests already falsify the scientific premises of the
review; CAP-10 adds no new numerical claim that warrants another test.

## Correlations

A typed parameter correlation is part of `Parameters`; its coefficients are
ordinary parameter identities. A published observation correlation is an
owner-admitted source representation with explicit outputs, units, domain, and
derivatives. Regression may integrate a real multivariate correlation only
when the owner descriptor defines its measure and outputs. The existing
one-target restriction is generalized only for such a concrete admitted case,
not preemptively.

## Implementation sequence

1. **Descriptor seam.** Define the internal owner descriptor/block protocol in
   the Engine without adding an installed ABI or second package boundary.
2. **Resolver simplification.** Make fit resolution query owner descriptors and
   preserve the one active-set identity/order. Delete Regression-owned
   thermodynamic capability registries.
3. **EOS fixed-state adapter.** Route all admitted fixed-state rows through the
   EOS generic block; remove binary/single-pair assumptions.
4. **Pure-saturation adapter.** Route all saturation rows through the EOS owner
   and share a solve for co-located observations.
5. **Native core reduction.** Leave only residual scaling/ties and Ceres-facing
   plumbing in Regression native code; delete EOS/saturation implementations
   and family/target switches.
6. **Result consolidation.** Preserve owner row identities, active parameter
   identities, values, residuals, Jacobians, and failures in `FitResult`.
7. **Dead-surface cleanup.** Remove old source-specific docs/tests and update the
   quickstart/science guide to the owner-driven model.
8. **Admission gate.** Completed by CAP-10 and Governance D-031: all current
   Equilibrium families are deferred and no dormant fit surface is added.

## Numerical verification

Compact Regression checks must establish:

- simultaneous fitting of several arbitrary `k_ij` pairs in a multicomponent
  fixed-state problem;
- a mixed-family active set where every owner-provided Jacobian column reaches
  Ceres in the same order;
- component reordering changes only declared row/coordinate permutations;
- owner exact Jacobians and residual-scaled Jacobians agree with independent
  central differences of the complete residual at interior points;
- parameter ties apply one optimizer coordinate to several exact identities;
- pure-saturation rows share one owner solve and preserve distinct outputs;
- non-evaluable trials remain typed diagnostics and do not become penalties;
- rank deficiency and active bounds remain visible; and
- accepted parameter transfer does not mutate the input document.

Tests may use compact manufactured or carefully selected numerical cases. Named
literature recovery, broad temperature tables, and workflow acceptance belong
in installed-wheel Validation.

## Issue-ready decomposition

Do not publish these units until the four canonical plans are approved.

| ID | Work unit | Depends on | Acceptance |
| --- | --- | --- | --- |
| REG-01 | Owner observable descriptor/block protocol | EOS-07 | One typed internal protocol carries values and ordered Jacobians |
| REG-02 | Descriptor-driven fit resolver | REG-01 | No Regression thermodynamic registry or target index remains |
| REG-03 | Generic EOS fixed-state residual block | REG-02 | Arbitrary supported active sets and component counts fit |
| REG-04 | EOS pure-saturation residual block | REG-02 | Co-located outputs share one solve and exact sensitivities |
| REG-05 | Native Ceres core reduction | REG-03, REG-04 | Native code owns optimization only |
| REG-06 | FitResult diagnostics and transfer | REG-05 | Complete identity/value/Jacobian/failure receipt is immutable |
| REG-07 | Regression dead-surface audit | REG-06 | Duplicated dispatch, source routes, and legacy docs/tests are absent |
| REG-08 | Equilibrium admission review | EQ-08 | Completed by Governance D-031: all current families deferred; no dormant code |

## Done definition

Regression contains no thermodynamic equation, phase solver, parameter-family
derivative implementation, native target switch, or source-specific route. It
optimizes any admitted EOS block over one `ActiveParameterSet`, reports exact
owner Jacobians transparently, and keeps Equilibrium outside the current fit
surface until a separately shaped implementation issue satisfies D-031.

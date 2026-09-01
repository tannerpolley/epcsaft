# Owner-Driven Parameter Fit Plan

- Status: completed owner contract; native Ceres cutover accepted by D-042
- Date: 2026-08-07
- Amended: 2026-09-01 for issue #153
- Owner: Regression
- Parent: [General scientific capability hard cutover](../../plans/2026-08-07-general-scientific-capability-hard-cutover.md)
- Decision: [ADR 0004](../../scientific/adr/0004-one-active-parameter-and-observable-system.md)

This document is the retained scientific owner contract and CAP provenance. It
is not a live execution queue; current work is selected from GitHub Issues.

## Outcome

Make Regression a small generic Ceres optimizer over opaque owner-provided
thermodynamic observable handles. Regression owns ordered observations,
caller-defined residual scaling, physical bounds, affine coordinates, solver
controls, diagnostics, and immutable result mapping. EOS and Equilibrium own
all scientific values, domains, solves, and exact selected-parameter
derivatives.

The hard-cutover boundary began with EOS fixed-state and pure-saturation
observations. Governance D-034 admits fixed-T,P homogeneous reactive
observations, and D-035 admits fixed-topology reactive bubble/VLE observations.
D-036 approves the generic GREPE evaluator target but admits no additional
Observation Family without separate evidence.

## What remains

- Public API: owner compilers plus `least_squares`, `curve_fit`, and
  `LeastSquaresResult`.
- Opaque immutable owner handles with sidecars carrying the complete row and
  parameter contract.
- Explicit ordered active parameter identities, physical values, bounds, and
  affine origin/scale.
- Ordered observations with optional scalar/vector diagonal discrepancy
  scales and linear loss.
- One private native Ceres core and complete fit diagnostics, including both
  physical and weighted optimizer-coordinate Jacobians.
- Owner-published closed typed catalogs for parameter correlations and admitted
  observable descriptors, queried read-only by Regression. Regression owns no
  duplicate thermodynamic registry or family switch.
- Immutable explicit candidate-parameter transfer remains a caller concern;
  Regression
  assigns no scientific acceptance or adoption status.

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

1. Validate one explicit owner compiler request and one `ActiveParameterSet`.
2. Validate each opaque sidecar's parameter contract, row order, units, bases,
   topology, and branch identity before entering Ceres.
3. Validate physical starts, bounds, affine coordinates, and caller-defined
   residual scales.
5. At each trial, materialize the physical active values and call each opaque
   owner handle in declared row order.
6. Form caller-defined scaled difference residuals; scales are not promoted to
   measurement uncertainty.
7. Route each owner's exact physical prediction Jacobian through the affine
   and residual chain into Ceres while retaining both Jacobian forms.
8. Record one-start failures, convergence, resource counters, active bounds,
   predictions, residuals, Jacobians, owner sidecars, and covariance status in
   one immutable `LeastSquaresResult`.

An owner evaluation may fail at a bounded trial. Regression records the typed
failure and allows Ceres to reject/backtrack that trial; it
never fabricates a scientific penalty value.

## Initial observation boundary

### Admitted

- EOS-owned fixed-state quantities with explicit state, composition, and basis.
- EOS-owned pure-saturation pressure, liquid density/volume, and vapor
  density/volume, with one shared saturation solve per state.
- Direct auxiliary EOS quantities only when EOS exposes a complete descriptor,
  domain, and exact selected-parameter derivative.
- Equilibrium-owned fixed-T,P homogeneous reactive rows with immutable local
  continuation anchors, exact total Jacobians, and diagnostic-only results.
- Equilibrium-owned fixed-topology reactive bubble/VLE rows with immutable
  local branch results, exact total Jacobians, and typed failures.

### Deferred

- unadmitted mixture VLE/LLE/flash quantities;
- phase appearance or phase-count observations other than the D-035 declared
  incipient-vapor family;
- reactive LLE and unadmitted nonreactive phase-equilibrium quantities;
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
its existing physical bounds and diagnostics ownership remains explicit; there
is no start list, basin selection, campaign stopping rule, or alternate
optimizer lifecycle.

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

1. **Owner seam.** Compile EOS and already-admitted D-034/D-035 Equilibrium
   requests into opaque handles with immutable sidecars.
2. **Native core.** Keep only owner-neutral residual evaluation, exact Jacobian
   transforms, Ceres controls, typed failures, cache accounting, and final
   validation in Regression.
3. **Public mapping.** Expose `least_squares` and `curve_fit` with one physical
   start, caller-defined diagonal scales, and a frozen `LeastSquaresResult`.
4. **Dead-surface cleanup.** Delete the thermodynamic/IPOPT lifecycle and its
   active callers, while retaining completed numerical evidence as provenance.
5. **Admission gate.** Preserve D-031's deferred discovery/globality limits;
   D-034 and D-035 remain the only admitted Equilibrium owner families.

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
- active bounds remain visible, while singular covariance truncation remains
  Ceres-owned without an Engine rank label; and
- candidate parameter export does not mutate the input document or assert
  scientific adoption.

Tests may use compact manufactured or carefully selected numerical cases. Named
literature recovery, broad temperature tables, and workflow acceptance belong
in installed-wheel analysis sets.

## Issue-ready decomposition

Do not publish these units until the four canonical plans are approved.

| ID | Work unit | Depends on | Acceptance |
| --- | --- | --- | --- |
| REG-01 | Owner observable descriptor/block protocol | EOS-07 | One typed internal protocol carries values and ordered Jacobians |
| REG-02 | Descriptor-driven fit resolver | REG-01 | No Regression thermodynamic registry or target index remains |
| REG-03 | Generic EOS fixed-state residual block | REG-02 | Arbitrary supported active sets and component counts fit |
| REG-04 | EOS pure-saturation residual block | REG-02 | Co-located outputs share one solve and exact sensitivities |
| REG-05 | Native Ceres core reduction | REG-03, REG-04 | Native code owns owner-neutral optimization only |
| REG-06 | LeastSquaresResult diagnostics | REG-05 | Complete identity/value/Jacobian/failure result is immutable |
| REG-07 | Regression dead-surface audit | REG-06 | Duplicated dispatch, source routes, and legacy docs/tests are absent |
| REG-08 | Equilibrium admission review | EQ-08 | D-031 deferred phase families; D-034 admits only fixed-T,P homogeneous reactive observations |

## Done definition

Regression contains no thermodynamic equation, phase solver, parameter-family
derivative implementation, native target switch, or source-specific route. It
optimizes admitted owner blocks over one `ActiveParameterSet`, reports exact
owner Jacobians transparently, and admits Equilibrium families only through a
separately shaped issue and decision such as issue #30/D-034.

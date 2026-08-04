# Literature Reproduction Contract

Status: authoritative Regression workflow and evidence-classification contract;
no parameter, prediction, or promotion authority.

Date: 2026-07-30

## Purpose

Regression reproduces a declared numerical contract. It does not promise that a
published parameter tuple can be recovered from a citation alone.

Literature often omits one or more of the exact rows, correlation sampling
grid, objective, residual normalization, property weights, loss function,
outlier treatment, bounds, starts, optimizer, tolerances, model version,
association topology, combining rule, or binary-parameter convention. Missing
choices must be reported as missing and may be replaced only by explicit
reconstruction decisions.

This contract applies to accepted paper-specific workflows and authority-neutral
general Regression. It complements
`docs/science/general-parameter-regression.md`; it does not broaden the
parameter-family or observation capabilities advertised by an installed EOS.

## Reproduction classes

Every result that claims literature reproduction uses exactly one class:

| Class | Required meaning |
|---|---|
| `EXACT_AUTHOR_METHOD_REPLAY` | The exact author data, model, objective, numerical contract, and reported comparison are available and executed without an undisclosed replacement. |
| `SOURCE_FAITHFUL_RECONSTRUCTION` | The governing model, source data, and intended objective meaning are retained, but at least one unavailable author implementation detail is replaced by an explicit, versioned, source-constrained reconstruction decision. |
| `PUBLISHED_TUPLE_PROPERTY_REPLAY` | A published parameter tuple is evaluated against declared reported or independently acquired properties; no parameter recovery is claimed. |
| `MODERN_REFIT` | Source-bound data are fitted after deliberately changing the model, row selection, estimand, objective, weights, loss, or another scientific fitting choice. Agreement or disagreement with a published tuple is descriptive. |

The class is an evidence description, not a quality ranking. A well-designed
`MODERN_REFIT` can be scientifically more useful than an incompletely reported
historical fit, but it is not an author-method replay.

Apply the classes in this order so they do not overlap:

1. An evaluation of a published tuple with no parameter fitting is
   `PUBLISHED_TUPLE_PROPERTY_REPLAY`.
2. A fit with the complete author contract and no replacement is
   `EXACT_AUTHOR_METHOD_REPLAY`.
3. A fit intended to recover the author method is
   `SOURCE_FAITHFUL_RECONSTRUCTION` only when each missing implementation
   detail is replaced by the narrowest source-constrained decision while the
   model, rows, and objective meaning remain unchanged.
4. Any newly selected scientific fitting choice, including an objective that
   cannot be reconstructed from the source's intended metric, makes the result
   `MODERN_REFIT`.

An unknown optimizer or start policy alone may be reconstructed under item 3;
changing what is fitted or how discrepancies are scientifically scored falls
under item 4. A fixed-parameter calculation without a published tuple is not
automatically tuple replay: it uses items 2 or 3 only when it replays or
reconstructs an author method, and otherwise remains an ordinary source-bound
calculation without a literature-reproduction claim.

## Complete replay contract

### Model identity

Retain:

- EOS package version, commit/tree, installed wheel/header/library hashes, and
  capability fingerprint;
- parameter identities, units, order, transforms, fixed values, and component
  order;
- formulation, phase convention, reference state, association scheme, site
  topology and multiplicities, and topology fingerprint;
- mixing and combining rules, including the sign and placement of `k_ij`,
  `l_ij`, resolved cross-association parameters, or any source-defined
  transform;
- treatment of polar, electrolyte, caloric, and temperature-correlation terms.

Discrete chemistry, formulation, and topology are specifications. Regression
does not infer or optimize them.

### Data identity

Retain:

- citation, durable row/table/figure locator, source artifact hash, license or
  use basis, and acquisition date where material;
- exact raw observations when available, including units, bases, uncertainty,
  reported precision, component order, and phase labels;
- whether each row is a direct measurement, digitized value, database record,
  author correlation, or newly reconstructed correlation value;
- for correlation-generated rows, the complete equation, coefficients, units,
  validity interval, sampling grid, and transformation record;
- duplicate, exclusion, critical-region, censoring, and outlier decisions;
- canonical transformed-dataset hash and immutable training, held-out, and
  stress partitions.

Correlation-generated points do not become independent experimental evidence
merely because a dense grid was sampled.

### Objective identity

Retain:

- residual equation for every observation family;
- absolute, relative, logarithmic, aggregate, or censored interpretation;
- property and row weights, scales, covariance assumptions, and uncertainty
  interpretation;
- squared, absolute, Huber, or other finite supported loss and every loss
  parameter;
- handling of failed, skipped, out-of-model, and nonfinite evaluations;
- the exact scalar objective used to rank candidate solutions.

If Regression cannot express the author objective through its finite typed
vocabulary, the result cannot be `EXACT_AUTHOR_METHOD_REPLAY`. A modern
least-squares or robust-loss replacement is a declared reconstruction or
modern refit, never a silent compatibility approximation.

### Numerical identity

Retain:

- physical bounds and affine scales;
- complete primary and confirmation start vectors;
- solver and linear-solver configuration, tolerances, iteration/evaluation
  limits, deterministic seeds, and multistart policy;
- nuisance-coordinate origins, starts, bounds, and closure tolerances;
- derivative mode and independent derivative-oracle configuration;
- termination, active-bound, rank, singular-value, conditioning, row-accounting,
  and confirmation diagnostics.

When the author numerical method is unknown, Regression may demonstrate
deterministic recovery under its declared method but must not attribute that
method to the source.

### Comparison identity

Retain:

- full-precision fitted values and the source's printed/rounded values
  separately;
- raw and scaled residuals by property family;
- source-reported summary metrics recomputed from the declared rows where
  possible;
- parameter differences using precision-compatible comparisons;
- whether `k_ij` or another mixture coordinate was zero, fixed, independently
  sourced, or refitted;
- solver, numerical, physical/workflow, scientific, predictive, and authority
  statuses as separate axes.

A different parameter tuple that reproduces fitted properties is a
property-equivalent result unless stronger identifiability evidence exists.

## Association-specific workflow

For a new species, use this order:

1. Classify the species and intended mixture domain as nonassociating,
   self-associating, cross-associating, or capable of induced association.
2. Select each candidate association scheme and site topology from source,
   chemical, spectroscopic, molecular-simulation, or separately declared
   engineering evidence. Treat each as a fixed model alternative.
3. Fit continuous pure parameters only within an installed advertised
   topology. Select typed ordinary-pure and component/site-pair slots from the
   immutable EOS descriptor. Energy-only, volume-only, joint, constrained
   subset, and explicit shared-parameter selections all use the same generic
   request/result path. A 2B five-parameter selection is one reference case,
   not a separate runtime capability.
4. Evaluate exact scaled sensitivities, nuisance-projected rank, conditioning,
   active bounds, multistart behavior, and nuisance-reoptimized association
   profiles. Solver convergence and local full rank alone do not establish
   practical identifiability.
5. If the accepted profile remains open, retain a property-equivalent family or
   use an explicitly constrained parameterization. A fixed association volume,
   externally estimated association energy, homologous-family relation, or
   other constraint must carry its own source and validity domain.
6. Validate mixture transfer with the `k_ij` treatment visible. A refitted
   `k_ij` is correlation evidence, not predictive transfer.
7. Treat induced and resolved cross association as fixed mixture-model
   alternatives requiring matching EOS capabilities and source-bound mixture
   observations. Compare, where advertised, (a) no induced association with
   declared `k_ij` treatment, (b) a predictive source-defined induced-
   association combining rule with no fitted mixture-association coordinate,
   (c) one resolved cross-association energy/volume pair, and (d) a
   dispersive-`k_ij`-only fit without induced association. They are not
   additional pure-component coordinates.
8. Reject a simultaneous `k_ij` plus resolved cross-association fit unless the
   declared data design, exact scaled Jacobian, nuisance-reoptimized profiles,
   and accepted-region evidence distinguish the complete active block.

Current Regression issue 28 owns practical identifiability for the retained
fixed-pure-2B literature case.
Issue 34 owns the future fixed induced/cross-association alternatives after EOS
support. Neither issue authorizes dynamic topology discovery.

Primary method anchors include:

- [Gross and Sadowski (2002)](https://doi.org/10.1021/ie010954d) for the joint
  five-parameter pure associating fit;
- [Clark et al. (2006)](https://doi.org/10.1080/00268970601081475) for
  association/dispersion parameter degeneracy and topology discrimination;
- [Kleiner and Sadowski (2007)](https://doi.org/10.1021/jp072640v) for one
  predictive induced-association construction;
- [Albers, Heilig, and Sadowski
  (2012)](https://doi.org/10.1016/j.fluid.2012.04.014) and [Fuenzalida et al.
  (2016)](https://doi.org/10.1016/j.fluid.2016.07.001) for conditional
  reduced-parameter strategies; and
- [Nikolaidis et al.
  (2024)](https://doi.org/10.1021/acs.jced.2c00781) for a broad comparison of
  induced association, fitted cross-association, and `k_ij` strategies.

These references support methods and limits; they do not supply universal
parameter defaults.

## Regression routine-test contract

Regression keeps small deterministic tests that continuously prove the package
can still execute its main generic numerical workflows and supported active
parameter variations. These tests run through public or installed-artifact
interfaces and remain suitable for routine CI.

The compact matrix covers:

1. accepted methane and ethane joint pure-saturation numerical results;
2. representative scalar component, pair, and direct-observable fits across
   each distinct native observation/derivative adapter;
3. joint pure `(m, sigma, epsilon/k)` recovery;
4. descriptor-driven fixed-topology association mechanics covering the 2B
   reference, a non-2B multi-pair topology, `N_active = 1, 2, >2`, constrained
   subsets, sharing, permutation rejection, and fail-closed identity checks;
5. constant `k_ij` and `l_ij` fixed-composition VLE mechanics;
6. the composed positive-observation evaluator transport; and
7. deterministic result serialization once the general record exists.

Each retained numerical sentinel:

- uses a compact manufactured analytic case or the smallest source-backed
  anchor that uniquely tests the workflow;
- verifies exact production derivatives against an independent analytic,
  automatic-differentiation, or bounded numerical oracle where needed;
- checks fitted values with
  `abs(value - reference) <= atol + rtol * abs(reference)`;
- justifies `atol` and `rtol` from parameter scale, oracle accuracy, source
  precision, and deterministic numerical behavior;
- checks solver termination, complete row evaluation, projected rank,
  conditioning policy, active bounds, and confirmation-start agreement as
  applicable; and
- never generates the reference value with the implementation path under test.

The matrix is parameterized where cases share one scientific invariant. It is
not expanded into one broad literature campaign per parameter family.

## Validation campaign contract

The sibling Validation repository owns durable installed-artifact evidence
whose breadth would make Regression tests slow, source-heavy, or
authority-confusing. Validation owns:

- full literature data grids and exact or reconstructed author objectives;
- broad temperature, pressure, composition, component-family, and property
  matrices;
- practical-identifiability profile surfaces, bootstrap/data-perturbation or
  parameter-distribution studies, and extrapolation uncertainty;
- competing topology, induced-association, cross-association, and `k_ij`
  strategy comparisons;
- mixture, caloric, electrolyte, reactive, held-out, and stress transfer;
- reproducibility receipts, independent review, and promotion evidence.

Every retained Validation campaign, profile surface, uncertainty distribution,
plot, and comparison table must carry an immutable repository/commit/path
reference and content hash together with the exact installed EOS and Regression
artifact identities. Regression result records retain those immutable
references and hashes when the artifacts apply; absence is explicit.

Validation consumes immutable installed EOS and Regression artifacts through
public contracts. It does not import either source checkout, replace routine
Regression mechanics tests, or turn a reconstructed literature campaign into
author-method evidence.

## Acceptance boundary

A literature-related fit is usable only when:

- its reproduction class and complete available contract are retained;
- every replacement decision is explicit;
- exact derivatives and all required rows are evaluated or the fit fails
  closed;
- numerical and physical diagnostics pass their declared gates; and
- scientific, predictive, and authority claims do not exceed the retained
  evidence.

No reproduction class alone persists parameters to the EOS catalog, establishes
global uniqueness, proves prediction, or grants production authority.

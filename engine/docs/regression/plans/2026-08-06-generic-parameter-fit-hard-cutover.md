# Generic ePC-SAFT Parameter Fit: Canonical Design and Hard-Cutover Plan

Status: canonical pre-implementation design

Date: 2026-08-06

Doctrine: Governance revision 12

Owner: Engine

Execution ledger:
[`2026-08-06-generic-parameter-fit-execution.md`](2026-08-06-generic-parameter-fit-execution.md)

## Authority and reading order

This is the current design and execution authority for replacing the existing
Regression implementation. It supersedes live design guidance in
`../science/general-parameter-regression.md`, the current
`../general-regression-quickstart.md`, and source-specific runtime designs.
Those files remain implementation history or current pre-cutover usage until
the final cutover; they do not constrain the replacement.

Read this plan with:

1. `../../../CONTEXT.md` for the canonical vocabulary;
2. `../../../../governance/GOVERNANCE.md` for Numerical
   Verification-Driven Development and monorepo ownership; and
3. `../../../ARCHITECTURE.yaml` for the one-wheel Engine boundary.

The implementation must not introduce a second active public Regression route,
a compatibility period, or a source-specific exception. The branch may build
private replacement internals incrementally, but the public change is one hard
cutover.

## 1. Outcome

Replace the present collection of generic-looking contracts plus specialized
Born, water-factor, aqueous-`k_ij`, pure-saturation, and source-named workflows
with one deep Parameter Fit module:

```python
from epcsaft.regression import FitResult, ParameterFit, fit

result = fit("fit.json")
fitted = result.to_parameters()
```

The same interface must fit any admitted continuous ePC-SAFT coefficient or
combination of coefficients when all of the following are declared:

- one canonical `Parameters` document;
- an exact Active Parameter Set;
- one or more typed Observation Families;
- source data, units, bases, domains, and provenance;
- residual definitions and optional weights;
- physical bounds, affine coordinate scales, and multistart controls; and
- an EOS capability that supplies the value and exact selected-parameter
  derivative.

The optimizer is generic. Components, papers, campaigns, and applications are
data and evidence, never runtime routes.

## 2. Non-negotiable decisions

### 2.1 Runtime scope

- `engine/` remains the sole runtime owner and builds one `epcsaft` wheel.
- `epcsaft.regression` owns fit resolution, residual assembly, Ceres execution,
  diagnostics, and immutable result transfer.
- EOS owns every thermodynamic value, density closure, pure-saturation solve,
  physical domain, and exact selected-parameter derivative.
- The first cutover admits EOS-owned fixed-state and pure-saturation
  observations only.
- Regression must not import or invoke `epcsaft.equilibrium` in this cutover.
- Mixture VLE, LLE, flash, reactive observations, reaction-equilibrium
  constants, and other reaction-model coefficients remain excluded until
  Equilibrium receives a separate interface and numerical redesign.
- Regression fits ePC-SAFT model coefficients only. It does not fit an
  application's reaction chemistry or process-model coefficients.

### 2.2 Public interface

The intended public Regression surface is:

```text
ParameterFit
FitResult
fit
```

`ParameterFit.from_mapping(mapping)` and `ParameterFit.from_json(path)` are two
loaders for the same strict, versioned schema. `fit(...)` may accept a resolved
`ParameterFit`, a matching mapping, or a JSON path and always resolves the same
model. No Python-only shorthand may alter behavior.

No current coordinate classes, observation classes, source-specific records,
workflow result classes, capability callbacks, or native transport details
remain public after cutover.

### 2.3 Naming and provenance

- Runtime symbols, files, parameter families, observations, capabilities, and
  tests use source-neutral physical or numerical names.
- Author and paper names remain valid in citations, source locators,
  Validation campaign paths, historical plans, and receipts.
- Established model names such as Born and SSM+DS remain valid scientific
  names; a `Figiel` solver, callback, parameter family, or Regression workflow
  does not.
- A literature value is evidence for a generic capability, not permission to
  encode the literature source in that capability.

### 2.4 No compatibility layer

This private project has no legacy API requirement. The final cutover deletes
displaced exports, implementations, packaged datasets, tests, and docs in the
same branch. It adds no aliases, wrappers, deprecation warnings, fallback
imports, or old-schema readers.

## 3. Scientific basis reviewed for this design

The following conclusions constrain the design; literature remains provenance
rather than runtime identity.

| Topic | Evidence reviewed | Design consequence |
|---|---|---|
| Pure nonassociating parameters | Gross and Sadowski (2001), section 3.1, reports `m`, `sigma`, and `epsilon/k` fitted to vapor pressure and liquid molar volume, with additional PVT data in some cases; equation 20 states the squared relative `P_sat`/volume objective used in the preceding universal-constant regression. `[verified]` | The standard pure baseline combines saturation pressure and saturated-liquid density or molar volume. |
| Homologous-series relation | Gross and Sadowski (2001), Appendix B, correlates already fitted n-alkane parameters against molar mass. `[verified]` | The Appendix-B relation is not a normal single-component Parameter Fit mode and is not a saturation-pressure polynomial. |
| Saturation source correlations | Existing project data is mainly tabulated; the MEA evidence cites DIPPR saturation-pressure and saturated-liquid-density correlations. `[verified]` | Observation sources may be explicit rows or typed published correlations with equations, coefficients, units, and domains. |
| Permanent moments | Gross (2005), Albers and Sadowski (2011), Albers et al. (2012), Esper et al. (2023), Kleiner and Sadowski (2007), and Vrabec and Gross (2009) use measured or quantum-calculated moments as fixed inputs in their primary parameter workflows. `[verified]` | Dipole and quadrupole moments are fixed Molecular Moments, not ordinary optimizer coordinates. |
| General association equations | Chapman et al. (1989), Huang and Radosz (1990), Veytsman (1998), Müller and Gubbins (2001), and Gross and Sadowski (2002) describe general site-based mass action plus named practical schemes. `[verified]` | Association uses one explicit graph; 1A, 2B, 3B, and 4C are presets, not solvers. |
| Association parameter fitting | Gross and Sadowski (2002) fit base and association parameters against pure vapor-pressure and liquid-density data. `[verified]` | Association energy and volume are ordinary selected coefficients on a fixed topology. |
| Relative permittivity | Rueben et al. (2024), Yang et al. (2024), and Zuber et al. (2014) treat relative permittivity as a separately parameterized constitutive input to primitive electrolyte EOS terms. `[verified]` | Permittivity is an explicit auxiliary EOS submodel with its own typed coefficients and observations. |
| Numerical identifiability | Albers and Sadowski (2011), Esper et al. (2023), and Fuenzalida et al. (2016) show that data choice and sensitivity structure materially affect fitted parameters. `[verified]` | FitResult must expose rank, singular values, conditioning, active bounds, sensitivities, and competing basins; objective reduction alone is insufficient. |

## 4. Canonical domain model

### 4.1 One Parameters document

One strict, versioned Parameters document contains every input needed to
construct the selected EOS model. It is not a flat dictionary. Its typed
sections distinguish:

- component and model specification: identities, molar masses, charges,
  Molecular Moments, fixed association topology, and discrete model-family
  choices;
- coefficient namespaces: base PC-SAFT, pair, association, polar,
  electrolyte, Born, and auxiliary constitutive model coefficients;
- validity domains and units; and
- provenance.

Only declared continuous model coefficients receive Parameter Identities.
Molar mass, charge, component identity, Molecular Moments, site count and
multiplicity, topology, model-family selection, citations, and provenance are
fixed specifications and cannot become optimizer coordinates merely because
they occur in the same serialized document.

`FitResult.to_parameters()` returns a new immutable Parameters value with the
accepted fitted coefficients and fit provenance. It never mutates the input
object or writes a file implicitly.

### 4.2 Model Family Namespaces

Base and auxiliary coefficients share the same Parameters document but remain
in separate source-neutral namespaces. Examples include:

```text
component/<component-id>/segment_count
component/<component-id>/segment_diameter/<correlation>/<coefficient>
component/<component-id>/dispersion_energy_over_k
pair/<component-a>/<component-b>/k_ij/<correlation>/<coefficient>
association/<endpoint-a>/<endpoint-b>/energy_over_k
association/<endpoint-a>/<endpoint-b>/volume
model/relative_permittivity/<family>/<coefficient>
component/<ion-id>/born_diameter
```

These are canonical identities, not JSON paths. Pair and association endpoints
use canonical ordering. The same identity is one coordinate everywhere unless
an explicit Parameter Tie makes distinct identities share a coordinate.

### 4.3 Typed Parameter Correlations

A parameter correlation has a fixed registry-owned family, independent
variables, units, validity domain, coefficient identities, and exact chain
rule. Any declared continuous coefficient may be active independently; the
correlation family itself is fixed during a fit.

The first registry contains only forms required by active models and data. It
does not accept arbitrary expressions or callbacks. Source-named current
coefficient families are renamed by their physical role while citations remain
in provenance.

### 4.4 Active selection

An Active Parameter Set is selected by exact Parameter Identities or typed
selectors over component, pair, association endpoint, and Model Family
Namespace. Preflight must:

1. expand every selector deterministically;
2. reject empty, ambiguous, duplicated, discrete, unsupported, or
   derivative-incomplete selections;
3. apply explicit Parameter Ties;
4. freeze the exact ordered identities in the Resolved Fit; and
5. record the expansion in FitResult.

Arbitrary JSONPath, attribute-path wildcards, stored-field allowlists, and
runtime registration are forbidden.

### 4.5 Pure and Mixture

The root EOS interface becomes explicit:

- `Pure(parameters)` requires exactly one component and owns pure-state and
  pure-saturation evaluation.
- `Mixture(parameters)` requires at least two components and requires an
  explicit composition for every state.
- `Pure.state(...)` and `Mixture.state(...)` return the same `State` type.
- `Pure.saturation(T=...)` returns the common saturation pressure and both
  coexisting densities from pressure and chemical-potential or fugacity
  equality.

The existing fingerprint-approved saturation route in Equilibrium is removed
from live ownership. Pure saturation must be generic for any EOS model that
satisfies its declared domain and numerical requirements.

## 5. Association model

### 5.1 Topology

Association topology is one fixed graph of:

- component/site identities;
- positive multiplicities;
- advisory `donor`, `acceptor`, or `effective` site roles; and
- active site-pair edges.

Named schemes such as 1A, 2B, 3B, and 4C expand into this graph before EOS
evaluation. They do not select implementations. The current physically inert
`site_class` field is replaced by `site_role`.

Roles validate presets and declared combining rules but do not create edges,
infer chemistry, or globally prohibit a cited effective same-role or self-pair
model. The explicit graph is authoritative.

### 5.2 Edge sources and combining rules

Every active edge has exactly one source for association energy and volume:

1. explicit values; or
2. one named, cited Association Combining Rule.

Preflight rejects missing sources, competing sources, precedence rules, and
project-wide automatic defaults. It materializes each derived cross value and
its exact derivative chain. FitResult retains the rule identity, inputs,
resolved values, and derivative mapping.

Induced association is represented by the resulting cross edges while the
specified self edges remain inactive. It is not a specialized EOS or
Regression route.

### 5.3 Derivative scope

Selected association energy and volume coefficients are fit-ready only when
their exact derivatives work for every supported fixed-topology EOS context
that uses them, including:

- pure and mixture models;
- neutral and electrolyte models;
- constant and temperature-dependent segment diameter; and
- arbitrary admitted site graphs, multiplicities, subsets, and ties.

The current neutral/constant-diameter derivative gate is therefore a deletion
target, not a public limitation to preserve.

## 6. Polar and permittivity policy

### 6.1 Molecular Moments

Permanent dipole and quadrupole moments are fixed measured or
quantum-calculated inputs. A polar component may still fit ordinary base,
pair, or association coefficients while its moment remains fixed. Fitting an
unlabeled effective moment is forbidden. A future effective-moment model would
need a separately named physical meaning, dataset, and decision.

### 6.2 Permittivity Submodel

Relative static permittivity is supplied by one explicitly selected auxiliary
constitutive model. EOS evaluates it, uses it in electrolyte contributions,
and may report the bulk value. The model's coefficients live in a typed Model
Family Namespace in the same Parameters document.

The Staged Calibration Baseline is the default workflow:

1. calibrate the selected permittivity submodel against direct permittivity
   observations;
2. hold those coefficients fixed while fitting downstream EOS coefficients;
3. permit explicit re-selection or a Joint Calibration when scientifically
   justified.

A Joint Calibration may combine direct permittivity and downstream EOS
observations, but they remain separate residual blocks with explicit block
weights and family-specific diagnostics. No base parameter fit activates
permittivity coefficients implicitly.

## 7. Observation and data model

### 7.1 Observation Families

An Observation Family declares:

- the target thermodynamic quantity;
- state inputs, components, composition, and reference convention;
- units and basis;
- valid EOS model domain;
- residual representation;
- required exact derivatives; and
- structured failure meanings.

First-cutover families may use EOS-owned fixed-state values, density closure,
reference-dependent electrolyte quantities with explicit conventions, and
Pure Saturation State values. Any family that requires a mixture phase split,
flash, or reaction solve is unavailable.

The first mandatory vertical slice is the Pure Parameter Baseline:

- saturation pressure; and
- saturated-liquid density or molar volume;

over an explicit temperature domain. Saturation-pressure-only fits remain
valid when explicitly requested but are labeled reduced-information fits.
Saturated-vapor density or volume may be added when the source reports it.

### 7.2 Source representations

A Fit Data Source is either:

1. typed observation rows; or
2. a Published Observation Correlation.

Both require citation, source locator or digest, units, basis, validity domain,
and transformation record. Files are locators, not scientific identities.
Regression contains no packaged paper tables and discovers no Data packet by
name.

A Published Observation Correlation belongs to a closed registry and declares
its exact equation, coefficients, independent/output transformations, units,
domain, and derivatives. “Polynomial” alone is invalid because it does not
identify whether the basis is `T`, `1/T`, reduced temperature, pressure, or
log-pressure. DIPPR, Antoine, Wagner, and an explicitly defined power-series
form may be admitted when active project data requires them.

For a continuous source correlation, preflight increases deterministic
Gaussian quadrature order until the normalized objective and exact Jacobian
converge, then freezes nodes and weights for Ceres. FitResult retains the
correlation identity, domain, measure, quadrature, convergence evidence,
predictions, residuals, and Jacobian rows. Explicit tabulated points remain the
right representation for direct measurements and selected operating states.

## 8. Residual and objective contract

- Positive quantities use a natural-log ratio residual by default:
  `log(predicted / observed)`.
- Signed quantities use a declared scaled difference:
  `(predicted - observed) / scale`.
- A reported row standard deviation may scale a residual only when its meaning
  and confidence conversion are explicit.
- Every scaled residual has equal weight by default.
- An observation block may declare an explicit multiplier; block identity and
  contribution remain visible in FitResult.
- The first cutover uses squared loss only.
- It has no robust loss, prior, regularization, hidden stabilization residual,
  covariance claim, or uncertainty interpretation.
- Missing uncertainty does not invent statistical uncertainty.
- Custom residual callbacks and executable expressions are forbidden.

The historical Gross--Sadowski squared relative-error objective is retained as
source evidence. The generic default uses log ratio for positive quantities
because it is symmetric with respect to reciprocal relative over- and
underprediction; a source-faithful reproduction may explicitly request its
published residual when that exact evidence class is required.

## 9. Exact derivative contract

### 9.1 ActiveParameterSet

EOS owns `ActiveParameterSet` beside `Parameters`. Every public EOS evaluation
used by Regression accepts the same ordered active set and returns values plus
exact selected-parameter derivatives. Regression never maps raw stored slots,
reimplements equations, or decides derivative support from record shape.

### 9.2 CppAD authority

CppAD remains the sole nonlinear production derivative authority. Production
finite differences, handwritten nonlinear derivatives, and differentiation of
solver iteration history are forbidden.

For solved states, including density closure, association mass action, and
pure saturation, EOS returns total derivatives obtained by exact implicit
linear algebra from CppAD partial tensors of the defining equations.

Every capability record binds:

- parameter identities and order;
- Observation Family and model domain;
- units and basis;
- value and derivative shapes;
- active topology and parameter fingerprints; and
- structured invalid-input, model-domain, numerical, and physical failures.

### 9.3 Derivative verification

Independent numerical or analytic checks may verify CppAD derivatives in
Engine evidence, but they are not runtime fallbacks. Verification must cover
representative base, pair, association, correlation, polar-with-fixed-moment,
and auxiliary-permittivity contexts before those families are admitted.

## 10. Optimizer and execution

### 10.1 Coordinates and bounds

- Ceres remains the optimizer.
- Every optimizer coordinate is an affine-scaled physical parameter.
- Bounds, scales, and final values remain expressible in physical units.
- No hidden log, logit, sign, or chemistry-derived transform is permitted.
- Discrete model choices and topology never enter Ceres.

### 10.2 Problem graph

- Use one Ceres parameter block per independent optimizer coordinate.
- Use one residual block per observation state or correlation quadrature node;
  vector observations may share one EOS evaluation.
- Cache only immutable resolved model structure and within-evaluation shared
  state; never cache a mutable fitted result across parameter points.
- Start with `DENSE_QR`. Adopt sparse linear algebra only after a retained
  benchmark demonstrates a material benefit for the actual block graph.

### 10.3 Multistart

- Use a declared deterministic low-discrepancy multistart design over physical
  bounds, optionally including the input parameter vector.
- Each Ceres solve is single-threaded.
- Independent starts may run in deterministic parallel workers.
- Retain every basin result, including rejected, non-evaluable, bound-active,
  rank-deficient, and higher-cost outcomes.
- Do not report the lowest objective as unique or globally optimal.

### 10.4 Non-evaluable trials

When an EOS or Pure Saturation State cannot return its declared value and
exact derivative at a bounded trial point, the residual evaluation returns a
structured non-evaluable status. It does not fabricate a residual, clip a
state, or substitute a fixed penalty. Ceres may reduce or reject the step. A
failed start remains in FitResult.

## 11. FitResult

`FitResult` is immutable and contains enough evidence to understand the fit
without rerunning it:

- Resolved Fit schema version and digest;
- input Parameters fingerprint and complete active identity expansion;
- ties, bounds, scales, starts, solver options, and deterministic start order;
- source and observation block identities;
- final fitted physical and affine coordinates;
- predictions, residuals, block contributions, and row accounting;
- complete exact Jacobian with row and column identities;
- singular values, rank, conditioning, and active-bound diagnostics;
- local parameter correlation and sensitivity diagnostics without labeling
  them as statistical uncertainty;
- physical and domain acceptance for each observation;
- all multistart basin outcomes;
- resolved association combining-rule values and derivative chains;
- source-correlation quadrature evidence where applicable; and
- artifact, package, and parameter fingerprints.

`FitResult.to_parameters()` creates the fitted Parameters value. Export to JSON
or another file is an explicit caller action. Regression does not mutate Data,
Validation, source packets, or the input parameter file.

Covariance matrices, standard errors, confidence intervals, bootstrap,
profile likelihood, cross-validation, scientific acceptance, and predictive
claims are not inferred by Engine. Validation owns those campaigns when a
concrete scientific question requires them.

## 12. Canonical schema shape

The exact JSON Schema is generated from the same closed registry used by the
runtime parser; it is not maintained independently and does not use Pydantic.
The following is the normative shape, with final field spelling frozen during
Phase 4:

```json
{
  "schema_version": 1,
  "parameters": {
    "path": "parameters.json",
    "sha256": "..."
  },
  "model": {
    "kind": "pure"
  },
  "active": {
    "selectors": [
      {
        "component": "methane",
        "families": [
          "segment_count",
          "segment_diameter",
          "dispersion_energy_over_k"
        ]
      }
    ],
    "coordinates": {
      "component/methane/segment_count": {
        "bounds": [0.8, 1.3],
        "scale": 0.1
      },
      "component/methane/segment_diameter/constant/value": {
        "bounds": [3.2, 4.2],
        "scale": 0.1,
        "unit": "angstrom"
      },
      "component/methane/dispersion_energy_over_k": {
        "bounds": [100.0, 220.0],
        "scale": 10.0,
        "unit": "kelvin"
      }
    },
    "ties": []
  },
  "observations": [
    {
      "family": "pure_saturation",
      "source": {
        "kind": "table",
        "path": "methane-saturation.csv",
        "sha256": "...",
        "citation": "...",
        "domain": {"temperature_k": [110.0, 180.0]}
      },
      "targets": [
        "saturation_pressure",
        "saturated_liquid_mass_density"
      ],
      "residuals": {
        "saturation_pressure": "log_ratio",
        "saturated_liquid_mass_density": "log_ratio"
      },
      "multiplier": 1.0
    }
  ],
  "solver": {
    "loss": "squared",
    "multistart": {
      "kind": "low_discrepancy",
      "count": 16,
      "include_input": true
    },
    "linear_solver": "dense_qr"
  }
}
```

Selector expansion never supplies bounds or scales. The `coordinates` mapping
must cover every resolved independent identity after ties. This keeps selector
convenience separate from the complete numerical contract.

## 13. Target module shape

The replacement should be a deep module with a small public interface and a
small number of private owners:

```text
engine/src/epcsaft/regression/
  __init__.py              # ParameterFit, FitResult, fit
  _model.py                # immutable schema records and strict loaders
  _resolve.py              # selectors, ties, sources, capabilities, preflight
  _result.py               # immutable result and explicit export
  native/
    fit.hpp
    fit.cpp                # generic residual blocks and EOS calls
    ceres_core.hpp
    ceres_core.cpp         # multistart, bounds, solve, matrix diagnostics
    module.cpp             # generic binding only
```

The exact private split may change if implementation locality improves, but
the public interface and ownership may not expand. There is no runtime
`regression/data/` directory and no source-specific native file.

## 14. Current-to-target deletion map

The following live surfaces are displaced, not wrapped:

| Current surface | Final action |
|---|---|
| `parameter_regression.py` public coordinate, observation, capability, and problem classes | Replace with strict private schema/resolution and the three-export public interface. |
| `usability.py` public preparation/preflight object graph | Fold into `ParameterFit.resolve()`/internal preflight; expose diagnostics through FitResult or a resolved read-only view only if a real caller needs it. |
| `workflow.py` source-specific fit functions and result classes | Delete. |
| `records.py` Figiel/Born/water-factor specifications and loaders | Delete. |
| `native/born_diameter_fit.*` | Delete after generic component-parameter fitting proves the same physical capability. |
| `native/figiel_kij_fit.*` | Delete after generic pair-parameter fitting proves the Figiel validation case. |
| `native/figiel_water_factor_fit.*` | Delete after generic component/correlation fitting proves the same capability. |
| `native/general_fit.*` 4,000-plus-line observation switch | Replace with registry-resolved generic blocks and EOS capability calls; do not layer a second dispatcher around it. |
| packaged source CSV files under `epcsaft.regression` | Remove; canonical reusable inputs belong in Data and literature campaigns in Validation. |
| source-named Engine Regression tests | Replace only with source-neutral numerical invariants; move literature breadth to Validation. |
| source-specific Regression science docs | Retain only as clearly marked provenance or move their live conclusions into this plan/current docs. |
| current quickstart | Rewrite completely against `ParameterFit`, `FitResult`, and `fit` after cutover. |

### Documentation migration matrix

The implementation branch must classify and update documentation as follows:

| Document or group | Cutover treatment |
|---|---|
| `engine/README.md` | Keep the canonical plan link now; replace pre-cutover usage with the three-export interface at Phase 10. |
| `engine/CONTEXT.md` | Keep the vocabulary and plan link authoritative; revise only if implementation reveals an approved vocabulary change. |
| `engine/ARCHITECTURE.yaml` | Replace `workflow.py` and current native-owner entries with the final generic owners and record the EOS-only Regression boundary. |
| `engine/docs/python-api.md` | Document `Pure`, strict `Mixture`, `ParameterFit`, `FitResult`, `fit`, and explicit export; remove current public Regression types. |
| `engine/docs/regression/general-regression-quickstart.md` | Rewrite from scratch using one pure baseline and one generic pair/correlation example; show JSON and mapping equivalence. |
| `engine/docs/regression/science/general-parameter-regression.md` | Preserve as marked pre-cutover implementation evidence or archive it; no installed-SDK, source-specific, or current-capability claim may remain live. |
| `engine/docs/regression/science/pure-saturation-regression.md` | Move literature data/results to Validation-facing evidence and replace live owner statements with `Pure.saturation` plus the generic pure baseline. |
| `engine/docs/regression/science/figiel-aqueous-kij-regression.md` and `figiel-born-diameter-tracer.md` | Retain only as source campaign evidence; replace any live specialized interface instruction with the generic parameter/observation identities. |
| `engine/docs/regression/science/literature-reproduction-contract.md` | Keep its evidence-class distinctions; update only public interface and artifact-path assumptions. |
| `engine/docs/regression/science/neutral-hydrocarbon-next-slice.md` | Recast any live implementation instruction as generic pair-parameter Validation or archive it if fully completed. |
| `engine/docs/regression/science/mea-coupled-regression-master-plan.md` | Mark blocked on the separate Equilibrium redesign; do not pull its reactive coordinates into this cutover. |
| `engine/docs/regression/research/**` | Preserve paper/source findings as provenance; update only explicit broken live links, never historical scientific wording. |
| `engine/docs/regression/superpowers/plans/**` | Preserve as historical plans and add no normalization edits. |
| `validation/campaigns/2025-figiel/README.md` and affected Regression campaign READMEs | Replace specialized Engine commands with installed-wheel generic fit inputs while preserving source names and results as campaign provenance. |
| reusable Data packet READMEs and manifests used by the migrated fits | Document the source rows/correlations and materialization contract, with no dependency on an Engine workflow name. |
| Governance plans, decisions, manifests, and accepted receipts | Preserve historical content. Update live topology pointers only if the implementation changes a governed owner; this Regression redesign alone does not require a doctrine revision. |

Pre-monorepo campaign lock files and general EOS handoff runners are not
compatibility targets or phase gates for this Regression rewrite. Preserve
completed receipts as provenance, but do not retrofit old Khudaida or EOS
handoff infrastructure merely to keep it runnable. Phase 9 builds the smallest
new installed-wheel campaigns needed to prove the generic Regression slices.
The broader Khudaida phase-equilibrium campaign may be rewritten separately
after this cutover and does not gate the Regression public API replacement.

Phase 0's inventory expands this matrix to every concrete reference before any
public deletion. A historical document is not rewritten merely to remove an old
name; a live instruction is not allowed to remain stale merely because it is
old.

## 15. Ordered implementation sequence

### Phase 0 — Freeze the replacement contract

1. Make this plan and `engine/CONTEXT.md` the only live design authority.
2. Inventory every public export, native binding, packaged dataset, test,
   Validation campaign, and documentation page that references the current
   Regression interface.
3. Classify each as retained evidence, migrated input, replacement target, or
   deletion target.
4. Freeze the first-cutover Observation Families and parameter families from
   the actual active data, not hypothetical future use.

Exit gate: complete machine-readable inventory with no unclassified current
surface.

### Phase 1 — Canonical Parameters and topology

1. Define the strict one-document schema and Model Family Namespaces.
2. Give every continuous coefficient a source-neutral Parameter Identity.
3. Separate fixed specification fields from selectable coordinates.
4. Replace `site_class` with advisory `site_role`.
5. Add named association-preset expansion and explicit edge-source records.
6. Add the closed Association Combining Rule registry.
7. Rename source-named live coefficient families by physical role.
8. Migrate all active parameter documents atomically; add no old-schema reader.

Exit gate: every active EOS model resolves from the new schema and the old
schema is absent from live inputs.

### Phase 2 — EOS object and pure-saturation ownership

1. Introduce strict `Pure` and strict `Mixture` construction.
2. Move generic pure saturation into EOS ownership as `Pure.saturation`.
3. Remove fingerprint-approved component routing.
4. Return saturation pressure and both coexisting densities with structured
   numerical and physical diagnostics.
5. Keep all mixture phase-equilibrium and reaction operations in Equilibrium,
   unavailable to Regression.

Exit gate: arbitrary admitted pure models solve and differentiate saturation
without an Equilibrium import or approved-fingerprint table.

### Phase 3 — ActiveParameterSet and exact derivatives

1. Add the EOS-owned ActiveParameterSet and capability records.
2. Make base, pair, correlation, association, Born, polar-with-fixed-moment,
   and selected permittivity coefficients derivative-ready where admitted.
3. Generalize association derivatives across the required EOS contexts.
4. Implement exact implicit sensitivities for density closure, association,
   and Pure Saturation State.
5. Remove regression-facing native SDK/callback concepts and raw slot matching.

Exit gate: a source-neutral derivative matrix passes independent numerical or
analytic checks for every admitted family and unsupported families fail before
optimization.

### Phase 4 — ParameterFit schema, sources, and resolution

1. Implement the strict mapping and JSON loaders from one registry.
2. Generate JSON Schema from that registry.
3. Implement exact identities, typed selectors, ties, deterministic expansion,
   and Resolved Fit hashing.
4. Implement row and Published Observation Correlation sources.
5. Implement units, basis, domain, reference convention, and source checks.
6. Implement pure-saturation integrated-correlation quadrature.
7. Reject unknown keys, unknown families, unsupported correlations, incomplete
   controls, and executable callbacks.

Exit gate: JSON and mapping inputs resolve byte-for-byte-equivalent fit records
and all scientific/numerical choices are explicit.

### Phase 5 — Generic residual engine and Ceres execution

1. Build one generic native residual owner around EOS capabilities.
2. Implement positive log-ratio and signed scaled-difference residuals.
3. Implement separate observation blocks and explicit block multipliers.
4. Implement one-coordinate parameter blocks and one-observation residual
   blocks.
5. Implement deterministic low-discrepancy multistart, single-threaded solves,
   and deterministic parallel start scheduling.
6. Implement non-evaluable trial backtracking without penalty residuals.
7. Implement exact Jacobian rank, singular-value, conditioning, active-bound,
   and basin diagnostics.

Exit gate: one native path executes every admitted EOS-only Observation Family
without family-specific Ceres code.

### Phase 6 — FitResult and transfer

1. Implement immutable FitResult with the complete evidence contract.
2. Implement `to_parameters()` as a new immutable value.
3. Implement explicit JSON export and canonical hashing.
4. Remove uncertainty language not supported by an uncertainty model.
5. Preserve failed and competing starts instead of summarizing only the winner.

Exit gate: the result can be independently audited and transferred without a
mutable parameter catalog or filesystem side effect.

### Phase 7 — Migrate the required scientific slices

Migrate through the generic path in this order:

1. nonassociating Pure Parameter Baseline: `m`, `sigma`, and `epsilon/k` from
   saturation pressure plus saturated-liquid density or volume;
2. associating pure components: selected base, association-energy, and
   association-volume coefficients on generic preset-expanded topology;
3. fixed-state mixture pair coefficients, including constant and admitted
   correlated `k_ij`, with no paper-specific route;
4. direct permittivity-submodel calibration and explicit staged/joint use;
5. EOS-owned electrolyte observables with explicit reference conventions,
   including the generic capabilities previously demonstrated by Born,
   water-factor, and aqueous-`k_ij` source workflows; and
6. polar components using fixed Molecular Moments and selected ordinary model
   coefficients.

Each slice must use the same public interface, schema, native solver, result,
and transfer path. A slice that needs a new Observation Family or parameter
family extends the closed registries, not the public interface.

Exit gate: all retained current scientific capabilities execute through the
generic path or are explicitly removed as unsupported.

### Phase 8 — Delete displaced implementation

1. Remove every source-specific Python export, result type, loader, packaged
   dataset, and native fit file.
2. Remove the old public coordinate/observation object graph.
3. Remove old Regression-facing SDK/header remnants and capability callbacks.
4. Remove contract, compatibility, serialization-shape, prose, and
   source-specific Engine tests.
5. Remove obsolete CMake sources and public documentation.
6. Run source-name, legacy-path, duplicate-owner, and dead-export audits.

Exit gate: no old route is importable, compiled, documented as current, or
needed by a retained test.

### Phase 9 — Scientific evidence and installed-wheel Validation

1. Retain compact Engine evidence for exact derivatives, implicit
   sensitivities, association graph/multiplicity behavior, saturation
   coexistence, and non-evaluable trial semantics.
2. Run the Gross (2001) pure baseline as an installed-wheel Validation
   campaign using paired saturation pressure and saturated-liquid density or
   volume.
3. Run the Gross (2002) association case through generic topology and
   association coordinates.
4. Run the Figiel pair-parameter case through generic `k_ij` selection and
   compare the recovered physical value to the cited value using the
   campaign's engineering tolerance; the previously approved `1e-1` parameter
   tolerance belongs in Validation, not an Engine unit test.
5. Run direct permittivity and electrolyte downstream blocks separately and,
   where useful, as one explicit Joint Calibration.
6. Report parameter rank, conditioning, active bounds, basin agreement,
   residual structure, source domain, and held-out or profile results without
   claiming global uniqueness.

Exit gate: installed-wheel campaigns prove the retained literature workflows
without importing source checkout internals or reviving source-specific APIs.

### Phase 10 — Atomic public cutover

1. Export only `ParameterFit`, `FitResult`, and `fit` from
   `epcsaft.regression`.
2. Rewrite the quickstart and Python interface docs.
3. Update `ARCHITECTURE.yaml`, Engine README, equation ownership navigation,
   and root affected-check routing.
4. Build one non-editable wheel and verify module origins.
5. Run focused Engine evidence, all retained Regression tests, applicable
   Validation campaigns, root affected checks, full checks, preflight, and the
   cleanup audit.
6. Review the complete deletion diff before any commit or push.

Exit gate: one public route, one native owner, one parameter artifact, one
wheel, no compatibility layer, and no unresolved current documentation.

## 16. Evidence allocation

| Claim or failure | Owner |
|---|---|
| CppAD selected-parameter derivative and implicit-sensitivity correctness | Compact parameterized Engine numerical evidence |
| Association graph, multiplicity, preset expansion, explicit edge source, and combining-rule chain rule | Compact Engine scientific evidence |
| Pure saturation pressure/fugacity equality and physical phase separation | Compact Engine scientific evidence |
| Bounds, rank, conditioning, basin retention, and structured non-evaluable trials | Compact Engine numerical evidence |
| Strict schema, manifests, hashes, forbidden names, and generated JSON Schema freshness | Readiness/preflight commands, not scientific tests |
| Gross, Figiel, permittivity, electrolyte, and broader component/property agreement | Installed-wheel Validation campaigns |
| Parameter uncertainty, profile likelihood, bootstrap, cross-validation, and predictive acceptance | Validation only when explicitly designed |
| Application-specific acceptance for MEA or lithium extraction | Their future end-to-end vertical slices after Engine and Equilibrium readiness |

Engine must not retain tests for imports, legacy aliases, exact exception text,
serialization snapshots, source wording, callback ABI shape, or source-specific
workflow names. Every persistent numerical test must identify a distinct
scientific failure.

## 17. Cutover acceptance checklist

The implementation is complete only when all items are true:

- [ ] `epcsaft.regression` exports only `ParameterFit`, `FitResult`, and `fit`.
- [ ] JSON and mapping loaders share one strict schema and behavior.
- [ ] One Parameters document owns base and auxiliary model inputs without a
      flat untyped map.
- [ ] Every active coefficient has one source-neutral Parameter Identity.
- [ ] Typed selectors expand deterministically and FitResult freezes the exact
      expansion.
- [ ] Pure and Mixture have strict cardinality and share State.
- [ ] Pure saturation is EOS-owned, generic, and exactly differentiated.
- [ ] No first-cutover Regression path imports Equilibrium or fits reaction
      constants.
- [ ] Molecular Moments are fixed inputs.
- [ ] Permittivity is an explicit auxiliary model with staged default and
      explicit Joint Calibration.
- [ ] Association presets expand to one graph; site roles are advisory; every
      edge has one explicit value or cited rule source.
- [ ] Association derivatives cover every supported fixed-topology EOS context.
- [ ] Any coefficient of an admitted typed parameter correlation can be active.
- [ ] Tabulated and published-correlation observation sources are both
      supported without executable callbacks.
- [ ] Positive and signed residual semantics are explicit and source-neutral.
- [ ] Ceres uses affine physical coordinates, physical bounds, exact
      derivatives, deterministic multistart, and retained basin outcomes.
- [ ] Non-evaluable trials never receive fabricated penalties.
- [ ] FitResult exposes full numerical and identifiability diagnostics without
      unsupported uncertainty claims.
- [ ] `to_parameters()` returns a new immutable Parameters value.
- [ ] No runtime file, symbol, public test, or live instruction is named for a
      paper, author, campaign, or application.
- [ ] No specialized Born, water-factor, Figiel, pure-saturation, or callback
      Regression implementation remains.
- [ ] No packaged literature data remains in `epcsaft.regression`.
- [ ] Engine tests contain only retained numerical/scientific claims and
      literature breadth runs from the installed wheel in Validation.
- [ ] One non-editable wheel passes focused, affected, full, preflight, and
      cleanup checks.

## 18. Explicitly deferred work

The following are not hidden TODOs inside the first cutover:

- redesigning Equilibrium and admitting its VLE, LLE, flash, or reactive
  capability records;
- fitting reaction-equilibrium constants or application chemistry;
- total caloric properties before an EOS-owned ideal-gas thermochemistry model
  exists;
- robust losses, priors, regularization, covariance, confidence intervals,
  Bayesian inference, bootstrap, profile likelihood, or experimental design;
- arbitrary user expressions, Python callbacks, C++ plugins, or mutable
  runtime registration;
- ordinary fitting of Molecular Moments;
- a homologous-series parameter model such as the Gross (2001) Appendix-B
  n-alkane relation; and
- sparse Ceres linear algebra without an actual benchmark demonstrating need.

Each deferred capability requires its own scientific question, owner,
interface, derivatives, data contract, and acceptance evidence. None may enter
through a compatibility helper or source-specific shortcut.

## 19. Principal risks and controls

| Risk | Control |
|---|---|
| Generic schema becomes an untyped expression language | Closed registries, strict unknown-key rejection, no callbacks, generated schema from runtime definitions |
| “Any parameter” becomes raw record-slot fitting | Canonical Parameter Identities and EOS-owned capability records |
| Saturation and association derivatives silently remain narrow | Explicit all-context derivative gates and independent numerical checks |
| Joint fits hide weak identifiability | Separate residual blocks, full Jacobian diagnostics, rank/conditioning, active bounds, and all basin outcomes |
| Source-specific behavior survives under generic wrappers | Replace-don't-layer implementation and final forbidden-name/dead-owner audits |
| Temporary dual paths become permanent | Private incremental internals followed by one atomic public cutover and same-branch deletion |
| Published correlations are used outside their meaning | Required equation family, units, transformations, validity domain, citation, and quadrature receipt |
| Equilibrium instability contaminates Regression | EOS-only first boundary and explicit future Equilibrium admission gate |
| Objective reduction is mistaken for scientific validity | Engine reports numerical/physical/identifiability evidence; Validation adjudicates scientific and predictive acceptance |

## 20. First execution action

When implementation is authorized, begin with Phase 0's machine-readable
inventory and Phase 1's canonical Parameters schema. Do not begin by wrapping
the current source-specific functions or extending `general_fit.cpp`. The
first implementation diff must establish the replacement model and its
deletion map so later phases cannot accidentally preserve the current public
surface.

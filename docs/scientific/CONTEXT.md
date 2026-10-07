# Scientific context

repository role: tool

## Agent role

agent role: chemical engineer specializing in electrolyte equations of state, phase and reaction equilibrium, and numerical thermodynamics

Maintain the Engine’s source-backed thermodynamic equations, shared native calculations, exact derivatives, equilibrium and fitting interfaces, and numerical qualification. Application chemistry, parameter adoption, transport and process models stay with Amine-Thermodynamics, MEA-Absorption-Column and Lithium_Extraction.

## Purpose and first acceptance

Build a general ePC-SAFT fluid thermodynamics platform: single-state properties,
phase and reaction equilibrium, exact derivative actions, and parameter fitting.
The three official applications are its first required acceptance uses, not
special cases encoded into the Engine. Generality means admitting physical
specifications and source-backed models within documented domains; it does not
promise every fluid, specification or derivative order before evidence exists.

Amine-Thermodynamics first requires reactive properties and fitting;
MEA-Absorption-Column requires consistent total caloric properties and the
solved-state actions its process equations consume; Lithium_Extraction requires
reactive liquid-liquid equilibrium and extraction observations. Their chemistry,
parameter adoption, transport and process models belong to the applications.
Solids remain a separately justified extension; kinetics and process simulation
are not current Engine responsibilities.

## The project map

```mermaid
flowchart TD
  M[Mathematics and algorithm requirements] --> E[Engine: shared native calculations]
  D[Data: explicit scientific inputs] --> E
  E --> W[One installed wheel]
  W --> A[Analyses: numerical evidence and interpretation]
  W --> U[Applications: chemistry and process calculations]
  A --> R[Review agreement and limitations]
  U --> R
  R --> M
```

The review arrow means a discrepancy prompts diagnosis or a justified proposed
scientific revision. Code behavior never silently changes the mathematics.
GitHub leaf issues record missing and unqualified capability separately from
this intended architecture.

| Owner | Responsibility |
|---|---|
| `docs/scientific/` | This map, accepted decisions, review evidence and writing guidance |
| `engine/docs/` | Canonical equations, conventions, formulation and algorithms |
| `data/` | Reusable scientific inputs, provenance and explicit materialization |
| `engine/native/` | Shared scientific calculations and numerical solves |
| `engine/src/epcsaft/` | Python input/output mapping to the native implementation |
| `engine/tests/` | Compact numerical falsifiers of Engine claims |
| `analyses/` | One-question calculations, retained results and interpretation using the installed wheel |
| Official application repositories | Application chemistry, parameter adoption, process equations and application acceptance |

The Pyomo/IDAES adapter ships `libepcsaft_asl.so` in the same non-editable Engine
wheel. Its neutral functions transport shared Helmholtz pressure derivatives;
the [neutral analysis](../../analyses/2026-idaes-asl-adapter/README.md) owns that
qualification. Its equilibrium function transports E5 phase-component amounts
and exact first/second feed actions at fixed T/P on a declared two-liquid branch.
The [Khudaida analysis](../../analyses/2026-khudaida-idaes-settler/README.md) owns
the nonreactive qualification. An explicit affine feed map, supports, scalar
reaction targets and phase-identity criterion are caller inputs, not application
selectors. Each downstream application must qualify its own chemistry, input
chart, nonlinear process chain and operating domain; neither reference qualifies
reactive Lithium by implication.

## Engine responsibilities

The [native conversion evidence](../../analyses/README.md#native-conversion-evidence)
records the completed mapped conversions and retirements, changed comparison
claims, and the remaining Wangler source-reference calculation.
The layout check's absence of `scipy.optimize` is a lexical result, not a physical validation claim.

Native dependencies run downward in `ARCHITECTURE.yaml`; the root
[native calculation rule](../../AGENTS.md#native-calculation-rule) owns execution.

| Layer | Owns | Boundary |
|---|---|---|
| core | Common constants, parameter coordinates, numerical policy and diagnostics | Small shared foundations; no phase or application selection |
| Contributions | Ideal-gas Helmholtz formulation and individual residual contributions, including association and electrolyte terms | Own thermal correlations and ideal mixing as physical inputs; preserve all derivative dependencies |
| Properties | Ideal, Residual and Total quantities derived from Contributions; state evaluation, density closure and reference conversions | One native calculation; governing definitions determine assembly and derivative chains |
| equilibrium | Specification compilation, shared rows, one equilibrium driver, continuation, local checks, explicit tangent-plane testing/discovery and solved-state actions (qualified cases in the retained analysis) | All equilibrium specifications use the same driver; no independent bubble, dew, reactive or paper-specific solvers |
| `regression` | Observation composition, weighted fitting and fit diagnostics | Consumes model/equilibrium values and actions; does not implement them again |
| surface | Physical input resolution and admission | Explicit parameters, units, choices and applicability |
| bindings | The Python/native boundary | One extension; follow the [native calculation rule](../../AGENTS.md#native-calculation-rule) |

Properties `ReferenceTaylor` means ordinary reference-contraction and volume
Taylor coefficients along a fixed-composition, pressure-closed base path.
`reference_partial` supplies raw source-pressure/T/one-coordinate partials;
`NeutralReferenceReactionDiagnostic::Derivative` holds each consumed C3/C4
column and both terminal changes. Their meanings and bounded inventory follow
[the single-insertion equation](../../engine/docs/equations.md#reference_single_insertion_taylor)
and [F8/F9 notebook](../../analyses/2026-reactive-caloric-actions/notebook.qmd).
These standalone Properties derivatives and all five lower reference partials
are numerically verified on the stated manufactured cases. The current-source
replay and its distinct wheel identity are retained in that notebook; Properties
reuses the resolved terminal States through one single-insertion Taylor owner.
The same diagnostic retains actual lower-column evaluation/admission causes
for chemical rows and caloric completion; the notebook retains the bounded
corrected-source failure case separately from the admitted-case replay.
The [current runtime review](../../analyses/2026-reactive-caloric-actions/results/current-source-review.md)
is Supported at `14c32d1a`, accepted for the bounded integrated reference and
prescribed-T/P caloric slice, without manuscript promotion. The same notebook
retains separate driven-T/P numerical
checks of anchored Gibbs-volume calorics, A1/A2 with one parameter-bearing
direction and feed-only A3. Referenced free-intensive combinations and
apparent-Cp fitting remain separate #225 qualification. The subsequent source
implementation refreshes referenced free-T/P row jets at their current point
and transports coefficient increments through the native law at T(t), including
free-T addition endpoints. The subsequent source admits finite Gibbs-conjugate
volume with finite nonzero pressure slope, recording stable-chart signs separately.
The temperature-scaling correction received Supported delivered re-review at
`e4342588`. The current-main replay at `2428781b`, integrating `d7739584`, retains
that correction, the density publication behavior, native duplication removals,
trace payload construction inside `Scope::observe` and the separate exact-inventory
incipient extension. The notebook reads its retained native values, action
comparisons, reference refinements, lower-admission ladder, binary128 pin,
independent moved-row evidence and paired work/timing. All shared movements from
the common native recording have independent comparisons before refresh.
This is numerical verification on the declared manufactured branches.
Independent delivered review of the current integration remains with the lead;
earlier scoped reviews do not establish its acceptance. The subsequent apparent-Cp
build uses Equilibrium extensive H jets and its shared action factorization,
including feed-normalization actions and exact fitting columns. The
[slice-four numerical record](../../analyses/2026-reactive-caloric-actions/results/slice4/provenance.json)
and the existing notebook retain independent H-value comparisons, both quadrature
identities, nonreacting reduction, rebuilt-prediction columns and synthetic recovery.
Delivered review of this build remains pending; it establishes no physical calorimetry.

### Residual, ideal and state responsibilities

Contributions owns the ideal-gas Helmholtz formulation and the individual residual
formulations: hard-chain, dispersion, association and applicable electrolyte terms.
Residual thermodynamics remains the scientific priority. Ideal thermal correlations
and mixing contain physical information; reference conventions set enthalpy/entropy
and other consistent offsets. Ideal thermodynamics is not a synonym for reference
conventions, and it does not need a separate parallel category.

Properties has three scientific groups:

- **Ideal**: quantities derived from the ideal formulation.
- **Residual**: quantities derived from the residual formulation.
- **Total**: assembled physical quantities, including experimental observables,
  where supported. Assembly follows each quantity's governing definition; it is
  not always the arithmetic sum of identically named Ideal and Residual outputs.

Preserve the connected differentiation chain and shared calculations across these
responsibilities. Grouping and naming do not establish implemented availability.

A caller supplies a parameter bundle once to a reusable mixture. Construction
validates and prepares species/site mappings and genuinely state-independent
quantities. Temperature-dependent parameters, association fractions, density
roots and properties require the state inputs and are evaluated there. Shared
intermediates and useful computed results are reused within the state; retained
mixture preparation must remain valid across changing temperature and composition.

State requests specify temperature, composition and either density or pressure.
A pressure request resolves the declared density branch and reuses that resolved
state. Support both a full bundle and selected properties, computing only the
required dependencies and derivative orders. Full output means supported
properties with explicit unavailability for missing required inputs, not invented
ideal thermal data. Missing ideal calorics must not disable valid residual work.
Python exposes these capabilities; native owners perform the science and reuse.

Contributions and Properties are the layer names; `engine/native/contributions/`
and `engine/native/properties/` implement their dependency boundary under the
same names. Callable names, output
container shape, default evaluation policy and pure-fluid/mixture interface
unification remain bounded design decisions. The plan distinguishes the accepted
Ideal/Residual/Total organization from currently implemented outputs.

**One equilibrium driver is an owner requirement.** It owns the solve lifecycle
for all admitted equilibrium specifications and uses one row assembly and one
thermodynamic implementation. This is stronger than a facade dispatching to
independent equilibrium solvers. Density and association remain lower-level
numerical kernels. Continuation and phase discovery control sequences of
problems and reuse the same equilibrium driver where equilibrium is solved;
their search operations do not duplicate it.

The owner requires a single Ipopt driver with exact curvature, including bubble
and dew specifications. This does not certify the current implementation or
promise convergence for every input. Qualify its pinned/free-intensive equations,
derivatives and cost; report a demonstrated limitation for an owner decision
rather than substitute another equilibrium method.

Runtime accepts caller-supplied physical inputs and one documented unit system;
source-parameter validity ranges are the caller's responsibility. Mathematical
input validity and numerical diagnostics remain Engine responsibilities.
Phase counts/types/supports are declared unless the caller explicitly requests
`Discover()` at fixed T/P. Native tangent-plane testing and phase discovery are
implemented under [#223](https://github.com/tannerpolley/ePC-SAFT/issues/223),
with retained numerical qualification on the stated neutral and admitted
ionic/reacting cases; the independent delivered review supports those cases
and stated limits.
Only an explicit `stability_test(...)` call tests a declared result, including
pinned/free-intensive snapshots at solved T/P with frozen actual standards.
Declared solves, continuation, parameter fitting and ASL incur no implicit trial
solves. [ALG-STABILITY](../../engine/docs/science/algorithms.md#alg-stability)
owns the finite native search through the shared driver. Mechanical/local checks
and limits on global claims remain required. Caller tolerances govern solve
requests; available candidates that miss them retain warnings and achieved
residuals. The required state/property derivative target is directional actions
through third order on admitted smooth branches, with correct implicit chains.
Fitting admits residual-EOS parameters and named source-basis reaction constants
(scalar ln K or stored a/b/c/d), with exact first derivatives, including
finite-addition enthalpy actions for referenced endpoints at fixed temperature;
reference declarations stay fixed. Referenced free-temperature finite-addition
endpoints consume the same H owner at each solved temperature; driven-T/P
apparent isobaric heat-capacity fitting consumes mixed thermal H actions.
The [reaction-constant fitting notebook](../../analyses/2026-reaction-constant-fitting/notebook.qmd)
retains the manufactured first-derivative and recovery comparisons; it makes
no application or physical-validation claim. Concrete input forms and
outstanding qualification are recorded in the owning GitHub issues.
Declared neutral PH/PS specifications use extensive caloric balances and exact
controls/actions in the existing driver; UV remains deferred. Caloric properties/actions
consumed by MEA and the absorber remain current requirements. Fitting uses
density, total molar isobaric heat capacity with fixed thermochemistry,
finite-addition enthalpy with fixed ideal-gas records and declared reference
states, and phase observations. Public state units
follow SI molar conventions; composition actions use explicit zero-sum
directions. Owning GitHub issues record precise implementation inputs and
remaining work.

### Accepted cross-cutting scientific requirements

- Callers supply species, parameters and physical configuration. The Engine
  validates mathematics and reports numerical failure; it neither guesses inputs
  nor polices source-published fit ranges.
- Public state quantities use K, Pa, mol, m³, mol/m³, J/mol and J/(mol·K), with
  dimensionless mole fractions. Source-unit and reaction-convention conversions
  are explicit and never silently reinterpret parameter records.
- The single Ipopt driver uses exact curvature for admitted TP, bubble/dew and
  vapor-fraction specifications at fixed T or P. Phase fractions and absolute
  phase amounts are distinct. Neutral PH/PS uses the same driver; UV remains future work; explicit native
  fixed-T/P discovery has retained numerical qualification on the stated cases.
- Requests declare phase counts, types and supports or explicitly select fixed-T/P
  discovery. Discover at initialization, then reuse its declared topology where
  that reuse is qualified. The retained 100 mol neutral case has an admitted
  discovery but rejected declared reuse; [#243](https://github.com/tannerpolley/ePC-SAFT/issues/243)
  owns that conditioning limit. Native finite searches admit charge-neutral trials and conserved
  donor splits with estimated uncertainty; they supply no global certificate.
  Mechanical and local checks remain available; a declared result implies no
  unperformed search. [ADR 0012](decisions/0012-native-tangent-plane-phase-discovery.md)
  records the adopted choices, not numerical support.
- Charged trial supports in explicit testing/discovery require a declared bulk
  permittivity form with the vacuum limit: Schreckenberg `DensityTemperature`
  or Rueben `DipolarPerturbation` (full EOS §6.2/§6.4). Density-independent
  `IonSuppression` and `ComponentMixing` describe condensed-phase electrolyte
  terms and raise `ParameterError`: "permittivity model has no low-density limit;
  declare phases". Neutral-only trial families and declared calculations are
  unchanged. The retained solvent-only gas-density descent is refused under
  this owner-adopted search rule; its numerical operands remain available.
  A vacuum limit establishes neither physical validation nor global stability.
  No density cutoff or optimizer-success filter is introduced.
- Caller tolerances govern solve requests. Return finite candidates with achieved
  residuals and explicit unmet-tolerance warnings; return named absence and
  diagnostics when no finite candidate exists.
- Reaction conventions require explicit stoichiometry, standard concentration or
  molality, solvent basis and ionic gauge where applicable. Native `K(T)` uses
  `ln K = a + b/T + c ln(T/Tref) + dT` with explicit `Tref`; inconsistent supplied
  thermochemistry and `K(T)` remain an error rather than a hidden preference.
- Native calorics use caller-supplied references and supported polynomial,
  Shomate, NASA7, NASA9 and DIPPR107 forms. Preserve source-owned piecewise
  endpoints, reject gaps or ambiguous overlaps, never smooth discontinuities or
  integrate an anchor through a gap, and require inputs only for consuming outputs.
- General directional actions cover admitted smooth observables through third
  order in composition, density, temperature, pressure and active parameters.
  Composition directions are caller-supplied and zero-sum; singular actions do
  not erase valid central values or other available actions.
- Built-in regression observations cover density, total molar isobaric heat
  capacity at fixed T/P/composition and branch, finite-addition enthalpy between
  two declared equilibrium states with an explicit feed state, direct
  source-basis ln K at fixed T for a named reaction, saturation/bubble/dew
  pressure or temperature, phase compositions, phase fractions and
  phase-component amount ratios. Fit residual-EOS parameters and named reaction
  constants, retain failed-trial causes, include complete
  implicit prediction Jacobians, distinguish absolute uncertainty from relative
  weights, and condition covariance claims on rank, bounds and stated assumptions.

### Electrolyte fitting quantities

The shared active `ParameterRecord<Scalar>` includes Born, Debye-Huckel and
packing diameters (Angstrom) and ion suppression (dimensionless charged slots).
Properties owns the exact pressure-closed fugacity/permittivity columns,
terminal-reference admission and electrolyte contractions. Parameter regression
owns five observations, the dependent hydration chart, unweighted final checks
and local SVD/covariance diagnostics; [ALG-FIT](../../engine/docs/science/algorithms.md#alg-fit)
owns execution doctrine. `CoordinateTarget` means one physical family and its
component or reciprocal-pair indices with a proportional coefficient, defined
by [the affine parameter map](../../engine/docs/equations.md#affine_parameter_coordinate_map).
`HydrationHold` means a single-ion solvation equality
with a dependent Born diameter. `ObservationEvaluator` means the shared native
point evaluator with caller-selected rows and exact columns. `MinimaxRequest`
means a scaled epigraph estimation problem with declared max/mean guards;
`MinimaxResult` means its local point, feasibility, multipliers and KKT diagnostics.
`MinimaxRow` labels a target/max/slack-link/mean inequality in the units defined
by [ALG-FIT](../../engine/docs/science/algorithms.md#alg-fit). These use the same
thermodynamic, coordinate, check and candidate owners as least squares.

The [electrolyte parameter-gradient analysis](../../analyses/2026-electrolyte-parameter-gradients/README.md)
retains installed-wheel numerical verification at 298.15 K, 1 bar and the
declared Figiel solvent/salt selections. This is no physical validation,
#209 refit, packet adoption or mixed-salt qualification. Methanol and ethanol
are separate solvent selections; the packet has no ethanol-methanol pair.

## How the setup evolves

1. Identify the physical question, requested observable, domain and existing owner.
2. Apply the current leaf issue's reuse decision and adoption gate. Inspect
   current and superseded owners/evidence, preserve justified behavior, and
   separate unresolved scientific decisions from specified implementation.
   Resolve the method/coordinate/domain questions and smallest falsifier before
   dependent coding; independent in-scope mapping work may proceed.
3. Change the smallest complete scientific slice, including required documentation
   and numerical evidence; preserve the one-driver and ownership requirements.
4. Review agreement and limitations, then update the owning issue and retained
   evidence against the actual code and installed wheel. Application acceptance
   remains distinct.

New capability must earn its equations, inputs, ownership, evidence and budget.
The code floor is an aspirational engineering target; soft and hard ceilings
constrain the current scope. General-platform ambition does not authorize
speculative frameworks or automatic budget growth.

Superseded setup documents, stage briefs and agent execution instructions live
only in Git history. Keep current decisions and scientifically useful evidence
with their owners. Consult old records only for an identified provenance or
recovery question. Historical stage completion is never current acceptance.

## Authority

| Subject | Sole owner |
|---|---|
| Equation definitions, stable equation IDs and source locators | [equations.md](../../engine/docs/equations.md); generated `equations_registry.yaml` is an index |
| EOS symbols, units, conventions and source validity | [full-eos.md](../../engine/docs/science/full-eos.md) |
| Equilibrium formulation and specification admission | [formulation.md](../../engine/docs/equilibrium/formulation.md) |
| Algorithms, derivative orders, initialization, termination and failures | [algorithms.md](../../engine/docs/science/algorithms.md) |
| Parameter regression observation and fit meaning | [parameter regression formulation](../../engine/docs/regression/science/pure-saturation-regression.md) and [reproduction classes](../../engine/docs/regression/science/literature-reproduction-contract.md) |
| Current implementation plan, prerequisites and acceptance | [GitHub issue roadmap](https://github.com/tannerpolley/ePC-SAFT/milestones) |
| Current code owners and required scientific evidence | This scientific context, canonical Engine documents, tests, analyses and dated reviews |
| Project structure: components, Engine layers, path ownership and file roles | `ARCHITECTURE.yaml` (`project_map`), shown live in the companion Map view |
| Runtime implementation and executable topology | `engine/native/`, `engine/src/epcsaft/`, `ARCHITECTURE.yaml` |
| Accepted decisions | [decisions/](decisions/README.md), sole sequence |
| Reusable source packets and source applicability | `data/packets/`, `data/docs/packet-contract.md` |
| Numerical claims and retained evidence | `engine/tests/`, one-question sets under `analyses/`, dated reviews |
| Numerical policy values | `engine/native/core/numerical_policy.hpp`; algorithm requirements explain their meaning |
| Runtime source-size budgets | `ARCHITECTURE.yaml:code_budgets`; [meaning and review](code-budgets.md) |
| Execution and check selection | root `AGENTS.md`, `tools/check-affected.sh`, `tools/check-all.sh` |

Mathematics defines required behavior. Code records actual behavior. Evidence
supports bounded claims about their agreement and physical adequacy. Neither
an architecture scan nor a successful solve proves model validity.

## Review rules

Reviewers apply these to every readiness and delivered review:

- **Falsifier tolerance.** A falsifier's tolerance is at least as tight as the
  claim it protects (a 1e-12 state-preservation claim needs a 1e-12 check, not
  an existing 1e-6 composition test).
- **Evaluation point.** Evidence is evaluated at the point the claim concerns:
  a solved-endpoint claim uses the Jacobian at the solved coordinates, not a
  generic-point rank; a returned-state claim compares returned states.
- **Causal claims.** "Caused by" needs an isolating comparison (one variable
  changed, the other inputs and the wheel identified). Otherwise report the
  difference and its conditions.
- **Numbers in prose** come from retained data or are computed inline from it.
- **Speed.** The owner's no-slowdown rule for ionic LLE and the MEA flash is
  judged on work counts (iterations, evaluations, Hessians); paired wall time is
  supporting evidence with its spread.

## Vocabulary and conventions

Project vocabulary: [GLOSSARY.md](../../GLOSSARY.md).

Use bubble, dew, saturation and flash as mathematical specifications or data
constructors. Shared solvers must not select physics by an application or paper
name. All equilibrium specifications use the one driver described above; lower-level
kernels and search controllers retain their distinct mathematical roles. ADR 0011 admits
only an optional Gross-Vrabec dipole-dipole residual term; other residual multipoles remain
retired by ADR 0006, and dipolar relative permittivity is a distinct input and equation.


`DensityDiagnostics` means a bounded fixed-temperature/composition report of
pressure crossings, spinodals ($P_\rho=0$), declared connected cells and
finite-evaluation brackets; [ALG-DENSITY](../../engine/docs/science/algorithms.md#alg-density)
owns these operations. `DensityDiagnosticPoint` records SI density, packing,
pressure and slope only for available evaluations; `cause` marks unavailable
samples. Root brackets carry pressure residual in Pa or spinodal slope divided
by $RT$. For `finite_domain`, only density bounds, stopping cause and
availability cause are defined; `evaluation_order` names the consumed request.
`DensityDiagnosticBracket` retains the original declared-cell bounds and its
two pressure residuals in Pa. Its count excludes adaptive discoveries. These
quantities do not assert complete root enumeration or physical validation.

## Implementation completion order

The [GitHub milestones](https://github.com/tannerpolley/ePC-SAFT/milestones)
and owning issues define the current completion order, status, scope,
acceptance, reuse decisions and estimates. Native issue relationships own
prerequisites; milestone numbering expresses delivery priority, not an additional
scientific dependency. Do not maintain a competing sequence here.
Finish and verify each prerequisite's current agreed scope before dependent
implementation; deferred UV and global certification do not block current work.
Later integration defects reopen the affected earlier issue. Design readiness is
scoped to the affected slice; the leaf issue owns its gates and execution plan.
Infrastructure needed to execute checks is a prerequisite. Existing higher
layers may be used to verify a lower layer without claiming their own completion.

## Source/data adoption

source/data adoption: zotero

## Analysis workflow choice

The owner confirmed deferred Snakemake adoption for these populated analyses
on 2026-09-29. Existing study commands remain in use; a later adoption needs
confirmation.

| Analysis | Snakemake adoption |
|---|---|
| `analyses/2026-association-native-comparison` | Deferred |
| `analyses/2026-binary-methane-ethane-bubble-point` | Deferred |
| `analyses/2026-binary-methane-ethane-tp-flash` | Deferred |
| `analyses/2026-density-closure-diagnosis` | Deferred |
| `analyses/2026-dielectric-alternatives` | Deferred |
| `analyses/2026-fitting-state-reuse` | Deferred |
| `analyses/2026-fixed-pressure-anchor` | Deferred |
| `analyses/2026-grepe-hard-cut` | Deferred |
| `analyses/2026-idaes-asl-adapter` | Deferred |
| `analyses/2026-initialization-diagnosis` | Deferred |
| `analyses/2026-khudaida-idaes-settler` | Deferred |
| `analyses/2026-khudaida-lineage-diagnosis` | Deferred |
| `analyses/2026-mea-association-diameter` | Deferred |
| `analyses/2026-mea-grepe-regression` | Deferred |
| `analyses/2026-mea-parameter-estimation` | Deferred |
| `analyses/2026-mea-reactive-fit-timing` | Deferred |
| `analyses/2026-mea-recovered-state-actions` | Deferred |
| `analyses/2026-neutral-pcsaft-equilibrium-comparison` | Deferred |
| `analyses/2026-open-source-equilibrium-comparison` | Deferred |
| `analyses/2026-polar-dipole` | Deferred |
| `analyses/2026-pure-saturation-regression` | Deferred |
| `analyses/2026-reacting-mixed-actions` | Deferred |
| `analyses/2026-reactive-electrolyte-driver` | Deferred |
| `analyses/2026-reactive-liquid-cold-start` | Deferred |
| `analyses/2026-shared-derivative-qualification` | Deferred |
| `analyses/2026-stability-noise-band` | Deferred |
| `analyses/2026-trace-action-diagnosis` | Deferred |
| `analyses/2026-wolbach-sandler-requalification` | Deferred |
| `analyses/phase-stability-cutover` | Deferred |
| `analyses/reproduction/1989-esteso` | Deferred |
| `analyses/reproduction/1998-chase` | Deferred |
| `analyses/reproduction/2001-gross` | Deferred |
| `analyses/reproduction/2002-gross` | Deferred |
| `analyses/reproduction/2004-glos` | Deferred |
| `analyses/reproduction/2005-cameretti` | Deferred |
| `analyses/reproduction/2008-held` | Deferred |
| `analyses/reproduction/2010-nasrifar-mdea` | Deferred |
| `analyses/reproduction/2012-held` | Deferred |
| `analyses/reproduction/2013-pahlavanzadeh-mdea` | Deferred |
| `analyses/reproduction/2014-held` | Deferred |
| `analyses/reproduction/2015-baygi` | Deferred |
| `analyses/reproduction/2015-may` | Deferred |
| `analyses/reproduction/2015-uyan` | Deferred |
| `analyses/reproduction/2018-iapws` | Deferred |
| `analyses/reproduction/2018-wangler` | Deferred |
| `analyses/reproduction/2019-bulow` | Deferred |
| `analyses/reproduction/2019-iapws` | Deferred |
| `analyses/reproduction/2020-bulow` | Deferred |
| `analyses/reproduction/2020-cleeton` | Deferred |
| `analyses/reproduction/2020-haslam` | Deferred |
| `analyses/reproduction/2021-bulow` | Deferred |
| `analyses/reproduction/2021-bulow-sour-gas` | Deferred |
| `analyses/reproduction/2022-ascani` | Deferred |
| `analyses/reproduction/2023-ascani` | Deferred |
| `analyses/reproduction/2023-belov` | Deferred |
| `analyses/reproduction/2023-schick` | Deferred |
| `analyses/reproduction/2024-hubach` | Deferred |
| `analyses/reproduction/2024-yu` | Deferred |
| `analyses/reproduction/2025-figiel` | Deferred |
| `analyses/reproduction/2025-figiel-engine-crosswalk` | Deferred |
| `analyses/reproduction/2025-pakravesh` | Deferred |
| `analyses/reproduction/2025-perdomo` | Deferred |
| `analyses/reproduction/2025-rezaee` | Deferred |
| `analyses/reproduction/2026-09-03-stage-05-repair` | Deferred |
| `analyses/reproduction/2026-eos-literature-handoff` | Deferred |
| `analyses/reproduction/2026-khudaida` | Deferred |
| `analyses/reproduction/2026-rezaee` | Deferred |

Declared equilibrium starts means caller-supplied initial physical guesses or full
compiled coordinates for one unchanged thermodynamic problem (ALG-EQUILIBRIUM).
`EquilibriumStart` carries those guesses; `EquilibriumWorkLimits` bounds attempts,
optimizer iterations and elapsed seconds; `EquilibriumStartsResult` records ordered
attempts, residual-bearing results, initialization causes, consumed work and the
first admitted index. `DriverLimits` passes the remaining work into the existing
Ipopt driver. These names describe inputs, quantities and operations, not new
thermodynamic models.

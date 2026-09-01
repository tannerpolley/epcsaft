# General Reactive Equilibrium and Phase Equilibrium

- Status: configured finite-phase GREPE implemented
- Date: 2026-08-14
- Owner: Equilibrium vertical-slice task
- Authority: ADR 0007 mapping of the retained D-036 decision
- Public owner: `epcsaft.equilibrium.phase_equilibrium`

## Outcome

General Reactive Equilibrium and Phase Equilibrium (GREPE) evolves the existing
public `phase_equilibrium` owner into one declarative equilibrium compiler. It
defines one target physical vocabulary for homogeneous reaction, reactive
bubble/dew/VLE, and configured-phase reactive equilibrium.

GREPE is an architecture, not an application route or a globality claim. MEA,
CO2, lithium, Ascani, source names, and operation labels never select runtime
behavior. The request's state roles, phase roles, species support, reactions,
and model capabilities determine the assembled problem.

The implementation routes fixed-T,P homogeneous reaction, fixed-topology
reactive bubble/VLE, and two-or-more configured finite phases with or without
reactions through the same public request/result and generic Regression boundary. The finite
multiphase formulation supports neutral and charged systems, electrochemical
phase transfer, source-standard-state conversion, exact selected-parameter
sensitivities, and arbitrary configured fluid-role labels. It attempts the
caller-provided configuration and reports mathematical, Provider, or solver
failure directly. Automatic phase discovery, incipient multiphase LLE, solved
temperature, solids, precipitation, and kinetics need additional formulations.

## Public API

The solve boundary is `phase_equilibrium(model, problem)`. The model supplies
the EOS and ordered components. The problem supplies state variables, phases,
reactions, starts, and requested outputs. Application and source names are
metadata only and never choose equations.

The model and condition together contain:

- one immutable resolved `Parameters` selection and optional ordered
  `ActiveParameterSet`, supplied by the separate model argument;
- canonical species identities, formula/conservation vectors, charges, and
  aliases resolved before the solver;
- one source-neutral resolved feed in amounts, mole fractions, mass fractions,
  or solvent-composition-plus-molality form; every conversion exposes its exact
  affine amount Jacobian before thermodynamic compilation;
- `Fixed(value)` or `Solved(initial, bounds)` temperature and pressure roles;
- one or more phases with stable IDs, fluid role, finite or incipient role,
  component support, and any deterministic label anchor;
- an optional ordered reaction stoichiometry matrix with eligible phases and
  typed standard-state records when reactions are present;
- explicit units and composition bases; and
- typed requested outputs.

The JSON and in-memory loaders produce the same condition. Unknown fields,
callbacks, executable expressions, ambiguous units, and implicit composition
bases fail before the numerical solve.

Thin chemical-engineering helpers such as `bubble_pressure`,
`bubble_temperature`, `dew_pressure`, and `dew_temperature` may
construct the canonical request. They do not own equations, result types, or
solver selection. The core has no bubble/LLE/application dispatch enum.

## One thermodynamic encoding

GREPE encodes each equilibrium condition exactly once.

### Finite phases

A finite phase has a positive phase amount. The compiler constructs one
constrained total-Gibbs formulation over phase/species amounts with the
declared conservation, charge, support, and state constraints. A generic sparse
phase/species direction row may describe a phase-local reaction, interphase
transfer, ion exchange, or any coupled combination. The native compiler proves
direction rank, global conservation, per-phase charge conservation, and the
orthogonal invariant complement before Ipopt sees the problem. For conservation
matrix `C`, phase reaction matrices `S_p`, and the support-restricted transfer
basis `T`, the permitted composition nullspace is the declared span of the
phase blocks `S_p^T` and `T`, with every column in `null(C)`. Preflight proves
equality with that declared admissible span and freezes every other direction,
so the minimization cannot invent an undeclared reaction. Reaction affinities and
interphase chemical-potential equalities over shared declared support arise
from its KKT conditions; they are not imposed a second time as independent
equations.

The solver establishes a local stationary equilibrium for the declared phases.
Nonconvexity means an Ipopt solution alone is not a global phase-stability
certificate. Every result therefore records `globality_claim =
not_established` unless a separate discovery/stability procedure
proves more.

### Incipient phases

An incipient phase has zero phase amount and contributes no material-balance
split. Because total Gibbs cannot determine a zero-amount phase composition,
the compiler lifts the finite phase's Gibbs KKT equations into the same Ipopt
feasibility problem with the incipient normalization and required chemical-
potential or fugacity equalities over the incipient/shared declared support.
It does not attach those equalities directly to the primal finite-phase
minimization, which would alter its stationarity. The coupled solve therefore
retains one physical encoding of the finite reacting phase while solving the
incipient composition and state variable.

A finite phase approaching zero, phase coalescence, support loss, label
exchange, active-inequality ambiguity, or a singular chart is a typed
Non-Evaluable Trial. GREPE never silently changes between finite and incipient
formulations during a differentiable evaluation.

### Reaction standards

A reaction standard-state record declares:

- eligible phase ID or phase role;
- dimensionless natural-log equilibrium convention or absolute standard
  chemical-potential representation;
- activity/concentration basis, reference pressure, units, independent
  variables, validity domain, and coefficient identities;
- temperature and pressure dependence; and
- the electrostatic reference when charged standard potentials require it.

Preflight verifies independent reaction rank relative to the permitted
nullspace of the conservation matrix, material conservation, charge
conservation, dimensions, domain, and consistency between any supplied
`ln(K)` and standard chemical potentials. Redundant reactions and incomplete
or extraneous reactive directions fail before solve construction. The compiler
converts accepted input to one internal potential representation and never
enforces both representations independently. A phase role is enabled for
reaction only when its standard-state convention is complete. Initial GREPE
slices keep reactions liquid-only unless equally complete vapor or second-
liquid evidence exists.

An empty reaction matrix compiles ordinary species conservation with no
reaction-affinity rows. It does not require a fictitious reaction or select a
different finite-phase solver.

### Electrochemical phases

The finite multiphase formulation uses phase-local electroneutrality,
electrochemical-potential equality over shared declared ionic support
`mu_tilde_i = mu_i + z_i F phi_phase`, and one explicit electrostatic-potential
gauge. Ordinary chemical-potential equality is insufficient for ions.

## Numerical method and derivatives

Ipopt solves one thermodynamic condition at a time. Its configured formulation
provides exact objective and constraint gradients and the exact Lagrangian
Hessian from the CppAD authority plus exact implicit linear algebra. Limited-
memory Hessians, numerical derivatives, and solver-history differentiation
cannot support the production claim.

For an incipient formulation, finite-phase stationarity is itself a lifted
constraint. Its exact Ipopt Lagrangian Hessian therefore requires the
corresponding third-order Gibbs/EOS tensors. That formulation is not yet
implemented; first-order output sensitivity is not sufficient.

The finite-phase formulation normalizes total feed to one mole internally. It
uses logarithmic amount and packing coordinates, dimensionless `G/RT`, and
feed-scaled conservation rows, then rescales amounts, volumes, free energy, and
their parameter derivatives at the public boundary. The conservation
certificate is the dimensionless feed-relative residual of those scaled rows;
pressure, affinity, transfer, and KKT certificates retain their existing
tolerances. Numerical scaling is otherwise internal to Ipopt. A bounded scalar
packing-coordinate Newton refinement tightens the direct fixed-pressure
residual after Ipopt without changing phase amounts or composition, after which
the equality multipliers are recomputed at the returned state.

Each certified fixed-branch result exposes requested values and exact total
Jacobians in `ActiveParameterSet` order. Finite differences and high-precision
oracles validate those derivatives only. Branch/support changes,
ill-conditioning, and unavailable derivative order fail closed.

Starts are deterministic functions of the request plus explicit hints. For two
configured finite phases, an omitted start set is generated from
conservation-preserving reaction-sign, neutral-component, and equal-feed
partitions with Provider fixed-pressure volume roots. The first completely
certified candidate selects the local branch. Explicit starts bypass generation
and select the intended local basin. Cached continuation states may reduce work
only for an equivalent request identity and must not make values, branches, or
derivatives depend on evaluation history.

## Results and failures

One result contains:

- ordered finite and incipient phase states, compositions, densities, and
  amounts where defined;
- reaction/speciation state and requested chemical potentials, fugacities,
  phase fractions, total pressure, or species partial pressures;
- exact selected-parameter Jacobians;
- solver convergence, numerical convergence, and physical-validity statuses;
- the local branch, requested phase set, support, condition estimates,
  residual certificates, and explicit globality claim; and
- exact request, parameter, data, software, solver-setting, and result-table
  identities/hashes.

Failure distinguishes invalid configuration, unavailable formulation,
nonconvergence, no incipient phase, requested-phase disappearance, unresolved
or exchanged branches, provider-domain rejection, support loss, singular
sensitivity, and failed physical certification. A useful candidate may be
returned with a non-success status; it is never an ordinary observation value.

The current proof targets two configured finite fluid phases. The underlying
finite-phase formulation accepts additional configured phases using the same
variables and equations. Solids and zero-amount incipient phases require
different state variables.

The installed numerical evidence now includes a neutral reactive LLE, a
nonreactive electrolyte liquid/vapor topology, exact finite-phase parameter
Jacobians, amount-scale invariance over six orders of magnitude, a generic
Regression recovery through the public GREPE owner, and a sparse lithium/sodium
exchange oracle. The electrolyte finite-phase proof currently keeps every
Provider component in every finite phase; exact liquid-only ion support in a
finite vapor and fixed-pressure incipient equilibrium with a solved generic feed
coordinate remain implementation work. The solved-bubble formulation already
supports a declared neutral-only incipient vapor.

## Generic Regression boundary

Regression consumes one owner-published evaluator protocol:

1. describe supported parameter identities, output selectors, units, bases,
   shapes, domains, and derivative availability;
2. evaluate each unique thermodynamic condition once for all requested rows;
3. return typed values, exact Jacobian rows, certificates, and failures.

Regression owns ordered differences, optional diagonal discrepancy scales,
physical bounds, affine coordinates, one start, native Ceres controls, and
optimizer accounting. It contains no species, reaction, phase, bubble, LLE,
application-specific branch, multistart, or Python iteration callback.

Only compiled EOS and anchored local D-034/D-035 Equilibrium owners enter the
current route. Unsupported formulations refuse before Ceres is constructed.
An invalid owner trial returns a typed failure; Regression does not fabricate
penalties or stale fitted values. A usable result exists only after a fresh,
cache-bypassing final owner value/Jacobian evaluation.

## Observation, fit-quality, and plotting contract

Each compiled owner sidecar retains ordered row identities, units, bases,
source and evaluator fingerprints, and topology/branch evidence. A
`LeastSquaresResult` retains ordered observations, predictions, raw and scaled
differences, physical and optimizer-coordinate Jacobians, final parameters,
bounds, solver options, failures, counters, and local covariance availability.
Diagonal residual scales remain caller-defined discrepancy weights unless the
caller independently establishes a stronger statistical meaning.

Campaign statistics, grouped fit quality, plots, multistart orchestration,
candidate document export, and parameter adoption belong to the calling
application or analysis. Regression makes no globality, uncertainty,
identifiability, parameter-selection, promotion, or scientific-adoption claim.

## Ownership

- Engine owns GREPE, generic Regression, exact derivatives, compact numerical
  evidence, and the one installed wheel.
- Data owns reusable canonical parameter and observation packets with source,
  units, domain, transformation, use basis, and hashes.
- Analysis sets own durable installed-wheel literature and cross-module evidence.
- A downstream MEA or lithium application owns its system configuration,
  selected rows and mappings, campaign controls, analysis, plots, fit-quality
  bands, and parameter-adoption decision. It pins immutable Data inputs and an
  immutable installed Engine wheel; it does not copy thermodynamic algorithms
  or canonical reusable packets.

The sub-ten-second warm developer gate retains compact Engine evidence. Slow
literature campaigns, broad regression runs, and publication plots run only by
explicit owner commands and do not enter that gate.

## Staged implementation

The generic GREPE owner and configured finite-phase Gibbs formulation remain
the Equilibrium authority. Native Regression compiles only the admitted,
anchored local D-034/D-035 slices; later formulations require their own
implementation issue and numerical evidence.

1. **Generic seam and schema.** Role-based request/result types, strict
   loaders, phase-aware reaction preflight, and the generic evaluator protocol
   are retained. Campaign tables and candidate export remain outside Engine
   Regression.
2. **Reactive bubble first.** Route the D-035 full nonideal-vapor bubble solve
   through the evolved public owner and generic evaluator. Prove identical
   physical values, certificates, and exact Jacobians. Cut MEA pressure work to
   `y_i P`; family-specific Regression dispatch is retired.
3. **Homogeneous reaction.** Route the D-034 fixed-T,P owner through the same
   finite-phase compiler and evaluator, preserve its source-reference and
   local-branch evidence, then remove its specialized public route.
4. **Configured reactive LLE.** Implement two or more finite phases
   using exact or
   manufactured Gibbs/KKT evidence and installed-wheel tie-line reproduction.
   Ascani/Senina-style endpoints without measured feed or phase fraction test
   conditional tie-line composition/equilibrium only. A documented midpoint
   feed is reproduction input, never independent phase-fraction validation or
   a fitted latent feed fraction.
5. **Electrochemical reactive LLE.** Implemented with exact per-phase
   electroneutral coordinates, Galvani-gauge-invariant transfer residuals,
   source-standard-state transforms, and exact selected-parameter derivatives.

Each slice uses one accountable Engine task across its exact Data, Engine, and
Validation files. A specialized implementation is removed only after reviewed
equivalence and new-claim evidence; the final public surface retains no
compatibility shim or duplicate equation owner.

## Numerical evidence

Before a slice can broaden current capability, Engine evidence must falsify:

- conservation, charge where applicable, normalization, KKT/potential, and
  phase-role certificates in original physical units;
- exact objective/constraint gradients, exact Lagrangian Hessians, and exact
  output Jacobians against independent derivative oracles;
- phase-label permutation and deterministic anchor behavior;
- reaction-rank compatibility with the permitted conservation nullspace;
- structural exclusion and trace-support behavior;
- phase disappearance, coalescence, branch ambiguity, singularity,
  nonconvergence, and unavailable-formulation failure semantics; and
- equivalence to every displaced route within physically justified
  tolerances.

Installed-wheel analysis sets then own bounded physical/literature comparisons,
exact plotted-data tables, and figures where required.

## Completion criteria

The admitted GREPE core covers homogeneous reaction, solved-pressure
bubble/VLE, imposed intensive boundaries, and declared-phase reactive LLE
through the one public `phase_equilibrium` owner and one generic Regression
evaluator. Temperature-solved Bubble-T and Dew-T observations remain issue
#42. Displaced specialized public routes, reduced tracers, family switches,
and stale live guidance are removed, and exact derivative and numerical
evidence passes through the installed wheel. Configured neutral and
electrochemical two-phase equilibrium is proven; automatic phase discovery and
solids remain separate numerical formulations.

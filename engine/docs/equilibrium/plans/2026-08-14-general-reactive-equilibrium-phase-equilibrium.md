# General Reactive Equilibrium and Phase Equilibrium

- Status: neutral reactive-LLE implementation candidate; formal admission pending
- Date: 2026-08-14
- Owner: Equilibrium vertical-slice task
- Authority: Governance D-036, doctrine revision 21
- Public owner: `epcsaft.equilibrium.phase_equilibrium`

## Outcome

General Reactive Equilibrium and Phase Equilibrium (GREPE) evolves the existing
public `phase_equilibrium` owner into one declarative equilibrium compiler. It
defines one target physical vocabulary for homogeneous reaction, reactive
bubble/dew/VLE, and declared-phase reactive LLE. A later independently admitted
slice may add electrochemical reactive LLE.

GREPE is an architecture, not an application route or a globality claim. MEA,
CO2, lithium, Ascani, source names, and operation labels never select runtime
behavior. The request's state roles, phase roles, species support, reactions,
and model capabilities determine the assembled problem.

Current admission remains Governance D-034 fixed-T,P homogeneous reactive and
D-035 fixed-topology reactive bubble/VLE observations. The feature branch now
implements the first two-finite-liquid neutral reactive-LLE candidate through
the same public request, result, receipt, and generic Regression boundary. Its
compact manufactured installed-wheel evidence proves the declared local
minimum, balances, reaction and transfer potentials, exact selected-parameter
sensitivities, label covariance, failure semantics, and Regression recovery.
Formal admission still requires review and an accepted hash-bound receipt.
Reactive phase discovery, electrochemical multiphase work, three-plus reacting
phases, solids, precipitation, kinetics, and reaction-correlation fitting remain
unavailable.

## Public contract

The canonical request contains:

- one immutable resolved `Parameters` document and optional ordered
  `ActiveParameterSet`;
- canonical species identities, formula/conservation vectors, charges, and
  aliases resolved before the solver;
- `Fixed(value)` or `Solved(initial, bounds)` temperature and pressure roles;
- an explicit composition role: overall feed, analytical/apparent totals, or
  known-phase composition;
- one or more phases with stable IDs, fluid role, finite or incipient role,
  allowed species, and any deterministic label anchor;
- an ordered reaction stoichiometry matrix, eligible phases, and typed
  standard-state records;
- explicit units and composition bases; and
- typed requested outputs.

The strict JSON and in-memory mapping loaders share one schema and produce the
same immutable request. Registered identifiers select EOS models, reaction
correlations, and residual transforms. Unknown fields, arbitrary import paths,
callbacks, executable expressions, ambiguous units, and implicit composition
bases fail preflight.

Thin chemical-engineering helpers such as `bubble_pressure`,
`bubble_temperature`, `dew_pressure`, and `dew_temperature` may
construct the canonical request. They do not own equations, result types, or
solver selection. The core has no bubble/LLE/application dispatch enum.

## One thermodynamic encoding

GREPE encodes each equilibrium condition exactly once.

### Finite phases

A finite phase has a positive phase amount. The compiler constructs one
constrained total-Gibbs formulation over phase/species amounts with the
declared conservation, charge, support, and state constraints. Phase-local
reaction extents restrict total-composition changes to the declared
stoichiometric subspace; species-transfer coordinates redistribute the same
species only across the participating phases' shared declared support and sum
to zero. For conservation matrix `C`, phase reaction matrices `S_p`, and the
support-restricted transfer basis `T`, the permitted composition nullspace is
the declared span of the phase blocks `S_p^T` and `T`, with every column in
`null(C)`. Preflight proves equality with that declared admissible span and
freezes every other direction, so the minimization cannot invent an undeclared
reaction. Reaction affinities and
interphase chemical-potential equalities over shared declared support arise
from its KKT conditions; they are not imposed a second time as independent
equations.

The solver establishes a local stationary equilibrium for the declared phases.
Nonconvexity means an Ipopt solution alone is not a global phase-stability
certificate. Every result therefore records `globality_claim =
not_established` unless a separate admitted discovery/stability procedure
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

### Later electrochemical phases

The future electrolyte reactive-LLE slice adds phase-local electroneutrality,
electrochemical-potential equality over shared declared ionic support
`mu_tilde_i = mu_i + z_i F phi_phase`, and one explicit electrostatic-potential
gauge. Ordinary chemical-potential equality is insufficient for ions. This
schema target does not constitute capability admission.

## Numerical and derivative contract

Ipopt solves one thermodynamic condition at a time. Its admitted formulation
provides exact objective and constraint gradients and the exact Lagrangian
Hessian from the CppAD authority plus exact implicit linear algebra. Limited-
memory Hessians, numerical derivatives, and solver-history differentiation
cannot support the production claim.

For an incipient formulation, finite-phase stationarity is itself a lifted
constraint. Its exact Ipopt Lagrangian Hessian therefore requires the
corresponding third-order Gibbs/EOS tensors. The slice remains unsupported
until CppAD supplies and independent oracles verify every derivative order the
formulation consumes; first-order output sensitivity is not sufficient.

Numerical scaling is internal to Ipopt. Certification evaluates original
unscaled conservation, charge, reaction/KKT, normalization, and phase-potential
closures at their unchanged physical tolerances. No post-Ipopt root polish is
part of the initial contract.

Each certified fixed-branch result exposes requested values and exact total
Jacobians in `ActiveParameterSet` order. Finite differences and high-precision
oracles validate those derivatives only. Branch/support changes,
ill-conditioning, and unavailable derivative order fail closed.

Starts are deterministic functions of the request plus explicit hints. Cached
continuation states may reduce work only for an equivalent request identity and
must not make values, branches, or derivatives depend on evaluation history.

## Result and failure contract

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

Failure distinguishes invalid configuration, unsupported capability,
nonconvergence, no incipient phase, requested-phase disappearance, unresolved
or exchanged branches, provider-domain rejection, support loss, singular
sensitivity, and failed physical certification. A useful candidate may be
returned with a non-success status; it is never an ordinary observation value.

First admission supports one or two fluid phases. Three-plus phases and solids
return typed `unsupported_topology`; declared problems never fall back to
discovery.

## Generic Regression boundary

Regression consumes one owner-published evaluator protocol:

1. describe supported parameter identities, output selectors, units, bases,
   shapes, domains, and derivative availability;
2. evaluate each unique thermodynamic condition once for all requested rows;
3. return typed values, exact Jacobian rows, certificates, and failures.

Regression owns observation mapping, residual transforms and covariance
whitening, parameter coordinates/bounds/ties, multistart, Ceres, optimizer
accounting, and descriptive/statistical diagnostics. It contains no species,
reaction, phase, bubble, LLE, or application-specific thermodynamic branch.
Current homogeneous and bubble family dispatch will be removed after the
generic protocol proves equivalent; no future LLE `isinstance` branch is
added.

The layering is nested and explicit: inner GREPE/Ipopt solves each condition;
outer Regression/Ceres changes parameters and minimizes the declared dataset
SSE using GREPE's exact total Jacobian. Parameters and all equilibrium states
are not combined into one monolithic campaign-wide Ipopt problem.

An invalid trial returns a typed failure so Ceres rejects or shrinks the step.
Regression does not fabricate penalties or drop failed rows. No usable
certified point produces a typed fit failure.

## Observation, fit-quality, and plotting contract

Every observation declares row/source/series identity, target selector,
physical unit and basis, uncertainty/covariance or preregistered residual
scale, partition, and grouping metadata. Structural zeros, measured zeros, and
censoring are distinct. Composition residuals use a closure-aware `N-1`
contrast while retaining every raw reported component.

`FitResult` stores an immutable evaluated-row table containing at least the row
and source IDs, selector, series metadata, observed and predicted values, raw
difference, scaled residual, unit/basis, evaluation status, candidate parameter
identity, and table hash. Generic plotting reads that exact table and never
re-evaluates thermodynamics. A CO2 pressure plot is therefore an ordinary
`vapor.partial_pressure(species_id)` output plot, not a CO2 code path.

Fit diagnostics include overall and grouped SSE, RMSE, AARD, valid-equilibrium
coverage, rank, conditioning, parameter correlation, active bounds, and
multistart clusters. Chi-square, covariance, and confidence intervals appear
only when declared uncertainty and rank assumptions justify them. A downstream
dataset may preregister descriptive fit-quality bands; Regression reports them
without deciding parameter adoption.

`FitResult` exposes one explicit candidate-export method. It creates a new
provenance-linked parameter document and never overwrites the source. An
exportable point requires a certified owner evaluation and exact Jacobian for
every required fit row; failed-row subsets never define a partial SSE. Among
complete usable starts, deterministic lowest objective and then stable start
identity select the candidate. Its immutable evaluated-row-table hash is bound
to the export. A budget-limited or nonconverged fit may export and plot that
complete certified point with its exact termination label. Partial evaluated
rows may remain diagnostics but cannot be exported as the fit candidate. No
complete certified point means no export. Regression has no parameter
acceptance, admission, promotion, or adoption authority.

## Ownership

- Engine owns GREPE, generic Regression, exact derivatives, compact numerical
  evidence, and the one installed wheel.
- Data owns reusable canonical parameter and observation packets with source,
  units, domain, transformation, use basis, and hashes.
- Validation owns durable installed-wheel literature and cross-module evidence.
- A downstream MEA or lithium application owns its system configuration,
  selected rows and mappings, campaign controls, analysis, plots, fit-quality
  bands, and parameter-adoption decision. It pins immutable Data inputs and an
  immutable installed Engine wheel; it does not copy thermodynamic algorithms
  or canonical reusable packets.

The sub-ten-second warm developer gate retains compact Engine evidence. Slow
literature campaigns, broad regression runs, and publication plots run only by
explicit owner commands and do not enter that gate.

## Staged implementation

Implementation checkpoint on this branch: stages 1--3 are routed through the
generic owner and evaluator, and stage 4 is implemented as a neutral-only,
exactly-two-finite-liquid candidate. The retained specialized homogeneous and
bubble entry points remain temporarily callable because downstream MEA work
still imports them; their thermodynamic and Regression ownership has already
moved behind GREPE. Removing those names is a coordinated downstream cutover,
not part of the neutral-LLE numerical claim.

1. **Generic seam and schema.** Add role-based request/result types, strict
   loaders, phase-aware reaction preflight, generic evaluator protocol, exact
   evaluated-row table, authority-neutral candidate export, and unambiguous
   `best_usable_start` optimizer vocabulary. Retain current admitted algorithms
   as evidence, not alternate public owners.
2. **Reactive bubble first.** Route the D-035 full nonideal-vapor bubble solve
   through the evolved public owner and generic evaluator. Prove identical
   physical values, certificates, and exact Jacobians. Cut MEA pressure work to
   `y_i P` and retire the reduced liquid-fugacity tracer guidance and family-
   specific Regression dispatch.
3. **Homogeneous reaction.** Route the D-034 fixed-T,P owner through the same
   finite-phase compiler and evaluator, preserve its source-reference and
   local-branch evidence, then remove its specialized public route.
4. **Neutral reactive LLE.** Implement and admit two finite neutral liquids
   using exact or
   manufactured Gibbs/KKT evidence and installed-wheel tie-line reproduction.
   Ascani/Senina-style endpoints without measured feed or phase fraction test
   conditional tie-line composition/equilibrium only. A documented midpoint
   feed is reproduction input, never independent phase-fraction validation or
   a fitted latent feed fraction.
5. **Electrochemical reactive LLE.** Add the charged formulation only after
   phase-aware ionic standard states, gauge handling, exact derivatives, and an
   independent source-complete two-liquid electrolyte case exist.

Each slice uses one accountable Engine task across its exact Data, Engine, and
Validation files. A specialized implementation is removed only after reviewed
equivalence and new-claim evidence; the final admitted surface retains no
compatibility shim or duplicate equation owner.

## Admission evidence

Before a slice can broaden current capability, Engine evidence must falsify:

- conservation, charge where applicable, normalization, KKT/potential, and
  phase-role certificates in original physical units;
- exact objective/constraint gradients, exact Lagrangian Hessians, and exact
  output Jacobians against independent derivative oracles;
- phase-label permutation and deterministic anchor behavior;
- reaction-rank compatibility with the permitted conservation nullspace;
- structural exclusion and trace-support behavior;
- phase disappearance, coalescence, branch ambiguity, singularity,
  nonconvergence, and unsupported topology failure semantics; and
- equivalence to every displaced admitted route within physically justified
  tolerances.

Installed-wheel Validation then owns the bounded physical/literature slice,
exact plotted-data table, figures where required, artifact hashes, and
independent review. A user-approved receipt admits only that demonstrated
slice. Design approval and low SSE alone never broaden runtime authority.

## Completion criteria

GREPE is complete only when all admitted homogeneous, bubble/dew/VLE, and
declared-phase reactive-LLE configurations use the one public
`phase_equilibrium` owner and one generic Regression evaluator; displaced
specialized public routes, reduced tracers, family switches, and stale live
guidance are removed; exact derivative and certificate evidence passes through
the installed wheel; and the applicable promotion receipts are accepted.

Electrochemical reactive LLE, phase discovery, three-plus phases, and solids
remain separately incomplete until their own criteria and receipts are met.

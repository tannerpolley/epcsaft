# Engine Context

`governance_doctrine_revision: 24`

`engine/` is the sole runtime owner inside `tannerpolley/ePC-SAFT`. It
builds one Python 3.13 distribution, `epcsaft`, containing EOS, equilibrium,
and regression. All native modules share one Meson graph and wheel.

Reusable packets remain owned by `data/` and are materialized with exact
commit, packet, fingerprint, and file hashes. Regression returns an immutable
`FitResult`. D-036 requires the promotion-gated transfer to become an explicit
candidate export that creates a new provenance-linked `Parameters` document
without accepting, admitting, promoting, or adopting it. No package, SDK,
filesystem, or legacy compatibility handoff exists between the three runtime
modules.

Application observations, including phase-equilibrium and extraction
observations, may be declared Regression targets when their units, basis,
topology, support, standard states, and residual meanings are explicit. A
successful fit produces a provenance-linked candidate parameter document for
the caller's stated domain. Engine numerical acceptance establishes that the
declared mathematical problem and its exact derivatives satisfy the stated
checks; it does not assert that fitted parameters are universal or
decide whether an application may use them. Application suitability,
calibration scope, surrogate construction, and flowsheet claims remain with
the downstream repository.

Engine retains compact numerical evidence. Installed-wheel physical and
literature campaigns belong in `validation/`; temporary experiments belong in
`lab/` and carry no durable evidence authority. Superseded live plans, designs,
implementation reviews, and source-specific runtime guidance are deleted.
Primary citations, research, source locators, raw data, receipts, results, and
still-relevant numerical evidence retain historical identities as provenance.

The canonical cross-module execution plan is
[`docs/plans/2026-08-07-general-scientific-capability-hard-cutover.md`](docs/plans/2026-08-07-general-scientific-capability-hard-cutover.md).
The approved GREPE target that evolves the sole public `phase_equilibrium`
owner is
[`docs/equilibrium/plans/2026-08-14-general-reactive-equilibrium-phase-equilibrium.md`](docs/equilibrium/plans/2026-08-14-general-reactive-equilibrium-phase-equilibrium.md).
Its full owner plans cover
[`EOS`](docs/eos/plans/2026-08-07-general-eos-parameter-derivatives.md),
[`Regression`](docs/regression/plans/2026-08-07-owner-driven-parameter-fit.md),
and [`Equilibrium`](docs/equilibrium/plans/2026-08-07-general-phase-equilibrium.md).

## Language

**Parameter Fit**:
A declared selection of ePC-SAFT parameters and observations to fit together
from one canonical parameter set.
_Avoid_: Regression problem, paper-specific fit, specialized solver

**Fit Capability**:
An Engine-advertised combination of parameter identity, observation family,
model domain, and derivative support that is eligible for a Parameter Fit.
_Avoid_: Author route, campaign solver, inferred support

**Observation Family**:
A source-neutral definition of a measured quantity, its thermodynamic inputs,
units or basis, residual meaning, and required derivatives.
_Avoid_: Paper target, source-specific evaluator

**Targetable Quantity**:
A public thermodynamic quantity with an explicit basis, unit, state domain,
data mapping, residual meaning, and CppAD selected-parameter derivative.
_Avoid_: Contribution breakdown, solver diagnostic, unqualified property value

**Source-Neutral Symbol**:
An Engine identity named for its physical model, equation, parameter,
observable, topology, or numerical operation rather than an author, paper,
campaign, or application. Established scientific model names such as Born
remain valid; citations and provenance may retain source names.
_Avoid_: Source-named class, function, capability, file, fixture, or test

**Parameter Identity**:
The unique ePC-SAFT parameter selected by family and model, component, pair, or
association scope. The same identity is one fitted coordinate everywhere.
_Avoid_: Native slot, paper coefficient, row-local parameter

**Active Parameter Set**:
The core Parameters-owned selection of Derivative-Ready Parameter identities
treated as independent CppAD inputs while every other parameter remains fixed.
EOS, Equilibrium, and Regression share this type without reversed module
dependencies.
_Avoid_: Capability block, parameter-family callback, fixed fit tuple

**Pair Interaction Coordinate**:
One unordered component-pair Parameter Identity applied reciprocally to its
declared cross-interaction rule. `k_ij` corrects cross dispersion, `l_ij`
corrects cross size, and a typed `k_HB,ij` corrects a selected
cross-association-energy combining rule.
_Avoid_: Directed pair slot, fixed binary callback, association-volume alias

**Cross-Association Correction**:
An unordered component-pair `k_HB,ij` resolved through a declared combining
rule onto eligible role-aware association edges with its exact derivative
chain. Site roles select edges; they do not make the numeric pair correction
directional.
_Avoid_: Generic association parameter, implicit cross edge, directed kHB

**Observable Derivative Block**:
An EOS- or Equilibrium-owned source-neutral set of observable row identities,
values, units/bases, domain, exact Jacobian columns in Active Parameter Set
order, and structured failures. It contains no observed targets or optimizer
policy.
_Avoid_: Regression target switch, native target index, duplicated evaluator

**Active Parameter Selection**:
An explicit ordered list of Parameter Identities resolved to one Active
Parameter Set. Preflight rejects missing, duplicate, discrete, or unsupported
identities; adding parameters to an input packet cannot silently expand a fit.
_Avoid_: Family wildcard, arbitrary JSONPath, implicit packet expansion

**Canonical Fit Schema**:
The single strict, versioned representation parsed by both persisted JSON and
in-memory Python mapping loaders. Both loaders perform identical validation and
produce the same Parameter Fit; neither has private shorthands or behavior.
_Avoid_: Dictionary-only shortcut, file-only feature, permissive key map

**Typed Parameter Correlation**:
A fixed, registry-owned functional family with explicit independent variables,
coefficient Parameter Identities, units, domain, and exact derivatives. Any
declared continuous coefficient may be selected, bounded, tied, and reported
independently without changing the correlation family during optimization.
_Avoid_: Whole-correlation activation, arbitrary expression, scalar-only fit

**Parameter Tie**:
An explicit equality relation that makes distinct Parameter Identities share
one fitted coordinate.
_Avoid_: Inferred sharing, arbitrary parameter expression

**Resolved Fit**:
A Parameter Fit whose supported identities, observation handling, starts,
bounds, scales, weights, and execution choices have been made explicit.
_Avoid_: Prepared problem, hidden defaults

**Fit Data Source**:
User-selected observation rows or a data file whose columns and units are
explicitly mapped to an Observation Family. Its path is a locator, not its
scientific identity.
_Avoid_: Required Data packet, untyped numeric array, path-based identity

**Fit Objective**:
A declarative residual representation, optional measurement weighting, and
loss applied to a Parameter Fit. Missing uncertainty may use a declared
Observation Family default but never implies statistical uncertainty.
_Avoid_: Custom callback, implicit weighting, raw mixed-unit error

**Derivative-Ready Parameter**:
A continuous, fixed-topology Parameter Identity with a valid domain and the
CppAD derivative required by its selected Observation Family.
_Avoid_: Stored-field allowlist, discrete fitted field, numerical fallback

**Molecular Moment**:
An independently measured or quantum-calculated permanent dipole or
quadrupole moment consumed by a polar EOS contribution. Molecular Moments are
fixed physical inputs, not ordinary Parameter Fit coordinates. An effective
fitted moment would require a future, explicitly named model and use case.
_Avoid_: Default fitted moment, unlabeled effective moment, inferred moment

**Permittivity Submodel**:
The auxiliary constitutive model that supplies state- and composition-dependent
relative static permittivity to primitive electrolyte EOS contributions. EOS
evaluates and may report its bulk value, while Parameter Fit calibrates its
coefficients only through an explicitly selected permittivity model family and
declared direct or downstream observations.
_Avoid_: PC-SAFT-generated permittivity, implicit participation in a base
parameter fit, untyped dielectric coefficient

**Model Family Namespace**:
A source-neutral, typed section of the single canonical `Parameters` document
owned by one base or auxiliary model family. It keeps coefficient identities,
domains, and model selection explicit without splitting parameter artifacts or
flattening unrelated coefficients into an inferred key map.
_Avoid_: Auxiliary parameter file, author-named field, unstructured flat map

**Staged Calibration Baseline**:
The normal sequence that first calibrates an auxiliary constitutive model to
its direct observations, then holds or explicitly reselects those coefficients
when fitting downstream EOS observations. It is a recommended workflow, not a
restriction on the generic optimizer.
_Avoid_: Implicit joint fit, permanently frozen auxiliary coefficient,
workflow-specific optimizer

**Joint Calibration**:
An explicit Parameter Fit whose Active Parameter Set spans more than one Model
Family Namespace. Direct auxiliary observations and downstream EOS observations
remain separate residual blocks with declared weights and family-specific
diagnostics so their influence and identifiability remain visible.
_Avoid_: Merged anonymous residual pool, automatic cross-family activation,
hidden staged result

**Association Topology**:
The fixed graph of component sites, advisory Association Site Roles, positive
multiplicities, and active site-pair interactions used by the general
association mass-action equations. Named schemes such as 1A, 2B, 3B, and 4C
are convenience presets that expand to this graph; they do not select solver
implementations.
_Avoid_: Scheme-specific solver, fitted site count, campaign-specific topology

**Association Site Role**:
A typed donor, acceptor, or effective annotation used to validate named
topology presets and the applicability of declared combining rules. The
explicit Association Topology remains authoritative: a role neither creates an
edge nor globally prohibits a cited effective same-role or self-pair model.
_Avoid_: Physically inert site class, inferred chemistry, universal edge ban

**Association Edge Source**:
The single declared origin of an active edge's association energy and volume:
either explicit values or one Association Combining Rule. Preflight rejects
missing, competing, and project-default sources, then materializes the resolved
values and their exact parameter derivative chain for evaluation and receipts.
_Avoid_: Source precedence, implicit default rule, unresolved cross edge

**Association Combining Rule**:
A cited, explicit rule that derives cross-association energy and volume from
declared component parameters on an Association Topology. Induced association
is represented by the resulting cross edges while self edges remain inactive;
it is not a specialized EOS or Regression route.
_Avoid_: Hidden cross parameter, campaign-specific induced-association path,
automatic chemistry inference

**Fit-Ready Association Parameter**:
An association energy or volume Parameter Identity whose exact selected-
parameter derivative is available for every supported pure and mixture,
neutral and electrolyte, constant- and temperature-dependent-diameter model
that uses its fixed Association Topology.
_Avoid_: Neutral-only fit slot, constant-diameter gate, pure-only derivative

**Implicit Sensitivity**:
A total derivative through a solved state assembled by exact linear algebra
from CppAD partial tensors of the defining equations.
_Avoid_: Solver-history differentiation, handwritten nonlinear derivative,
production finite difference

**Phase-Equilibrium Discovery**:
A nonreactive phase-equilibrium solve in which phase count and liquid/vapor
roles are unknown and the owner may certify one stable phase or a two-/three-
phase split.
_Avoid_: Source-case route, assumed two-phase flash, local solve called global

**Constrained Phase Equilibrium**:
A nonreactive solve for one to three caller-declared liquid/vapor phase roles.
Its certificate covers only that topology unless discovery separately proves a
global phase claim.
_Avoid_: Hint treated as constraint, constrained success called discovery

**Phase Hint**:
Optional initialization or continuation information that may reorder work but
cannot alter equations, required certificates, or the scientific claim.
_Avoid_: Hidden source bound, fingerprint gate, acceptance override

**Caloric Thermochemistry Boundary**:
The current Engine and Parameter Fit scope retains residual caloric quantities
but has no species ideal-gas caloric model and exposes no total caloric target.
Ideal density and mixing terms required for pressure, fugacity, and phase
equilibrium remain part of the Helmholtz model; they are not a phase-selected
caloric reference. Total internal energy, enthalpy, entropy, Gibbs energy, and
heat capacities require a future explicit EOS-owned thermochemistry design.
_Avoid_: Residual-only EOS, phase-inferred caloric reference, implicit total
caloric property

**Property Capability Ownership**:
EOS owns fixed-state and pure-fluid property values, domains, and exact
selected-parameter derivatives. Equilibrium owns mixture and reactive solved
property values, domains, and exact selected-parameter derivatives. Regression
may discover only capabilities admitted by the current Regression Cutover
Boundary, maps observations, forms scaled residuals, and routes their exact
Jacobians through its reduced-space Ipopt owner.
_Avoid_: Regression property implementation, Regression-owned thermodynamic
derivative, duplicated observation solver

**Regression Cutover Boundary**:
The hard cutover admits EOS-owned fixed-state and Pure Saturation State
Observation Families plus the D-034 fixed-T,P homogeneous reactive and D-035
fixed-topology reactive bubble/VLE families. Every admitted Equilibrium family
owns its descriptor, values, exact total Jacobians, branch certificate, and
typed Non-Evaluable Trial result. D-036 approves GREPE as the next target
architecture but admits no reactive LLE, electrochemical multiphase, reaction-
coefficient fit, or other future capability without its own executable
evidence and receipt.
_Avoid_: Family-specific Regression thermodynamics, reaction fit coordinate,
target architecture presented as current admission

**Non-Evaluable Trial**:
A bounded optimizer point at which a required owner evaluation cannot produce
its declared value, certificate, and exact derivative. The evaluation returns
structured failure without a fabricated residual, Ipopt may backtrack, and the
affected multistart outcome remains in FitResult diagnostics.
_Avoid_: Fixed scientific penalty, silent clipping, whole-fit first-failure abort

**Pure Saturation State**:
The EOS-owned pure-component vapor--liquid boundary at one temperature,
containing the common saturation pressure and both coexisting densities. It is
resolved from pressure and chemical-potential or fugacity equality, not from a
single arbitrary-density pressure evaluation. Any observed subset of pressure,
liquid density or volume, and vapor density or volume may share one solve.
_Avoid_: Fingerprint-approved saturation route, pressure-at-arbitrary-density,
mixture VLE owner

**Pure Parameter Baseline**:
The standard nonassociating pure-component Parameter Fit recipe combining
saturation pressure with saturated-liquid density or molar volume over an
explicit temperature domain. A saturation-pressure-only fit remains valid when
declared, but is reported as a reduced-information fit rather than as the
baseline recipe. Tabulated observations and Published Observation Correlations
are alternative source representations of the same Observation Families.
_Avoid_: Implicit pressure-only baseline, polynomial-as-EOS model, mandatory
homologous-series correlation

**Pure**:
The public EOS model constructed from exactly one selected component. It owns
pure-state evaluation and Pure Saturation State evaluation without a caller
composition. `Pure.state` and `Mixture.state` return the same `State` type.
_Avoid_: One-component Mixture, duplicate PureState, separate pure EOS

**Mixture**:
The public EOS model constructed from at least two selected components. Every
state evaluation requires an explicit normalized composition; mixture VLE,
LLE, flash, and reactive solves remain owned by Equilibrium.
_Avoid_: Implicit composition, pure-fluid saturation owner

**Reference-Dependent Quantity**:
An EOS quantity such as activity, osmotic coefficient, solvation free energy,
or transfer free energy that is valid only under an explicit named solution
reference convention. It does not require ideal-gas caloric thermochemistry.
_Avoid_: Hidden aqueous default, unqualified activity, implicit standard state

**Published Observation Correlation**:
A citation-backed property equation admitted to a closed registry with explicit
coefficients, units, validity domain, outputs, and exact derivatives. It may
materialize observation targets when a real fit requires that published form;
it is not a user callback or an ePC-SAFT parameterization mode.
_Avoid_: Arbitrary expression, uncited formula, homologous-series parameter fit

**Integrated Correlation Objective**:
The normalized continuous-domain squared residual for a Published Observation
Correlation. Preflight increases a deterministic Gaussian quadrature order
until the objective and exact Jacobian converge, then freezes the nodes and
weights for Regression. The fit receipt retains the correlation identity, domain,
measure, quadrature, convergence evidence, predictions, residuals, and Jacobian
rows. Explicit points remain valid for tabulated sources or selected operating
conditions.
_Avoid_: Hidden uniform grid, quadrature adaptation during optimization, unnormalized
row-count weighting, summary-only receipt

**Optimizer Coordinate**:
An affine-scaled physical ePC-SAFT parameter with explicit transformed bounds.
Starts, bounds, active-bound diagnostics, and final values remain expressed as
direct physical parameters; fit files do not select log or logit transforms.
_Avoid_: Hidden parameter transform, unscaled physical coordinate, per-fit
coordinate system

**Fit Prior Boundary**:
The current Parameter Fit has no prior, regularization, or hidden stabilization
residual. Rank deficiency, poor conditioning, and competing basins remain
visible. A future explicit typed prior residual may be designed for a concrete
underidentified scientific case but has no dormant implementation now.
_Avoid_: Automatic weak prior, start-value penalty, implied identifiability

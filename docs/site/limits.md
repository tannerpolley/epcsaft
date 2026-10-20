# Limits

The Engine evaluates the declared model and reports mathematical admission,
values and numerical causes. Accuracy for an application also depends on the
parameter sources, reference choices and physical domain. A converged solve
alone does not validate that model. The
[EOS domain rules](theory/eos.md),
[equilibrium formulation](theory/equilibrium.md) and
[termination algorithm](theory/algorithms.md#alg-termination) define these
boundaries; the [guide](guide/index.md) connects them to runnable cases.

## Source domain and physical validation

Mathematical input admission is distinct from a source's fit range. The Engine
checks required physical inputs, units and the requested mathematical domain;
the caller must establish source applicability and state extrapolation.
Combining published pure-component records does not create a published fitted
mixture parameter set. Unknown uncertainty remains unknown.

The [water notebook](examples/parameters-and-water-state.ipynb) verifies native
pressure closure without a measured-density comparison. The
[ammonia notebook](examples/ammonia-synthesis-equilibrium.ipynb) distinguishes
its low-pressure ideal-limit verification from its extrapolated dense-gas
prediction. The [propane notebook](examples/propane-parameter-estimation.ipynb)
calibrates parameters against observations; calibration is not validation.

## Equilibrium specifications and calorics

Declared solves restrict the phases and supports considered. They perform no
implicit phase-discovery search. Pure stoichiometric solids use declared
formulae and references; this is not a general solid-solution model. The
[halite notebook](examples/halite-solubility.ipynb) supports its one-salt,
one-temperature comparison with separate complementarity and balance checks.

Neutral PH/PS specifications require fixed pressure, free temperature,
declared free phase amounts, one extensive enthalpy or entropy target and
complete consuming thermal records. Charged or source-referenced PH/PS is
unavailable. UV remains deferred. A reacting thermal calculation additionally
requires consistent reaction and caloric standards; supplying contradictory
standards is an error. The
[caloric algorithm](theory/algorithms.md#alg-calorics) and
[declared-equilibrium guide](guide/declared-equilibrium.md) explain admission.
The ammonia notebook is isothermal and does not qualify reacting calorics.

Without ideal thermal records, residual properties can remain available while
Ideal and Total calorics are unavailable. An available endpoint value does not
promise an action across a discontinuous source interval. Preserve the cause
rather than substituting residual or zero values. The
[methane caloric notebook](examples/methane-caloric-properties.ipynb) shows
thermal records and their independent numerical checks.

## Derivatives and support boundaries

Derivative actions describe a smooth fixed branch and declared coordinate
chart. A singular Jacobian, critical/coalesced phases, a singular pressure
root or association problem can make an action unavailable. A central state
may remain useful even when an action does not exist or cannot be resolved.

Zero-fraction insertion and changing support have additional restrictions.
Charged insertion on an ion-free support is not an ordinary smooth action;
entropy insertion and optional polar terms have their own limits. Admitting
an incipient central equilibrium does not qualify every insertion derivative.
Read [ALG-ACTIONS](theory/algorithms.md#alg-actions) and the per-action cause.
The [solved-state notebook](examples/water-isobutanol-solved-state-actions.ipynb)
verifies its stated smooth binary branch, not a derivative through a phase
change. Actions retain observable units; normalize by the physical increments
before reporting partial-derivative units.

## Finite stability searches

Local stability searches are not global proofs. Finite tangent-plane trials,
positive mechanical slopes and local curvature cannot exclude every distant
lower-energy phase. Record unresolved trials, tested supports, work budgets
and estimated uncertainty. Charged trial admission depends on the selected
permittivity model's physical domain; an unavailable trial does not prove
that a phase is physically absent.

The [discovery notebook](examples/water-isobutanol-discovery.ipynb) demonstrates
a neutral finite-search comparison. Its scope and the
[stability algorithm](theory/algorithms.md#alg-stability) do not establish a
global certificate, universal phase enumeration or arbitrary reacting search
accuracy.

## Estimation and process qualification

Least-squares covariance is local and conditional on the supplied error model,
rank and bounds. A local minimax result has separate feasibility and
stationarity diagnostics and no covariance or global-optimality claim.
Measurement reconciliation uses the explicitly declared measurement groups
and uncertainty model. See [fitting](theory/algorithms.md#alg-fit) and
[parameter estimation](guide/parameter-estimation.md).

The [methane ASL notebook](examples/methane-asl-properties.ipynb) checks neutral
pressure-property transport. The
[IDAES settler notebook](examples/idaes-two-liquid-settler.ipynb) checks its
fixed-T/P ionic two-liquid feed mapping and outer process equations. Neither
establishes arbitrary ionic/reacting thermal process derivatives or physical
process validation. Application acceptance requires its own source-backed
inputs and independent evidence for the intended operating domain.

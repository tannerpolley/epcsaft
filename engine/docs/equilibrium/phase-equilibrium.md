# Phase-Equilibrium Documentation Authority

Status: public nonreactive owner and configured two-phase GREPE

`epcsaft.equilibrium.solve` is the only public equilibrium evaluator.
Governance D-030 and the
[CAP-09 promotion receipt](../../../validation/campaigns/2026-khudaida/results/receipts/cap-09-phase-equilibrium-promotion.json)
complete the protected transition defined by
[ADR 0002](../scientific/adr/0002-protect-held-until-phase-equilibrium-promotion.md).
`ARCHITECTURE.yaml` records the machine-readable boundary and the
[general plan](plans/2026-08-07-general-phase-equilibrium.md) records the target
design.

Historical research, source locators, results, and receipts retain their
original HELD, HELD2, and `tp_flash` names as provenance. The HELD/HELD2 native
controller also remains the internal discovery and start-generation algorithm;
the retired name is no longer a Python API.

## Ownership

| Problem | Public owner | Boundary |
| --- | --- | --- |
| Pure-fluid saturation | EOS `Pure.saturation` | Pure vapor--liquid boundary |
| Nonreactive phase equilibrium | Equilibrium `solve` | Fixed-T,P retained discovery or declared 1--3 phases; fixed-T two-phase coexistence |
| General reactive equilibrium | Equilibrium `solve` | Homogeneous, finite--incipient, and configured finite multiphase equilibrium |
| Solid equilibrium | Excluded | Requires another physical model |

## Fixed-T,P feed

`FixedTPFeed` owns positive component amounts at fixed temperature and pressure.
Its `topology` is optional:

- omitted topology runs the retained HELD/HELD2 finite discovery controller;
- `DeclaredTopology` solves exactly one to three ordered all-component phases,
  with at most one vapor, and certifies only that topology.

Discovery retains the incumbent algorithm's implemented binary-neutral and
one-/two-phase charged search boundary and reports `not_established` globality.
A declared solve never falls back to discovery after failure and never turns a
local certificate into a global claim. One phase selects a mechanically
strict-stable pressure root. Two or three phases use the HELD2 Step-8 material,
composition, volume, and phase-fraction formulation.

`PhaseHint` supplies non-authoritative composition and density ordering advice
for an explicitly declared phase index. It may reorder the generated finite
start set but does not populate nonlinear variables. A
`DeclaredPhaseInitialGuess` instead supplies the actual composition, density,
and phase fraction for every phase of a fixed-T,P multiphase topology. Its
phase indices must cover the topology exactly; compositions and phase
fractions are normalized, positive, finite, and representation-valid. Hints
and actual initial guesses are mutually exclusive. Neither can define phase
count or type, change equations or tolerances, or weaken acceptance.

## Fixed-T coexistence

`FixedTCoexistence` fixes one reference-phase composition and returns the
two-phase boundary pressure. Its result uses `IncipientBoundary`; an incipient
boundary has no thermodynamically determined phase fractions. Fixed-T,P feed
results instead use `FeedSplit` with phase amounts, fractions, and total reduced
free energy.

Every `PhaseSpecification` declares a liquid/vapor role and component support.
The reference phase uses `AllComponents`. A `NeutralComponents` partner is
built from every neutral model component and excludes every charged component
structurally. It does not insert zero or epsilon ion coordinates and contains
no water- or source-specific rule. Perdomo Table 3 therefore has a water-only
vapor because water is the model's only neutral component; a multisolvent
electrolyte keeps every neutral solvent.

## Results and certificates

Results retain declared phase order, supported component identities,
composition, molar volume and density, packing fraction, pressure, supported
chemical potentials, mechanical class, parameter fingerprint, route-native
diagnostics, and exact available first-order selected-parameter sensitivity.
Unsupported, singular, ill-conditioned, coalesced, nonconverged, or uncertified
states fail explicitly.

Physical acceptance is independent of the promotion comparison. It requires
the applicable material and charge balances, pressure stationarity, KKT and
bound conditions, phase distinctness, native mechanical stability and ordered
volume-branch identity, and ordinary or charge-constrained modified-potential
equality. Fixed-T coexistence uses all configured neutral-potential equalities or
the charged HELD2 `C-1` modified-potential basis plus common pressure. Exact
implicit sensitivities differentiate those defining equations in the selected
parameter order.

After the declared fixed-T,P Stage-III terminal and pressure update, a state
whose modified-potential mismatch exceeds the existing `1e-6` certificate
enters an interior
exact-derivative Newton correction of the equality-constrained KKT system with
bounded backtracking. The corrected state is used only when its residual
reaches `1e-10`; the complete existing physical certificates are then
reevaluated without relaxed tolerances. Diagnostics report whether correction
was attempted and accepted, its iteration count, and its pre/post KKT
residual. A candidate that reaches phase construction but fails a later role
or modified-potential certificate is retained on
`PhaseEquilibriumError.best_candidate` for diagnosis and has no accepted-state
or globality claim.

## Evidence

The user-approved promotion tolerance is `1e-4` for candidate-versus-incumbent
historical state lineage only. It does not relax conservation, charge,
pressure, KKT, modified-potential, or stability checks. The public installed-
wheel campaign retains two HELD2 claims over all eight Khudaida Figure 2 feeds:
declared two-liquid continuation and omitted-topology discovery. Both compare
the packet-v5 state lineage and physical certificates without asserting private
routes, work counts, predictive cutoffs, or globality.

Engine evidence separately covers discovery amount scaling and parameter
order, strict-stable role selection, neutral and charged one-/two-/three-phase
declared solves, component and declared-phase permutation, non-authoritative
hints, structural neutral-only partner support, and exact implicit
sensitivities against centered re-solves. Perdomo source packets, results, and
receipts remain historical cross-EOS provenance; they are not a current
installed-wheel campaign because the paper used SAFT-gamma-Mie.

The unrelated scratch Gibbs solver remains removed. The private HELD/HELD2
sources and their algorithmic certification evidence remain because they are
the live discovery core, not a second public route. Homogeneous reaction is a
GREPE formulation, not a separate public evaluator.

## General reactive equilibrium

`GeneralReactiveEquilibriumProblem` is the declarative GREPE condition. Fixed
temperature and pressure, conserved totals, phase roles, component support,
reaction eligibility, standard-state records, starts or immutable continuation
anchors, active parameters, and requested outputs determine the compiled
problem. The compiler binds that condition to the selected Parameters and
ordered components.
Application and source names never select an algorithm.

For two or more configured finite phases at fixed temperature and pressure, GREPE
minimizes one total Gibbs objective in positive phase/species amounts and phase
packing variables subject to global conservation and exact per-phase charge
neutrality. The configured phases may use any fluid-role labels. Reaction
affinities and neutral or electrochemical interphase-potential equalities arise
from the minimization KKT conditions. Source-standard-state equilibrium
constants are transformed to the EOS basis before the solve, including
their selected-parameter derivatives. The result describes the configured
local minimum; it does not claim that no unconfigured phase is more stable.

Finite-phase starts are optional. When they are omitted for two configured
phases, GREPE constructs deterministic positive candidates from the declared
reaction feed. Reaction-product/reactant partitions are tried first, followed
by neutral-component partitions and an equal split. Every partition preserves
the total amount of each species; charged species retain an equal split so
each phase remains electroneutral. A fixed-pressure EOS state supplies
each candidate phase volume. Candidates are attempted in stable order until
one passes the complete numerical and physical checks. Explicit starts bypass
this generation and remain the way to select a particular known local branch.

The finite-phase NLP uses logarithmic positive-amount coordinates and a
logarithmic packing coordinate. Its total Gibbs objective is dimensionless
(`G/RT`), while each conservation row is normalized by its conserved total
and feed contribution. IPOPT receives exact gradients and the exact
Lagrangian Hessian. Final acceptance rechecks the original unscaled material,
charge, pressure, reaction-affinity, interphase-potential, KKT, phase-
distinctness, and local-curvature quantities.

For the Ascani 2023 hypothetical case, 33 of 34 deliberately varied feasible
starts reached the same two-phase minimum. Successful starts required 7--61
IPOPT iterations (median 14); the only rejected start drove one configured
phase to effectively zero amount. The automatic reaction-sign candidate
reached the same phase compositions in 11 iterations and 286 EOS
evaluations before callback caching, and 122 afterward, without using the
digitized tie-line compositions.

The same owner returns generic output rows and exact total Jacobians directly
to Regression. Current finite-phase outputs include system or phase pressure,
phase amount, mole fraction, molar density, and dimensionless chemical
potential. For one declared finite and one declared incipient liquid/vapor
phase, solved temperature with fixed pressure compiles to fixed-pressure
coexistence. Either phase can be the fixed-composition reference, so the same
public request represents Bubble-T or Dew-T. The caller supplies the local
temperature bracket and optional start; exact implicit total derivatives flow
through the same output rows and Regression owner receipt. Every configuration
is attempted according to its declared mathematical variables; a missing
formulation or failed solve is returned directly and is never decided by an
application or source name.

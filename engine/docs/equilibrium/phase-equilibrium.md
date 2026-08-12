# Phase-Equilibrium Documentation Authority

Status: protected HELD transition; candidate implemented but not promoted

This document describes the implemented nonreactive phase-equilibrium boundary
and the evidence required before the incumbent may be retired.

## Authority

The active transition decision is
[protect HELD until the phase-equilibrium promotion gate passes](adr/0002-protect-held-until-phase-equilibrium-promotion.md).
The destination remains
[one general nonreactive phase-equilibrium owner](adr/0001-one-general-phase-equilibrium-owner.md),
but its atomic-retirement sequence is superseded. The complete target design is
in the
[General Phase-Equilibrium Plan](plans/2026-08-07-general-phase-equilibrium.md).
`CONTEXT.md` defines shared language and `ARCHITECTURE.yaml` records the current
machine-readable boundary.

Historical research, source locators, results, and receipts keep their original
HELD/HELD2 names as provenance. During this protected transition the same names
also identify the exact incumbent implementation whose evidence remains live.

## Ownership

| Problem | Current route | Boundary |
| --- | --- | --- |
| Pure-fluid saturation | EOS `Pure.saturation` | Pure vapor--liquid boundary and exact selected-parameter sensitivities |
| Nonreactive incumbent | Equilibrium `tp_flash` | Exact retained HELD/HELD2 discovery workflow and certificates |
| Nonreactive candidate | Equilibrium `phase_equilibrium` | Protected discovery, declared fixed-T,P feed equilibrium, and declared fixed-T coexistence |
| Homogeneous chemical equilibrium | Equilibrium `chemical_equilibrium` | One homogeneous reacting state over explicit species, reactions, balances, and charge |
| Reactive multiphase equilibrium | Deferred | Shared state/phase/parameter/sensitivity seam only; no runtime route |
| Solid equilibrium | Excluded | Requires a separate physical model and plan |

`tp_flash` and `phase_equilibrium` remain two explicit interfaces. Neither is
an alias or fallback, and neither silently substitutes the other. Protected
discovery uses the exact incumbent HELD/HELD2 controller. Declared feed and
coexistence specifications reuse its phase coordinates, pressure envelopes,
Step-8 local-equilibrium equations, tolerances, EOS provider, and physical
certificate basis. They do not recreate the removed unrelated Gibbs solver.

## Typed specifications

`FixedTPFeed` owns a positive amount feed at fixed temperature and pressure.
`DiscoverPhases()` runs the exact protected `tp_flash` controller and retains
its admitted one-/two-phase finite-search and `not_established` globality
scope. It does not generalize the incumbent's binary-neutral or two-phase
charged discovery boundary. `DeclaredTopology` instead solves one to three
ordered all-component phases, with at most one vapor, and certifies only that
declared topology. One phase selects a mechanically strict-stable pressure-root
branch. Two or three phases use the HELD2 Step-8 material-balance,
composition, volume, and phase-fraction equations. A local declared solve never
claims global discovery.

The candidate is presently blocked from promotion by Khudaida Figure 2 feed 6:
its first returned phase amount is `0.5365460177251793 mol`, while protected
`tp_flash` returns `0.536547300021282 mol`. The `1.2822961027e-6 mol`
difference exceeds the unchanged `1e-7 mol` HELD2 lineage tolerance. The other
seven feed rows and the additional candidate claims pass on the same installed
wheel, but ADR 0002 requires the complete sweep. This failed comparison is
retained in the CAP-09 blocked-promotion receipt and keeps `tp_flash` live.

`FixedTCoexistence` fixes one reference-phase composition and solves a two-phase
boundary with pressure as an output. Phase fractions are undefined at an
incipient boundary, so the result carries `IncipientBoundary` rather than fake
amounts. A feed result carries `FeedSplit` with extensive amounts, fractions,
and total reduced free energy.

Every `PhaseSpecification` declares its role and support. Fixed-T,P feed solves
currently require `AllComponents`; an unsupported structural request fails
before native execution. In fixed-T coexistence, `NeutralComponents`
constructs the partner phase from all neutral model components and excludes
every charged species structurally; it does not insert zero or epsilon ion
coordinates. Thus the Perdomo Table 3 vapor contains only water because water
is that model's sole neutral component, not because the route recognizes water
or that source. A multisolvent electrolyte keeps every neutral solvent on that
partner support.

`PhaseHint` is indexed by declared phase. In declared fixed-T,P work it reorders
the retained candidate set; in coexistence it reorders the fixed neutral-
enrichment start matrix. It never contributes a new candidate or initializes
an authoritative direct solve. A support-restricted phase hint uses that
phase's own component coordinates. Discovery does not consume hints and
rejects them explicitly. Hints do not change equations, bounds, tolerances,
certificates, selected candidate set, or claim scope.

## Inputs and results

The incumbent retains its exact normalized-feed input, source-domain admission,
typed `PhaseState`, HELD/HELD2 diagnostics, certificates, and sensitivity
availability. These are protected scientific behavior, not a compatibility
facade. In particular, its charged-multiphase sensitivity refusal remains
unchanged.

The candidate accepts one typed specification, optional indexed hints, and an
optional ordered `ActiveParameterSet`. Results preserve declared phase order,
component identities, support, mole fractions, molar volume and density,
packing fraction, pressure, supported chemical potentials, mechanical class,
parameter fingerprint, certificate diagnostics, and exact available
first-order sensitivity. Unsupported, singular, ill-conditioned, coalesced, or
uncertified states fail explicitly.

Candidate diagnostics count Provider callbacks owned by the candidate route.
For protected discovery this is the post-controller native role audit; the
incumbent controller keeps its exact existing HELD/HELD2 timing and evaluation
evidence rather than duplicating those counts into the candidate field.

## Numerical boundary

The restored neutral HELD and strong-electrolyte HELD2 incumbent sources remain
the protected pre-CAP-09 implementation, including the audited source-domain
gate and HELD2 ten-step controller. Candidate-only native entry points add a
neutral ordinary-simplex chart to the HELD2 local phase block without changing
the charged chart. Fixed-T,P feed solves use the incumbent candidate set when
available or generic declared-topology starts otherwise. Hints only reorder
that finite set. Candidate neighborhoods do not become physical bounds.

Role evidence is native mechanical stability plus ordered molar-volume branch
identity. Protected discovery adds the incumbent compressibility branch only
after a native strict-stability audit; an inconclusive gap fails closed. When
only one stable root exists, the same explicit compressibility split is used
solely as a fail-closed one-phase role check. Multiphase acceptance also
requires material and charge balance, pressure stationarity, KKT and bound
conditions, distinct phases, and ordinary or charge-constrained modified
potential equality as applicable.

Fixed-T coexistence solves partner composition and both phase volumes over a
finite, deterministic neutral-component enrichment matrix, then selects the
lowest-pressure certified root found in that declared scope. Neutral support
equates every admitted neutral chemical potential. A charged
all-component partner uses the HELD2 `C-1` modified-potential basis plus common
pressure, not redundant raw single-ion equalities. Exact implicit sensitivities
differentiate these defining equations and preserve selected-parameter order.

The unrelated native CAP-09 scratch phase solver remains removed. The exact
`tp_flash`, HELD/HELD2 diagnostic payloads, certificate audits, and historical
evidence remain protected until promotion.

## Promotion evidence

Engine compact falsifiers cover protected discovery parity, native role
auditing, fail-closed role selection, neutral one- and three-liquid solves, a
six-component multisolvent/multi-ion liquid/vapor root pair under component
permutation, two-neutral structural vapor support, hint-set invariance, and
exact neutral two-/three-phase plus charged model-parameter sensitivities
against centered re-solves. Validation replays all eight accepted Khudaida
Figure 2 feeds, verifies the charged declared split and its exact sensitivity,
and solves a charged fixed-T coexistence boundary in modified coordinates.

The Perdomo installed-wheel campaign fixes the published Table 3 liquid
composition, excludes charged species structurally from the vapor, certifies
the returned ePC-SAFT boundary, and checks its exact sensitivity. The published
2508 Pa remains a descriptive cross-EOS comparison because the paper used
SAFT-gamma-Mie rather than this parameterized ePC-SAFT model.

This evidence supports review; it does not itself retire the incumbent. ADR
0002 still requires successful affected/full checks, the remaining-tree audit,
independent phase-science and solver-architecture acceptance, one promotion
receipt, and explicit user approval before `tp_flash` or the transitional extra
interface may be removed.

`chemical_equilibrium` remains separate. Its arbitrary
species/reaction/balance/charge compiler and solver are unchanged by this
transition.

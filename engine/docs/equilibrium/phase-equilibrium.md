# Phase-Equilibrium Documentation Authority

Status: promoted sole public nonreactive owner

`epcsaft.equilibrium.phase_equilibrium` is the only public nonreactive
phase-equilibrium function. Governance D-030 and the
[CAP-09 promotion receipt](../../../validation/campaigns/2026-khudaida/results/receipts/cap-09-phase-equilibrium-promotion.json)
complete the protected transition defined by
[ADR 0002](adr/0002-protect-held-until-phase-equilibrium-promotion.md).
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
| Nonreactive phase equilibrium | Equilibrium `phase_equilibrium` | Fixed-T,P retained discovery or declared 1--3 phases; fixed-T two-phase coexistence |
| Homogeneous chemical equilibrium | Equilibrium `chemical_equilibrium` | One homogeneous reacting state |
| Reactive multiphase equilibrium | Deferred | No runtime route |
| Solid equilibrium | Excluded | Requires another physical model |

## Fixed-T,P feed

`FixedTPFeed` owns positive component amounts at fixed temperature and pressure.
Its `topology` is optional:

- omitted topology runs the retained HELD/HELD2 finite discovery controller;
- `DeclaredTopology` solves exactly one to three ordered all-component phases,
  with at most one vapor, and certifies only that topology.

Discovery retains the incumbent algorithm's admitted binary-neutral and
one-/two-phase charged search boundary and reports `not_established` globality.
A declared solve never falls back to discovery after failure and never turns a
local certificate into a global claim. One phase selects a mechanically
strict-stable pressure root. Two or three phases use the HELD2 Step-8 material,
composition, volume, and phase-fraction formulation.

`PhaseHint` supplies optional composition and density starts for an explicitly
declared phase index. Hints may reorder the finite start set but cannot define
phase count or type, change equations or tolerances, or alter the accepted
physical result. Hints without a declared topology are rejected.

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
equality. Fixed-T coexistence uses all admitted neutral-potential equalities or
the charged HELD2 `C-1` modified-potential basis plus common pressure. Exact
implicit sensitivities differentiate those defining equations in the selected
parameter order.

## Evidence

The user-approved promotion tolerance is `1e-4` for candidate-versus-incumbent
historical state lineage only. It does not relax conservation, charge,
pressure, KKT, modified-potential, stability, or derivative checks. The
installed-wheel promotion run passed all eight Khudaida Figure 2 feeds, both
selected-parameter derivative checks, and the other retained campaign evidence
(`12 passed, 4 deselected in 437.88s`).

Engine evidence separately covers discovery amount scaling and parameter
order, strict-stable role selection, neutral and charged one-/two-/three-phase
declared solves, component and declared-phase permutation, non-authoritative
hints, structural neutral-only partner support, and exact implicit
sensitivities against centered re-solves. The installed Perdomo campaign fixes
the published Table 3 liquid composition, structurally excludes ions from the
vapor, certifies the ePC-SAFT boundary, and checks exact sensitivity. Its
published pressure remains a descriptive cross-EOS comparison because the
paper used SAFT-gamma-Mie.

The unrelated scratch Gibbs solver remains removed. The private HELD/HELD2
sources and their algorithmic certification evidence remain because they are
the live discovery core, not a second public route. `chemical_equilibrium`
remains a separate homogeneous reacting-state capability.

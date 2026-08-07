# General Phase-Equilibrium Plan

- Status: canonical owner plan; implementation not started
- Date: 2026-08-07
- Owner: Equilibrium
- Parent: [General scientific capability hard cutover](../../plans/2026-08-07-general-scientific-capability-hard-cutover.md)
- Decision: [one general phase-equilibrium owner](../adr/0001-one-general-phase-equilibrium-owner.md)

## Outcome

Replace the separate neutral and electrolyte source-derived controllers with
one nonreactive `phase_equilibrium` owner. It accepts explicit components,
species, charge, feed, thermodynamic specification, and optional phase
information. The same physical algorithm handles neutral and charged systems;
neutral systems naturally skip charge/electrolyte operations. Source cases,
fingerprints, and component names never select runtime behavior.

The first owner supports one to three liquid/vapor phases. Solids are excluded.
Homogeneous chemical equilibrium remains a separate owner. This plan defines a
shared state/phase/parameter/sensitivity seam for future reactive phase
equilibrium but does not implement that coupled problem.

## Public problem modes

### Discovery

The caller does not declare phase count/types. The owner performs stability and
phase finding and may certify one stable phase or a two-/three-phase split using
liquid and vapor phases.

### Constrained specification

The caller declares one to three expected liquid/vapor phase roles. The owner
solves and certifies only that topology. Success does not claim that another
topology is globally absent unless a separate discovery certificate proves it.

### Hints

The caller may provide compositions, densities, phase labels, or continuation
states that reorder initialization/search work. A hint cannot change the
mathematical claim, skip required certification, or turn failure into success.

The result always states the mode, phase roles/count, conserved balances,
stability/globality claim actually established, tolerances, and failure stage.
Constrained success certifies only the declared topology. Discovery reports the
finite candidate/search domain, budgets, exclusions, and completeness criterion
behind its certificate; it may claim globality only when that declared
criterion establishes it. A locally converged or search-exhausted but
uncertified candidate fails closed rather than being promoted to a stable or
global result.

## Current defects to retire

### Priority 0

- Separate neutral HELD1 and electrolyte HELD2 runtime routing.
- Exactly-binary neutral flash data and result names such as `x_methane`.
- `_FLASH_FINGERPRINT` and source T/P/feed bounds.
- Provider active-parameter transport limited to segment diameter and
  dispersion-energy component coordinates while Python accepts broader shapes.
- Parallel provider/native mappings for the same model and parameter identities.

### Priority 1

- Hidden equal-solvent neutral reference construction in native transport.
- Author/paper-named controller, progress, tolerance, and diagnostic strings.
- Fixed `{solvent,cation,anion}` order or selected species tuple.
- Stage names that encode a publication rather than a physical or numerical
  operation.

## Retained physical distinctions

- Neutral and charged formulations may require different algebraic operations,
  coordinates, or constraints inside one owner.
- Specialized numerical algorithms are allowed when selected by typed physical
  capability and problem conditioning.
- A selected cation/anion formula unit remains valid when pair identity and
  stoichiometry are explicit.
- Pure saturation remains EOS-owned with its natural three unknowns.
- The existing homogeneous chemical-equilibrium compiler/solver remains because
  it is already organized around arbitrary species, reactions, balances, and
  charge; it is not rewritten merely to share a facade.

## Canonical problem model

`phase_equilibrium` receives:

- immutable resolved `Parameters` and the exact EOS model topology;
- temperature plus pressure/volume/other supported ensemble specification;
- ordered feed amounts and explicit species/component identities;
- charges and electroneutrality/formula-unit information where applicable;
- mode: discovery or constrained specification;
- zero to three typed phase roles (`liquid` or `vapor`) for constrained mode;
- optional non-authoritative hints; and
- optional `ActiveParameterSet` for exact solved-state sensitivities.

It returns ordered phase states, fractions/amounts, compositions, densities,
chemical potentials/fugacities as applicable, balance and stationarity
residuals, stability/globality certificates, numerical diagnostics, and exact
selected-parameter sensitivities when requested.

## One physical algorithm

The owner constructs the admissible phase/state coordinates from topology and
declared phase roles. It enforces material balance, positivity, normalization,
phase equilibrium, and charge constraints where charge is present. Discovery
adds phase-stability search and candidate management; constrained mode solves
only the declared topology. Charged systems add electroneutrality and the
required modified-potential/reference treatment. Neutral systems omit those
rows because their physical topology has no charge, not because they route to a
different named controller.

Internal stages use physical/numerical names such as `reference_state`,
`stability_search`, `candidate_generation`, `local_equilibrium`,
`globality_check`, and `certificate`. Literature algorithms can justify an
implementation but do not become public symbols.

## Active parameters and sensitivities

The provider receives the same ordered `ActiveParameterSet` resolved by
`Parameters`. It transports component, unordered pair, association,
Born/solvation, permittivity/model, and typed correlation coefficients wherever
EOS advertises partial derivatives. No Equilibrium enum remaps the identities.

For solved equations `F(y, p) = 0`, Equilibrium computes total sensitivities
from EOS-provided exact partials:

```text
dy/dp = -(dF/dy)^(-1) (dF/dp)
```

The first-order contract requires exact `F_y` and `F_p` blocks with declared
row/column identities. Any certificate, KKT, continuation, or exact-Hessian
route must separately declare and obtain every state-state,
state-parameter, and parameter-parameter CppAD partial it consumes. First-order
implicit sensitivity never implies Hessian support; unavailable orders or
nondifferentiable branch changes are typed pre-solve or result outcomes.

The implementation uses a rank-revealing linear solve with conditioning and
singularity diagnostics; it does not differentiate solver history. Phase
appearance/disappearance and nonunique solutions are explicit nondifferentiable
or ill-conditioned outcomes. These sensitivities are made generic now, but no
Equilibrium observation family is admitted to Regression by this plan.

## Implementation sequence

1. **Problem/result types.** Introduce source-neutral phase specification,
   phase role, mode, hint, certificate, and failure types. Remove binary result
   field names and fingerprint gates.
2. **Provider unification.** Consolidate Python/native model mapping and pass the
   canonical active set without family-specific transport slots.
3. **Reference explicitness.** Replace equal-solvent and positional conventions
   with explicit solvent composition, ion pair/stoichiometry, phase, and
   reference convention.
4. **Kernel extraction.** Identify reusable physical operations from current
   neutral/electrolyte implementations; keep specialized numerical kernels only
   behind typed physical applicability.
5. **General owner.** Implement discovery and constrained solves for one to
   three liquid/vapor phases in one controller.
6. **Neutral reduction.** Demonstrate that zero-charge topologies skip charged
   operations and reproduce neutral numerical invariants without a route split.
7. **Solved-state sensitivities.** Add exact implicit sensitivities in canonical
   active-set order with rank/conditioning diagnostics.
8. **Public cutover.** Make `phase_equilibrium` the sole nonreactive public
   owner; remove HELD1/HELD2, fixed flash, source bounds, old schemas, aliases,
   and compatibility tests atomically.
9. **Validation relocation.** Keep named literature and application campaigns
   in Validation against the installed wheel; Engine retains only compact
   algorithmic/numerical evidence.
10. **Future seam review.** Document how homogeneous chemical equilibrium
    states/reactions could join the phase seam later, without adding a dormant
    reactive-phase controller now.

## Numerical verification

Persistent Equilibrium checks must cover:

- one stable phase and two-/three-phase liquid/vapor results;
- discovery and constrained modes with claims distinguished in results;
- bounded discovery certificates recording candidate/search scope, budgets,
  completeness criterion, and any unsearched topology;
- neutral, ternary electrolyte, multisolvent, and multi-ion topologies through
  the same owner;
- arbitrary component order and no fixed solvent/cation/anion positions;
- material, charge, normalization, pressure, and potential residuals;
- phase-label permutation invariance and deterministic canonical result order;
- hints changing work order but not the converged claim;
- active component, several pair, association, and model coefficients where EOS
  derivatives exist;
- exact implicit sensitivities against independent central differences away
  from phase boundaries and singular Jacobians; and
- explicit failure for instability, unavailable topology, singularity,
  nonconvergence, and uncertified globality.

Named publication matrices and broad phase envelopes belong in Validation. A
looser literature acceptance tolerance cannot weaken local balance,
stationarity, derivative, or certificate tolerances.

## Failure semantics

Failure is typed by physical/numerical stage and includes the supplied mode,
declared topology, best candidate, residuals, conditioning, search/certificate
status, and budget exhaustion. A constrained success never implies global
phase discovery. A local optimum never implies global stability. No result is
accepted solely because it resembles a source case. Exhausted search without
the declared completeness criterion is an uncertified outcome, not a global
certificate.

## Issue-ready decomposition

Do not publish these units until the four canonical plans are approved.

| ID | Work unit | Depends on | Acceptance |
| --- | --- | --- | --- |
| EQ-01 | Source-neutral phase problem/result model | None | Modes, roles, hints, and claims are explicit for 1–3 phases |
| EQ-02 | Canonical provider model and active-set transport | EOS-03 | No duplicated parameter enum/mapping or family-limited transport remains |
| EQ-03 | Explicit electrolyte reference inputs | EQ-02, EOS-05 | No positional or equal-solvent assumption remains |
| EQ-04 | General physical kernels | EQ-01, EQ-02 | Kernels select by topology/conditioning, never source identity |
| EQ-05 | Constrained 1–3 phase solve | EQ-04 | Declared liquid/vapor topologies solve and certify only their claim |
| EQ-06 | Discovery 1–3 phase solve | EQ-05 | Stable one-phase and two-/three-phase outcomes carry evidence |
| EQ-07 | Exact solved-state sensitivities | EQ-05, EOS-07 | Active-set ordered Jacobians pass independent differences |
| EQ-08 | Atomic public cutover and dead-surface audit | EQ-06, EQ-07 | Only `phase_equilibrium` remains for nonreactive phase equilibrium |
| EQ-09 | Reactive-phase seam note | EQ-08 | Future state/phase/sensitivity boundary is explicit; no dormant solver exists |

## Done definition

Equilibrium has one source-neutral nonreactive phase owner for one to three
liquid/vapor phases, explicit discovery/constrained semantics, generic active
parameter transport, exact solved-state sensitivity diagnostics, and no
neutral/electrolyte author-route split. Homogeneous chemical equilibrium still
works independently, and reactive phase equilibrium remains consciously
deferred.

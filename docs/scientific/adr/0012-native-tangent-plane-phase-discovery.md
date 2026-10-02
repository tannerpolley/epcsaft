---
status: accepted
date: 2026-10-01
---

# Native tangent-plane testing and fluid phase discovery

Kind: thermodynamic method and computational ownership.

The owner [adopted D1–D10 in #223](https://github.com/tannerpolley/ePC-SAFT/issues/223#issuecomment-5943111636).
The subsequent [Ready review](https://github.com/tannerpolley/ePC-SAFT/issues/223#issuecomment-5943167683)
accepts the design for Build; it is not delivered numerical qualification.
Local reduced-Hessian positivity cannot exclude a remote lower-G fluid split.
Michelsen1982 Part I, Baker–Pierce–Luks1982, Perdomo2025 §3.3 and Ascani2022
§3.3 support tangent-plane testing; the conserved reacting-plane and donor
derivations are owned by [equations.md](../../../engine/docs/equations.md).

## Decision

Use native explicit fixed-T/P tangent-plane trials through the existing exact
Ipopt driver and shared thermodynamics. Python maps physical declarations and
results. [ALG-STABILITY](../../../engine/docs/science/algorithms.md#alg-stability)
owns search, admission, uncertainty, topology changes and work counts;
[formulation](../../../engine/docs/equilibrium/formulation.md#discovery-and-global-claims)
owns the physical meaning. This extends ADR 0009's one-driver requirement.
Only explicit searches own a native work capture. Its recorded calls supply
cost diagnostics; they never determine a thermodynamic or solver decision.

- **D1/D2:** one issue, neutral qualification followed by ionic/reacting
  qualification; measured physical-line counts must fit the current owner
  budgets. Numeric ceilings remain in `ARCHITECTURE.yaml:code_budgets`.
- **D3/D4:** run all enumerated distinct trials. Default supports are union
  Liquid and neutral-only Liquid/Vapor; caller support/label entries replace
  defaults. Charge-neutral reconstruction admits ionic Liquid only. Declaration
  solves, trial invocations and trial rounds have separate finite budgets.
- **D5:** unresolved resolved descent makes discovery unsuccessful, with
  retained diagnostics. Replacement at the balance-rank phase limit and
  generalized multi-donor/reactive withdrawal are separately scoped.
- **D6/D7:** `Discover()` and `stability_test(...)` are explicit only. Test
  pinned/free-intensive declared results at solved T/P with frozen actual
  standards and positive donors. Declared solves, continuation, parameter
  fitting and ASL perform no implicit search; their numerical values and work
  counts must be unchanged. Discover once, then use the discovered declaration.
- **D8:** refined FeOs neutral comparisons are the required independent
  implementation comparison; Clapeyron remains optional. Achieved reference
  refinement is required before composition/density comparisons.
- **D9/D10:** finite conserved descent with propagated estimated uncertainty is
  numerical evidence, distinct from a proved instability witness requiring
  justified total bounds. Finite no-witness search never certifies global
  stability. This is the sole ADR for these method choices.

Normalized reacting trials are search compositions; mole-changing reactions
need not have zero affinity there. Full added declarations restore original
reaction/reference inputs and affinity rows. Frozen plane coefficients retain
log-chart curvature; transformed trial balance multipliers use their complete
right-hand-side contraction. Addition uses a positive donor conserved split;
removal merges or transfers finite inventory using the existing topology events.

## Alternatives and limits

Reject restoration of HELD2/GREPE and its separate Python/SciPy searches.
Reuse retained numerical cases and lessons in `analyses/phase-stability-cutover`,
not the deleted implementation. Reject a separate tangent-plane minimizer and
implicit stability work in known-topology calculations. Random/Powell search,
global-bound methods, rigorous error enclosures, free-intensive discovery and
continuation across phase appearance are not adopted.

The explicit native path is implemented pending numerical qualification and
independent delivered review. Same-wheel comparisons check integration and
topology; FeOs agreement checks independent implementations. Neither constitutes
physical validation or universal stability. The owning issue and retained
analysis hold required evidence and remaining work; this ADR is not a parallel
implementation plan.

The reacting domain diagnosis is retained in the qualification analysis.
Existing pressure-root admission does not classify a sole root as a physical
liquid, and solvent-only permittivity supplies no sourced density exclusion.
A stricter branch-label search would be a new numerical-domain choice rather
than proof of electrolyte-model applicability. That restriction and dielectric
model redesign are not adopted here; the issue's domain addendum owns the
unresolved scientific decision.

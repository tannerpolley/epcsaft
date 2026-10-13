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

The explicit native path has retained numerical qualification on the stated
neutral and admitted density-dependent ionic/reacting cases; the independent
delivered review supports those cases and their stated limits. Same-wheel
comparisons check integration and topology; FeOs agreement checks independent implementations. Neither constitutes
physical validation or universal stability. The owning issue and retained
analysis hold required evidence and remaining work; this ADR is not a parallel
implementation plan.

## Addendum: charged-trial admission by permittivity form

The owner adopted this condition on 2026-10-02 in #223. Charged trial phases in
testing and discovery require the declared bulk-permittivity form to reach
$\epsilon_r\to1$ as $\rho_m\to0$. Schreckenberg `DensityTemperature` and Rueben
`DipolarPerturbation` satisfy the condition by full EOS §6.2 and §6.4;
`IonSuppression` and `ComponentMixing` are density independent. Their electrolyte
terms are defined only for condensed phases, so they receive the typed refusal
"permittivity model has no low-density limit; declare phases" for charged trials.

One native rule in shared trial-support admission serves both explicit entry
points. It rejects a family containing refused charged supports before trial
roots/optimizers and before discovery declaration solves. Caller-supplied
neutral-only families remain available. There is no density cutoff, failed-trial
filter or dielectric substitution. Existing declared equilibrium, continuation,
fitting and ASL remain unchanged. The vacuum limit is necessary for this search,
without establishing physical applicability everywhere it searches.

The reproduced solvent-only descent at 318.15 K, 5918.5286 Pa and
2.3214665 mol/m³ is preserved in the analysis with its calculated decrease and
estimated uncertainty; the current search refuses it with the stated reason.
Case 6's ion-fraction-suppression selection likewise receives a charged-search
refusal while its declared result remains a no-overhead reference. Reacting
qualification uses an explicitly identified Engine-owned density-dependent
selection and manufactured reaction reference, separate from the MEA packet.

Reacting discovery prepares its initial union composition from the existing
positive charge-neutral base before the unchanged ideal-speciation sweeps.
The manufactured reference exposed a zero-fed participant whose generic seed
lost chloride and generated an invalid log coordinate. The bounded correction
retains every required species without changing declared initialization, trial
starts, search budgets, standard terms or tolerances.

The delivered review identified a dimensional defect in the addition acceptance
margin: present-phase molar energy estimates must be multiplied by phase amounts
before comparison with extensive G/RT in mol. The corrected controller retains
that extensive estimate. A 100 mol rescaling of the neutral two-liquid reference
checks the amount-weighted margin and unchanged intensive equilibrium, using
the existing tolerances without a new numerical method.

The 100 mol default-discovery margin has relative error zero against the
amount-weighted reference; the old unweighted expression has relative error
0.99 and fails the unchanged 1e-13 criterion. The same discovery preserves
the neutral intensive split and extensive Gibbs scaling. Cold declared and
supplied-result public reuse at this feed scale then coalesced; that
pre-existing conditioning limit was removed by
[#243](https://github.com/tannerpolley/ePC-SAFT/issues/243) (PR #400). Original qualified
reuse cases and frozen declared work remain unchanged. No driver equilibration
is adopted here, and no universal reuse guarantee is claimed.

## Addendum: explicit cooling-onset dew

The owner adopted this operation on 2026-10-10 in
[#391](https://github.com/tannerpolley/ePC-SAFT/issues/391), after a Ready-with-changes
readiness review. `cooling_onset(compiled)` is a third explicit search beside
`Discover()` and `stability_test(...)`, with free temperature at fixed pressure and
vapor composition. Its principle: the dew is the first liquid to appear on cooling the
vapor. It attempts the ordinary declared solve once, checks an ordinary success with the
liquid trial family at its temperature, and searches after a resolved descent there or a
local-instability or infeasibility failure; [ALG-EQUILIBRIUM](../../../engine/docs/science/algorithms.md#cooling-onset-dew)
owns its procedure, bounds and refusals. D6/D7 stand: declared solves, continuation,
fitting and ASL perform no implicit search, and an empty report list means no trial
family ran. Its coverage is finite: the sampled temperatures and enumerated liquid trials
cannot exclude a narrow negative-distance interval between samples or a missed
composition minimum, and a refusal means unresolved selection, not an absent dew.
It is numerical recovery, not a global-stability certificate.

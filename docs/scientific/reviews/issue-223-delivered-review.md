# #223 delivered review — initial verdict

Read-only reviewer: Claude Opus 5.5, high effort; task `77f29c43-6d91-4695-9fb2-4440a0655b40`, ordinal 1. Inspected commit `798cb53550229f2cd12016b8b57378fb315773da` against `443d828908f4ad5127ab14384f97635dca30bd30`. The following delivered review and correction are retained verbatim. Required dimensional finding accepted; final acceptance remains open until the corrected tree and evidence are rechecked.

# Delivered review of #223 at `798cb535` (base `443d8289`)

**Verdict: Not supported**, because of one required defect. The defect is narrow and does not change any retained qualification number. Once it is corrected and the affected lens rechecked, I expect this to be Supported within the limits below.

I did not edit, build, install or run anything in the Engine. I inspected the source, CSVs, hashes and logs read-only.

## Required defect

**The test for accepting an added phase mixes units** (`engine/native/equilibrium/gibbs_minimizer.cpp:1847–1849`).
- The code accepts an added phase only if the new $G/RT$ is lower than the candidate's by more than `error_estimate`.
- The objective is extensive $G/RT$, in mol.
- `error_estimate` adds up `energy_estimate(q)` for each present phase. Those are per-mole estimates from `phase_energy` (unit-mole composition), and they are not multiplied by phase amounts.
- So the margin has the wrong units. For a 1 mol feed it is conservative. For larger feeds it is too small by roughly the feed total; the 86.7 mol MEA basis is an example of the size applications use.
- This contradicts ALG-STABILITY's "uncertainty-qualified lower Gibbs energy". The finite-split check itself is weighted correctly (`remaining·ε + β·ε + N_p·ε`).
- **Smallest correction:** weight each present-phase estimate by its amount:
  `for (std::size_t p = 0; p < current.present_phases.size(); ++p) error_estimate += best.amounts[p] * energy_estimate(current.present_phases[p]);`
- **Effect on the evidence:**
  - Every qualified addition (cases 1, 3, 4 at 10 MPa, FeOs, case-6 neutral-only) is on a 1 mol feed and was accepted on its first attempt. The margin is never used in the solves.
  - With the fix the margin stays the same or gets smaller. Every accepted addition stays accepted, and no solve's numbers change.
  - Still needed: a rebuild, the gates, a confirming rerun of discovery qualification, and a focused re-review of this branch.

## What the evidence supports

- **Thermodynamics matches the design.** I re-derived these from the source:
  - The donor slope is $d_p = g(w,v_t) - w^T(\mu_p/RT+s) + v_p(P_p-P)/RT$.
  - The present-phase identity is $t = x^T r - r_P$.
  - The ionic reconstruction uses $\sqrt{Q_-/Q_+}$ scaling, then the two-pivot correction.
  - The reduced trial balance basis is the row space of $[e;z]$ on the support.
  - The frozen snapshot standards and the positive-donor rule are as designed.
- **One driver, explicit calls only.** Only `solve_equilibrium` with `request.discovery` and the bound `stability_test` reach the new code. Continuation, fitting, ASL and `solve_compiled` have no caller of it.
- **Permittivity rule.** The four-line guard in `stability_supports` sits after the support, charge and label checks.
  - In `stability_test` it runs before any work.
  - In Discover it runs before the first declaration solve.
  - The owner's exact `ParameterError` reason was observed for ion-fraction suppression, ion-specific suppression and component mixing. Solvent-only was covered through the MEA packet, ion-fraction suppression also through Khudaida.
  - Both eligible forms match their defining equations exactly (6 comparisons, error 0 against 1e-13).
- **Indeterminate trials are not fixable failures.** The 26/4/3 named trials (36 rows including repeats) all have |TPD| ≤ 7.2e-15 and |slope|/error ≤ 0.18. They are trivial roots at the present phase. No trial is "unresolved" or "negative plane without a donor".
- **Integrity checks.**
  - The wheel SHA, the three gate log hashes, all five cost-CSV hashes and the four before-CSV hashes match the committed blobs; the `-text` rule is active.
  - Qualification has 44,836 rows.
  - The render source, included files and HTML hashes match.
  - Acceptance 10 shows `identical: true` for 1,841/838/217/1,478 rows, plus the 396 ASL work rows.
- **Budgets.** `diffstat --budgets` reproduces the issue's figures (Equilibrium 6919/7575, total 21108/22000). The final runtime amendment is exactly eight lines.

## Limits to keep with the verdict

- **What kind of evidence this is.** It is numerical verification of the supplied models and finite trial sets. Descent comes with an estimated uncertainty, not a proved witness. A finite search without a witness is not a global stability certificate. None of this is physical validation.
- **No charged split is qualified.** On an admitted permittivity form, no qualified case has a charged trial producing a resolved descent and an ionic phase addition. The reacting fixture stays one phase.
- **Case-6 neutral-only discovery finds a different, higher-G state.** It returns an ion-free second liquid with $G/RT = -8.662995$ mol. The declared ionic split is $-8.663321$ mol. "No witness" holds only for the restricted trial family; the notebook should state this difference.
- **Zero work for a refused Discover is shown by source, not by the trace.**
  - `discover_phases` keeps its own work trace, so the caller's trace can't see inside it.
  - `compile_problem`, including a Jacobian at a generic point, runs before the guard.
  - The source shows no driver, root, trial, iteration or Hessian work before the guard, but the ~0.1–1 ms wall time is compile work.
  - This is allowed by the owner rule. The "zero evaluations" columns should not be read as "no EOS work".
- **Untested paths and FeOs.**
  - The removal loop inside the controller is not exercised by any discovery case. Only the native fixtures exercise removal.
  - FeOs agreement is between independent implementations, using a stopping test that allows a negative tangent-plane distance.
- **Not portable.** The ASL internal-count evidence is specific to this Linux/glibc build.

## Optional, not required

- Record in the notebook the case-6 neutral-only G difference and that Discover's refusal includes compile work.
- In `tools/check-engine-architecture.py`, the rule against solvers naming problem kinds no longer checks strings, `finite` was dropped from its word list, and the `Capture` read-back check was removed. These look justified, but they are weaker guards.

**Next stage:** Build, for the one-line correction. Then rebuild, rerun the gates and discovery qualification, and send the focused re-review. After Supported, the lead can stop for this issue; no Write is requested.

Correction: the MEA feed basis is about 89.2 mol (CO₂ 2.487 + MEA 9.729 + water 76.975), not 86.7 mol. The conclusion still holds: for large feeds the acceptance margin is too small by roughly the feed total.

## Scope clarification, ordinal 2

+# Scope advice: should #223 fix the 100 mol declared cold start and reuse?

**Classification: this is a separately scoped limitation of declared-phase conditioning. It is not required within #223.** #223 should record it as a named limit on its intended use. It should not repair the shared driver.

## Why it is out of scope

- **No owner case covers it.** Every owner-enumerated falsifier uses 1 mol, except the 2 mol reacting case. The returned-declaration recompile, test and re-solve check was required only for those cases, and it passes. Your 100 mol cold-start and reuse checks are extra probes you added after review, not accepted acceptance criteria.
- **Nothing in #223 caused it.** The declared path is unchanged in source: `run_driver` only gained optional null pointers. All 4,374 frozen rows and the 396 ASL work rows are bit-identical.
  - The cold 100 mol two-liquid declaration ends with Ipopt reporting SUCCESS and rows ≈1e-15, but the declared driver then rejects it as coincident phases, at $G/RT = -54.0584$ mol. That is 100 × the one-liquid value, i.e. both phases sit on the one-liquid composition.
  - This is how the existing declared initialization behaves at that feed scale.
- **Fixing it here would break an owner requirement.** Your untested hypothesis (rescaling the amount coordinate as N/feed) would change the shared declared driver. That would put the mandatory bit-identical values and work counts, and the owner's no-overhead rule, at risk. That is a scope decision for the owner, not a correction within #223.
- **Discovery itself scales correctly at 100 mol.** The default Discover run succeeds with two liquids. Against 100 × the unit-feed reference:
  - composition differs by at most 1.67e-15;
  - amounts by at most 2.42e-13 mol;
  - relative $G$ by 6.70e-15;
  - inventory by 7.1e-15 mol.

  Exercising the corrected addition margin at 100 mol is therefore legitimately #223 work.

## What #223 must still do

1. **Capture the 100 mol Discover margin values** you requested, with no new starts, models or tolerances. These are the operands that verify my required units defect. The current diagnosis JSON stopped before recording the margin, so this is still outstanding.
2. **Bound the reuse claim.** "Discover once, then declare" is qualified only for the original reuse cases. At a 100 mol feed, re-solving the returned declaration through the public declared path coalesces (11 iterations, 18 evaluations, 12 Hessian calls; residual 6.85e-15). It does so even with amount, composition and density carried over. State this as a known limitation of the intended use, caused by the shared declared driver and not proven to be a feed-scale law.
   - Keep the diagnosis files: they are currently only in `/tmp` and need to move into the analysis's results.
   - Do not call the result a physical finding, and do not make it a persistent failure test.
3. **Route the root cause to its owner.** Open a separate issue for declared-driver feed-scale conditioning, carrying the Sol diagnosis (d3, ordinal 5) and the untested amount-coordinate hypothesis. Under the fix-before-negative rule, that issue owns the root-cause attempt. #223 does not need to close it.

## Not a new requirement

The existing `solve_equilibrium(..., start=<returned result>)` warm start may be the reuse route that works at 100 mol. If you already have that evidence, you can bound the claim to it. I am not asking for it.

I'll recheck the complete affected lens once you supply the final frozen tree and the margin values.

Lead disposition: the amount-weighting defect is corrected natively and the
100 mol default discovery margin has been measured. The separate declared
conditioning finding is routed to [#243](https://github.com/tannerpolley/ePC-SAFT/issues/243),
without implementation authority. The supplied-result public reuse route
mentioned as a possibility above was already probed and coalesced; no successful
100 mol reuse route is inferred. Final affected-lens recheck remains pending.

# Permittivity admission readiness

Independent reviewer: GPT-6.1 Sol, xhigh, Codex, 2026-10-02.
Task: `22a45c46-c782-447e-89f7-a46ef55e8306`.
This record is design readiness only.

**Ready. No required design defects found** in the new charged-trial addendum or the related formulation, ALG-STABILITY and ADR edits. Reviewed the live #223 body and all seven comments against `ac7078d8` plus the three doctrine edits.

- **Admission matches the equations and owner decision.** `DensityTemperature` and `DipolarPerturbation` approach εᵣ=1 at vanishing density; `IonSuppression` and `ComponentMixing` lack that vacuum behavior. The edits preserve the condition’s limited meaning: search admission, without claiming physical validation throughout the search domain. [Full EOS §6](/home/tnnrpolley21/.t3/worktrees/ePC-SAFT/t3code-223-build-tpd-stability-test-and-phase-35ebff8f/engine/docs/science/full-eos.md:393), [owner decision](https://github.com/tannerpolley/ePC-SAFT/issues/223#issuecomment-5947638210).
- **One native owner is sufficient.** Both explicit entry points already reach `stability_supports` before search work; discovery calls it before its first declaration solve. Existing `std::invalid_argument` mapping supplies `ParameterError` with **“permittivity model has no low-density limit; declare phases”**. Neutral-only families and declared consumers remain outside the proposed restriction. [Admission owner](/home/tnnrpolley21/.t3/worktrees/ePC-SAFT/t3code-223-build-tpd-stability-test-and-phase-35ebff8f/engine/native/equilibrium/gibbs_minimizer.cpp:1202), [exception mapping](/home/tnnrpolley21/.t3/worktrees/ePC-SAFT/t3code-223-build-tpd-stability-test-and-phase-35ebff8f/engine/native/bindings/module.cpp:19).
- **The reacting references are defined consistently before execution.** Species order, charges, frozen standards, reaction targets and expected amounts agree with the manufactured construction. The first reference satisfies both stated moieties; the pinned variant fixes total phase amount at `100/99` mol through the feed-relative `Pinned(50/99)`. Tolerances are explicit, while phase count and descent remain measured outcomes. [Reference construction](/home/tnnrpolley21/.t3/worktrees/ePC-SAFT/t3code-223-build-tpd-stability-test-and-phase-35ebff8f/engine/tests/native/equilibrium_checks.cpp:2907), [addendum](/tmp/issue223-permittivity-addendum.md:1).

Case 6’s default charged families correctly become refusals, with declared values/work preserved. The retained MEA decrease remains numerical evidence for the original model and inputs, with estimated uncertainty—not a converged physical phase or proved instability. [Retained diagnosis](/home/tnnrpolley21/.t3/worktrees/ePC-SAFT/t3code-223-build-tpd-stability-test-and-phase-35ebff8f/analyses/2026-tpd-phase-discovery/results/diagnostics/phase-domain/diagnosis.md:1).

No Engine calls, edits, builds, installations or Git/GitHub mutations were performed. Runtime qualification, unchanged declared-work verification and independent delivered review remain required.

**Ready remains unchanged; no required defect.** The [Ascani packet](/home/tnnrpolley21/.t3/worktrees/ePC-SAFT/t3code-223-build-tpd-stability-test-and-phase-35ebff8f/data/packets/ascani-2022-case-study-2-electrolyte-lle/3/parameters/parameters.json) matches the supplied SHA-256 and declares `component-permittivity-mixing` at 298.15 K and 1 bar. It supplies water/Na⁺/Cl⁻ permittivity records and all three subset pairs, so the charged refusal check at x=(0.98,0.01,0.01) needs no synthetic parameters.

This remains design readiness only. No code was built, no Engine calls were made, and no files were changed.

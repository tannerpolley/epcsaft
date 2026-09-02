# HELD2 stages I–II: ranked findings

Read against \cite{Pereira2012} Table 1 and §3 (Markdown transcription in Zotero `Y2J6BEAZ`). Files read in full: `held2_step1.cpp` (1126), `held2_step2.cpp` (746), `held2_step3.cpp` (174), `held2_step4.cpp` (257), `held2_step5.cpp` (1984), `held2_step6.cpp` (103), `held2_algorithm.cpp` (372), `held2_certificates.cpp` (474), `held2_tolerances.hpp` (283), `held2.hpp` (273), the step headers, `flash.cpp` (evaluator wiring), `Held2AlgorithmCache` (`held2_step8.hpp`). Measurements: instrumentation branch `perf/instrumentation` at `0db4d46` plus one commit exporting `step_timings` (`perf: export HELD2 per-step timings`); scripts and raw output in `oracle/held2_steps.*`, `oracle/held2_sweep.*`, `oracle/held2_khudaida.*`. Kernel document: `engine/docs/equilibrium/phase-stability.md` on branch `docs/equilibrium-formulation` (PR #194).

Tags follow the CSE audit vocabulary.

## Findings, most consequential first

1. **Stage III dominates unstable feeds and is entered too early.** `replace`. On the Khudaida LLE feeds (5 species, 2 charged) Problem 67 (step 8) takes 42 %, 34 %, 14 % of evaluations (points 1, 4, 8) and the Eq. 12 free-energy test (step 9) rejects 15/16, 6/7, 9/10 of its solutions, sending the loop back to Stage II. Two causes: the candidate gap `ε_b` is `5e-2` (`kHeld2PaperStep6Gap`) where the paper uses `1e-2`, so candidates qualify before the bounds have closed; and Stage III re-solves from scratch each time. Major iterations 38, 31, 18 against the paper's 29 for its ternary LLE, with 4–8× the paper's evaluations (10,287 / 8,092 / 5,366 vs 1,241). Target: paper's `ε_b`, and the Gibbs minimizer on the candidate set (#188).

2. **Stage I never converges; it exhausts a budget.** `replace`. `run_direct_search` sets `maxeval = 6500/2 = 3,250` and DIRECT-L stops only at a witness or at `MAXEVAL_REACHED`. Every stable feed of the aqueous electrolyte (0.1 %–18 % ions; 1 kPa, 3 kPa, 1 bar) costs 3,700–4,630 evaluations and 0.6–0.7 s, and the same one-phase claim. The paper stops after `10·nc` tunneling minimizations. Stage I alone costs 3–8× the paper's whole flash. Target: fixed start list, local `TPD` minimization from the best samples, `search_budget_exhausted` as a named stop.

3. **Pressure envelopes cost about 1,900 evaluations per flash for a root the density owner finds in tens.** `delete`. `evaluate_held2_pressure_envelope` (`held2.cpp`, 450 lines): 64 intervals × 3 points, subdivision to depth 8, up to 120 bisections per bracket; about 270 evaluations per call. Step 3 calls it once per bounding point (`2·d` = 6 calls, 1,635–1,647 evaluations, 16–30 % of each LLE flash) and Step 2 once for the reference plus once per witness refinement. `density.hpp` owns a bracketed Newton root already. Target: reference and bounding points at the density owner's root; Hessians only at accepted roots.

4. **Two electroneutral charts.** `dedupe`. `make_held2_coordinates` + transforms/lift/Jacobian/potentials/cube map/anchor (`held2_step1.cpp`, about 850 lines) and `make_amount_chart` (`chemical_equilibrium.cpp`, registered as `electroneutral_amount_chart`) implement one idea: eliminate one ion by electroneutrality. The modified-fraction chart is the amount chart on `Σ n = 1`. Both neutral-chart defects (#188 review comment; `phase-stability.md` §5) come from the HELD2 chart treating "no charged species" as a sentinel (`eliminated_index == C`) that six consumers must test and two forget. Target: one chart; the neutral case has no eliminated coordinate because it has no charged species.

5. **Step 5 is 1,984 lines for a 150-line inner problem.** `delete` 1,300 of them. `Step5Tnlp` (290 lines) is Ipopt on `g + λ·(x̃⁰ − x̃)` with polytope rows; that is the Gibbs minimizer with one phase and fixed multipliers (`formulation.md` §2). The rest is scaffolding the paper does not have: `certify_step5_terminal` (570 lines) re-derives KKT residuals in physical coordinates and re-verifies the Step-4 cut snapshot, multipliers, and active-cut ids by recomputation; `audit_held2_step5_kkt` (250 lines in `held2_certificates.cpp`) including a Gaussian-elimination rank of the active constraint rows; a dilute-face restart with `mu_init = 1e-16` (100 lines); a 12-iteration Newton pressure polish with a dense solve and line search (170 lines); `same_coordinates` (40 lines) comparing every chart field between three copies of the same struct. The minimizer's own reduced-Hessian and stationarity residuals (`formulation.md` §5) are the certificate.

6. **The paper's start rotation is used only in recovery.** `replace`. `make_step5_start` returns `random_interior` for every start unless `state.step5_requires_new_member`; the paper always alternates random / shifted-previous / near-pure starts and credits the shifted-previous family with up to 50 % cost reduction. Measured: 2.9 inner solves per major iteration at 27 evaluations each (29–31 % of each LLE flash). Target: the three families in rotation from the first start.

7. **Step 2's second search and trace scan are unexercised insurance.** `delete`. After the joint DIRECT-L finds a witness that is not a strict pressure root, a second DIRECT-L over composition alone runs a 65-point volume scan plus up to 64 bisections per trial (`prepare_composition`, 160 lines), and a 129-point volume scan at the trace ion composition precedes both (`kTraceVolumeIntervals`). On every measured feed the first search's witness was accepted after `refine_stationary_witness`; neither the preparation search nor the trace scan produced the witness. Target: witness = local `TPD` minimizer from the best sample, at the density root.

8. **Step 4's LP certificate restates HiGHS.** `delete` 90 lines. `solve_upper_highs` recomputes primal, dual, and complementarity residuals of the LP solution and audits them against `kHeld2Lp*`; HiGHS's `kOptimal` with feasible primal and dual status is the same statement. The LP itself (60 lines of model assembly) follows problem OP exactly, feed cut included, and carries over.

9. **Per-step result structs and reason strings outnumber the decisions.** `dedupe`. Stages I–II produce `Held2Step{1,2,3,4,5,6,7}Result`, `Held2LocalCertificate` (one per Ipopt start, kept in a vector), `Held2Step5KktCertificate` (30 fields), `Held2DiluteFaceRestartEvidence`, `Held2LpCertificate`, `Held2Step4CutEvidence`, `Held2CandidateDecision`, plus a `Held2StepTiming` per invocation, and 30 distinct reason strings. `held2_discovery_to_dict` exports six scalars from all of it. Target: one `Discovery` result (witness, candidates, counts, named stop).

10. **Tolerance drift from the paper.** `prevent`. `ε_b` `1e-2 → 5e-2` (finding 1); composition floor `1e-8 → 1e-10·f_k` (`kHeld2ModifiedLowerScale`), which puts inner-problem solutions closer to the bound than the paper's Step 9 trace treatment assumes; `ε_λ = 0.5`, `ε_x = ε_η = 1e-3` unchanged. 24 of the 40 `Held2Tolerance` entries are consulted in stages I–II; seven of those (`root_*`, `mechanical_margin`, `stable_objective_tie`) are density-root tolerances that belong to the density owner. Target: the paper's values as category-3 entries in the shared numerical policy (#186).

11. **Stage II is not in the gate.** `prevent`. The only gate tests that touch HELD2 are the one-phase electrolyte discovery (`test_discovery.py`, Stage I only) and two synthetic certificate unit tests (`test_internal_certificates.py`, deleted in the gut prototype). Steps 3–10 with the real EOS run only in the Khudaida campaign (`pytest -m campaign`, `tools/run-slice.sh khudaida`). Target: one Khudaida feed in the equilibrium data tests (#190's data-backed set) so the kernel's two-liquid result is a gate claim.

12. **A one-phase verdict at 1 kPa that nobody checks.** `investigate`. At 298.15 K and 1 kPa the 2 % NaCl feed is declared one phase with molar volume 2.477 m³/mol, i.e. a vapour carrying 4 % ions (`RT/P = 2.479`). Whether any liquid within the ion ceiling (0.38) lowers the free energy is not tested by the search's 3,250 samples in a way the result can show, and no test asks. Not a demonstrated defect; a physical claim without a data anchor.

13. **Every HELD2 evaluation records a fresh tape and pays `2n−3` Hessian sweeps.** `replace` (owned by #189). `recordings == evaluations` on every run (3,625 of 3,625; 9,803 of 9,803); Hessian calls per evaluation 3 for `n = 3`, 7 for `n = 5`; association 16 and 42 iterations per evaluation. DIRECT-L, the envelope scans, and the trace scan use only the value and the pressure, so 90 % of Stage I's Hessians are discarded. The memoising `Held2AlgorithmCache` keys on exact coordinates and, by the counters, almost never hits.

14. **Bookkeeping split across files.** `move`. Steps 3, 6, 7 (307 lines over 6 files) and the timing lambdas repeated in every step (`finish` closures, about 40 lines each) fold into the controller.

## Stage table

See `phase-stability.md` §3 for the per-step verdicts (follows / deviates deliberately / drifted) and §4 for the evaluation attribution.

## Lines

Stages I–II plus shared HELD2 code today: `held2_step1` 1,126; `held2_step2` 746; `held2_step3` 174; `held2_step4` 257; `held2_step5` 1,984; `held2_step6` 103; `held2_step7` 30; `held2.cpp/.hpp` 751; `held2_algorithm.cpp/.hpp` 442; `held2_certificates` (Step-5 part) about 300; `held2_tolerances.hpp` 283; step headers 284. Total about 6,480 (of 10,538 HELD2 lines; steps 8–10 with `step8_nlp` are the rest).

Removed by the findings above: chart duplication 850 → 60; envelope 450 → 0; Step 5 scaffolding 1,300 → 0; Step 2 second search and trace scan 330 → 0; LP audit 90 → 0; result structs, timings, reason plumbing about 500 → 60; tolerances 283 → policy entries; step files 3/6/7 → controller. What remains is written fresh at about 1,000 lines plus a 100-line header (`phase-stability.md` §7).

## Number for #188

**1,100 lines** (range 950–1,300) for `phase_stability.cpp` + `.hpp`, replacing the earlier 1,400 and the `investigate` tag. The two open decisions that move it: keeping the DIRECT-L start list (+80 for the NLopt adapter) and reusing the minimizer's state cache (−100).

## What is not established

- Stage attribution for the neutral binary (4,632 evaluations) is HELD v1 and has no per-stage counters; HELD v1 is scheduled for deletion, so it was not instrumented.
- The nine-species MEA model was not run: the feed built from the parameter file failed the electroneutrality check because the charge field is not exposed on `Mixture`; not pursued.
- Cache hit rates of `Held2AlgorithmCache` were not counted directly; the equality `recordings == eos states` is indirect evidence that hits are rare.

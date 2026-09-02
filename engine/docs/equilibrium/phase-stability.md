# Phase stability

Status: kernel document for phase discovery (`find_phases` in the target; HELD2 stages I–II today). It serves §6 of `formulation.md`: the tangent-plane criterion and its Lagrangian dual. Equation identifiers are registered in `engine/docs/equations_registry.yaml`; runtime source is the executable authority.

Source: \cite{Pereira2012}, "The HELD algorithm for multicomponent, multiphase equilibrium calculations with generic equations of state" (HELD: Helmholtz free Energy Lagrangian Dual), Table 1 and §3; the charged-species chart follows Perdomo et al. 2025 (`Perdomo2025`) §3; the dual extremum principle is \cite{MitsosBarton2007}. Where this document says "the paper" it means \cite{Pereira2012}.

## 1. The problem the kernel solves

Given `T`, `P`, and a feed composition `x⁰`, decide whether one phase is stable and, if not, produce a set of candidate phases for the Gibbs minimizer of `formulation.md` §1. The kernel does not itself deliver an equilibrium; it delivers a witness of instability and candidates whose free energies bracket the dual value.

The paper's objects, in the chart of §2:

| paper | meaning | `EqID` |
|---|---|---|
| `g(x̃, v) = (A + P V)/(R T)` | one-phase Gibbs energy in Helmholtz coordinates | `system_gibbs_helmholtz_coordinates` (one phase) |
| `TPD(x̃, v)` (Eq. 9) | tangent-plane distance of a trial from the feed's tangent plane | `tangent_plane_distance` |
| `L(x̃, v, λ) = g + λ·(x̃⁰ − x̃)` (Eq. 6) | Lagrangian of the one-phase minimization with the balance rows dualised | `helmholtz_lagrangian_dual` |
| `G^D = max_λ min_{x̃,v} L` (problem D) | the dual; its value is the equilibrium Gibbs energy, its inner minimizers the stable phases | `helmholtz_lagrangian_dual` |
| inner problem `IP(λ)` | `min_{x̃,v} L(x̃, v, λ)` at fixed `λ`; a lower value `L(λ)` | `dual_inner_problem` |
| outer problem `OP` | LP over the retained set `M`: an upper bound `UBD` and the next `λ` | `outer_approximation_bound` |
| Eq. 7 | at an inner solution `λ_k = μ̃_k − μ̃_d` | `multiplier_potential_link` |
| Eq. 8 | at an inner solution `−∂A/∂V = P` with no explicit row | `pressure_equality` |
| system (PP) | which retained points are candidate stable phases | `candidate_phase_criteria` |

The stability test (Eq. 9) is the inner problem at `λ = ∇g(x̃⁰)` offset by `g(x̃⁰)`; the paper states this equivalence in §3.1. A negative `TPD` anywhere is a proof of instability; a nonnegative minimum is a claim only as strong as the search that produced it. `globality_certificate` stays `not_guaranteed` (crosswalk of 2026-07-31).

## 2. Coordinates

**Neutral phase.** Independent mole fractions `x̃_k`, `k = 1 … C−1`, the last species dependent; volume as `ln v` (the paper uses packing fraction `η`; `ln v` bounded by the packing-fraction limits is the same interval).

**Charged phase.** Eliminate one ion `e` (largest `|q|`, chosen so that every factor below is positive) by electroneutrality and rescale: `EqID: modified_fraction_chart`

```text
x̃_i = f_i x_i,   f_i = 1 − q_i / q_e   (i ≠ e),   Σ_{i≠e} x̃_i = 1,
x_e = −(1/q_e) Σ_{i≠e} q_i x_i,   x_i = x̃_i / f_i.
```

The modified fractions sum to one by construction, so the chart has `C−2` independent coordinates plus one dependent neutral species. Chemical potentials transform contragrediently, `EqID: multiplier_potential_link`:

```text
μ̃_i = (μ_i − (q_i/q_e) μ_e) / f_i,   Σ_{i≠e} μ̃_i x̃_i = Σ_i μ_i x_i,   λ_k = μ̃_k − μ̃_d.
```

so the free energy is unchanged by the chart and the reduced gradient of `g` in the independent coordinates is `μ̃_k − μ̃_d`. The domain is a polytope: `x̃^l f_k ≤ x̃_k ≤ f_k · q_opp/(|q_k| + q_opp)` for ions (`q_opp` the largest opposite charge), the dependent species `≥ x̃^l f_d`, the eliminated ion nonnegative, and the EOS ion ceiling `Σ_ions x_i ≤ x_ion,max` when the model declares one.

**One chart, not two.** The target keeps one electroneutral chart. `formulation.md` §2 registers `electroneutral_amount_chart` (amounts, one ion eliminated, used by the Gibbs minimizer); the modified-fraction chart above is that chart on the normalised simplex `Σ n = 1`. The kernel therefore uses the amount chart with one normalisation row and needs no modified factors, no `eliminated_index` sentinel, and no radial map from a unit cube onto the polytope (§5, defect 1). Today `make_held2_coordinates` (`held2_step1.cpp`) and `make_amount_chart` (`chemical_equilibrium.cpp`) are two implementations of one idea.

## 3. Stages and steps: paper against implementation

Paper Table 1 has nine steps in three stages. The implementation numbers its steps 1–10 and the mapping is not one-to-one. For each paper step: the current owner, whether the implementation follows the paper, deviates deliberately, or has drifted, and what the target keeps.

| paper step | current owner | verdict | target |
|---|---|---|---|
| Stage I, Step 1: stability test, tunneling from `10·nc` random starts, stop at first negative `TPD` | `run_held2_step1` (chart, domain, feed lift; no paper counterpart) and `run_held2_step2` (`held2_step2.cpp`): reference pressure envelope at the feed, a 129-point trace-composition scan, DIRECT-L over `(ln v, cube)` with `maxeval = 6500/2`, a second "preparation" DIRECT-L over composition alone with a 65-point volume scan per trial, then re-certification and stationary refinement of the witness | **Deviates deliberately** on the search: DIRECT-L global sampling replaces tunneling (`tunneling_transform` is Documentation-only). **Drifted** on termination: the paper stops after `10·nc` local minimizations; the implementation stops only at a witness or when the fixed sample budget is exhausted, so a stable feed always costs the full budget (§4). The `TPD` itself (`evaluate_held2_tpd`) follows Eq. 9 exactly. | `search_tangent_plane`: fixed start list, `TPD` at the feed's density root, local minimization from the best samples, stop at the first witness or after the start list; named budget stop. |
| Stage I, Step 2: `k = 0`, `UBD = G^P`, initialise `M` per Pereira et al. 2010 | `run_held2_step3` (`held2_step3.cpp`): `M = {feed} ∪ {2·d bounding points, each half-way to the polytope boundary along one coordinate, at its own pressure root} ∪ {stage-I witness}` | **Follows** the 2010 strategy the paper cites. **Drifted** on cost: every bounding point pays a full pressure envelope (about 270 evaluations), 1,640 evaluations per flash before the first outer problem (§4). | Same set; each point at the density owner's root (tens of evaluations, not hundreds). |
| Stage II, Step 3: outer problem `OP` | `run_held2_step4` (`held2_step4.cpp`): HiGHS LP, one cut per `M` point plus the feed cut `v ≤ G^P`; primal, dual, and complementarity residuals recomputed and audited | **Follows** exactly. The residual audit (about 90 lines) restates HiGHS's own optimality status. | `outer_approximation_bound`: the LP and its multipliers; HiGHS status is the certificate. |
| Stage II, Step 4: inner problem `IP(λ)` by multistart until `L ≤ UBD`; starts: random, previous solution shifted by 0.01, near a pure component; alternate mole and mass fractions | `run_held2_step5` (`held2_step5.cpp`, 1,984 lines): Ipopt on `L` in `(x̃, ln v)` with the polytope rows; start families `random_interior`, `boundary_aware`, `shifted_retained`, `space_filling` (Halton); stop at the first start whose value is `≤ UBD + 1e-8` and whose point is new in `M` | **Follows** the acceptance rule and the start families. **Deviates**: the structured families are used only in recovery mode (`step5_requires_new_member`); in the nominal pass every start is random. **Scaffolding**: `certify_step5_terminal` and `audit_held2_step5_kkt` (about 820 lines) re-derive KKT residuals in physical coordinates, verify the step-4 cut snapshot and multipliers by re-computation, and add a dilute-face restart and a Newton pressure polish, none of which the paper has and all of which the Gibbs minimizer's own stationarity check covers (`formulation.md` §5). | `minimize_lagrangian`: the Gibbs minimizer with one phase, fixed `λ`, no balance rows, the paper's three start families always in rotation. |
| Stage II, Step 5: system (PP) | `run_held2_step6` (`held2_step6.cpp`): gap, derivative agreement away from the lower bound, pairwise distinctness in packing fraction or composition | **Follows**. **Drifted** tolerance: gap `ε_b` is `5e-2` (`kHeld2PaperStep6Gap`) against the paper's `1e-2`; `ε_λ = 0.5`, `ε_x = ε_η = 1e-3` as in the paper; the composition floor is `1e-10·f_k` against the paper's `1e-8`. | `distinct_candidates` with the paper's values as category-3 tolerances. |
| Stage II, Step 6: `k ← k+1`, go to Step 3 | `run_held2_step7` (30 lines) and the loop in `run_held2_algorithm` | **Follows**; cap 128 major iterations. | controller loop with a named stop. |
| Stage III, Steps 7–9 | `run_held2_step8` (Problem 67, `held2_step8_nlp.cpp`), `run_held2_step9` (Eq. 12–13 tests), `run_held2_step10` (trace refinement) | Out of scope here; measured in §4 because it drives the loop count. | Replaced by the Gibbs minimizer on the candidate set (`formulation.md` §1, issue #188). |

Two things the paper has that the implementation does not: the alternation between mole-fraction and mass-fraction coordinates in Step 4, and the `nc_pp` cap on converged components for `nc > 6`. Two things the implementation has that the paper does not: the pressure envelope (a scan of 64 log-volume intervals subdivided to depth 8 with up to 120 bisections per bracket, `held2.cpp`) used to pick the lowest-free-energy mechanically stable root, and the memoising `Held2AlgorithmCache` of states, volume bounds, and packing fractions keyed on exact coordinates.

## 4. Evaluation accounting

Measured on the instrumentation branch (`perf/instrumentation`, counters in `instrumentation.hpp`; `OMP_NUM_THREADS=1`; one EOS evaluation is one `evaluate_phase_state` call, which records a fresh tape and takes `2n−3` Hessian sweeps for `n` species).

**Stable feed, three species (water, Na⁺, Cl⁻; 298.15 K, 1 bar).** Outcome `one_phase_no_negative_witness_detected` for every feed from 0.1 % to 18 % ion fraction and at 1 kPa and 3 kPa:

| item | evaluations |
|---|---|
| DIRECT-L samples (`maxeval = 6500/2`) | 3,250 |
| trace-composition volume scan | 129 |
| reference pressure envelope | about 270 |
| volume-bound calls (composition changes) | 980 |
| total (`eos_evaluation_count`) | 3,700–4,630 |

The count does not depend on the feed: the search ends only at the budget. Wall 0.6–0.7 s; 16 association iterations per evaluation.

**Unstable feeds, five species (Khudaida water + ethanol + isobutanol + NaCl LLE, 293.15 K, 1 bar; `analyses/reproduction/2026-khudaida`, Figure 2 points 1, 4, 8).** Two liquids accepted at every point:

| step | point 1 | point 4 | point 8 | what the evaluations buy |
|---|---|---|---|---|
| 2 stability | 1,250 (12 %) | 1,207 (15 %) | 1,282 (24 %) | witness found after about 400 DIRECT samples; the rest is the reference envelope, the trace scan, and witness refinement |
| 3 initial `M` | 1,639 (16 %) | 1,647 (20 %) | 1,635 (30 %) | six bounding points × one envelope each |
| 4 outer LP | 0 | 0 | 0 | 38, 31, 18 solves |
| 5 inner problem | 2,994 (29 %) | 2,452 (30 %) | 1,682 (31 %) | 109, 91, 56 Ipopt solves; 27 evaluations per solve; 2.9 solves per major iteration |
| 6 candidates | 38 | 34 | 17 | packing fractions |
| 8 Problem 67 | 4,361 (42 %) | 2,747 (34 %) | 745 (14 %) | 16, 7, 10 completed solves; the Step-9 free-energy test (Eq. 12) rejected 15 of 16, 6 of 7, 9 of 10 |
| total | 10,287 | 8,092 | 5,366 | major iterations 38, 31, 18; wall 17.5 s, 13.0 s, 8.7 s; 42 association iterations per evaluation |

The paper's Table 2 reports 1,241 evaluations and 29 major iterations on average for its ternary LLE (system IVa, SAFT-VR) and 481 for the binary. The implementation spends 4–8× that on a five-species LLE. Attribution of the excess, in order: Stage III entered on candidates that the loosened gap admits early and rejected by Eq. 12 almost every time (34–42 %); pressure envelopes in Steps 2 and 3 (about 1,900 per flash, 18–35 %); inner-problem solves at 27 evaluations each with Hessians on every callback (29–31 %).

The neutral binary discovery in the gate (4,632 evaluations) runs HELD v1, not this kernel, and has no per-stage counters; it is not attributed here.

## 5. The two neutral-chart defects, in context

The neutral chart is HELD2's degenerate case: `charged_count == 0` in `make_held2_coordinates` sets `eliminated_index = C` as a sentinel and every consumer must remember to test it. Four do (`held2_lift_modified_fractions`, `held2_physical_fraction_jacobian`, `held2_transform_modified_potentials`, `evaluate_held2_phase_block`); two do not.

1. `polytope_anchor` (`held2_step1.cpp:143-187`) reserves `ion_total = 0.5` when the ion ceiling is undefined and gives the neutral species `(1 − 0.5)/neutral_count`; with no ions the anchor sums to 0.5 and `held2_transform_physical_fractions` throws at `:498`. The anchor exists only to map DIRECT's unit cube onto the polytope radially (`held2_map_unit_cube_to_independent_fractions`). The paper has no such map: it samples `X` directly. With the amount chart on the normalised simplex the map disappears and so does the anchor.
2. `physical_total_ion_coefficients` (`held2_step1.cpp:128`) and `solved_state_sensitivity.cpp:235` index `charges[eliminated_index]` without the sentinel test; both are reached only when the ceiling is defined, which no neutral model declares, so the read is latent. The fix is not a guard: it is a chart whose neutral case has no eliminated ion because it has no charged species, which is the shape `make_amount_chart` already has.

Both are consequences of building the neutral chart as a special case of the charged one instead of the charged chart as the neutral chart plus one eliminated coordinate.

## 6. Named stops and tolerances

The kernel reports one of: `witness_found` (Stage I negative `TPD`), `search_budget_exhausted` (Stage I start list consumed with no witness; the one-phase claim is only as strong as the list), `candidate_set` (Stage II produced ≥ 2 candidates), `major_iteration_cap`, `inner_problem_start_cap`, `evaluation_failure` (named EOS failure). Today these are spread over `Held2Step2Outcome`, `Held2Step5Result.reason`, `Held2Step6Result.reason`, and `Held2AlgorithmResult.failure_stage/failure_reason` with 30 distinct reason strings across stages I–II.

Tolerances the kernel needs (category numbers per `formulation.md` and issue #186): `TPD` negativity margin (3), inner acceptance gap `L ≤ UBD + ε` (3), (PP) gap `ε_b`, derivative `ε_λ`, distinctness `ε_x`, `ε_η` (3, the paper's values), composition floor `x̃^l` (2), LP feasibility (Ipopt/HiGHS solver targets, 4). Today `held2_tolerances.hpp` names 40 `Held2Tolerance` entries of which stages I–II consult 24; the seven root-finding and topology entries (`root_*`, `mechanical_margin`, `stable_objective_tie`) belong to the density owner.

## 7. What the target keeps

Load-bearing, rewritten small: the chart of §2 as the amount chart plus normalisation; `TPD` with the reference at the feed's density root; a fixed start list and local `TPD` minimization; the LP with the feed cut; the inner problem as the Gibbs minimizer with fixed multipliers and the paper's three start families; the (PP) filter with the paper's tolerances; `M` retention with representation-equivalence de-duplication; a state cache keyed on coordinates.

Carried over close to verbatim: `evaluate_held2_tpd` (30 lines), the HiGHS model assembly in `solve_upper_highs` (60), `representation_equivalent` and `retain_held2_m_point` (40), the Halton start generator (45), the (PP) loop (60).

Not carried: the pressure envelope, the unit-cube radial map and anchor, the trace-composition scan, the preparation search, `certify_step5_terminal` and the Step-5 KKT audit, the dilute-face restart, the Newton pressure polish, the cut-snapshot re-verification, per-step result structs and timings, and steps 8–10 (issue #188 hands candidates to the Gibbs minimizer).

Budget for `phase_stability.cpp` (fresh lines): chart and polytope 60, `TPD` and search 230, LP 120, inner minimization 150, candidates 60, controller with cache, stops, and hand-off 370 → about 1,000 lines plus a 100-line header; range 950–1,300 depending on whether the DIRECT-L start list is kept (issue #188, owner decision) and whether the state cache is reused from the minimizer.

## 8. Registry entries

Registered by this document: `modified_fraction_chart` (`held2_transform_physical_fractions`), `multiplier_potential_link` (`held2_transform_modified_potentials`), `outer_approximation_bound` (`run_held2_step4`), `candidate_phase_criteria` (`run_held2_step6`), `helmholtz_lagrangian_dual` (Documentation-only; realized by the Step 4/5 alternation), `tunneling_transform` (Documentation-only; not implemented). Already registered by `formulation.md`: `tangent_plane_distance` (`evaluate_held2_tpd`), `dual_inner_problem` (`run_held2_step5`), `pressure_equality`, `electroneutral_amount_chart`. The rewrite moves the refs to `phase_stability.cpp`; the identifiers do not change.

# Performance floor, derived and measured

Measured 2026-09-02 on branch `perf/instrumentation` (`f0abace`, `main @9eccbfd` plus counters only), wheel `287ad813144acc856381d64ca43177e061531f28cec605b7a036a3943b34f57b`, `OMP_NUM_THREADS=1`, gcc 13.3.0, AMD Ryzen 5 5500, Python 3.13, `ldd` in `floor-measurements.json`. Script: `measure_floor.py`. Counters: tape recordings and recording time, tape operations (`ADFun::size_op`), Jacobian/Hessian/Forward calls, sweep time, association solves and iterations, exposed as `_core._instrumentation()`.

## 0. The premise, tested first

The estimate this replaces assumed re-taping dominated the 33 ms reactive evaluation. It does not.

| case | wall | recording | sweeps (Jacobian + Hessians + Forward) | Hessian calls | assoc. iterations |
|---|---|---|---|---|---|
| neutral binary, `Mixture.state(rho)` | 660 µs | 17.8 µs (2.7%) | 49 µs (7.5%) | 1 | 0 |
| associating binary, `Mixture.state(rho)` | 1480 µs | 37 µs (2.5%) | 791 µs (53%) | 5 | 126 |
| aqueous electrolyte (3), `Mixture.state(rho)` | 1006 µs | 31 µs (3.1%) | 397 µs (39%) | 3 | 14 |
| nine-species electrolyte, `Mixture.state(rho)` | 8494 µs | 100 µs (1.2%) | 7315 µs (86%) | 7 | 44 |
| nine-species, GREPE per-evaluation path (`evaluate_parameterized_phase`) | 4815 µs | 124 µs (2.6%) | 4572 µs (95%) | 8 (two tapes) | 44 |
| reactive bubble point, whole solve (4 cases) | 10.9–11.4 s | 160–180 ms (1.5%) | 10.6–11.0 s (97.5%) | 11,041–11,547 | 69,649–74,020 over ~1,600 solves |
| neutral binary stability decision | 671 ms | 68 ms (10%) | ~185 ms (est., small-tape sites not timed) | 4,632 | 0 |

Two consequences. First, one cached tape per topology is worth 1–10% of an evaluation, not the factor the estimate assumed; the premise was wrong. Second, the reactive solve reports 183–357 "EOS evaluations" but performs 1,700–1,830 tape recordings and 11,000–11,500 Hessian sweeps: about 10 tapes and 63 dense Hessians per counted evaluation. The cost is the Hessian count, and the count is set by two multipliers stacked: `(S + 1)` Hessians per tape to assemble the implicit association second derivative, and roughly 10 tapes per counted evaluation from re-evaluating the same state across the objective, gradient, Hessian, vapour-closure, and derivative-audit callbacks.

The `Mixture.state` path adds about 590 µs of Python and pybind marshalling on top of 67 µs of native work for the neutral binary; that is the 50× the earlier note attributed correctly to marshalling.

## 1. Arithmetic floor for one residual Helmholtz value

The operation count is measured, not estimated: CppAD records one tape operation per elementary scalar operation, and `size_op` reports it.

| case | tape ops (value + mass-action residuals) | transcendental ops (counted from the kernels) | association fixed-point iterations | double-precision floor |
|---|---|---|---|---|
| n = 2 neutral | 301 | 2 `exp` (diameters), 4 `sqrt` (pair energies), 3 `log` (hard sphere, 2 chain contacts) | 0 | 0.4–0.8 µs |
| n = 2 associating, S = 4 | 601 | as above + 16 `exp` (site-pair strengths) + 4 `log` | 126 (damped substitution, 0.5) | 1.5–4 µs incl. the solve |
| n = 3 electrolyte, S = 2 | 572 | + 1 `sqrt` (κ), 2 `log` (χ), permittivity correlation | 14 | 1–2 µs |
| n = 9 electrolyte, S = 6 | 2164 | 9 `exp`, 81 `sqrt`, 11 `log`, 36 `exp` (site pairs), 6 `log` (χ) | 44 (1,584 site-pair MACs per iteration) | 5–12 µs |

Scaling read off the kernels (`helmholtz.hpp:215-357`): the pair loop is `n²` with a `sqrt` per pair and about 20 flops; moments and chain terms are `n` with a `log` per component; the seven-term integrals and `C1` are constant (~120 flops); association strengths are `S²` with an `exp` each and the residual is `S²` per iteration; Debye–Hückel is a `sqrt` plus a `log` per ion; Born is `O(n_ion)`.

Throughput assumptions (line items):

| assumption | value used | sensitivity |
|---|---|---|
| scalar double throughput, compiled, no vectorisation (the loops are short and data-dependent) | 1 flop / 0.5–1 ns | linear in the floor |
| transcendental cost (`exp`, `log`, `sqrt` via libm) | 15–40 ns each | the 9-species value has ~140 transcendentals: 2–6 µs of its 5–12 µs; dominant at n = 9 |
| cache | everything (parameters, ~50 KB tape) resident in L1/L2 | if not, ×2–3 on the upper end |
| association solve | Newton with `H_X` reaches 1e-13 in 3–5 iterations from the previous point (today 14–126 damped-substitution iterations) | today's 126 iterations at S = 4 are ~15 µs on their own; the floor assumes ≤5 |

CppAD's interpreted forward sweep runs the same 2164 ops at 14 ns per op measured (below), so a value through the tape costs ~30 µs against 5–12 µs compiled.

## 2. Derivative multiplier

Measured per-op cost on the sweep path: nine species, 7 Hessians in 7,315 µs over a 2164-op tape of width `w = 2 + n + S = 17` (`T, ρ, x, X_sites`). A dense CppAD Hessian is one forward-0 plus `w` (forward-1 + reverse-2) pairs, `1 + 2w = 35` sweeps; 7 × 35 = 245 sweeps in 7,315 µs → 30 µs per sweep → **14 ns per tape op per sweep**. Recording: 100 µs / 2164 ops → 46 ns per op recorded. Neutral binary (301 ops, w = 4): 9 sweeps ≈ 38 µs, consistent with the 49 µs measured for Jacobian + Hessian + Forward.

Relative to one value sweep (30 µs at n = 9):

| quantity | sweeps | multiple of the value | today at n = 9, S = 6 |
|---|---|---|---|
| value | 1 | 1 | 30 µs |
| gradient (reverse, one dependent) | 1 + 1 | 2 | ~60 µs; today the Jacobian has `1 + S` dependents |
| dense Hessian, `w` independents | `1 + 2w` | 23 at `w = 11` (sites eliminated), 35 at `w = 17` | 1.0 ms each, 7 of them: 7.3 ms |
| parameter Jacobian, P active | included in the same Hessian when parameters are tape inputs | +2P sweeps | — |

The implicit association second derivative today: `(S + 1)` dense Hessians of the full tape, then `X_qr` assembled by linear solves. With one Hessian of the Lagrangian `a + λ·H` (λ from one linear solve of `H_X`) the same tensor costs one Hessian. Multiplier: `(S + 1)` → 1 (7× at S = 6, 5× at S = 4), and if the sites leave the tape as independents (they are eliminated by the Lagrangian), `w` drops from 17 to 11: 35 → 23 sweeps, another 1.5×. Combined: **245 sweeps → 23 sweeps per derivative-rich evaluation, about 10×**, from 7.3 ms to ~0.7 ms at the same 14 ns/op.

Floor for value + gradient + Hessian + P-parameter Jacobian, nine species:

- **0.7 ms** with the tape interpreted by CppAD as today (23 + 2P sweeps × 30 µs);
- **60–100 µs** if the tape is compiled to machine code (CppADCodeGen or equivalent), where a sweep runs at ~1.5 ns per op;
- textbook "3–5× the value" applies to the gradient only; a dense Hessian of `w` variables is `O(w)` sweeps and cannot be brought below that without sparsity, and this tape is dense.

## 3. Tape cost, and whether caching is the whole story

Measured: recording is 1.2–3.1% of a single evaluation, 1.5% of the reactive solve, 10% of the stability decision (small tapes, cheap sweeps). Caching one tape per topology removes at most that. The premise that caching is the main lever is falsified by these numbers. The main levers, in order:

1. Hessian count per evaluation, `(S + 1)` → 1 (7×), and site elimination (1.5×).
2. Tapes per counted evaluation, ~10 → 1, by evaluating each state once per Ipopt iteration and serving objective, gradient, constraint Jacobian, and Hessian callbacks from it (today each callback and the incipient-vapour closure re-evaluate). Worth ~5–10× on the reactive solve on its own, with no change to the derivative code.
3. Evaluation count in discovery (section 4).
4. Recording (≤10%).
5. Python marshalling on the `Mixture.state` path only (590 µs per call), irrelevant to native solves.

## 4. Evaluation-count floor for a binary stability decision

A tangent-plane stability decision for a binary at fixed T, P in the `(x, v)` space: minimise `TPD(x, v)` from a start set that covers every basin the EOS can have. For a binary there are at most two composition basins per density branch; a defensible start set is the two pure-component corners on each of the liquid and vapour branches plus the feed on the other branch, 4–6 starts. Newton with the exact 2×2 Hessian from the tape converges in 6–12 iterations from those starts. The dual bound loop adds one LP (no EOS work) per major iteration and 1–3 major iterations for a one-phase conclusion. Floor: **40–150 evaluations** (lower end: 4 starts × 8 Newton steps plus 8 for bound verification; upper end: 6 starts × 12 steps, 3 major iterations, and a 20-point coarse sweep to seed the starts).

Measured: 4,631 evaluations, 42 optimizer solves. The excess is spent on (a) the DIRECT-L global search with a work budget of 6,500 (`held2_step2.cpp`, `step2_search_work_budget`), which samples the domain rather than descending from a start set, and (b) the pressure-envelope scan (64 log-volume intervals × 3 points, subdivided to depth 8, up to 120 bisections per bracket) used to find density roots at each trial composition, where a bracketed Newton from the previous root needs 3–6. Both are algorithmic, not implementation, choices. Each of those 4,631 evaluations also records a fresh 301-op tape and takes a dense Hessian it never uses for a root scan.

## 5. Whole-problem floor

Composition: evaluations × per-evaluation cost + KKT linear algebra. A 30-variable KKT system factorises in ~9,000 flops (≈10 µs); over 30 Ipopt iterations that is 0.3 ms and never sets the floor.

**(a) Nine-species reactive bubble point** (finite electrolyte liquid, incipient neutral vapour, two reactions). Unknowns: 9 amounts + 1 packing coordinate for the liquid, 3 for the incipient vapour, P; with multipliers ~25–30. Ipopt iterations: 15–30 with an exact Hessian (today the finite solve reaches its answer in 11–61 iterations per the GREPE plan; the bubble closure adds an outer Newton of 5–18 steps today, which folds into the same system when P is a row). Evaluations: 100–200 liquid evaluations with Hessian plus 50–100 vapour evaluations with gradient.

| endpoint | evaluations | per evaluation | total |
|---|---|---|---|
| lower: compiled tape, one Hessian, 100 + 50 evaluations, 15 Ipopt iterations | 150 | 60–100 µs liquid, 10 µs vapour | **10–15 ms** |
| upper: interpreted CppAD tape cached, one Hessian, 200 + 100 evaluations, 30 iterations | 300 | 0.7 ms liquid, 0.1 ms vapour | **150–220 ms** |
| today | 183–357 counted, 1,700–1,830 tapes, 11,000–11,500 Hessians | 33 ms per counted evaluation | 7.7–11.4 s |

The range is 50–1,000× below today. The upper end needs only the two derivative changes (one Hessian, one evaluation per state per iteration); the lower end needs a compiled tape.

**(b) Neutral binary stability decision.**

| endpoint | evaluations | per evaluation | total |
|---|---|---|---|
| lower: 40 evaluations, compiled tape, Hessian only at Newton steps (w = 4, 9 sweeps at 1.5 ns/op) | 40 | ~5 µs | **0.2 ms** |
| upper: 150 evaluations, interpreted tape (38 µs Hessian + 5 µs value), plus density roots by bracketed Newton at ~4 value-only sweeps | 150 | ~60 µs | **9 ms** |
| today | 4,631 | 145 µs incl. recording and an unused Hessian | 671 ms |

## 6. Assumptions and what falsifies them

| assumption | used | falsified by |
|---|---|---|
| 14 ns per tape op per sweep (interpreted) | measured at n = 9 and n = 2 | a sweep counter that shows more sweeps per Hessian than `1 + 2w` (CppAD may use more for `Hessian(x, l)` with many dependents); measure with `size_forward_bool`/timing per call |
| 1.5 ns per op compiled | literature for CppADCodeGen on this size of tape; not measured here | compile one 2164-op tape and time a sweep; if > 4 ns/op the lower endpoints double |
| `(S + 1)` → 1 Hessian yields the identical tensor | algebra of the implicit function theorem | the #189 equivalence table (elementwise 1e-10) |
| 10 tapes per counted evaluation are redundant re-evaluations of one state | inferred from 1,759 recordings / 183 evaluations | instrument the callback sequence: if the extra tapes are distinct states (line search, vapour closure) the 5–10× from state caching shrinks to 2–3× |
| Newton association in ≤ 5 iterations | standard for a well-scaled mass-action system from a warm start | the iteration counter after the change; today 14–126 |
| 40–150 evaluations suffice for a binary stability decision | algorithm analysis, section 4 | the eight LLE feeds and the neutral test discovering the wrong phase count from the reduced start set |
| Ipopt 15–30 iterations for the bubble point | GREPE plan's measured 7–61 with median 14, plus the pressure row | the iteration counter on the eight cases after the boundary owner lands |
| KKT algebra negligible | 30³/3 flops | a KKT larger than ~200 variables, which no admitted problem has |

## 7. Acceptance criteria, as written into #189 and #188

**#189 (Helmholtz tape and Hessian), quantity: derivative-rich evaluation cost.** Case: the nine-species MEA model (`analyses/2026-mea-parameter-estimation/inputs/neutral-selected-parameters.json`) at 313.15 K, ρ = 3.0e4 mol/m³, the composition in `measure_floor.py`, through the native per-evaluation entry with one active parameter. Threshold: value + gradient + Hessian + parameter Jacobian in ≤ 1.0 ms (interpreted tape; the upper endpoint with 30% margin) and ≤ 23 + 2P sweeps and ≤ 1 Hessian call per evaluation on the counters; the associating binary at 376.88 K in ≤ 5 association iterations per solve. Measurement: `OMP_NUM_THREADS=1 EPCSAFT_WORKTREE=<worktree> python measure_floor.py <out>` (the script in `oracle/`; the counters must remain available as `_core._instrumentation()`), median of 20 calls after 3 warm-ups. Environment record: wheel SHA-256, commit, `ldd` of `_core*.so`, compiler, CPU model, `OMP_NUM_THREADS`, as `measure_floor.py` writes them. A compiled-tape implementation is judged against 100 µs on the same case.

**#188 (phase stability), quantity: evaluations and wall time per stability decision.** Case: the neutral binary fixture (`analyses/2026-fixed-pressure-anchor/inputs/fixtures/neutral_binary.json`) at 300 K, 1e5 Pa, feed (0.5, 0.5), which is one phase. Threshold: ≤ 150 EOS evaluations and ≤ 150 tape recordings per decision (today 4,631), wall ≤ 10 ms on the interpreted tape; the eight water/ethanol/isobutanol/NaCl feeds each in ≤ 600 evaluations. Measurement: `discovery_case()` in `measure_floor.py` and the eight-feed comparison, same environment record. A decision that meets the count but changes the phase count or compositions beyond the Table 3 uncertainties fails.

**Reactive bubble point (composite, #187 + #189):** the four `vle_obs_*` cases in `measure_floor.py` in ≤ 250 ms each (upper endpoint) with ≤ 300 recordings and ≤ 300 Hessian calls per solve.

## 8. Does the floor justify the greenfield rebuild?

On performance alone, no. The measurements say that three bounded changes inside the existing tree reach most of the gain: one Lagrangian Hessian instead of `(S + 1)` (7×), one state evaluation per Ipopt iteration shared across callbacks instead of ~10 (5–10× on the reactive solve, pending the callback instrumentation named above), and Newton-from-starts instead of a 6,500-work global sample plus envelope scan in discovery (30–100× on evaluation count). Tape caching, the change the plan had at the front, is worth 1–10%. Together these take the bubble point from 11 s to the 150–220 ms upper endpoint and the stability decision from 671 ms to under 10 ms without a new tree; the remaining 10× to the lower endpoints needs a compiled tape, which is orthogonal to the rebuild. The rebuild is still justified by what the audit found elsewhere: six problem classes, fifteen definitions of R, the Python numerics, and 46k removable lines. #189 should therefore be reordered so the Hessian change and the per-iteration state cache come first and tape caching last, and its acceptance is the table above, not "one cached tape".

---

# Revision 2 (2026-09-02): every public entry point

Sweep of 39 entry points on `perf/instrumentation` `0db4d46` (adds tape width, range, and density-iteration counters), wheel `324eac4f…`, `OMP_NUM_THREADS=1`, gcc 13.3.0, Ryzen 5 5500. Full table with every counter: `entry-point-profile.md` / `.json`; script `profile_entry_points.py`. Counters are reset before each call; state calls are medians of 5–20, solves single calls after one warm-up.

Counter gap, stated: sweep time is measured only at the two full Helmholtz tapes (`thermodynamic_derivatives`, `fixed_topology_parameterized_mixture_phase_block`); the five small tapes (pure phase, packing, inverse packing, mixture phase) report recordings and Hessian calls but not sweep time, so "sweeps 0%" for saturation, discovery, and the declared routes means unmeasured, not zero. Ipopt iteration counts for saturation and discovery are not exposed. Python marshalling is obtained only by subtraction.

## The four mechanisms, per call

Shares are of wall time; `rec` = tape recordings; `H/rec` = Hessian calls per recording; `assoc/solve` = successive-substitution iterations per association solve; `dens` = density-root iterations.

| entry point | wall ms | rec | rec % | sweep % | H/rec | assoc/solve | dens | dominant mechanism |
|---|---|---|---|---|---|---|---|---|
| `Mixture.state` ρ, neutral | 0.46 | 1 | 3 | 8 | 1 | 0 | 0 | Python marshalling (~85% of wall) |
| `Mixture.state` ρ, associating | 0.99 | 1 | 3 | 53 | 5 | 126 | 0 | (S+1) Hessians; association 126 iterations |
| `Mixture.state` ρ, electrolyte 3 | 0.69 | 1 | 3 | 39 | 3 | 14 | 0 | Hessians |
| `Mixture.state` ρ, nine species | 5.25 | 1 | 1 | 87 | 7 | 44 | 0 | (S+1) Hessians, width 17 |
| `Mixture.state` P liquid, neutral | 7.5 | 158 | 19 | 67 | 1 | 0 | 141 | Hessian per density iteration; 141 iterations; re-taping 19% |
| `Mixture.state` P vapor, neutral | 13.4 | 321 | 21 | 70 | 1 | 0 | 304 | 360-interval scan; re-taping 21% |
| `Mixture.state` P liquid, associating | 42.5 | 80 | 3 | 89 | 5 | 381 | 64 | 5 Hessians × 64 iterations; 30,515 association iterations |
| `Mixture.state` P vapor, associating | 150 | 321 | 4 | 95 | 5 | 27 | 304 | scan × 5 Hessians |
| `Mixture.state` P liquid, nine species | 382 | 77 | 2 | 98 | 7 | 44 | 61 | 7 Hessians of width 17 per density iteration |
| `Mixture.state` P vapor, nine species | 1,691 | 322 | 2 | 98 | 7 | 28 | 306 | scan × 7 Hessians |
| `Mixture.reference`, electrolyte 3 | 108 | 387 | 6 | 92 | 3 | 13 | 307 | four fixed-P states, each a scan |
| `ElectrolyteReference.state` | 47.6 | 158 | 6 | 88 | 3 | 13 | 124 | one fixed-P state |
| `observable_block` fixed-state, neutral | 0.72 | 2 | 4 | 11 | 1 | 0 | 0 | marshalling |
| `observable_block` fixed-pressure-fugacity, neutral | 7.1 | 159 | 18 | 64 | 1 | 0 | 141 | as fixed-P state |
| `Pure.saturation`, propane | 454 | 5,422 | 11 | 15 (partial) | 1 | 0 | 0 | 21 Ipopt starts × ~258 callbacks; ~75% outside tape work |
| `Pure.saturation` with start | 1.3 | 6 | 8 | 8 | 1 | 0 | 0 | already near floor |
| `Pure.saturation`, associating pure | 531 | 1,700 | 4 | 24 (partial) | 3 | 22 | 0 | starts × Hessians |
| declared fixed-T,P one phase | 4.6 | 126 | 24 | n/m | 1 | 0 | 0 | 48-interval pressure scan; re-taping 24% |
| declared fixed-T,P three phase | 10.1 | 45 | 7 | n/m | 1 | 0 | 0 | near floor for its route |
| fixed-T coexistence, associating | 98 | 240 | 5 | 90 | 5 | 212 | 0 | 5 Hessians; 50,840 association iterations |
| fixed-P coexistence, associating | 1,559 (7,166 on the first run; variance unexplained) | 4,082 | 5 | 90 | 5 | 211 | 0 | outer Python bisection × inner solve; 863,025 association iterations |
| fixed-T coexistence, charged 4-comp | 63 | 260 | 8 | 87 | 3 | 23 | 0 | Hessians |
| discovery, neutral binary | 478 | 4,632 | 10 | n/m | 1 | 0 | 0 | evaluation count (DIRECT-L 6,500 budget + envelope) |
| discovery, charged electrolyte 3 | 870 | 3,626 | 8 | 87 | 3 | 16 | 0 | evaluation count × 3 Hessians |
| GREPE finite LLE, ternary | 99 | 1,430 | 17 | 67 | 1 | 0 | 1,028 | 1,028 density iterations inside the start construction |
| GREPE finite + incipient, solved P | 438 | 2,213 | 7 | 78 | 2.3 | 20 | 527 | Hessians; density iterations in vapour closure |
| GREPE homogeneous reaction | 13 | 0 | 0 | 0 | — | — | 0 | Python only (no EOS evaluation) |
| GREPE solved-T bubble, binary | 1,678 | 4,082 | 5 | 90 | 5 | 211 | 0 | same path as fixed-P coexistence |
| GREPE reactive bubble, nine species | 7,416 | 1,759 | 2 | 97 | 6.5 | 46 | 572 | 11,485 Hessians of width 16 |
| `regression.least_squares`, saturation 3 params | 11,705 | 140,562 | 11 | 15 (partial) | 1 | 0 | 0 | 26 cold saturation solves × 5,422 |
| `regression.curve_fit`, same | 880 | 10,848 | 10 | 15 (partial) | 1 | 0 | 0 | 2 cold solves |
| `regression.least_squares`, equilibrium owner | 227 | 795 | 5 | 52 | 1.9 | 15 | 0 | 8 solved-P evaluations |

**Mechanism 1, re-taping per iteration.** Recordings per unit of useful work: 158–322 per fixed-pressure state, 387 per reference, 5,422 per cold saturation, 4,632 per binary stability decision, 4,082 per fixed-P boundary, 140,562 per three-parameter saturation fit. Share of wall: 1–4% on wide tapes (nine species, associating), 6–8% on the electrolyte-3 tape, 18–25% on the small neutral tapes (width 2–4). So tape caching is worth up to a quarter on the neutral fixed-P, declared, and discovery paths and 1–3% on the wide-tape paths. It moves up for the EOS entry points but is nowhere the largest term.

**Mechanism 2, Hessians where a first derivative suffices.** The shared entry is `NativeModel::evaluate` → `thermodynamic_derivatives` (`derivatives.hpp:62`), which always builds the Jacobian and `(1+S)` dense Hessians. `solve_density` (`density.hpp:109`) and `polished_pressure_density` (`eos.cpp:105`) call it once per iteration to read `P` and `dP/dρ`, so every fixed-pressure state pays 61–306 full Hessian sets: 67–98% of its wall. A value-plus-`dP/dρ` path costs one forward-0 sweep, one reverse-1 sweep for the gradient, and one forward-1/reverse-2 pair for `∂²a/∂ρ²`: 4 sweeps against `1 + 2w` per Hessian times `(1+S)`; at width 17 with S = 6 that is 4 against 245, 60×; at width 4 neutral, 4 against 9, 2×. The same entry serves `Pure.saturation` through `pure_phase_block` with an extra third-order sweep per callback.

**Mechanism 3, association convergence.** Iterations per solve today: 126 (associating state at fixed ρ), 381 (fixed-P liquid associating), 212 (coexistence), 44–46 (nine species), 13–27 (electrolyte 3), 22 (associating pure saturation). Distribution: the fixed-P associating liquid ran 30,515 iterations over 80 solves; the fixed-P boundary 863,025 over 4,082. Newton on the site fractions with `H_X` converges in 3–6 from the previous state (the residual is 2e-13). Cost share: at S = 4 an iteration is ~150 flops, so 863k iterations are ~0.1–0.2 s of the boundary's 1.6–7 s (5–10%); on `Mixture.state` at fixed ρ the 126 iterations are ~20 µs of 0.99 ms. This mechanism is real, isolated, and second-order in wall time; it is being fixed separately and is assumed landed below.

**Mechanism 4, evaluation count per solve.** Saturation: 5,422 recordings for a three-variable problem, from 21 Ipopt starts (7 vapour × 3 liquid seeds, `pure_saturation.cpp:317-330`) each running ~258 callbacks (max_iter 50, no callback caching); with a start it takes 6. Floor without a start: 2–3 seeds from a Clausius–Clapeyron estimate × 6–8 Newton steps = 15–25. Discovery: 4,632 against 40–150 (section 4 above). Fixed-P boundary: 4,082 recordings from a Python temperature bisection (64 iterations cap) around an inner fixed-T solve; with T as an unknown row, one Newton of 10–15 steps over two phases: ~30 evaluations. Fixed-pressure states: 141–306 density iterations from the 360-interval scan; with an anchor or a bracket from the previous root, 5–8 (the `density_anchor` path already does this: `Pure.saturation` with a start is 6 recordings). GREPE finite LLE: 1,028 density iterations in start construction against ~40. Regression on saturation: 140,562 recordings for 26 residual evaluations because each is a cold 21-start solve; with continuation, 26 × 6.

## Per-entry-point floors

Lower endpoint: compiled tape (1.5 ns/op/sweep), one Hessian, Newton counts; upper endpoint: interpreted CppAD tape cached (14 ns/op/sweep), one Hessian, Newton counts; both after association Newton. Marshalling floor for a `py::class_` binding is taken as 20–50 µs per call.

| entry point | today | floor (compiled … interpreted) | what sets the endpoints |
|---|---|---|---|
| `Mixture.state` ρ, neutral | 0.46 ms | 30–110 µs | marshalling 20–50 µs + 9 sweeps × 301 ops |
| `Mixture.state` ρ, nine species | 5.25 ms | 0.1–0.7 ms | one Hessian at width 11 (23 sweeps × 2164 ops) |
| `Mixture.state` P, neutral, with anchor | 7.5 ms | 0.1–0.2 ms | 6 iterations × 4 sweeps |
| `Mixture.state` P, neutral, no anchor | 7.5–13 ms | 0.3–0.5 ms | 20–40 value-only scan points + Newton |
| `Mixture.state` P, nine species | 0.38–1.7 s | 1–5 ms | 6–40 iterations × (value + 3 sweeps) |
| `Mixture.reference` | 108 ms | 2–5 ms | four fixed-P states with anchors |
| `Pure.saturation`, cold | 454 ms | 0.3–1 ms | 15–25 evaluations, 3-variable Newton, no Ipopt |
| `Pure.saturation`, associating pure | 531 ms | 1–3 ms | as above with one Hessian |
| declared fixed-T,P three phase | 10 ms | 2–5 ms | 20–40 evaluations at width 4 |
| fixed-T boundary, associating | 98 ms | 2–5 ms | 10–15 Newton steps × 2 phases, one Hessian |
| fixed-P boundary | 1.6–7.2 s | 3–8 ms | T as an unknown row; no outer bisection |
| discovery, neutral binary | 478 ms | 0.2–9 ms | 40–150 evaluations |
| discovery, charged 3 | 870 ms | 15–30 ms | 150 evaluations × one Hessian at width 6 |
| GREPE finite LLE, ternary | 99 ms | 3–6 ms | 30–60 evaluations × 0.1 ms; anchored density |
| GREPE solved P, 4-comp charged | 438 ms | 20–50 ms | 150 evaluations × 0.3 ms |
| GREPE reactive bubble, nine species | 7.4 s | 10–15 ms … 150–220 ms | section 5 |
| `regression.least_squares`, saturation | 11.7 s | 30–50 ms | 26 evaluations × warm saturation |
| `regression.least_squares`, equilibrium owner | 227 ms | 10–20 ms | 8 evaluations × 2 ms |

## Revised ordering for #189 (per family, not one list)

- **EOS fixed-pressure states, reference, and every route that solves density** (`Mixture.state(P)`, `observable_block` fixed-pressure, `Mixture.reference`, GREPE start construction): (1) a derivative entry with an order argument so density iterations get value and `dP/dρ` from 4 sweeps, not `(1+S)` Hessians (2–60×); (2) Newton from an anchor or the previous root instead of the 360-interval scan (20–50× on iterations); (3) tape caching (10–25% on width ≤ 4, 1–3% otherwise); (4) one Hessian where a Hessian is wanted.
- **Saturation and saturation fits**: (1) warm continuation seeds and a 3-variable Newton in place of 21 Ipopt starts (5,422 → 15–25 evaluations; the fit from 140,562 to ~160); (2) drop the third-order sweep from callbacks that do not need it; (3) tape caching (11%).
- **Associating and electrolyte states, coexistence, discovery (charged)**: (1) one Lagrangian Hessian instead of `(1+S)` (3–7×); (2) evaluation count (discovery 30–100×; boundaries 100×); (3) tape caching (4–10%).
- **Reactive bubble point**: unchanged from revision 1: Lagrangian Hessian, one state evaluation per Ipopt iteration shared across callbacks, tape caching last.

## Acceptance criteria (revised; posted on #189)

Held to on `profile_entry_points.py` with the environment record it writes, medians of the stated repetitions, after association Newton lands:

| case | quantity | threshold |
|---|---|---|
| `Mixture.state` P liquid, nine species, with a density anchor | recordings ≤ 8, Hessian calls ≤ 1, wall ≤ 5 ms | today 77 / 539 / 382 ms |
| `Mixture.state` P vapor, neutral, no anchor | recordings ≤ 40, wall ≤ 0.5 ms | today 321 / 13 ms |
| `Mixture.reference`, electrolyte 3 | recordings ≤ 40, wall ≤ 5 ms | today 387 / 108 ms |
| `Pure.saturation`, propane 300 K, cold | recordings ≤ 30, wall ≤ 1 ms | today 5,422 / 454 ms |
| `regression.least_squares`, saturation, 3 parameters | recordings ≤ 300, wall ≤ 50 ms | today 140,562 / 11.7 s |
| `Mixture.state` ρ, nine species | Hessian calls = 1, sweeps ≤ 23 + 2P, wall ≤ 1 ms | today 7 / 245 / 5.3 ms |
| GREPE reactive bubble, nine species | ≤ 300 recordings, ≤ 300 Hessian calls, ≤ 250 ms | today 1,759 / 11,485 / 7.4 s |
| discovery, neutral binary | ≤ 150 recordings, ≤ 10 ms (#188) | today 4,632 / 478 ms |

## Verdict, from the complete measurement

Unchanged in conclusion, sharpened in detail. Every entry point is dominated by one of four bounded defects inside the current tree: a derivative entry that always builds `(1+S)` Hessians (fixed-pressure states, reference, GREPE start construction, saturation callbacks), a 21-start Ipopt saturation with no continuation (saturation and every fit on it), evaluation-count policies in discovery and the fixed-P boundary, and `(1+S)` Hessians per evaluation in the reactive solve. Tape caching is the fifth lever and worth 1–25% depending on tape width; it is not the plan's first item on any path. None of the four requires a new tree; the rebuild stands on ownership and size, and #189's acceptance is now the per-entry table above rather than one number.

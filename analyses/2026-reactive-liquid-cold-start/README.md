# Cold starts of the reacting MEA liquid (#154, #157)

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The completed
cold-start diagnosis depends on the original MEA and absorber checkouts and pinned
wheels; its frozen application loaders are excluded from current declaration checks.

Question: why did cold one-liquid solves of the nine-species reacting MEA liquid exhaust Ipopt's
200-iteration budget, and what fixes them?
- #154: Böttinger 042/046 and Matin 016 fail at the pre-refit record.
- #157: about 1 % of absorber states fail, and a 1 Pa pressure change flips the outcome.

The two issues share one cause.

**Inputs.** All application clones are read only (`scripts/cases.py`).
- **Absorber (#157).** Clone of MEA-Absorption-Column `4f2d12e9` (the #148 merge), with record `868a5018…`, its Engine reactions and neutral reference, and its ideal-gas records.
  - The problem is declared as `EnginePhase._problem` declares it: one liquid, an `Amounts` feed of apparent CO2/MEA/water, and no composition guess.
- **MEA (#154).** Checkout `calibration-misfit` at `e522dc4`, through #152's loader `analyses/2026-mea-reactive-fit-timing/scripts/run.py`, which uses MEA's own request mapping, its packet guesses, the pre-refit record and refit C's start.
- **Wheels.** Every JSON file records the SHA-256 of its wheel.
  - Before: `684b213a…` (main `62f136ab`).
  - After: `4105e2a4…` (this branch).

## Diagnosis

The cold start seeds the one liquid from the caller's guess or from the feed.
- An unguessed ion gets 1/|support|, so each ion starts at 6.7 % in the absorber.
- A caller's zero water ions are raised to the same value.

An Ipopt file log of the stalled absorber state (329.39 K, 110900 Pa) shows the optimizer wandering through compositions far from equilibrium:
- HCO3⁻ = 1.0 at iteration 27;
- CO3²⁻ = 1.0 at iteration 31;
- iterations 44–72 in a basin where H3O⁺ ≈ 0.028 carries the positive charge while MEAH⁺ sits on the 1e-8 box.

The log also shows repeated restoration phases, inertia corrections up to 1e4, and near-feasible full steps of size one that never contract.

The mechanism is the log-composition chart:
- The objective gradient for species i is (N x_i/N_feed)(μ_i − (Bᵀλ)_i).
- An under-populated species therefore has almost no restoring force, even when its chemical-potential defect is ten or more.
- Its Lagrangian Hessian diagonal, N x_i (1 + defect), is negative when the defect is below −1.

The equilibria themselves are unique. Every route that converges agrees to 1e-12 relative (below), so the cause is the start, not a hard state.

## Fix (algorithms.md ALG-EQUILIBRIUM)

A reacting declaration with one free-amount phase now starts from the ideal-solution speciation of the feed balances. This also covers #152's finite-phase sub-solve.
- **Frozen terms.** Each chemical row's residual potential, −ln v and standard term are frozen at the seed.
- **Iteration.** The start minimizes Σ n_j (c_j + ln n_j/n) subject to Bn = b with the Gordon–McBride element-potential Newton iteration (NASA RP-1311).
- **Sweeps.** Two sweeps refresh c and the density root at the new composition. The last sweep's root seeds the start's packing, so the start costs one density-root solve more than before.
- **Warm starts.** A previous result whose phase `start_point` copies skips the speciation.

The density root must follow the speciation. A variant that held the molar volume at the unspeciated feed's root reintroduced absorber stalls in a probe that was not retained. So each sweep continues the previous root, and solves it anew when speciation moves it beyond the anchored window.

## Results

**Reproduction states** (`states-{before,after}.csv`). The other route is:
- for #157, native continuation from a converged absorber state (318.15 K, loading 0.2);
- for #154, the solve anchored on the other parameter set's solution.

| State | before: cold (iterations) | after: cold (iterations) | after: cold vs other route |
|---|---|---|---|
| #157 329.39 K, 110900 Pa | stall, residual 6.2e-4 (200) | converged (13) | ≤ 5.9e-13 relative |
| #157 329.39 K, 110901 Pa | converged (145) | converged (13) | " |
| #157 318.15 K, 111100 / 111110 / 111200 Pa | converged (58 / 140 / 58) | converged (13–14) | " |
| #157 318.15 K, 110900 Pa, CO2 2.550636 | converged (154) | converged (13) | " |
| #154 Böttinger 042 / 046, Matin 016, pre-refit | stall (200) on all three | converged (13–15) | " |
| #154 the same states, refit C | converged (24 / 16 / 19) | converged (13–15) | " |

The issue's feeds are rounded. At those values the first row reproduces its stall and residual exactly; the other five rows converge before the fix, but only after 58–154 iterations, which is the same slow phase.

**#157 random absorber sweep** (`sweep-{before,after}.csv`): seeds 1 and 2, 350 states, 310–362 K, 1.05–1.12 bar, loading 0.03–0.56.

| | failures | iterations: median / p90 / max |
|---|---|---|
| before | 3 (budget) | 22 / 35 / 200 |
| after | 0 | 13 / 14 / 17 |

**123-state MEA packet timing** (`timing-*.csv`), #152's protocol, in seconds (median / p90 / max). Passes were interleaved on one host that was running an unrelated game (load average 1.9–5.6), so compare within this table first.

| | cold solve | anchored (trial protocol) | MEA forked child | cold failures |
|---|---|---|---|---|
| #152 table (quiet host) | 0.26 / 0.36 / 2.2 | 0.26 / 0.31 / 0.37 | 0.38 / 0.42 / 0.45 | 3 |
| before, eight passes | 0.254–0.275 / 0.352–0.387 / 2.14–2.48 | 0.253–0.273 / 0.295–0.327 / 0.326–0.354 | 0.360–0.386 / 0.402–0.437 / 0.434–0.638 | 3 |
| after, two passes | 0.243–0.247 / 0.281–0.283 / 0.302–0.337 | 0.256–0.260 / 0.294–0.298 / 0.312–0.335 | 0.364–0.368 / 0.403–0.405 / 0.439–0.440 | 0 |

- **Cold solves.** Faster at every statistic, and the three budget failures are gone. Median cold iterations on the 44 speciation states fall from 47.5 to 13.
- **Anchored solves.**
  - Iterations are unchanged on 115 of 123 states: anchored guesses are already speciated. The three #154 states save 3–9 iterations; five others move by one.
  - The paired per-state median difference against the interleaved pass `before-8` is +0.2 to +1.3 ms. The anchored median rises 0.3–2 %; p90 falls in both passes and the maximum in one.
  - Every statistic is within #152's table.

**Persistent falsifier.** `engine/tests/equilibrium/test_trace_reliability.py::test_reacting_liquid_cold_start_recovers_manufactured_speciation`.
- **Setup.** Packet-3 EOS parameters (identical to `868a5018…`), with Böttinger 046 and Matin 016 as manufactured speciations. ln K_r = Σ ν_ri (ln(x_i ρ) + μ_i^res/RT) at x, so x is exact by construction. The feed and unspeciated guess are MEA's own.
- **After.** The test recovers x to rtol 1e-8.
- **Before.** On 684b213a both states exhaust the budget.

## Limits

- **Evidence coverage.** Evidence covers the MEA nine-species liquid only. Other one-phase reacting declarations, for example in Lithium, now start the same way. Two-phase declarations, including Lithium's reactive liquid–liquid system, are unchanged.
- **Start sensitivity.** Ipopt stays sensitive to start quality on this system: the frozen-volume variant (probe not retained) stalled again. The fix gives an accurate start; it does not make the barrier path in the log chart robust from a poor one. That is a separate question for the driver.
- **Anchored-solve cost.** Anchored one-phase solves run one extra density-root solve and two row evaluations; the paired per-state median cost is +0.2 to +1.3 ms.
- **Timing conditions.** Timings are single-threaded on a loaded 2-CPU host. They are performance evidence for equivalent work, not a certificate.

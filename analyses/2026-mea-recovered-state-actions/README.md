# MEA solved states: intensive residuals, consumed actions and calorics

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The retained
actions, sweep and calorimetry require the identified MEA application revisions and
wheels; these frozen application loaders are excluded from current declaration checks.

Since #308 `scripts/run.py` makes no outcome-dependent choice: the packet start
is the only cold start (the native #152 speciated seed replaces the Python
restart), anchors come from the declared temperature and loading, and both
continuation protocols run on every request. The retained `results/*.csv` and
`summary.json` are the #61 record of the earlier driver; the
[same-wheel comparison](results/issue-308-native-comparison.json) shows the two
drivers agree (see "Native conversion" below). The #331 rerun with the phase-pressure row is
in `results/issue-331-*.json` (see "Phase-pressure row" below).

Questions (#131, extended by #61): on the adopted MEA model,
1. do the production anchored path and declared-branch feed continuation reach the
   cold branch of all 79 solved-pressure requests, and does every accepted state meet
   species/moiety and charge balance, intensive (unscaled) row, interphase fugacity,
   reaction-affinity, root and support criteria, trace species included?
2. are the solved-state actions MEA consumes, including the CO2–H2O k_ij direction and
   the caloric actions of its heat and Cp stages, correct on recovered states?
3. does MEA's rebuilt calorimetry stage evaluate every interval on this wheel, and does it
   reproduce MEA's pinned-wheel (`b66c7b96`) values?

Inputs are read from an MEA checkout at `04a9328` (sweep) or `799b201` (actions, heat: the rebuilt heat stage; `shared_evaluation` and the inputs are unchanged) (record `868a5018…`, packet
`86f60041…`), mapped to problems by MEA's production
`shared_evaluation._problem_from_request`. The Engine is the installed non-editable
wheel `ad12b19f…` (Engine main `f5ab8102`); `results/summary.json` records both per mode.

The retired command annotations recorded 23 s for actions, 521 s for the
sweep and 50 s for cold heat evaluation; the state cache reduced repeat heat
evaluations to seconds. `scripts/density_shift.py` supplied the #331 falsifier
on MEA `04a9328`. These timings belong to the original wheel.

## Frozen protocol

- **States for actions.** `B050`: Böttinger 2008 state 050 (313.15 K, 101325 Pa), the #142
  k_ij state. `R360`: the same liquid at 360 K, started from the cold-solved 345 K liquid.
  VLE 0128, 0145, 0227: the packet start, which the Engine seeds from the liquid speciated
  alone at the packet pressure seed (#152); the #61 record used a Python two-stage start.
  Each state also has a liquid-only TP twin at its solved pressure, solved with MEA's
  neutral ideal-gas records (`verify_reference_calorics.record`), for the caloric rows.
- **Actions** (first order, native units): every liquid mole fraction and residual chemical
  potential over RT, liquid ln a_CO2, amount and volume, plus pressure and vapor CO2 fraction
  for VLE. Directions T (K), P (Pa, fixed-P states only), CO2 feed (mol) and CO2–H2O k_ij.
  Caloric rows: liquid total enthalpy H, amount and ln a_CO2 along T and CO2 feed.
- **Reference.** External central differences of full re-solves seeded from the solved
  state, steps h and h/2 (T 0.04/0.02 K, P 200/100 Pa, feed 2e-4/1e-4 mol, k_ij 1e-3/5e-4,
  the k_ij re-solves rebuilding the model). Richardson value R = (4D(h/2) − D(h))/3,
  resolution e = |D(h/2) − D(h)|/3. Pass: |a − R| ≤ 1e-5|R| + e; a reference with
  e > 1e-3|R| is unresolved. VLE 0227 sits at the 393.15 K upper bound of R2's source
  correlation, so its T reference uses the second-order backward difference.
- **Heat identity.** Absorption heat h_CO2^ig − ∂(nH)/∂n_CO2 from the enthalpy feed action
  equals RT² ∂ln a_CO2/∂T|P + RT (Gibbs–Helmholtz), relative ≤ 1e-6.
- **Sweep.** Reference branch: the MEA packet start. Anchored attempts start once from the
  nearest-loading neighbour at the same temperature and from the nearest cross-temperature
  neighbour, chosen from declared temperature and loading; an unsolved neighbour is recorded
  as `neighbour_unsolved`. Continuation runs from the same-temperature neighbour's state with
  step 0.25, minimum 1/64, 64 steps and 4 retries, and again with minimum 1/512, 512 steps and
  8 retries. The #61 record restarted a failed cold solve from a speciated liquid, kept only
  solved anchors and ran the finer protocol only after a stop.
- **Acceptance.** Engine success, MEA material and charge balance ≤ 1e-7, pressure inside
  the packet bounds; branch identity with the cold state |ΔP|/P, |ΔpCO2|/pCO2 ≤ 1e-8.
- **State criteria (#61), on every accepted state.**
  - Balance: n_out − n_feed must be a reaction extent (every conserved moiety), and
    Σ z_i n_i = 0, both ≤ 1e-10 relative to the feed.
  - Every unscaled chemical and pressure row ≤ 1e-10 (the Engine's default caller row
    criterion), not the amount-weighted rows; trace species (x < 1e-6) are reported
    separately.
  - Reaction affinity |Σ_i ν_ri r_i|/Σ|ν_ri| ≤ 1e-10 for R1–R5 (the moiety multipliers cancel).
  - Interphase neutral fugacity: each phase is re-evaluated independently at its solved
    (T, ρ, x); |Δ(μ_res/RT + ln x_i ρ)| ≤ 1e-9 for CO2, MEA, H2O, and each phase pressure
    as the Engine's pressure row, |P_phase − P|/(ρ_phase R T) ≤ 1e-10 (since #331; the #61
    record used |P_phase/P − 1| ≤ 1e-9).
  - Roots and support: ∂P/∂ρ > 0 in every phase; liquid packing fraction > 0.2 and incipient
    vapor < 0.01; every liquid species present; vapor support exactly CO2, MEA, H2O.

## Result

**Sweep** (`results/sweep.csv`). Every accepted state on every path meets every state
criterion; the largest values over all 316 accepted states are below.

| Path | Accepted, on cold branch | Iterations | Largest balance / charge | Largest unscaled row (trace) | Largest affinity | Largest Δ ln f |
|---|---|---|---|---|---|---|
| Cold packet start | 79 / 79 | 212 | 1.8e-16 / 1.5e-16 | 1.4e-13 (1.4e-13) | 6.4e-14 | 1.4e-13 |
| Same-temperature anchor | 79 / 79 | 196 | 1.8e-16 / 2.0e-16 | 2.0e-13 (1.7e-13) | 8.5e-14 | 1.7e-13 |
| Cross-temperature anchor | 79 / 79 | 211 | 1.7e-16 / 1.7e-16 | 1.7e-13 (1.7e-13) | 7.7e-14 | 1.2e-13 |
| Feed continuation, frozen budget | 76 / 79 | 2309 | 1.5e-16 / 2.6e-16 | 1.7e-13 (1.7e-13) | 7.1e-14 | 1.9e-13 |
| Feed continuation, finer steps | 3 / 3 | 496 | 9.8e-17 / 9.4e-17 | 7.1e-14 (5.7e-14) | 5.0e-14 | 4.2e-14 |

- The trace species are H3O+ (316 phase states), OH− (160) and molecular CO2 at low
  loading (104); the smallest liquid mole fraction is 6.2e-12. Their unscaled rows are at
  the same rounding level as the major species.
- Independent phase pressures agree with the solved pressure to ≤ 7.6e-10 relative. For
  the nearly incompressible liquid (Z ≈ 5e-5) this is the double-precision floor of the
  re-evaluation, not a tolerance: one double-precision step of liquid density moves the
  relative pressure by 8.6e-11, so 7.6e-10 is about nine such steps (#331). The
  largest incipient-vapor packing fraction is 0.0077.
- All 79 packet starts now converge without the speciated-liquid restart (212 iterations,
  against 1809 on the #131 wheel, after #154/#157).
- Every pressure equals the #131 wheel `2a8f0904` value to ≤ 2.2e-13 relative (312 paired
  states), so the adopted-model states are unchanged by the Engine changes since.
- The three frozen-budget continuation stops (0201→0202, 0209→0208, 0232→0233) are the
  #131 step-floor stops: the first step leaves the 1 % anchored density window
  (`DENSITY_ANCHOR_LOG_WINDOW`), cause `branch_lost`, 0 accepted steps. With the finer
  steps all three reach the target on the cold branch. None is unsupported or physically
  rejected. Every accepted state agrees with the cold state, so one root was found per
  request and no multiple-root selection rule is needed.

**Actions** (`results/actions.csv`). All 405 rows are Available and pass; no reference is
unresolved; the largest error/criterion is 0.81.

| Direction | Rows | ≤ 1e-7 relative | Largest relative error |
|---|---|---|---|
| T | 111 | 111 | 8.4e-8 |
| P (B050, R360) | 42 | 32 | 1.0e-5 |
| CO2 feed | 111 | 111 | 5.2e-8 |
| CO2–H2O k_ij | 111 | 78 | 4.3e-6 |
| Caloric, T (equilibrium Cp) | 15 | 15 | 9.9e-8 |
| Caloric, CO2 feed (partial molar enthalpy) | 15 | 15 | 7.0e-8 |

- The ten P rows above 1e-7 belong to the nearly incompressible fixed-P liquids; a 100 Pa
  re-solve moves them by 3e-10 to 8e-7 relative, so the reference, not the action, bounds
  the agreement.
- The B050 x(MEACOO−) k_ij action is −9.3287292e-6, the #142 value; the reference here is
  −9.3287258e-6 (3.6e-7). The k_ij rows above 1e-7 have references resolved only to
  1e-6–9e-5 of their value (k_ij re-solves at the row tolerance), and each passes.
- Heat identity 5/5, largest 3.2e-13. On the liquid-only TP twins the absorption heat is
  72.84 (B050, 313 K), 73.67 (R360), 82.29 (0128), 83.48 (0145) and 96.09 kJ/mol (0227,
  393 K); the closed-system equilibrium Cp, (n ∂h/∂T + h ∂n/∂T)/n, is 73.5–79.0 J/(mol K) of liquid.

**Heat** (`results/heat.csv`). MEA's `evaluate_direct_absorption_heat.score` on this wheel
evaluates all 113 calorimetry intervals. RMSE 11.50 kJ/mol CO2 (bias +0.3) for the 66
Kim–Svendsen 40/80 °C calibration intervals, 29.28 (bias −15.8) for the 20 Kim–Svendsen
120 °C holdout intervals and 10.61 for the 27 Kim 2014 rows. The 109 intervals MEA's pinned
wheel evaluates agree with it to ≤ 4.0e-9 kJ/mol. The other four start at loading 0.003 or
0.037 at 40 °C, where the pinned wheel's bubble solve is infeasible from every start; this
wheel solves them from the packet start in three iterations (the #154/#157 cold-start change; `heat.cold_packet_start` in `summary.json`) and scores them
81.4–81.6 kJ/mol against 82.1–84.1 measured.

## Native conversion (#308)

The earlier and converted drivers ran on one wheel (`4181cdf3…`, Engine main `32b0f48d`)
with MEA `04a9328` (sweep) and `799b201` (actions); criterion
|a − b| ≤ 1e-10 max(|a|, |b|) + floor, floors 1e-12 Pa and 1e-14
([comparison](results/issue-308-native-comparison.json)).

- The earlier driver took the packet start on all 79 requests; its restart never fired.
- Cold, anchored and frozen-continuation values are bit-identical, with the same neighbours
  and success flags. The three frozen-protocol stops reproduce `step_too_small`, cause
  `branch_lost`, 0 accepted steps in 5 attempts.
- The finer protocol now runs on all 79 requests (395 rows instead of 319). Its three
  earlier rows are bit-identical; the 76 added rows follow their frozen rows' step history
  and are bit-identical to them. The frozen continuations took 4 to 64 steps, so the
  finer controls change nothing unless the step floor, step budget or retry limit is reached.
- VLE 0128/0145/0227 from the packet start agree with the two-stage start to 1.3e-12
  relative; all 405 actions agree to 2.0e-12, are Available and pass; heat identity 5/5.
- On this wheel the independent phase pressure of the 0188 frozen continuation (and its
  finer twin) differs from the solved pressure by 1.31e-9 relative, above the #61 1e-9
  criterion, in both drivers. The #61 record's largest was 7.6e-10; the change predates
  #308 and is diagnosed in #331 (below).

## Phase-pressure row (#331)

The 0188 miss was rounding, not solve error. The liquid has Z = 5.1e-5 and
ρ(∂P/∂ρ)/P = 5.9e5, so one double-precision step of its density is 8.6e-11 in relative
pressure, and the independent re-evaluation scatters by 1.3e-10. #325 moved the 0188
continuation density by about ten such steps. The check therefore uses the quantity the
Engine's pressure row solves, (P_phase − P)/(ρ_phase R T), with the Engine's row tolerance
1e-10. For the vapor (Z ≈ 1) it is the relative check, tightened from 1e-9 to 1e-10. For
this liquid it is looser than #61's: it allows a relative pressure difference of about 2e-6
(0.013 Pa at 6.7 kPa), a density error near 3.3e-12, where #61 allowed 1e-9 (a density
error near 1.7e-15, about 12 rounding steps). #61's liquid bound sat inside the
re-evaluation's own scatter, and the solver converges this row only to 1e-10, so an
independent recomputation is held to the solver's standard.

Rerun on wheel `c921e015…` (Engine main `645c58a4`, after #322), MEA `04a9328` (sweep) and
`799b201` (actions), retained with the comparison against the #308 outputs (wheel
`4181cdf3…`, main `32b0f48d`) in [`issue-331-rerun.json`](results/issue-331-rerun.json), a
row-by-row bookkeeping of the two drivers' `sweep.csv` and `actions.csv` outputs:

- **Sweep.** 392 of 395 rows are accepted, on the cold branch, and meet every criterion; the
  three rejected rows are the frozen-budget `branch_lost` stops, which the finer protocol
  completes. The largest phase-pressure row is 5.7e-14 for a liquid (0125, cross-temperature
  anchor) and 8.8e-13 for a vapor (0130, same-temperature anchor). The 0188 continuation
  rows read 3.9e-14 (liquid) and 2.7e-16 (vapor), against 1.31e-9 relative before.
- **Against #308.** Only the two 0188 continuation rows change outcome (criteria now met). No
  other criterion changes pass or fail. Solved values move by rounding: largest 1.7e-13
  relative in pressure, 2.5e-13 in the smallest liquid mole fraction and 1.9e-13 in vapor
  packing fraction. Work is identical except the 0210 → 0209 continuations (262 → 261
  evaluations). Iteration totals per path are unchanged (212, 196, 211, 2309, 2797).
- **Actions** (the five states use the same criteria). All 405 rows are Available and pass,
  the largest error/criterion is 0.44, heat identity 5/5. Values and actions move from #308
  by at most 8.5e-13 and 1.3e-12 relative.
- **Falsifier** ([`issue-331-density-shift.json`](results/issue-331-density-shift.json),
  `scripts/density_shift.py`). Shifting the 0188 continuation liquid density by ±1e-11
  relative gives rows of 3.01e-10 and 3.02e-10 and fails. A shift of ±1e-12 gives 3.01e-11
  and 3.02e-11 and passes. The predicted shift contribution is 3.015e-10 and 3.015e-11; the
  unshifted row is 3.9e-14.

## Limits

- The evidence is bounded to this model, packet, wheel and the listed starts. It is
  numerical verification, not physical agreement with measurements.
- Actions are first order, which is what MEA consumes. Finite differences are external
  references only; the Engine has no finite-difference runtime fallback. Reaction-coefficient
  actions do not exist (#61 limit).
- Caloric rows are on liquid-only TP states; pressure-moving caloric actions on the bubble
  path are #147's scope and are not consumed by MEA.
- Continuation along temperature or with free temperature is not exercised beyond the
  Engine test's referenced ionic T/P path.

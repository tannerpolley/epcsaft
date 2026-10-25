# Final application acceptance on one Engine wheel (#64)

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The replay
requires the original pinned MEA and Lithium application checkouts and non-editable
wheel; its frozen application loaders are not current-main qualification commands.

Question: do the accepted MEA-Thermodynamics and Lithium_Extraction results hold on
one identified, non-editable Engine wheel? Checks invalidated by Engine changes since
each application's pin are replayed on that wheel. Numerical correctness, physical
comparison error, surrogate accuracy and process results are reported separately.

## Accepted basis (owner decision, 2026-09-29)

- **MEA:** adopted record `868a5018…`, packet `86f60041…` (#61), from MEA checkout
  `799b201`. The record omits both Born coefficients; MEA's loader
  (`shared_evaluation.py:parameter_mapping`) fills `c_shell = c_dielectric = 1.0`, so the
  SSM+DS (1,1) Born form is evaluated (water has f_solv = 1.5; MEA PR #125). The Engine's
  own loader rejects omitted coefficients.
- **Lithium:** the unbuffered one-stage calibration, ln K_exchange 48.056433 and
  ln K_Na −0.445812 (Lithium PRs #147 `38af594` and #161 `78f8140`), from Lithium main
  `53f1714`.
- A new MEA record (model D) or a new Lithium basis is replayed under its own issue.

## Final wheel

`28181e72e429c6e87fc6361082af1a7b21c8747e76bda65a1c30abb4a97402f2`, Engine `1a905ff8`
(the environment wheel for that source and toolchain fingerprint, restored from the shared
cache; the same SHA-256 as the qualified wheel kept for the #129 fit audit). The full
replay first ran on `b5566590…` (Engine `80cf018b`) and was repeated on this wheel; the two
runs differ by at most 9.1e-13 relative (MEA heat RMSE) and not at all for Lithium. The
retained evidence it replaces as the common identity:

| Evidence | Wheel | Engine |
|---|---|---|
| MEA Engine record (#61: sweep, actions, heat) | `ad12b19f…` | `f5ab8102` |
| MEA application scoring (pCO2, speciation) | `b66c7b96…` (MEA pin) | `443a9da4` |
| Lithium MS-2 chain (#147/#161) | `18c7812f…` (Lithium pin) | `1303c119` |

## Engine changes since the pins

| Merge | Runtime change | Effect on consumed results |
|---|---|---|
| #141–#167 (`443a9da4` → `f5ab8102`) | Gibbs minimizer, Helmholtz model and tape, association, phase boundary, infinite-dilution reference, regression | Affect MEA's scoring on its pin (`b66c7b96`); replayed |
| #170 | Mixed T second actions and ∂h/∂P on reacting problems; changes the infinite-dilution reference and request compilation | Replay shows consumed values unchanged |
| #174 | Inventory-capped composition box, moiety echelon balance basis | Changes optimizer paths: MEA replay required |
| #177 | Compensated action refinement under an ASL host | Rounding level (≤ 1e-13) on actions |
| #183 | Build only (NLopt/Boost removed; neither was included by any source) | None |
| #184 | No state after fatal Ipopt returns; `FEASIBLE_POINT_FOUND` polished; strict mapping/Reference admission | Admission only; replay confirms no consumed solve changed class |
| #197 | Fitting-state reuse and bound-set refresh in the Engine fitter; diagnostic probes (active only in a diagnostic scope) in the infinite-dilution reference | No replayed check calls the Engine fitter; MEA moves ≤ 9.1e-13 relative from the `b5566590` run, Lithium not at all |
| #168, #169, #171, #178, #185, #189, #194, #195, #196 | None (docs, analyses, tooling) | None |

MEA is affected by all of the runtime changes; Lithium by #183 and #184 only. Every consumed
check was therefore replayed.

## Method

`replay.sh` runs one step at a time, single-threaded, once the load average is below 5
(about 26 min total). The Lithium checkout's `.venv` holds the final wheel. MEA steps use
the #61 driver `analyses/2026-mea-recovered-state-actions/scripts/run.py` (copied to
scratch so the #61 record is kept) and `scripts/mea_scores.py`, which runs MEA's own
`generate_figure_data` and `score_observation_coverage` with the wheel pin overridden.
Lithium steps use Lithium's own drivers with `--engine-wheel/--engine-sha256/--engine-commit`.
Both checkouts are disposable detached worktrees; their in-place outputs are compared
with the committed files.

**Criteria.** Each check applies its stored criteria unchanged. In addition, every solved
value (state, action, prediction) must agree with its retained value to 1e-8 relative,
the #61 branch-identity criterion. Convergence diagnostics (residuals, KKT, timings,
iterations) are gated by the stored criteria, not compared.

## Result

All 11 replayed checks pass (`results/replay.csv`):

| Application | Check | Values compared (bit-identical) | Largest relative change | Stored criterion | Result |
|---|---|---|---|---|---|
| MEA | VLE sweep, 79 requests × cold / two anchors / continuation | 948 (35) | 2.3e-13 | 316/319 path states accepted by the pre-#308 driver, the same set as #61; all state criteria met, cold branch | pass |
| MEA | First-order actions vs Richardson re-solves (T, P, CO2 feed, k_ij, calorics) | 810 (61) | 1.6e-11 | 405/405 Available and within \|a − R\| ≤ 1e-5\|R\| + e | pass |
| MEA | Calorimetry intervals | 113 (20) | 1.2e-12 | 113/113 evaluated | pass |
| MEA | pCO2 (161 states) and packet speciation predictions | 292 (1) | 2.1e-11 | 0 failed states | pass |
| MEA | Observation-coverage speciation | 332 (1) | 8.9e-13 | 274/274 positive rows scored, 58 nondetections | pass |
| Lithium | One-stage refit, five starts | 289 (289) | 0 | 5/5 starts converge; strict calibration state accepted | pass |
| Lithium | Zhang 2017 holdouts (Fig. 5, Table 4, Table 2 cascade) | 160 (160) | 0 | 11 states strict | pass |
| Lithium | Direct-state volume scaling | not compared | — | Scales 1–1e5 evaluated; invariance ≤ 8.7e-14 (limit 1e-10) | pass |
| Lithium | Live IDAES campaign | 824 (824) | 0 | 13/13 units accepted, cohort holds | pass |
| Lithium | Consumed derivatives | 740 (740) | 0 | 84/84 ladders pass; 50-digit A1 action 2.2e-16 | pass |
| Lithium | v3 surrogate states | 29,980 (29,980) | 0 | 5,996/5,996 strict; refitted networks identical | pass |

- **MEA** values move at rounding level only, as expected from #174's new optimizer path.
  The three frozen-budget continuation stops (0201→0202, 0209→0208, 0232→0233,
  `branch_lost`) are unchanged and reach the target with finer steps. Cold starts still
  take 212 iterations in total. The four 40 °C calorimetry endpoints that MEA's pinned wheel
  cannot solve evaluate on the final wheel, as on #61's.
- **#308** converted the MEA driver after this replay: every request now declares both
  continuation protocols (395 rows), and `compare.py` requires each retained row to keep
  its success flag and every added row to be accepted. This table was not rerun. On wheel
  `4181cdf3…` the 0188 continuation missed the #61 phase-pressure check (1.31e-9 relative
  against 1e-9), a double-precision floor for this liquid, not a solve error. Since #331 the
  driver checks the Engine's pressure row |P_phase − P|/(ρ_phase R T) ≤ 1e-10. All 392
  accepted rows meet every criterion on the current main wheel (`c921e015…`, main
  `645c58a4`). The same-wheel comparison is in
  [`issue-308-native-comparison.json`](../2026-mea-recovered-state-actions/results/issue-308-native-comparison.json),
  the #331 rerun in [`issue-331-rerun.json`](../2026-mea-recovered-state-actions/results/issue-331-rerun.json).
- **Lithium** is bit-identical in every solved value, as expected: #183 and #184 change no
  computation after a non-error Ipopt stop, and #197 changes no calculation they call. The regenerated `loaded_state_surrogate.json` differs from
  the retained file only in its wheel, commit and native-core hashes and fit timings.
- **Volume scaling** was rerun, not replayed: the retained run predates the #147 fit. The
  3e5 scale is again rejected (balance 1.9e-9 mol against 1e-9), the recorded #63 limit.

## Physical comparison and process results (`results/quantities.csv`)

They are model/data comparisons, not numerical checks, and carry no acceptance band
beyond the application's own. The pCO2, speciation and Lithium values are unchanged by
the wheel. The heat RMSE differs from MEA's pinned-wheel record only because four
40 °C intervals, which that wheel cannot solve, now evaluate. The 109 common intervals
agree with the pin to 4.9e-11 relative.

- **MEA pCO2** (calibration, 143 rows): RMS ln 0.544, AARD 53 %. Xu 2011 holdout (18 rows):
  RMS ln 0.440, AARD 49 %. The model's pCO2 error is roughly a factor of 1.5–1.7; it is
  a model/reference issue under MEA #107/#121, not solver error.
- **MEA speciation** RMS ln: Böttinger packet 0.23, Böttinger outside packet 0.41,
  Matin 0.42, Jakobsen 0.99.
- **MEA heat of absorption** RMSE (kJ/mol CO2): 11.50 calibration (66 intervals; 11.77 on
  the pin's 63), 29.28 120 °C holdout (20), 10.61 Kim 2014 (27; 10.80 on the pin's 26).
- **Lithium holdouts:** Fig. 5 max |log10(D_pred/D_obs)| 0.256; Table 4 R = 70 % D_Li 17.7
  against 7.77; Table 2 raffinate 0.0072 against 0.045 g/L.
- **Lithium live designs** (unbuffered fresh-organic MS-2 basis): 1 and 2 stages reach at
  most 75.4 % and 94.4 % Li; 3 stages 95.00 % at $420.2/kg Li (O/A 1.080); 4 stages
  98.82 % at $417.6/kg Li (O/A 1.107).

## Shared capabilities and fitting evidence

`bash tools/check-all.sh` on this branch runs the persistent numerical checks of every
#81 row on the final wheel. It covers the Python suite and the native `contributions`,
`properties`, `diagnostics`, `specification`, `surface_example` and `regression`
(R1, D86) checks; it passed on this branch (107 Python checks). Three native programs it
omits cover #81 rows and were run once on the same build, all passing: `equilibrium`
(E5 solved-state actions, including the #130 forward-bound and #177 ASL-host falsifiers;
22 s), `boundary` (S3–S15, E4; 13 s) and `regression_fit` (R2, R3; 22 s). `check-affected`
now runs them when a lower layer changes (#190). `check-all` also omits the native
`properties_floor` and `association-comparison` checks; they are not run here.

| #81 row | Owner | Evidence | Final-wheel check |
|---|---|---|---|
| Reference and thermal chains | #84, #138, #140, #147 | `2026-reactive-electrolyte-driver`, `2026-reacting-mixed-actions`; MEA caloric actions | MEA actions/heat replay; check-all neutral-reference and calorics tests; native `equilibrium` |
| Active EOS parameter families | #85 | Reviewed no-code inventory (seven families) | check-all `test_fit.py`, native `regression` |
| Density prediction actions | #86 | Native D86 rows | check-all native `regression` D86 |
| Fitting: observations, fits, failures, real data | #45–#49, #152 | `2026-pure-saturation-regression`, `2026-mea-reactive-fit-timing`; Lithium ln K refit | check-all `test_fit.py`, native `regression_fit`; Lithium refit replay |
| MEA candidate and qualification | #82, #61 | `2026-mea-recovered-state-actions` | Replayed above |
| Lithium reconciliation and qualification | #83, #63 | Lithium PRs #147/#161 | Replayed above |
| Trace species and consumed small actions | #87, #130, #172 | `2026-trace-action-diagnosis`, `2026-stability-noise-band` | check-all `test_trace_reliability.py`; native `equilibrium` E5-8 |
| Density closure and branch selection | #89 | `2026-density-closure-diagnosis` | check-all native `properties` P2/M8 |
| Initialization and continuation | #90, #131, #157 | `2026-initialization-diagnosis`, `2026-reactive-liquid-cold-start`; MEA sweep | MEA sweep replay; check-all continuation tests; native `boundary` |

The application fits (MEA record `868a5018`, Lithium ln K) are fixed inputs here. The
Lithium ln K refit is reproduced bit-identically. MEA's record is not refitted:
parameter regression is unnecessary for a fixed-input replay.

## Documentation

Live views now separate historical evidence from current capability:
- `analyses/reproduction/2026-khudaida/README.md` labels its HELD2 discovery campaign as
  historical on its own wheels. `run-slice.sh khudaida` runs the declared two-liquid
  settler.
- `analyses/reproduction/2018-wangler/README.md` labels its calculations historical
  (#186): they use retired interfaces and do not run on the current wheel.
- `docs/scientific/CONTEXT.md` lists local stability checks and marks phase discovery as
  deferred.

No live view claims that a B1–B4 capability is missing.

## Limits

- **Numerical scope.** These are local, declared-phase results for the listed states, starts
  and inputs. There is no phase-discovery or global-stability claim. Agreement with the
  retained values shows that the final wheel reproduces them; it does not validate them
  physically.
- **MEA physical accuracy.** pCO2 and speciation errors are calibration residuals, and the
  Jakobsen and high-loading Böttinger rows are previously accessed data. The 120 °C heat
  holdout is under-predicted by 15.8 kJ/mol on average. Liquid Cp and loaded density are
  open in MEA #124/#123. In the record's SSM+DS (1,1) Born form the Born shell diameter
  depends on f_mix, which averages over every component including the ions, so the Born
  term adds to ln a_w and MEA's check E3 fails (MEA PR #125).
- **Lithium physical caveat.** Lithium #158 (PR #160) showed that fresh organic on MS-2
  brine self-acidifies it to pH 2.6–3.1 and extracts at most 4 % per stage, not the 41–75 %
  of this basis. The fresh-organic MS-2 chain is numerically qualified but physically
  superseded. A saponified alternative flowsheet is separate work. Physical MS-2 transfer
  accuracy is unvalidated, and the model does not represent brine alkalinity. MS-2's Ca
  (36.9 g/L), Mg, K and Sr are not modelled, and β-diketones extract Mg and Ca before Li
  (Lithium PR #160), so the live designs' Li recovery and cost assume divalent removal
  upstream.
- **Surrogate.** Accuracy is the stored p95/max limit (0.03/0.08) against Engine states,
  not against measurements. Lithium's surrogate consumers refuse a different native core by
  design. Moving them to this wheel means updating the hashes in `loaded_state_surrogate.json`, which is a Lithium
  decision; the states and networks are shown identical here.
- **Pins.** On 2026-09-29, MEA main pins `f66d972c` (Engine `1303c119`); its own
  qualification on `28181e72` found its 84-state, 160-target initial point bit-identical.
  Lithium main pins `b5566590` by its own decision. Re-pinning is an application decision.
  This record shows both applications' #64 results reproduce on `b5566590` and
  `28181e72`.
- **Known numerical limits kept.** The MEA frozen-budget continuation stops and the Lithium
  3e5 volume scale, both above.

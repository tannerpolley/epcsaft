# Shared-state FeOs / Clapeyron cross-check

**The shared solved values agree within the cross-library value targets for C1–C4.** FeOs acetone saturation misses the 1e-8 pressure/density target (about 1.2e-7 relative). The large motivating NaCl gap is reproduced and decomposed: screening diameter explains part, zeroing the water–chloride interaction reproduces another part, and the finite reference explains a much smaller remainder. It is not all a diameter effect.

Strict common-state agreement is **not universal**. Dispersion coefficient precision, dipole-unit conversion, the Clapeyron DH rational approximation and floating-point cancellation are retained separately from raw differences. After those computed corrections, most Helmholtz/Z differences are below their targets, but dense-liquid log-fugacity misses remain; C4 also has a Clapeyron derivative remainder. The corrected Engine meets the unchanged 1e-10 relative-pressure target for all 12 primary and 14 ladder C3 rows (worst 6.96e-11/8.46e-11); the original 3.88e-10 miss remains under its old wheel. Clapeyron C4 pressure/chemical closure still reaches about 2.88e-9. Passing value comparisons do not make those closure checks pass. This is numerical verification on these inputs, not experimental validation or a library ranking.

[Research notebook](notebook.qmd) contains the conditions, equations, diagnostics, interpretation and limits. [Issue #222](https://github.com/tannerpolley/ePC-SAFT/issues/222) owns the accepted design and amendments.

## Agreement

Maxima over selected solved rows; relative differences use Engine as denominator. Comparisons exclude deliberately unmatched NaCl convention probes.

| Case | FeOs versus Engine | Clapeyron versus Engine | Meaning |
|---|---|---|---|
| C1 propane | P 5.68e-10; density 5.74e-10 relative | P/density at most 4.83e-15 relative | Solved target 1e-8 met |
| C2 ethanol/butane | P 1.65e-10; density 2.26e-10 relative; vapor x 1.94e-11 absolute | P/density at most 3.42e-13 relative | Solved targets met |
| C3 NaCl, matched 0.88σ | density 7.99e-11 relative; Δlnγ 1.18e-6; ΔΦ 2.48e-8 | density 5.05e-11; Δlnγ 1.18e-6; ΔΦ 2.27e-8 | Electrolyte target 1e-5 met; FeOs zero-ion versus Engine finite reference distinguished |
| C3 NaCl, matched σ | Not an available FeOs screening selection | density 1.03e-10 relative; Δlnγ 1.18e-6; ΔΦ 4.62e-8 | Unchanged Held screening convention |
| C4 water/pentanol | density 7.77e-11 relative; x 7.40e-11; water-rich mole fraction 9.40e-11 absolute | density at most 5.03e-12 relative | Nontrivial split; common water feed remains 0.5 |
| C5 acetone | P/vapor density 1.21e-7/1.22e-7 relative: **miss** | Solved target met | Raw and matched dipole conversion rows remain separate |
| C1/C2 pressure gradients | maximum relative differences 1.92e-10 / 2.24e-10 | Not selected by adopted D6 | Every selected coordinate meets mixed 1e-6 relative + 1e-8 absolute rule |

The exact rows, signs, units and targets are in [values.csv](results/values.csv). [source-adjusted.csv](results/source-adjusted.csv) does not relabel raw misses as passes. Negative-pressure association/acetone falsifiers retain valid energies and derivatives; real log fugacity is outside its domain there.

## Original optimized timing campaign

Milliseconds, median of three fresh-process medians; each process has five warmups and 25 retained samples. Complete public case-output costs, **not equal algorithms or derivative work**. Engine's general equilibrium driver and T3 duplicated prediction are included. FeOs C3 includes Julia contraction and transport; Clapeyron C3 includes both native helpers and their repeated reference solves.

| Block | Engine† | FeOs | Clapeyron |
|---|---:|---:|---:|
| C1 saturation, 300 K | 8.761 | 0.04949 | 0.05638 |
| C2 bubble, xEtOH=0.5 | 21.990 | 1.148 | 0.3292 |
| C4 LLE, 298.15 K | 22.282 | 43.109 | 0.8908 |
| C4 LLE, 323.15 K | 22.776 | 11.183 | 1.027 |
| C3 0.88σ, reused Engine reference / complete comparator helpers | 2.226 | 0.3511 | 0.07103 |
| C3 0.88σ, fresh Engine reference | 9.157 | — | — |
| C1 value + pressure gradient | 18.583 | 0.01043 | — |
| C2 value + pressure gradient | 55.888 | 1.123 | — |
| C1 fixed-state aʳ/Z/lnφ | 0.02657 | 0.06554 | 0.005811 |
| C2 fixed-state aʳ/Z/lnφ | 0.07899 | 0.07391 | 0.003857 |
| C3 0.88σ fixed-state aʳ/Z/lnφ | 0.1423 | 0.07989 | 0.01046 |

AMD Ryzen 5 5500, CPU 0 affinity, qualification-time performance governor snapshot, `nice -n 10`, all thread pools limited to one. The qualification-time host snapshot records performance; governor changes during the timing processes were not logged. Normal optimized Julia (`compile_enabled=1`, optimization level 2), **not `--compile=min`**. Measured cold work summed across imports/model construction/first calls was 604–696 s per Julia process, versus 0.230–0.334 s Engine and 0.603–0.709 s FeOs. These are measured block sums, not whole process startup times. Total measured warm work including warmups was 22.955 s. Scatter and every sample are retained, including the 124.8 ms Clapeyron warm outlier; no trimming.

† The original Engine timing processes retained no numerical return records. Their costs remain historical measurements with unqualified returned values. The separate a0 replay below does not qualify those absent returns.

See [timing-statistics.csv](results/timing-statistics.csv), [timings.csv](results/timings.csv) and [budgets.csv](results/budgets.csv). The 642 retained comparator fields reproduce their implementation's qualification baseline (306 FeOs, 336 Clapeyron). The original Engine timing files contain no returned numerical records, so this check does not yet qualify Engine timed outputs. Neither comparator file retains an independent numerical snapshot of every warm sample.

## Separate Engine replay

The 2026-10-03 a0 replay retains 14 selected request/selection records and 350 warm samples per process, over three serial processes with five warmups/25 samples per selection. All 183 Engine advertised fields reproduce the old qualification baseline exactly. [timed-output-checks.csv](results/timed-output-checks.csv) contains those fields plus 150 FeOs/147 Clapeyron advertised first-call fields from the original campaign, also unchanged; the original 642 qualification checks remain in [timed-output-checks-original.csv](results/timed-output-checks-original.csv). Only the last measured Engine return is retained, so individual warm-sample equality is not established.

The [notebook](notebook.qmd#separate-engine-replay-returned-values-and-observed-costs) gives the separate cost table and process-median ranges. [Replay statistics](results/engine-replay-statistics.csv) preserve per-process p10/p90, [samples](results/engine-timing-replay.csv) preserve all measurements, and [budgets](results/engine-replay-budgets.csv) record 10.2104 s of warm work and a 57.4224 ms maximum sample. Actual niceness was 19 after an inherited positive niceness plus the declared +10 adjustment; CPU affinity was [0], governor endpoint observations were performance, and one-minute load was 1.395–1.524. No quiet-host or paired comparator speed claim follows. The historical strict closure/source misses remain unresolved; this selected-value replay does not reassess or remedy them.

[Replay identity](results/engine-replay-summary.json) and [raw returns/logs](results/engine-replay-raw/Engine-timing-replay-1.json.gz) separate the a0 wheel from case-driver checkout `b5a13cea` and historical qualification. The installed a0 bytes were verified against all 14 wheel RECORD hashes before replay; the corrected #231 wheel was restored and verified after it.

## Native density correction

The owner-approved revision 1 uses the bracket endpoint with the smaller pressure error for native Newton and publication, and returns the sole stable anchored local-scan root with composed counts. All 26 cold C3 rows retain two stable roots, positive slope and liquid branch. Iterations decrease 8 → 7 and evaluations 53–54 → 52–53. All five pinned/near-window anchor cases meet pressure closure; evaluations remain 11/14/15 while reported iterations include the existing local refinement. The existing associating-binary fixture at 353.15 K and 1 bar fails the old closure target (2.71894e-10) and passes after correction (8.83010e-12), earning one value check in the existing P2 group. The ethane probe already passed and remains retained.

[Corrected density rows](results/density-closure-corrected.csv), [reference rows](results/density-closure-reference.csv), [revised benchmark values](results/density-closure-values.csv) and [notebook](notebook.qmd#native-density-publication-and-corrected-engine-replay) preserve wheel identities, both endpoints and native traces. All 99 Engine requests returned without error. 790 original Engine value rows change across 41 requests; revised differences use the original comparator returns, with no outcome flips. Frozen common-state values and C5 do not change, so comparator execution was not repeated.

The first corrected Engine timing replay preserves five warmups/25 samples/three processes: 1050 warm/75 cold rows and 183 exact final-return matches to the new qualification baseline. Original comparator first-call fields add 297 unchanged matches. The notebook reports all 14 cost deltas and limits. Khudaida's native eight-point median block decreases from 252.338 to 226.117 ms (10.39%); corrected 25-sample totals are 5.7844/5.7882/5.6636 s versus reference 6.6212/6.1099/6.5140 s. Campaigns span reboot and do not support a paired causal estimate or new hard ceiling.

Native source `55e3b3a4` produced wheel SHA-256 `6039c149e7df951236fea278995ff9f29078ea19c58f354e92679c6da546cb57`; reference density/Khudaida rows use `e28d258df06eac39711c59a79f9ce2b69fc8ac7565637bfba0dbef21863d8dc7`. The initial full gate has seven mismatches to old retained native values, and the original native recorder refuses E5-3 reference uncertainty. The disposable larger-node probe passes every unchanged criterion; the owner approved these nodes for E5-3 only. All 268 unchanged criteria and allowances pass and are retained in [the reference-resolution rows](results/density-closure-E5-3-reference-criteria.csv); the refresh comprises 252 density-explained changes and one separately proven four-increment stale pressure prediction from the earlier accepted-warm publication correction (`33359588`). Both causes and old/new rows are retained in the notebook. The canonical Khudaida slice passed in 7.15 s. [Summary](results/density-closure-summary.json) retains failures, identities and raw-output hashes. The final merged-wheel results follow below; #222 closed when PR #273 merged (2026-10-06).


Final qualification on tree `81cd8582`, wheel `835d9bc27b243acb6d2288f4911828be973f07dc0d20766f397538065a9368b0`, reproduces all first corrected-wheel numerical values exactly. All 26 C3 rows, five anchors, both reference/activity traces, the P2 group, twelve recorder selectors and 480 timed-field comparisons pass. `check-all` passes (116 Python tests, one deselected, six native programs; 30.64 s), and the Khudaida slice passes (6.14 s). The native refresh changes 252 density-explained rows plus the one separately proven stale accepted-warm prediction; 225 remain unchanged.

The earlier final Khudaida campaign is retained as unpaired host timing: its median 252.337534 → 260.564580 ms (+3.26%) and aggregate +6.14% do not isolate density cost. The controlled comparison below owns the current cost interpretation.

## Paired ionic LLE and MEA cost

**No wall-time slowdown is resolved beyond six-pair spread.** Cold ionic point 7 has a median paired increase of 2.25% (range −21.50% to +21.31%; median absolute deviation 3.06 percentage points). It reproducibly needs four additional polish trial evaluations near the raw-row rounding floor, rather than additional density pressure evaluations. This does not prove zero overhead. MEA cold counts remain 43/57/44; its median paired change is −10.68%, with a −27.62% to +3.98% range, so no portable speedup is claimed.

The two wheels use merged tree `86b1ef14` including main `8f5f5685`, with identical inputs except the density owner. A restores pre-`55e3b3a4` density publication and its two policy-text lines; B keeps the pressure-closest endpoint and local-scan anchored root. Their SHA-256 values are `da6897394f85e0fd2a487899ad2296e209f2ee18c591cdfcc5aa7f68cf3ff382` and `f10774c1ff665861e80e3449a9ad3bb5bb5557d002546e16596db108663b811d`. The default #238 cases are unchanged: Khudaida ionic point 7 at 293.15 K/100000 Pa, and nine-species MEA reacting flash at 318.15 K/5918.5286 Pa. Each short lease installs A then B, with CPU 0, actual niceness 19, one thread per pool, five warmups and 25 samples. Cold and accepted-start costs are separate.

| Workload | A median / ms | B median / ms | Median paired change / % | Paired range / % | Paired MAD / percentage points | Iterations/evaluations/Hessians, A → B |
|---|---:|---:|---:|---|---:|---|
| Khudaida eight-point LLE block | 271.735197 | 276.684836 | -4.31 | -12.82 to +6.27 | 5.02 | 66/118/74 → 66/119/74 |
| Ionic point 7, cold | 55.688299 | 55.393255 | +2.25 | -21.50 to +21.31 | 3.06 | 16/49/17 → 16/53/17 |
| Ionic point 7, accepted start | 30.586151 | 30.085726 | -9.67 | -25.32 to +9.16 | 11.77 | 7/14/8 → 7/14/8 |
| MEA flash, cold | 114.361308 | 96.865381 | -10.68 | -27.62 to +3.98 | 6.34 | 43/57/44 → 43/57/44 |
| MEA flash, accepted start | 65.724829 | 66.113982 | -4.19 | -21.00 to +36.12 | 14.41 | 18/23/19 → 18/23/19 |

The wheel columns are medians of six process medians; the difference column is the median of within-pair percentage differences, which can differ in sign from the ratio of the two marginal medians. [Paired rows](results/density-paired-pairs.csv), [samples](results/density-paired-samples.csv), [host observations](results/density-paired-host-observations.csv), [polish trials](results/density-paired-polish-trials.csv) and [9360 output checks](results/density-paired-output-checks.csv) retain the result and spread. All states preserve two phases and native row criteria; maximum density/composition/amount changes are 1.65e-14 relative, 1.38e-14 absolute and 2.32e-14 relative. Load varies 4.03–11.75; this is numerical verification and host timing, not physical validation or a sub-percent cost bound.

Wheel B passes the full gate once (119 Python tests, one deselected and six fast native programs). [Provenance](results/density-paired-provenance.json) records the wheel identities of each timing run. [The notebook](notebook.qmd#paired-ionic-lle-and-mea-reacting-flash-cost) gives the cause and limits.

## Reproduction

Use the sole canonical non-editable Engine installation for new work. The retained original calculation uses the historical wheel described below; its collector snapshot is local commit `c1c29e94`. Its measured wheel SHA-256 is `da56e55eab2dac8b75c2497da07e4895d838d82d2d78b305e4a30b81bfef3950`; numerical Engine owners are from merged main `efcb50fc`. The original benchmark changed no Engine runtime; the separately authorized density revision below changes native endpoint publication with unchanged numerical criteria.

Comparator environments are task-local, not Engine dependencies:

- FeOs **0.10.1**, with [requirements-feos.lock](inputs/requirements-feos.lock) and [native controls/source locators](inputs/feos-controls.md).
- Julia **1.13.1**, Clapeyron **0.6.29**, with [Project.toml](inputs/Project.toml), [Julia dependency versions](inputs/Manifest.toml), and [archive hash/source controls](inputs/julia-source-controls.md).
- Install FeOs under `${XDG_CACHE_HOME:-~/.cache}/epcsaft-comparators/222/feos/`; extract the verified Julia archive as `julia-1.13.1/` and instantiate the locked project with `JULIA_DEPOT_PATH` pointing to this cache's `depot/`.
- `systemd-run --user --scope` supplies Julia's 6 GB resident-memory cap. Do not use a virtual-memory `ulimit`.

The historical a0 replay used wheel SHA-256 `398bcf8b90185b68cc796b5db258a9bfa3e919b9fb74bbb29171e003a03b3cc7`. The density correction uses the wheel and reference identities above. For current reproduction, install this checkout's non-editable wheel first (the first command below): the collector records its `build/environment-wheel/preparation.json`. Rerunning `bash .codex/environments/setup.sh` restores the editable install.

Within that lease and heavy scope, the existing Engine collector declares the full 99-request grid (including the 26 C3 states); a repetition argument selects the unchanged timing blocks. Write new output names to preserve previous campaigns.

```bash
bash tools/build-engine.sh --install-wheel
timeout 180 nice -n 10 taskset -c 0 .venv/bin/python \
  analyses/2026-competitor-cross-check/scripts/engine.py \
  analyses/2026-competitor-cross-check/inputs/cases.json build/Engine-current.json
for rep in 1 2 3; do
  timeout 180 nice -n 10 taskset -c 0 .venv/bin/python \
    analyses/2026-competitor-cross-check/scripts/engine.py \
    analyses/2026-competitor-cross-check/inputs/cases.json "build/Engine-current-timing-$rep.json" "$rep"
done
```

The density qualification additionally captures native `state_at_pressure` traces on the selected liquid requests and on the pinned density multiplied by 1/1.001/0.999/1.009/0.991. It reads the final `density_root_bracket` and native `density_diagnostics` endpoint pressure/slope values to compare endpoint selection, residual, radius, branch and count; it performs no Python density refinement. Native reference/activity spans are retained at 1 mol/kg under both screening selections. The temporary qualification declarations and retention readout remain in `build/density-closure-qualification/`; retained outputs and input hashes provide the durable numerical evidence. Notebook rendering uses retained rows with execution disabled.

Each replay subprocess has a 180 s limit, one allowed CPU, one thread per pool, five warmups and 25 samples. All three processes run serially. Replay outputs, logs, samples, statistics and budgets use separate names; the original campaign remains unchanged. The last measured return from each Engine block/selection is retained after its clock stops, with actual interpreter/preparation, case-driver hash, UTC times, affinity, load and governor observations. Qualification compares each campaign’s returned values to its identified retained baseline under unchanged targets. It cannot qualify missing historical returns or every individual warm sample.

The existing `qualify`, `diagnose` and `timing` modes are the complete-study producers; they rewrite primary outputs and require separately authorized reproduction of that study. The retained original study's wheel, libraries, locks, inputs and reproduction-driver hashes remain its scientific basis. Current scripts use the invoking interpreter's installed Engine and are not authority to relabel that historical snapshot.

The original 300 s Julia first-compilation cap was exhausted; the lead authorized **1800 s**, retaining normal optimized compilation. Warm per-case limit remains 180 s and total warm campaign 1200 s. That original campaign retired its explicit heavy lock. Current heavy runs use `agent-heavy` inside `/tmp/epc-shared-env.lock`, acquired first; never acquire the environment lock while holding the heavy lock. Runs require the accountable lead's exclusive heavy-work lease and are serial except FeOs plus its contraction helper, pinned together to the same CPU (at most two computational processes). Thermodynamic solves/gradients use native Engine C++, FeOs Rust or Julia; source automatic differentiation, diagnostic density refinements and reference fits run in Julia; Python only declares inputs, transports results and compares/plots them.

## Evidence ownership and limits

- [cases.json](inputs/cases.json) freezes source packets, component/site/pair order, solver selections, reference ladder and physical gradient coordinates; [states.csv](inputs/states.csv) distinguishes own and common states.
- [diagnostics.csv](results/diagnostics.csv) retains observable closures and available native diagnostics. Library success flags are not accuracy bounds; absent internals are explicit in [summary.json](results/summary.json).
- [density-closure-probe.csv](results/density-closure-probe.csv) preserves three historical Engine anchor outputs from commit `08c9ef8e`; its generating command is unavailable, so these rows do not qualify a remedy or replace primary outputs.
- [source-errors.csv](results/source-errors.csv) retains 256/512-bit source-function/derivative/refinement and ladder evidence; [decomposition.csv](results/decomposition.csv) and [screening-contrast.csv](results/screening-contrast.csv) retain ordered and same-density contrasts.
- C3 Born is off, εr is constant solvent-only, ion packing is 0.88σ in both screening selections. This does not adopt a new Held packet, fit parameters, establish global phase stability, test multiple-solvent permittivity, or establish Clapeyron parameter-gradient capability.
- Compiler/cold costs, interface extra work, three-process scatter and the specific processor prevent portable speed rankings.
- The 19 gate, build-restore, render and prose-check logs (`results/density-closure-*.log.gz`, `results/density-paired-*.log[.gz]`) are at record commit `21d7a112`.

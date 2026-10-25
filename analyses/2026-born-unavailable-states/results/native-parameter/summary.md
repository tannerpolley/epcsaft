# Native straight EOS/Born parameter continuation

The native declared-branch controller follows physical EOS/Born parameters with
exact first parameter actions and the existing full-coordinate Ipopt corrector.
The two copied reacting bubble requests retain their original source references,
feeds, supports and 313.15 K temperature. Pressure is free. No parameter fitting,
physical validation, global-stability search or turning-point method is claimed.

The shared family view handles scalar and vector coefficients; Born
`c_dielectric` and `ionic_region_relative_permittivity` now use the active scalar
record. Trial models are rebuilt from immutable base records with fixed
thermochemistry and remain alive through all compiled evaluations and anchored
root checks. The implicit tangent solves
$F_y y_s=-F_\theta(\theta_1-\theta_0)$ on
$\theta(s)=\theta_0+s(\theta_1-\theta_0)$, with $s\in[0,1]$.
The EOS/Born, packing, association and source-reference parameter chains are
exact. A partial prefix is not endpoint success.
The mathematical owners remain ALG-EQUILIBRIUM, ALG-ACTIONS, ALG-DENSITY and
ALG-TERMINATION, and EqIDs `implicit_sensitivity`, `born_radial_charging_work`
and `born_residual_helmholtz`. Their physics is unchanged.

## Current qualification

The main `8f5f5685` refresh reproduces the approved numbers on merged source `4d6118b0`.
The [research notebook](../../notebook.qmd#native-physical-parameter-continuation-after-the-main-refresh)
records the unchanged prefixes, falsifiers and limits;
[current-main-refresh.json](current-main-refresh.json) identifies its wheel and
retained replay hashes. The performance measurements below retain their stated
older source and wheel identities.

## Numerical comparisons

Tables report scalar falsifier maxima and finite last attempts. A plotted
boundary would not clarify these checks and could imply a fold location not
established by the calculations. Historical chart/turning tables are not used.

| Falsifier | Measured error or residual | Required criterion |
|---|---:|---:|
| 24 frozen-state Born values/columns, both copied anchor compositions; c=0,1 and epsilon=2,8 | Maximum scaled error 1.83e-16 | Relative 1e-12; absolute 1e-12 at zero |
| Manufactured six-identical-ethane path, 191.42→201.42 K dispersion energy, 200 K and 5 MPa | Phase amount error 8.88e-16 mol; composition error 4.44e-16 | Absolute 1e-9 |
| Manufactured amount/composition tangent per unit s | Maximum absolute action 4.25e-17 | Absolute 1e-9 |
| All accepted Born raw and classified rows, six paths | Maximum absolute residual 2.84e-13 | Absolute 1e-10 |
| 4,500 common historical-point quantity comparisons | Maximum comparison error 3.34e-13 | Relative 1e-7 for P/rho and nonzero Born/DH/permittivity; absolute 1e-7 for x, 1e-10 for zero terms |

`born-columns.csv`, `manufactured-exchange.csv` and
`parameter-path-comparison.csv` own those comparisons. Manufactured compositions
and amounts come from the analytic exchange laws, not a second fitted solution.
Every manufactured accepted state also meets the original-row 1e-10 check;
the maximum retained raw residual is 1.11e-14.

Centered native re-solves verify full row and solved-state first actions at the
smooth accepted c=0.5 and epsilon=3 states. Guesses come from retained historical
states, and the native driver recomputes the central state before differentiation.
Differences are verification only; no Python solver or corrector is present.

| Physical step | Case 1 maximum scaled action error | Case 2 maximum scaled action error |
|---|---:|---:|
| 1e-3 | 3.184e-5 | 5.635e-7 |
| 1e-4 | 3.184e-7 | 5.249e-9 |
| 1e-5 | 4.699e-9 | 3.707e-9 |

The required 1e-5 criterion is met on the refined rungs. Both central systems
have rank 20. Native condition estimates are 269.834 and 316.552; maximum linear
backward residuals are 1.425e-17 and 1.809e-17. Central residuals are 5.684e-14.
The native propagated inner-evaluation row floors are 1.801e-10 and 1.833e-10;
these are numerical estimates, not rigorous bounds on the solved-state derivative
error. `parameter-actions.csv` and `parameter-action-diagnostics.csv` retain the
full values and diagnostics. No physical parameter uncertainty is propagated.

## Finite continuation limits

The six paths use the requested coarse/fine/refined physical spacings through
one native continuation call per declaration. Native anchored-density admission
can halve a step; no Python sequence of cold restarts replaces continuation.
All final failed attempts return native `infeasible` causes and retain their rows.

| Path | Last accepted physical parameter | Final attempted parameter | Achieved failed residual |
|---|---:|---:|---:|
| Case 1 coarse | 0.331250 | 0.3296875 | 1.159e-5 |
| Case 1 fine | 0.330000 | 0.329375 | 3.847e-5 |
| Case 1 refined | 0.330000 | 0.329000 | 7.069e-5 |
| Case 2 coarse | 2.0302734375 | 2.030029296875 | 2.927e-6 |
| Case 2 fine | 2.030625 | 2.030000 | 5.212e-6 |
| Case 2 refined | 2.030500 | 2.030000 | 5.211e-6 |

Compared with the historical first-failed-step procedure, native predictors and
retries retain additional accepted epsilon points below 2.031. Both refined
historical accepted points remain available; both next historical targets fail
the unchanged 1e-10 criterion. The old failures and state tables are preserved.
The native failure residual need not match a historical cold/warm corrector's
rejected iterate. Different rejected iterates are not different physical roots.
Neither path reaches its requested c=0 or epsilon=2 endpoint. The returned
prefixes and bounded retries do not identify a mathematical fold or exclude
disconnected solutions. Target existence remains undetermined, as in the
historical diagnosis. Held-composition and turning-point calculations remain
retired and are not completion prerequisites.

Each `coarse/`, `fine/` and `refined/` directory retains state, coordinate, raw and
classified row, Born/DH/permittivity and trace CSVs, with compact path causes and
work counts in JSON. The native trace export reads the span list once; rebuilding
it for every event caused two export timeouts, not slow equilibrium correction.
The post-merge path calls took approximately 9–177 s under shared host load.

## Declared-solve performance

The two frozen #223 declarations use exactly the same inputs and initial guesses
on main and this branch. Solves run at `nice -n 10` with OMP and OpenBLAS threads
set to one. The final comparison consists of three independent main/branch
pairs on CPU 11, ordered main/branch, branch/main, main/branch. Each process runs
one solve per declaration (`--runs 1`); model construction is outside the timed
interval. Wall times below are the best of the three solves. All repetitions
have the reported iteration counts.

| Frozen declaration | Main Ipopt iterations | Branch Ipopt iterations | Main best wall time (s) | Branch best wall time (s) |
|---|---:|---:|---:|---:|
| Ionic LLE, Figure 2 point 7 | 16 | 16 | 0.0560403 | 0.0568194 |
| Reacting MEA flash, 318.15 K / 5918.5286 Pa | 43 | 43 | 0.0997766 | 0.100138 |

Inputs, returned states, iteration/evaluation/Hessian work are unchanged. Against
main `875526d0`, the final matched best-of-three wall times are 1.39% higher for
LLE and 0.362% higher for MEA. These small increases are reported, not discarded
or classified as zero. Three shared-host samples do not establish whether they
are a systematic regression; a strict no-wall-time-increase claim is not
supported. No speedup or universal non-regression guarantee is claimed. Earlier sequential cohorts
also include slower branch measurements; none has been discarded. Main/branch
best times on CPU 0 were 0.0859414/0.0923852 s (LLE) and
0.125359/0.174775 s (MEA); the subsequent CPU-11 sequential cohort gave
0.0572574/0.0615260 s and 0.108177/0.112173 s. This large between-cohort
variation motivated the paired, alternating-order comparison rather than a
comparison against an older host-load sample. Shared-host load and CPU frequency
remain measurement limitations. `performance/*-runs.csv` retains every cohort,
including the initial and unpaired samples, with residuals and native work.
The final paired tables also retain process CPU times. The JSON records retain
the frozen inputs, reaction constants, states, packet hashes and wheel identities;
identical per-pair JSON duplicates are removed after equality verification.

The current main-reference wheel was built from archived `875526d0` inside this
worktree. It was installed only sequentially into the same environment; the
branch non-editable wheel was restored and verified after the comparison.
Main wheel SHA-256 is `08b49ecf7ea2334e49740cf7089f3f464ebdfa6a22238bc1716adafcb3f5abf6`.
The branch wheel SHA-256 for that performance cohort is `1d7ae6add23f2c272476408fa7e2a769fd28dfb45579d85457f8b766033912c5`.
All six paths, manufactured rows, full first actions, strict Born columns,
global scalar chart/work identity and gates were rerun after merging this main.
Earlier `f593e803` timing cohorts remain retained; their wheel identities are in
`performance/wheel-snapshots.json`. `wheel.json` identifies that older qualification
wheel; the current source/wheel identity is in `current-main-refresh.json`; `inputs/manifest.json` above this result set owns the exact copied inputs.

## Scope and delivery

The existing native controller replaces the Python continuation loop in
`reproduce.py`; its allowlist entry is removed. Original frozen inputs and all
historical numerical records remain. #234 batch 15 owns consolidation of the
historical README/notebook prose; this result record does not revive its retired
chart scripts or establish their current execution. No main-checkout edit or
second Engine installation was made.

Runtime growth against `875526d0` is 103 gross added, 60 removed, +43 net, giving
21,732 lines against the recorded 24,000 envelope. The code-budget closeout
records the per-layer counts, required exact chains and reuse. No ceilings were
edited. On merged source commit `dd5a013d`, `bash tools/check-all.sh` passed
(115 Python checks and six native programs). `bash tools/check-affected.sh
origin/main` passed, including 115 Python checks and native Equilibrium,
Boundary and Parameter-regression-fit. Its base was `875526d0`. The Boundary
program retains the T/P/feed continuation falsifiers. The first cold post-merge
build reached all 253 linked targets but exceeded its 600-second cap; the warm
resume installed the current wheel and completed both gates and all replays. Explicit Born and manufactured falsifiers passed. Global Born
fitting coordinates now have no component targets, and native code owns their
physical affine chart and exact map. Seven retained comparisons in
`global-born-chart.csv` verify the coefficients/map (absolute 1e-15) and the
fixed-state Born correction ratio 16/35 (relative 1e-12); maximum absolute error
is 5.55e-17. A nonzero shell width is declared and the c=0 contribution is
subtracted. Born-column records enforce relative 1e-12 for every nonzero
reference, with an absolute 1e-12 criterion only at zero. The no-movement diagnostic
keeps its original intensive/feed meaning and names a stationary parameter path
separately; this correction restored the existing invalid-path check without
changing a numerical result.
[Independent delivered review](review.md) returned **Supported — scientific
#240 capability** for merged candidate `853574a5` against `875526d0`.
Both findings (global scalar semantics and literal Born tolerances) are closed.
Strict no-wall-time-increase is not supported; the measured timing increases and
shared-host limits remain explicit. This is numerical
verification for the stated declarations, not physical validation or manuscript
promotion.

## Reproduction

Run from the task worktree, with its installed `.venv/bin/python` and one thread:

```bash
out=analyses/2026-born-unavailable-states/results/native-parameter
script=analyses/2026-born-unavailable-states/scripts/reproduce.py
# Six independent declarations; each invocation contains one native path.
.venv/bin/python "$script" --case 1 --continuation --out "$out/coarse"
.venv/bin/python "$script" --case 2 --continuation --out "$out/coarse"
.venv/bin/python "$script" --case 1 --continuation --from-value 0.35 --step 0.01 --out "$out/fine"
.venv/bin/python "$script" --case 1 --continuation --from-value 0.33 --step 0.001 --out "$out/refined"
.venv/bin/python "$script" --case 2 --continuation --from-value 2.25 --step 0.01 --out "$out/fine"
.venv/bin/python "$script" --case 2 --continuation --from-value 2.04 --step 0.001 --out "$out/refined"
.venv/bin/python analyses/2026-born-unavailable-states/scripts/verify_parameter_continuation.py
```

The executed commands had explicit timeouts, `nice -n 10`,
`OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`. Builds and installed-wheel
gates used `/tmp/t3-engine-build-slot.sh timeout <cap> nice -n 10 env
OMP_NUM_THREADS=1 EPCSAFT_BUILD_JOBS=2 <command>`. The performance script
`benchmark_declared.py --label <main|branch> --wheel <identified-wheel>` reads
the same frozen #223 declarations and records three native solves per case by
default. For the final comparison, use `--runs 1` in each of three paired
main/branch processes, alternating their order as stated above, and retain every
sample. Process CPU time accompanies wall time; neither timer defines a
thermodynamic acceptance criterion.

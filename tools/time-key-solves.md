# `time-key-solves`: native work and paired solve time

`tools/time-key-solves <prepared-A> <prepared-B> [--pairs N] [--expect-change] --out <new dir>`
compares deterministic native work before reporting wall time. Each input selects
an existing Python interpreter, virtual environment, or prepared build root with
a `python` marker. Prepare both environments before invoking the command; it
verifies their loaded Engine identities without installing packages or building
wheels. Use the same prepared environment twice for an A=A control.

**Cases.** The declarations are reused, not selected here: the ordinary cases
from `analyses/2026-tpd-phase-discovery/scripts/baseline.py`, the eight ordered
feeds and previous-result sequence from
`analyses/2026-khudaida-idaes-settler/scripts/run_fixed_settler.py`, and the
unloaded-MEA declaration in `engine/tests/equilibrium/test_trace_reliability.py`.
Ionic liquid–liquid equilibrium is Khudaida Figure 2 point 7 (293.15 K,
100000 Pa); the reacting MEA flash is at 318.15 K and 5918.5286 Pa. Each runs
cold and from a common accepted coordinate vector taken from environment A and given
to both environments as one declared native start. The Khudaida sequence starts cold
at point 1 and passes each environment's previous result to its next solve through
point 8. Unloaded MEA is at 353.15 K and 101325 Pa with its existing composition
guess, reactions and neutral reference.

**Timing.** Each ordinary case has one complete warmup per environment. Each pair holds an AB
and a BA block, alternating which comes first. Workers share one CPU
(`taskset`), one numerical thread and niceness 10 (inherited niceness above 10
is rejected). Times sum equilibrium-call wall times only, excluding imports,
declaration and start preparation, messages and output reading. A block's
change is `100 (B/A - 1)` percent; the two blocks of a pair are averaged, and
the tool reports median, range and median absolute deviation. Negative means B
is faster. Timing on a shared host is supporting evidence, not a causal speedup.

**Outputs.** `work.csv`: iterations, native evaluations (including exact polish
work, not every equation-of-state evaluation), Hessian calls and achieved
raw/classified residuals per solve and sample. `states.csv`: returned T (K),
P (Pa), phase amounts (mol), mole fractions and molar densities (mol/m³).
`starts.csv`: the common accepted coordinates. `pairs.csv`: block times, load
averages and scaled state differences. `summary.json`: prepared identities under
`environments` (source revision, configuration and loaded-image paths and hashes),
input hashes and the comparison summary. Analyses retain `summary.json` only.

**Comparison and exit.** Zero exit is the no-slowdown gate. The tool exits
nonzero when any iteration, evaluation or Hessian-call count of ionic LLE
(case-6/figure-2-point-7) or the MEA flash (case-7/flash-5918.5286-Pa), cold or
accepted start, differs between A and B or between repeated samples, printing
each `WORK CHANGE` with both environments' counts. `--expect-change` prints the same
lines and does not fail on them, for a slice that intends a work change; the
other two cases' work changes are always printed only. It also exits nonzero
when any returned state value differs by more than 1e-10 of that quantity's
largest magnitude in the same returned state (both environments), in the declared
phase order: 1e-10 relative at the quantity's scale, with that scale as the
absolute floor for near-zero entries such as trace mole fractions. The worst
row is printed as `WORST STATE ROW`; `--expect-change` does not waive it.
Array-length or point-order differences stop the comparison. Paired wall times,
the median change and `zero_within_spread` (the absolute median within the
largest deviation from it) are information only. When an environment's Engine
lacks a type, enum member or field that a frozen fixture sets, that case is
printed as `NOT COMPARABLE <case>: <A|B> lacks …`, recorded under
`not_comparable` with no rows from either environment, and is a gate failure
(with its `MISSING STATE` line); the other cases still run. `summary.json` records
`gate_failures`; `tools/refresh-main-wheel` records a merged change, including a
`NOT COMPARABLE` line, in its README instead of stopping.

The harness self-check (an A=A run on the then-current main wheel) was retained
in `analyses/2026-key-solve-timing/` until #341 moved this manual here; it
remains at record commit `21d7a112`.

**Reacting F1 evaluation (#497).** The seventh case is the frozen seven-coordinate
MEA F1 packet in `tools/fixtures/mea-f1.json.gz`: 47 pressure rows, 94 species
rows and eight finite-addition enthalpies. Python reconstructs these declarations;
the Engine calculates their predictions, residuals and exact Jacobian together
at the packet's starting point. No application checkout is read. The identity
file records the original source hashes, compressed hash and uncompressed packet
SHA-256 `862ecf1e5d2756b7639493ea9d08444e890c66af6f427e5d9ee436910783a8db`.

F1 runs serially in A/B pairs, outside the ordinary warmup/AB/BA timing repeats.
The same `work.csv`, `states.csv`, comparison and `summary.json` retain it.
Its four work counts are reference computations after cache lookup (misses),
computed `reference_column` calls, `solve_density` calls, and the sum of
driver iteration ordinals from `ipopt_iteration` events, summed across every
attempt (including failed attempts and restoration). Missing iteration evidence
fails. A decrease passes; `--expect-change` waives only an increase from A to B.
Repeated counts within either environment must agree, including timer-on/off
samples, even with `--expect-change`. A missing
or incomplete native trace, unavailable row or nonfinite output stops the run.
Predictions and derivatives use the same `STATE_TOLERANCE`: entrywise relative
error with the packet's observation scales as absolute floors; physical
parameter derivatives divide that floor by the coordinate scale. Native scaled
residual derivatives use a unit floor. This prevents pressure magnitudes from
hiding changed species or heat outputs.

Serial CPU is process CPU around the native `evaluate_observations` call with
non-diagnostic trace capture. Declaration decoding, trace reduction, process
startup and output writing are excluded. Candidate B must use at most 1.25 times
the fixed CPU value in the fixture identity file; `--expect-change` cannot waive
this ceiling. The frozen value does not follow later main timings upward.
A baseline A, or either arm of an A=A control, measures the host only when its
four work counts equal the identity file's frozen counts; an arm doing other work,
such as a baseline from before a work change, is not compared with the ceiling (#595).
A host-measuring arm over the ceiling invalidates the host measurement. The command repeats
the serial A/B pair up to three times, sequentially, to obtain a valid baseline.
If all three controls exceed the ceiling, it exits nonzero with
`invalid host measurement: rerun on a quieter host`. This is not a ceiling
failure or a pass. With a valid control, candidate B must still meet the same
absolute ceiling. Every attempt's work and states remain subject to the work
and numerical comparisons; `serial_attempts` retains both CPU times per pair and fixture.
The historical 21.79 s measurement on `36485a88` is **wall time**, not CPU.

**Reacting Cp evaluation (#595).** The eighth case is the frozen
`tools/fixtures/mea-cp.json.gz`: the Hilliard (2008) equilibrium isobaric heat capacity
row at 353.15 K, 101325 Pa and CO₂ loading 0.358 (feed MEA 1, water 7.929776437398854,
CO₂ 0.358 mol; liquid only, five reactions) at the #170 Route 2b larger-dV start of 16
coordinates. Its value and exact derivatives need the solved state's third temperature
and mixed temperature-parameter reference columns. The packet carries the CO₂/MEA
`no_cross_association` declaration inside its hash; the identity file records both
source hashes. It runs the F1 serial A/B evaluation, work counts, numerical comparison,
host rule and CPU ceiling, without a fit arm; an A=A `--fit` control checks its serial
evaluations only. Its 20 physical columns (17 active parameters and three reaction-law
directions) do not map one-to-one to coordinates, so a physical derivative's floor is
the row's largest physical derivative. Frozen work is **3 / 280 / 35 / 13** in the F1
count order; serial CPU **1.978993 s** (ceiling **2.473741 s**) on the reference host.

**Full fit.** `--fit` uses the same prepared environments and worker protocol
for one complete fit per arm, with 40 iterations maximum, a 120 s time budget,
and native default workers (32 on the reference host). The ordinary samples
and F1 work measurement run with the timer off. Fresh workers then enable
`EPCSAFT_FIT_CPU_PROFILE` and repeat F1 once to check identical work and numerical
values before fitting. Each fit must converge and return usable, finite numerical arrays
and available observations without trial failures, and finish within 120 s
wall. The numerical comparison and wall ceiling cannot be waived.

For A=A, invoke `--fit` without a control to retain two same-image fits. For
A/B, supply that passing summary as `--control <summary.json>`; its host, A
images, packet and worker count must match. The observed A=A wall range is
printed alongside the two absolute-ceiling decisions. A control that itself
breaks the ceiling, lacks either serial or fit arm, has nonfinite required
times, or reports a failed fit is not eligible; host contention never waives
the limit.

`fit-cpu-A/` and `fit-cpu-B/` contain native per-process rows:
`evaluation,worker,jacobian,cpu_seconds`. The timer is off by default, with
one branch in `compute` before its unchanged calculation. When enabled, it
uses `CLOCK_PROCESS_CPUTIME_ID` around that calculation; file opening and output
writing lie outside the interval. The #405/#498 quantity is
$C_{\mathrm{fit}}=\sum_e\max_w t_{e,w}$, the sum of each evaluation's slowest-worker
compute CPU, not the slowest worker's whole-fit total. Missing, duplicated or
nonfinite worker/evaluation records stop the run. This excludes controller,
startup, communication and merging outside `compute`; total fit CPU including
child processes and total compute CPU are reported separately. No instruction
counts or physical-validation claim follow from these measurements.

**Frozen reference on main `5ad73fb9`, host `tnnrpolley21-MS-7C95`.**

| Quantity | Issue baseline `36485a88` | Refrozen current main |
|---|---:|---:|
| Serial evaluation wall, s | 21.79 | 16.282463 |
| Serial evaluation CPU, s | not recorded | 16.279413 |
| Central Taylor evaluations | 664 | 675 |
| Fit wall, s | 46.4 | 31.264311–31.566555 |
| Whole-fit CPU, s | 833 | 620.181655–631.669772 |
| Fit iterations | 14 | 14 |
| Residual / Jacobian evaluations | 15 / 11 | 15 / 12 |
| Critical-path compute CPU, s | not recorded for this F1 packet | 28.646385–29.049397 |

Both fits used 32 workers and converged. The fixed serial CPU ceiling is **20.349267 s**. Current work per evaluation is **295 / 3388 / 2627 / 1426** in the four-count order above, refrozen by #595, which stopped computing unread reference columns; the serial CPU value is unchanged. The input identity file owns these values and their source and measurement identities.

The main A=A command with `--fit` took 143.596 s; the default main-versus-candidate command took 39.391 s. Full fitting is explicit because it repeats the fixed-point work throughout the optimizer, while the default command already checks the four work counts and serial ceiling on every runtime change. The Gates table requires full fits for changes that can change this repeated work or worker distribution.

These are performance verification measurements on a shared desktop. The two current fits give an observed range, not a confidence interval. All intervening main changes separate the historical and current values; the difference is not attributed to this observation-only change. The original six-case gate also compared unmodified main with the instrumented baseline, returning identical work counts and states with the timer off.

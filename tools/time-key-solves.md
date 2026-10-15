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

**Timing.** Every case has one complete warmup per environment. Each pair holds an AB
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
largest deviation from it) are information only. `summary.json` records
`gate_failures`; `tools/refresh-main-wheel` records a merged change in its
README instead of stopping.

The harness self-check (an A=A run on the then-current main wheel) was retained
in `analyses/2026-key-solve-timing/` until #341 moved this manual here; it
remains at record commit `21d7a112`.

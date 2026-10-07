# `time-key-solves`: native work and paired solve time

`tools/time-key-solves wheel-A.whl wheel-B.whl --pairs N --out <new dir>`
compares deterministic native work before reporting wall time. Each wheel is
installed with `uv` into a throwaway environment under `build/time-key-solves/`;
the checkout's `.venv` supplies only the Python version and is never modified.

**Cases.** The declarations are reused, not selected here: the ordinary cases
from `analyses/2026-tpd-phase-discovery/scripts/baseline.py`, the eight ordered
feeds and previous-result sequence from
`analyses/2026-khudaida-idaes-settler/scripts/run_fixed_settler.py`, and the
unloaded-MEA declaration in `engine/tests/equilibrium/test_trace_reliability.py`.
Ionic liquid–liquid equilibrium is Khudaida Figure 2 point 7 (293.15 K,
100000 Pa); the reacting MEA flash is at 318.15 K and 5918.5286 Pa. Each runs
cold and from a common accepted coordinate vector taken from wheel A and given
to both wheels as one declared native start. The Khudaida sequence starts cold
at point 1 and passes each wheel's previous result to its next solve through
point 8. Unloaded MEA is at 353.15 K and 101325 Pa with its existing composition
guess, reactions and neutral reference.

**Timing.** Every case has one complete warmup per wheel. Each pair holds an AB
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
averages and state differences. `summary.json`: wheel and input hashes and the
comparison summary. Analyses retain `summary.json` only.

**Comparison.** State differences use `abs(A-B)/max(abs(A),abs(B))` in the
declared phase order (two zeros differ by zero); array-length or point-order
differences stop the comparison. Every work-count difference is marked,
including between repeated samples. For identical wheels the tool succeeds only
when work is identical, the largest state difference is exactly zero, and the
absolute median timing change is within the largest observed deviation from
that median in each case (a spread check, not a confidence interval).

The harness self-check (an A=A run on the then-current main wheel) was retained
in `analyses/2026-key-solve-timing/` until #341 moved this manual here; it
remains at record commit `21d7a112`.

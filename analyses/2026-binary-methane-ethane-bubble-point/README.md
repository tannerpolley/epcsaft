# Binary methane/ethane PC-SAFT bubble-point comparison

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The original
cross-package measurements require their original Engine and comparator environments;
external packages were not rerun by the later Engine refresh.

This folder preserves the original cross-package run. The
[current Engine refresh](../2026-open-source-equilibrium-comparison/README.md)
rechecks the same case on wheel `4268f839…`: pressure is 1596859.1113027686 Pa,
with a 13.035 ms solve median (p10–p90 9.575–15.685 ms; 25 samples).
The tables and figures below are historical, not measurements of that wheel.
The noisy fresh timing does not establish a speedup or regression against the
unpaired historical 10.021 ms result. External packages were not rerun.

This analysis uses the Engine's source-backed `gross-2001-methane-ethane`
packet as a more demanding comparison than pure-fluid saturation:

- neutral binary PC-SAFT, no association, polarity, or electrolyte terms;
- methane/ethane parameters from Gross–Sadowski (2001);
- `k_ij = 0` as retained in the Engine packet;
- methane `(m, sigma, epsilon/k) = (1.0, 3.7039 Å, 150.03 K)` and ethane
  `(1.6069, 3.5206 Å, 191.42 K)`;
- liquid feed composition `x = z = [0.5 methane, 0.5 ethane]`;
- isothermal bubble point at `T = 180 K`;
- outputs: bubble pressure, incipient vapor composition, and both phase
  molar densities.

The case is solved by ePC-SAFT, FeOs, and ThermoPack from the local forks.
The external package comparison is intentionally limited to packages that
expose a direct multicomponent bubble-point calculation for this formulation.

## Performance protocol

Each package is imported once in its own process. Model construction is timed
separately. The equilibrium call is warmed up five times and then timed for 25
repetitions using `time.perf_counter_ns()`. The reported solve statistics are
wall-clock milliseconds on the same host, with `FEOS_MAX_THREADS=1` and
`OMP_NUM_THREADS=1`. Import/interpreter startup is excluded from both setup and
solve timings.

## Results

The exact machine-readable result is in `results/values.csv`; provenance,
source commits, commands, and timing samples are in `results/provenance.json`.
The figures are rendered only from that retained CSV.

At this state, the three implementations agree closely:

| implementation | bubble pressure [MPa] | vapor y(CH4) | liquid density [kmol m^-3] | vapor density [kmol m^-3] |
| --- | ---: | ---: | ---: | ---: |
| ePC-SAFT | 1.596859111303 | 0.956360242961 | 19.6985002963 | 1.28845507053 |
| FeOs | 1.596859111206 | 0.956360242976 | 19.6985002978 | 1.28845507044 |
| ThermoPack | 1.596859111303 | 0.956360242961 | 19.6985002963 | 1.28845507053 |

The largest cross-implementation pressure difference is about `6.1e-11` in
relative terms; the Engine equilibrium solve reports a maximum absolute
residual of `3.67e-14`.

The measured warm-solve median and model-setup times were:

| implementation | warm solve median [ms] (p10–p90) | setup [ms] |
| --- | ---: | ---: |
| ePC-SAFT | 10.021 (9.744–12.010) | 0.577 |
| FeOs | 0.215 (0.207–0.220) | 6.915 |
| ThermoPack | 0.487 (0.454–0.621) | 5.308 |

These timings compare the exposed case-level calls, not identical internal
algorithms: ePC-SAFT uses the general `bubble_point` plus
`solve_equilibrium` driver, while FeOs and ThermoPack use their specialized
bubble-point routines. They are useful for this workflow on this machine, but
not a portable package benchmark.

## Limits

This is a declared two-phase bubble-point calculation, not a global phase
stability search. It tests binary phase splitting and incipient-vapor
composition, but not association, reaction equilibrium, electrolyte chemistry,
or LLE. Runtime values are environment-specific performance evidence, not a
portable benchmark claim.

# Neutral methane/ethane PC-SAFT TP-flash comparison

The integrated Engine's canonical cold median is **14.30 ms versus 41.09 ms**
on the repaired baseline; previous-result reuse is **6.62 ms versus 21.42 ms**.
All 270 integration controls pass with the repaired branch and multiplier
checks retained. See the [integration record](notebook.qmd#integration-with-repaired-admission)
for the actual wheel identities, numerical changes and limited input domain.

The [performance notebook](notebook.qmd) retains the original local performance
study and its wheel identities. Its 6.80 ms versus 29.32 ms cold-solve comparison
predates the density-branch and ionic-multiplier repairs; it is historical
evidence, not a timing claim for the integrated Engine. The cross-package
timings below are also retained historical measurements.

This is a true isothermal-isobaric flash, not a bubble or dew calculation:

- Gross–Sadowski (2001) neutral PC-SAFT, methane/ethane, `k_ij = 0`;
- `T = 180 K`, `P = 1.0 MPa`, and overall feed `z = (0.5, 0.5)`;
- both phase compositions and the vapor fraction are solved or reconstructed
  from the solved two-phase state;
- component order is `(methane, ethane)`.

The case is inside the two-phase region. The ePC-SAFT result has distinct
liquid and vapor compositions and `beta_V = 0.3261546790`, so it is not a
boundary-point calculation.

## Values

| implementation | vapor fraction | liquid x(CH4) | vapor y(CH4) | liquid density [mol m^-3] | vapor density [mol m^-3] |
| --- | ---: | ---: | ---: | ---: | ---: |
| ePC-SAFT | 0.326154679039 | 0.294537595547 | 0.924491472211 | 19234.211443 | 749.100424 |
| Clapeyron.jl | 0.326154679072 | 0.294537595516 | 0.924491472211 | not captured | not captured |
| FeOs | 0.326154679020 | 0.294537595549 | 0.924491472245 | 19234.211445 | 749.100424 |
| teqp | 0.326154679039 | 0.294537595547 | 0.924491472211 | 19234.211443 | 749.100424 |
| ThermoPack | 0.326154677926 | 0.294537596588 | 0.924491472211 | 19234.211446 | 749.100424 |

The four implementations with captured densities agree to better than about
`1.1e-9` in methane composition and `1.5e-6 mol m^-3` in liquid density.
Clapeyron's completed flash agrees at the same level for the quantities it
returned. The complete machine-readable table is in
[`results/values.csv`](results/values.csv), with provenance in
[`results/provenance.json`](results/provenance.json).

## Runtime

Warm solve medians used 3 warm-ups and 15 timed calls in a single process per
package. Imports and model construction were excluded. Each timed process was
pinned to one CPU with nice level 19 and BLAS/OpenMP thread counts set to one.

| implementation | solver call | median [ms] | p10–p90 [ms] |
| --- | --- | ---: | ---: |
| ePC-SAFT | general `solve_equilibrium` TP driver | 30.799 | 28.223–34.001 |
| FeOs | `PhaseEquilibrium.tp_flash` | 0.412 | 0.284–0.452 |
| teqp | `mix_VLE_Tp` plus lever rule | 0.199 | 0.162–0.236 |
| ThermoPack | `two_phase_tpflash` | 0.800 | 0.640–1.035 |
| Clapeyron.jl | `tp_flash` | not measured | bounded first specialization exceeded safe memory cap |

These are exposed case-level timings, not an algorithm-neutral benchmark.
ePC-SAFT uses the general derivative-aware equilibrium driver, whereas FeOs and
ThermoPack use specialized TP-flash paths. teqp's Python API does not expose a
single feed-aware TP-flash call in this checkout: its native `mix_VLE_Tp`
solves the two coexisting phases at fixed `T,P`, and the vapor fraction is then
obtained from the material balance for the specified `z`. That is still a
two-phase TP flash, but the distinction matters when interpreting runtime.

teqp's Newton result had a maximum phase residual of `4.82e-11`; its Python
wrapper reported `success=false`/`return_code=unset` because the Newton branch
does not set those status fields even when the returned residual is finite and
stationary. This is recorded rather than silently converting it to a success
flag.

## Resource guard and limits

No compiler or Julia process is intentionally left running. Clapeyron was not
retried without a memory guard after a bounded first method-specialization run
hit `OutOfMemoryError` with a 1.2-GB virtual-memory limit and 861616 kB maximum
resident set. Its state values are retained from a completed, converged flash;
its density extraction and steady-state runtime remain unmeasured.

This is one neutral binary PC-SAFT case. It does not qualify association,
electrolyte, reactive, LLE, global-stability, or derivative-action behavior.

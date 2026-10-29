### Current main against the merged candidate

Source A: `3cf744eed72b9b896c5a77150cf4ed0a96b8fb3d`; source B: `16ccdf4a69ecd19f193d36035617d389ec48ae2b`. Reporting edits leave the Engine source unchanged.

The timing verdict is **unresolved**: the median saving does not exceed the same-image control envelope, so no speed improvement is claimed. The change is retained as a numerically neutral cleanup.

| Prepared whole-call quantity | Result |
|---|---:|
| R | 0.881904095042 |
| negative_log_R | 0.125671964668 |
| C | 0.278960856467 |
| saving_median_us | 4.1995 |
| saving_min_us | -4.73 |
| saving_max_us | 15.559 |
| maximum_selected_load1 | 3.81 |

Six paired medians use three repetitions, two workers per image, five warmups and 25 samples per arm/block on logical CPU 31 at niceness 10. The middle ordering is reversed. Selected boundaries have load1 below four and no compiler observations. Whole disturbed operations are rejected by the unchanged reducer; none is replaced because of speed.

[Pairs and decision arithmetic](results/issue-865-r2/interim-summary.json), [block medians](results/issue-865-r2/blocks.csv), [selection](results/issue-865-r2/selection.csv), [comparator ratios](results/issue-865-r2/paired-ratios.csv), and image identities retain the measurements. All raw samples and deferrals remain in the main checkout under `build/issue-865/build-r2/cost-analysis/results/`.

| Numerical preservation | Result |
|---|---:|
| regular_states | 251 |
| regular_values | 10040 |
| states | 256 |
| values | 10135 |
| admission_cases | 26 |

The regular values are bit-identical, with zero quantity failures, diagnostic changes or admission differences. All admission decisions, terminal indicators, messages and limit counters agree in both phase orders. The successful unbounded path uses one cold search and both anchors; bounded calls retain the original sequence. Sole-root and injected-refusal cases preserve finite-search classification and failure precedence. The [comparison](results/issue-865-r2/equivalence-summary.json) and admission/root tables retain these operands.

All equilibrium corpora have zero differences: Cameretti 1800/1800, Gross 721/724, neutral 63/70, matrix 42/42, replay 4381 rows and MEA 73/73. Existing refusals remain refusals.

| Key solve quantity | Candidate | Ceiling |
|---|---:|---:|
| mea-f1/residual-jacobian serial CPU / s | 4.39442061 | 20.349267 |
| mea-f1/residual-jacobian computations/columns/density solves/iterations | 295/3388/1223/1426 | unchanged |
| mea-cp/hilliard-353.15-K-0.358 serial CPU / s | 0.607212338 | 2.473741 |
| mea-cp/hilliard-353.15-K-0.358 computations/columns/density solves/iterations | 3/280/11/13 | unchanged |
| Full fit A wall / s | 7.15662439 | 120 |
| Full fit B wall / s | 7.4412407 | 120 |

A=A fit wall range: 7.25955522–7.29538203 s. Fits retain full CPU affinity and 32 native workers; serial samples use CPU 31. Work and states agree, with zero state differences and zero central Taylor evaluations. [A=A control](results/issue-865-r2/key-control-summary.json.gz) and [A/B comparison](results/issue-865-r2/key-comparison-summary.json.gz) retain the operands.

At the time of this measurement, both the affected and the full checks stopped at the total runtime-line ceiling: 30059 / 30025, over by 34. Equilibrium was 11884 / 11950. The total ceiling has since been set to the measured count, 30059, and the current budget check passes: total 30059 / 30059, with every layer within its ceiling. Supplemental architecture, 238 Python cases, Ruff and strict type checking pass. The explicit native result and final frozen-head native qualification are reported in `build/issue-865/build-r2/build-summary.md`.

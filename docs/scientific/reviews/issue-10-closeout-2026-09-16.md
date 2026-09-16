# Issue #10 neutral IDAES prototype closeout - 2026-09-16

## Result

The complete nonequilibrium prototype is supported on current main base
`7ac882e0`. Clean source commit `9486fccb` built the one non-editable Engine wheel
and produced the retained campaign in `analyses/2026-idaes-asl-adapter/results/`.
The adapter changes no equilibrium code and performs no phase or reaction solve.

The wheel SHA-256 is
`2fc1e01ccb6953ceb93815036c6bf88aaf8a696ab0546b9019dc2ad975beb170` and
the installed ASL library SHA-256 is
`525ad4bca9116963315581eee06a5c189f5a2839430bdd42d079afb5c1743f0e`.
The library exports the named neutral residual-Helmholtz and pressure callbacks;
the process uses the retained propane packet fingerprint
`4593b9dbd813004096c7e396fc65df6ae3f4c1c65933c36918a6d316111326e2`.

## Numerical evidence

The 17.26-second campaign passed 54 conditioning states, seven independent
references, 72 scalar derivative ladders, 216 whole-NLP directional ladders,
27 explicit Hessian matrices and four optimal KKT cases. The minimum sampled
pressure slope divided by `RT` is 0.9498 against the 0.8 admission limit. Maximum
whole-NLP scaled Hessian error is `1.033e-9`; maximum KKT stationarity is
`1.72e-16` and maximum complementarity is `9.10e-11`.

The stable domain, derivative, whole-NLP and independent-reference tables are
identical to the pre-rebase accepted tables after line-ending normalization.
The changed wheel and ASL hashes reflect the current-main native archive, including
E1 work; the bounded adapter outputs did not change.

Four second-virial reference rows at 400 K and 0.001–1 mol/m³ are deliberately
native-only dilute-limit checks. They lie outside the adapter's sampled
5–110 mol/m³ domain, do not call ASL and do not expand the adapter claim.

## Presentation

`scripts/plot_results.py` renders five figures from retained CSV/JSON only:
the sampled pressure/conditioning envelope, derivative convergence and cost,
process convergence/KKT results, the actual throughput flowsheet and the compiled
adapter data path. `notebook.qmd` assembles those figures and retained tables into
a self-contained HTML report without rerunning physics.

## Limits

This closes the pure neutral, steady, isothermal, nonequilibrium example only.
It does not qualify association, electrolytes, mixtures, equilibrium, reactive
loading, phase stability, Lithium/PROMMIS behavior, equipment cost, heat duty or
compressor work. Those boundaries are now separated into issue #66 and formal
Lithium qualification issue #63.

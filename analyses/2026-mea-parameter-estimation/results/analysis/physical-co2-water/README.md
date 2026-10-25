# Corrected CO2-water physical-tier evidence

`scripts/physical_co2_water_run.py` (leaf `scripts/`) now declares native constant and reciprocal-temperature fits. The
[current result](result.json) and [comparison](batch2-comparison.csv) retain
objective and pressure agreement, exact columns, rank and bounds. A historical
unresolved coexistence trial now closes locally; global uniqueness and the
two-sided nuisance-refitted nonlinear profile remain unestablished.

The earlier fit used
wheel SHA-256
`d59abedb354fbd163ea918d3cc96b607f6e8ce8ec7c1f730d12b9653a1541e5d`,
Engine `0def56adc38e133f5fcd1ad031649f41529a8888`; its result was recorded in
commit `1ed0d0bc1af5915b910997f3b4ddff0ba97c1e1e` (#115).
Current MEA qualification uses the adopted application record and
[`2026-mea-recovered-state-actions`](../../../../2026-mea-recovered-state-actions/README.md)
(#61). Retaining a fitted coefficient in later input packets does not requalify
this older fit, source-reference treatment or uncertainty calculation.

## Historical selection

For that wheel and parameter family, the selected constant binary interaction was

\[
k_{\mathrm{CO_2,H_2O}} = 0.0132628791769.
\]

It uses the corrected Gross-Sadowski CO2 segment count `m = 2.0729`. The fitted
point is interior to the retained `[-1, 1]` bound, the exact sensitivity has rank
1/1, and all 39 retained Kiepe rows solve without a phase failure. Its nonlinear
one-dimensional 95% scaled-profile crossings are `0.0113745020` and
`0.0152445662`. The scale is `SSE/(n-p)` for unweighted log-pressure residuals;
these crossings are diagnostic model-fit bounds, not source-measurement confidence
limits.

The reciprocal-temperature form is retained only as a candidate:

\[
k(T) = 0.0215945542225 - 30.8969651510\,\mathrm{K}
\left(\frac{1}{T}-\frac{1}{353.11\,\mathrm{K}}\right).
\]

It is interior to the bounds (`k` in `[-1, 1]`, slope in
`[-1000, 1000] K`), has rank 2/2, and materially lowers AICc and residual error.
However, a reproduced off-optimum nuisance-profile trial at `k=0.02` and slope
`-28 K` reached `multiple_unresolved_coexistence_branches` at
`kiepe-313.20-13`. The retained
quadratic curves are local Jacobian/covariance uncertainty only; they are not
finite nonlinear profile bounds and are not used to promote the temperature
dependence.

## Deterministic results

| Model | Parameters | log-pressure RMSE | mean log bias | max absolute log error | final phase failures |
|---|---:|---:|---:|---:|---:|
| corrected `m=2.0729`, fixed `k=0` | `k=0` | 0.3620541961 | -0.3364744096 | 0.6067268150 | 0/39 |
| selected constant interaction | `k=0.0132628792` | 0.1853632802 | -0.0942751487 | 0.3989135365 | 0/39 |
| reciprocal-temperature candidate | `k_ref=0.0215945542`, slope `-30.8969651510 K` | 0.03596102365 | 0.00606557191 | 0.1133468770 | 0/39 |

For the selected constant interaction, log-pressure RMSE by temperature is
`0.1361456280` at 313.20 K, `0.1247887781` at 353.11 K, `0.2088356938` at
373.26 K, and `0.2751049435` at 393.17 K. By observed-pressure band it is
`0.1264751926` below 1 MPa, `0.1877355718` from 1 to 5 MPa, and
`0.2501811957` at or above 5 MPa. By liquid CO2 mole-fraction band it is
`0.1280790571` below 0.001, `0.1894507787` from 0.001 to 0.01, and
`0.2112638098` at or above 0.01. Exact row predictions and every grouped metric
are retained in `predictions.csv` and `result.json`.

The data-supported domain is 313.20-393.17 K, 95.09-9257.6 kPa observed
pressure, and liquid CO2 mole fraction 0.00003-0.02481. The six 393.17 K rows
are 0.02 K above the active candidate ceiling of 393.15 K. They were executed
under an analysis-only boundary extension and remain boundary evidence; this
packet does not expand the active model domain.

## Inputs and method

The historical calculation imported `epcsaft 0.2.0.dev0` from the non-editable installed
post-#103 wheel and hard-gated these inputs:

- wheel SHA-256: `d59abedb354fbd163ea918d3cc96b607f6e8ce8ec7c1f730d12b9653a1541e5d`
- materialized parameter document SHA-256: `12e380a27eb615e25bae580f8c1ec852020879917b08fa16a5821ae76311550d`
- retained Kiepe table SHA-256: `ec28659dc696643c49bca771b65c79c7852b4c2e24296c9aa79b555ffdc7b9bd`

The binary selection retains the corrected CO2 216-304 K fit range together
with unrelated reaction/model-family records qualified only at 313.15 K, so its
source domains have an empty intersection before candidate qualification.
Following the campaign's established replay path, the run script records every
candidate-provenance domain rewrite; physical coefficients and association topology
were unchanged. The historical calculation solved fixed-temperature incipient
coexistence through `PhaseEquilibriumProblem` and fitted unweighted
`log(P_calc/P_obs)` using exact public sensitivities.

The historical `m = 2.079` model was not executed and remains historical only.
No loaded-MEA observations or holdout data are read, and this folder does not
create or modify an active bundle.

## Retained files and data-only rendering

The retained predictions can be rendered from the repository root:

```bash
.venv/bin/python analyses/2026-mea-parameter-estimation/scripts/physical_co2_water_render.py
```

The original input/provenance record was removed in commit `9a7fc51a`.
It remains readable at
`9a7fc51a^:analyses/2026-mea-parameter-estimation/results/analysis/physical-co2-water/input-receipt.json`.
The native producer retains current input identities in `input-result.json`.

- `scripts/physical_co2_water_run.py`: native fit declarations and input identities.
- `result.json`: parameters, metrics, ranks, bounds, grouped evidence, profiles,
  status, and supported domain.
- `predictions.csv`: exact row-level plotted and audited values.
- `profiles.csv`: local quadratic two-parameter candidate curves; not nonlinear
  profile evidence.
- `pressure-evidence.svg` and `.png`: byte-deterministic parity and
  selected-model residual views.
- `scripts/physical_co2_water_render.py`: rendering from the retained tables, without Engine evaluation.

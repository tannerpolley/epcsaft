# Generic Pure-Saturation Parameter Fit

Pure-component parameter fitting uses the sole public Regression workflow:
`ParameterFit` resolves selected parameter identities, `fit` executes the
generic native residual engine, and `FitResult` retains every start and the
exact fitted parameter transfer.

The EOS owns coexistence through `Pure.saturation` and its exact selected-
parameter derivatives. Regression supplies no saturation solver, equation,
target index, paper-specific route, packaged literature table, or alternate
parameter record; it consumes the EOS-owned observable block.

For each admitted temperature, the generic `pure_saturation` observation may
compare saturation pressure and saturated liquid or vapor density. Positive
quantities use declared log-ratio residuals and explicit block multipliers.
Source units, basis, reference convention, transformations, validity domain,
and SHA-256 identity are required by the strict `ParameterFit` document.

The compact Engine tests falsify coexistence, derivative, deterministic-
multistart, non-evaluable-trial, rank, conditioning, active-bound, and transfer
failures. A future literature or predictive claim belongs in an explicit
installed-wheel Validation campaign. The retained
`validation/campaigns/2026-pure-saturation-regression/` request and receipt are
historical provenance: the final wheel produced no comparable quantitative fit
output, so no current campaign test or held-out predictive claim is admitted.

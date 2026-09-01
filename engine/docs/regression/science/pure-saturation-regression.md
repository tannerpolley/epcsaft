# Pure-saturation owner regression

Pure saturation is compiled by `epcsaft.eos.compile_regression_evaluator` with
an explicit temperature, ordered saturation row identities, and one
`ActiveParameterSet`. The EOS owner supplies the common pressure, selected
phase densities, and exact physical selected-parameter derivatives. Regression
only supplies observations, discrepancy weights, bounds, affine coordinates,
and Ceres controls.

The returned `LeastSquaresResult` keeps physical owner Jacobians separate from
the weighted optimizer-coordinate Jacobian, along with native failures,
resource counters, covariance availability, and solver options. No result
promotes or adopts parameters. Quantitative literature or held-out claims
remain installed-wheel analysis evidence.

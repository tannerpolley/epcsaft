# Parameter estimation

Declare observations and physical coordinates, then let the native Engine
evaluate predictions, exact first derivatives and the estimate. The
[propane notebook](../examples/propane-parameter-estimation.ipynb) estimates
three pure parameters from published saturation pressure and liquid-density
observations. The [estimation formulation](../theory/parameter-estimation.md)
and [fitting algorithm](../theory/algorithms.md#alg-fit) define the objective,
coordinate maps and diagnostic meaning.

## Observations and residuals

Specify the observation kind, temperature in K, pressure in Pa where fixed,
composition/feed basis and phase branch. Give the measured value in the
observation's declared units, including the mass or molar density basis.
A saturation pressure observation requires a boundary problem; a homogeneous
density observation requires a pressure-closed state. These are different
physical predictions even at the same nominal temperature.

Supported observation constructors include density, fixed-composition total
isobaric heat capacity, finite-addition enthalpy, source-basis reaction ln K,
and declared phase observations. Electrolyte constructors also expose selected
activity and osmotic quantities. Each kind has its own admission and reference
requirements; see [observation calls](../api/regression.rst) and the
[native estimation records](../api/native/parameter-estimation.rst).

Distinguish a difference residual from a log-ratio residual. A scale makes a
residual dimensionless according to its declared definition. An absolute
measurement standard deviation supports a different statistical interpretation
from a relative weight chosen to balance observations. See [propane calibration](../examples/index.md#propane-calibration)
for the example's weighting assumption.

## Coordinates, bounds and holds

A coordinate declares its physical family and component or reciprocal-pair
targets, origin, scale and bounds. The optimizer coordinate is dimensionless;
the scale maps it to native physical units. Reciprocal pair slots belong to
one unlike interaction coordinate. Charges, masses and discrete topology are
fixed inputs rather than fitting coordinates.

Reaction coordinates can select a named reaction's scalar ln K or supported
stored correlation coefficients while retaining the reference declaration.
Do not refit reference conventions implicitly. `HydrationHold` imposes a
single-ion solvation equality using a dependent Born diameter with its declared
bracket and tolerance. Its eliminated-variable chain and conditioning remain
part of the returned derivative evidence.

## Run a bounded native estimate

Set `FitControls` iteration and elapsed-time budgets, tolerances and declared
workers as appropriate to the calculation. Python supplies the request;
the Engine owns prediction solves, exact Jacobians and optimization. No Python
density search, equilibrium loop or numerical optimizer belongs around it.

See [propane calibration](../examples/index.md#propane-calibration)
for the bounded estimate and its interpretation.

## Inspect the point and its uncertainty assumptions

Read physical parameters, predicted values, residuals, objective and trial
failures, together with stopping reason and consumed work. A failed trial
retains its thermodynamic cause; do not replace it with an invented physical
prediction. Re-evaluate the final point using `fitted_parameters` and the
original observations.

Rank and singular values describe local sensitivity in the declared chart.
A covariance estimate is conditional on the stated error model, rank, active
bounds and local linearization. Full rank alone does not prove global
identifiability; a small fitted error does not give out-of-sample uncertainty.
Calibration is not validation. Independent observations and an intended-use
domain are needed for a physical validation claim.

## Minimax and measurement reconciliation

`MinimaxRequest` declares a local scaled epigraph estimate with target,
maximum-error and mean-error guards. Inspect feasibility, physical and scaled
inequalities, multipliers, stationarity and complementarity separately.
The implementation uses exact first columns with a limited-memory optimizer
Hessian; it does not promise exact parameter Hessians or a global optimum.
Minimax supplies no covariance interpretation.

Errors-in-variables estimation declares uncertain measured inputs and output
standard deviations in one Gaussian measurement objective. Shared measurement
groups and optional correlations must be supplied explicitly. Zero input sigma
fixes that input; the Engine does not infer measurement sharing from coincident
state values. Inspect reconciled inputs, standardized corrections and residuals
alongside the fitted parameters.

These advanced declarations are specified by
[ALG-FIT](../theory/algorithms.md#alg-fit) and the
[estimation reference](../api/native/parameter-estimation.rst). The propane
notebook demonstrates ordinary least squares; it does not qualify hydration
holds, minimax or measurement reconciliation for arbitrary application data.

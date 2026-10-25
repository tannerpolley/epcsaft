# Solved-state derivatives

Differentiate a named observable on a solved, regular branch. Compile the
problem once, obtain an available central solution for that same declaration,
then supply physical directions. The
[water/isobutanol action notebook](../examples/water-isobutanol-solved-state-actions.ipynb)
shows exact first and mixed actions and an independent finite-difference
verification. The [action algorithm](../theory/algorithms.md#alg-actions)
defines the implicit chains and availability conditions.

## Prepare one central problem

`compile_problem(model, problem)` prepares the native equations and keeps the
model alive. Use the corresponding `solve_equilibrium` result as the central
solution for `solved_state_actions`. Do not combine coordinates or a solution
from a different feed, support, source reference or parameter set.

A central candidate must satisfy the defining equations closely enough for
the requested action. Even an available central value can have unavailable
derivatives because the local defining Jacobian is singular or the requested
observable lacks its required energy order. See the
[equilibrium calls](../api/equilibrium.rst) and
[native action records](../api/native/equilibrium.rst).

## Declare physical directions

`SolvedStateActionDirection` supplies temperature increments in K, pressure
increments in Pa, feed-amount increments in mol and active-parameter increments
in their native units. Caloric control increments use J or J/K where admitted.
Parameter vectors follow the model's active-coordinate order; feed vectors
follow component order. A direction is not an arbitrary normalized optimizer
coordinate.

Choose `SolvedStateObservable` entries for quantities such as phase amount,
phase composition, log activity, phase property or system enthalpy. Distinguish
phase fraction from molar amount, and fixed-composition Cp from equilibrium
isobaric heat capacity. The latter includes admitted reaction and phase
responses and has additional thermal/reference requirements.

## Actions retain the observable unit

Each path variable is dimensionless. For a temperature increment ΔT,
an enthalpy first action is ΔT times the partial derivative with respect to T,
so the action is J/mol. Divide by ΔT in K to compare with Cp in J/(mol K).
The [methane caloric notebook](../examples/methane-caloric-properties.ipynb)
performs this normalization explicitly.

For dimensionless ln activity, both first and mixed actions are dimensionless.
With a temperature increment and an association-energy-over-k increment,
both increments are in K. Divide by their product to label the mixed partial
in K⁻². A unit numerical increment still carries its physical unit.

The native general action request forms the mixed contraction across its
declared directions. `solved_state_first_actions` instead returns separate
first-order batches for the supplied directions. Separate first-order results
do not constitute a mixed derivative. Available orders depend on the
observable and smooth branch; no order is implied merely by central success.

## Inspect availability and conditioning

Read each result's `value`, `action`, `unit`, `order` and named cause. Keep the
batch diagnostics: equation dimension, rank, condition estimate, central
residual/floor and linear-solve residual. A missing action must not be replaced
with zero. One unavailable derivative need not invalidate a central value or
another available action.

Support changes, ionic insertion, coalesced phases, a singular pressure root,
association singularity, missing thermal records and source-interval boundaries
can prevent the requested action. Do not differentiate through a changed
phase count by reusing the previous branch's derivative. The
[limits](../limits.md) state these restrictions.

## Verify without replacing the native calculation

The action notebook re-solves at two temperature step sizes solely to check
the exact native first action. Those finite differences are numerical
verification, not the derivative supplied to fitting or process calculations.
Agreement at that smooth branch does not establish accuracy across a phase
transition, source extrapolation or experimental validation.

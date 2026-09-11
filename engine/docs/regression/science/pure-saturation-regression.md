# Regression formulation

Required general fitting mathematics; pure saturation is one observation
construction. Installed availability is recorded in the
[transfer plan](../../../../docs/scientific/implementation-plan.md).

Use the plan's adoption gate to adapt justified existing native fitting and
evaluator logic. Implement only observation forms consumed by current scope;
the broader algebra below defines meaning rather than requiring unused nodes.

An observation combines a state/equilibrium specification, an expression, a
reported value and an explicit discrepancy/uncertainty model. Expressions use
four node kinds: Leaf (property, coordinate or constant), weighted Sum,
Product with declared exponents, and Log. This algebra covers partial pressure,
finite-dose heat, extraction fractions and ratios of selectivities without
application-specific property code. Every operation must respect its domain:
positive log arguments, nonzero denominators and valid real powers.

Fixed-coordinate leaves are valid constants with zero derivative. Observable
and parameter order, units, bases and chain rules must remain explicit through
physical and optimizer coordinates. The common state/equilibrium owner supplies
values and required actions; Regression owns no duplicate EOS or equilibrium.
The [action table](../../science/algorithms.md#alg-actions) determines required
energy order; first residual Jacobians can require third energy derivatives
when the observed quantity is a heat capacity.

Fit the declared discrepancy under stated bounds, transforms, starts, loss and
numerical controls. Retain raw predictions, physical Jacobians, weighted
optimizer-coordinate Jacobians, failed row identities, conditioning and
termination. A failed or unsupported state must not become a successful zero
residual or silently disappear from the objective.

Covariance interpretation depends on the data model. With known absolute
independent Gaussian uncertainty and a regular locally linear fit, the inverse
weighted information matrix estimates parameter covariance. With an unknown
common variance estimated from residuals, a separately estimated variance
factor is required. Correlated errors, constraints, rank deficiency, nonlinear
losses and model discrepancy require their stated treatment; no universal
variance multiplier or unconditional covariance promise applies.

A manufactured identifiable case can verify parameter recovery. Agreement with
published property data does not imply unique recovery of a published physical
parameter tuple. Record identifiability, fitted-domain and holdout evidence;
parameter adoption remains an application decision. Use the retained
[literature reproduction classes](literature-reproduction-contract.md), including
the distinction between measured values, evaluated references and source model
curves. Historical Ceres callable names and stage line/time budgets are not
current interfaces or universal acceptance criteria.

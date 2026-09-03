# Regression — a shell around Ceres, with a four-node expression

Stage 6's design document. Written 2026-09-03 from the planner's layering draft after the
independent review and a reading of the three downstream applications; the owner's
decisions are recorded in `EQUILIBRIUM.md` §12–§13 and `EOS.md` "Properties", and this
document does not restate them. It depends on `CORE.md` §3.5 (the fit floor) and §7 (the
`diagnostics` block), `EQUILIBRIUM.md` §11–§12 (the `Problem` and the sensitivity), and
`SURFACE.md` §2.5 (the two declared documents).

## 1. What regression is, in one rule

`fit(model, observations)` is a shell around Ceres with no property-specific code. An
observation is data plus a selector plus a target; the selector is an expression over the
solved or evaluated state; the Jacobian is the chain rule per node composed with the
equilibrium sensitivity. No property name is known to `regression/`: adding a field to
`State` adds it to regression with no change.

    Observation { target: Problem | StateSpec; select: Expression; value; uncertainty }
    Expression =
        Leaf(Property(phase, field) | Coordinate(lever or amount n_{p,i}) | Constant)
      | Sum(w_j · e_j)
      | Product(e_j ^ p_j)
      | Log(e)

Why four nodes and not a weighted sum plus a ratio: the applications fit `pCO2 = y_CO2 · P`
scored as `log10(pred/obs)`, `net_na = na_feed − na_aq`, an extraction fraction `li_org /
(li_org + li_aq)`, and a selectivity that is a ratio of ratios. Four nodes cover all of them
and finite-dose heat; the tree is smaller than a combinator list plus a residual-form flag.

## 2. Files

| file | lines | owns |
|---|---|---|
| `native/regression/least_squares.{hpp,cpp}` | ~490 | the Ceres kernel: one `ResidualBlock` per observation, `PredictionCache` holding the last result per observation and per Ceres iterate, covariance, the fit floor `κ_sys · residual` per observation from the solve's own `Convergence` record |
| `native/regression/expression.{hpp,cpp}` | ~80 | the four nodes, value and Jacobian; the one bound evaluator a caller may use outside a fit |
| `src/epcsaft/regression.py` | ~200 | mapping only: observations in, `Fit` out |

Rule 3 lets `regression/` include `equilibrium/problem.hpp`; rule 4 forbids it from
reaching `phase_stability.hpp`, `find_phases`, or `Discover`, so a fit can only move
parameters where the declared physics still holds.

## 3. Rules

- `StateSpec` declares its `DensityBranch`; `Unspecified` is refused inside a fit.
- `Coordinate(T)` or `Coordinate(P)` on a target where that lever is fixed is rejected by
  name at observation compile, never a silent zero Jacobian.
- No `unit` field natively; the boundary is SI (`SURFACE.md` §3.4). Weights such as
  `1/(n₂ − n₁)` are caller data or computed natively; Python does none of that arithmetic.
- A solve that fails the declared-topology contract, or leaves the domain, returns a
  failure with no state; the cost functor returns `false`, Ceres shrinks the trust region,
  and nothing is thrown across the fit or reported as a zero residual.
- No third derivatives, anywhere. Ceres never forms the residual Hessian; the implicit
  sensitivity needs second derivatives in `(y, θ)`, which the tape has.
- The active parameters are the six EOS families plus the `K(T)` correlation coefficients
  per reaction (`EQUILIBRIUM.md` §12); `∂g°/∂θ` comes from the model's
  `standard_gibbs(T, order)`.

## 4. Gate

From `TRANSFER.md` Stage 6, unchanged: NIST parameter recovery on the
`light_alkane_saturation` rows; the Jacobian against centred differences per node kind;
permutation invariance; the three-parameter saturation fit in **≤ 400 ms** (today 11,705)
at **≤ 200 recordings** (today 140,562). The saturation continuation from Stage 5 is what
makes it affordable.

Checks to name, each with its mutation in the brief: Jacobian versus centred difference per
node kind (`Property`, `Coordinate`, `Sum`, `Product`, `Log`); rule 4 with `Discover`; a
fixed-lever `Coordinate` rejected by name; a failed residual shrinking the trust region
rather than scoring zero; the A1 falsifiers for the sensitivity (`EOS.md` Properties §5).

### 4.1 Named checks

Each is a function `check_r1` … `check_r8` in `tests/native/regression_checks.cpp`, called
from `main`, non-aborting, every solve under `core::Capture`; `tools/check-equilibrium-claims.py`
reads this table as it reads `EQUILIBRIUM.md` §4.

| # | claim | asserts | mutation (must go red) |
|---|---|---|---|
| R1 | Jacobian per node kind | value and `∂/∂θ` of a `Property` leaf on a `StateSpec`, a `Coordinate` leaf on a solved `Problem`, a `Sum`, a `Product`, and a `Log`, each against a centred difference in every active parameter, to 1e-6 relative | `chain(outer, inner)` returns `outer` alone |
| R2 | several targets | an observation over two solves (a `Sum` of the same residual field of two solved liquids at two feeds, minus a constant) evaluates, and its Jacobian matches the centred difference through both solves | the evaluator visits only the first target |
| R3 | fixed lever rejected by name | `Coordinate(P)` on a fixed-`P` target fails at observation compile with `fixed_lever` in the message; no fit runs | the name changes |
| R4 | failed residual | an observation whose declared `{V,L}` solve fails (the ternary at 30 bar) inside a fit is reported failed by name, the functor returned `false`, and the fit's cost carries no zero residual for it | the failure path returns a zero residual as success |
| R5 | manufactured recovery | synthetic saturation rows from ethane's three parameters `θ*`; a fit from a 20 % perturbed start recovers `θ*` to 1e-6 relative with `Stop::converged` | `Coordinate(P)` reads the temperature column |
| R6 | the NIST rows | the ten `light_alkane_saturation` rows, pressure and liquid density as `Coordinate(P)` and `Property(L, molar_density)`, fit for the three ethane parameters: stop `converged`, residual `≤ max(criterion, floor)` from the record; AARD in `P` and `ρ`, the deviation from the published tuple, recordings and wall time **reported** | shares R5's |
| R7 | permutation invariance | relabelling the two components of a binary observation set permutes the fitted parameters and leaves the final cost unchanged to 1e-10 | the parameter offset ignores the index |
| R8 | covariance | the returned covariance equals `σ² (JᵀJ)⁻¹` recomputed from the returned Jacobian, symmetric, positive semidefinite | `σ²` replaced by one |

Two Stage 5 deviations, folded in with their own mutations: **S6-skip**, a failed existence
certificate reports its dependent cells skipped by name (harness function in
`equilibrium_support.hpp`); **S9-zero**, on the reactive fixture the feed sensitivity along a
reaction direction (`B d = 0`) is exactly zero.

## 5. What is deliberately not here

An observable enum. A residual-form flag. Any extensivity table (`Σ_p N_p P_p` is the
caller's error to make). Unit handling. A Python loop over observations: the path from the
one Python call to the one Python result is native throughout.

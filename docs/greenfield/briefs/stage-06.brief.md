# Stage 06 — Regression

Branch `stage-06`, based on `greenfield@be3f329` (the coordinator re-pins if canon moves
before the worktree is cut). Machine-readable claims are in `briefs/stage-06.brief.json`; that
file governs where the two disagree.

**`REGRESSION.md` is the specification, not background**, together with `EQUILIBRIUM.md` §12
(the observation, the expression, the levers, the rules) and `EOS.md` "Properties" §4 (the
order flag on `State`). Where this brief and those differ, they win, and that is a
`source_contradicted` stop, not a judgement call.

## What this stage is

`fit(model, observations)`: a shell around Ceres with no property-specific code. An
observation is data plus a selector plus one or more targets; the selector is a four-node
expression over the solved or evaluated states; the Jacobian is the chain rule per node
composed with the equilibrium sensitivity. No property name is known to `regression/`.

    Observation { targets: [Problem | StateSpec]; select: Expression; value; uncertainty }
    Expression =
        Leaf(Property(target, phase, field) | Coordinate(target, lever or amount n_{p,i}) | Constant)
      | Sum(w_j · e_j)
      | Product(e_j ^ p_j)
      | Log(e)

Every leaf names its target by index, so an observation over two solves (finite-dose heat is
`Sum` of a property of solve 2 minus the same of solve 1 and a constant) is one observation
with one residual. That is the multi-target form `HANDOFF.md` asked for.

**Scope decision, stated so it can be overruled.** This stage writes the regression shell
and the `PropertyOrder` derivative block on `State` that a `Property` leaf's Jacobian needs.
It does **not** write `Thermochemistry`, `physics/ideal.hpp`, the `State` totals,
`standard_gibbs(T, order)` or the seventh active family: those are EOS-category writes that
the NIST gate never touches, and they form a Properties stage that follows this one. Until
that stage lands, a `Property` leaf may name any *residual* field of `State` and the
active parameters are the six EOS families. Nothing in this stage may pre-empt those
designs by writing a partial version of them.

## The order it is written in

1. **`PropertyOrder` on `State`** (`model/helmholtz_model.{hpp,cpp}`): `state(T, ρ, x, order)`
   and `state_at_pressure(T, P, x, branch, anchor, order)`. At `FirstDerivatives` every
   first-derivative field (`pressure`, `compressibility_factor`,
   `residual_chemical_potential_over_rt`, `log_fugacity_coefficient`, `residual_entropy`,
   `residual_internal_energy`, `residual_enthalpy`, `molar_density`) carries `∂/∂T`, `∂/∂ρ`,
   `∂/∂x_i`, `∂/∂θ_k` at fixed `(T, ρ, x)` from the tape, and `∂/∂T|_P`, `∂/∂x_i|_P`,
   `∂/∂θ_k|_P` through the density root. The three fields that are already second
   derivatives of `A` carry none and say so by name (`MissingDeclaration`). No third
   derivative anywhere; the block is assembled from the Hessian `state()` already computes.
2. **`expression.{hpp,cpp}`**: the four nodes, value and Jacobian in the fit's parameter
   vector `θ`; one shared `chain(outer, inner)` helper every non-leaf node uses; leaves take
   their target's `EquilibriumResult` plus `dy/dθ`, or a `State` at `FirstDerivatives`; a
   `LeverIndex` enum `{Temperature = 0, Pressure = 1, Feed = 2}` matching
   `lever_jacobian`'s columns. One bound evaluator (value plus Jacobian) for a caller outside
   a fit.
3. **`least_squares.{hpp,cpp}`**: `fit(model, observations)`; one Ceres `CostFunction` per
   observation with the analytic Jacobian; `PredictionCache` holding the last
   `EquilibriumResult` per observation as the next solve's `start`; a failed solve returns
   `false` from the functor (the trust region shrinks, nothing scores zero, the observation
   is reported failed by name); covariance `σ² (JᵀJ)⁻¹` at the solution; the fit's
   `Convergence` record (`CORE.md` §7 row `fit`): residual is the gradient norm at
   termination against `FIT_GRADIENT_TOLERANCE`, floor is `Σ_obs ‖∂obs/∂u‖ κ_sys ·
   residual_obs` from each observation's own record, amplification is the condition
   estimate of `JᵀJ` at the solution (the number the covariance divides by), iterations
   against `FIT_MAX_ITERATIONS`. Every bound function opens one `Scope`.
4. Bindings, `regression.py` (mapping only, no arithmetic), `_core.pyi`, the wheel verifier's
   entry-point names, `meson.build`.
5. The checks, named below, and the two Stage 5 deviations.

## Traps

- **No property name in `regression/`.** A selector names a `State` field by string and
  the evaluator reads it through one accessor; a `switch` over field names in
  `regression/` is the observable enum the design removed.
- **No extensivity table.** `Σ_p N_p · field` is what the caller writes as a `Sum`; the
  shell offers no "total" helper.
- **`StateSpec` declares its branch.** `DensityBranch::Unspecified` inside a fit is refused
  by name, as undeclared topology is.
- **A fixed lever is not a zero.** `Coordinate(T)` on a target with `T` fixed is rejected
  by name at observation compile, never a silent zero column.
- **Weights are data.** `1/(n₂ − n₁)` and the like arrive in the observation; Python does
  none of that arithmetic (`python_does_no_arithmetic` covers `regression.py`).
- **Counts and times are reported, never asserted** (`TESTING`, rule 8). The saturation
  fit's recordings and wall time are printed by R6 against `TRANSFER.md`'s targets (≤ 200
  recordings, ≤ 400 ms) and read at acceptance; the check asserts convergence and the
  record, not the numbers.
- **No golden parameters.** R5 recovers the parameters that generated its own synthetic
  rows; R6 fits the NIST rows and *reports* the deviation from the published tuple.
- **Read the tree through Serena before writing.** `get_symbols_overview` on
  `equilibrium/problem.hpp`, `equilibrium/gibbs_minimizer.hpp`, `model/helmholtz_model.hpp`;
  `find_symbol` for `solve_equilibrium`, `sensitivity`, `lever_jacobian`,
  `state_at_pressure`, `ActiveParameter`; `find_referencing_symbols` for `sensitivity` and
  `lever_jacobian`. Record in the report what it found. Raw reads follow the overview,
  never replace it.

## What you inherit and may not change

- `equilibrium/`: `solve_equilibrium(model, problem, start)`, `compile_problem`,
  `sensitivity(compiled, point, lever)` with `lever` a weight vector over
  `[T, P, feed_1..feed_C]`, `lever_jacobian`. The shell calls them; it writes no row, no
  solver, no start.
- `CORE.md` §5.2: a solver never names a residual stop; call `classify`. Rule 14.
- Rule 3: `regression/` includes `equilibrium/problem.hpp` and `model/`; rule 4: it never
  reaches `phase_stability.hpp`, `find_phases` or `Discover`.
- Rule 9: a Ceres loop reads `FIT_MAX_ITERATIONS` and therefore carries `core::Scope` and
  `Convergence`. Rule 13: every bound function opens exactly one `Scope`.
- `nlp_scaling_method = none`, the tolerances and `classify` are not touched.

## Budgets

Ceilings only (owner's ruling, 2026-09-03, `CANON.md` "Budget"): a ceiling bounds the
stage, a minimum proves nothing, and completion is the named-claims gate, the mutations and
the reading. A breach is `budget_ceiling`, a stop with a number, never compaction: the file
is measured as first written.

## The checks are the stage, and they are named

Each is a function `check_r1` … `check_r8` in `tests/native/regression_checks.cpp`, called
from `main`, non-aborting (the harness of `equilibrium_support.hpp`), every solve under
`core::Capture`, a failing check printing its trace. `tools/check-equilibrium-claims.py` is
extended to read `REGRESSION.md`'s claim table too; a claim with no function is a missing
check, a function never called is worse.

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

## Stopping

`budget_ceiling`, `gate_impossible`, `source_contradicted`, `environment_broken`,
`mutation_invalid`. Stop and report rather than adapting if: a `must_contain` marker is
absent; a write would pass its ceiling; a gate green at base goes red for a reason other
than the write it names; the work seems to require a property name in `regression/`, an
extensivity helper, a second sensitivity, a third derivative, or a change to a tolerance or
to `classify`; or `REGRESSION.md` does not decide something the code must decide — report
the node or the rule and stop under `gate_impossible`. Do not accept your own work: the
coordinator rebuilds from a clean directory, runs every gate, applies every mutation, and
reads the checks.

Model and effort are the owner's per attempt; never downgrade.

# Coupled-solver tolerances: the floor is measured, the factor is a hypothesis

How the greenfield engine keeps the tolerances of nested solvers mutually consistent, so
that no stopping test demands more accuracy than its inputs can carry. Companion to
`numerical_policy.hpp` (the nodes of the tolerance graph) and checker rule 6 (no literal
outside it). This document adds the edges and the runtime check on them.

Evidence classes used throughout: **measured** (read from a file or a run), **derived**
(follows from the equations or the code by an argument given here), **guessed** (a number
or claim with no evidence yet; each one is marked).

Revised 2026-09-02 against `SPECIFICATION-SPACE.md` and the twelve-stage order: the
boundary driver is a Newton on the compiled rows, not Ipopt (§3.4); `Convergence` and
`Stop` live in the policy header, not in `equilibrium/` (§3.8); `derived()` needs
`constexpr` entries (§4); the mechanism is stated for the neutral and the electrolyte model
at once (§3.9). Nothing in §2–§6 is derived from the superseded tree; that tree appears
only in §1 as the measured failure and in Appendix A as the place the failure still runs.

## 1. The failure, in one sentence

In the tree this replaces, an outer Newton on the bubble pressure tests
`|residual| <= 2.0e-10` (`engine/cpp/equilibrium/grepe.cpp:2303`), where the residual is
`mu_vapour,ref - mu_liquid,ref` built from two separately converged inner solves: the
reactive liquid through Ipopt at `tol = 1.0e-10` (`reactive_finite_phase.cpp:27`,
`nlp_scaling_method none`) and the incipient vapour through its own Newton at
`final_norm <= 2.0e-10` (`grepe.cpp:1815`). The outer residual therefore carries at least
the sum of the two inner errors before any amplification, and the outer test sits below
that sum by construction. *Measured:* 31 of 161 rows of the MEA workload fail
(`MEA-Thermodynamics/analyses/mea_parameter_bundle/results/figure-calculation-failures.csv`);
28 carry `bubble Newton iteration limit reached`, 27 of them at 120 °C, one at 60 °C; the
other three, all at 40 °C, are `reactive_finite_phase_certificate_failed` and are a
different failure. *Measured:* the failing path builds a bare `GrepeFailure` with no
residual (`grepe.cpp:2366`), so whether the 18 closures plateaued at 3e-10 or at 3e-3 is
unrecorded.

Two corrections already established and not repeated here: Ipopt's `acceptable_tol` never
feeds the outer loop (any non-`Solve_Succeeded` status aborts at
`reactive_finite_phase.cpp:3318`), and loosening all nine `2.0e-10` sites at once proved
nothing because it moved the inner vapour Newton along with the outer test.

## 2. The principle

Write the outer residual as a function of the outer unknowns `p` and the exact inner
solutions `u_i*(p)`:

    r_out(p) = R(p, u_1*(p), ..., u_m*(p))

What is computed is `R(p, ũ_1, ..., ũ_m)` with inner residuals `‖F_i(ũ_i)‖ <= τ_i`. The
standard backward-to-forward bound gives `‖ũ_i - u_i*‖ <= κ_i τ_i`, with `κ_i = ‖J_i^{-1}‖`
evaluated at the inner solution: the conditioning of the inner system *at the point being
solved*. To first order,

    |r̃_out - r_out|  <=  Σ_i ‖∂R/∂u_i‖ · κ_i · τ_i  +  ε_eval  =:  φ      (the floor)

where `ε_eval` is the roundoff of evaluating `R` itself. Three consequences, each checkable:

1. **An outer tolerance `τ_out < φ` is unreachable by construction.** Because each inner
   solve stops on a tolerance, its error is a piecewise-constant-then-jumping function of
   `p` (the iteration count changes as `p` moves), so `r̃_out` is `r_out` plus jumps of
   size up to `φ`. A Newton iteration on such a function makes progress down to scale `φ`
   and then wanders. This is why the old tree's eighteen closures can all succeed and
   the outer test still never pass: not an inner failure and not a wrong seed.

2. **The inners sum and each is amplified through its own map.** "Outer looser than inner"
   is the special case `m = 1`, `‖∂R/∂u‖ = 1`, `κ = 1`. It holds when the inner and outer
   tests are the *same norm on the same variables*, which is why SUNDIALS can fix the
   ratios 0.1 and 0.05 between its integrator, nonlinear and linear tolerances statically
   (§4). Across a map with conditioning it does not hold, and that is the case here.

3. **The operative quantity is the achieved residual, tested against `max(τ_out, c·φ)`.**
   The declared `τ_out` is a wish; `φ` is what the inputs allow at this point; the achieved
   residual is what happened. A result that reports all three has said everything a caller
   needs. A result that reports only "iteration limit reached" has said nothing.

**Knowable a priori:** the `τ_i` (declared in the policy header), which inner feeds which
outer (`∂R/∂u_i ≠ 0` is a structural fact), and the *nominal* floor `Σ τ_i` with all
amplifications set to one. That nominal floor is the most a static check can compare
against, and it already catches the old tree (`3e-10 > 2e-10`).

**Knowable only at runtime:** `κ_i` and `‖∂R/∂u_i‖`, both point-dependent by orders of
magnitude across a sweep (an ionic liquid at 120 °C is not the same problem as at 40 °C,
and the old tree's failures cluster at 120 °C), and `ε_eval`.

Two kinds of inner solve, which the brief's framing did not separate:

- **Re-evaluated every outer iteration** (association sites inside every Helmholtz
  evaluation; the density root inside a pressure-specified state). Their error is *noise*
  on the outer residual and sets the floor of the outer stopping test.
- **Solved once, before the outer problem is posed** (the infinite-dilution reference
  that fixes `g°`). Their error is a *bias* on the problem, not on its solution process:
  the outer solve converges exactly to the solution of a slightly wrong problem. A bias
  does not enter the stopping floor. It enters the *accuracy claim*, and is reported as
  such (§3.6).

Only the first kind is a tolerance-consistency problem. Mixing the two is how a policy
header ends up with a reference threshold of `5e-5` beside a minimizer tolerance of
`1e-10` and no one able to say whether that is a contradiction. It is not: the `5e-5` is a
bias on `ln φ^∞`, which caps the accuracy of any activity claim at about that level
whatever the minimizer achieves, and an oracle tolerance on activities tighter than it
would test noise.

## 3. The mechanism

Every solver returns, beside its answer, one record:

```cpp
struct Convergence final {
    double residual;        // achieved, in the test's own norm
    double tolerance;       // the policy entry it was tested against
    double floor;           // Σ amplified inner tolerances + roundoff, at the solution
    double amplification;   // the measured κ that produced `floor`
    std::size_t iterations;
    Stop stop;              // §5
};
```

`residual / tolerance` says whether the wish was met; `residual / floor` says whether the
inputs could have allowed more; `floor / (Σ τ_i)` is the **measured amplification**, and is
what makes the static declaration of §4 auditable. Three doubles per solve, no new
allocation, and every number in the record is one the solver already had in hand.

### 3.1 Association Newton (physics, exists)

Residual `‖H‖∞ <= ASSOCIATION_RESIDUAL` (`2e-13`, dimensionless O(1)). The solver already
factorises the Jacobian with full pivoting and reads the pivot ratio
(`association.hpp`, `solve_linear`); the pivot ratio is an estimate of `κ_X` within a
small factor for `S <= 10`. Reported as `amplification`. Cost: zero, already computed.
Floor: roundoff only, `ε_mach · S`. This solve is the leaf; it has no inner.

### 3.2 The tape's gradient floor (model, one line in the Hessian path)

The minimizer's dual infeasibility is built from `Φ_y`, and `Φ_y` is evaluated at the
association solution `X̃` with residual `H(X̃)`. By the envelope property the *value* error
is second order in `δX` and can be ignored; the *gradient* error is first order:

    δ(Φ_y) = Φ_yX · δX,      δX ≈ -Φ_XX^{-1} · F(X̃),      F_a = -(n_tot m_a / X_a) · H_a

so the gradient's floor at this point is

    floor_grad = ‖ Φ_yX · Φ_XX^{-1} · diag(n_tot m_a / X_a) ‖∞ · ASSOCIATION_RESIDUAL

`Φ_yX` is `cross` and `Φ_XX^{-1}` is `factorization` in
`helmholtz_model.cpp:419-451`; the product `factorization.solve(cross)` is already formed
for the Schur complement. The equilibrium assembly evaluates every phase's tape at unit
amount (`SPECIFICATION-SPACE.md` §2.2), so there `n_tot = 1` and the floor is intensive,
one per phase per iteration; a pinned phase's floor is as real as a free one's, because
its chemical-potential rows are evaluated at the same point. One more row-sum norm of a `w × S` matrix, `O(wS)` flops against
`1 + 2w` tape sweeps: unmeasurable. Add `double gradient_floor` to `HelmholtzDerivatives`,
set on `Hessian` and `PressureSlope` orders (for the latter only the `V` row of `cross` is
in hand, which is exactly the row the density root needs).

This is the coordinator's reading made concrete: the amplification is not a number
someone declares, it is the Schur factor the tape computes anyway.

### 3.3 Density root (model, exists)

Residual `|r| <= DENSITY_ROOT_RESIDUAL` (`1e-12`, reduced units). The floor on `r` from the
association noise is `V · floor_grad[V]` (since `Z = 1 - V Φ_V`), one multiplication. The
amplification of `r` into the answer is `|dr/dρ|^{-1}`, which the solver already computes
as `derivative`; near a spinodal it diverges, and that is the existing mechanical-stability
rejection restated as a floor. Report `residual, tolerance, floor, amplification =
1/|dr/dρ|`. If `floor > tolerance`, the root test cannot be met and the stop is
`converged_to_floor` (§5), not "did not converge".

### 3.4 The two equilibrium drivers on one assembly (Stages 7, 8, 9)

`SPECIFICATION-SPACE.md` §2.4 fixes the shape: one row assembly `F(y; spec)` with its
exact Jacobian `F_y` (Stage 6, `request_compile.cpp`), consumed by a damped Newton for
`k ≥ 1` (saturation at Stage 7, every boundary at Stage 9) and by Ipopt for `k = 0` (the
flash, Stage 8). The floor is a property of the *assembly*, computed once per iteration
and read by whichever driver is running:

- **chemical-potential rows** of phase `p`: that phase's `gradient_floor` (§3.2), at unit
  amount;
- **pressure-equality rows**: `V_p · floor_grad[V]` of that phase, as in §3.3;
- **material-balance, charge and normalisation rows**: linear in the unknowns; roundoff
  only, `ε_mach · Σ_p N_p` or `ε_mach · C`.

The assembly returns the floor as one number, the maximum over rows, beside `F` and
`F_y`. Per-row floors cost nothing more but S1 tests the ∞-norm, so one number is what a
result carries.

**The boundary Newton (`k ≥ 1`).** `residual = ‖F‖∞`, `tolerance` the policy entry,
`floor` from the assembly, `amplification = κ(F_y)` from the LU the step already takes
(the pivot ratio, as in §3.1; `implicit_sensitivity.cpp` inverts the same matrix and can
report the exact smallest singular value at the solution). Classification is §5 directly:
this driver owns its loop, so the stall window and the budget are its own counters and no
post-hoc reading of a foreign status is needed. This is the driver the first equilibrium
stage runs, and the regression fit's inner solve (§3.5), so it is the one the tolerance
core must serve before Ipopt exists.

**Ipopt (`k = 0`).** Ipopt tests the scaled NLP error against `tol` and the unscaled parts
against `dual_inf_tol`, `constr_viol_tol`, `compl_inf_tol`
(https://coin-or.github.io/Ipopt/OPTIONS.html). Set `nlp_scaling_method = none`, because
the floor is computed in the assembly's own units and a scaled test would need the same
factors applied to the floor; if scaling is ever enabled, divide the floor by them.

The floor on dual infeasibility is the running maximum of `gradient_floor` over the
evaluations Ipopt made (one `std::max` per evaluation, inside the one-state-per-iteration
cache Stage 8 already plans). The floor on primal infeasibility is roundoff on the linear
conservation and charge rows, `ε_mach · C · n_tot`, which is why a gate of `1e-10` on those
rows is safe: it sits five decades above the floor.

Ipopt has no hook to move `tol` mid-run, and none is wanted. Run with the declared
tolerance. On any status other than `Solve_Succeeded`, read the final unscaled `inf_du`,
`inf_pr` (available from the `intermediate_callback` and `finalize_solution`), and
classify against the floor (§5). On `Solve_Succeeded`, still report the floor: a solve
that met `1e-10` with a floor of `9e-11` is one whose next point may not, and a sweep
should see it coming.

At the solution, `implicit_sensitivity.cpp` factorises `F_y` (the KKT matrix for `k = 0`,
the boundary Jacobian for `k ≥ 1`) anyway; its smallest singular value gives `κ_sys`, the
amplification of the achieved residual into the *solution's* error, which is the floor
the next layer up needs (§3.5). The dimension is `C·Π + rows`, twenty-ish; one dense SVD
at the solution costs microseconds and runs once per solve.

### 3.5 Regression (Stage 10, Ceres)

The outer here is a least-squares fit over solved equilibria. Its inner solve is
whichever driver the observation's specification compiles to: the boundary Newton for
every saturation and bubble-point observation, Ipopt for a flash observation. Each
residual is an observable of a solution with error `κ_sys · residual` read from that
solve's own `Convergence` record, mapped through `∂obs/∂u`, and the Jacobian comes
through the same implicit sensitivity. Ceres has no floor hook either, so the same
post-hoc rule applies: the fit's gradient norm at termination is compared with the sum
over observations of `‖∂obs/∂u‖ · κ_sys · residual`, and `function_tolerance` /
`gradient_tolerance` are `derived()` from the driver's convergence entry with that
structure (§4). Because a pure fit and a mixture fit share one assembly and one policy
entry (`TRANSFER.md` §"Saturation is equilibrium, fully"), they share one floor, and a
parameter set that converges cleanly at `C = 1` and sits at its floor at `C = 2` says so
in the fit result rather than in a silent difference between two solvers.
Ceres's defaults (`1e-6`, `1e-10`, `1e-8`) are loose relative to the minimizer's floor in
the cases measured so far, so this is expected to classify as `converged` and the value
of doing it is that the fit result *says so* with numbers rather than by silence.

### 3.6 Infinite-dilution reference (model, exists)

A once-solved bias, not a floor (§2). Its ladder convergence
(`REFERENCE_CONVERGENCE_TOLERANCE`, `5e-5` on `Δ ln φ`) is reported on the reference
result as `bias`, and any consumer's accuracy claim on activities or reaction constants
is capped at that value. It does not enter any stopping test.

### 3.7 What this costs

*Derived:* per Helmholtz Hessian one `w × S` row-sum norm; per density root two
multiplications; per assembly one `max` over rows; per driver iteration one `max`; per
solve one small SVD that Stage 7 needs regardless. No new evaluations, no new sweeps, no
new allocations beyond three doubles per result. *Guessed:* under one percent of any
solve's wall time; to be read off `model_floor.cpp` when the gradient floor lands.

### 3.8 Where it lives

`Convergence` and `Stop` are consumed first by `DensityRoot` in `model/` at Stage 5, and
`physics` is the only layer everything may include, so they are defined in the policy
header beside the tolerances they report on. `problem.hpp` carries one, it does not
define one.

| layer | file | stage | adds |
|---|---|---|---|
| physics | `numerical_policy.hpp` | 5 | `Convergence`, `Stop`, `derived()` (§4), `FLOOR_MARGIN`, `STALL_WINDOW` |
| physics | `association.hpp` | 5 | returns the pivot ratio it already computes |
| model | `helmholtz_model.hpp` | 5 | `HelmholtzDerivatives::gradient_floor` |
| model | `density.cpp` | 5 | `DensityRoot::convergence` |
| equilibrium | `request_compile.cpp` | 6 | the assembly's floor: max over rows of the per-phase floors and the linear rows' roundoff |
| equilibrium | `phase_boundary.cpp` | 7, 9 | the Newton's `Convergence`; `κ(F_y)` from its LU; the classification |
| equilibrium | `implicit_sensitivity.cpp` | 7 | `κ_sys` from the factorisation it owns |
| equilibrium | `gibbs_minimizer.cpp` | 8 | the running max over Ipopt's evaluations; the post-hoc classification |
| regression | `least_squares.cpp` | 10 | the aggregate floor, the classification |

No structure here is carried from the superseded tree. Its results had a name–value list
of the tolerances in force and a condition number nothing read; neither is copied, and the
record above is written from §2, not from either.

### 3.9 One model, both parameterisations

There is one Helmholtz model; the neutral and the electrolyte cases are one code path with
terms switched by the parameter record (`STAGE-LOG.md`, Stage 2: one permittivity
equation, four forms). The floor depends only on which *inner solves* exist, and the
electrolyte terms add none: Debye–Hückel, Born and every permittivity form are closed
forms on the tape. So:

- **Neutral, non-associating** (`S = 0`): `gradient_floor` is zero and every floor is
  roundoff. `Convergence` still carries the four numbers; they are small.
- **Neutral, associating**: §3.2 in full.
- **Electrolyte**: §3.2 unchanged, plus three things that are not floors and must not be
  mistaken for them. The charge rows are linear (roundoff, §3.4). The finite-size series
  switch (`SCREENING_SERIES_SWITCH`) is a recording-time predicate: both forms are exact
  and the value is continuous across it to the cancellation floor, so a re-record inside
  one solve moves the residual by roundoff, not by a tolerance. The infinite-dilution
  reference is a bias (§3.6), and it is the electrolyte case's only bias.

The mechanism is therefore parameter-set-agnostic by construction, which is the same
property the checker's rule 7 demands of the solvers: nothing in it names a system kind.

## 4. The static guard, sized honestly

```cpp
// A tolerance that cannot be tighter than what feeds it. `factor` is a hypothesis
// about amplification, tested at runtime by Convergence::floor / upstream.value.
consteval Tolerance derived(const Tolerance& upstream, double factor, const char* why) {
    if (!(factor >= 1.0)) throw "derived tolerance tighter than its input";
    return {upstream.value * factor, Category::AlgorithmTolerance, why};
}
```

Two conditions on the header for this to compile, both cheap and both Stage 5's:

- **The entries become `inline constexpr`.** Today all 27 are `inline const`
  (*measured*), which is not a constant expression, and a `consteval` call on one is
  ill-formed. `Tolerance` is already a literal type. Checker rule 5 then sees every entry
  by name; the names are unique, so it stays green, and it now also refuses a second
  definition of any tolerance, which is the drift rule 5 was written for.
- **One `derived()`, two callers.** `TRANSFER.md` §"Oracle" uses the same name for the
  closeness constructor a brief builds from a parameter set's reported deviation. It is
  the same function: the brief constructs the upstream `Tolerance` from the source's
  deviation with that source as `source`, and `derived()` refuses to test the engine
  tighter than the parameter set was ever claimed to be. One name, one meaning: a
  tolerance that cannot be tighter than what feeds it.

What it catches, and did catch in the old tree's shape:

- **Inversion in the same-norm case.** The old tree's outer `2e-10` below the vapour
  Newton's `2e-10` plus the liquid's `1e-10` is caught by the nominal floor, factor one,
  with no knowledge of conditioning at all.
- **Propagation.** Twelve literals, five meanings, nine sites for one quantity: the failure
  was that no one *could* change one without changing the others, because nothing recorded
  that they were related. Rule 6 makes each tolerance one node; `derived()` makes each
  dependency one edge, greppable, and a change to an upstream value moves every downstream
  one by construction.

What it cannot catch, which is the case that actually failed:

- **A wrong factor.** Declare `10.0` where the point's conditioning is `1e4` and the check
  passes with the tolerance still unreachable. The coordinator's reading is correct: on
  its own, `derived()` moves the unjustified number from the comparison site to the
  declaration site.

Why it still earns its place: the runtime record reports `floor / Σ τ_i`, the measured
amplification at the solved point. A declared factor of `10.0` beside a measured `1e4` is
a contradiction the result can print. The static declaration is therefore not a check; it
is a **hypothesis that the runtime falsifies**, and a declaration that is never falsified
across the oracle rows is one that has been tested against real conditioning rather than
asserted. That closes the loop the coordinator's objection left open. Six lines; keep it.

The field's practice agrees with this division. SUNDIALS fixes ratios statically and is
right to, because its three levels share one weighted norm on one set of variables
(https://sundials.readthedocs.io/en/latest/cvode/Mathematics_link.html: nonlinear test
`R‖δ_m‖ < 0.1 ε`, linear residual `< 0.05 · 0.1 ε`). PETSc's Eisenstat–Walker sets the
inner tolerance from the outer residual *at runtime* (`-snes_ksp_ew`,
https://petsc.org/release/manualpages/SNES/SNESKSPSetUseEW/; Eisenstat & Walker, SIAM J.
Sci. Comput. 17 (1996) 16–32, https://doi.org/10.1137/0917003), which is the dual of the
mechanism here: they adapt `τ_in` to `r_out`, this engine fixes `τ_in` and measures `φ`.
Eisenstat–Walker is the better choice when the inner solve is expensive and its accuracy
is wasted early; the association Newton costs 3–8 steps to `2e-13` from a warm start, so
there is nothing to save, and the fixed-inner, measured-floor form is the right one here.

## 5. Named stops

`iteration limit reached` conflates four situations a caller must tell apart. The enum,
with the numbers each carries in `Convergence`:

| stop | meaning | carries | caller's move |
|---|---|---|---|
| `converged` | `residual <= tolerance` | all three | none |
| `converged_to_floor` | `tolerance < residual <= FLOOR_MARGIN · floor` | residual, tolerance, floor, amplification | accept with the degraded claim; the tolerance was unreachable at this point, not the problem unsolvable |
| `stalled` | `residual > FLOOR_MARGIN · floor` and no contraction over the last `STALL_WINDOW` steps | residual, floor, best residual, contraction ratio | a different seed or a continuation step; the basin, not the budget |
| `budget_exhausted` | still contracting at the iteration budget | residual, contraction ratio, iterations | a larger budget or a closer seed; distinct from `stalled` |
| `left_domain` | an iterate left the physical domain | which variable, its value | the seed or the bounds |
| `singular` | Jacobian or KKT rank-deficient at the iterate | pivot ratio or rank, dimension | the specification (well-posedness), not the numerics |
| `infeasible` | Ipopt restoration failed or `inf_pr` cannot fall | `inf_pr`, its floor | the specification |

Two policy entries, both category 3 and both honest about being declared:
`FLOOR_MARGIN` (a decade: the first-order bound's slack across the norm inequalities;
*guessed*, to be tightened once measured ratios exist across the oracle) and
`STALL_WINDOW` (five steps: below Ipopt's `acceptable_iter` default of 15 because a stall
at the floor shows in far fewer). Neither changes an answer. Both only choose a name.

`converged_to_floor` is a **success with a stated accuracy**, and that is the point of the
whole design: the old tree's eighteen successful closures were this outcome, misreported
as a failure with no number attached. It is not a pass on a test row: `TESTING.md` gates
compare against measured data or closed forms with their own tolerances, and a gate
constructed as `derived(...)` below the floor tests the floor rather than the equations,
which is exactly what TRANSFER.md's Stage 7 and 8 gates were rewritten to avoid.

One outcome is deliberately **not** a stop: a declared phase set that converges to fewer
phases than declared, or to two coincident phases. That iterate has met its tolerance;
what it has failed is the specification's contract (`SPECIFICATION-SPACE.md` §3.4). It is
returned as a failure with no equilibrium state, named `declared_phase_absent` or
`phases_coincide`, carrying the phase and either Michelsen's deficit `1 − Σ` or the phase
distance. Putting it in the `Stop` enum would let a caller read amounts off a result whose
phases are not the ones asked for, which is the trivial solution reported as an answer.

## 6. Which couplings survive the rebuild

| coupling | old tree | greenfield | relation |
|---|---|---|---|
| bubble pressure ↔ liquid + vapour solves | outer Newton over `P` around two inner solves; the failure of §1 | **gone**: `P` (or `T`) and the pinned phase's composition are unknowns of the one compiled row system (`k = f = 1`), solved by the boundary Newton on the whole system at once; one residual, one tolerance | none to state; the instance no longer exists |
| association sites → tape → either driver | `2e-13` under `(S+1)` Hessians, no accounting | survives, for the flash (Ipopt) and every boundary and saturation (Newton) alike | `floor = max over phases and rows` of §3.4; `‖Φ_yX Φ_XX^{-1} D‖∞ · 2e-13` per phase (§3.2) |
| association → density root | as above under the 360-interval scan | survives for `state_at_pressure`, the reference ladder, and the `C = 1` seed's root count | `floor_r = V · floor_grad[V]`; `κ = 1/|dr/dρ|` (§3.3) |
| flash solution or previous boundary point → boundary seed | 21-start grid | survives as a seed only | a seed carries no tolerance; nothing to relate |
| single-phase flash ↔ density root (S13) | separate solvers | the pressure row *is* the density equation at `Π = 1` | agreement tested within the larger of the two `max(tolerance, floor)`; a bare number here would test the floor |
| either driver → implicit sensitivity → regression | owner compilation | survives; the inner is the Newton for every saturation and boundary observation | `floor_fit = Σ_obs ‖∂obs/∂u‖ κ_sys · residual` (§3.5) |
| infinite-dilution reference → `g°` → minimizer | four fixed-P scans | survives | a bias, `5e-5` on `ln φ^∞`; reported, never a floor (§3.6) |
| stability inner minimizations → dual bound loop | HELD2 stages I–II | survives at Stage 11 | inner is §3.4; the gap test `ε_b = 1e-2` sits eight decades above any floor and needs no relation |
| Ipopt's own linear solves → its KKT test | Ipopt-internal | unchanged | same norm, same variables: Ipopt's iterative refinement owns it |

The first row is a design commitment for Stages 7 and 9 and the reason the instance
disappears. Under `SPECIFICATION-SPACE.md` §2 a pinned phase keeps its composition and
volume as unknowns and loses only its amount; `P` (or `T`) is freed to pay for the pin;
the chemical-potential and pressure rows of every phase are rows of one system; and that
system is solved *as one system* by a damped Newton, never as a minimisation and never as
an outer loop over `P` around `solve_equilibrium`. If a brief writes the outer loop, the
coupling returns with the old tree's floor, and the brief must then carry §3.4's floor
through it. Prefer the rows. teqp makes the same move for its VLE by taking density as a
variable of the outer system rather than an inner solve, and has no density-loop coupling
to report.

## 7. Prior art, surveyed

For each: where a tolerance is declared and who owns it; whether nested solvers relate
their tolerances; what scaling precedes the test; what a failure returns. All *measured*
from the cited source on 2026-09-02 unless marked.

| package | declared where | nested relation | scaling before test | failure returns | verdict |
|---|---|---|---|---|---|
| **FeOs** ([bubble_dew.rs:18-24](https://github.com/feos-org/feos/blob/main/crates/feos-core/src/phase_equilibria/bubble_dew.rs), [density_iteration.rs:159](https://github.com/feos-org/feos/blob/main/crates/feos-core/src/density_iteration.rs)) | module `const`s: `TOL_INNER 1e-9`, `TOL_OUTER 1e-10`, density `abstol 1e-12, reltol 1e-14`; `SolverOptions{tol}` overrides | none; the outer `1e-10` is *tighter* than the inner `1e-9` it consumes, the same inversion as §1 | K-value residual `Σ|K x/y - 1|`, unscaled; density test `|p - p_spec| < max(abstol, ρ·reltol)` | `FeosError::NotConverged("bubble-dew-iteration")`, a string, no residual | same bug |
| **teqp** ([iteration.hpp](https://github.com/usnistgov/teqp/blob/main/include/teqp/algorithms/iteration.hpp), [VLE.hpp](https://github.com/usnistgov/teqp/blob/main/include/teqp/algorithms/VLE.hpp)) | caller-supplied `StoppingCondition` objects; `MixVLETpFlags{atol, reltol, axtol, relxtol, maxiter}` | none needed: density is a variable of the VLE system, so there is no inner density loop | pressure rows `/p`; per-row `atol + reltol·|r|` | `VLE_return_code` enum plus `MixVLEReturn{r, initial_r, num_iter, num_fev, success}` | avoids nesting by formulation; the best failure return in the set |
| **Clapeyron.jl** ([bubble_chempot.jl:52-57](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/master/src/methods/property_solvers/multicomponent/bubble_point/bubble_chempot.jl), [volume.jl:153](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/master/src/methods/property_solvers/volume.jl)) | method struct `atol 1e-8, rtol 1e-12, max_iters 1000`; volume solver literal `1e-12` on `log V` | none | volume step scaled by `1/(V dp/dV)`; outer unscaled | `NaN`, silently | same bug, worst reporting |
| **Reaktoro** ([EquilibriumOptions.hpp](https://github.com/reaktoro/reaktoro/blob/main/Reaktoro/Equilibrium/EquilibriumOptions.hpp), [Optima Result.hpp](https://github.com/reaktoro/optima/blob/main/Optima/Result.hpp)) | `epsilon 1e-16` amount floor; `Optima::Options::convergence` | none needed: activity models are closed-form, no inner solve | Optima-internal | `Result{succeeded, failure_reason, error, error_optimality, error_feasibility, iterations}` | not the same problem; returns numbers |
| **IDAES / Pyomo** ([DiagnosticsToolbox](https://idaes-pse.readthedocs.io/en/stable/reference_guides/core/util/diagnostics/diagnostics_toolbox.html)) | Ipopt options plus toolbox thresholds (`constraint_residual_tolerance 1e-5`, `jacobian_large_value_warning 1e8`, ...) | none needed: equation-oriented, one KKT; nested solves are a modelling error the toolbox hunts | user's job; toolbox reports extreme Jacobian entries and the Jacobian condition number | Pyomo `TerminationCondition` plus Ipopt status | avoids nesting by formulation; the only one that reports κ, but a posteriori and for a human |
| **Cantera ChemEquil** ([ChemEquil.h](https://github.com/Cantera/cantera/blob/main/include/cantera/equil/ChemEquil.h), [ChemEquil.cpp](https://github.com/Cantera/cantera/blob/main/src/equil/ChemEquil.cpp)) | `EquilOpt{relTolerance 1e-8, absElemTol 1e-70, maxIterations 1000}` | none | element rows `elMoles·rel + abs` | `CanteraError("no convergence in N iterations")`, no residual | same bug |
| **thermo** (Bell) ([flash_utils.py:266-395](https://github.com/CalebBell/thermo/blob/master/thermo/flash/flash_utils.py)) | `sequential_substitution_2P(tol=1e-13, maxiter=1000, trivial_solution_tol=1e-5)` | none; the phase volume solve has its own | `Σ(K x/y - 1)²`, unscaled | `UnconvergedError("End of SS without convergence")`, no numbers | same bug |
| **phasepy** ([bubble.py:139-178](https://github.com/gustavochm/phasepy/blob/master/phasepy/equilibrium/bubble.py)) | literals `tol 1e-8, itmax 10` | none | `|f|`, unscaled | no exception; `error` in `full_output` | same bug; at least returns the number |
| **CoolProp** ([VLERoutines.cpp:647-805](https://github.com/CoolProp/CoolProp/blob/master/src/Backends/Helmholtz/VLERoutines.cpp)) | literals `1e-9`, `1e-10`, `1e-11`, `50`/`100` iterations, a `1e-3` post-check on `p` | none | per-row residual norms | `SolutionError` strings, some with the current error ("current error is %Lg"; "solver was good, but went bad") | same bug; the "went bad" branch is a floor symptom handled as an exception |
| **DWSIM** ([NestedLoops.vb:321-324, 647-671](https://github.com/DanWBR/dwsim/blob/windows/DWSIM.Thermodynamics/FlashAlgorithms/NestedLoops.vb)) | `FlashSettings`: `PTFlash_External_Loop_Tolerance`, `..._Internal_Loop_Tolerance`, both iteration limits | the two loops' tolerances are *named as a pair* in user settings; no relation enforced | Rachford–Rice `|F|`, unscaled | `Exception` with T, P, z in the message | same bug; the only package that makes the nesting explicit as data |

**Plainly:** none of the ten thermodynamic packages relates a nested tolerance to what
feeds it, and none computes a floor. Three (teqp, IDAES, Reaktoro) do not have the
problem because their formulation has no inner solve, which is the greenfield's own move
for pressure (§6, row 1) and the strongest form of the fix where it is available. Two
(teqp, Reaktoro) return numbers on failure; DWSIM names the two loops; IDAES reports the
condition number. The full mechanism, a floor measured from the conditioning at the solved
point and tested against the achieved residual, appears in none of them. It does appear in
general numerics, in two halves that no thermodynamic package has put together:
Eisenstat–Walker for the runtime linkage and SUNDIALS for the same-norm static ratio.
The brief's prior was right. The static guard is a copy; the runtime floor is a
contribution.

## 8. Where the brief's framing was wrong, and where it was right

- **Right:** runtime measurement is the mechanism, the static check is hygiene. §4 goes
  one step further: the static declaration is a hypothesis the runtime tests, which is
  what makes it worth keeping at all.
- **Right:** "everyone hardcodes it" is the finding. §7.
- **Incomplete:** "the inners sum, amplified" treats every inner as noise. Once-solved
  inners are bias (§2), and the policy header already contains one (`5e-5`) that would
  otherwise look like a contradiction beside `1e-10`.
- **Understated:** the pressure coupling does not merely "disappear" as a consequence of
  the KKT formulation; it disappears only if the boundary is written as rows of the one
  compiled system and not as an outer Newton. §6 makes that a commitment rather than an
  observation.
- **Rejected:** any version of the mechanism that adapts the inner tolerance to the outer
  progress (Eisenstat–Walker). Correct in its domain, wrong here: the inner solve is
  already cheap at full accuracy, and adapting it would trade a measured floor for a
  moving one.

## Appendix A. The superseded tree: evidence, and a retention patch

This appendix is the only part of the document that reads the superseded tree, and it
reads it as the place where the measured failure of §1 still runs today. Nothing in it is
a source for §2–§6, and nothing in it is to be carried.

Values unchanged. Nothing below moves a tolerance. The order is forced: the residual is not
recorded, so no experiment can be read until it is.

### A.1 Evidence retention on the failing path (do this first)

At `grepe.cpp:2295`, the eighteen-iteration loop, and `:2365-2368`, the bare failure.
Keep the residual trajectory and the last state; the `evidence` vector on `GrepeResult`
and `incipient_liquid_partial` already exist, so this is additive:

```cpp
// grepe.cpp, inside solve_incipient_pressure, before the `for (int iteration ...)` loop:
std::vector<double> trajectory;
// after `++closure_count;`:
trajectory.push_back(current->residual);
// replace the bare failure at :2365-2368 with:
GrepeResult partial = current.has_value()
    ? incipient_liquid_partial(request, current->liquid, current->pressure_pa)
    : GrepeResult{};
partial.evidence.push_back({"bubble_log_closure_abs", current ? std::abs(current->residual) : NAN});
partial.evidence.push_back({"bubble_closure_slope", current ? current->residual_log_pressure_derivative : NAN});
partial.evidence.push_back({"bubble_liquid_kkt_condition", current ? current->liquid.solved.equilibrium.condition_number_inf : NAN});
for (std::size_t k = 0; k < trajectory.size(); ++k)
    partial.evidence.push_back({"bubble_residual_" + std::to_string(k), trajectory[k]});
failures.push_back(failure.value_or(GrepeFailure(
    "numerical_convergence_failure", "bubble Newton iteration limit reached",
    std::move(partial), start_index)));
```

Twelve lines, no behaviour change on any converging row, and the partial already flows to
the caller through `GrepeFailure::partial_result()`. Then re-run the 28 rows and read
three numbers per row: the final residual, the ratio of the last two residuals, and the
liquid's KKT condition number.

### A.2 What the numbers will say

*Derived, not yet measured:*

- Final residual in `[2e-10, ~1e-8]` with ratio near one: a plateau at the floor.
  §2 predicts this, and `converged_to_floor` is the honest classification. The rows are
  solved; the claim is degraded to the plateau value.
- Final residual `≫ 1e-8` with ratio near one: `stalled`, a seed or basin problem,
  and no tolerance change would help.
- Ratio well below one at iteration 18: `budget_exhausted`; eighteen is a literal too.
- A condition number that is decades larger at 120 °C than at 40 °C explains the
  clustering without any further hypothesis. If it is not, the clustering needs one.

### A.3 Only after A.2, and only if it plateaus

The minimal reversible change is to *one* site, `grepe.cpp:2303`, and it is not a new
tolerance: accept when `|residual| <= max(2.0e-10, margin · floor)` with `floor` the
measured plateau per row, and record `converged_to_floor` with the residual on the
result. The inner two tolerances stay where they are, which is what the earlier
all-nine experiment failed to do. `TESTING.md` is not violated because no test row is
made to pass by this; the MEA workload's rows are reported, never asserted, and the
report gains a number it did not have.

What not to do: loosen `2.0e-10` to `1e-8` at the outer site and call the 28 rows
converged. A row that converges because a tolerance moved has not converged; it has been
re-described. The floor and the residual on the result are what make the description
honest.

### A.4 Out of scope here

The three 40 °C rows fail the liquid certificate
(`reactive_finite_phase_certificate_failed`), which tests the *inner* solve against
`kChemicalKktTolerance = 1e-7` and its siblings (`chemical_equilibrium.hpp:13-23`).
Different mechanism, same shape: a certificate threshold with no relation to Ipopt's
`1e-10`. Not this document's problem, and not to be touched with it.

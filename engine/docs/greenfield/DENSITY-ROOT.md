# The density root: three defects and what they cost

Status: findings, 2026-09-02, from a measured sweep of 630 cases (5 fixtures × 7 temperatures
× 6 pressures × 3 branches) plus three real packets. Independently verified on the two claims
that carry the rest. Evidence is **measured**, **derived** or **guessed**.

## 1. The central defect: `PressureSlope` returns a frozen-site derivative

`DerivativeOrder::PressureSlope` sets

```cpp
result.volume_second_derivative = second[2 * 1 + 1];        // helmholtz_model.cpp:385
```

which is the raw tape second derivative at **frozen site fractions**, while the `Hessian`
order sets

```cpp
result.volume_second_derivative = result.hessian[1 * total_coordinates + 1];   // :460
```

which is the Schur-corrected total derivative `base(1,1) − cross_Vᵀ H_X⁻¹ cross_V`, including
`dX/dV`. **Verified by reading both branches.** Central differences on `P` confirm which is
which — associating binary at 313.15 K, ρ = 33000: finite difference 70463.61706, Hessian
path 70463.61705, `PressureSlope` path 40494.85915.

`solve_density` consumes `PressureSlope` for both its Newton derivative and its
mechanical-stability filter, so with association active it is iterating on the wrong slope.

**Consequences, measured:**

- The Newton degenerates to a linear iteration at rate `|1 − exact/frozen|`. Predicted 0.7358
  for the 313.15 K liquid root; the trace measures 0.735 per step over 60 steps. Steps to
  reach 1e-12: 90 associating, 99 aqueous, 118 nine-species — against a
  `DENSITY_MAX_ITERATIONS` budget of 60.
- **102 of 630 sweep cases throw** `density root refinement did not converge`. **Verified on a
  real packet through the public API:** `figiel-2025-reference-electrolytes/2`, NaCl at
  280 K and 5 × 10⁵ Pa, raises `EosError: density root refinement did not converge`.
- The stability filter tests a quantity up to 368× wrong. Its *sign* never disagreed across
  471 converged cases, so no returned density is currently wrong from this — luck, not a
  guarantee.
- With the exact slope the same brackets converge in 4–8 evaluations.

**Why no gate caught it.** `model_checks.cpp:227` finite-difference-checks
`State::pressure_density_derivative`, which comes from the *Hessian* path and is correct. The
only fixed-pressure native test (`:276`) uses `model_support::ethane()` — **non-associating**,
where the frozen and exact slopes are identical by construction. `PressureSlope`'s
`volume_second_derivative` has no finite-difference check anywhere. The suite is not weak; it
tests the path that works and never exercises the one that does not.

## 2. It can miss the liquid root

At **ethane, 280 K, 10⁴ Pa** the reference (1200-interval bisection) finds three roots:
4.299 (dP/dρ = +2324), 9567.8 (−461), 11163.5 (+645). The two dense roots are a factor 1.167
apart and one scan interval is 1.985× wide, so both fall inside a single interval, the two
sign changes cancel, and neither is seen. `phase="liquid"` returns the **vapour** root with
`stable_root_count = 1` and no diagnostic. 200 intervals are needed to find it.

`DENSITY_BRACKET_INTERVALS = 40` is inherited rather than tuned; the miss count against the
reference is non-monotone in granularity — 10 → 26 misses, 20 → 2, 30 → 4, 40 → 1, 200 → 0 —
and polish cost is flat above 20 intervals, so a finer scan buys root coverage and nothing
else.

## 3. The anchor path bypasses branch selection

Requesting `Liquid` with a vapour-side anchor returned the **vapour** root at **57 of 57**
two-root conditions, labelled `branch = Liquid` with `stable_root_count = 1` — a count it
never computed. Example: associating binary, 450 K, 10⁵ Pa, liquid root 29013, anchored
answer 27.05.

This is why threading more anchors through the callers (§5) must wait for the fix: more warm
starts would mean more silently wrong branches.

## 4. Cost

Per cold solve, `OMP_NUM_THREADS=1`:

| fixture | evaluations | scan | polish | wall |
|---|---|---|---|---|
| ethane | 46 | 41 | 5 | 273 µs |
| neutral binary | 46 | 41 | 5 | 392 µs |
| associating binary | 99 | 41 | 58 | 2007 µs |
| aqueous electrolyte | 90 | 41 | 50 | 1915 µs |
| nine species | 86 | 41 | 45 | 6437 µs |

**This corrects `performance_floor.md` Mechanism 4**, which attributes the cost to the scan.
On neutral fluids polishing is 11 % of evaluations; on associating fluids it is **60–66 %**.
The scan is not the dominant term on the fluids that matter — the broken Newton is.

Through Python: nine-species cold `Mixture.state` is 8024 µs at 108 density evaluations;
with an exact anchor, 460 µs at 1 evaluation. `Reference.activity(1.0)` is 4.04 ms — one
fully cold solve per activity call.

## 5. Ranked, with what each changes

| # | change | saving | changes an answer? |
|---|---|---|---|
| 1 | Give `PressureSlope` the true `dP/dρ` — the same Schur correction the `Hessian` order already forms | polish 27–79 → 4–8 evaluations; nine-species solve 7.6 → 3.2 ms | **Yes.** Turns 102 throws into answers and changes what the stability filter accepts |
| 2 | Scan at `DerivativeOrder::Gradient` — the scan reads only `residual`, which needs 2 sweeps, and asks for 4 | 1.36 ms per solve; a further 43 % after #1. Answers bit-identical on 20 of 20 rows | No |
| 3 | Make the anchor path respect `branch` and stop asserting a stable-root count it never computed | correctness only | **Yes**, explicitly |
| 4 | Thread anchors through `Pure.state`, `electrolyte_activity`, and the two cold `state_at_pressure` sites in `infinite_dilution_reference.cpp` | 2.9–4.7× end to end | No — **but only after #3** |
| 5 | Seeded two-branch Newton (`P/RT` for vapour, 0.6 × the hard-sphere bound for liquid) with the scan as fallback | 41 + 27–79 → about 9–15 evaluations; 75 of 75 at two-root conditions | **Yes** — the scan enumerates roots and can report ambiguity; two seeds cannot prove there are only two. Argue separately |

Changes 1–4 need no new policy entry and no literal. Change 5 needs one seed fraction,
expressible as a multiple of `CLOSE_PACKING_FRACTION` through the existing
`molar_density_upper_bound`; it is a seed, not a tolerance, so it does not sit under
`DENSITY_ROOT_RESIDUAL` in the nesting order.

## 6. What this says about `TOLERANCES.md` §3.3

The `Convergence` record sketched there wants `amplification = 1/|dr/dρ|` — **which is exactly
the quantity defect 1 corrupts.** It cannot be built correctly until that is fixed.

And the floor is **not** currently the binding constraint: with a correct Newton the residual
reaches 1e-16 to 5e-13 everywhere tested, and over the 471 converged cases the achieved
reduced residual is p50 1.17 × 10⁻¹³, max 9.97 × 10⁻¹³. `converged_to_floor` would never
fire here. The 102 failures are a real non-convergence, not a mislabelled floor.

`DENSITY_MAX_ITERATIONS`'s own source string — "5-8 steps from an anchor, 20-40 for a cold
bracket" — is falsified by measurement today and becomes true after change 1. The value needs
no change.

## 7. Delivered accuracy, stated

The tolerance is on the **reduced** residual, so the delivered relative pressure error is
`|r|/z`, and on a dense liquid root `z ≈ 7e-5…7e-4`. A reduced 1e-12 therefore buys only about
1e-8 relative. Measured worst cases: 1.23 × 10⁻⁸ (MEA nine-species), 1.70 × 10⁻⁸ (Figiel
NaCl), 6.6 × 10⁻¹³ (methane–ethane). Neutral fluids sit at p50 3 × 10⁻¹⁶; associating and
electrolyte fluids at p50 3.3–3.7 × 10⁻¹². This is worth stating in the policy header's
rationale, because "1e-12" reads tighter than what a caller actually receives.

## Not measured

- The association-noise contribution to the residual floor (`TOLERANCES.md` §3.3's
  `V · floor_grad[V]`): `solve_site_fractions` does not expose its converged residual, so the
  floor was bounded from above but not measured.
- The isolated wall cost of the S×S correction in change 1 — bounded below one sweep of noise
  by difference, not timed directly.
- Composition. The sweep varies temperature and pressure at one composition per fixture, so
  the two-in-98 real-packet failure rate is a **lower bound**.

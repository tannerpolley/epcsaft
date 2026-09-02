# The specification space: how the equilibrium routes prove they are general

Status: design, 2026-09-02. Stage briefs 6, 7, 8, 9 and 11 are written against it. Nothing
here is code yet; every line estimate is a guess until `STAGE-LOG.md` records the actual.

Evidence is tagged **measured** (read off the tree or run), **derived** (follows from the
equations in `formulation.md`), or **guessed** (a prediction).

## 1. The claim, and why rows cannot test it

Sixty-three oracle rows prove that the routes those rows take work. The owner's diagnosis
is exact: a test that says a method call works for one application never says the route is
general. The old tree is the proof. `engine/cpp/equilibrium/grepe.cpp:3540-3580` on `main`
dispatches to six solvers on `request.temperature.solved`, `request.boundary.has_value()`,
`incipient_count == 1` and `request.pressure.solved` — it names no problem kind in an
identifier, so checker rule 7 (`tools/check-engine-architecture.py:190`) passes it, and it
is nonetheless a union of special cases. Each of the six had its own tests and each test
passed. **measured**

What has to be established instead: that one row system, compiled from a specification by
one rule, is what every equilibrium problem becomes, and that every point of the
specification space that *should* be solvable compiles to a square, full-rank system while
every point that should not is rejected by name. That is a property of the compile step and
of the rows, and it can be checked for thousands of specifications without solving one.

Two things this document does not claim, stated up front (§8): the *absence* of a hidden
route is not observable from outside the code, and a lattice of cells is a longer list, not
a proof. What the space buys is that a special case now has to survive every cell its axes
generate, inside a 2,250-line budget, instead of the cells someone remembered to write.

## 2. One system, one activation rule

### 2.1 Variables

Per declared phase `p` with support `S_p`:

| coordinate | meaning | status the specification may give it |
|---|---|---|
| `x_p` (as `ln x_{p,i}`, `i ∈ S_p`) | composition | always free |
| `N_p` | phase amount, mol | **free**, or **pinned** to `β_p N_feed` (`β_p = 0` is an incipient phase) |
| `V_p` (as `ln η_p`) | phase volume | always free |

System-wide: `T` and `P`, each **fixed** or **free**; feed `n_feed` (fixed, in any basis);
the balance matrix `B` (identity when no reactions; otherwise the conserved moieties, whose
nullspace is the reaction span); charges `z`; standard-state offsets `g°`.

This differs from `formulation.md` §2 in one place, and the difference is load-bearing.
The formulation carries `ln n_{p,i}` per phase, and its §3 then needs a *second variable
kind*, `x_inc`, for a phase of zero amount, because `ln 0` is not a coordinate. That second
kind is where `TRANSFER.md` §"Stage 6 versus Stage 7" concludes that an incipient phase
"cannot be a row of the minimizer". It is right about the coordinate and wrong about the
formulation: with the phase amount `N_p` separated from composition, `n_{p,i} = N_p x_{p,i}`
is a linear change of variables in log space (the tape does not care), and pinning `N_p`
is a specification that removes a variable, not a limit that breaks one. **derived**

### 2.2 Rows

With multipliers `λ` on the balance and `φ_p` on each charge row, the rows are the KKT
conditions of `G = Σ_p [A_p + P V_p] + Σ g° n` (`EqID: system_gibbs_helmholtz_coordinates`),
written **intensively** — the stationarity row in `n_{p,i}` divided by `N_p`:

| row | equation | count | EqID |
|---|---|---|---|
| chemical potential | `μ_{p,i}(T, V_p/N_p, x_p) − λ·B_i − z_i φ_p = 0`, `i ∈ S_p` | `Σ_p |S_p|` | `chemical_potential_stationarity` |
| pressure equality | `−∂A_p/∂V_p − P = 0` | one per phase | `pressure_equality` |
| material balance | `Σ_p N_p B x_p − B n_feed = 0` | `rank B` | `material_balance` |
| electroneutrality | `z · x_p = 0` | one per phase whose support carries charge | `phase_electroneutrality` |
| normalisation | `Σ_i x_{p,i} − 1 = 0` | one per **pinned** phase | `incipient_normalization` |

For a free phase, normalisation is implied by the definition `N_p = Σ_i n_{p,i}` and is not
a row. For a pinned phase it is the row that replaces the missing amount variable. The
chemical-potential row is the *same equation* for a free and a pinned phase because `μ` is
intensive; only the intensive normalisation makes that literal. Euler's relation then
gives, for every phase, `g_p = λ·B x_p + φ_p z·x_p` as a consequence — which for a pinned
phase is the statement that the phase is at the threshold of existing (`TPD = 0`), and for
a free phase is the usual equilibrium condition. Nothing is imposed twice. **derived**

Reaction affinity is not a row: `Σ_i ν_i μ_i = λ·(Bν) + φ z·ν = 0` because `Bν = 0` and
`z·ν = 0` (`EqID: reaction_affinity`). It is a consequence, and §5 uses it as a check.

### 2.3 The count rule

Let `k` be the number of pinned phases and `f` the number of free intensives (`T`, `P`).
Counting unknowns against rows in §2.1–2.2 (the multipliers pay for the constraint rows
exactly, so they drop out):

    rows − unknowns = k − f

The system is square **iff `f = k`**. Necessary, not sufficient: the rank of the active
Jacobian at the start point is the sufficient check (§4.1), because a square system can
still be singular by construction — all phases pinned, a redundant reaction, a redundant
charge row. **derived**

Read off the rule, with no route named anywhere:

| `k` | `f` | what it is |
|---|---|---|
| 0 | 0 | the flash: fixed `T`, `P`, feed; `Π` phases free |
| 1 | 1 | a boundary: one phase pinned at `β = 0`, `T` or `P` free — bubble, dew, cloud point, saturation (`C = 1`) are one cell family, distinguished only by *which* phase is pinned |
| 1 | 1 | the same with `β > 0`: the fixed-fraction flash (`T` at which a fraction `β` has vaporised at fixed `P`). No route in the old tree; the same rows |
| 2 | 2 | two phases pinned, `T` and `P` free: the three-phase point at a fixed reference composition |
| 0 | 1 | under-specified: rejected, "free `T` with no pinned phase" |
| 1 | 0 | over-specified: rejected, "phase `q` pinned with `T` and `P` fixed: free one" |

The owner's insight holds exactly under §2.1: bubble point = the flash's specification with
`N_V` pinned to zero and `P` freed; dew point = `N_L` pinned instead. The framing in the
brief, "a dew point pins it to one", is redundant — pinning one phase's amount to the whole
feed is pinning the other to zero, and the rule counts pins, so it is stated once.

Rachford–Rice is the balance rows reduced to a scalar in `N_V` *after* the
chemical-potential rows have been replaced by fixed K-values. In this system it is a start
generator at most; if `gibbs_minimizer.cpp` or `phase_boundary.cpp` contains a scalar solve
in a phase fraction, that is a second route. Rule 7 cannot see it; §4.3's checks can.

### 2.4 Two drivers, one assembly

A specification with `k = 0` is a minimisation and goes to Ipopt: objective `G`, equality
constraints = balance, charge and normalisation rows, bounds on `ln η_p`, exact Lagrangian
Hessian from the tape's `∂²A_p` per phase. A specification with `k ≥ 1` is **not** a
minimisation — `G` is not minimised over `P`, and the pinned-multiplier condition is an
equation, not an optimality condition — so it is a square Newton on the same rows with the
same Jacobian. It cannot be an Ipopt NLP with zero objective either: Ipopt's Hessian of
`Σ σ_i F_i` needs third derivatives of `A`, which the tape does not carry. **derived**

So "one solver" means one row assembly `F(y; spec)` with its exact Jacobian `F_y`, consumed
three ways: Ipopt's objective/constraint/Jacobian/Hessian callbacks are views of it; the
boundary Newton iterates on it; `implicit_sensitivity.cpp` inverts it. The Hessian is
consumed in all three. §4.3 has the check that the three see one matrix.

## 3. The space

### 3.1 Axes

| axis | values | what it changes |
|---|---|---|
| phase set, by identity | `{L}`, `{V}`, `{V,L}`, `{L,L}`, `{V,L,L}` | a phase's identity is (support, branch hint, amount status). `V` versus `L` is a **start** attribute — the density-branch hint for the initial `η_p` — and changes no row. `{L,L}` is in the space precisely because it is where a vapour-like reference root would be assumed if any part of the start or the seed assumed one |
| amount status per phase | free; pinned `β = 0`; pinned `0 < β < 1` | `k`, the normalisation row, the balance's view of the phase |
| `T`, `P` | each fixed or free | `f` |
| reactions | none; one; two independent; two plus a dependent third; one that does not conserve charge | `rank B`, `g°`, the rejection |
| ions | absent; present in a support with both signs; a charged feed; a support with one ion sign; an incipient neutral phase beside a charged finite one | the charge rows, the rejection |
| feed basis | mole fractions; amounts; mass fractions; molalities | **nothing in the compiled system** — that is the check |
| topology | declared; omitted | declared is the main path; omitted is rejected by name until Stage 11, then accepted with `k = f = 0` |
| feed coordinate free with a target row | (later; the old tree's "affine feed boundary") | the same count rule with a target row on the right: a freed feed coordinate is paid for by one target equation. Recorded so the rule is seen to extend; not built before Stage 12 needs it |

Solids and critical points are outside the space, as `formulation.md` §7 says: different
state variables, and second-order conditions as equations.

### 3.2 Unreachable cells, and what each must say

These carry the weight. A compile step that accepts them silently is the old tree's
"it never defaults" claim without its test.

| cell | rejection names |
|---|---|
| `f < k` | the pinned phase and "free `T` or `P`" |
| `f > k` | the free intensive and "fix it or pin a phase amount" |
| every phase pinned (`k = Π`) | "no free phase carries the feed" — square by the count, singular by rank; this cell is why the rank check is not redundant |
| dependent reaction set | the reaction that lies in the span of the others |
| `z·ν ≠ 0` | the reaction and "does not conserve charge" |
| `z·n_feed ≠ 0` | "feed carries net charge" |
| a charged support with one sign | the phase and "support cannot be electroneutral" |
| a species in no support | the species |
| an empty support | the phase |
| topology omitted, before Stage 11 | "phase count unknown; declare phases" |
| topology omitted with `k ≥ 1` or `f ≥ 1`, at Stage 11 | "discovery is a fixed-`T`,`P` problem" |

One cell is not a rejection but a decision the compile step must make and say: a species
with zero feed that no reaction produces has zero amount in every phase, and `ln 0` is not
a coordinate. Compile removes it from every support and reports the removal.

### 3.3 The compile-tier lattice is exhaustive, not sampled

The core axes multiply to about 3,200 specifications (5 phase sets × ≤ 8 amount-status
patterns × 4 intensive patterns × 5 reaction values × 4 ion values); each compiles in
microseconds and the rank of a matrix under 40 × 40 costs the same. Feed basis is an
invariance applied to every cell (four compiles, compared), not a fifth axis. Nested loops
over four enums are the generator; there is no property-based framework and no random
specification — the space is small enough to walk. The one place randomness would have
mattered, the point at which the Jacobian rank is taken, uses a fixed generic point
(`x_i ∝ 1 + 0.37 i`, `η_p = 0.25 + 0.05 p`) so that no symmetry of an "equal amounts"
point manufactures a rank drop. **guessed** on the count; **derived** on the cost.

## 4. What each check asserts, and the mutation that makes it red

`TESTING.md` binds hardest here. A synthetic system has no measured answer, so no check
below asserts a value. Each asserts a property that holds for every correct implementation
of the equations, and each carries the mutation the brief must apply at acceptance.

Every quantity asserted is coordinate-free: **nullity** (unknowns minus rank), never the
NLP's variable count — a reimplementation that eliminates the density inline has fewer
variables and the same nullity. Row *labels* by `EqID` are asserted only where the
formulation registers the row.

### 4.1 Compile tier — no solve; lands with Stage 6

| # | claim | asserts | mutation (must go red) |
|---|---|---|---|
| C1 | count rule | for every lattice cell, compile's verdict (square / over / under) equals `f − k` computed in the test from the specification alone | in the compile step, count a pinned phase as free |
| C2 | reject by name | every unreachable cell throws, and the message contains the offending phase / reaction / species name from §3.2 | swap any two messages |
| C3 | rank | for every square cell, `F_y` at the generic point has nullity 0; for `k = Π`, nullity ≥ 1 and the rejection names it | drop the normalisation row of pinned phases |
| C4 | feed-basis invariance | the same feed as mole fractions, amounts, mass fractions and molalities compiles to identical rows and identical `B n_feed` to 1e-12 | molality conversion omits the solvent mass |
| C5 | reaction-basis invariance | reactions `{ν₁, ν₂}` and `{ν₁ + ν₂, ν₁ − ν₂}` (and any invertible mixing) give the same `B`, the same nullity, the same `g°` transform | build `B` from the reaction list in order instead of from its span |
| C6 | pin identity | the boundary specification's residual vector at a point equals the flash specification's residual vector at the same point with `N_q = 0`, on every shared row — the owner's "same rows, one entry differs", made numerical | the balance keeps the pinned phase at its start amount |
| C7 | permutation | relabelling components, or reordering phases, permutes the rows and leaves nullity and residual norms unchanged to 1e-14 | index a phase where a component was meant in any row |
| C8 | chart-or-row | an electrolyte flash compiles with nullity 0 — this catches the double count in `formulation.md` (§2 eliminates an ion by the electroneutral chart *and* §4 lists the charge row; both together are a redundant row) | add the eliminated coordinate back |

C6 is the check the brief asked for. It has a limit worth stating: if `phase_boundary.cpp`
assembles rows through the same function as the minimizer, no mutation of that function can
make C6 red, because both sides move together. C6's mutation therefore lives in the
specification handling, and the checks that catch a *bypass* of compile are S3–S5 below.

### 4.2 Solve tier — the minimal input set

Five parameter sets, four of which already exist in `engine/tests/native/model_support.hpp`
(**measured**: `ethane()` :42, `neutral_binary()` :51, `associating_binary()` :60,
`aqueous_electrolyte()` :87, `nine_species()` :116). One carries from `main`:
`engine/tests/eos/_cases.py:373-397`, the immiscible ternary — three components with
`m = 1.5`, `σ = 3.5, 3.7, 3.9 Å`, `ε/k = 300 K`, `k_ij = 0.20` on every pair. It is the LL
fixture, and with reactions declared among its own components (`α ⇌ β`, `β ⇌ γ`, and
`α ⇌ γ` as the dependent third) it is the reactive fixture too. One input, three axes.

Each fixture's cell must *exist* at the chosen state, and no measurement says so. The
existence certificate is the engine's own, one-sided: a two-phase solve with `G` below the
one-phase `G` at the same feed proves the one-phase state unstable (S6). A cell whose
certificate fails is a fixture problem, recorded in `STAGE-LOG.md`, not a code failure.

| id | fixture | cells | stage |
|---|---|---|---|
| ethane | `ethane()` | `{V,L}`, `V` pinned, `P` free; then `T` free | 7 |
| binary | `neutral_binary()` | `{V,L}` flash inside the envelope; bubble; dew; `β`-flash; both intensives | 8, 9 |
| ternary | immiscible ternary | `{L,L}` flash; `{V,L,L}` flash at low `P`; `{L,L}` with `L₂` pinned, `T` free; `{V,L,L}` with `V`, `L₂` pinned, `T` and `P` free; one phase with `α ⇌ β`; `{L,L}` with both reactions | 8, 9 |
| electrolyte | `aqueous_electrolyte()` | `{L}`; `{V,L}` with `V` support `{water}` and `V` pinned, `P` free — the MEA shape at three species | 8, 9 |
| nine species | `nine_species()` | compile tier only: the composite cell's rows, rank, nullity, and reject-by-name with the bundle's two reactions. The solve is Stage 12's gate | 6 |

| # | claim | asserts | mutation |
|---|---|---|---|
| S1 | rows satisfied | every solved cell's residual vector is within `derived(minimizer_convergence_tolerance)` of zero on every row — including rows the driver did not iterate on (affinity for a reactive cell; Euler's `g_p = λ·Bx_p` for every phase) | drop one balance row from the Ipopt constraint set: affinity goes nonzero |
| S2 | local minimum | for `k = 0`, the reduced Hessian `Zᵀ ∇²L Z` is positive definite (`EqID: reduced_hessian_stability`) | flip the sign of the `P V` term in the objective |
| S3 | dew ↔ bubble | for the binary at `T`: bubble from liquid `x` gives `(P*, y*)`; dew from `y*` must give `(P*, x)` to the convergence tolerance. An exact symmetry of the formulation, and it fails on any asymmetry between "which phase is pinned" | start the pinned phase on the vapour branch regardless of its declared identity |
| S4 | `β → 0` | the fixed-fraction flash at `β = 10⁻², 10⁻⁴, 10⁻⁶` converges to the bubble point's `T` with the expected first-order rate | the balance omits the pinned phase for `β > 0` |
| S5 | round trip | solve `P` at fixed `T`; solve `T` at that `P`; recover `T` | freeze `T` in the Jacobian when `T` is the unknown |
| S6 | instability certificate | `G₂ < G₁` at the ternary feed; `G₃ ≤ G₂ + tol` for the three-phase cell | return the one-phase start as the result |
| S7 | lever rule | two feeds on the `{L,L}` tie line give the same two phases to 1e-8 | seed each phase from its own feed instead of from the tie line |
| S8 | extensivity | every solved cell with the feed doubled: amounts double, every intensive unchanged | an absolute trace floor in the start |
| S9 | sensitivity | `dy/dθ` from `implicit_sensitivity.cpp` against a centred difference of the solve in `θ` to 1e-6 (the Stage 3 precedent for the FD floor) | flip the sign of `F_θ` |
| S10 | one Jacobian | at the binary flash solution, Ipopt's `eval_h` with `σ = 1` and the converged multipliers equals the stationarity block of the boundary driver's `F_y` and of the sensitivity's `F_y`, elementwise to 1e-12 | scale `eval_h` by `0.5` |
| S11 | ideal-gas limit | the reactive one-phase cell at `P → 0` approaches the closed form `K = exp(−Δg°/RT)` in mole fractions (a limit; the residual terms vanish, which Stage 2 already checks) | perturb the `g°` transform |
| S12 | discovery agrees (Stage 11) | topology omitted at the ternary feed finds the phase count of the declared `{L,L}` flash and the same `G` | accept the first sample as a witness |

**Reported, never asserted**, per cell: Ipopt iterations, Newton iterations, evaluations,
`hessian_calls`, recordings, wall. `hessian_calls ≥ 1` per iteration on every cell is how
"the exact Hessian is consumed everywhere" is *seen*; it is not a pass condition. One
architecture rule is cheap enough to add to the checker: `gibbs_minimizer.cpp` defines
`eval_h` and contains no `limited-memory` string. Lexical and weak, and it says what the
design intends.

### 4.3 What the tiers establish together

C1–C8 prove that the compile step is one rule over the whole space and that the rejections
work. S1, S2, S8 prove each solved point is what it claims. S3–S5 are the checks that see a
*bypass*: a second bubble-point path that never enters compile still has to be symmetric
under pin exchange, continuous in `β`, and invertible in `T`/`P`, and a route written for
one application will usually fail one of them. S10 is the only check on the derivative
path itself. The artifact's "seven genericity equivalences" (round trip, reaction basis,
extensivity, a problem stated two ways) are not seven checks here; they are C4, C5, S5, S8
applied to every cell they fit, which is what makes them coverage rather than spot checks.

## 5. From specification to NLP

    Problem  ──compile──►  CompiledSystem  ──►  Ipopt (k = 0)  ──►  EquilibriumResult
                                │                     │
                                ├──►  Newton (k ≥ 1)  ┤   both: F(y), F_y from one assembly
                                │                     │
                                └──►  implicit sensitivity: −F_y⁻¹ F_θ

1. **`Problem`** (`problem.hpp`) is data: phases with support, branch hint and amount
   status; `T`, `P` each fixed or free; feed in a declared basis; reactions; active
   parameters. It has no method that solves.
2. **compile** (`request_compile.cpp`) resolves the feed basis to amounts (native — see §7,
   item 4), builds `B` from the reaction span by SVD, applies the count rule, assembles the
   active row list with its `EqID`s, takes the rank at the start point, and either returns
   a `CompiledSystem` or throws naming the deficiency. It is the one file allowed to branch,
   and only on status flags.
3. **assembly** (`gibbs_minimizer.cpp`): `F(y)` and `F_y` from one tape evaluation per
   phase per iteration — value, gradient, Hessian of `A_p` in `(T, V_p, n_p, θ)` — plus the
   linear constraint blocks. The Hessian of `A_p` is the tape's reduced Hessian; the
   site-fraction elimination happened in Stage 3 and nothing here sees sites.
4. **drivers**: Ipopt for `k = 0` with bounds on `ln η_p` and `ln x_{p,i}`, exact Hessian;
   damped Newton on `F` for `k ≥ 1`, seeded from the flash solution or the previous
   boundary point when the caller is sweeping (`TRANSFER.md` Stage 7).
5. **result**: amounts, volumes, intensives, multipliers, every residual row, the reduced
   Hessian's smallest eigenvalue, counters, a named status. Sensitivity when asked.

The coverage exercises this path and not a bypass of it because C1–C8 act on step 2's
output, S1 reads step 5's residuals for rows step 4 did not iterate on, and S10 holds the
three consumers of step 3 to one matrix.

## 6. What of `formulation.md` should land, and how

The branch `docs/equilibrium-formulation` (PR #194, old repository) was merged at
`c24da84` and reverted at `363a739` with no reason in the message. **measured**: the
registry sync (`engine/scripts/docs/sync_equation_registry.py --check`) already fails on
`greenfield` with a traceback, because every `implementation_refs` entry points into
`src/epcsaft/_native/`, which Stage 1 deleted; the merge would not have changed that.
Whatever the reason for the revert, landing needs the following.

**Keep, as equations:** §1 (one problem), §2 (free energy and coordinates, with §2.1's
change: `N_p` separate from composition; the `x_inc` variable kind goes), §4's row table,
§5 (KKT and the two consequences), §6 (stability, verbatim), §8 (implicit sensitivity). The
21 `Equilibrium`-section identifiers (**measured** on the branch's registry; 15 registered
by `formulation.md` §9 and 6 by `phase-stability.md` §8) keep their names and statuses.

**Replace §7.** Its nine-row "specification patterns" table is the taxonomy in prose — "T,
reference-phase composition, one incipient phase" is a route with a name. §2.3's count rule
and §3.2's rejection table replace it. `formulation.md` §7's two "fails" rows survive as the
`f ≠ k` lines.

**Resolve §2 against §4.** §2 says a charged phase uses the electroneutral chart (`C − 2`
coordinates); §4 lists `phase_electroneutrality` as a row. Both is a redundant row (C8
catches it). Recommendation: the row, with its multiplier `φ_p`, and no elimination in the
minimizer — the charged chart is then the neutral chart plus one row, which is smaller
than "plus one elimination" and has no `eliminated_index`. `phase-stability.md` §2 keeps
the chart for the TPD search, where trial compositions must be neutral by construction;
`electroneutral_amount_chart` stays registered for that owner.

**Drop:** §4's sentence naming `compile_phase_direction_system` and the two Python rank
functions (deleted files); §9 entirely (owners that no longer exist).

**Registry:** land the 21 entries with empty `implementation_refs`; extend the sync script
so that a non-`Documentation-only` entry with no owner is a *warning* until the stage that
annotates it, and its `MATH_OWNER_FILES` lists `native/**` instead of the 21 deleted paths.
The `// EqID:` annotation system is worth keeping exactly as `CARRY-OVER.md` says; 58 of
them are already in `native/` (**measured**) and the script cannot currently read them.

## 7. Where the brief's framing is corrected

1. **"An incipient phase cannot be a row of the minimizer"** (`TRANSFER.md`, Stage 6 versus
   Stage 7). It cannot be an *Ipopt variable* of the minimizer, and the boundary is not a
   minimisation. It is a row of the same system, and Stage 9 is a driver on that system,
   not a second row set. If Stage 9's `phase_boundary.cpp` has its own chemical-potential
   residuals, the design has failed in exactly the way the owner described.
2. **"Visible from the compiled row set alone."** Only if the boundary passes through
   compile. A path that bypasses compile is invisible to C6 and visible to S3–S5.
3. **"A dew point pins the phase fraction to one."** Redundant under a per-phase amount
   status; stated once as "pins the other phase to zero" so the count rule has one case.
4. **Feed-basis conversion in Python.** `TRANSFER.md`'s target layout gives `loader.py`
   250 lines for "mass-fraction and molality feeds". That is arithmetic, which the same
   document forbids in Python. It moves into compile (native), where C4 can run under the
   `CHECK` harness and Python stays a mapping layer with no exception clause.
5. **Rule 7 is lexical.** Agreed, and it stays: it is cheap and it names the intent. The
   non-lexical replacement is C1–C8 plus S3–S5, and §8 says what even those miss.
6. **"Seven structural equivalences."** Spot checks on path pairs, as the brief says. Here
   they are properties applied per cell (§4.3).
7. **The test budget.** `engine/tests/native/` holds 1,187 lines after two closed stages
   (**measured**), against a 1,650 total for the whole engine. This design adds about 560
   (§9). With the oracle-row tests still to write, the engine lands near 2,400 lines of
   tests, not 1,650. `TESTING.md`'s rule is that a stage over its share must show its
   checks pass the diagnostic, not that the budget grows; every check above states the
   property it tests and the mutation that fails it, and the number is recorded rather than
   trimmed. `model_floor.cpp` (186) is a measurement harness, not a test, and arguably sits
   outside the budget; that does not close the gap.

## 8. What completeness does not establish

- **Starts.** A cell solved from one start proves one basin. Continuation across a sweep,
  the shifted-previous start family, and the fixed start list of Stage 11 are not covered by
  any check above; a cell can pass every property and fail from the next caller's start.
- **Budgets.** Iteration caps and the named stops that follow them are reported, never
  asserted. A cell that converges in 400 iterations passes S1.
- **Shared physics.** Every cell inherits one `μ`. A wrong `μ` — a sign in the Born term,
  a permittivity form filled wrongly by the model layer — makes every cell consistently
  wrong, and every row-satisfied check passes, because the rows are satisfied. Only the
  oracle rows see it. The two are complementary and neither substitutes for the other.
- **Absence is not observable.** No external check proves there is no hidden route. S3–S5
  catch the routes that behave differently; a hidden route that behaves identically is,
  by that token, not a defect anyone can find, and the line budget is the only bound on it.
- **The lattice is a list.** Once it exists, code will be written against it. It is a
  better list — generated, so a new axis value regenerates thousands of cells — but it is
  the same kind of comfort one level up. The count rule and the rejections are the part
  that is not a list; they are a closed form.
- **Existence is self-certified.** A synthetic cell "exists" because the engine's own
  `G₂ < G₁` says so. That is one-sided: it proves instability, not that the two phases
  found are the right two.
- **Tolerances define satisfaction.** "Satisfies its own rows" means within Ipopt's and the
  policy header's tolerances. A tolerance loosened by mistake weakens every S-check at
  once; S1's mutation is the guard, and it is one mutation.
- **The discovery cell asserts a rejection until Stage 11.** Between Stage 6 and Stage 11,
  "topology omitted" is proven only to fail by name.

## 9. Cost and placement

All native, under the existing `CHECK` harness (`engine/tests/native/physics_checks.cpp:18`);
no framework, no fixture library, no Python until Stage 4's bindings exist, which under the
revised order (`TRANSFER.md` §"Execution order") is before Stage 6 anyway. **guessed**:

| file | lines | ceiling | stage | holds |
|---|---|---|---|---|
| `tests/native/equilibrium_support.hpp` | 130 | 180 | 6 | includes `model_support.hpp`; adds the immiscible ternary, its reactions and `g°`, specification builders, the four-enum lattice walk, the generic point, the count-rule predictor |
| `tests/native/specification_checks.cpp` | 180 | 240 | 6 | C1–C8 |
| `tests/native/equilibrium_checks.cpp` | 130 | 180 | 8, then +15 at 11 | S1, S2, S6, S7, S8, S11 on the flash cells; S12 at Stage 11 |
| `tests/native/boundary_checks.cpp` | 120 | 160 | 7 opens it (ethane, S5, S9 at `C = 1`, ~40); 9 completes it | S3, S4, S5, S9, S10, the `{L,L}` and `{V,L,L}` boundary cells |
| `tools/check-engine-architecture.py` | +12 | — | 8 | rule 8: `eval_h` present, no `limited-memory` |
| **total** | **~570** | **~770** | | |

Order inside Stage 6: `equilibrium_support.hpp` and C1–C3 first, before any solver exists —
compile and the count rule are the first thing the stage writes, and C1–C3 are green on
the compile step alone. That is also the cheapest possible preflight for the stage brief:
if C1's predictor and the compile step disagree on a cell, one of them misreads §2.3, and
that is found before a line of `gibbs_minimizer.cpp` is written.

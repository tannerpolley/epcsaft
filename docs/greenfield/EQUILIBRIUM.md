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
| `N_p` | phase amount, mol | **free**, or **pinned** to `β_p N_feed`, where `N_feed` is the feed's total mole number as given — a specification, whether or not the reactions conserve total moles (`β_p = 0` is an incipient phase) |
| `V_p` (as `ln η_p`) | phase volume | always free, bounded below |

The driver's unknowns are `ln x_{p,i}`, `ln η_p`, `N_p`, and the free intensives. Bounds:
`ln x_{p,i} ∈ [ln STABILITY_COMPOSITION_FLOOR, 0]`, `ln η_p ∈ [ln 10⁻¹², ln PACKING_FRACTION_CEILING]`
with the ceiling `0.74` (`policy::PACKING_FRACTION_CEILING`, an `AlgorithmTolerance`: close
packing, above which no fluid root exists and `validate_point` refuses anyway),
`N_p ≥ 0` with no upper box bound (with reactions the mole total is not conserved, mass is, and the balance rows carry it; a bounded KKT system needs no box above, and for a non-reactive cell `N_p ≤ Σ n_feed` follows from the same rows), `T, P > 0`. A plain-volume coordinate with no ceiling admits the
`v → ∞`, `P → 0` solution of the pressure row, which attempt 3 returned as a converged bubble
point. This is held by reading the driver's `get_bounds_info`, like "one assembly": no check
can distinguish a bounded coordinate from an unbounded one that happened not to run off.

The sensitivity levers are `T`, `P` and the species feed `n_feed,i`, one column per species
(§12): `dy/dθ` is one call with one column per lever, and the model's active parameters are
the same call with more columns.

System-wide: `T` and `P`, each **fixed** or **free**; feed `n_feed` (fixed, in any basis);
the balance matrix `B` (identity when no reactions; otherwise the conserved moieties, whose
nullspace is the reaction span); charges `z`; standard-state offsets `g°`.

This differs from `formulation.md` §2 in one place, and the difference is load-bearing.
The formulation carries `ln n_{p,i}` per phase, and its §3 then needs a *second variable
kind*, `x_inc`, for a phase of zero amount, because `ln 0` is not a coordinate. That second
kind is where `CORE.md` §6 concludes that an incipient phase is "not a phase of the
minimizer". It is right about the coordinate and wrong about the
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
| electroneutrality | `z · x_p = 0` | one per ion-containing phase **except one phase of nonzero amount**; all of them if every ion-containing phase is pinned at zero | `phase_electroneutrality` |
| normalisation | `Σ_i x_{p,i} − 1 = 0` | one per phase | `phase_normalization` (new; see §6) |

**The rows are dimensionless.** `nlp_scaling_method = none` is a commitment (`CORE.md` §3.4),
so the rows scale themselves. The chemical-potential rows are already in units of `RT`. The
pressure row is written `Z_p(T, v_p, x_p) − P v_p / (RT) = 0` (`EqID: pressure_equality`, the
same equation multiplied by `v_p / RT`), which is `O(1)` on every branch. The balance rows are
divided by `Σ_i n_{feed,i}`; the amounts `N_p` stay extensive and `n_feed` is never normalised,
so S8's mutation is unchanged. Normalisation and charge rows are `O(1)` as written.

The normalisation row is active for **every** phase, so the coordinates `(N_p, x_p, V_p)`
are the same for every phase and a pinned phase differs from a free one by exactly one
removed variable — nothing else changes. (Writing a free phase in `n_{p,i}` and dropping
its normalisation row is a legitimate implementation choice with the same nullity; it is
not the same *statement*, and the statement is what the checks assert.) The
chemical-potential row is the *same equation* for a free and a pinned phase because `μ` is
intensive; only the intensive normalisation makes that literal.

**How the intensive row is evaluated — the rule that makes `N_p = 0` a plain zero.** The
tape takes `(T, V, n)` extensively, and `V_p / N_p` is `0/0` for a pinned phase. So the
assembly evaluates every phase's tape at **unit amount**: `n = x_p`, `V = v_p` (the molar
volume, which `ln η_p` already is up to the packing coefficient), and multiplies extensive
results by `N_p` afterwards. `μ`, `P` and `g` are homogeneous of degree zero and come out
unchanged; `A_p`, `G` and the Hessian blocks scale by `N_p` — in the `(N_p, x_p, v_p)`
chart the assembly works in. A driver that works in extensive `(n_p, V_p)` coordinates
sees the inverse scaling, and the driver states which chart it uses. A pinned phase then
contributes `0 · (…)` to the balance and to `G`, its rows are evaluated at the same point
as any other phase's, and no `ln 0` is ever formed. This is the one implementation
decision the count rule depends on; an assembly that evaluates at `(T, V_p, N_p x_p)`
cannot represent a pinned phase and will grow a second path. **derived**

In Ipopt's own KKT system for `k = 0` the composition rows carry the factor `N_p` and the
normalisation multiplier `ψ_p`; at the solution `ψ_p = 0` follows from the `N_p`
stationarity row, so the intensive rows above hold there as consequences. S1 asserts them
as such. The two drivers therefore see rows that differ by a scaling of `N_p` per phase,
which is why §4.2's S10 is a finite-difference agreement and not an elementwise equality
between drivers.

Euler's relation gives, for every phase, `g_p = λ·B x_p + φ_p z·x_p` as a consequence of
the intensive rows — which for a pinned phase is the statement that the phase is at the
threshold of existing (`TPD = 0`), and for a free phase is the usual equilibrium
condition. Nothing is imposed twice. **derived**

**Why one charge row fewer than the ion-containing phases.** The feed is neutral and
every admitted reaction conserves charge, so the charge vector `z` lies in the row space
of `B` and overall neutrality `Σ_p N_p z·x_p = 0` is already a consequence of the balance
rows. With one ion-containing phase of nonzero amount its neutrality follows; with `π_ch`
such phases only `π_ch − 1` rows are independent, and their multipliers `φ_p` are the
electrostatic potential differences between the phases, the row-less phase's potential
being arbitrary. **The row-less phase must have nonzero amount**: the balance sum weights
each phase by `N_p`, so a phase pinned at zero drops out of it and its neutrality follows
from nothing. The rule is therefore "every ion-containing phase except one of nonzero
amount", and if every ion-containing phase is pinned, all of them carry the row. The
three sources below treat finite-phase flashes only and are silent on this case;
S1 asserts `|z·x_p|` within the floor for every ion-containing phase so that a wrong
assignment is red at the solve tier.
This is Tsanas et al. 2019 §2.2.1 and 2022 §2.4, and Ascani, Sadowski and Held 2022 eq. 5
(`π − 1` rows "due to the electroneutrality of the whole system"). Writing one row per
charged phase, as the first draft of this table did, adds a row that is a linear
combination of the balance rows: the compile step's rank check (C3, C8) rejects a valid
electrolyte flash, and an implementer who did not know why would loosen the rank
tolerance. **derived**, three sources.

Reaction affinity is not a row: `Σ_i ν_i μ_i = λ·(Bν) + φ z·ν = 0` because `Bν = 0` and
`z·ν = 0` (`EqID: reaction_affinity`). It is a consequence, and §5 uses it as a check.
The multipliers `λ` are the chemical potentials of the conserved moieties, and
`G_min = λ·B n_feed` (Tsanas et al. 2017 eq. 29–30); that identity is one more consequence
S1 can read off a solution.

### 2.3 The count rule

Let `k` be the number of pinned phases and `f` the number of free intensives (`T`, `P`).
Counting unknowns against rows in §2.1–2.2 (the multipliers pay for the constraint rows
exactly, so they drop out):

    rows − unknowns = k − f

The system is square **iff `f = k`**. Necessary, not sufficient: the rank of the active
Jacobian at the start point is the sufficient check (§4.1), because a square system can
still be singular by construction — all phases pinned, a redundant reaction, a charge row
written for a phase whose neutrality the balance already implies. **derived**

One more condition is the phase rule in rank form, and it is what makes the fixed-`T`,`P`
three-phase flash of a binary unreachable. `N_p` appears only in the balance rows, so the
`Π − k` free-amount columns can be independent only if `Π − k ≤ rank B`: free phases at
most the number of independent moieties. The `C = 1` `{V,L}` flash at fixed `T`, `P`
(`2 > 1`) and the binary `{V,L,L}` flash (`3 > 2`) are square by the count and singular
at every point; pinning a phase or freeing an intensive restores them (the three-phase
point of §2.3 is `k = 2`). A reacting ternary with one reaction has `rank B = 2`, so its
three-phase flash is unreachable too. Compile rejects these by name (§3.2). **derived**

Read off the rule, with no route named anywhere:

| `k` | `f` | what it is |
|---|---|---|
| 0 | 0 | the flash: fixed `T`, `P`, feed; `Π` phases free. **`Π = 1` is a member**: the single-phase flash, whose solution is the feed composition and the density root — see §2.5 |
| 1 | 1 | a boundary: one phase pinned at `β = 0`, `T` or `P` free — bubble, dew, cloud point are one cell family, distinguished only by *which* phase is pinned. **`C = 1` is a member**: pure saturation — see §2.5 |
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

### 2.4 One driver, one assembly

Every cell is solved by Ipopt on one `TNLP` class whose two inputs come from compile as
data: the **objective** and the **constraint set**.

- `k = 0` is a minimisation: objective `G`, constraints = balance, charge and normalisation
  rows, bounds on `ln η_p`, `ln x_{p,i}` and `N_p ≥ 0`, exact Lagrangian Hessian from the
  tape's `∂²A_p` per phase. Minimising `G` rather than solving its KKT rows is what tells
  the trivial solution and every saddle apart: they satisfy the rows, and only descent on
  an objective refuses them.
- `k ≥ 1` is **not** a minimisation — `G` is not bounded below in `P` — but it is square by
  the count rule, and a square feasibility problem is one Ipopt solves as Newton. The
  objective is zero and *every* row of §2.2 is a constraint. Ipopt's step solves
  `[W + δI, Aᵀ; A, 0] [dx; dλ] = −[∇f + Aᵀλ; g]`, and with `A` square and nonsingular the
  second block row alone fixes `dx = −A⁻¹ g`: the exact Newton step on the residuals,
  whatever the Hessian `W`. The multipliers of a square feasibility problem vanish at the
  solution, so the exact Lagrangian Hessian there is zero, and the Hessian callback
  returns zero for `k ≥ 1`. No third derivative of `A` is needed anywhere; the earlier
  draft's claim that a boundary "cannot be an Ipopt NLP" conflated the constraint
  Hessian with the primal step. **derived**

What one driver buys: bounds, line search, regularisation, `left_domain`, and one
`Convergence` record are Ipopt's for every cell, and no second solver with its own loop,
tolerance and stop names is written — which is the disease in another form. The "pick"
is compile filling two fields of the compiled system, not a branch in a driver.

So "one solver" means one row assembly `F(y; spec)` with its exact Jacobian `F_y`, consumed
two ways: the `TNLP` callbacks are views of it, and `implicit_sensitivity.cpp` inverts it.
S10 holds each consumer to finite differences of its own rows; that there is one assembly
is a placement decision (§5), not a check.

### 2.5 Three levels, one row set — and the two cells the design might quietly skip

    flash (feed, T, P)        every N_p free                      k = 0, f = 0
      bubble / dew            one N_p pinned at 0, T or P freed   k = 1, f = 1
        saturation            the same, at C = 1                  k = 1, f = 1

The nesting is a claim about rows, so it is checked as one. Count the `C = 1` boundary
cell — one component, `{V, L}`, `V` pinned at zero, `P` free — against §2.1–2.2 with no
special case anywhere: unknowns `x_V, V_V, x_L, N_L, V_L, λ, P` (seven); rows: two
chemical-potential, two pressure-equality, one balance (`B = I₁`), two normalisation
(seven). Square. The normalisation rows make both `x` equal to one, the balance makes
`N_L` the feed, the two chemical-potential rows together say `μ_V = μ_L`, and the two
pressure rows say `P_V = P_L = P`: three equations of substance in `V_V, V_L, P`. That is
pure saturation, reached with **no reaction span, no charge row, no stability, no phase
count above two** — the active subset is seven rows of five registered kinds and no
charge or reaction rows, which is the evidence that a regression fit of pure-component
parameters depends on `equilibrium/` lightly (one compile, one assembly, one Ipopt call)
rather than structurally. **derived**; S14 asserts the count, the nullity and the absence
of reaction and charge rows, reports the row list (the multiset itself is not asserted,
because a legitimate implementation may drop a free phase's normalisation row), and
requires that pin exchange at `C = 1` (pin `L` instead of `V`) give the same `P*`, which
is the sharpest available test that "which phase is pinned" is a specification entry and
not a route. If the engine needs a distinct saturation solver, it is visible here without
solving: the compiled system for ethane would not be the general rows at `C = 1`.

The single-phase flash, `Π = 1`, `k = f = 0`, is the other cell a two-phase habit skips.
Unknowns `x (C), N, V, λ (C)`: `2C + 2`; rows: `C` chemical-potential, one pressure, one
normalisation, `C` balance: `2C + 2`. Square, and its solution is forced — `x` is the feed,
`N` is the feed amount, `V` is the density root of the pressure row, `λ` is `μ`. A single
phase is a density solve plus stationarity, and anything in a start, a distinctness test, a
K-value seed or a branch hint that indexes a second phase is a genericity defect of the
bubble-point family. With reactions the cell is not trivial (`x` moves in the reaction
span) and is S11's ideal-gas limit. S13 exercises `Π = 1` on every fixture. **derived**

**Where the layer line falls.** `model/` answers a state: `(T, V or ρ, n)` in, properties
and exact derivatives out; its one iteration is the density root, over one phase.
`equilibrium/` answers coexistence: anything where two states satisfy a shared condition.
Saturation needs two states, so it is equilibrium by the rule, not by judgement, and
nothing saturation-shaped may exist in `model/` — no seed, no vapour-pressure estimate, no
pure-component coexistence path. The Clausius–Clapeyron seed `TRANSFER.md` Stage 5 names
lives in `phase_boundary.cpp`. The scientific reason outranks the tidiness one: a pure fit
and a mixture fit then share one convergence behaviour and one tolerance entry, so a
parameter set cannot converge cleanly at `C = 1` under a private solver and fail at
`C = 2` under the general one with no signal. S5 and S9 running through the same check
function at `C = 1` and `C = 2`, against the same policy entries, is that claim as a test.
`regression/` including `equilibrium/problem.hpp` is already what checker rule 3 permits
and rule 4 bounds; the owner's discomfort with the dependency is answered by its size, and
S14's compiled row list for the one-component cell is where the size is measured.

## 3. The space

### 3.1 Axes

| axis | values | what it changes |
|---|---|---|
| phase set, by identity | `{L}`, `{V}`, `{V,L}`, `{L,L}`, `{V,L,L}`, `{L,L,L}` | a phase's identity is (support, branch hint, amount status). Three-liquid systems exist with ePC-SAFT (Ascani et al. 2022, water + ionic liquid + alkane) and cost nothing at the compile tier. `V` versus `L` is a **start** attribute — the density-branch hint for the initial `η_p` — and changes no row. `{L,L}` is in the space precisely because it is where a vapour-like reference root would be assumed if any part of the start or the seed assumed one |
| amount status per phase | free; pinned `β = 0`; pinned `0 < β < 1` | `k`, and the balance's view of the phase; no row is added or removed |
| `T`, `P` | each fixed or free | `f` |
| reactions | none; one; two independent; two plus a dependent third; one that does not conserve charge | `rank B`, `g°`, the rejection |
| ions | absent; present in one phase's support with both signs; present in two liquids (`π_ch = 2`: one charge row, one interphase potential); an incipient charged phase beside a charged finite one (the row goes to the pinned phase); a charged feed; a support with one ion sign; an incipient neutral phase beside a charged finite one | the charge rows and which phase carries none, the rejection |
| feed basis | mole fractions; amounts; mass fractions; molalities | **nothing in the compiled system** — that is the check |
| topology | declared; omitted | declared is the main path; omitted is rejected by name until Stage 7, then accepted with `k = f = 0` |
| feed coordinate free with a target row | (later; the old tree's "affine feed boundary") | the same count rule with a target row on the right: a freed feed coordinate is paid for by one target equation. Recorded so the rule is seen to extend; not built before Stage 8 needs it |

Solids and critical points are outside the space, as `formulation.md` §7 says: different
state variables, and second-order conditions as equations.

**Two counts that this document never conflates.** *Binary, ternary, quaternary* count
components (`C`); *`{L,L}`, `{V,L,L}`* count phases and name their identities. "The
immiscible ternary" is a three-*component* fixture used for a two-*phase* cell. Many
components with two phases is the ordinary case; the component count adds columns and
indices, never rows of a new kind, and the checks that see index errors (C7, S1, S8)
already run at nine components.

### 3.2 Unreachable cells, and what each must say

These carry the weight. A compile step that accepts them silently is the old tree's
"it never defaults" claim without its test.

| cell | rejection names |
|---|---|
| `f < k` | the pinned phase and "free `T` or `P`" |
| `f > k` | the free intensive and "fix it or pin a phase amount" |
| every phase pinned (`k = Π`) | "no free phase carries the feed" — square by the count, singular by rank; this cell is why the rank check is not redundant |
| `Π − k > rank B` | the counts, and "pin a phase or free an intensive" — the phase rule in rank form (§2.3); the pure-component two-phase flash and the binary three-phase flash at fixed `T`, `P` land here |
| dependent reaction set | the reaction that lies in the span of the others |
| `z·ν ≠ 0` | the reaction and "does not conserve charge" |
| `z·n_feed ≠ 0` | "feed carries net charge" |
| a charged support with one sign | the phase and "support cannot be electroneutral" |
| a species in no support | the species |
| an empty support | the phase |
| topology omitted, before Stage 7 | "phase count unknown; declare phases" |
| topology omitted with `k ≥ 1` or `f ≥ 1`, at Stage 7 | "discovery is a fixed-`T`,`P` problem" |

One cell is not a rejection but a decision the compile step must make and say: a species
with zero feed that no reaction produces has zero amount in every phase, and `ln 0` is not
a coordinate. Compile removes it from every support and reports the removal.

### 3.3 The compile-tier lattice is exhaustive, not sampled

The core axes multiply to several thousand specifications (6 phase sets × ≤ 8
amount-status patterns — free or pinned-at-zero per phase, plus one `β > 0` pattern per
two-phase set — × 4 intensive patterns × 5 reaction values × 6 ion values, about 5,800
before the generator discards combinations a fixture cannot express); each compiles in
microseconds
and the rank of a matrix under 40 × 40 costs the same. Feed basis is an invariance applied
to every cell (four compiles, compared), not a fifth axis. Nested loops over four enums
are the generator; there is no property-based framework and no random specification — the
space is small enough to walk. The generator picks the fixture from the (reaction, ion)
value — the binary for neither, the reactive ternary of §4.2 for reactions, the
electrolyte for ions, the nine-species set for both — and assigns supports by one rule:
liquids carry the full support, a vapour carries the neutral species. **guessed** on the count; **derived**
on the cost.

The point at which the Jacobian rank is taken is fixed, not random, and it must give
**every phase a different composition**: `x_{p,i} ∝ 1 + 0.37 i + 0.11 p`,
`η_p = 0.25 + 0.05 p`. The reason is not aesthetic. When two phases share a composition,
the columns of `N_1` and `N_2` in the balance rows are identical and vanish in every other
row (the chemical-potential and pressure rows are intensive), so the Jacobian has nullity
at least one *at every such point* — this is the trivial solution of the flash, which is a
genuine stationary point of every multiphase cell, and a rank check taken there would
reject every square cell. The rule concerns free phases: at `C = 1` every phase has
`x = 1`, and the boundary cell is fine because only one phase has an amount column.
**derived**

### 3.4 Declared topology: what it buys, and what it must refuse

A declared phase set is a specification, and it is what lets a solve skip stability
entirely: the rows are the equalities of §2.2 for every declared phase, including the
normalisation `Σ x_p = 1` as an *equality* for each. That is a claim about the answer —
every declared phase is present and distinct — and the engine either delivers exactly that
or delivers nothing. The rules:

- A solution is accepted only if every declared `N_p` is above the composition floor,
  every pair of phases is distinct (composition distance or packing fraction above the
  category-3 entries Stage 4 already carries from the paper), and the reduced Hessian is
  positive definite.
- Otherwise the call **fails with no equilibrium state**: no amounts, compositions,
  volumes or `G` are returned. What is returned is the failure name, the phase, and the
  quantity that condemned it. For a vanished phase that quantity is Michelsen's criterion
  (1994, p. 545: an absent phase satisfies `Σ_i Y_{i,m} < 1` where a present one satisfies
  `= 1`), written on these rows as `Σ_i x_{p,i} exp(−r_{p,i})` (the rows are dimensionless, so no `RT`) with `r_{p,i}` the
  phase's chemical-potential row residual `μ_{p,i} − λ·B_i − z_i φ_p` at the converged
  multipliers; the deficit `1 − Σ` is the report, and a positive deficit says the phase is
  not stable to add. For two phases that collapsed onto each other it is their distance,
  `d(p, q) = max( max_i |x_{p,i} − x_{q,i}|, |ln η_p − ln η_q| )`, and the refusal is
  `d ≤ PHASE_COINCIDENCE_DISTANCE` (`policy`, `AlgorithmTolerance`, `10⁻⁶`). The scale is
  derived, not chosen: at distance `d` the two phases' amount columns in the balance rows
  differ by `O(d)`, so the KKT matrix's conditioning is `O(1/d)`, and below `10⁻⁶` a row
  residual of `EQUILIBRIUM_ROW_RESIDUAL` cannot separate the pair in double precision. The
  test runs on **every multiphase cell, `k = 0` and `k ≥ 1` alike, over every pair, pinned
  or free**. A pinned phase that has slid onto the free phase's branch is the trivial
  solution of the boundary cell exactly as a coincident pair is the trivial solution of the
  flash; attempt 3 guarded the test with `pinned_count == 0` and returned that solution as
  converged at 15.8 bar, 149 Pa and 3.8 Pa for the same problem. This is the owner's rule
  stated as a mechanism: the trivial
  solution is never an answer, and neither is a two-phase result with one phase empty.
- Declared topology does not make the solution unique. A declared set may admit several
  solutions (Coatléven and Michel 2022 §2.3: "the choice of present phases leading to the
  existence of a solution is not necessarily unique"; Othmer's up to three chemical-
  equilibrium solutions in a strongly non-ideal reactive ternary, cited in Ascani et al.
  2023 §2.1). The reduced Hessian certifies a local minimum; when the driver runs from more
  than one start the lowest `G` is taken; a global claim is never made (§8).
- The start for a declared multi-liquid set must split composition (§3.3), by the three
  rules of §5 step 4 and no other: without discovery the split is the caller's
  `composition_guess`, and a same-kind pair with no guess is refused `start_symmetric`.
  For a reactive set each start phase is first brought to chemical equilibrium along its
  own stoichiometric line, which is Ascani et al. 2023's step 3 and costs one single-phase
  reactive solve per phase.

Discovery (Stage 7) is the complement: phases omitted, Michelsen's inequality replaces
the equality, and the phase count is an output. The two are one row system with the
normalisation row as equality or as complementarity; nothing else differs.

### 3.5 The space the engine is built to solve, cell by cell

The owner's enumeration — neutral or electrolyte, non-reactive or reactive, any component
count, one phase of either kind, two phases, three if it is cheap — against what covers
each cell and at which tier. *Solve* means a solve-tier check or an oracle gate exercises
the cell; *compile* means the rows, rank and rejections are checked and nothing is solved. `C` runs
from one to nine across the table; no cell has a component-count branch to exercise, and
C7 and S8 at nine components are the checks that see an index used where a component was
meant.

| system | `Π = 1` (`{L}` or `{V}`) | `Π = 2`, `{V,L}` | `Π = 2`, `{L,L}` | `Π = 3` (`{V,L,L}`, `{L,L,L}`) |
|---|---|---|---|---|
| neutral, non-reactive | solve: `C = 1` ethane, S13; `C = 2` binary, S13 | solve: `C = 1` saturation, S5 S9 S14 and the NIST oracle; `C = 2` flash, bubble, dew, fixed-fraction flash, S3–S5, and the Cai 1996 oracle on the associating binary | solve: `C = 3` the reactive ternary with its reaction omitted, S6 S7; cloud point; `C = 4` Schmitt and Hasse 2005 Table 13 quaternary tie lines at 298.15 K (rows to transcribe) | compile tier: the fixed-`T`,`P` three-phase flash is rejected by the phase rule wherever `Π − k > rank B` (the binary), and accepted on the three-component fixture; the `k = 2` three-phase point is square on the binary; opportunistic solve on the ternary |
| neutral, reactive | solve: `C = 3` isomerisation as `{V}`, S11 | solve: a reacting liquid under a pinned vapour, on the reactive ternary | solve: `C = 3` Ascani 2023 fixture, S16; `C = 4` Senina 2022 Table 3 reactive tie lines judged against Ascani 2023's predicted tie lines, catalyst caveat (rows to transcribe); Schmitt and Hasse 2006 supplies `K(T)` | compile tier; the reacting ternary's fixed-`T`,`P` three-phase flash is rejected (`rank B = 2`) |
| electrolyte, non-reactive | solve: `C = 3` aqueous NaCl, S13; the Hamer–Wu activity oracle is single-phase states on the EOS surface, not an equilibrium cell | solve: `C = 3` aqueous NaCl, neutral vapour pinned, `P` free | solve: `C = 5` Khudaida 2026 Tables 3–4 tie lines under declared `{L,L}`: `π_ch = 2`, one charge row, one interphase potential, within the paper's Table 9 AAD; then Ascani et al. 2022 case study 2 from its refs 82–84; both under the Stage 4 declared-variant precondition | compile tier; Ascani et al. 2022 case study 1, three-liquid ionic-liquid systems, is the only published `Π = 3` electrolyte example in the library |
| electrolyte, reactive | solve: `C = 9` MEA speciation oracle | solve: `C = 9` MEA bubble point, neutral vapour pinned, at the composite gate | **compile tier only**: the nine-species parameters under a declared `{L,L}` exercise the rows with one charge row; the system itself has one liquid. No published measured case: Tsanas et al. 2022's demixing-amine example (water, DEEA, MAPA, CO2, two ion-containing liquids) is eNRTL with no ePC-SAFT parameters. Solve tier is filled by **engine-internal** evidence: the calibrated HBTA/TOPO reactive network in Lithium_Extraction (`analyses/hbta_topo_calibrated_surrogate`, engine commit 90688d1, seven conditions with their consequence norms), a reproduction target labelled internal rather than published | compile tier |

Two cells carry compile-tier evidence only, and they are named so that coverage is not
mistaken for proof: reactive electrolyte two-liquid, and every three-phase cell except
the opportunistic neutral one. The structure supports them because they add no row kind;
what they lack is a system known to exhibit them under this model.

## 4. What each check asserts, and the mutation that makes it red

`CANON.md` binds hardest here. A synthetic system has no measured answer, so no check
below asserts a value. Each asserts a property that holds for every correct implementation
of the equations, and each carries the mutation the brief must apply at acceptance.

Every quantity asserted is coordinate-free: **nullity** (unknowns minus rank), never the
NLP's variable count — a reimplementation that eliminates the density inline has fewer
variables and the same nullity. Row *labels* by `EqID` are asserted only where the
formulation registers the row.

A check binary records every failed assertion and continues; each named check prints one
line, `check_s3 pass` or `check_s3 FAIL <file:line: condition>`, and `main` exits nonzero if
any failed. Every solve inside a check runs under `core::Capture`, and a failing check prints
its trace: spans, the per-iteration events of §4.4, and the record. Two reds that hide ten
unknowns are what attempt 3's abort-first macro produced.

### 4.1 Compile tier — no solve; lands with Stage 5

| # | claim | asserts | mutation (must go red) |
|---|---|---|---|
| C1 | count rule | for every lattice cell, compile's verdict (square / over / under) equals `f − k` computed in the test from the specification alone | in the compile step, count a pinned phase as free |
| C2 | reject by name | every unreachable cell throws, and the message contains the offending phase / reaction / species name from §3.2 | swap any two messages |
| C3 | rank | for every cell that is square by the count, `F_y` at the generic point has nullity 0 **iff** `Π − k ≤ rank B` and not every phase is pinned; otherwise nullity ≥ 1 and the rejection names the condition | drop the normalisation row of pinned phases |
| C4 | feed-basis invariance | the same feed as mole fractions, amounts, mass fractions and molalities compiles to the same row list and an identical feed *composition* `n_feed / Σ n_feed` to 1e-12 — the total implied by each basis differs (one mole, one kilogram, one kilogram of solvent), so `B n_feed` differs by that total, and S8 is what makes it irrelevant | molality conversion omits the solvent mass |
| C5 | reaction-basis invariance | reactions `{ν₁, ν₂}` and `{ν₁ + ν₂, ν₁ − ν₂}` (and any invertible mixing) give the same row-space projector `BᵀB` to 1e-12, the same nullity, and the same reaction affinities `Σ_i ν_i g°_i` — not the same `B`, which is a basis, nor the same per-species `g°`, which depends on the reference-species choice | build `B` from the reaction list in order instead of from its span |
| C6 | pin identity | the boundary specification's residual vector at a point equals the flash specification's residual vector at the same point with `N_q = 0`, on every shared row — the owner's "same rows, one entry differs", made numerical | the balance keeps the pinned phase at its start amount |
| C7 | permutation | relabelling components, or reordering phases, permutes the rows and leaves nullity and residual norms unchanged to 1e-14 | index a phase where a component was meant in any row |
| C8 | charge rows | an electrolyte flash with one ion-containing phase compiles with nullity 0 and no charge row; with two, nullity 0 and one row; with the charged phase pinned, the row is on the pinned phase. This catches both the double count in `formulation.md` (§2's chart *and* §4's row) and the one-row-per-charged-phase error of §2.2 | write a charge row, with its multiplier, for every ion-containing phase — the constraint Jacobian loses row rank and the KKT matrix is singular |
| C9 | one phase is a member | `{L}` alone and `{V}` alone, at every component count the fixtures span and with and without reactions and ions, compile square with nullity 0 and `k = f = 0`; a single pinned phase is rejected as `k = Π`. Structurally implied by C1 and C3, named here so that a `Π ≥ 2` assumption anywhere in compile has a check with its own mutation | compile requires at least two declared phases |

C6 is the check the brief asked for. It has a limit worth stating: if `phase_boundary.cpp`
assembles rows through the same function as the minimizer, no mutation of that function can
make C6 red, because both sides move together. C6's mutation therefore lives in the
specification handling, and the checks that catch a *bypass* of compile are S3–S5 below.

### 4.2 Solve tier — the minimal input set

Parameter sets: five already exist in `engine/tests/native/model_support.hpp`
(**measured**: `ethane()` :42, `neutral_binary()` :51, `associating_binary()` :60 — the
Baygi 2015 MEA–water model the Cai 1996 oracle needs — `aqueous_electrolyte()` :87,
`nine_species()` :116); one is published; two arrive with oracle packets (Khudaida 2026,
and the esterification systems). The published one and the packets come from the
literature, not from the old tree:

- **The reactive two-liquid fixture is published.** Ascani et al. 2023 (*Molecules* 28,
  1768, §3.2 Table 3) define a hypothetical `A + B ⇌ C` with PC-SAFT parameters
  (`m = 2.40, 1.08, 2.80`; `σ = 3.20, 3.00, 3.80 Å`; `ε/k = 200, 400, 280 K`; `A` and `B`
  two-site associating with `ε^{AB}/k = 2500 K`, `κ = 0.05`; `k_AB = −0.045`,
  `k_BC = −0.025`, `k_AC = 0.045`) and show (their Figure 3) that at `K_a = 2.25` the
  chemical-equilibrium curve is cut by the miscibility gap into two disconnected branches
  with exactly two reactive tie lines. Existence is therefore *documented*, not
  self-certified, and the fixture is a two-phase, three-component, one-reaction cell built
  for the purpose by people who solved it. It replaces the old tree's immiscible ternary
  for every reactive cell.
- **The non-reactive two-liquid fixture** is the same three components with the reaction
  omitted, which keeps one parameter set for two axes. If its gap at the chosen state is
  inconvenient, a symmetric binary with one large positive `k_ij` is the fallback; the
  existence certificate below decides.

Each non-reactive synthetic cell must *exist* at the chosen state, and no measurement
says so. The existence certificate is the engine's own, one-sided: a two-phase solve with
`G` below the one-phase `G` at the same feed proves the one-phase state unstable (S6). A
cell whose certificate fails is a fixture problem, recorded in `STAGE-LOG.md`, not a code
failure. The reactive cell needs no certificate; the paper is the certificate.

**A real reactive two-liquid oracle exists and is in the library.** Senina et al. 2022
(acetic acid + 1-pentanol + pentyl acetate + water, 318.15 K, 1.013 bar, nine quaternary
tie lines at chemical equilibrium) and Schmitt et al. 2006 (the hexyl acetate analogue,
293–403 K), with the PC-SAFT parameters and `K_a` Ascani et al. 2023 tabulate (Tables 1, 2,
4). That is a four-component, two-liquid, one-reaction system with measured tie lines: a
published-model reproduction in exactly the sense `TRANSFER.md` §"Oracle" gives Baygi
2015. Recommended as an oracle case for Stage 5 (declared `{L,L}`), which puts reactive
LLE under measured data before any discovery stage exists. The rows are to be transcribed
from the papers already attached in Zotero; nothing needs acquiring.

| id | fixture | cells |
|---|---|---|
| ethane | `ethane()` | `{L}` alone and `{V}` alone at fixed `T`, `P`; `{V,L}` with `V` pinned, `P` free; the same with `L` pinned (must give the same `P*`); then `T` free. The Stage 5 file list is the record that these ran with no `gibbs_minimizer.cpp`, no `phase_stability.cpp`, and `B = I` |
| binary | `neutral_binary()` | `{V,L}` flash inside the envelope; bubble; dew; `β`-flash; both intensives |
| reactive ternary | Ascani 2023 `A + B ⇌ C` | **required:** `{L,L}` with the reaction, from three feeds on one tie line's moiety segment (S16); `{V}` with an isomerisation of the same components (S11); `{L,L}` without the reaction (S6, S7); `{L,L}` with `L₂` pinned, `T` free, with and without the reaction; `{V,L}` with the reaction and `V` pinned, `P` free — a reacting liquid under a pinned vapour. **Opportunistic:** `{V,L,L}` without the reaction and the three-phase point, if the fixture admits them at some state; recorded either way |
| electrolyte | `aqueous_electrolyte()` | `{L}`; `{V,L}` with `V` support `{water}` and `V` pinned, `P` free — the MEA shape at three species |
| electrolyte two liquids | oracle `aqueous_alcohol_salt_lle` (Khudaida et al. 2026, JCED 71, 708; packet revision 6 carries parameters only) | declared `{L,L}` from the tie-line midpoint feed; `π_ch = 2`, one charge row, one interphase potential. Rows: Tables 3 and 4, full quaternary tie lines with NaCl in both phases at 5 and 10 wt % over 293.15–313.15 K (a transcription task). Acceptance: within the paper's own ePC-SAFT AAD per component and phase (Table 9: grand 0.016 at 5 wt %), not the row uncertainty. Parameters: Tables 5–7 (Born diameters 3.445 and 4.100 Å, `f_k`, `k_ij`). **Precondition:** the greenfield model supports the variant the paper declares (Born diameters and `f_k` per Figiel 2025), a Stage 4 declared-choice question to check; the row is out of scope until answered. Second oracle: Ascani, Sadowski and Held 2022 case study 2 (water + 1-butanol or 1-propanol + two salts, 298.15 K) once its cited tie-line sources (its refs 82–84) are acquired and transcribed, with the same precondition (mixed dielectric by solvent mass fraction and ion mole fraction, ion `ε_r = 8`); its case study 1 is the only published `Π = 3` electrolyte example in the library. Withdrawn as oracles: Hubach et al. 2024 and Yu et al. 2024, extraction efficiencies only on switch-laden parameterisations. Discovery later reaches the same rows with topology omitted |
| reactive two liquids, real | Senina et al. 2022 Table 3 (acetic acid + n-amyl alcohol + n-amyl acetate + water, 318.15 K, nine reactive tie lines, `u(x) = 0.002`); parameters and `K_a` from Ascani 2023 Tables 1, 2, 4 | declared `{L,L}` on the nine tie lines; judged against Ascani 2023's **predicted** tie lines (Figure 9, right) rather than the data, because the paper states its tie-line prediction is qualitative and worse with acetic acid (dimerisation not represented). Caveat recorded: an acid catalyst is present in the data (Ascani 2023 §3.1: HCl(aq) below 2 wt %), which no model in the library represents. Schmitt and Hasse 2006 is the `K(T)` source for the hexyl system; Schmitt and Hasse 2005 Table 13 is a separate, non-reactive LLE at 25 °C. Two records once the rows are transcribed: Senina's measured tie lines as an oracle case in their own right, and the comparison against Ascani's predicted curve as reproduction evidence, never as an oracle |
| neutral two liquids, real | Schmitt and Hasse 2005 Table 13 (water + acetic acid + 1-hexanol + hexyl acetate, five quaternary tie lines at 298.15 K, both phases, no catalyst); parameters from Ascani 2023 Tables 1–2 | declared `{L,L}` on the five tie lines; the measured neutral quaternary oracle, acetic-acid caveat recorded. The engine-certified fixture stays the Ascani ternary with its reaction off. An oracle case (measured rows) once the rows are transcribed |
| reactive electrolyte two liquids, internal | the calibrated HBTA/TOPO reactive network in Lithium_Extraction (`analyses/hbta_topo_calibrated_surrogate`, engine commit 90688d1) | declared `{L,L}` with a cross-phase reaction and one ion-bearing phase at the seven evaluated conditions; **engine-internal evidence, not published**: a reproduction target with its consequence norms (balance, affinity, KKT and reduced-Hessian norms in `condition_results.csv`). The compile-tier nine-species exercise stays as is |
| nine species | `nine_species()` | compile tier: the composite cell's rows, rank, nullity, and reject-by-name with the bundle's two reactions; and its parameters under a declared `{L,L}` compiled with one charge row, which exercises the row assembly for that shape and says nothing about the MEA system, whose liquid is one aqueous phase. Its solves are oracle gates: speciation as `{L}`, and the bubble point at the composite gate |

Transcription tasks, each with its locator: Khudaida 2026 Tables 3–4; Senina 2022 Table 3;
Schmitt and Hasse 2005 Table 13; Ascani et al. 2022 refs 82–84 after acquisition. Yu et al.
2025 is not in the library and is not cited.

Every synthetic state below is certified by S6 (`G₂ < G₁ − limit`) **before** any cell at that
state runs; a failed certificate is a fixture finding for `STAGE-LOG.md`, and the cells that
depend on it are reported skipped by name, neither green nor aborted. Pressures in parentheses
are Raoult estimates for placing the state; the engine's own `P*` is the record. Attempt 3
copied fixtures without states: `associating_binary()` at 280 K and 2.5 kPa is subcooled
liquid, and `neutral_binary()` at 280 K and 50/50 is above the mixture critical point.

| fixture | file | state | cells |
|---|---|---|---|
| `neutral_ternary()` — methane, ethane, propane, Gross–Sadowski 2001 pure parameters (`m = 1.0000, 1.6069, 2.0020`; `σ = 3.7039, 3.5206, 3.6184`; `ε/k = 150.03, 191.42, 208.11`), every `k_ij = 0` (methane–ethane is 0 by Gross's statement; the other two pairs are absent from the packet, so the fixture is synthetic, not an oracle) | `model_support.hpp` | `T = 180 K`, feed `1/3` each; all three subcritical (trap 13) | bubble, `P` free (~11 bar); dew from `y*`, cold; `{V,L}` flash at 1 bar (S2, S6 certificate, S8, S10 flash half, S13 as `{L}` and `{V}`); `β`-flash (S4); `T` free at `P*` (S5); `{V,L}` at 30 bar for S15 (subcooled) |
| `neutral_binary()` — methane, ethane | exists | `T = 180 K`, feed `0.5/0.5` (bubble ~16 bar, dew ~0.16 bar) | S3 dew ↔ bubble (cold, trap 13); the `C = 2` compile cells |
| `ethane()` | exists | `T = 250 K`, `P` free | S14 |
| Ascani ternary, reaction omitted | `equilibrium_support.hpp` | `T = 298.15 K`, `P = 1 bar`, feed on the Figure 3 tie line, `composition_guess` from the two branches | `{L,L}` (S6 certificate, S7); `{L,L}` with `L₂` pinned, `T` free |
| Ascani ternary, `A + B ⇌ C`, `K_a = 2.25` | exists | same `T`, `P`, three feeds on the moiety segment | S16; `{V}` isomerisation (S11); `{V,L}` reacting under a pinned vapour |
| `aqueous_electrolyte()` | exists | `T = 298.15 K`, `P = 1 bar`, 1 molal | `{L}` (S1, S13); `{V,L}` with `V` support `{water}`, `V` pinned, `P` free |
| MEA + water, `associating_binary()` | exists | **oracle**: Cai et al. 1996 isobaric VLE at 101.3 kPa, `T` free (bubble `T` at the measured `x`), rows transcribed into `oracle.json` | published-model reproduction; never a synthetic state |
| Khudaida 2026 electrolyte two liquids | packet revision 6 | 293.15–313.15 K, 0.1 MPa, 5 and 10 wt % NaCl, tie-line midpoint feeds; the observation rows (Tables 3–4) are **not in the packet** and must be transcribed into `oracle.json`; acceptance is the paper's Table 9 AAD | S-tier oracle for `{L,L}` with one charge row, under the Stage 4 declared-variant precondition |

| # | claim | asserts | mutation |
|---|---|---|---|
| S1 | rows satisfied | every solved cell's residual vector is within `max(tolerance, floor)` of zero on every row, the two numbers read from the result's own `Convergence` record (`CORE.md` §2, consequence 3) — including rows the driver did not iterate on (affinity for a reactive cell; Euler's `g_p = λ·Bx_p` for every phase; `G_min = λ·B n_feed`; and `|z·x_p|` for **every** ion-containing phase, whether or not it carried a row) | drop one balance row from the Ipopt constraint set: affinity goes nonzero |
| S2 | local minimum | for `k = 0`, the reduced Hessian `Zᵀ ∇²L Z` is positive definite (`EqID: reduced_hessian_stability`) | negate one phase's Hessian block in the assembly |
| S3 | dew ↔ bubble | for the binary at `T`: bubble from liquid `x` gives `(P*, y*)`; dew from `y*` must give `(P*, x)` to the convergence tolerance. An exact symmetry of the formulation, and it fails on any asymmetry between "which phase is pinned" | start the pinned phase on the vapour branch regardless of its declared identity |
| S4 | `β → 0` | the fixed-fraction flash at `β = 10⁻², 10⁻⁴, 10⁻⁶` has `|T(β) − T_bubble|` decreasing and bounded by a constant times `β` (each step of a hundred in `β` shrinks the gap by at least ten) | the balance omits the pinned phase for `β > 0` |
| S5 | round trip | solve `P` at fixed `T`; solve `T` at that `P`; recover `T` | freeze `T` in the Jacobian when `T` is the unknown |
| S6 | instability certificate | `G₂ < G₁` at the ternary feed; `G₃ ≤ G₂ + tol` for the three-phase cell | return the one-phase start as the result |
| S7 | lever rule | two feeds on the `{L,L}` tie line give the same two phases to 1e-8. This is S1 plus isolation of the solution, so it has no mutation of its own: S1's reddens it | shares S1's |
| S8 | extensivity | every solved cell with the feed doubled: amounts double, every intensive unchanged | compile normalises `n_feed` to unit total before the balance rows |
| S9 | sensitivity | `dy/dθ` from `implicit_sensitivity.cpp`, `θ` one of three levers: a fixed intensive (`T` or `P`) or one species-feed coordinate `n_feed,i` that carries a moiety (`B[:, i] ≠ 0`) (the model's active parameters are the same call with more columns and land with regression), against a centred difference of two full solves in `θ` to 1e-6 relative; and, on the reactive fixture, the derivative along a reaction direction `d` with `B d = 0`, formed as `Σ_i d_i · dy/dn_feed,i`, is exactly zero. `F_θ` is the assembly's: the tape's Hessian column in `T` plus the standard term's own `∂/∂T`, `−v_p/(RT)` on the pressure row (written as `−c_p/(η_p RT)` with `c_p` the packing coefficient, the derivative of `Z_p − P v_p/(RT)` in `P`), and on the balance rows `−B[:, i] / Σ_j n_feed,j + Σ_{q pinned} β_q (B x_q)_r / N_feed` for the feed; the test supplies no derivative. Cells: the ternary flash (`T`, `P`, one feed column), the fixed-fraction flash at `β = 10⁻²` with the feed column (the pinned-amount term), a reactive free-`T` cell with `±1,000 J/mol` standard-state offsets (the standard term's derivative), and an active-parameter column on the binary flash against two full solves with the parameter perturbed | flip the sign of `F_θ` |
| S10 | Jacobian agreement | the assembly's `F_y`, analytic from the tape's Hessian and the log-coordinate chain rule, against a centred difference of `F` in every unknown, at a boundary cell and at a flash cell, to 1e-6; and Ipopt's `eval_h` at the flash solution with the converged multipliers against a centred difference of the Lagrangian gradient (`eval_grad_f + Σ λ_j ∇c_j`) to 1e-6. Nothing about call counts. The drivers' rows differ by `N_p` per phase (§2.2), so no elementwise equality between consumers is asserted; "one assembly" is held by §5's placement and rule 7, not by a check | scale the pressure term of `evaluate`'s pressure row by 2: an `F_y` that is a difference of `evaluate` follows it and stays green, and a green mutant rejects the attempt; and scale the tape's Hessian blocks by 0.5 for the `eval_h` half |
| S11 | ideal-gas limit | the reactive one-phase cell **as `{V}`** — the vapour branch is the one that reaches the ideal-gas limit; a liquid hint stays on its own branch — with the isomerisation `α ⇌ β` (so `Σ ν_i = 0` and the `ρ_ref` term cancels) at `P → 0` approaches `x_β / x_α = exp(−Δg°/RT)` (a limit; the residual terms vanish, which Stage 4 already checks). A reaction with `Σ ν ≠ 0` has a density factor in its limit and is not used here | perturb the `g°` transform |
| S12 | discovery agrees (with the stability kernel) | topology omitted at the ternary feed finds the phase count of the declared `{L,L}` flash and the same `G` | accept the first sample as a witness |
| S13 | single phase | `Π = 1` on every non-reactive fixture at fixed `T`, `P`: the solved composition equals the feed within the balance rows' `max(tolerance, floor)` (the balance is an Ipopt equality, so 1e-14 is not reachable), and the solved volume equals `model/`'s `solve_density` root at the same `(T, P, x)` within the larger of the two solves' `max(tolerance, floor)` — the pressure row *is* the density equation when there is one phase. Reactive single phases are S11's | seed the start from phase index 1 |
| S14 | `C = 1` is the general row set | the ethane boundary compiles with `rows − unknowns = 0` and nullity 0 with no reaction span and no charge row, and pinning `L` instead of `V` gives the same `P*` to the convergence tolerance. The row list by `EqID` is *reported*, not asserted: a legitimate implementation may drop a free phase's normalisation row (§2.2), and `CANON.md`'s diagnostic forbids asserting the multiset | a `component_count == 1` branch in compile that omits the balance row |
| S15 | a declared phase that vanishes | the binary `{V,L}` flash at a state outside its envelope, and the `{L,L}` flash at a feed outside the gap, **fail with no equilibrium state** (§3.4): the failure names the phase and carries Michelsen's deficit `1 − Σ` or the phase distance; no amounts, compositions or `G` are returned, and a symmetric start on any two-phase cell is refused the same way; and (c) the binary `{V,L}` boundary cell with `V` pinned and `P` free at 280 K, feed 0.5/0.5, a state above the mixture critical temperature where only the trivial solution exists (attempt 3 returned it as converged at 3.8 Pa), fails with no state and its failure names `phases_coincide` with the distance; it lives in `boundary_checks.cpp`. The start rule re-seeds a pinned phase on its declared branch, so no caller can place it on the free phase's branch by hand; the trivial solution is reached physically | (a, b) return the converged iterate when one `N_p` sits at its bound; (c) skip the coincidence test when `k ≥ 1` |
| S16 | reactive lever rule | on the `A + B ⇌ C` fixture at fixed `T`, `P`, the declared `{L,L}` reactive flash from three feeds whose moiety vectors `B n_feed` lie on the segment between the two phases' moiety vectors of a first solve returns the same two phase compositions to the convergence tolerance, and affinity vanishes in both phases with one `λ`. Gibbs' rule for the reactive system gives `F = 2 − 2 + 3 − 1 = 2`, both consumed by fixing `T` and `P`, so reactive tie lines are *isolated*, not unique — Ascani et al. 2023's Figure 3 has two, and "unique reactive tie-line" is their term of art for one that appears in a ternary — which is why the three feeds are placed on one tie line rather than anywhere in the gap. A route that solved the reaction in one phase and partitioned into the other would move the phases with the feed | drop the affinity consequence by solving the balance with `B = I` |

**Reported, never asserted**, per cell: Ipopt iterations, evaluations,
`hessian_calls`, recordings, wall. `hessian_calls ≥ 1` per iteration on every cell is how
"the exact Hessian is consumed everywhere" is *seen*; it is not a pass condition. One
architecture rule is cheap enough to add to the checker: `gibbs_minimizer.cpp` defines
`eval_h` and contains no `limited-memory` string. Lexical and weak, it tests the driver
rather than the physics, and it says what the design intends.

### 4.3 What the tiers establish together

C1–C8 prove that the compile step is one rule over the whole space and that the rejections
work. S1, S2, S8 prove each solved point is what it claims. S3–S5 are the checks that see a
*bypass*: a second bubble-point path that never enters compile still has to be symmetric
under pin exchange, continuous in `β`, and invertible in `T`/`P`, and a route written for
one application will usually fail one of them. S10 is the only check on the derivative
path itself, and it is a derivative-agreement check of the kind `CANON.md` already
admits; that the two drivers share one assembly is enforced by where the assembly lives
(§5) and by rule 7 covering that file, not by a test. The artifact's "seven genericity equivalences" (round trip, reaction basis,
extensivity, a problem stated two ways) are not seven checks here; they are C4, C5, S5, S8
applied to every cell they fit, which is what makes them coverage rather than spot checks.

### 4.4 Diagnostics the brief carries

The `CORE.md` §7 block, one entry per solver or assembly this stage writes. The brief's
`diagnostics` array is this table verbatim, and the preflight refuses a brief without it.

| span | events | convergence | check |
|---|---|---|---|
| `solve_equilibrium` | `ipopt_iteration` per iteration from `intermediate_callback`: `inf_pr, inf_du, mu, d_norm, regularization_size, alpha_pr, alpha_du, ls_trials`; `phase_distance` per pair; `michelsen_deficit` per phase; `reduced_hessian_minimum`; `start_rule` (1, 2 or 3 of §5 step 4) | residual: infinity norm of the dimensionless rows at the returned point; criterion `EQUILIBRIUM_ROW_RESIDUAL`; floor: largest phase `gradient_floor`; amplification: Ipopt's last `regularization_size` | `record.residual` equals `max |row|` recomputed from the returned state through `CompiledSystem::evaluate` |
| `equilibrium_evaluate` | `row_values`, `gradient_floor` per phase | — | — |
| `boundary_start` | `k_factor` per component, `branch_density` per phase | — | — |

## 5. From specification to NLP

    Problem  ──compile──►  CompiledSystem  ──►  Ipopt: objective G (k = 0) or 0 (k ≥ 1),
                                │               constraints = the rows compile marks  ──►  EquilibriumResult
                                │                     F(y), F_y from one assembly
                                └──►  implicit sensitivity: −F_y⁻¹ F_θ

1. **`Problem`** (`problem.hpp`) is data: phases with support, branch hint and amount
   status; `T`, `P` each fixed or free; feed in a declared basis; reactions; active
   parameters. It has no method that solves.
2. **compile** (`request_compile.cpp`, Stage 5) resolves the feed basis to amounts
   (native — see §7, item 4; a molality feed needs the solvent set, which is a `Problem`
   field, not an inference), builds `B` from the reaction span by SVD, applies the count
   rule, assembles the active row list with its `EqID`s, takes the rank at the generic
   point of §3.3, and either returns a `CompiledSystem` or throws naming the deficiency.
   It is the one file allowed to branch, and only on status flags.
3. **assembly** — `CompiledSystem::evaluate(y)` giving `F(y)` and `F_y`, **also Stage 5,
   in `request_compile.cpp`**, not in `gibbs_minimizer.cpp`: the compile-tier checks
   C1–C9 and `implicit_sensitivity.cpp` consume it with no driver present, the driver is
   written *against* it, and a driver file that owned the rows would be the place a
   second row set grows, which is the failure this document exists to prevent. One tape evaluation per phase per iteration at unit amount (§2.2) — value,
   gradient, Hessian of `A_p` in `(T, v_p, x_p, θ)` — scaled by `N_p`, plus the linear
   constraint blocks. The Hessian of `A_p` is the tape's reduced Hessian; the site-fraction
   elimination happened in Stage 4 and nothing here sees sites. Because the assembly moves
   here, `equilibrium/request_compile.cpp` joins rule 7's `SOLVER_FILES`
   (`tools/check-engine-architecture.py:163`): it may branch on status flags, and it may
   not name a problem kind — the comment at `:154-156` is updated to say so.
4. **driver** (`gibbs_minimizer.cpp`, Stage 5): one Ipopt `TNLP` with bounds on `ln η_p`,
   `ln x_{p,i}` and `N_p ≥ 0`, `nlp_scaling_method = none` (`CORE.md` §3.4),
   objective and constraint set read from the `CompiledSystem` (§2.4), exact Hessian for
   `k = 0` and zero for `k ≥ 1`. Seeds and continuation for boundary cells live in
   `phase_boundary.cpp` (Stage 5): a `k ≥ 1` solve is seeded from the flash
   solution or the previous boundary point when the caller is sweeping, and there is never
   an outer loop around `solve_equilibrium` (`CORE.md` §6). Every multiphase start must give the phases *different*
   compositions, by one of three rules and no other:
   1. **The caller's `composition_guess`** for every phase that carries one.
   2. **A one-step K-factor split** for a pair whose declared kinds differ (`V` and `L`),
      from the model alone: evaluate `ln φ_i` at the feed composition on the vapour and on
      the liquid density branch (`state_at_pressure`, one call per branch), set
      `K_i = exp(ln φ_i^L − ln φ_i^V)`, `y_i ∝ z_i K_i`, `x_i ∝ z_i / K_i`, renormalise, and
      take each phase's `η` from its own branch root. No critical constants, no Wilson
      correlation, general in `C`.
   3. **Refusal**, named `start_symmetric` with the two phase names, for a same-kind pair
      (`L`, `L`) with no guess. The design does not invent a liquid–liquid split; the
      stability kernel supplies one when discovery exists, and until then the caller does.

   For `k ≥ 1` with a `start` given, compositions and `η` come from the start phase by
   phase; the pinned phase's `η` is re-seeded on its **declared** branch from the density
   root at the start's `T`, `P` and its own composition, never copied from the free phase.
   Amounts from a start are kept, not overwritten with an equal split (attempt 3's
   `start_point` overwrote them). The trivial solution is a stationary point (§3.3) and a
   start on the symmetric subspace, equal composition and equal `η`, converges to it; a
   composition-symmetric start with branch-distinct `η` is not on that subspace.
5. **result**: amounts, volumes, intensives, multipliers, every residual row, the reduced
   Hessian's smallest eigenvalue, counters, a named status. Sensitivity when asked. Or,
   under §3.4, no state at all: a failure name, the phase, and the number that condemned
   it.

The coverage exercises this path and not a bypass of it because C1–C8 act on step 2's
output, S1 reads step 5's residuals for rows step 4 did not iterate on, and S10 holds the
three consumers of step 3 to one matrix.

**Tolerance and budget.** Four policy entries, in `CORE.md` §3's categories:

| entry | value | category | rationale |
|---|---|---|---|
| `EQUILIBRIUM_ROW_RESIDUAL` | `1e-10` | AlgorithmTolerance | infinity norm of the dimensionless rows (§2.2); exact Newton on the analytic `F_y` reaches it in a handful of steps; two decades above the association floor `ASSOCIATION_RESIDUAL = 2e-13` that enters the record as `floor` |
| `EQUILIBRIUM_MAX_ITERATIONS` | `200` | Budget | rule 9's predicate: a budget makes the loop a solver; attempt 3 used Ipopt's default 3000 without a policy budget and rule 9 never saw it |
| `PHASE_COINCIDENCE_DISTANCE` | `1e-6` | AlgorithmTolerance | §3.4 |
| `PACKING_FRACTION_CEILING` | `0.74` | AlgorithmTolerance | §2.1 |

Ipopt reads `tol`, `constr_viol_tol`, `dual_inf_tol` from `EQUILIBRIUM_ROW_RESIDUAL`,
`max_iter` from the budget, `acceptable_iter = 0`. The record's `criterion` points at
`EQUILIBRIUM_ROW_RESIDUAL`; `floor` is the largest phase `gradient_floor`. Nothing under
`equilibrium/` names `DENSITY_ROOT_RESIDUAL`: a density Newton tolerance on rows of another
unit was attempt 3's `5e-12` chase.

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
coordinates); §4 lists `phase_electroneutrality` as a row "for every phase with charged
support". Both together is a redundant row, and one row per charged phase is *already* one
too many (§2.2: `π_ch − 1` rows, three sources). Recommendation: the row form, with its
multiplier `φ_p` as the interphase potential, `π_ch − 1` of them, and no elimination in
the minimizer — the charged chart is then the neutral chart plus at most one row per
additional ion-containing phase, which is smaller than "plus one elimination" and has no
`eliminated_index`. `phase-stability.md` §2 keeps the chart for the TPD search, where trial
compositions must be neutral by construction; `electroneutral_amount_chart` stays
registered for that owner.

**Standard states.** `stoichiometry.cpp`'s "least certain" part is the transform from
equilibrium constants to per-species `g°`. Tsanas et al. 2017 §2.4 gives the recipe:
choose `R` reference species each in at least one reaction, set the others' offsets to
zero, and solve `ν̂ᵀ g° = −RT ln K` for the `R` unknowns. The offsets are fictitious and
consistent, and every affinity is reproduced. Ascani et al. 2023 §3.3 fixes `K_a` from one
measured equilibrium composition through the model's own activity coefficients, which is
the only route when formation energies are absent, and is the route the esterification
oracle needs.

**Drop:** §4's sentence naming `compile_phase_direction_system` and the two Python rank
functions (deleted files); §9 entirely (owners that no longer exist).

**One new identifier, and a rule about labels.** `incipient_normalization` is registered
as "one per incipient phase"; §2.2 makes the row universal, so the every-phase row is
registered as `phase_normalization` and the old identifier is retired rather than
redefined. Identifiers are cited in `// EqID:` comments only, which the checker strips;
the row labels a `CompiledSystem` carries at runtime are an enum whose names contain no
rule-7 token (`normalization`, `balance`, `charge`, `pressure`, `potential`), because the
string `incipient` inside `request_compile.cpp` would otherwise fail rule 7 once that file
joins its list (§5).

**Registry:** land the 21 entries with empty `implementation_refs`; extend the sync script
so that a non-`Documentation-only` entry with no owner is a *warning* until the stage that
annotates it, and its `MATH_OWNER_FILES` lists `native/**` instead of the 21 deleted paths.
The `// EqID:` annotation system is worth keeping exactly as `CARRY-OVER.md` says; 58 of
them are already in `native/` (**measured**) and the script cannot currently read them.

## 7. Where the brief's framing is corrected

1. **"An incipient phase is not a phase of the minimizer"** (`CORE.md` §6). It is
   not an *Ipopt variable* of the minimizer, and the boundary is not a
   minimisation. It is a row of the same system, solved by the same driver on that
   system, not a second row set. If `phase_boundary.cpp` has its own chemical-potential
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
7. **Saturation is equilibrium, fully.** The owner's tension is real and the resolution
   holds under test: §2.5 counts the `C = 1` cell against the general rows with no special
   case and gets saturation, so a private pure-component solver would be a second
   implementation of seven rows that already exist. **measured** residue check: nothing
   saturation-shaped exists under `native/` today — no seed, no estimate, no coexistence
   path, no policy entry — and `infinite_dilution_reference.cpp` is a sequence of
   single-phase states, correctly placed. Two planned items were residue and
   `target-layout.md` is revised: the `Saturation` `py::class_` leaves `eos_bindings.cpp`
   (the call is `solve_equilibrium`, bound in `equilibrium_bindings.cpp`), and `eos.py`'s
   `Pure.saturation` is a thin call onto that binding, kept where users reach for it,
   because the Python surface is organised by capability and was never meant to mirror
   the layers.
8. **The test budget.** `engine/tests/native/` holds 1,188 lines after three closed stages
   (**measured** 2026-09-02), against a 1,650 total for the whole engine. This design adds
   about 640 (§9). With the oracle-row tests still to write, the engine lands near 2,400 lines of
   tests, not 1,650. `CANON.md`'s rule is that a stage over its share must show its
   checks pass the diagnostic, not that the budget grows; every check above states the
   property it tests and the mutation that fails it, and the number is recorded rather than
   trimmed. `model_floor.cpp` (186) is a measurement harness, not a test, and arguably sits
   outside the budget; that does not close the gap.

## 8. What completeness does not establish

- **Starts.** A cell solved from one start proves one basin. Continuation across a sweep,
  the shifted-previous start family, and the fixed start list of Stage 7 are not covered by
  any check above; a cell can pass every property and fail from the next caller's start.
  The trivial solution is the standing example: every multiphase cell has it, a symmetric
  start finds it, and only S6 and S15 notice.
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
- **Existence is self-certified for the non-reactive synthetic cells.** They "exist"
  because the engine's own `G₂ < G₁` says so. That is one-sided: it proves instability,
  not that the two phases found are the right two. The reactive fixture and the
  esterification oracle are exempt; a paper and a data set are the certificates.
- **Declared topology guarantees a local answer, not the answer.** Several declared
  solutions may exist (§3.4). A cell that returns one and passes S1, S2 and S16 has still
  not been shown to be the global minimum; only discovery can say that, and §6 of the
  formulation says even discovery records `global = None` unless the bound loop closes.
- **Tolerances define satisfaction.** "Satisfies its own rows" means within Ipopt's and the
  policy header's tolerances. A tolerance loosened by mistake weakens every S-check at
  once; S1's mutation is the guard, and it is one mutation.
- **The discovery cell asserts a rejection until Stage 7.** Between Stage 5 and Stage 7,
  "topology omitted" is proven only to fail by name.
- **A single-phase cell proves consistency, not stability.** S13 shows the rows and the
  density closure agree; whether that one phase is the stable state at `(T, P)` is Stage
  11's question, and until then a `{L}` result carries no stability claim.

## 9. Cost and placement

All native, under the existing `CHECK` harness (`engine/tests/native/physics_checks.cpp:18`);
no framework, no fixture library, no Python until Stage 3's bindings exist, which under the
revised order (`TRANSFER.md` §"Execution order") is before Stage 5 anyway. **guessed**:

| file | lines | ceiling | lands with | holds |
|---|---|---|---|---|
| `tests/native/equilibrium_support.hpp` | 140 | 190 | compile | includes `model_support.hpp`; adds the Ascani 2023 reactive ternary, its reaction and `g°`, specification builders, the four-enum lattice walk, the generic point, the count-rule predictor |
| `tests/native/specification_checks.cpp` | 180 | 240 | compile | C1–C9 |
| `tests/native/equilibrium_checks.cpp` | 170 | 220 | the driver and starts; +15 with discovery | S1, S2, S6, S7, S8, S10's `eval_h` half, S11, S13 (multicomponent), S15, S16 on the flash cells; S12 at Stage 7 |
| `tests/native/boundary_checks.cpp` | 140 | 180 | the one-component cell opens it (ethane: S5, S9, S10's `F_y` half, S13 at `C = 1`, S14, ~60); the multicomponent cells complete it | S3, S4, S5, S9, S10's `F_y` half, the `{L,L}` and `{V,L,L}` boundary cells |
| `tools/check-engine-architecture.py` | +2, then +12 | — | compile, then the driver | `request_compile.cpp` joins rule 7's file list with compile; rule 8 (`eval_h` present, no `limited-memory`) with the driver |
| **total** | **~640** | **~840** | | |

The "light dependency" claim of §2.5 is measured rather than argued by the one-component
cell itself: S14 records its compiled row list, which carries no reaction span and no
charge row, and its solve reaches no stability kernel. `STAGE-LOG.md` records the row list
the first time the cell runs.

Order inside Stage 5: `equilibrium_support.hpp` and C1–C3 first, before any solver exists —
compile and the count rule are the first thing the stage writes, and C1–C3 are green on
the compile step alone. That is also the cheapest possible preflight for the stage brief:
if C1's predictor and the compile step disagree on a cell, one of them misreads §2.3, and
that is found before a line of `gibbs_minimizer.cpp` is written.

## 10. Traps, one line each — the brief copies this list

Every item below is a place an implementer would otherwise guess, or a place where a check
as first written could only be passed by fudging. Found on a re-read of this design against
`CORE.md`, the model layer and the checker; the first four are corrections to this
document's own first draft.

1. Evaluate every phase's tape at unit amount and scale by `N_p` afterwards (§2.2), or a
   pinned phase is `0/0` and a second path grows.
2. The rank point must give every phase a different composition (§3.3), or every square
   multiphase cell is rejected at the trivial solution.
3. S10 compares against centred differences; the two drivers' rows differ by `N_p` per
   phase and are not elementwise equal.
4. C4 compares feed *composition*; the totals implied by the four bases differ.
5. The assembly `F`, `F_y` lives in `request_compile.cpp`, and that file joins rule 7's
   list; the driver and the sensitivity must find it there, never write their own.
6. Multiphase starts split composition by one of the three rules of §5 step 4 (the caller's
   guess, the K-factor split, or the `start_symmetric` refusal) and no other. The symmetric
   subspace is equal composition **and** equal `η`: distinct branch densities alone break it
   for a vapour–liquid pair, which is why rule 2 gives both a split and a branch per phase,
   and why a liquid–liquid pair with no guess is refused (two liquid branches share `η`).
   The trivial solution is a stationary point and S15 requires it reported by name.
7. S1 and S13 test against `max(tolerance, floor)` from the result's `Convergence` record,
   never against a bare policy number — `CORE.md` §2 says why a bare number can be
   unreachable.
8. S11 uses the isomerisation only: with `Σ ν ≠ 0` the ideal-gas limit carries a density
   factor and the closed form the check names is wrong.
9. A molality feed needs a declared solvent set on `Problem`; compile does not infer one.
10. A support restriction is a zero amount passed to the tape; `validate_point`
    (`helmholtz_model.cpp:250`) accepts zeros, and a pure-water vapour inside an
    electrolyte model takes the no-ion tape branch — Stage 4's predicate-signature cache
    key exists for exactly this, and the assembly must not bypass it.
11. `nlp_scaling_method = none`, bounds on `N_p ≥ 0`, and no outer Newton around
    `solve_equilibrium` — three commitments from `CORE.md` that S1, S15 and the
    Stage 5 gate inherit.
12. The `C = 1` seed needs a pressure at which both density roots exist; the density
    solver's `stable_root_count` says whether it does, and a seed where it does not is a
    named stop, not a retry loop.
13. S3's dew is seeded cold: `solve_equilibrium(model, dew)` with no `start`, or the
    symmetry is trivial; the state is chosen subcritical for both components so the dew
    pressure is unique. Attempt 3 passed the bubble result as the start and still failed,
    which is how the trivial solution was found.
14. Nothing asserts a variable count, an iteration count, or which driver ran. Counters
    are printed per cell; the design intent that the exact Hessian is consumed is visible
    in `hessian_calls`, and rule 8 is the only thing that turns it into a failure.
15. One electroneutrality row per ion-containing phase **except one phase of nonzero
    amount** (§2.2); one per charged phase is a dependent row and the rank check rejects
    valid electrolyte flashes, and leaving a *pinned* charged phase row-less lets it
    converge non-neutral.
16. A declared phase set is delivered whole or not at all (§3.4): a vanished phase, a
    collapsed pair or the trivial solution fails with no state, carrying the phase name and
    Michelsen's deficit or the distance.
17. The reactive two-liquid fixture is Ascani et al. 2023 Table 3 with `K_a = 2.25`, not
    the old tree's immiscible ternary; its existence is the paper's Figure 3.
18. At fixed `T`, `P` a one-reaction ternary has `F = 2`, both consumed, so reactive tie
    lines are isolated, not unique; S16 places its feeds on one tie line's moiety segment
    and asserts the same phases (the reactive lever rule). Feed dependence along that
    segment is a route solving the reaction in one phase.
19. Start phases of a declared reactive multi-liquid set are first brought to chemical
    equilibrium along their own stoichiometric lines (Ascani et al. 2023 step 3), one
    single-phase reactive solve each; the feed itself stays on the stoichiometric line and
    is *not* on the chemical-equilibrium surface after the split.
20. `Π − k ≤ rank B` (§2.3) is the phase rule in rank form; the pure-component two-phase
    flash and the binary three-phase flash at fixed `T`, `P` are square by the count and
    singular at every point, and compile rejects them by name.
21. No identifier under `native/` is named after saturation: `solve_saturation` in
    `phase_boundary.cpp` fails rule 7 (`saturation` is a case token). The call is
    `solve_equilibrium` on a one-component problem; only Python's `Pure.saturation`
    carries the word.
22. The assembly's Hessian blocks scale by `N_p` in the `(N_p, x_p, v_p)` chart (§2.2); a
    driver in extensive coordinates sees `1/N_p`, and says which chart it uses.
23. Michelsen's deficit is `1 − Σ_i x_{p,i} exp(−r_{p,i}/RT)` with `r` the phase's
    chemical-potential row residuals at the converged multipliers (§3.4).
24. A mass-fraction feed needs molar masses, which `Parameters` carries; a molality feed
    needs the solvent set, which only the `Problem` can supply.
25. A boundary seed with no coexisting density roots is a named failure before the solve;
    an iterate that leaves the domain is Ipopt's `left_domain` stop, for every `k`.

## 11. The front surface: how a caller states the levers

The levers are the activation of §2: per phase an identity, a support and an amount
status; `T` and `P` each fixed or free; a feed in a stated basis; reactions; topology
declared or omitted. The question is how a Python call carries them to `Problem`.

**Decision: typed keyword arguments on one `Problem` object; nothing assumed; a mapping
form that is the same object serialised; convenience constructors that build a `Problem`
and show it.** Rejected: a free-form dictionary as the primary intake, a configuration
file at the call site, and defaults that pick a problem kind.

Why each:

- **Typed keywords, not a dictionary.** A dictionary is a vocabulary, and the loader that
  reads it branches on keys — the old tree's `_grepe_io.py` was 956 lines of exactly that,
  and its "request families" were the routes. A typed `Problem(phases=…, T=…, P=…,
  feed=…, reactions=…)` has one shape; Python does no validation beyond types and hands
  it to the native constructor, where compile rejects by name (§3.2). One validation, in
  one language, which is `TRANSFER.md`'s rule for parameters applied to problems.
- **Nothing assumed.** `T`, `P`, `phases` and `feed` have no defaults: a missing lever is
  rejected naming the lever, the same way an over-specified one is. The one default is the
  amount status, which is *free*, because the flash is the parent problem and every other
  cell is a pin on it. A support defaults to every species. Discovery is stated, not
  inferred from silence: `phases=Discover()`; an empty phase list is an error.
- **Free is a value, and it may carry a seed.** `P=Free()` or `P=Free(guess=2.0e5)`;
  `T=313.15` is fixed. The seed is a start, never a lever: it changes no row.
- **A feed names its basis by type.** `MoleFractions({...})`, `Amounts({...})`,
  `MassFractions({...})`, `Molalities({...}, solvent=[...])`. Python carries the numbers
  and the basis; the conversion is native (§7.4), and C4 is what makes the four
  indistinguishable to the rows.
- **The mapping form is the same object.** `problem_from_mapping` reads a JSON or YAML
  document whose keys are the constructor's keywords and nothing else, for campaigns and
  for the oracle rows. It adds no vocabulary; a key that is not a keyword is rejected.
- **Convenience constructors build a `Problem` and are not routes.** `bubble_point(
  liquid=..., T=...)` returns the `Problem` with the vapour pinned at zero and `P` free;
  `Pure.saturation(T)` returns the one-component form of the same. They live in Python,
  are three lines each, and every one of them is checked to give the same result as the
  explicit `Problem` it builds (one Python test per constructor, at the stage that lands
  it). A caller can always write the explicit form; the constructor only saves typing.
- **One native entry point.** `solve_equilibrium(model, problem)` — the name the brief
  fixed — compiles and runs the one Ipopt driver with the objective and constraint set
  compile chose (§2.4). Those are compile's data, not a caller's choice; a second public
  entry named for the boundary (`solve_phase_boundary`) would be the taxonomy on the
  surface, so it is not bound and `phase_boundary.cpp` holds seeds and continuation only. `find_phases` stays exposed for a caller who wants the discovery output alone.
- **The `Problem` shows its activation.** `problem.levers()` prints the phases with their
  statuses, which intensives are free, `k`, `f`, the reaction rank and the charge-row
  count, before any solve — the compile step's view, readable from Python, so that "what
  did I just ask for" has an answer that is not the solver's.

The MEA bubble point, stated in full:

```python
from epcsaft.equilibrium import Problem, Phase, Free, Pinned, Molalities, solve_equilibrium

problem = Problem(
    phases=[Phase("L"), Phase("V", support=["water", "CO2"], amount=Pinned(0.0))],
    T=313.15,
    P=Free(guess=2.0e5),
    feed=Molalities({"MEA": 5.0, "CO2": 2.0}, solvent=["water"]),
    reactions=bundle.reactions,
)
problem.levers()   # L free · V pinned 0 · T fixed · P free · k = 1, f = 1 · rank B = 7 · charge rows 0
result = solve_equilibrium(model, problem)
```

and pure saturation, which is the same call with one component:

```python
problem = Problem(phases=[Phase("L"), Phase("V", amount=Pinned(0.0))],
                  T=280.0, P=Free(), feed=Amounts({"ethane": 1.0}))
```

`Pure.saturation(280.0)` builds exactly that and nothing else. What this rules out is any
method on the EOS surface that solves coexistence without a `Problem` passing through
compile; the one-line test per constructor is what keeps it ruled out.

## 12. What regression needs from equilibrium: a shell around Ceres

Owner decisions, 2026-09-03, after the independent review and a reading of the three
downstream applications: equilibrium stays strictly generic; regression is a shell around
Ceres with no property-specific code; `State` is one batch object with an order flag
(`EOS.md`, "Properties"); ePC-SAFT parameters and pure-component thermochemistry are two
declared objects (`SURFACE.md` §2.5).

`EquilibriumResult` carries per phase `(T, v_p, x_p, N_p)`, the multipliers, rows, record,
and the sensitivity when asked. No property. A property of phase `p` is
`model.state(T, 1/v_p, x_p, order)`; a system total is `Σ_p N_p · property_p`, and this is
generic only if no extensivity table ever appears in `regression/`: the rule admits
`Σ_p N_p P_p`, which is the caller's error to make.

**Sensitivity levers: `T`, `P`, and the species feed** `n_feed,i`, one column per species.
`F_θ` for column `i` is `−B[:, i] / Σ_j n_feed,j + Σ_{q pinned} β_q (B x_q)_r / N_feed` on balance row `r` (the second term the same for every feed column, zero when every pinned amount is zero) and zero elsewhere
(balance row `r` is `(Σ_p N_p (B x_p)_r − (B n_feed)_r) / N_feed`, so its derivative in
`n_feed,i` is `−B[r,i]/N_feed − row_r/N_feed`, and the second term vanishes at a converged
point). The conserved amounts `b = B n_feed` are not a lever: `b` has rank-`B` entries and no
unique preimage. A caller who wants a derivative along any feed direction `d`, a conserved
coordinate included, forms `Σ_i d_i · dy/dn_feed,i` from the `C` columns, which is what the
absorber's loading-to-feed transform does; a direction in the reaction span (`B d = 0`) gives
exactly zero, because the solve is invariant to feed moved along a reaction. Trivial to assemble, and it is the
derivative the absorption film finite-differences today (`Reactive_Film.py:280-330`,
`d(log x)/d(log loading)` by `np.gradient`) and the derivative the superseded engine computes
natively but filters from its public rows (`reactive_finite_phase.cpp:2077`, thread A1). S9
carries a feed column. The reacting heat capacity `(∂H_eq/∂T)_{P,b}`, which the absorption
energy balance needs and which frozen-composition `cp` is not, is then the selector's chain
rule: `∂H/∂y · dy/dT + ∂H/∂T`, with no new code.

**The seventh active-parameter family.** The six EOS families (`helmholtz_tape.hpp:61-68`)
gain the `K(T)` correlation coefficients per reaction; `∂g°/∂θ` is already listed in `F_θ`
(`formulation.md:204`) and comes from the model's `standard_gibbs(T, order)` (`EOS.md`,
"Properties" §3). Without it neither application's fit is analytic: the amine fit shifts
`reaction:R*:correlation:a, b_k`; the lithium fit moves `ln_k_exchange`, `ln_k_na_uptake`.

The seventh family, `standard_gibbs(T, order)` and the totals land in the Properties stage after
Stage 6 (owner's ruling, 2026-09-03; `REGRESSION.md` §3, Staging), not in Stage 6.

An observation is data plus a selector plus a target:

    Observation { target: Problem | StateSpec; select: Expression; value; uncertainty }
    Expression =
        Leaf(Property(phase, field) | Coordinate(lever or amount n_{p,i}) | Constant)
      | Sum(w_j · e_j)
      | Product(e_j ^ p_j)
      | Log(e)

A weighted sum plus one ratio cannot express what the applications fit. Verified needs:
`pCO2 = y_CO2 · P` (a product) scored as `log10(pred/obs)` (a log); `net_na = na_feed −
na_aq` (a constant in a sum); extraction fraction `li_org / (li_org + li_aq)` (ratio of a
sum); selectivity (ratio of ratios). Four nodes cover all of them and finite-dose heat, and
the tree is smaller than a combinator list plus a residual-form flag.

`fit(model, observations)` builds one Ceres residual block per observation; the model output
is the expression evaluated on the solved or evaluated state; the Jacobian is the chain rule
per node: `Property` from `State` at `FirstDerivatives` composed with `dy/dθ`, `Coordinate`
from `dy/dθ`, `Sum`/`Product`/`Log` by their rules. No property name is known to
`regression/`; adding a field to `State` adds it to regression with no change.

Rules that follow:

- `StateSpec` declares its `DensityBranch`; `Unspecified` is refused inside a fit, as
  undeclared topology is.
- `Coordinate(T)` or `Coordinate(P)` on a target where that lever is fixed is rejected by
  name at observation compile, never a silent zero Jacobian.
- No `unit` field natively; the boundary is SI (`SURFACE.md` §3.4).
- Weights such as `1/(n₂ − n₁)` are caller data or computed natively; Python does none of
  that arithmetic.
- Rule 4 gains `Discover` in `STABILITY_SYMBOLS` so a fit cannot omit phases.
- The expression evaluator (value plus Jacobian) is one bound function, so a caller that
  wants a sensitivity outside a fit (the lithium script reads one by hand,
  `independent_network.py:855-858`) does not re-derive it.
- **No third derivatives, anywhere.** Ceres is Gauss–Newton or Levenberg–Marquardt on the
  residual Jacobian; it never forms the residual Hessian. Ipopt on a square feasibility
  problem needs none (§2.4). The implicit sensitivity needs second derivatives in `(y, θ)`,
  which the tape has.
- **Every regression observation declares its topology.** A parameter step that makes a
  declared two-phase point cease to exist is a *failed residual*, not a silently different
  problem: the cost functor returns `false`, Ceres shrinks the trust region, and nothing is
  thrown across the fit or reported as a zero residual.
- **Pure and mixture fits are one path.** A saturation observation is a `C = 1`, `k = 1`
  `Problem` through the same driver and the same sensitivity as a two-liquid tie line.

Budget: `TRANSFER.md` Stage 6 listed `eos_observable.cpp` (150) and
`equilibrium_observable.cpp` (250); with no property names left, both go, and the expression
tree is about 80 native lines net. `REGRESSION.md` is Stage 6's design document and carries
the Ceres kernel and the gate.

## 13. What the three applications need, and what the answer is

Read on 2026-09-03: MEA-Thermodynamics (the reaction-temperature fit running under the
superseded wheel at commit 8438ce5), Lithium_Extraction (the calibrated HBTA/TOPO reactive
two-liquid network at commit 90688d1), MEA-Absorption-Column (pinned at 8007a70, three
worktrees, and the two Codex threads: the absorber's three-plan diagnosis and the engine
thread that merged PRs #199–#201 and is designing A1 as issue #202). These citations are the
provenance of every row below.

| application needs | answered by | status |
|---|---|---|
| reactive VLE bubble pressure, vapour composition, speciation at fixed `T`, feed (MEA) | the `k = f = 1` cell with a vapour support restriction | §2–§3 |
| total enthalpy of a solved liquid; finite-dose heat (MEA calorimetry) | `State` totals + `Sum`/`Product` over two solves and `ideal_enthalpy` | `EOS.md` Properties §4; §12 |
| reacting `(∂H_eq/∂T)_{P,b}`, not frozen `cp` (absorber energy balance) | selector chain rule with `dy/dT` | §12 |
| equilibrium-path derivative `dy/db` at fixed `T, P` (absorber film tangent; removes its finite difference) | feed as a sensitivity lever | §12, S9, **new** |
| fixed-`(T, P)` chemical-potential and density derivatives in constrained log-composition coordinates (absorber film, `EpcsaftLiquidTransportState`) | `State` at `FirstDerivatives`, fixed-`P` block | `EOS.md` Properties §4, **new** |
| typed `K(T)` correlations with `∂ ln K/∂T`, van 't Hoff consistency, standard-state transfer with `∂/∂T` (all three) | `EOS.md` Properties §3 | **new** in the model |
| reaction constants and correlation coefficients as fitted parameters (MEA, lithium) | the seventh active family | §12, **new** |
| distribution ratios, extraction fractions, selectivity as fit outputs (lithium) | `Product`/`Sum`/`Log` nodes | §12 |
| two liquids with disjoint supports, cross-phase reaction, one ion-bearing phase (lithium) | supports, whole-system stoichiometry, C8 | §2–§4 |
| higher-order directional derivatives of the solved map for CasADi outer Hessians (absorber Plan A2) | **not provided**: third order is excluded by design; the absorber's own second candidate, explicit chemical-potential variables, needs only second derivatives of value maps, which the tape has | stated limit |
| hydraulics, mobility, loading transforms, column balances, collocation | not the engine's, by the absorber's own plan | — |

What the superseded engine added in the week of 2026-09-01 and the greenfield must not lose:
explicit reference thermochemistry with complete `h`, `u`, `cp` and equilibrium total
enthalpy (PR #199), reference-pressure consistency (#200), reaction slopes retained at
incipient boundaries (#201), and the A1 design (#202): certified state plus exact first
derivatives in `T`, `P` and the complete conserved feed coordinates, with the falsifiers it
lists (ideal reactive mixture with closed-form extent; same-branch centred ladders in `T`,
`P`, feed and mixed directions; `h = u + Pv` and `dh/dT|_{P,x} = cp` separately from the
reacting `dH_eq/dT`). Those falsifiers are the solve-tier checks for the totals and the
sensitivity and are adopted as written. The absorber still pins 8007a70 and supplies no
thermochemistry, which its own Plan 2 calls an adoption gap, not an engine gap.

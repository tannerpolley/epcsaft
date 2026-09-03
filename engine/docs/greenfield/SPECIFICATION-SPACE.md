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
| electroneutrality | `z · x_p = 0` | one per ion-containing phase **beyond the first** | `phase_electroneutrality` |
| normalisation | `Σ_i x_{p,i} − 1 = 0` | one per phase | `incipient_normalization` |

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
unchanged; `A_p`, `G` and the Hessian blocks scale by `N_p`. A pinned phase then
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
rows. With one ion-containing phase its neutrality follows; with `π_ch` such phases only
`π_ch − 1` rows are independent, and their multipliers `φ_p` are the electrostatic
potential differences between the phases, the first phase's potential being arbitrary.
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
still be singular by construction — all phases pinned, a redundant reaction, a redundant
charge row. **derived**

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
consumed in all three. S10 holds each consumer to finite differences of its own rows;
that there is one assembly is a placement decision (§5), not a check.

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
count above two** — the active subset is exactly seven registered rows and no others,
which is the evidence that a regression fit of pure-component parameters depends on
`equilibrium/` lightly (one compile, one assembly, one Newton) rather than structurally.
**derived**; S14 asserts the subset, and pin exchange at `C = 1` (pin `L` instead of `V`)
must give the same `P*`, which is the sharpest available test that "which phase is pinned"
is a specification entry and not a route. If the engine needs a distinct saturation solver,
it is visible here without solving: the compiled row multiset for ethane would differ from
the multiset the rule predicts.

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
pure-component coexistence path. The Clausius–Clapeyron seed `TRANSFER.md` Stage 7 names
lives in `phase_boundary.cpp`. The scientific reason outranks the tidiness one: a pure fit
and a mixture fit then share one convergence behaviour and one tolerance entry, so a
parameter set cannot converge cleanly at `C = 1` under a private solver and fail at
`C = 2` under the general one with no signal. S5 and S9 running through the same check
function at `C = 1` and `C = 2`, against the same policy entries, is that claim as a test.
`regression/` including `equilibrium/problem.hpp` is already what checker rule 3 permits
and rule 4 bounds; the owner's discomfort with the dependency is answered by its size, and
§4.2's Stage 7 file list is where the size is measured.

## 3. The space

### 3.1 Axes

| axis | values | what it changes |
|---|---|---|
| phase set, by identity | `{L}`, `{V}`, `{V,L}`, `{L,L}`, `{V,L,L}`, `{L,L,L}` | a phase's identity is (support, branch hint, amount status). Three-liquid systems exist with ePC-SAFT (Ascani et al. 2022, water + ionic liquid + alkane) and cost nothing at the compile tier. `V` versus `L` is a **start** attribute — the density-branch hint for the initial `η_p` — and changes no row. `{L,L}` is in the space precisely because it is where a vapour-like reference root would be assumed if any part of the start or the seed assumed one |
| amount status per phase | free; pinned `β = 0`; pinned `0 < β < 1` | `k`, and the balance's view of the phase; no row is added or removed |
| `T`, `P` | each fixed or free | `f` |
| reactions | none; one; two independent; two plus a dependent third; one that does not conserve charge | `rank B`, `g°`, the rejection |
| ions | absent; present in one phase's support with both signs; present in two liquids (`π_ch = 2`: one charge row, one interphase potential); a charged feed; a support with one ion sign; an incipient neutral phase beside a charged finite one | the charge rows (`π_ch − 1` of them), the rejection |
| feed basis | mole fractions; amounts; mass fractions; molalities | **nothing in the compiled system** — that is the check |
| topology | declared; omitted | declared is the main path; omitted is rejected by name until Stage 11, then accepted with `k = f = 0` |
| feed coordinate free with a target row | (later; the old tree's "affine feed boundary") | the same count rule with a target row on the right: a freed feed coordinate is paid for by one target equation. Recorded so the rule is seen to extend; not built before Stage 12 needs it |

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

The core axes multiply to about 3,800 specifications (6 phase sets × ≤ 8 amount-status
patterns — free or pinned-at-zero per phase, plus one `β > 0` pattern per two-phase set —
× 4 intensive patterns × 5 reaction values × 4 ion values); each compiles in microseconds
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
reject every square cell. **derived**

### 3.4 Declared topology: what it buys, and what it must refuse

A declared phase set is a specification, and it is what lets a solve skip stability
entirely: the rows are the equalities of §2.2 for every declared phase, including the
normalisation `Σ x_p = 1` as an *equality* for each. That is a claim about the answer —
every declared phase is present and distinct — and the engine either delivers exactly that
or delivers nothing. The rules:

- A solution is accepted only if every declared `N_p` is above the composition floor,
  every pair of phases is distinct (composition distance or packing fraction above the
  category-3 entries Stage 2 already carries from the paper), and the reduced Hessian is
  positive definite.
- Otherwise the call **fails with no equilibrium state**: no amounts, compositions,
  volumes or `G` are returned. What is returned is the failure name, the phase, and the
  quantity that condemned it. For a vanished phase that quantity is Michelsen's criterion
  (1994): a phase absent at equilibrium has `Σ_i x_i exp(…) < 1` where a present one has
  `= 1`, so the deficit `1 − Σ` is the report. For two phases that collapsed onto each
  other it is their distance. This is the owner's rule stated as a mechanism: the trivial
  solution is never an answer, and neither is a two-phase result with one phase empty.
- Declared topology does not make the solution unique. A declared set may admit several
  solutions (Coatléven and Michel 2022 §2.3: "the choice of present phases leading to the
  existence of a solution is not necessarily unique"; Othmer's up to three chemical-
  equilibrium solutions in a strongly non-ideal reactive ternary, cited in Ascani et al.
  2023 §2.1). The reduced Hessian certifies a local minimum; when the driver runs from more
  than one start the lowest `G` is taken; a global claim is never made (§8).
- The start for a declared multi-liquid set must split composition (§3.3). Without
  discovery the split direction is the caller's, or the default: one phase at the feed,
  the other displaced toward the pure component with the largest positive pair deviation.
  For a reactive set each start phase is first brought to chemical equilibrium along its
  own stoichiometric line, which is Ascani et al. 2023's step 3 and costs one single-phase
  reactive solve per phase.

Discovery (Stage 11) is the complement: phases omitted, Michelsen's inequality replaces
the equality, and the phase count is an output. The two are one row system with the
normalisation row as equality or as complementarity; nothing else differs.

### 3.5 The space the engine is built to solve, cell by cell

The owner's enumeration — neutral or electrolyte, non-reactive or reactive, any component
count, one phase of either kind, two phases, three if it is cheap — against what covers
each cell and at which stage. A stage number means solve-tier evidence exists; "compile
tier" means the rows, rank and rejections are checked and nothing is solved. `C` runs
from one to nine across the table; no cell has a component-count branch to exercise, and
C7 and S8 at nine components are the checks that see an index used where a component was
meant.

| system | `Π = 1` (`{L}` or `{V}`) | `Π = 2`, `{V,L}` | `Π = 2`, `{L,L}` | `Π = 3` (`{V,L,L}`, `{L,L,L}`) |
|---|---|---|---|---|
| neutral, non-reactive | `C = 1` ethane, S13 (7); `C = 2` binary, S13 (8) | `C = 1` saturation, S5 S9 S14 and the NIST oracle (7); `C = 2` flash, bubble, dew, fixed-fraction flash, S3–S5 and the Cai 1996 oracle (8, 9) | `C = 3` the reactive ternary with its reaction omitted, S6 S7 (8); cloud point (9) | compile tier on every fixture; opportunistic solve on the ternary |
| neutral, reactive | `C = 3` isomerisation, S11 (8) | a reacting liquid under a pinned vapour: compile tier (6); solve on the reactive ternary at Stage 9 | `C = 3` Ascani 2023 fixture, S16 (8); `C = 4` esterification tie lines, Senina 2022 and Schmitt 2006 (8, rows to transcribe) | compile tier |
| electrolyte, non-reactive | `C = 3` aqueous NaCl, S13 (8); the Hamer–Wu activity oracle is single-phase states, not an equilibrium cell (4) | `C = 3` aqueous NaCl, neutral vapour pinned, `P` free (9) | `C = 5` Khudaida 2026 tie lines under declared `{L,L}` (8): `π_ch = 2`, one charge row, one interphase potential | compile tier; the three-liquid ionic-liquid systems of Ascani et al. 2022 are the literature case if one is ever wanted |
| electrolyte, reactive | `C = 9` MEA speciation oracle (8) | `C = 9` MEA bubble point, neutral vapour pinned (9, 12) | **compile tier only**: nine species under declared `{L,L}` — the same rows with one charge row, but no fixture and no data | compile tier |

Two cells carry compile-tier evidence only, and they are named so that coverage is not
mistaken for proof: reactive electrolyte two-liquid, and every three-phase cell except
the opportunistic neutral one. The structure supports them because they add no row kind;
what they lack is a system known to exhibit them under this model.

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
| C4 | feed-basis invariance | the same feed as mole fractions, amounts, mass fractions and molalities compiles to identical rows and an identical feed *composition* `n_feed / Σ n_feed` to 1e-12 — the total implied by each basis differs (one mole, one kilogram, one kilogram of solvent) and S8 is what makes that irrelevant | molality conversion omits the solvent mass |
| C5 | reaction-basis invariance | reactions `{ν₁, ν₂}` and `{ν₁ + ν₂, ν₁ − ν₂}` (and any invertible mixing) give the same `B`, the same nullity, the same `g°` transform | build `B` from the reaction list in order instead of from its span |
| C6 | pin identity | the boundary specification's residual vector at a point equals the flash specification's residual vector at the same point with `N_q = 0`, on every shared row — the owner's "same rows, one entry differs", made numerical | the balance keeps the pinned phase at its start amount |
| C7 | permutation | relabelling components, or reordering phases, permutes the rows and leaves nullity and residual norms unchanged to 1e-14 | index a phase where a component was meant in any row |
| C8 | chart-or-row | an electrolyte flash compiles with nullity 0 — this catches the double count in `formulation.md` (§2 eliminates an ion by the electroneutral chart *and* §4 lists the charge row; both together are a redundant row) | add the eliminated coordinate back |

C6 is the check the brief asked for. It has a limit worth stating: if `phase_boundary.cpp`
assembles rows through the same function as the minimizer, no mutation of that function can
make C6 red, because both sides move together. C6's mutation therefore lives in the
specification handling, and the checks that catch a *bypass* of compile are S3–S5 below.

### 4.2 Solve tier — the minimal input set

Six parameter sets. Four already exist in `engine/tests/native/model_support.hpp`
(**measured**: `ethane()` :42, `neutral_binary()` :51, `associating_binary()` :60,
`aqueous_electrolyte()` :87, `nine_species()` :116). Two come from the literature, not
from the old tree:

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
2015. Recommended as an oracle case for Stage 8 (declared `{L,L}`), which puts reactive
LLE under measured data before any discovery stage exists. The rows are to be transcribed
from the papers already attached in Zotero; nothing needs acquiring.

| id | fixture | cells | stage |
|---|---|---|---|
| ethane | `ethane()` | `{L}` alone and `{V}` alone at fixed `T`, `P`; `{V,L}` with `V` pinned, `P` free; the same with `L` pinned (must give the same `P*`); then `T` free. The Stage 7 file list is the record that these ran with no `gibbs_minimizer.cpp`, no `phase_stability.cpp`, and `B = I` | 7 |
| binary | `neutral_binary()` | `{V,L}` flash inside the envelope; bubble; dew; `β`-flash; both intensives | 8, 9 |
| reactive ternary | Ascani 2023 `A + B ⇌ C` | **required:** `{L,L}` with the reaction, from three feeds in the gap (S16); `{L}` with the reaction (S11 uses an isomerisation variant of the same components); `{L,L}` without the reaction (S6, S7); `{L,L}` with `L₂` pinned, `T` free, with and without the reaction. **Opportunistic:** `{V,L,L}` and the three-phase point, if the fixture admits them at some state; recorded either way | 8, 9 |
| electrolyte | `aqueous_electrolyte()` | `{L}`; `{V,L}` with `V` support `{water}` and `V` pinned, `P` free — the MEA shape at three species | 8, 9 |
| electrolyte two liquids | oracle `aqueous_alcohol_salt_lle` (Khudaida 2026, packet revision 2) | declared `{L,L}` from the tie-line midpoint feed; `π_ch = 2`, one charge row, one interphase potential; both endpoint compositions within the row uncertainties. Measured data on the electrolyte two-liquid cell before any discovery stage exists; Stage 11 later reaches the same rows with topology omitted | 8 |
| reactive two liquids, real | Senina 2022, Schmitt 2006; parameters and `K_a` from Ascani 2023 Tables 1, 2, 4 | declared `{L,L}` on the nine quaternary tie lines; a published-model reproduction. A new oracle case once the rows are transcribed | 8 |
| nine species | `nine_species()` | compile tier only: the composite cell's rows, rank, nullity, and reject-by-name with the bundle's two reactions; and declared `{L,L}` compiled with one charge row, which is the reactive electrolyte two-liquid cell's only evidence. The solve is Stage 12's gate | 6 |

| # | claim | asserts | mutation |
|---|---|---|---|
| S1 | rows satisfied | every solved cell's residual vector is within `max(tolerance, floor)` of zero on every row, the two numbers read from the result's own `Convergence` record (`TOLERANCES.md` §2, consequence 3) — including rows the driver did not iterate on (affinity for a reactive cell; Euler's `g_p = λ·Bx_p` for every phase) | drop one balance row from the Ipopt constraint set: affinity goes nonzero |
| S2 | local minimum | for `k = 0`, the reduced Hessian `Zᵀ ∇²L Z` is positive definite (`EqID: reduced_hessian_stability`) | flip the sign of the `P V` term in the objective |
| S3 | dew ↔ bubble | for the binary at `T`: bubble from liquid `x` gives `(P*, y*)`; dew from `y*` must give `(P*, x)` to the convergence tolerance. An exact symmetry of the formulation, and it fails on any asymmetry between "which phase is pinned" | start the pinned phase on the vapour branch regardless of its declared identity |
| S4 | `β → 0` | the fixed-fraction flash at `β = 10⁻², 10⁻⁴, 10⁻⁶` has `|T(β) − T_bubble|` decreasing and bounded by a constant times `β` (each step of a hundred in `β` shrinks the gap by at least ten) | the balance omits the pinned phase for `β > 0` |
| S5 | round trip | solve `P` at fixed `T`; solve `T` at that `P`; recover `T` | freeze `T` in the Jacobian when `T` is the unknown |
| S6 | instability certificate | `G₂ < G₁` at the ternary feed; `G₃ ≤ G₂ + tol` for the three-phase cell | return the one-phase start as the result |
| S7 | lever rule | two feeds on the `{L,L}` tie line give the same two phases to 1e-8 | seed each phase from its own feed instead of from the tie line |
| S8 | extensivity | every solved cell with the feed doubled: amounts double, every intensive unchanged | an absolute trace floor in the start |
| S9 | sensitivity | `dy/dθ` from `implicit_sensitivity.cpp` against a centred difference of the solve in `θ` to 1e-6 (the Stage 3 precedent for the FD floor) | flip the sign of `F_θ` |
| S10 | Jacobian agreement | the assembly's `F_y` against a centred difference of `F` in every unknown, at a boundary cell and at a flash cell, to 1e-6; and Ipopt's `eval_h` at the flash solution against a centred difference of the Lagrangian gradient (`eval_grad_f + Σ λ_j ∇c_j`) to 1e-6. The earlier form of this check, elementwise equality between the two drivers' matrices, was false: the drivers' rows differ by the factor `N_p` per phase (§2.2), so an implementer could only pass it by fudging. "One assembly" is a design statement, held by §5's file placement and rule 7's file list, not by a check | scale one Hessian block by `0.5` |
| S11 | ideal-gas limit | the reactive one-phase cell with the isomerisation `α ⇌ β` (so `Σ ν_i = 0` and the `ρ_ref` term cancels) at `P → 0` approaches `x_β / x_α = exp(−Δg°/RT)` (a limit; the residual terms vanish, which Stage 2 already checks). A reaction with `Σ ν ≠ 0` has a density factor in its limit and is not used here | perturb the `g°` transform |
| S12 | discovery agrees (Stage 11) | topology omitted at the ternary feed finds the phase count of the declared `{L,L}` flash and the same `G` | accept the first sample as a witness |
| S13 | single phase | `Π = 1` on every fixture at fixed `T`, `P`: the solved composition equals the feed to 1e-14, and the solved volume equals `model/`'s `solve_density` root at the same `(T, P, x)` within the larger of the two solves' `max(tolerance, floor)` — the pressure row *is* the density equation when there is one phase | seed the start from phase index 1 |
| S14 | `C = 1` is the general row set | the ethane boundary's compiled rows, counted by `EqID`, are exactly two `chemical_potential_stationarity`, two `pressure_equality`, one `material_balance`, two `incipient_normalization`, and nothing else; pinning `L` instead of `V` gives the same `P*` to the convergence tolerance | a `component_count == 1` branch in compile that omits the balance row |
| S15 | a declared phase that vanishes | the binary `{V,L}` flash at a state outside its envelope, and the `{L,L}` flash at a feed outside the gap, **fail with no equilibrium state** (§3.4): the failure names the phase and carries Michelsen's deficit `1 − Σ` or the phase distance; no amounts, compositions or `G` are returned, and a symmetric start on any two-phase cell is refused the same way | return the converged iterate when one `N_p` sits at its bound |
| S16 | reactive tie line is feed-invariant | on the `A + B ⇌ C` fixture at fixed `T`, `P`, the declared `{L,L}` reactive flash from three different feeds inside the gap returns the same two phase compositions to the convergence tolerance, and affinity vanishes in both phases with one `λ`. Gibbs' rule for a reactive system, `F = 2 − π + C − R = 0` here, says the tie line is unique (Ascani et al. 2023 §2.1, "unique reactive tie-line"); a route that solved the reaction in one phase and partitioned into the other would move the tie line with the feed | drop the affinity consequence by solving the balance with `B = I` |

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
path itself, and it is a derivative-agreement check of the kind `TESTING.md` already
admits; that the two drivers share one assembly is enforced by where the assembly lives
(§5) and by rule 7 covering that file, not by a test. The artifact's "seven genericity equivalences" (round trip, reaction basis,
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
2. **compile** (`request_compile.cpp`, Stage 6) resolves the feed basis to amounts
   (native — see §7, item 4; a molality feed needs the solvent set, which is a `Problem`
   field, not an inference), builds `B` from the reaction span by SVD, applies the count
   rule, assembles the active row list with its `EqID`s, takes the rank at the generic
   point of §3.3, and either returns a `CompiledSystem` or throws naming the deficiency.
   It is the one file allowed to branch, and only on status flags.
3. **assembly** — `CompiledSystem::evaluate(y)` giving `F(y)` and `F_y`, **also Stage 6,
   in `request_compile.cpp`**, not in `gibbs_minimizer.cpp`: Stage 7's saturation and
   `implicit_sensitivity.cpp` need it before the Ipopt driver exists, and a Stage 7 that
   cannot find it will write its own, which is the failure this document exists to
   prevent. One tape evaluation per phase per iteration at unit amount (§2.2) — value,
   gradient, Hessian of `A_p` in `(T, v_p, x_p, θ)` — scaled by `N_p`, plus the linear
   constraint blocks. The Hessian of `A_p` is the tape's reduced Hessian; the site-fraction
   elimination happened in Stage 3 and nothing here sees sites. Because the assembly moves
   here, `equilibrium/request_compile.cpp` joins rule 7's `SOLVER_FILES`
   (`tools/check-engine-architecture.py:163`): it may branch on status flags, and it may
   not name a problem kind — the comment at `:154-156` is updated to say so.
4. **drivers**: Ipopt for `k = 0` with bounds on `ln η_p`, `ln x_{p,i}` and `N_p ≥ 0`,
   exact Hessian, `nlp_scaling_method = none` (`TOLERANCES.md` §3.4), objective `G` and
   constraint rows read from the `CompiledSystem`; damped Newton on `F` for `k ≥ 1`, seeded
   from the flash solution or the previous boundary point when the caller is sweeping
   (`TRANSFER.md` Stage 7), and never an outer loop around `solve_equilibrium`
   (`TOLERANCES.md` §6). Every multiphase start must give the phases *different*
   compositions — the trivial solution is a stationary point (§3.3) and a symmetric start
   converges to it.
5. **result**: amounts, volumes, intensives, multipliers, every residual row, the reduced
   Hessian's smallest eigenvalue, counters, a named status. Sensitivity when asked. Or,
   under §3.4, no state at all: a failure name, the phase, and the number that condemned
   it.

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
7. **Saturation is equilibrium, fully.** The owner's tension is real and the resolution
   holds under test: §2.5 counts the `C = 1` cell against the general rows with no special
   case and gets saturation, so a private pure-component solver would be a second
   implementation of seven rows that already exist. **measured** residue check: nothing
   saturation-shaped exists under `native/` today — no seed, no estimate, no coexistence
   path, no policy entry — and `infinite_dilution_reference.cpp` is a sequence of
   single-phase states, correctly placed. Two planned items were residue and
   `target-layout.md` is revised: the `Saturation` `py::class_` leaves `eos_bindings.cpp`
   (the binding is `solve_phase_boundary`'s, in `equilibrium_bindings.cpp`), and `eos.py`'s
   `Pure.saturation` is a thin call onto that binding, kept where users reach for it,
   because the Python surface is organised by capability and was never meant to mirror
   the layers.
8. **The test budget.** `engine/tests/native/` holds 1,187 lines after two closed stages
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
- **The discovery cell asserts a rejection until Stage 11.** Between Stage 6 and Stage 11,
  "topology omitted" is proven only to fail by name.
- **A single-phase cell proves consistency, not stability.** S13 shows the rows and the
  density closure agree; whether that one phase is the stable state at `(T, P)` is Stage
  11's question, and until then a `{L}` result carries no stability claim.

## 9. Cost and placement

All native, under the existing `CHECK` harness (`engine/tests/native/physics_checks.cpp:18`);
no framework, no fixture library, no Python until Stage 4's bindings exist, which under the
revised order (`TRANSFER.md` §"Execution order") is before Stage 6 anyway. **guessed**:

| file | lines | ceiling | stage | holds |
|---|---|---|---|---|
| `tests/native/equilibrium_support.hpp` | 140 | 190 | 6 | includes `model_support.hpp`; adds the Ascani 2023 reactive ternary, its reaction and `g°`, specification builders, the four-enum lattice walk, the generic point, the count-rule predictor |
| `tests/native/specification_checks.cpp` | 180 | 240 | 6 | C1–C8 |
| `tests/native/equilibrium_checks.cpp` | 170 | 220 | 8, then +15 at 11 | S1, S2, S6, S7, S8, S11, S13, S15, S16 on the flash cells; S12 at Stage 11 |
| `tests/native/boundary_checks.cpp` | 140 | 180 | 7 opens it (ethane: S5, S9, S13, S14 at `C = 1`, ~60); 9 completes it | S3, S4, S5, S9, S10, the `{L,L}` and `{V,L,L}` boundary cells |
| `tools/check-engine-architecture.py` | +2, then +12 | — | 6, then 8 | `request_compile.cpp` joins rule 7's file list at Stage 6; rule 8 (`eval_h` present, no `limited-memory`) at Stage 8 |
| **total** | **~640** | **~840** | | |

Stage 7 is where the "light dependency" claim of §2.5 is measured rather than argued: its
brief's file list is `phase_boundary.cpp`, `implicit_sensitivity.cpp` and
`boundary_checks.cpp`, and its gate runs with `gibbs_minimizer.cpp` and
`phase_stability.cpp` absent from the tree. `STAGE-LOG.md` records the list.

Order inside Stage 6: `equilibrium_support.hpp` and C1–C3 first, before any solver exists —
compile and the count rule are the first thing the stage writes, and C1–C3 are green on
the compile step alone. That is also the cheapest possible preflight for the stage brief:
if C1's predictor and the compile step disagree on a cell, one of them misreads §2.3, and
that is found before a line of `gibbs_minimizer.cpp` is written.

## 10. Traps, one line each — the brief copies this list

Every item below is a place an implementer would otherwise guess, or a place where a check
as first written could only be passed by fudging. Found on a re-read of this design against
`TOLERANCES.md`, the model layer and the checker; the first four are corrections to this
document's own first draft.

1. Evaluate every phase's tape at unit amount and scale by `N_p` afterwards (§2.2), or a
   pinned phase is `0/0` and a second path grows.
2. The rank point must give every phase a different composition (§3.3), or every square
   multiphase cell is rejected at the trivial solution.
3. S10 compares against centred differences; the two drivers' rows differ by `N_p` per
   phase and are not elementwise equal.
4. C4 compares feed *composition*; the totals implied by the four bases differ.
5. The assembly `F`, `F_y` lives in Stage 6's `request_compile.cpp`, and that file joins
   rule 7's list; Stage 7 must find it there, not write it.
6. Multiphase starts must split composition; the trivial solution is a stationary point and
   S15 requires it to be reported by name, never returned as `converged`.
7. S1 and S13 test against `max(tolerance, floor)` from the result's `Convergence` record,
   never against a bare policy number — `TOLERANCES.md` §2 says why a bare number can be
   unreachable.
8. S11 uses the isomerisation only: with `Σ ν ≠ 0` the ideal-gas limit carries a density
   factor and the closed form the check names is wrong.
9. A molality feed needs a declared solvent set on `Problem`; compile does not infer one.
10. A support restriction is a zero amount passed to the tape; `validate_point`
    (`helmholtz_model.cpp:250`) accepts zeros, and a pure-water vapour inside an
    electrolyte model takes the no-ion tape branch — Stage 3's predicate-signature cache
    key exists for exactly this, and the assembly must not bypass it.
11. `nlp_scaling_method = none`, bounds on `N_p ≥ 0`, and no outer Newton around
    `solve_equilibrium` — three commitments from `TOLERANCES.md` that S1, S15 and the
    Stage 9 gate inherit.
12. The `C = 1` seed needs a pressure at which both density roots exist; the density
    solver's `stable_root_count` says whether it does, and a seed where it does not is a
    named stop, not a retry loop.
13. S3's dew is seeded cold, not from the bubble solution, or the symmetry is trivial; the
    state is chosen subcritical for both components so the dew pressure is unique.
14. Nothing asserts a variable count, an iteration count, or which driver ran. Counters
    are printed per cell; the design intent that the exact Hessian is consumed is visible
    in `hessian_calls`, and rule 8 is the only thing that turns it into a failure.
15. One electroneutrality row per ion-containing phase **beyond the first** (§2.2); one
    per charged phase is a dependent row and the rank check will reject valid electrolyte
    flashes.
16. A declared phase set is delivered whole or not at all (§3.4): a vanished phase, a
    collapsed pair or the trivial solution fails with no state, carrying the phase name and
    Michelsen's deficit or the distance.
17. The reactive two-liquid fixture is Ascani et al. 2023 Table 3 with `K_a = 2.25`, not
    the old tree's immiscible ternary; its existence is the paper's Figure 3.
18. At fixed `T`, `P` a reactive tie line in a one-reaction ternary has `F = 0` and is
    feed-invariant (S16). Any feed dependence is a route solving the reaction in one phase.
19. Start phases of a declared reactive multi-liquid set are first brought to chemical
    equilibrium along their own stoichiometric lines (Ascani et al. 2023 step 3), one
    single-phase reactive solve each; the feed itself stays on the stoichiometric line and
    is *not* on the chemical-equilibrium surface after the split.

# Greenfield transfer: build order and per-stage gates

Branch `greenfield`, worktree `~/Workspaces/Engineering/ePC-SAFT-greenfield`, cut from `main @9eccbfd`.
Companion: `CARRY-OVER.md` (what is copied, what is written fresh, what has no successor).

## Why the issues change shape

Issues #185–#196 were written against the existing tree, so most of their content is
deletion: "delete HELD v1 and steps 8–10", "delete coexistence and declared routes",
"replace `(S+1)` Hessians", "strip Python numerics". In a greenfield those clauses are
not work. They are the absence of work.

What survives from each issue is only its positive half — the thing to write — plus its
acceptance threshold, which carries over unchanged because it was always stated against
data and counters, never against the old code.

| issue | migration content | survives as | writing work |
|---|---|---|---|
| #195 layout and gate | move `_native/` → `native/`, reconcile three archives to one, fix the checker's contradicted module ban | the tree shape and the four checker rules | Stage 0, ~180 lines (checker + meson) |
| #186 numerical policy | find 15 `R` definitions, 12 duplicated names, unify | `physical_constants.hpp` + `numerical_policy.hpp` written before any consumer | Stage 1, 280 lines |
| #189 EOS derivatives | replace `(S+1)` Hessians, add tape cache, Newton association | one Lagrangian Hessian, site elimination, cached `ADFun` | Stage 2, 350 lines |
| #190 EOS Python | strip numerics from Python, cut exports to 11, replace `py::dict` round trips | Python as a pure mapping layer, `py::class_` bindings | Stage 3, 1,100 lines |
| #191 Problem/Result | collapse 6 specification types and 2 request families | one degrees-of-freedom `Problem`, one `EquilibriumResult` | Stage 4, 1,050 lines |
| #187 GREPE ownership | delete `coexistence.cpp`, the declared routes, `solve_incipient_pressure` | one `GibbsMinimizer`, `solve_phase_boundary` as its rows | Stage 5, 1,650 lines |
| #188 phase stability | delete HELD v1 (3,051) and HELD2 steps 8–10 (3,338); fix the neutral-chart sentinel | `find_phases` written to Pereira 2012 | Stage 6, 1,100 lines |
| #192 regression | remove owner compilation from the fit path | `fit(model, observations)` on Ceres | Stage 7, 890 lines |
| #196 consumers | migrate analyses and two downstream repos, rename result fields | the same migration — this one does not change | Stage 8, unchanged scope |
| #181 policy-edit gate | reject unreviewed solver-policy edits | folds into Stage 1's duplicate-constant gate | — |
| #178 build retention | bound build-tree and wheel-cache growth | orthogonal to the tree; stays on `main` | — |

Two issues are wholly absorbed (#181 into #186's gate, #178 unaffected). Six lose their
deletion halves entirely. One (#196) is unchanged, because consumers are outside the tree.

## Scale

| | today | target | cut |
|---|---|---|---|
| native C++ | 43,910 (80 files, `_native/` + `cpp/`) | 10,020 (29 files) | 77% |
| Python | 15,743 | 2,430 | 85% |
| tests | 4,404 | 1,650 | 63% |
| **total** | **64,057** | **14,220** | **78%** |

The 78% is the audit's number reached by construction rather than by deletion. It is a
prediction, not a budget: a stage that needs more lines to be correct takes them, and the
overrun is recorded here rather than absorbed silently.

## Execution order: MEA first, stability last

Revised 2026-09-02 on the owner's priority: the MEA reactive workload comes first, phase
stability goes last. Stage *numbers* are kept as names — they are in commit messages,
`STAGE-LOG.md`, and the checker — and only the execution order changes. Two stages split.

Nothing on MEA's path calls `find_phases`: its case is a finite electrolyte liquid against
an incipient neutral vapour, topology declared. Stage 6 is not merely deprioritised, it is
off the critical path.

| run | stage | why MEA needs it |
|---|---|---|
| 1 | **2** `model/` | the Lagrangian Hessian and site elimination, 7–10× per evaluation — MEA's single biggest lever. Also `infinite_dilution_reference.cpp`, which its molality standard state requires |
| 2 | **3** bindings + Python EOS | every name on MEA's surface |
| 3 | **4** minimizer | fixed-T,P finite reactive electrolyte liquid: the speciation states |
| 4 | **5a** `phase_boundary.cpp`, P unknown | the bubble point itself. **MEA switches here** |
| 5 | **8a** MEA cutover | the 39→11 and 53→9 contraction |
| 6 | **7** `regression/` | MEA fits only the neutral CO₂–water binary, an EOS-observable path |
| 7 | **5b** saturation, `implicit_sensitivity.cpp` | owed to Stage 7's saturation gate and to Lithium |
| 8 | **6** `phase_stability.cpp` | nothing |
| 9 | **8b** composite gate, analyses, Lithium, merge | — |

**Stage 4 versus Stage 5, disambiguated.** §"Problem as a degrees-of-freedom specification"
lists "fixed T with an incipient partner solving P" among the patterns that fall out without
a branch, and Stage 4's gate *accepts* that pattern in the well-posedness check. But it
cannot be a row of the minimizer: an incipient phase has zero amount, and `ln n_p` is `ln 0`.
Its rows are boundary rows. So the MEA path is 2 → 3 → 4 → **5a**, and the reactive bubble
point moves from Stage 8 to 5a's gate, where the code that must satisfy it is written.

`implicit_sensitivity.cpp` (250) leaves Stage 5 for 5b; `phase_boundary.cpp` (300) stays.
Stage 5b keeps the saturation-fit gate explicitly, so Stage 7 running earlier does not
silently lose it. Stage 8b re-runs 5a's gate as its closing step, so something still
exercises every layer at once after 6 and 7 land.

Reactive LLE remains absent from the oracle. MEA does not need it; it belongs with Stage 6.

## Correction: the Stage 8 reactive gate is a golden file

`TESTING.md` bans golden files of computed output — "a number that came out of this code is
not evidence that this code is right." The Stage 8 gate as first written pins the liquid
molar densities from `reactive_vle_baseline.json` (46,425–47,148 mol/m³), which are exactly
that: this engine's own output. The rule catches the plan that introduced it.

Replacement rows, all real measured data already transcribed with source hashes in
`MEA-Thermodynamics/data/reference/MEA/`:

| oracle case | source | gates |
|---|---|---|
| `aqueous_amine_binary_vle` | Cai 1996, isobaric T-x-y at 101.33 kPa | 5a — boundary rows on a neutral associating binary, no ions, no reactions; carries the Baygi reproduction |
| `aqueous_amine_carbon_dioxide_speciation` | Jakobsen 2005 / Böttinger 2008 NMR, 293–353 K | 4 |
| `aqueous_amine_carbon_dioxide_reactive_vle` | Jou 1995, Hilliard 2008, Xu 2011, 313–393 K | 5a |
| `loaded_amine_density` | Amundsen 2009, ±0.0005 g/cm³ | 3, 4 — the measured counterpart of the computed density the baseline pinned |

**These rows test the parameter set as much as the code**, and the current bundle is 2.26×
off in rmse-factor on pCO₂, so they cannot be pass/fail on agreement. Within `TESTING.md`
they can gate three things: that every row *solves*; that the solved point satisfies its own
equations (equal chemical potentials, conservation, charge, Σy = 1 to 1e-10); and a
**published-model reproduction** — a fact about the world in the same sense as the
Khudaida tie-lines.

*Corrected 2026-09-02.* This first named Uyan, Sadowski & Held 2015. That paper is
**MDEA**–CO₂–H₂O, not MEA, and the campaign that attempted it
(`codex/bulow-2021-fit-reproduction:validation/campaigns/2015-uyan/analysis.yaml`) records
its own failure: `grepe_ard_percent: 118.98`, `promotion_status: not-promotable`. It fails
as a gate on both the system and the evidence.

The replacement is **Baygi & Pahlavanzadeh 2015** (PC-SAFT, MEA–H₂O; MEA 3B, water 4C,
`k_ij = −0.0520`; Table 3 reports AAD 3.54% in x_MEA and 1.54% in y_MEA against Cai 1996).
Its parameters are *already* the Stage 2 associating-binary fixture — `2383.4744`,
`0.118488`, `k_ij = −0.052` in `strongly_associating_binary.json` **is** Baygi's model — so
the reproduction costs no new data. One item to read before asserting: whether Baygi's
Eq. 11 cross-association rule matches the fixture's Wolbach–Sandler volume
`0.1540467874409226`.

## Build order

The order is forced by the include direction (`physics ← model ← equilibrium ← regression
← bindings`) except that bindings and the Python mapping layer land early, at Stage 3, so
that every later stage is testable from Python against the oracle. Nothing is written
without a way to run it.

Each stage is one PR onto `greenfield`. A stage closes when its gate passes; no stage
begins before the previous one's gate passes.

---

### Stage 0 — Tree, build, checker

Write: `engine/native/{physics,model,equilibrium,regression,bindings}/` (empty),
`meson.build` producing one static archive `epcsaft_native` and one extension `_core`,
and `tools/check-engine-architecture.py` with five rules:

1. exactly one static archive and one extension in the Meson products;
2. `nm` on the archive shows no `PyInit`/pybind11 symbols;
3. every `#include "..."` under `native/<layer>/` resolves to the same or a lower layer;
4. `regression/` does not reach phase stability;
5. no `constexpr` name defined in two translation units, and no category-1 constant
   outside `physics/physical_constants.hpp`.

Rule 4 replaces the module-level regression ban that the old layout contradicted, with a
symbol rule: `regression/` may include `equilibrium/problem.hpp` and call
`solve_equilibrium` / `solve_phase_boundary`; it may not include `phase_stability.hpp`,
call `find_phases`, or construct a problem with phases omitted. A fit needs
differentiable residuals, and a discovered phase count is not one.

Delete in the same commit: `engine/src/epcsaft/_native/`, `engine/cpp/`,
`engine/tests/*/numerical/`, `engine/tests/*/behavior/`, and every Python module except
the three that survive (`CARRY-OVER.md` §4). The branch does not build after this commit
and is not expected to until Stage 3.

**Gate:** the checker reports all five rules green over the emptied tree, and each rule
is shown to fail on a deliberate violation. `meson setup` does not configure until Stage 1
supplies the first sources — a Meson archive with no source files is an error — so the
build gate begins at Stage 1.

---

### Stage 1 — `physics/` (950 lines)

`physical_constants.hpp` (30), `numerical_policy.hpp` (250), `helmholtz.hpp` (200),
`association.hpp` (150), `electrolyte.hpp` (260), `born.hpp` (60).

Mostly carry-over (`CARRY-OVER.md` §1): the Helmholtz terms, the association Newton
solver, Debye–Hückel, Born, and four of the five permittivity forms are correct and
validated, and are copied with their `// EqID:` annotations intact. The written work is
`numerical_policy.hpp`: every tolerance in the tree as a `Tolerance` entry with its
category (1 physical constant / 2 thermodynamic constraint / 3 algorithm tolerance /
4 unjustified legacy). Category-4 entries are not carried; each one is either justified
into category 3 with a source or replaced with the paper's value.

The HELD2 tolerance drift is resolved here, not in Stage 6: the candidate gap `ε_b`
enters as the paper's `1e-2` (not the drifted `5e-2`), and the composition floor as
`1e-8` (not `1e-10·f_k`). The seven density-root tolerances that HELD2 owned
(`root_*`, `mechanical_margin`, `stable_objective_tie`) enter as `model/density` entries.

**Gate:** closed-form, identity, limit, and permutation tests on the terms
(`tests/eos/`, ~200 lines). Checker rule 4 passes. No file in `physics/` includes CppAD
or a solver.

---

### Stage 2 — `model/` (1,810 lines)

`helmholtz_tape.hpp` (350), `helmholtz_model.{hpp,cpp}` (750), `density.cpp` (250),
`infinite_dilution_reference.cpp` (450).

This is #189, written rather than migrated. Three properties are designed in:

- **One Hessian, not `(S+1)`.** The implicit association second derivative comes from
  one Hessian of the Lagrangian `a + λ·H`, with `λ` from one linear solve of `H_X`.
  Measured multiplier at `S = 6`: 7×.
- **Sites are not tape independents.** They are eliminated by the Lagrangian, so the tape
  width at nine species is 11, not 17: `1 + 2w` drops 35 → 23 sweeps. Another 1.5×.
- **One `ADFun` per topology, cached.** Worth 1–10%, and last in importance — the audit's
  original ordering had this first, which the measurements falsified.

`density.cpp` is a bracketed scan plus Newton on the pressure row, seeded from an anchor
when the caller has one. The 360-interval scan and the 64×depth-8 pressure envelope are
not written; a value-plus-`dP/dρ` path (4 sweeps) serves the root iteration, not a full
Hessian set.

**Gate, correctness:** the `(S+1)`-Hessian tensor and the one-Lagrangian-Hessian tensor
agree elementwise to **1e-10** on the nine-species case — two exact AD formulations, so the
threshold is real. This licenses the 7× and is the only thing in Stage 2 that can
invalidate the floor. Derivatives against centred finite differences agree to **1e-6**,
not 1e-10: centred FD in double precision bottoms out near 1e-8, so the tighter figure this
gate first carried was unreachable by construction.

Three conditions the floor depends on, which the original gate stated as bare numbers:

- **`≤ 23 + 2P` sweeps holds only if `CppAD::ADFun::Hessian` is not used.** It discards the
  `Forward(1)` outputs, forcing `S` extra reverse sweeps for `Φ_Xy` (31 at nine species).
  The count is met only by driving `Forward(1)`/`Reverse(2)` directly and reading `Φ_Xy`
  off the forward sweeps.
- **"One `ADFun` per topology" is an insufficient cache key.** Two predicates in
  `physics/electrolyte.hpp` are evaluated from `double` inputs *while recording* and baked
  into the tape: the `finite_size_chi` series switch per ion (`:332`) and the
  `reference_charge_square_sum == 0.0` branch (`:242`). Both flip inside Stage 2's own
  scope — the dilute end of the molality ladder takes the series, a pure-solvent state
  inside an electrolyte model takes the no-ion branch. A tape keyed on topology alone is
  silently wrong across them, with no error raised. Key on (active set, predicate
  signature) and re-record when it changes.
- **Association iteration counts need a three-line physics edit** — an optional
  `std::size_t*` out-parameter on `solve_site_fractions`. Stage 1 removed the kernel
  counters on the grounds that the tape owns evaluation accounting; without this the row is
  reported as unmeasured.

**Gate, performance (reported, not pass/fail — `TESTING.md`):** the nine-species MEA model
at 313.15 K, ρ = 3.0e4 with one active parameter, value + gradient + Hessian + parameter
Jacobian in **≤ 1.0 ms**, at **≤ 23 + 2P sweeps** and **≤ 1 Hessian call** per evaluation;
the associating binary at 376.88 K in **≤ 5 association iterations** per solve (today 126).

**Measurement is a native harness, `tests/native/model_floor.cpp`, not
`oracle/measure_floor.py`.** The Python script opens `import epcsaft` and reads its
counters from `_core._instrumentation()` (`measure_floor.py:12,21`), and `_core` is built
from `bindings/` at Stage 3 — so the gate as first written could not run at the stage that
makes the claim. The harness times and counts the tape directly, on the same case, median
of 20 after 3 warm-ups, `OMP_NUM_THREADS=1`, recording commit, compiler and CPU. Stage 3
re-runs the same case through `measure_floor.py` and the two must agree; a divergence is
marshalling, which is Stage 3's to account for.

*This defect was found by executing the plan's own claims against the tree before Stage 2
began, which is the mechanism described in `ORCHESTRATION.md`.*

---

### Stage 3 — `bindings/` and Python EOS (1,100 + 1,240 lines)

`bindings/{module,eos_bindings}.cpp` (550), and `src/epcsaft/` `__init__.py`, `units.py`,
`parameters.py`, `eos.py`, `reference.py`, `_core.pyi`.

This is #190. Python maps inputs to `_core` and `_core` results to types, and does no
arithmetic — the rule the user set: *any Python frontend method call is completely
agnostic to what it is solving; it takes the inputs it is given, gives them to C++, gets
the answer back, and returns the output.* Bindings are `py::class_`, not `py::dict` round
trips. Eleven exports. Parameters are validated once, in the native constructor; the
three Python validation passes are not written.

**Gate:** oracle case `aqueous_alkali_halide_activity` (24 rows, LiCl/NaCl/KCl mean ionic
activity coefficients) reproduced through `Mixture.reference` within its recorded
convergence errors; `Mixture.state` ρ neutral in ≤ 110 µs (today 460), nine species in ≤ Stage 2's measured native case-A wall + 50 µs (today 5,250 µs; the
≤ 0.7 ms this first said is below Stage 2's own ≤ 1.0 ms native floor — Python cannot
beat the native it wraps); `Mixture.state` P nine species in ≤ 5 ms (today 380–1,690).
`grep` finds no `numpy` arithmetic in `src/epcsaft/*.py` outside array marshalling.

---

### Stage 4 — `equilibrium/` specification and minimizer (2,250 lines)

`problem.hpp` (300), `stoichiometry.{hpp,cpp}` (920), `request_compile.cpp` (600),
`gibbs_minimizer.{hpp,cpp}` (1,350) — and `problem.py` (450), `loader.py` (250).

This is #191 plus the core of #187. One objective:

    G = Σ_p [A_p(T, V_p, n_p) + P V_p] + Σ g° n     over (ln n_p, ln η_p)

Phase equilibrium, chemical equilibrium, and the coupled problem are the same KKT system
with `B = I`, one phase, or both; the affinity vanishes because `Bν = 0`. Pressure
equality is a stationarity row, not an inner density solve.

`Problem` is the activation matrix the user specified, not a union of route classes: the
caller declares which variables are specified, and the equation set is forced.
`compile_phase_direction_system` is the well-posedness check — unknowns against
independent equations, failing by naming the deficiency. Fixed T,P; fixed T with an
incipient partner solving P; fixed P solving T; reactive; nonreactive as `B = I`;
electrolyte as the charge row; omitted topology as `Π` unknown — all fall out without a
branch. Six `Ipopt::TNLP` subclasses become one `GibbsMinimizer`.

One state evaluation per Ipopt iteration, served to the objective, gradient, constraint
Jacobian, and Hessian callbacks. Today each callback re-evaluates: ~10 tapes per counted
evaluation. This is worth 5–10× on the reactive solve on its own.

**Gate:** the well-posedness check accepts each of the seven patterns above and rejects a
deficient one by name; conservation and charge residuals ≤ 1e-10; the homogeneous
analytic reaction case exactly.

---

### Stage 5 — `phase_boundary` and saturation (550 lines)

`phase_boundary.cpp` (300), `implicit_sensitivity.cpp` (250).

Saturation is not a file. `solve_saturation` is the one-component, two-phase call of
`solve_phase_boundary`; `Pure.saturation` in Python calls that binding. The three Newton
systems in the old tree (`density.hpp`, `pure_saturation.cpp`, `incipient_closure` +
`coexistence.cpp`) are rows of the one system, selected by which variables are unknown.
The one exception kept as distinct numerics: the 1-D density root keeps its bracketed
scalar solver, because multiple roots need branch selection.

The 21-start grid (7 vapour × 3 liquid) is not written. Two or three seeds from a
Clausius–Clapeyron estimate, and continuation from the previous point when the caller is
sweeping — the change that took the existing fit from 3,509 ms to 330 ms with identical
fitted parameters.

**Gate:** oracle case `light_alkane_saturation` (10 rows: ethane, propane, n-butane
against NIST) to the recorded pressures and molar volumes; `Pure.saturation` cold in
**≤ 1 ms** (today 454) at **≤ 25 evaluations** (today 5,422); bubble/dew Jacobians
against centred differences.

---

### Stage 6 — `phase_stability.cpp` (1,100 lines, range 950–1,300)

`phase_stability.hpp` (150), `phase_stability.cpp` (1,100). This is #188, and the
findings in `held2-stage-i-ii-findings.md` are its specification. Written to Pereira 2012
(HELD = *Helmholtz free Energy Lagrangian Dual*), replacing 10,538 lines of HELD2 plus
3,051 of HELD v1.

The findings, as build instructions:

- **One electroneutral chart** (finding 4). The number of dependent amount coordinates
  per phase is one closure plus one per charge row; the charged chart is the neutral
  chart plus one elimination. The neutral case has no eliminated coordinate because it
  has no charged species — so there is no `eliminated_index == C` sentinel, and the two
  defects that came from six consumers testing it and two forgetting cannot occur.
  850 lines → 60.
- **No pressure envelope** (finding 3). Reference and bounding points come from the
  density owner's root. 450 lines → 0, and 1,900 evaluations per flash → tens.
- **Stage I descends, it does not sample** (finding 2). A fixed start list and local TPD
  minimization from the best samples, with `search_budget_exhausted` as a named stop
  rather than `MAXEVAL_REACHED` standing in for convergence.
- **The inner problem is the Gibbs minimizer** (finding 5) with one phase and fixed
  multipliers, plus `−λ·n`. Its own reduced-Hessian and stationarity residuals are the
  certificate; the 570-line re-derivation, the 250-line KKT audit, the Gaussian-
  elimination rank, the dilute-face restart, the 12-iteration pressure polish, and
  `same_coordinates` are not written. 1,984 lines → 150.
- **The paper's start rotation from the first start** (finding 6): random / shifted-
  previous / near-pure alternating, not `random_interior` with the rotation reachable
  only in recovery. The paper credits shifted-previous with up to 50% cost reduction.
- **No second search, no trace scan** (finding 7): on every measured feed the first
  search's witness was accepted. 330 lines → 0.
- **HiGHS's `kOptimal` is the LP certificate** (finding 8). 90 lines → 0. The LP model
  assembly follows problem OP exactly and carries over.
- **One `Discovery` result** (finding 9): witness, candidates, counts, named stop.
  Fourteen result structs, 30 reason strings, and a per-invocation timing struct — from
  which six scalars were exported — become one type.
- **Candidates go to the Gibbs minimizer** rather than a from-scratch Stage III re-solve
  (finding 1), and `ε_b` is the paper's `1e-2`, set in Stage 1. On the Khudaida feeds
  Stage III took 14–42% of evaluations and had 60–90% of its solutions rejected by the
  Eq. 12 test.

**Gate:** the neutral binary at 300 K, 1e5 Pa, feed (0.5, 0.5) in **≤ 150 evaluations**
and **≤ 150 recordings** (today 4,631) and **≤ 10 ms** wall (today 671); the eight
water/ethanol/isobutanol/NaCl feeds each in **≤ 600 evaluations**; the Khudaida LLE feeds
reproducing their phase count and compositions within Table 3 uncertainties, in **≤ 29
major iterations** (the paper's ternary figure; today 18–38 with 4–8× its evaluations).
A decision that meets the count but changes the phase count fails.

One Khudaida feed enters the default gate (finding 11): today the only gated HELD2 test
is Stage I on a one-phase electrolyte, so the kernel's two-liquid result has never been a
gate claim.

---

### Stage 7 — `regression/` (890 + 200 lines)

`least_squares.{hpp,cpp}` (490), `eos_observable.cpp` (150),
`equilibrium_observable.cpp` (250), `regression.py` (200).

This is #192: `fit(model, observations)` with covariance, no owner compilation. Ceres
kernel, `ResidualBlock`, `PredictionCache`. The saturation continuation from Stage 5 is
what makes a parameter fit affordable: 140,562 recordings for a three-parameter
saturation fit today, because each of 26 residual evaluations is a cold 21-start solve.

**Gate:** NIST parameter recovery on the `light_alkane_saturation` rows; the Jacobian
against centred differences; permutation invariance; the three-parameter saturation fit
in **≤ 400 ms** (today 11,705) at **≤ 200 recordings** (today 140,562).

---

### Stage 8 — Composite gate and cutover

The reactive case, which touches every layer, is the last gate — it is the one problem
that exercises finite electrolyte liquid, incipient neutral vapour, two reactions, and a
solved pressure at once.

**Gate:** the four `vle_obs_*` reactive bubble points in **≤ 250 ms each** (today
7,390–11,400) at **≤ 300 recordings** and **≤ 300 Hessian calls** per solve (today
1,700–1,830 and 11,000–11,500), reproducing `oracle/reactive_vle_baseline.json` liquid
molar densities (46,425–47,148 mol/m³) to the recorded precision.

Then #196: migrate `analyses/`, the MEA and Lithium_Extraction repositories, result field
names, and packet names. Both downstream repos are already unrunnable against `main`
(MEA: 123 requests through the retired Provider; Lithium_Extraction: `ProviderModel` gone
and `dipole_moment` rejected per ADR 0006), so this stage repairs them rather than
preserving them — the migration was owed regardless of the rebuild.

Cutover is one merge of `greenfield` into `main`.

## Units: C++ resolves them, Python only reads the file

The whole vocabulary is 14 strings across every parameter document in the project:
`dimensionless` (1152), `kelvin` (503), `angstrom` (444), `elementary-charge` (113),
`kg / mol` (80), `kilogram / mole` (42), `1 / kelvin` (28), `gram / mole` (25), `pascal`
(20), `bar` (14), `1 / kelvin ** 2` (4), `megapascal` (2), `natural_log`,
`absolute_mole_fraction`. Three spellings of molar mass, two of pressure — the failure mode
is real (a document in g/mol read as kg/mol is wrong by 1000) and the table is tiny.

**The split: Python does structure, C++ does numbers.** Python walks the JSON — I/O, not
arithmetic — and emits `(identity, magnitude, unit_string)`. The native parameter builder
resolves each unit against the field's expected dimension, rejects an incompatible one by
name, and fills the array.

Why not Pint in Python, and why not a Python-side table:

- **It is validation, and validation happens once, natively.** A wrong unit is the same
  class of error as a negative segment count, which the Stage 2 constructor already rejects
  by field name. The tree this replaces ran three validation passes in two languages; one
  pass in one language is the point.
- **Python stays arithmetic-free with no exception clause.** A unit multiply in
  `parameters.py` is small, but it is the first crack in a rule whose value is that it has
  none.
- **It is testable without Python.** The same `CHECK` harness as the rest of parameter
  validation, and Stage 3's are the first gates that need a Python install at all.
- **Pint leaves the hot path.** `Mixture.state(T=313.15*ureg.kelvin, ...)` unwraps a
  quantity per call today, inside the measured ~590 µs of Python overhead on a call whose
  native work is 67 µs. The public API takes and returns SI floats: K, Pa, mol/m³, mole
  fractions, J/mol.

Cost: about 40 lines of table and dimension check in `model/`, and `Pint` leaves
`engine/pyproject.toml`. It is not installed in the greenfield venv today, so any Stage 3
gate importing `epcsaft` would have failed on `units.py` regardless.

## The MEA parameter set is the living bundle, not an estimation input

Stage 2's nine-species fixture was transcribed from
`analyses/2026-mea-parameter-estimation/inputs/neutral-selected-parameters.json`, an older
estimation input. The current set is
`MEA-Thermodynamics/analyses/mea_parameter_bundle/results/selected-current-best-parameters.json`,
tabulated in that analysis's `results/notebook.pdf` page 3 with a provenance class per
value (1 confirmed foundation, 2 literature-backed candidate, 3 local fit awaiting
promotion, 4 missing or default). Three values differ:

| | Stage 2 fixture | current bundle |
|---|---|---|
| CO₂ ε/k | 169.21 | **173.44025** |
| CO₂–water k_ij | absent | **0.013262879176919628** |
| permittivity formulation | ion suppression, per-ion α | **solvent-only** (α = 0 everywhere) |

Stage 2 is unaffected: its checks are closed forms, exact symmetries, limits and derivative
agreement, which hold for any parameter set — which is what `TESTING.md`'s rule buys. But
every MEA claim from Stage 3 onward reads the living bundle, and the fixture is refreshed
when Stage 3 lands.

**Two different roles, not two candidates.** Baygi & Pahlavanzadeh 2015 is a *published
model reproduced*: PC-SAFT, MEA–H₂O, no electrolyte, its own reported AAD, and therefore a
pass/fail gate — a fact about the world. The living bundle is the *model under development*:
nine species with electrolyte and reactions, currently 2.26× off in rmse-factor on pCO₂, so
its deviations are reported and never asserted. Baygi gates the binary boundary rows at 5a;
the bundle is what every reactive case is evaluated with.

## Order of the performance work, corrected

The audit put tape caching first. The measurements put it last. Written into the stages
above, in the order the levers actually matter:

1. Hessian count per evaluation, `(S+1)` → 1, plus site elimination — Stage 2, ~10×.
2. Tapes per counted evaluation, ~10 → 1, by one state evaluation per Ipopt iteration —
   Stage 4, 5–10× on the reactive solve, no derivative code involved.
3. Evaluation count in discovery, 4,631 → ≤150 — Stage 6, 30–100×.
4. Evaluation count in saturation, 5,422 → ≤25, and continuation in fits — Stage 5.
5. Recording — Stage 2, ≤10%.
6. Python marshalling, `Mixture.state` only, 590 µs per call — Stage 3.

The remaining ~10× to the compiled-tape endpoints (CppADCodeGen) is orthogonal to the
rebuild and is not in scope. It stays available because the tape is one object with one
owner.

## Naming

No proper nouns in identifiers. HELD, Held2, GREPE, Provider, and the step numbering do
not appear in the greenfield tree; `find_phases`, `solve_equilibrium`,
`solve_phase_boundary`, `GibbsMinimizer`, `Problem`, `EquilibriumResult` name what the
code does. Source names appear as citations in docs and `// EqID:` annotations, never as
identifiers. HELD in particular reads as the wrong author inside this repository's own
bibliography, where Christoph Held is a parameter source (`Ji2012`, `Figiel2025`,
`AscaniHeld2021`) and the algorithm is Pereira et al. 2012.

## What could invalidate this plan

- **The one-Hessian equivalence.** If `a + λ·H` does not reproduce the `(S+1)` tensor
  elementwise to 1e-10, Stage 2's 7× is not available and the whole floor moves.
  Falsifiable at Stage 2's gate, before anything depends on it.
- **The ~10 tapes per evaluation.** Inferred from 1,759 recordings / 183 evaluations. If
  the extra tapes are distinct states (line search, vapour closure) rather than
  re-evaluations of one, Stage 4's 5–10× shrinks to 2–3×. Instrument the callback
  sequence at Stage 4 before relying on it.
- **The reduced start set.** If 40–150 evaluations cannot decide a binary, or the eight
  LLE feeds find the wrong phase count from a fixed start list, Stage 6 keeps DIRECT-L
  (+80 lines for the NLopt adapter) and its budget becomes a named, justified stop rather
  than an accident.
- **`stoichiometry.cpp`'s standard-state transform** reads EOS reference quantities
  today; its line budget (600–900) is the least certain in the tree.

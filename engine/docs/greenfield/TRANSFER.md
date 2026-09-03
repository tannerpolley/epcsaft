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
| #195 layout and gate | move `_native/` → `native/`, reconcile three archives to one, fix the checker's contradicted module ban | the tree shape and the four checker rules | Stage 1, ~180 lines (checker + meson) |
| #186 numerical policy | find 15 `R` definitions, 12 duplicated names, unify | `physical_constants.hpp` + `numerical_policy.hpp` written before any consumer | Stage 4, 280 lines |
| #189 EOS derivatives | replace `(S+1)` Hessians, add tape cache, Newton association | one Lagrangian Hessian, site elimination, cached `ADFun` | Stage 4, 350 lines |
| #190 EOS Python | strip numerics from Python, cut exports to 11, replace `py::dict` round trips | Python as a pure mapping layer, `py::class_` bindings | Stage 3, 1,100 lines |
| #191 Problem/Result | collapse 6 specification types and 2 request families | one degrees-of-freedom `Problem`, one `EquilibriumResult` | Stage 5, 1,050 lines |
| #187 GREPE ownership | delete `coexistence.cpp`, the declared routes, `solve_incipient_pressure` | one `GibbsMinimizer`, the boundary Newton as its rows | Stage 5, 1,650 lines |
| #188 phase stability | delete HELD v1 (3,051) and HELD2 steps 8–10 (3,338); fix the neutral-chart sentinel | `find_phases` written to Pereira 2012 | Stage 7, 1,100 lines |
| #192 regression | remove owner compilation from the fit path | `fit(model, observations)` on Ceres | Stage 6, 890 lines |
| #196 consumers | migrate analyses and two downstream repos, rename result fields | the same migration — this one does not change | Stage 8, unchanged scope |
| #181 policy-edit gate | reject unreviewed solver-policy edits | folds into Stage 4's duplicate-constant gate | — |
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

## Execution order

One number per stage. There is no separate "run" and "stage" numbering any more: the
sequence below *is* the identifier, and the old numbering (which kept `TRANSFER.md`'s first
draft order as names) is retired. §"Crosswalk" maps the old labels for the three closed
stages, whose commit messages still carry them.

The order is the include direction (`physics ← model ← equilibrium ← regression ←
bindings`) with two deliberate departures: bindings land early so every later stage is
testable from Python against the oracle, and saturation is the first equilibrium stage —
the first thing built on the rows, before the flash. No workload is prioritised; the engine
is built inside-out.

Every stage owns exactly one brief pair, `briefs/stage-NN.brief.{md,json}`, written by the
planner and preflighted before a worker starts. Where a stage needs a design document
larger than a brief, it is named here and is the brief's specification.

| # | stage | category | writes | design document | state |
|---|---|---|---|---|---|
| 1 | **Canon** | canon | the layer tree, one archive, the checker's rules, and what may be tested | `CANON.md` | **closed** `e4d964c` |
| 2 | **Core** | core | `native/core/` — seven headers below everything: the four `EqID`-free tenants move in, and the declared-choice contract, the convergence record and the trace are written | `CORE.md` | designed — ~1,300 lines, ~+560 net |
| 3 | **Surface** | surface | `native/surface/` the parameter builder; `native/bindings/` the only pybind11 includes; `src/epcsaft/` as a pure mapping layer, units resolved in C++ | `SURFACE.md` | **reopened** — merged, then three defects |
| 4 | **EOS** | engine | the equation kernels and the policy they read; one tape per topology, sites eliminated, one Lagrangian Hessian, the density closure | `EOS.md` | **closed** `8bc6419`, `a2a2ac5` — reopened with 3 |
| 5 | **Equilibrium** | engine | the specification and count rule; one Ipopt driver for every `k`; seeds, continuation, sensitivity; starts and the topology contract | `EQUILIBRIUM.md` | designed — ~3,040 lines |
| 6 | **Regression** | engine | `fit(model, observations)` on Ceres, with covariance | `REGRESSION.md` | planned |
| 7 | **Phase stability** | engine | `find_phases` written to Pereira 2012, replacing 13,589 lines | `PHASE-STABILITY.md` | planned |
| 8 | **Final check** | canon | every layer exercised at once; one merge of `greenfield` into `main` | `FINAL-CHECK.md` | planned |

**Stages are named for what they write, not for what they solve.** Corrected 2026-09-02
after the owner caught it: the table called Stage 5 *saturation*, Stage 5 *the flash* and
Stage 5 *phase boundary*. Every one of those is a problem kind — the `k = f = 1` cell at one
component, the `k = f = 0` cell, and the `k ≥ 1` family — and naming a stage after a cell
reintroduces in the plan the special case the design removed from the code. Stage 5's own
body already said it: *`saturation` is a rule-7 case token, so nothing under `native/` is
named after it* — printed directly beneath a heading that named the stage after it.

There is one driver and one equilibrium. A cell is a specification, never a stage.

**Equilibrium is one stage, decided by the owner 2026-09-02.** It was four — a specification
stage, a driver stage, a starts stage and a boundary stage. Under one Ipopt driver serving
every `k` (objective `G` at `k = 0`, a zero objective with every row as a constraint at
`k ≥ 1`, because Ipopt's primal step on a square feasibility problem *is* the Newton step
whatever Hessian it is given) there are no rows for a second stage to write: compile emits
them once and the driver solves them. The multicomponent seed is a start rule, and
continuation across a sweep is one function that does not change with component count, so
neither is a second capability. Every split was route thinking surviving in the plan after
it had been removed from the code.

Budget ≈ 3,040 lines, ceiling ≈ 3,870 — specification 1,000, driver and one-component cell
1,330, starts and topology 710. The parts are the equilibrium planner's estimates, marked
**guessed** with per-file reasoning; `STAGE-LOG.md` records the actual.

**Stage 3 is reopened, and Stage 4 has a defect.** Stage 3's nine gates pass and every
mutation was shown red; the disagreement it surfaced is not its own.
`analyses/reproduction/2025-figiel-engine-crosswalk` localised it to the like-ion dispersion
pairs and then settled it against the literature: ePC-SAFT sets the dispersion energy
between cation–cation and anion–anion pairs to zero (Held 2008 excludes all ion–ion
dispersion; Held 2014 adds back only cation–anion; Figiel 2025 states `u_ii = 0` outright),
and `native/physics/helmholtz.hpp` applies the full combining rule to every pair.

**This is a Stage 4 defect** — the dispersion kernel — surfaced by Stage 3, and the engine is
internally inconsistent about it: it honours `u_ii = 0` in the ionic packing diameter, whose
`0.88 σ` is *derived from* that rule, while ignoring the rule in the dispersion energy.
Stage 4's gates could not have caught it: closed forms, symmetries, limits and derivative
agreement all hold for any `u_ij` matrix, including a wrong one. A modelling convention is
seen only by a comparison against the model that defined it, which is an argument for the
oracle rather than against the property tests.

Stage 2 does not start until the kernel is fixed and both stages re-accepted.

**Saturation before the flash, and why it is Stage 5 and not Stage 3.** Pure saturation is
the `k = 1`, `f = 1` cell of the one row system: one phase pinned at zero amount, `P` free
(`EQUILIBRIUM.md` §2.3). It is two rows over the density closure Stage 4 already
provides and it needs no flash seed — a Clausius–Clapeyron estimate is enough. What it does
need is the row assembly and the compile step, because writing its own chemical-potential
residuals would be the second route the whole rebuild exists to remove
(`EQUILIBRIUM.md` §7.1). So it runs as early as it can: immediately after the
specification stage, before the minimizer, and Stage 5 adds the multicomponent rows to the
same file. Saturation is never reimplemented; it is the case the general rows reduce to.

**Why the specification and the minimizer are two stages.** They were one 2,250-line stage.
`EQUILIBRIUM.md` §9 splits them at a seam it identifies itself: the compile tier's
checks C1–C8 are green on the compile step alone, before any solver exists, and it
recommends writing them first for exactly that reason. Splitting there gives Stage 5 a gate
that needs no driver, and lets saturation sit between the two.

**Stage 2 is a core, not a layer.** `CORE.md` is its specification. It sits between 4
and 6 because everything it adds to `physics/` and `model/` can be written and
mutation-tested against the tree that exists after Stage 3, and everything Stage 5, 9 and
10 must then do — the running maximum of the gradient floor, `κ_KKT` from the sensitivity
factorisation, the aggregate fit floor — is a commitment those briefs inherit rather than a
mechanism they invent. Its gate is the one the document names: a `derived()` tighter than
its input does not compile; every solver result carries `residual`, `tolerance`, `floor`,
`amplification` and a named stop; the old tree's `3e-10 > 2e-10` inversion is rejected by
the nominal floor at factor one.

### Outside the engine

Work that needs a brief, a gate and acceptance exactly like a stage, but is not a layer of
the engine. Lettered, not numbered, so the engine sequence stays unbroken.

| # | stage | writes | design document | gate |
|---|---|---|---|---|
| P | packet contract cut | one provenance form (`origin: {repo, commit, path, sha256}` per file); the lock manifest removed; the manifest's duplicated header dropped where the data file carries it; `validate_packets.py` and `data/docs/packet-contract.md` rewritten; the 14 packets migrated | `PACKET-CONTRACT-CUT.md` (to write) | every packet validates; every hash the old form recorded still resolves; the in-place edit found on 2026-09-02 is rejected by the new check |
| C | consumers (#196) | `analyses/`, the MEA and Lithium repositories, result field names, packet names | — | each consumer's own suite green on the greenfield wheel |

### Crosswalk

Stages are ordered by **category**: canon, then core, then surface, then the engine, then the
final check. That is dependency order — the core depends on nothing, the surface needs the
contract types, and every engine layer needs the core — and it is *not* the order the closed
work happened in. The number is the intended build order; the Status column carries what is
actually built.

| old | new | note |
|---|---|---|
| 0 | 1 | Canon |
| 1, 2 | 4 | physics and model are one thing: the EOS |
| 3 | 3 | Surface; the brief keeps its `stage-04` filename |
| T, 5 | 2 | Core; it was never numbered before the sequence existed |
| 4, 5a, 5b, 6 | 5 | Equilibrium is one stage |
| 7 | 6 | Regression |
| 6 | 7 | Phase stability |
| 8 | 8 | Final check |

## The line between `model/` and `equilibrium/`

Stated once, as a rule rather than a judgement call:

- **`model/` answers a state.** Give it `T`, `V` or `ρ`, and amounts, and it returns
  properties and exact derivatives. One state in, one answer out. The density root is the
  one iteration it performs, and it iterates over a single phase.
- **`equilibrium/` answers coexistence.** Anything requiring two states to satisfy a shared
  condition — equal chemical potential, equal pressure, a closed balance — lives here,
  whatever its component count.

**If it needs two states, it is not an EOS property.**

### Saturation is equilibrium, fully

Pure-component saturation is what pure parameters are regressed against, so it reads like an
EOS property, and routing it through `equilibrium/` makes `regression/` depend on the
equilibrium layer. Both observations are fair and the decision is still that saturation moves
fully into `equilibrium/`, as the one-component rows of `phase_boundary.cpp`. **Nothing
saturation-shaped stays in `model/`** — no seeds, no Clausius–Clapeyron estimate, no
vapour-pressure helper, no "quick" pure-component coexistence path. The seed generator
Stage 5 writes lives in `phase_boundary.cpp`, beside the rows it seeds.

1. **Physically it is coexistence, not a state evaluation.** Density, enthalpy and fugacity
   are answered from a single state. Saturation requires equal chemical potential and equal
   pressure between two phases. It resembles an EOS property only in what it is used for.
2. **A separate solver would be a second implementation of rows that already exist** — its
   own Newton, seeds, tolerance and failure modes, drifting from the general path. That is
   the special-case pattern rule 7 forbids, reached by a different route.
3. **The scientific argument, which is the strongest.** Shared rows mean a pure-component fit
   and a mixture fit share convergence behaviour and tolerances. A separate saturation solver
   lets a pure fit succeed under conditions a mixture fit would not, producing parameters that
   work only in the regime they were fitted in, with nothing in the result to signal it.

### The dependency is lighter than it reads

"Regression depends on equilibrium" overstates it. A pure-component saturation activates two
phases, one component, chemical-potential equality and pressure equality — and activates **no
reactions, no charge row, no stoichiometry, no stability**. `regression/` depends on
`problem.hpp` and a Newton with almost every row inactive. That is the degrees-of-freedom
design doing its job: the specification selects the rows, so a small problem compiles to a
small system rather than paying for a general one.

### This is why saturation runs at Stage 5

The order is not arbitrary and must not be "simplified" back. Saturation runs directly after
the specification stage and before the flash, so `phase_boundary.cpp` exists early with only
its one-component rows exercised, and Stage 6 gets what it needs without the rest of the
equilibrium stack ever having been built. The layering objection and the build order are
answered by the same choice.

### Single-phase is a legitimate specification

A degrees-of-freedom design must not assume `Π ≥ 2`: `{V}` alone and `{L}` alone are ordinary
specifications, a single-phase flash being a density solve plus stationarity. Verified against
`EQUILIBRIUM.md` on 2026-09-02: `{L}` and `{V}` are values of the phase-set axis
(§3.1), `{L}` is a solve-tier cell of the electrolyte fixture (§4.2), the count rule's
`k = f = 0` row says "`Π` phases free" with no minimum (§2.3), and S11 solves a one-phase
reactive cell. For a single free phase the balance rows fix the composition and the amount, so
the volume is what remains, and minimising `G` over it *is* the pressure-equality row — the
design degrades to a density solve exactly as it should.

A `Π ≥ 2` assumption anywhere in compile therefore turns C1 and C3 red on the `{L}` and `{V}`
cells, because the lattice walk emits them. That is structural rather than named, and a
named claim with its own mutation would be better; raised with the planner.

## Oracle

`engine/tests/oracle/oracle.json`, 71 rows in seven cases, measured columns only. Every
row is a fact about the world; no implementation can be right and fail one. Each case
names the analysis leaf that owns its evaluation (`analyses/reproduction/<year>-<author>`,
the repository's convention; a standard with no author keeps the standard's name).
Tolerances come in three kinds: a **row uncertainty** the source reports, a
**`derived()` constructor** the stage brief builds from the parameter set's own reported
deviation (rule 6 forbids the literal), or **reported, never asserted** where the parameter
set is the model under development and a threshold would test the fit, not the code.

| case | rows | system and property | conditions | closeness test | data | analysis |
|---|---|---|---|---|---|---|
| `light_alkane_saturation` | 10 | ethane; saturation pressure and saturated-liquid mass density, NIST SRD 69 | 100–280 K | `derived(saturation_reproduction_tolerance)` | `data/packets/nist-srd69-pure-saturation-properties/1` | `analyses/reproduction/nist-srd69` |
| `aqueous_alkali_halide_activity` | 24 | LiCl, LiBr, KCl, KBr in water; mean ionic activity coefficient (Hamer & Wu 1972) | 298.15 K, 1 bar, 0.001–2 mol/kg | per-salt AARD of γ± ≤ 10 % over the salt's six rows — the parameter source's own adequacy criterion (Figiel 2025 SI §S2); superseded engine gives 2.4–4.1 % | observations `data/packets/hamer-wu-1972-aqueous-alkali-halides/1`; parameters `data/packets/figiel-2025-reference-electrolytes/2` | `analyses/reproduction/1972-hamer-wu` |
| `aqueous_amine_binary_vle` | 6 | MEA–water; isobaric T–x–y (Cai 1996) | 373–443 K at 66.66 and 101.33 kPa | row uncertainties: ±133 Pa, ±0.1 K, ±0.001 x; carries the Baygi 2015 published-model reproduction, AAD 3.54 % x / 1.54 % y | row `source_file` + sha256 | `analyses/reproduction/1996-cai`, `analyses/reproduction/2015-baygi` |
| `loaded_amine_density` | 6 | 30 wt % MEA–water, CO₂-loaded; mass density (Amundsen 2009) | 298–353 K, 0–0.4 mol/mol | row uncertainty ±0.5 kg/m³ | row `source_file` + sha256 | `analyses/reproduction/2009-amundsen` |
| `aqueous_amine_carbon_dioxide_speciation` | 7 | MEA–CO₂–water, nine species; liquid speciation by NMR (Jakobsen 2005, Böttinger 2008) | 293–353 K, 0.11–0.99 mol/mol | reported, not asserted; structural gate: solves, conservation and charge ≤ policy | row `source_file` + sha256 | `analyses/reproduction/2005-jakobsen` |
| `aqueous_amine_carbon_dioxide_reactive_vle` | 10 | MEA–CO₂–water, nine species; CO₂ partial pressure at fixed loading (Jou 1995, Hilliard 2008, Xu 2011) | 313–393 K, 0.2–0.64 mol/mol | reported, not asserted; structural gate: converges two-phase, residuals ≤ policy, phase count | row `source_file` + sha256 | `analyses/reproduction/1995-jou` |
| `aqueous_alcohol_salt_lle` | 8 | water–ethanol–isobutanol–NaCl; both liquid-phase compositions on eight tie lines (Khudaida 2026) | 293.15 K, 5 wt % salt | row uncertainties per species: 0.0016 / 0.0007 / 0.0015 / 0.0002 mole fraction; phase count must match | `data/packets/khudaida-2026-figure-2-electrolyte-lle/2/observations/` | `analyses/reproduction/2026-khudaida` |

The rows do not name stages and stages do not name rows: a stage's gate says which
*capability* it makes available, and whichever oracle cases exercise that capability run
from then on. A case that a stage makes solvable stays in the run for every later stage.
Stage 3's brief resolved the activity tolerance above; the saturation constructor is owed
by Stage 5's brief and the analytic-reaction constructor by Stage 5's.

## Correction: the Stage 8 reactive gate is a golden file

`CANON.md` bans golden files of computed output — "a number that came out of this code is
not evidence that this code is right." The Stage 8 gate as first written pins the liquid
molar densities from `reactive_vle_baseline.json` (46,425–47,148 mol/m³), which are exactly
that: this engine's own output. The rule catches the plan that introduced it.

Replacement rows, all real measured data already transcribed with source hashes in
`MEA-Thermodynamics/data/reference/MEA/`:

| oracle case | source | gates |
|---|---|---|
| `aqueous_amine_binary_vle` | Cai 1996, isobaric T-x-y at 101.33 kPa | 9 — boundary rows on a neutral associating binary, no ions, no reactions; carries the Baygi reproduction |
| `aqueous_amine_carbon_dioxide_speciation` | Jakobsen 2005 / Böttinger 2008 NMR, 293–353 K | 4 |
| `aqueous_amine_carbon_dioxide_reactive_vle` | Jou 1995, Hilliard 2008, Xu 2011, 313–393 K | 9 |
| `loaded_amine_density` | Amundsen 2009, ±0.0005 g/cm³ | 3, 4 — the measured counterpart of the computed density the baseline pinned |

**These rows test the parameter set as much as the code**, and the current bundle is 2.26×
off in rmse-factor on pCO₂, so they cannot be pass/fail on agreement. Within `CANON.md`
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
Its parameters are *already* the Stage 4 associating-binary fixture — `2383.4744`,
`0.118488`, `k_ij = −0.052` in `strongly_associating_binary.json` **is** Baygi's model — so
the reproduction costs no new data. One item to read before asserting: whether Baygi's
Eq. 11 cross-association rule matches the fixture's Wolbach–Sandler volume
`0.1540467874409226`.

## Build order

Detail per stage, in the run order of §"Execution order". Each stage is one PR onto
`greenfield`, owns one brief pair `briefs/stage-NN.brief.{md,json}`, and closes when its
gate passes and every mutation has been shown to make it red. No stage begins before the
previous one's gate passes. Line counts are predictions; `STAGE-LOG.md` records the actual.

---

### Stage 1 — Canon

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
`solve_equilibrium`; it may not include `phase_stability.hpp`,
call `find_phases`, or construct a problem with phases omitted. A fit needs
differentiable residuals, and a discovered phase count is not one.

Delete in the same commit: `engine/src/epcsaft/_native/`, `engine/cpp/`,
`engine/tests/*/numerical/`, `engine/tests/*/behavior/`, and every Python module except
the three that survive (`CARRY-OVER.md` §4). The branch does not build after this commit
and is not expected to until Stage 3.

**Gate:** the checker reports every rule green over the emptied tree, and each rule
is shown to fail on a deliberate violation. Stage 1 landed five; rules 6 (no threshold
literals) and 7 (solvers are case-blind) landed later and are held to the same standard. `meson setup` does not configure until Stage 4
supplies the first sources — a Meson archive with no source files is an error — so the
build gate begins at Stage 4.

---

### Stage 2 — Core (~1,300 lines written, ~+560 net in `native/`)

`native/core/` — seven headers, no source file: `physical_constants.hpp`,
`numerical_policy.hpp`, `parameters.hpp` and `scalar.hpp` move in from `physics/`, and
`declared.hpp`, `convergence.hpp` and `trace.hpp` are written. `native/surface/` takes
`parameter_builder.{hpp,cpp}` from `model/`. `association.hpp`, `helmholtz_model.hpp` and
`density.cpp` gain their records; `epcsaft.trace()` reaches Python; the checker gains its
new rules and the layer ranks.

`CORE.md` is the specification, and `CORE.md` is superseded by it (§9 there is the disposition table). A nested stopping test cannot demand more accuracy
than its inputs carry: the outer residual's **floor** is the sum of the inner tolerances
amplified through their own conditioning at the solved point. The old tree tested a bubble
residual against `2e-10` while its two inner solves carried `3e-10` between them, so
eighteen *successful* closures were reported as iteration-limit failures.

Three parts. `derived(upstream, factor, why)` is `consteval` and refuses a tolerance
tighter than its input, making the inversion unwritable rather than discouraged. Every
solver returns one `Convergence` record — achieved `residual`, the `tolerance` it was tested
against, the `floor` its inputs allow, and the `amplification` that produced it — so a
declared factor is a hypothesis the runtime can falsify. `Stop` replaces "iteration limit
reached" with seven named outcomes, of which `converged_to_floor` is a success with a
stated accuracy.

The three floors this stage can already measure are all quantities the tree computes
anyway: the association Newton's pivot ratio, the tape's gradient floor from the Schur
product `Φ_yX Φ_XX⁻¹` already formed for the Hessian, and the density root's `1/|dr/dρ|`.
Stage 5 inherit the rest as commitments.

**Gate:** a `derived()` call tighter than its input fails to compile; the old tree's
`3e-10 > 2e-10` shape is rejected by the nominal floor at factor one; every solver result
carries the four numbers and a named stop; a solve driven below its floor reports
`converged_to_floor` rather than a failure.

---

### Stage 3 — Surface

`bindings/{module,eos_bindings}.cpp` (550), and `src/epcsaft/` `__init__.py`, `units.py`,
`parameters.py`, `eos.py`, `reference.py`, `_core.pyi`.

This is #190. Python maps inputs to `_core` and `_core` results to types, and does no
arithmetic — the rule the user set: *any Python frontend method call is completely
agnostic to what it is solving; it takes the inputs it is given, gives them to C++, gets
the answer back, and returns the output.* Bindings are `py::class_`, not `py::dict` round
trips. Eleven exports. Parameters are validated once, in the native constructor; the
three Python validation passes are not written.

**Gate:** oracle case `aqueous_alkali_halide_activity` (24 rows, LiCl/NaCl/KCl mean ionic
activity coefficients) reproduced through `Mixture.reference` within
`derived(activity_reproduction_tolerance)` — the rows carry `measured_...` and the
conditions, and nothing else; the "recorded convergence errors" this first named were a
computed column the oracle rebuild deleted. The brief constructs the factor from the
parameter set's own reported deviation on these salts; `Mixture.state` ρ neutral in ≤ 110 µs (today 460), nine species in ≤ Stage 4's measured native case-A wall + 50 µs (today 5,250 µs; the
≤ 0.7 ms this first said is below Stage 4's own ≤ 1.0 ms native floor — Python cannot
beat the native it wraps); `Mixture.state` P nine species in ≤ 5 ms (today 380–1,690).
`grep` finds no `numpy` arithmetic in `src/epcsaft/*.py` outside array marshalling.

---

### Stage 4 — EOS (2,760 lines closed; reopened on three defects)

#### 4.1 The equation kernels (950 lines)

`physical_constants.hpp` (30), `numerical_policy.hpp` (250), `helmholtz.hpp` (200),
`association.hpp` (150), `electrolyte.hpp` (260), `born.hpp` (60).

Mostly carry-over (`CARRY-OVER.md` §1): the Helmholtz terms, the association Newton
solver, Debye–Hückel, Born, and four of the five permittivity forms are correct and
validated, and are copied with their `// EqID:` annotations intact. The written work is
`numerical_policy.hpp`: every tolerance in the tree as a `Tolerance` entry with its
category (1 physical constant / 2 thermodynamic constraint / 3 algorithm tolerance /
4 unjustified legacy). Category-4 entries are not carried; each one is either justified
into category 3 with a source or replaced with the paper's value.

The HELD2 tolerance drift is resolved here, not in Stage 7: the candidate gap `ε_b`
enters as the paper's `1e-2` (not the drifted `5e-2`), and the composition floor as
`1e-8` (not `1e-10·f_k`). The seven density-root tolerances that HELD2 owned
(`root_*`, `mechanical_margin`, `stable_objective_tie`) enter as `model/density` entries.

**Gate:** closed-form, identity, limit, and permutation tests on the terms
(`tests/eos/`, ~200 lines). Checker rule 4 passes. No file in `physics/` includes CppAD
or a solver.

#### 4.2 The tape, the state and the density closure (1,810 lines)

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
threshold is real. This licenses the 7× and is the only thing in Stage 4 that can
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
  `reference_charge_square_sum == 0.0` branch (`:242`). Both flip inside Stage 4's own
  scope — the dilute end of the molality ladder takes the series, a pure-solvent state
  inside an electrolyte model takes the no-ion branch. A tape keyed on topology alone is
  silently wrong across them, with no error raised. Key on (active set, predicate
  signature) and re-record when it changes.
- **Association iteration counts need a three-line physics edit** — an optional
  `std::size_t*` out-parameter on `solve_site_fractions`. Stage 4 removed the kernel
  counters on the grounds that the tape owns evaluation accounting; without this the row is
  reported as unmeasured.

**Gate, performance (reported, not pass/fail — `CANON.md`):** the nine-species MEA model
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

*This defect was found by executing the plan's own claims against the tree before Stage 4
began, which is the mechanism described in `ORCHESTRATION.md`.*

---

### Stage 5 — Equilibrium (~3,040 lines, ceiling ~3,870)

One equilibrium, one driver, **one stage**. This was four — a specification stage, a driver
stage, a starts stage and a boundary stage — and every split was an artifact of route
thinking rather than a seam in the work. There are no rows for a second stage to write:
compile emits them once and the driver solves them for every `k`. The subsections below are
the order the stage is written in, not four gates.

`problem.hpp` (300), `stoichiometry.{hpp,cpp}` (600–900), `request_compile.cpp` (600),
`gibbs_minimizer.{hpp,cpp}` (870), `phase_boundary.cpp` (150), `implicit_sensitivity.cpp`
(250), starts and the topology contract (~340) — and `problem.py` (450).

#### 6.1 The specification (~1,000 lines)

`problem.hpp` (300), `stoichiometry.{hpp,cpp}` (920 → range 600–900), `request_compile.cpp`
(600) — and `problem.py` (450). No driver, no solve.

This is #191. `EQUILIBRIUM.md` is the specification. One row system, compiled from a
specification by one rule:

    rows − unknowns = k − f

with `k` pinned phases and `f` free intensives; square iff `f = k`, and rank at a generic
point is the sufficient check. The flash is `k = f = 0`; a boundary is `k = f = 1`; a
three-phase point is `k = f = 2`. No route is named anywhere — bubble and dew differ only
in *which* phase is pinned.

`Problem` is data with no method that solves. `request_compile.cpp` is the one file allowed
to branch, and only on status flags: it resolves the feed basis to amounts (natively —
mass-fraction and molality conversion is arithmetic, which Python may not do), builds `B`
from the reaction span by SVD, applies the count rule, assembles the active row list with
its `EqID`s, takes the rank, and either returns a `CompiledSystem` or throws naming the
deficiency.

**Gate:** the compile-tier lattice of `EQUILIBRIUM.md` §4.1, walked exhaustively —
several thousand specifications from four enums, none sampled. C1 count rule against a predictor
built from the specification alone; C2 every unreachable cell rejected by a message naming
the offending phase, reaction or species; C3 nullity 0 on every square cell and ≥ 1 where
all phases are pinned; C4 feed-basis invariance across four bases to 1e-12; C5
reaction-basis invariance under any invertible mixing; C6 the pinned specification's
residual equals the flash's at the same point with `N_q = 0`; C7 permutation to 1e-14; C8
an electrolyte flash compiles with nullity 0. Each with the named mutation from that table.

#### 6.2 The driver, and the one-component cell (~1,330 lines)

`gibbs_minimizer.{hpp,cpp}` — the one Ipopt driver, for every `k` (`EQUILIBRIUM.md`
§2.4: objective `G` or zero, constraint set from compile, exact Hessian for `k = 0` and
zero for `k ≥ 1`); `phase_boundary.cpp` opens with the one-component seeds and
continuation; `implicit_sensitivity.cpp` (250) lands with it. Revised 2026-09-02: the
separate boundary Newton is not written, because Ipopt's step on a square feasibility
problem *is* the Newton step, and a second solver loop with its own tolerance and stops
was the disease in another form. Budgets are in the section heading; `STAGE-LOG.md` records the actual.

Saturation is not a file, not a route, and not an identifier: `saturation` is a rule-7
case token, so nothing under `native/` is named after it and only Python's
`Pure.saturation` carries the word. Pure saturation is the `C = 1`, `k = 1`,
`f = 1` cell of `solve_equilibrium`, the one public entry: one phase pinned at zero amount, `P` free, over the
rows Stage 5 compiled. `Pure.saturation` in Python calls that binding. The three Newton
systems in the old tree (`density.hpp`, `pure_saturation.cpp`, `incipient_closure` +
`coexistence.cpp`) are rows of the one system, selected by which variables are unknown. The
one exception kept as distinct numerics: the 1-D density root keeps its bracketed scalar
solver, because multiple roots need branch selection.

The 21-start grid (7 vapour × 3 liquid) is not written. Two or three seeds from a
Clausius–Clapeyron estimate, and continuation from the previous point when the caller is
sweeping — the change that took the existing fit from 3,509 ms to 330 ms with identical
fitted parameters. **The seed generator lives in `phase_boundary.cpp`**, beside the rows it
seeds; nothing saturation-shaped is added to `model/`, including a vapour-pressure estimate
(§"The line between `model/` and `equilibrium/`").

If this stage writes its own chemical-potential residuals rather than calling Stage 5's
assembly, the rebuild has failed in the exact way it exists to prevent
(`EQUILIBRIUM.md` §7.1). The gate's S5 and S9 are what see it.

**Gate:** oracle case `light_alkane_saturation` (10 rows: ethane against NIST SRD 69) to
`measured_pressure_pa` and `measured_liquid_density_kg_m3` — the two columns the rows
actually carry — within `derived(saturation_reproduction_tolerance)`, constructed in the
brief from the parameter set's own reported deviation, not from what this engine happens to
return. S5 round trip: solve `P` at fixed `T`, then `T` at that `P`, and recover `T`. S9
sensitivity against centred differences to 1e-6. **Reported, not pass/fail:**
`Pure.saturation` cold in ≤ 1 ms (today 454) at ≤ 25 evaluations (today 5,422).

#### 6.3 Starts, and the topology contract (~710 lines)

`gibbs_minimizer.{hpp,cpp}` (1,350) — and `loader.py` shrinks, because the feed-basis
arithmetic it was budgeted for moved into Stage 5's compile.

This is the core of #187. One objective:

    G = Σ_p [A_p(T, V_p, n_p) + P V_p] + Σ g° n

Six `Ipopt::TNLP` subclasses become one `GibbsMinimizer`, the driver for `k = 0`. One state
evaluation per Ipopt iteration, served to the objective, gradient, constraint Jacobian and
Hessian callbacks; today each callback re-evaluates, about ten tapes per counted evaluation,
which is worth 5–10× on the reactive solve on its own.

The driver is Stage 5's; this stage writes the multiphase starts, the declared-topology
contract (§3.4 of the specification: whole or nothing) and the flash cells. The assembly is
Stage 5's, in `request_compile.cpp`, which rule 7 covers from Stage 5; that placement, not a
check, is what holds the driver and the implicit sensitivity to one matrix
(`EQUILIBRIUM.md` §5). S10 is a finite-difference agreement on each.

**Gate:** `EQUILIBRIUM.md` §4.2 solve tier on the flash cells. S1 every solved
cell's residual within `derived(minimizer_convergence_tolerance)` on every row, *including
rows the driver did not iterate on* — reaction affinity, and Euler's `g_p = λ·Bx_p` per
phase. S2 the reduced Hessian positive definite. S6 the instability certificate `G₂ < G₁`.
S7 lever rule to 1e-8. S8 extensivity under feed doubling. S10 the assembly's `F_y` against
centred differences of `F`, and Ipopt's `eval_h` against centred differences of the
Lagrangian gradient, both to 1e-6. S11 the ideal-gas limit `x_β / x_α = exp(−Δg°/RT)` on
the isomerisation. S13 the single-phase cell on every fixture. S15 a declared phase that
vanishes, or two that coincide, fails with no equilibrium state (`EQUILIBRIUM.md`
§3.4). S16 the reactive tie line on the Ascani 2023 `A + B ⇌ C` fixture is feed-invariant.
Oracle case `aqueous_alcohol_salt_lle` under **declared** `{L,L}` from each tie line's
midpoint feed, both endpoint compositions within the row uncertainties — the electrolyte
two-liquid cell under measured data before discovery exists (Stage 7 reaches the same
rows with topology omitted). The esterification tie lines of Senina 2022 and Schmitt 2006
(`EQUILIBRIUM.md` §4.2) become a further oracle case once transcribed.
Conservation and charge residuals ≤ 1e-10. Each with the named mutation from that table.

#### 6.4 The multicomponent cells — folded into 6.3

`phase_boundary.cpp` completes: seeds and continuation for `k ≥ 1` at more than one component. There are no new rows — compile emits them and Stage 5's driver solves them.

The same file and the same assembly Stage 5 opened, now with more than one component and
more than two phases: bubble, dew, cloud point and the fixed-fraction flash are one cell
family distinguished only by which phase is pinned and which intensive is free. Seeded from
the flash solution Stage 5 provides, or from the previous boundary point when sweeping.

**Commitment carried from `CORE.md` §6.** The boundary is written as *rows of the
KKT system*, never as an outer Newton over `P` around `solve_equilibrium`. That is the one
design choice that makes the old tree's failure mode disappear rather than reappear: an
outer Newton reintroduces exactly the nested-tolerance coupling of §1, and this brief would
then have to carry the minimizer's floor through it.

**Gate:** `EQUILIBRIUM.md` §4.2 on the boundary cells. S3 dew ↔ bubble: bubble from
liquid `x` at `T` gives `(P*, y*)`; dew from `y*` must return `(P*, x)`. S4 the
fixed-fraction flash at `β = 10⁻², 10⁻⁴, 10⁻⁶` approaches the bubble point's `T` at the
expected first-order rate. S5 and S9 extended to the multicomponent cells. Oracle case
`aqueous_amine_binary_vle` within its row uncertainties, carrying the published-model
reproduction. Each with the named mutation.

---

### Stage 6 — Regression (890 + 200 lines)

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

### Stage 7 — Phase stability (1,100 lines, range 950–1,300)

`phase_stability.hpp` (150), `phase_stability.cpp` (1,100). This is #188, and the
findings in `PHASE-STABILITY.md` are its specification. Written to Pereira 2012
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
  (finding 1), and `ε_b` is the paper's `1e-2`, set in Stage 4. On the Khudaida feeds
  Stage III took 14–42% of evaluations and had 60–90% of its solutions rejected by the
  Eq. 12 test.

**Gate:** the neutral binary at 300 K, 1e5 Pa, feed (0.5, 0.5) in **≤ 150 evaluations**
and **≤ 150 recordings** (today 4,631) and **≤ 10 ms** wall (today 671); the eight
water/ethanol/isobutanol/NaCl feeds each in **≤ 600 evaluations**; oracle case
`aqueous_alcohol_salt_lle` (8 tie lines) reproducing the phase count and both endpoint
compositions within the `composition_uncertainty_mole_fraction` each row carries, in
**≤ 29 major iterations** (the paper's ternary figure; today 18–38 with 4–8× its
evaluations). Each row's feed is the midpoint of its two measured endpoints, so the gate
depends on no digitized overall composition: every interior point of a tie line has the
same two equilibrium phases.
A decision that meets the count but changes the phase count fails.

One Khudaida feed enters the default gate (finding 11): today the only gated HELD2 test
is Stage I on a one-phase electrolyte, so the kernel's two-liquid result has never been a
gate claim.

---

### Stage 8 — Final check

The reactive case, which touches every layer, is the last gate — it is the one problem
that exercises finite electrolyte liquid, incipient neutral vapour, two reactions, and a
solved pressure at once.

**Gate:** the four `vle_obs_*` reactive bubble points each converge to a two-phase
solution with conservation and charge residuals ≤ 1e-10 and the phase count the chemistry
requires. **Reported, not pass/fail:** **≤ 250 ms each** (today 7,390–11,400) at **≤ 300
recordings** and **≤ 300 Hessian calls** per solve (today 1,700–1,830 and
11,000–11,500), and the liquid molar densities beside
`oracle/reactive_vle_baseline.json`'s 46,425–47,148 mol/m³.

Those densities are `measured_2026_09_02` — what the tree being replaced returned on one
day, not a fact about the world. Pinning them would make the old engine's answer the
definition of a right answer, which is the golden-file ban in `CANON.md`. They are a
change indicator: a rebuilt engine landing far outside that band is a conversation, not a
failure. The measured observable on these rows is `co2_partial_pressure_pa`, and the
bundle is the model under development, so its deviation there is reported too.

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
  class of error as a negative segment count, which the Stage 4 constructor already rejects
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

Stage 4's nine-species fixture was transcribed from
`analyses/2026-mea-parameter-estimation/inputs/neutral-selected-parameters.json`, an older
estimation input. The current set is
`MEA-Thermodynamics/analyses/mea_parameter_bundle/results/selected-current-best-parameters.json`,
tabulated in that analysis's `results/notebook.pdf` page 3 with a provenance class per
value (1 confirmed foundation, 2 literature-backed candidate, 3 local fit awaiting
promotion, 4 missing or default). Three values differ:

| | Stage 4 fixture | current bundle |
|---|---|---|
| CO₂ ε/k | 169.21 | **173.44025** |
| CO₂–water k_ij | absent | **0.013262879176919628** |
| permittivity formulation | ion suppression, per-ion α | **solvent-only** (α = 0 everywhere) |

Stage 4 is unaffected: its checks are closed forms, exact symmetries, limits and derivative
agreement, which hold for any parameter set — which is what `CANON.md`'s rule buys. But
every MEA claim from Stage 3 onward reads the living bundle, and the fixture is refreshed
when Stage 3 lands.

**Two different roles, not two candidates.** Baygi & Pahlavanzadeh 2015 is a *published
model reproduced*: PC-SAFT, MEA–H₂O, no electrolyte, its own reported AAD, and therefore a
pass/fail gate — a fact about the world. The living bundle is the *model under development*:
nine species with electrolyte and reactions, currently 2.26× off in rmse-factor on pCO₂, so
its deviations are reported and never asserted. Baygi gates the binary boundary rows at Stage 5;
the bundle is what every reactive case is evaluated with.

## Order of the performance work, corrected

The audit put tape caching first. The measurements put it last. Written into the stages
above, in the order the levers actually matter:

1. Hessian count per evaluation, `(S+1)` → 1, plus site elimination — Stage 4, ~10×.
2. Tapes per counted evaluation, ~10 → 1, by one state evaluation per Ipopt iteration —
   Stage 5, 5–10× on the reactive solve, no derivative code involved.
3. Evaluation count in discovery, 4,631 → ≤150 — Stage 7, 30–100×.
4. Evaluation count in saturation, 5,422 → ≤25, and continuation in fits — Stage 5.
5. Recording — Stage 4, ≤10%.
6. Python marshalling, `Mixture.state` only, 590 µs per call — Stage 3.

The remaining ~10× to the compiled-tape endpoints (CppADCodeGen) is orthogonal to the
rebuild and is not in scope. It stays available because the tape is one object with one
owner.

## Naming

No proper nouns in identifiers. HELD, Held2, GREPE, Provider, and the step numbering do
not appear in the greenfield tree; `find_phases`, `solve_equilibrium`,
`GibbsMinimizer`, `Problem`, `EquilibriumResult` name what the
code does. Source names appear as citations in docs and `// EqID:` annotations, never as
identifiers. HELD in particular reads as the wrong author inside this repository's own
bibliography, where Christoph Held is a parameter source (`Ji2012`, `Figiel2025`,
`AscaniHeld2021`) and the algorithm is Pereira et al. 2012.

## What could invalidate this plan

- **The one-Hessian equivalence.** If `a + λ·H` does not reproduce the `(S+1)` tensor
  elementwise to 1e-10, Stage 4's 7× is not available and the whole floor moves.
  Falsifiable at Stage 4's gate, before anything depends on it.
- **The ~10 tapes per evaluation.** Inferred from 1,759 recordings / 183 evaluations. If
  the extra tapes are distinct states (line search, vapour closure) rather than
  re-evaluations of one, Stage 5's 5–10× shrinks to 2–3×. Instrument the callback
  sequence at Stage 5 before relying on it.
- **The reduced start set.** If 40–150 evaluations cannot decide a binary, or the eight
  LLE feeds find the wrong phase count from a fixed start list, Stage 7 keeps DIRECT-L
  (+80 lines for the NLopt adapter) and its budget becomes a named, justified stop rather
  than an accident.
- **`stoichiometry.cpp`'s standard-state transform** reads EOS reference quantities
  today; its line budget (600–900) is the least certain in the tree.

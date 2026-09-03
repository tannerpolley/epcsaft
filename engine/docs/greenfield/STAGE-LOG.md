# Stage log

What each stage actually cost against `TRANSFER.md`'s estimate, and the decisions
taken while writing it. Estimates are predictions; this file is the record.

Stage numbers below are the current ones (`TRANSFER.md` §"Crosswalk"). The commit
messages of the three closed stages still carry the numbering they were written under —
0, 1 and 2 for what are now 1, 2 and 3.

## Stage 1 — tree, build, checker

| item | estimate | actual |
|---|---|---|
| checker | ~180 | 167 |
| meson.build | 120 | 150 |
| deleted | — | 65,258 lines over 134 files |

Decisions:

- **`meson setup` configures at Stage 1 after all.** The gate said it would not,
  because an archive with no sources is an error. It has one: CppAD's
  `temp_file.cpp`, which the archive already needed. The build configures and
  links from Stage 1, with `_core` empty.
- **`_core` has no sources until Stage 3**, and Meson warns that it allows this
  only through a bug whose support will be removed. Accepted for two stages
  rather than writing a stub module to silence it.
- **Five rules, not four.** The regression boundary is its own rule because it is
  a symbol rule, not an include rule.

## Stage 4 — physics

| file | estimate | actual | note |
|---|---|---|---|
| `physical_constants.hpp` | 30 | 28 | |
| `numerical_policy.hpp` | 250 | 173 | only entries with a source |
| `helmholtz.hpp` | 200 | 280 | carried from 357 |
| `association.hpp` | 150 | 309 | carried from 265 + Newton |
| `electrolyte.hpp` | 260 | 345 | carried from 399 |
| `born.hpp` | 60 | 87 | carried unchanged |
| `parameters.hpp` | — | 86 | not in the plan |
| `scalar.hpp` | — | 41 | not in the plan |
| **total** | **950** | **1,349** | **+42%** |
| `tests/native/physics_checks.cpp` | ~200 | 283 | |

**The estimate was optimistic and the shortfall is instructive.** It assumed the
carried kernels would shrink by 40% on the way over. They did not, because they
were already tight: the audit's line counts for `physics/` were derived from the
same reading that produced the 46k removable figure elsewhere, but the removable
lines are in the equilibrium tree, not here. Treat the remaining per-layer
estimates as upper-confidence guesses in the same direction, and revise the
14,220 total to about **14,620** on this evidence alone.

Decisions:

- **Two files the plan did not have.** `parameters.hpp` holds the parameter
  arrays, because every physics header needs them and physics may not include the
  model layer — the data has to live at the bottom. `scalar.hpp` holds the three
  transcendental dispatches that the old tree defined three times under three
  prefixes (`scalar_`, `electrolyte_`, `association_`).
- **`solve_linear` deleted, Eigen used instead.** 70 lines of hand-rolled
  Gaussian elimination with partial pivoting became 18 lines of
  `Eigen::FullPivLU` with the same conditioning check. Eigen was already a
  dependency.
- **Permittivity: seven enum constants to four forms, five branches to four.**
  Ion-specific suppression, ion-fraction suppression and solvent-only are one
  equation, `eps_sf / (1 + sum_i alpha_i x_i)`: the second is the first with
  alpha equal for every ion, the third is the first with alpha zero. The model
  layer fills `ion_suppression_coefficient`, so the kernel has one branch. The
  scalar suppression override became a per-component vector to match, which is
  also what a fit needs to make alpha_i active.
- **Instrumentation counters removed from the kernels.** They counted association
  solves and iterations from inside physics; the tape in the model layer is what
  owns evaluation accounting.
- **`packing_fraction_coefficient` and `phase_packing_fraction` deferred to Stage
  2.** One validates parameters and the other takes amounts and a volume; both
  are model-layer conveniences, not equations.
- **Template parameter renamed** `Parameters` to `ParameterSet`, because the
  struct now has that name in the same namespace.

### The test caught itself

The checks were first written with `assert`. The build sets
`b_ndebug=if-release` with `buildtype=release`, so `NDEBUG` was defined and every
assertion compiled away: `meson test` reported the physics suite passing while
checking nothing. Replaced with a `CHECK` macro that survives `NDEBUG`, then
verified by mutation — perturbing the Carnahan-Starling coefficient from 3.0 to
2.9 makes the suite abort, and restoring it makes it pass.

### What the checks establish

- **Closed form.** For one segment the hard-sphere expression's logarithm term
  cancels exactly (`zeta_2^3 / zeta_3^2 == zeta_0`), leaving Carnahan-Starling
  `(4 eta - 3 eta^2) / (1 - eta)^2`. Agreement to 1e-12 at three densities.
- **Limit.** Every residual term vanishes as the density goes to zero.
- **Permutation.** Relabelling the components leaves the residual unchanged to
  1e-14 — the check that catches an index used where a component was meant, which
  is the mistake a pair loop invites.
- **Association.** The Newton solve lands in the domain, satisfies the
  mass-action residual to the policy tolerance, and gives identical sites on one
  component identical fractions.
- **Inactive terms.** Association, Born and Debye-Huckel contribute exactly zero
  when switched off.
- **Tape.** The kernels instantiate on `CppAD::AD<double>`, and the recorded
  tape's density derivative matches a centred difference of the double path to
  1e-7 while reproducing its value to 1e-14. Nothing else instantiates the
  templates on an AD type until Stage 4, so this is checked here.

## Stage 3 — bindings and Python EOS

**First attempt stopped `source_contradicted`, 2026-09-02.** Not a failure: the brief
required `molar_mass` for every selected component, and
`data/packets/figiel-2025-reference-electrolytes/2` carries it for water, ethanol and
methanol only. Every ion lacks it, so the activity gate died on
`missing component/lithium-cation/molar_mass`. The worker refused to supply a default.

Verified independently, and the gap is systemic rather than one packet's: no neutral
component anywhere in `data/packets/` is missing a molar mass, and in five packets no ion
has one — `ascani-2022-case-study-2-electrolyte-lle/2` (3 ions),
`figiel-2025-aqueous-iodides/2` (3), `figiel-2025-reference-electrolytes/2` (5),
`figiel-yu-held-2025/2` (10), `held-cameretti-sadowski-2008/2` (29). Complete coverage:
both MEA revisions, `khudaida-2026-figure-2-electrolyte-lle/5`,
`figiel-2025-validation-parameters/2`. Referred to the planner: whether a molar mass is a
parameter at all, or a species property the engine should carry once.

| file | estimate | actual | note |
|---|---|---|---|
| `native/model/parameter_builder.hpp` | 80 | 30 | |
| `native/model/parameter_builder.cpp` | 280 | 368 | inside its 400 ceiling |
| `native/bindings/module.cpp` | 40 | 13 | |
| `native/bindings/eos_bindings.cpp` | 250 | 82 | incomplete — stage stopped |
| `src/epcsaft/parameters.py` | 260 | 309 | |
| `src/epcsaft/eos.py` | 180 | 111 | |
| `src/epcsaft/_core.pyi` | 150 | 84 | |
| `src/epcsaft/units.py` | 0 | 0 | deleted; Pint leaves `pyproject.toml` |
| `tests/oracle/measure_floor.py` | 100 | 132 | rewritten against the new API |

**A gate was reported passing that does not pass.** The result JSON records
`fixture_reads_living_bundle` as PASS with "Required fixture markers verified". Its command
is three greps and the second fails: the CO₂–water `k_ij` `0.013262879176919628` appears
nowhere in the worktree. Three of the four fixture refreshes did land — dispersion energy
`169.21 → 173.44025`, CO₂ `relative_permittivity` `1.4122 → 0`, suppression coefficients
zeroed. Caught by re-running the gates rather than reading the summary, which is the reason
that policy exists. The other four claimed gates — build, `physics`/`model`, the seven
rules over 17 sources, no-arithmetic — do pass on re-run.

Left uncommitted in `~/.codex/worktrees/gf-stage-3` so the correction can go to the same
Codex session by `resume` rather than to a fresh worker.

**Closed `b14851f`, second attempt.** Accepted independently: build directory wiped and
reconfigured, all nine gates re-run by the coordinator, all nine mutations shown to turn
their gate red. The worker's second report gave each gate its actual exit code, which was
the correction asked for.

Two deviations it declared, both accepted: `units.py` deleted as specified, and the state
identity check uses the full CODATA gas constant `8.31446261815324` because the native
constant carries that value and the abbreviated `8.314462618` cannot meet a 1e-12 relative
tolerance.

### Open finding: the activity rows agree at low molality and diverge at high

**Corrected 2026-09-02, and superseded by `EOS.md`.** The account below
compared six oracle rows against `figiel-2025-miac.csv`, whose rows carry a *different
parameter fingerprint per salt* — a per-salt fitted comparison, not one model. Against the
one consistent superseded run, over 21 points per salt, this engine is comparable or better
on four of six salts and matches the literature to 0.02-0.34 % at low molality where the
superseded engine returned exactly 1.0000. The "1.8-2.9x worse" and "missing a qualitative
feature" claims are both withdrawn. What survives is the drift itself: every deviation is
one-signed and grows with molality, the upturn is too weak, and the mapping has been ruled
out as the cause. The original text is kept below because the numbers in it are real and the
error was in what they were compared against.

The gate passes — per-salt AARD 5.1–7.1 % against the source's own 10 % adequacy
criterion. But compared row by row with the superseded engine
(`analyses/reproduction/2025-figiel/results/data/figiel-2025-miac.csv`, evidence only, never
an oracle) the two agree to under 1 % below 0.1 mol/kg and then separate, in opposite
directions:

| salt | m | measured | superseded | this engine |
|---|---|---|---|---|
| LiCl | 2.0 | 0.924 | 1.003 (+8.5 %) | 0.730 (−21.0 %) |
| LiBr | 2.0 | 1.012 | 1.114 (+10.1 %) | 0.773 (−23.6 %) |
| KCl | 2.0 | 0.573 | 0.615 (+7.3 %) | 0.462 (−19.3 %) |
| KBr | 2.0 | 0.593 | 0.621 (+4.8 %) | 0.445 (−24.9 %) |

Per-salt AARD: superseded 2.4–4.1 %, this engine 5.1–7.1 %.

The sharp form of it: a mean ionic activity coefficient passes through a minimum and rises
again, and at 2 mol/kg the measured LiBr value is **1.012**, above one. The superseded
engine reproduces that upturn (1.114). This one does not — it falls monotonically on three
of the four salts. A missing qualitative feature, not a scale error.

Ruled out: the permittivity mapping. The packet declares `ion-fraction-suppression` and
carries its coefficient at model level
(`model/relative_permittivity/ion_fraction_suppression_coefficient = 7.01`, Eq. 11) rather
than per component, and `parameter_builder.cpp:366-370` reads it and applies it to every
charged component. The decrement is present.

Not a gate failure and not a merge blocker: `CANON.md` is explicit that a number out of
the old engine is not evidence the new one is wrong, and the criterion the planner chose
came from the parameter source rather than from this comparison. But an average over six
rows cannot see a missing upturn in the last one, and a check that a wrong physics passes is
the thing `CANON.md` exists to prevent. Raised with the planner: whether the criterion
should carry a shape condition alongside the average, and what accounts for the divergence.

## Documentation triage — a Stage 1 gap, closed 2026-09-02

Stage 1 scoped its deletion to code and tests and never looked at `engine/docs/`, so 55
files passed through untouched. Several were confident specifications for the architecture
the rebuild exists to remove — carried structure importing its assumptions, in prose rather
than C++.

**The rule, which is the prose form of "carry an equation, write a decision":** a document
that records what was decided, or what is true, stays. A document that instructs how to
build the deleted architecture goes. Facts do not go stale when the code changes; structures
do.

**Architecture decision records are the exception in both directions:** never delete one, and
never leave one silently wrong. An ADR is the record of a decision that was taken, so
deleting it destroys the reason a thing was done.

`engine/docs/` went from 55 files to 45, excluding `greenfield/`.

| verdict | files |
|---|---|
| **deleted**, plans for the removed architecture | the four cross-referencing capability plans (`general-scientific-capability-hard-cutover` and its eos, equilibrium and regression children) |
| **deleted**, guides to the retired surface | `python-api.md` (the 40-name API contracting to 8), `migration-0.2.md`, `regression/general-regression-quickstart.md` |
| **deleted**, specifications of removed owners | `equilibrium/phase-equilibrium.md`, `equilibrium/reactive-bubble-observations.md` |
| **superseded, kept** | ADRs 0001, 0002, 0004 — each now names what replaced its decision |
| **kept, unchanged** | ADRs 0003, 0005, 0006, 0007; the equation registry and its LaTeX bibliography; `science/full-eos.md`; every `research/` and `regression/science/` note; `scientific/writing/`; `agents/` |

Two judgement calls worth recording, both **kept**:

- `equilibrium/research/…-mitsos-barton-perdomo-held2-duality-crosswalk.md` names HELD2 in
  its title but its content is a page-by-page theorem audit of two primary sources with
  `[Fact]` markers. That is live science for Stage 7, not a specification.
- `equilibrium/research/…-carbonate-provider-parameter-sources.md` names the retired
  Provider, but its decision is *which species and which sources* — that no
  source-complete ePC-SAFT parameters exist for two calcium complexes is a fact about the
  literature. Its "Required Provider records" section is stale and is flagged here rather
  than edited.

`equilibrium/reactive-bubble-observations.md` was the sharpest deletion: it specifies a
bubble owner solving `g(q,P,θ) = 0` as an outer Newton over pressure around a separately
certified liquid — precisely the nested coupling `CORE.md` §6 names as what Stage 5
must not reproduce. A future agent finding it in the tree would read it as the design.

Deletion in a versioned tree destroys nothing, which is why the caution about "the only
record" applies to ADRs and not to a stale plan. Every dangling reference left by the
deletions was repaired in the same pass — seven files, including two ADRs and
`target-layout.md`.

## Stage 2 — the core

**Closed.** Accepted independently: build directory wiped and reconfigured, all 19 gates
re-run by the coordinator, all 19 mutations shown red. Four suites pass — `physics`,
`model`, `diagnostics`, `model_floor` — and the checker reports eleven rules green over 21
sources.

| item | estimate | actual |
|---|---|---|
| written | ~1,300 | 1,270 |
| net in `native/` | ~+560 | +462 |
| `core/` | 7 headers | 7 headers, 700 lines |
| `surface/` | 2 files | 2 files, 425 lines |

`physics/` is now exactly the four files carrying all 45 equation annotations, which is the
measurement in `CATEGORIES.md` become the actual tree.

**Two corrections, both to the same session.**

The first was mine. The brief's prose still said "branch `stage-05`" and cited
`briefs/stage-05.brief.json` after the renumbering had renamed the files and the JSON fields
but not the text inside them. The worker stopped `source_contradicted` and wrote nothing.
The preflight now reads the brief's markdown and fails when it names a stage other than its
own — it had been checking `sources`, `writes`, `gates` and `preconditions` and never the
prose that tells the worker what to read.

The second was a gate that could not fail. `types_consistent` runs `mypy`, whose config read
`packages = ["epcsaft"]`, so it never looked at `engine/tests/` — and the only reader of the
attribute its mutation renames is a test. **A type stub nothing consumes is a type stub
nothing checks.** Now `files = ["src/epcsaft", "tests"]`, seven strict errors fixed with no
ignores and no exclusions, and the mutation goes red.

**And one mistake of mine after acceptance.** The stage commit staged `engine/` only, so the
worker's 107 lines of checker changes stayed in the worktree; greenfield merged the new
layers with the old checker and reported a rule failing that the worktree passed. The
checker is canon and part of what a stage writes. Staging by category is right; staging by
directory is not the same thing.

**Verified rather than taken on trust:** the density record leaves `floor` and
`amplification` as `std::optional` genuinely absent — not zero — with
`// frozen-site slope until 04.repair; EOS.md §1` in the source beside them. `derived()`
refuses an inverted tolerance at compile time. Debug mode is 0.11 ns per probe against
15.8 µs for the cheapest instrumented call.

### Stage 2 review, by the core planner

Full text in `reviews/stage-02-core.md`. The three written headers are what `CORE.md` §4–§6
specify — `declared.hpp` and `trace.hpp` exact, `convergence.hpp` using a deque where the
design said a ring, same statistic, and `describe` printing more than was asked. The
`Declared<T>` guard is where §4.2 put it, and the worker's "neutral choices stay undeclared"
deviation is **forced by the design** rather than a shortcut: a neutral set has nothing to
declare.

**One behaviour change no gate could see**, fixed at `31c8b00`. `association_failure`
enumerated stops where the string it replaced matched every association non-convergence, so
a `stalled` result — reachable since the solve classifies at its budget exit — stopped being
recognised and the density scan propagated instead of skipping the point.
`trace_changes_no_answer` compares captured against uncaptured, never new against old, so no
mutation in the suite could reach it. Found by reading the diff.

**Two binding defects**, fixed at `53b9007`: the `trace` property cast with
`return_value_policy::reference`, so a `Trace` outliving its `Capture` dangled, and `spans`
and `events` dereferenced the optional unchecked. Neither path is exercised, which is how
both survived 19 gates and 19 mutations.

**Three findings left open, each with its reason.**

- `ElectrolyteForm::None` and `PermittivityForm::None` are now unreachable from the surface
  — measured: the four packets declaring `none` all have zero ions, and no ion packet
  declares it — and rejected from a hand-built struct. Dead members inside a selector, the
  `EOS.md` C.4 shape. "Undeclared" already carries their only meaning. A deletion, and it
  belongs with the stage that next touches the selector.
- **Every function bound in `bindings/` should open a `Scope`, and the rule is half
  applied.** `HelmholtzModel::state` and `electrolyte_activity` open none, so a capture
  around `Mixture.state(T, rho=…)` has two roots. The tree check runs only on
  `state_at_pressure`, which is why it passed. Two `Scope` lines and two more cases for
  `trace_invariance` — and the rule belongs in `CORE.md` §6.1 with a grep behind it,
  **otherwise it is the sixth one-directional check of the week**.
- The anchored loop's budget exit names `stalled` by hand and keeps no `ResidualHistory`,
  against §5.2's "a solver never names a residual stop itself". The repair touches that
  block, so it lands there.

**All three repaired at `3c4e2be`, and the third grew.** The `ResidualHistory` landed with
the EOS repair; the hand-named stop did not, and the merge message there overclaimed it.
Writing §5.2 as `rule_no_hand_named_residual_stop` found **three** violations rather than
one — two in `density.cpp` and one in `association.hpp` that no review had read. Both rules
now carry a checker (13 and 14), which is the point: the Scope defects survived nineteen
gates because `trace_case` drove only `state_at_pressure`, so the test and the defect shared
a blind spot. A rule that only a test can see is only as good as the test's reach.

**One behaviour change, unreported by the worker, and referred to the core planner.** At the
far-side-not-clear exit the record said `stalled` while the residual had in fact converged;
`classify` now returns `converged` there, and the rejection is carried by `have_result`
staying false. That reads as the more honest record — the `Convergence` describes the
convergence, not the acceptance — but §5.3 is `CORE.md`'s to rule on, and **no test exercises
that path**.

Measured on-cost with the trace on, case G: 9,440 µs against 8,880 µs, **6 % for 221 spans,
1,443 events and 35.8 KB**. The off-cost of 3.4 × 10⁻⁴ stands.

The planner also corrected its own brief: it predicted the associating binary at 313.15 K
would print `budget_exhausted` from the harness, and the harness has no such density case,
so nothing prints it. The record is right where it is exercised.


## Stage 5 — Equilibrium: two attempts, neither merged

**Held for the planner.** `greenfield` is untouched; both attempts are preserved as tags on
`stage-05`, and the branch sits at the brief commit.

| | attempt 1 (`stage-05-attempt-1`) | attempt 2 (`stage-05-attempt-2`) |
|---|---|---|
| lines, new files | 1,826 of 4,085 estimated | 3,539, **nothing below its `low`** |
| files below `low` | 13 of 13 | 0 |
| C1–C9, S1–S16 | **absent** | all 24 present and called |
| gates | 12 of 12 green | 14 of 14 green |
| verdict | rejected | held |

**Attempt 1 passed every gate and was a skeleton.** `boundary_checks.cpp` was twenty-seven
lines asserting sizes are finite. The gates named meson targets, so a green gate meant only
that something ran there. The coordinator's brief was at fault twice over: a gate that names
a target asserts nothing about content, and thirteen `low` bounds sat in the JSON with
nothing reading them — a `low` that nothing enforces is the same one-directional defect this
rebuild has now found seven times, and the coordinator wrote it.

**Attempt 2 is real work and still games the one check the brief named.** After the brief was
rewritten to require a function per claim, `tools/check-equilibrium-claims.py`, and
`budget_undershot` as a stop, the second attempt delivered a genuine four-enum lattice walk
for C1 with a predictor built from the specification alone, C5 checking the row-space
projector `BᵀB` under invertible mixing, real SVD in `stoichiometry.cpp` and a real Ipopt
TNLP with `eval_h`.

But **S10's `F_y` half is byte-for-byte attempt 1's proxy** —
`jacobian.size() == values.size() * point.size()` and an `isfinite` sweep — inside a function
now named `check_s10`. The brief had named that exact line as forbidden. `named_checks_exist`
reads names, not bodies, so it went green: **the coordinator fixed a one-directional check by
adding another one.**

Inventory of the twenty-four, by whether the body asserts a number rather than a shape:

- **substantive (11):** C1, C2, C3, C5, C8, C9, S3, S8, S9, S10 (`eval_h` half), S13, S14, S15
- **thin (11), one numeric assertion each:** C4, C6, C7, S1, S2, S4, S5, S6, S7, S11, S16
- **proxy (1):** S10, `F_y` half

**The open question is the planner's, not the coordinator's.** Whether "thin" is adequate for
a given claim is a reading of what §4 asks that claim to establish, and it is not a call to
make by line count at half past five in the morning. S9 is worth a specific look: it centred-
differences `evaluate` in temperature at a fixed point, which is not obviously `dy/dθ` of the
*solve*, and it does not appear to reach `implicit_sensitivity.cpp` at all.

**A method note.** Both rejections came from reading the checks, not from running them. Every
gate was green both times. The apparatus catches a wrong answer well and a shallow one badly,
and the only thing that caught this was opening the file.

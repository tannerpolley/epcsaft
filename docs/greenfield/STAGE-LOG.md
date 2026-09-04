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
| lines, new files | 1,826 of 4,085 estimated | 3,539 as the worker reported its writes (by `git diff --numstat 88b5a2a`, 4,708 lines added: native 2,558, tests 778), **nothing below its `low`** |
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

## Records moved to the root, 2026-09-03

`engine/docs/scientific/` and `engine/docs/greenfield/` moved to `docs/scientific/` and
`docs/greenfield/` with `git mv`, by the owner's decision: the CSE record owner belongs at
the repository root, and the rebuild's record spans `data/`, `analyses/` and `tools/`, so it
was never an engine document. Engine-only documents (the equation registry, the bibliography
snapshot, the research notes) stay under `engine/docs/`.

Deleted in the same commit: `formulation.tex`, `methods.tex`, `evidence.tex` (every path
they pointed at had been deleted by the rebuild) and `adr/templates/` (generic starters,
never used). `docs/scientific/CONTEXT.md` is rewritten to the shape CSE issue 55 specifies:
question, owner map, glossary with sources, avoid table, default deliverable. Five terms
from `engine/CONTEXT.md`'s old glossary are listed without a source; retiring them belongs to
the stage that rewrites that file.

Every unpinned path was rewritten. Pinned refs of the form `greenfield@<sha>:engine/docs/...`
in closed briefs are left as written: they resolve at their commit and are provenance.
Repaired on the way: two dangling `plan:` entries in `engine/ARCHITECTURE.yaml`, the
doubled triage sentences in `engine/CONTEXT.md`, and the root `CONTEXT.md` still naming the
deleted Hard-Cutover Plan as the current program. The CSE test-gate hooks were not
installed: they block a session from stopping until its test fingerprint is approved, and
every stage worker runs `codex exec` unattended.

## Rule 13 loses its allowlist, and the wheel fingerprint was blind to `native/`, 2026-09-03

**Rule 13.** The owner ruled for `SURFACE.md` §4.2: every function bound in
`native/bindings/` opens a `core::Scope`, with no allowlist. `NON_SCOPED_BOUND_FUNCTIONS`
is an empty set; `build_parameters`, `counters` and `reset_counters` open a span. Shown red
with the allowlist emptied and the spans absent (three findings), green with them present,
and the spans are visible from Python: a captured `Parameters.from_json` has one root named
`build_parameters`. The coordinator wrote those three lines and one include on the owner's
direct instruction; it is the smallest engine edit made outside a stage.

**The fingerprint.** Proving the spans from Python failed at first: the installed wheel had
no `build_parameters` root. `tools/ensure-engine-wheel.sh` fingerprinted `engine/cpp` and
`engine/src`, and `engine/cpp` was deleted at Stage 1, so an edit under `engine/native/`
alone never changed the fingerprint and the cached wheel was reinstalled. Every stage so far
also edited `meson.build` or `src/`, which are in the fingerprint, so stage acceptance rebuilt;
a native-only edit did not. Now `engine/native` is hashed. `check-affected.sh` and
`check-analysis-layout.py` carried the same dead path and follow. This is the one-directional
shape again, in the build: a cache key that names a directory nothing checks exists.

## Stage 5 attempt 3 stopped, and the control plane reset, 2026-09-03

**Attempt 3.** The worker stopped itself `gate_impossible` after 47 minutes rather than
weaken a check. Reproduced by the coordinator in the worktree: the local-minimum check
(`equilibrium_checks.cpp:73`, reduced Hessian positive definite) fails because the
two-phase flash fixture collapses to one phase, and the dew–bubble symmetry check
(`boundary_checks.cpp:45`) fails on the recovered pressure. Every compile-tier check,
the architecture checker (15 rules over 30 sources), the 24 named claims and the wheel
build are green; `check-affected` is red by design until Stage 6. About 4,300 new lines
plus 221 edited. Not tagged, merged or resumed; the owner sent the two failures to the
Equilibrium planner for a ruling on design versus implementation before any resume.

**The reset.** An independent read-only audit the same day found the engine clean against
every checker rule and placed the drift in the transfer documents. The owner agreed and
had it applied before Stage 5 continues, as one canon change: `TRANSFER.md` no longer says
Stage 2 has not started, no longer argues for the two-stage split the crosswalk already
retired, and numbers the downstream stages 5, 6 and 7; its test figure is the prediction
beside the measurement (1,926 after Stage 4, about 2,400 expected). `CANON.md` says a line
count is never completion evidence in either direction. `CARRY-OVER.md`'s stage column
drops the retired numbers 7, 10 and 12. `target-layout.md` is marked historical and kept
because closed briefs pin it by sha. ADR 0005 is marked superseded in part, naming the
three decisions the rebuild reverses and what survives; ADR 0007's CAP-14 no longer names
a deleted test file. `engine/ARCHITECTURE.yaml` drops its stage status and the Pint
dependency, lists all seven layers the checker ranks, and points the constants home at
`core/`. The paired briefs and the preflight are untouched: they are the mechanism that
caught two of Stage 4's three plan defects, and they retire at cutover, not before.

## Attempt 3 tagged and rejected; the planner's revisions for attempt 4, 2026-09-03

Tagged `stage-05-attempt-3` at `d02f3c9` on `stage-05`. The planner reproduced both reds
from the attempt's own wheel and read them as one defect: the coincidence refusal guarded by
`pinned_count == 0`, so bubble, dew and cold dew on methane–ethane at 180 K all converged to
`x_L = y_V = 0.5/0.5` at 15.8 bar, 149 Pa and 15.8 bar. Behind it, four design gaps (no
coincidence scale or `k ≥ 1` wording, no start rule, no equilibrium row tolerance, fixtures
without states) and three procedural ones (abort-first check macro, a warm-seeded dew against
trap 13, no `diagnostics` array). The owner ruled: tag and reject; attempt 4 from the tag.

Applied to `EQUILIBRIUM.md` by the coordinator from the planner's draft, mapped by heading:
§2.1 bounds, §2.2 dimensionless rows, §3.4 the coincidence mechanism, §5 step 4 the three-rule
start and the tolerance-and-budget table, §4 preamble non-aborting checks, §4.2 the fixture
state table and S15 clause (c), §4.4 diagnostics, trap 6 and trap 13. One sentence outside
the draft was changed to agree with it: §3.4's "or the default: one phase at the feed, the
other displaced toward the pure component with the largest positive pair deviation" is now
the `start_symmetric` refusal, because rule 3 of the start rule removes that default. The
brief JSON gains the tokens, the moved budgets, `numerical_policy.hpp` as a write, five gates
and the `diagnostics` array; `preflight-brief.py` now refuses a native-writing brief without
one, which `CORE.md` §7 asked for and which would have refused attempts 1–3.

For the log, outside the stage: `Parameters.from_json` on the nine-species packet with
`components=("monoethanolamine", "water")` builds, then `HelmholtzModel` throws
"invalid ionic_dispersion". A neutral subset of an electrolyte packet cannot make a model.

## The layering draft: properties, thermochemistry, the regression shell, 2026-09-03

The planner's layering draft (revision 2, reviewed by an independent Fable 5.1 agent against
the tree and the three downstream applications) was approved by the owner and applied by the
coordinator. `EOS.md` gains a "Properties" section: `State` with a `PropertyOrder` flag,
totals absent by name without a declared `Thermochemistry`, the fixed-pressure derivative
block, `physics/ideal.hpp` on one basis, the standard-state transfer with a temperature
derivative and a cache, `K(T) → g°(T)` in the model, consistency reported not asserted.
`SURFACE.md` §2.5 makes `Thermochemistry` a second declared document with its own builder
and packet schema; rule 11 reaches `core/thermochemistry.hpp` once it exists. `EQUILIBRIUM.md`
§12 is the regression shell (four-node expression, no property names in `regression/`), the
levers gain the feed with S9's feed column, the seventh active family is named, and §13 is
the applications table with its provenance. `REGRESSION.md` is written; `TRANSFER.md`'s
Stage 6 drops the two observable files. Rule 4 gains `Discover`. The Stage 5 brief's S9 and
its sensitivity gate say `T`, `P` and the feed; preflight clean. `chemicals` 1.5.2 is a dev
dependency for the thermochemistry transcription tool.

**The registry is Markdown-only** (owner ruling the same day: "latex files are for
manuscripts"). Schema 3: `math` replaces `latex`, `label` is dropped, `source` is a plain
citation string (`\cite{Key}` became `[Key]`, `\eqref` became the identifier), and the
check no longer reads `references.bib`; a LaTeX-only command in `source` or `math` is
refused. Eleven entries added from the planner's bodies: the ideal-gas reference term, four
heat-capacity kinds, two `K(T)` kinds, total enthalpy and entropy, the transfer's temperature
derivative, and the reaction-derived thermochemistry. The reactive electrolyte two-liquid
oracle row is still pending the owner's ruling (candidate: the lithium extraction network).

## Attempt 4 launched, 2026-09-03

Owner's call: proceed, on `gpt-5.6-sol` at reasoning high rather than the luna model of
attempts 1–3 ("this needs to be done correctly"). The stage worktree was brought to the
attempt-3 tag merged with `greenfield@78f4bc7` (`d404a27`, clean merge; both sides had
edited the checker) so the worker reads the revised `EQUILIBRIUM.md` and brief from its
own tree. `uv sync`, checker green on 30 sources, preflight clean: 8 sources, 24 gates,
11 preconditions, 7 mutations deferred to acceptance. Launched with the verified recipe,
stdin closed, venv first on `PATH`; session `01a0690e-d72e-70b0-9d6e-5d299c9ef0d5`,
pid 1228759, log `attempt-4.log`.

**Interrupted and resumed, 14 minutes in.** The live digest showed "Identifying budget ceiling
breach" then "code compaction", "combining lines", "removing blank lines": `request_compile.cpp`
at 784 against a 780 ceiling. The planner ruled 600 / 750 / 900 (716 lines of real content at
the tag plus the feed column, the dimensionless rows and the `standard_gibbs` call) and that
the compaction is a rejection finding whatever the total; the owner chose to interrupt and
resume rather than let it run. Brief corrected at `7d6fe4d`, merged into the worktree, session
resumed with a prompt naming the ceiling, the rule, and a restore step. Two lessons for the
launch recipe: `$!` after `setsid nohup` is the wrapper's pid, not codex's, so the first kill
missed; and a killed run holds its thread-writer lock until the codex process itself is gone.
The interrupted run had already reached the boundary checks once, with S3 red on a mole
fraction; that is the resumed worker's to finish.

**First honest ceiling stop, and the second resume.** Twelve minutes after the resume the
worker stopped `budget_ceiling`: `equilibrium_checks.cpp` at 333 against 330 while writing S1's
four consequences per cell; no compaction, and it said which additions were untested. Every
gate it ran was green on its tree, S1–S16 (S12 excepted) included, and the compacted row
compiler was restored as first written (848 of 900). The planner re-budgeted the three test
files at 1.2× measured-plus-remaining (`1061a7d`); the owner chose to resume the same session.
The planner's note for acceptance: run `coincidence_every_cell` and `start_splits_composition`
before reading anything else, since those two would have caught attempt 3.

**Attempt 4 complete on the worker's tree; first acceptance run.** The second resume finished
`complete`: 24 gates green in brief order after its own clean rebuild. Committed as `6778be5`
on `stage-05`. Coordinator's acceptance (`accept.sh`, wipe and rebuild, every gate, every
mutation with the two attempt-3 catchers first): clean build green, 24 base gates green,
**16 of 24 mutations red**. The eight that did not go red: one acceptance-tooling artifact (the
abort gate's pre-fix command; red with the corrected gate); five stale seds whose anchors the
rewrite moved (`solve_boundary`, `compiled.F_y`, a global `pinned` rename that does not
compile, a first-occurrence `RowKind::Moiety` inside the residual classifier, a permutation
edit on one of two sites under the specification binary), re-anchored by the planner and
verified to change the file; and two real findings. **S15 clause (c) was not written**, so
with the coincidence floor zeroed nothing drives a boundary cell to the trivial solution: the
planner placed it as the binary bubble cell at 280 K, above the mixture critical temperature,
where only the trivial solution exists, in `boundary_checks.cpp`. **The composition-splitting
mutant did not reach the symmetric subspace**: with `K = 1` the compositions start at the feed
but each phase's `η` comes from its own branch; the subspace is equal composition *and* equal
`η` (trap 6 and §5 step 4 amended; the mutation now also forces one branch). Provenance note:
`request_compile.cpp` was first written at 838 lines, compacted to 781, and restored to 790;
its density (2.4 % blank, 0.64 statements per line) matches the tag version, so no penalty,
and the 48-line gap is recorded as unexplained. Verdict: not accepted as it stood; one
correction resume for clause (c), then the eight mutations re-run.

## Stage 5 — Equilibrium: attempt 4 accepted and merged, 2026-09-03

**Third resume, clause (c).** S15(c) added to `boundary_checks.cpp` (13 lines, 210 of 270):
the binary bubble cell at 280 K reached the `phases_coincide` refusal, not a budget stop.
Tree `d3f3392`, tagged `stage-05-attempt-4`.

**Acceptance run 2** (`accept.sh`, wipe and rebuild, every gate, every mutation with the
attempt-3 catchers first): clean build green, **24 gates green, 24 mutations red**, every
mutant built, worktree clean after restore. Merged into `greenfield` at `b40d5a7` as one
commit crossing engine and canon, as a stage does. Whole stage against `dd81eb5`: native
+3,236, bindings and Python +721, tests +1,067; test source now 2,993 lines, recorded.

**Read against `EQUILIBRIUM.md` §4.2.** The eight checks attempt 2 faked are real: feed-basis
invariance compares every row's kind, phase, index and identifier across four bases;
permutation builds the permuted point through the layout and relabels through a second
fixture; the residual verifier runs six solved cells and asserts Euler's relation per phase,
neutrality per ion phase, affinity per reaction per phase, the objective equal to the feed
dual, and the record's residual equal to the recomputed one; the local minimum on both
flashes; certificate, lever rule, isomerisation and reactive tie line as their rows say;
three levers against centred differences; both Jacobian halves; the dew cold-seeded.

**Recorded deviations, carried into Stage 6's brief, not gates:** no "skipped by name"
reporting when a certificate fails (the three-phase clause is absent rather than reported);
the reaction-direction zero on the feed sensitivity post-dates the brief the worker ran
against. Provenance: `request_compile.cpp` first written at 838, compacted to 781, restored
to 790, density equal to the tag's; 48 lines unexplained.

**Serena.** Owner's instruction: the coordinator, the planner and the workers read the tree
through Serena. It was pinned to the superseded checkout. Now: `.serena/project.yml` in the
greenfield and stage worktrees (gitignored), the Claude Code server's `--project` at the
greenfield worktree, a `serena` entry for the codex workers at the stage worktree. Takes
effect at the next session and the next launch; the next brief tells the worker to use it.

**Lessons that cost something today.** `$!` after `setsid nohup` is the wrapper's pid; find
codex with `ps` and kill that. A killed run holds its thread-writer lock until the codex
process is gone. `git merge` runs in whichever worktree you are in. A monitor's completion
line can be missed: check the log by hand when the expected duration passes.

## Stage 6 — Regression: brief, and attempt 1 launched, 2026-09-03

The planner's brief: twenty writes (estimate 2,076, ceiling 3,043, no lows per the owner's
ruling), twenty-two gates each with a mutation, nine sources pinned at `6bdc05f` with every
marker verified, two diagnostics entries, the five stops, a Serena reading step reported as
`serena_reading`. `REGRESSION.md` §4.1 carries the eight named checks R1–R8 and the two
Stage 5 deviations with their own mutations. Scope, owner's ruling: this stage writes the
regression shell and the `PropertyOrder` first-derivative block on `State`; thermochemistry,
the ideal side, the totals and the seventh active family form a Properties stage after it.
Worktree cut from `fb47283`, subprojects copied, fresh venv, Serena project file, the codex
Serena entry re-pointed. Preflight clean: 9 sources, 22 gates, 12 preconditions. Launched on
`gpt-5.6-sol` at high, session `01a069b3-a577-7fd1-94e8-f037660b4535`, pid 91635.

## Stage 5 repair: the second audit's findings, and Stage 6's two honest stops, 2026-09-03

**Stage 6 attempt 1** stopped twice before writing a line, both correctly. First
`source_contradicted`: the owner's scope split (Properties stage after Stage 6) had reached
`HANDOFF.md` and `TRANSFER.md` but not `REGRESSION.md` §3, `EQUILIBRIUM.md` §12 or `EOS.md`
Properties, which still required the seventh family and `standard_gibbs`; the coordinator wrote the
staging ruling into all three (`be3f329`), re-pinned, resumed. Then `gate_impossible`: a Problem
target's Jacobian for the six EOS families needs `dy/dθ = −F_y⁻¹ F_θ`, and the assembly's
`lever_jacobian` has only `T`, `P` and feed columns although §12 promised the active-parameter
columns as Stage 5's; Stage 6 may not change the assembly. Held until the repair merges.

**The audit** (read-only, `b74ae3b → f6efa0b`) found seven defects, all verified by reading, two
reproduced numerically by the reviewer: the feed sensitivity omits the pinned amount's dependence
on the feed total (22 % at `β = 0.2`); the temperature column omits the standard term's derivative
and compile evaluates `g°` once at the temperature guess; the Michelsen deficit divides dimensionless
rows by `RT`; `result.convergence` never reaches `scope.finish`; `compile_problem` has no
`keep_alive`; the Python mapping drops unknown species and turns an unknown amount status into
free; `Reaction(equilibrium_constant=1.0)` and `solvent_mass_kg=1.0` are physical defaults the
rule-12 checker missed because it scans top-level files only. An eighth: rule 13's allowlist had
returned through the Stage 5 merge, which kept the worker's copy of the checker; removed at
`be3f329`. Canon edits applied: S9's pressure-row derivative sentence (`−v_p/(RT)`), the
pinned-amount term in §12's `F_θ` sentence, the deficit without `RT` in §3.4, and S9's new cells.
The repair brief adds the EOS-parameter `F_θ` columns as finding 8 at the coordinator's request.
Lesson: a merge that "was clean" can still carry a policy reversal when both sides edited the same
tool; diff the enforcers after every merge.

**Repair launched.** Finding 8 (the EOS-parameter `F_θ` columns from the tape's mixed block,
`lever_jacobian` widened to `2 + C + A`, S9's active-parameter column) spliced into the brief;
§2.1 and §12 name the active parameters as a lever. Preflight clean at `be6ddf0`: 10 sources, 17
gates, 11 preconditions, 8 mutations deferred to acceptance. Launched on `gpt-5.6-sol` at high,
session `01a069d6-7868-7ec3-a66e-308fdd5ce3d6`, pid 252018, while Stage 6 waits.

**Repair, two honest stops.** First `source_contradicted` before writing: rule 12 is `SURFACE.md`'s
and governs every default under the subpackage, so widening the checker flags twelve defaults in
`problem.py`, not the two the brief named; and `sensitivity()` lives in `implicit_sensitivity.cpp`,
unlisted. Brief amended (`1fc76b6`): physical defaults required, structural empties `None` with an
explicit replacement, the file listed; then the planner's caveat (`363ff5d`): a replacement is an
empty structure, never a number, so `Free(guess=None)` and `packing_guess=None` cross as absent and
the start rule seeds them; the worker was interrupted before its first write to carry it. Second
`budget_ceiling`: `boundary_checks.cpp` 327 against 300 with three new centred-difference cells
(117 lines against 96); the planner's 96 was the estimate of two cells before finding 8 added the
third; re-budgeted 130 / 340 (`8cd5b58`), resumed with the file as first written and no refactor
of accepted cells. Note for the Stage 6 closeout: one helper for every sensitivity check, S9 and
R1 alike, if the reading finds verbatim copies of the centred-difference loop.

**Third stop, the coordinator's.** The re-budget put `ceiling_total: 340` under `lines` while the
brief keeps it at the write's top level, where 300 remained; the worker found two answers and
stopped `source_contradicted`, correctly. Fixed (`769efc5`), and the preflight now refuses a write
with `ceiling_total` under `lines`. Resumed a fourth time.

**Fourth stop, the fixture's.** Build and checker green on the repaired tree (30 sources), then
the boundary binary red at the new S9 reactive free-temperature cell: `{L,L}` with `L₂` pinned at
one percent, `T` free from 298.15 K, offsets `{+1000, −1000, 0}` J/mol; 200 iterations exhausted at
residual 7.3e-4. The worker's own diagnostic shows the assembled Jacobian's `T` column matching
centred differences at the start to 2.2e-11, so finding 2's assembly is right and the cell does
not converge. Stopped `gate_impossible` without touching fixture, start rule, driver or budget.
The planner is asked for the cell: the falsifier (the `T` column against centred differences of
the row evaluation) holds at any point, converged or not, so it need not ride on a free-`T` solve;
whether that free-`T` `{L,L}` solve should converge is a separate question, possibly its own
finding. Environment note: the repair worktree had subprojects but no configured build directory,
and the worker configured one itself and reported it as a deviation; the coordinator now runs
`meson setup` in every new worktree before launch (done for stage-06 as well).

**Finding for its own brief, not the repair's (planner's ruling).** The reactive `{L,L}` cell with
`L₂` pinned at one percent, `T` free at 1 bar, offsets `±1,000 J/mol`, exhausts
`EQUILIBRIUM_MAX_ITERATIONS` at residual 7.3e-4 with the assembly verified at the start. Either the
cell has no solution near 298 K with those offsets (a fixture-state question; certify it as the
state table requires) or the free-`T` boundary path needs its own continuation start. The trace
of that run (`gf-stage-05-repair/attempt-1.resume4.log`) is the evidence to keep. The repair's
S9 standard-term check becomes a derivative comparison at the start point, no solve.

**Fifth stop, three lines.** The two-part point comparison that replaced the free-`T` solve made
`boundary_checks.cpp` 343 against 340 (133 added against 130). Coordinator's re-budget, 145 per
write and 355 total, on the ground that the planner's own replacement is the whole overrun;
recorded here so the planner can object. No compaction, as before.

**Sixth stop, the coordinator's wording.** The build gate failed on an undeclared identifier in
the worker's own S9-T-standard replacement (`reactive_compiled` where it defined
`reactive_free`), and the resume prompt had said to keep the file "exactly as written", so the
worker stopped `source_contradicted` rather than fix its bug. The rule was always "no compaction
to fit a count", never "no correction"; the next prompt says so in those words. Held there on
purpose: a third audit (through `ecdccd7`) found that finding 8's column rule omits the chain rule
through the packing coordinate, since `v = c(T, x, θ)/η` with `c` built from segment count and
diameter (verified in `packing_terms()`), so the next resume carries the corrected derivation, a
stronger `reject_by_name` mutation, the `FreeValue` absent state the worker had not made, and the
bug-fix clarification, in one amendment. Lesson, the audit's and fair: reconcile the mathematics,
the checks, the permitted files and the budget before a launch; six of seven stops on this repair
were brief or prompt defects, not code.

**Seventh stop, a permitted-files gap.** The corrected finding 8 differentiates the packing
coefficient per active column and needs each active parameter's family and indices, which
`HelmholtzModel` keeps private; no model-layer write was listed. `helmholtz_model.hpp` is now a
write for one read-only accessor (`29d2fc7`). The worker had already fixed its identifier bug and
made the native absent state within their ceilings. Eighth resume.

**Eighth stop, a name collision the checker caught.** Rule 10 refuses any unqualified `active(`
outside `core/` and `bindings/`, since that token is the trace's accessor; the coordinator had
named the model's accessor `active()`. Renamed `active_parameters()` (`1509e3e`). Ninth resume.

**Ninth stop, an ill-posed cell.** S9-T-standard (ii) as written fixed `T` with `L₂` pinned and
`P` fixed, which the count rule refuses (`k = 1` needs `f = 1`). The coordinator read the intent,
"the same problem with `T` fixed", with `P` free, and resumed; the planner is told. On the same
run S4, S5, S10 and S14 were green on the repaired tree.

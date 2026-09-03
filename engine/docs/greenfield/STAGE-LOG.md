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
- **`_core` has no sources until Stage 4**, and Meson warns that it allows this
  only through a bug whose support will be removed. Accepted for two stages
  rather than writing a stub module to silence it.
- **Five rules, not four.** The regression boundary is its own rule because it is
  a symbol rule, not an include rule.

## Stage 2 — physics

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
  templates on an AD type until Stage 3, so this is checked here.

## Stage 4 — bindings and Python EOS

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

**Corrected 2026-09-02, and superseded by `DIELECTRIC-SUPPRESSION.md`.** The account below
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

Not a gate failure and not a merge blocker: `TESTING.md` is explicit that a number out of
the old engine is not evidence the new one is wrong, and the criterion the planner chose
came from the parameter source rather than from this comparison. But an average over six
rows cannot see a missing upturn in the last one, and a check that a wrong physics passes is
the thing `TESTING.md` exists to prevent. Raised with the planner: whether the criterion
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
  `[Fact]` markers. That is live science for Stage 11, not a specification.
- `equilibrium/research/…-carbonate-provider-parameter-sources.md` names the retired
  Provider, but its decision is *which species and which sources* — that no
  source-complete ePC-SAFT parameters exist for two calcium complexes is a fact about the
  literature. Its "Required Provider records" section is stale and is flagged here rather
  than edited.

`equilibrium/reactive-bubble-observations.md` was the sharpest deletion: it specifies a
bubble owner solving `g(q,P,θ) = 0` as an outer Newton over pressure around a separately
certified liquid — precisely the nested coupling `TOLERANCES.md` §6 names as what Stage 9
must not reproduce. A future agent finding it in the tree would read it as the design.

Deletion in a versioned tree destroys nothing, which is why the caution about "the only
record" applies to ADRs and not to a stale plan. Every dangling reference left by the
deletions was repaired in the same pass — seven files, including two ADRs and
`target-layout.md`.


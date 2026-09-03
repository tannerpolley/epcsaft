# Stage 5 — the core

You are implementing Stage 5 of a greenfield rebuild in this worktree (branch `stage-05`,
cut from `greenfield@21b8695`). You have no one to ask. Read, in order, before writing
anything:

1. `engine/docs/greenfield/TESTING.md` — a hard restriction on what may be tested
2. `engine/docs/greenfield/CORE.md` — the specification of this stage, completely
3. `engine/docs/greenfield/CATEGORIES.md` — why four files leave `physics/` and one leaves `model/`
4. `engine/docs/greenfield/TOLERANCES.md` §2 and §5 — kept by reference
5. `engine/docs/greenfield/DENSITY-ROOT.md` §1, §6, §7 — why the density record carries no floor yet, and the accuracy sentence
6. `engine/native/physics/*.hpp`; `engine/native/model/*.hpp,*.cpp`; `engine/native/bindings/*.cpp`;
   `engine/tests/native/*.cpp,*.hpp`; `engine/src/epcsaft/*.py,*.pyi`
7. `tools/check-engine-architecture.py` — you extend it; every rule must stay green
8. `engine/docs/greenfield/briefs/stage-05.brief.json` — sources, budgets, gates, mutations

Verify every source marker and every precondition in the JSON against the tree **before**
writing. If a claim is false, stop and say which.

## 0. Non-negotiable

**Layers.** `core <- physics <- model <- equilibrium <- regression <- surface <- bindings`.
This stage creates the bottom layer and the surface layer. A `core/` file includes only
`core/` and the standard library; `scalar.hpp` is the one core file that may include CppAD.
Includes are rooted at `engine/native`, never relative. Only `native/bindings/` includes
pybind11.

**Nothing in the engine reads the trace.** `core::Capture` and `core::active()` are named
only under `native/core/` and `native/bindings/`. Engine layers use `core::Scope`, which is
write-only. A floor crosses a layer as a return value.

**No answer changes and no value moves.** This stage adds records, spans and declarations;
it moves no tolerance value and alters no iteration. A captured call returns the same bits as
an uncaptured one. `DENSITY_ROOT_RESIDUAL` becomes a `derived()` edge with the same value.

**No hidden default.** A selector is declared or the call is refused by name — natively for a
`Parameters` struct (the model constructor), natively for a packet (the builder). No enum
member, struct default, builder branch, Python `.get(..., default)` or signature default may
supply one. Rule 11 forbids `Unspecified` in the parameter record; rule 12 forbids any
signature default but `None` on the surface.

**Constants and tolerances.** No numeric literal acting as a threshold in `engine/native/`
outside `native/core/numerical_policy.hpp`. Two entries are added (`FLOOR_MARGIN`,
`STALL_WINDOW`).

**Python does no arithmetic.** The three ratios on `Convergence` are native properties.

**Naming.** No proper nouns in identifiers. Event names are registered `EqID`s, policy entry
names, or field names of the result the function returns — never a route, solver, author or
stage. Rule 7 scans string literals in solver files.

**Testing.** Four properties and nothing else (`CORE.md` §8): the classifier's partition and
the window's closed form; the two guards (`derived()` inversion refuses to compile,
`Declared::require` refuses by name); the record's residual recomputed from the returned
answer; a captured call bit-identical to an uncaptured one with a trace that is a tree. No
check names a span, an event, a count or an iteration. Everything else is reported.

**Environment is not yours to repair.** You may edit `engine/meson.build` only to move the
builder's source line into a `surface_sources` list, to add one test executable and one
`test()` line. Not the dependency block, `default_options`, the virtual environment,
`uv.lock`, `engine/build/subprojects`, or how you were invoked. `PATH` leads with the
worktree's `.venv/bin`; `$HOME/.cache` is writable. Stop with `environment_broken`.

**Physics is open for exactly these edits**: include paths; the selector reads rewritten as
`.is(...)` (§2.2); the `solve_site_fractions` out-parameter and the `solve_linear` pivot ratio
(§2.4). No kernel arithmetic, tolerance value or other signature changes.

## 1. Stop conditions — stopping is a first-class outcome

- `source_contradicted` — a source lacks its stated marker, or says something other than this
  brief claims. Name source and marker.
- `gate_impossible` — a gate cannot run as written, or cannot be made to fail by its named
  mutation. Name the gate and why.
- `budget_ceiling` — a file would exceed its ceiling. Stop *before* exceeding it.
- `environment_broken` — configure, compile, install, or a precondition fails for a reason
  outside the tree. Quote the command and its first error line. Do not work around it.

An implementer that quietly absorbs a wrong plan is the failure. A correct stop is not.

## 2. Design

`CORE.md` is the specification; this section fixes the lines the gates mutate.

### 2.1 The moves, and `numerical_policy.hpp` made `constexpr` and self-naming

`git mv` four files from `native/physics/` to `native/core/`: `physical_constants.hpp`,
`numerical_policy.hpp`, `parameters.hpp`, `scalar.hpp`; and `native/model/parameter_builder.{hpp,cpp}`
to `native/surface/`. Then every `#include "physics/<one of the four>"` becomes
`"core/<file>"` (19 lines under `native/`, *measured*), `"model/parameter_builder.hpp"` becomes
`"surface/parameter_builder.hpp"` in `bindings/eos_bindings.cpp`, and `meson.build`'s
`model_sources` loses the builder to a new `surface_sources` list passed to the archive beside
the others. `grep -rn "physics/parameters\|physics/physical_constants\|physics/scalar\|physics/numerical_policy\|model/parameter_builder" engine/native engine/tests tools engine/meson.build`
must then find nothing.

In `core/numerical_policy.hpp`:

- `struct Tolerance final { const char* name; double value; Category category; const char* source;
  const Tolerance* upstream = nullptr; double factor = 1.0; };`
- every entry: `inline const Tolerance NAME{` → `inline constexpr Tolerance NAME{"NAME",` —
  27 entries, values and sources unchanged, the value on its own line as today.
- `derived()` exactly as `CORE.md` §3.2, with the guard line written exactly as

  ```cpp
    if (!(factor >= 1.0)) throw "derived tolerance tighter than its input";
  ```

- `DENSITY_ROOT_RESIDUAL` becomes `derived(ASSOCIATION_RESIDUAL, 5.0, "DENSITY_ROOT_RESIDUAL", "<why>")`;
  `2e-13 × 5 = 1e-12`, unchanged. Its `why` carries `CORE.md` §3.2's hypothesis sentence and
  `DENSITY-ROOT.md` §7's accuracy sentence: *the tolerance is on the reduced residual; the
  delivered relative pressure error is |r|/z, about 1e-8 on a dense liquid root, measured worst
  1.7e-8*.
- add `FLOOR_MARGIN{"FLOOR_MARGIN", 10.0, ...}` and `STALL_WINDOW{"STALL_WINDOW", 5.0, ...}` under
  `// --- convergence`, rationale from `TOLERANCES.md` §5, marked guessed.

### 2.2 `core/declared.hpp` and the parameter record

`MissingDeclaration` and `Declared<T>` exactly as `CORE.md` §4.1. In `core/parameters.hpp`,
`ElectrolyteForm electrolyte_form = ElectrolyteForm::None;` becomes
`core::Declared<ElectrolyteForm> electrolyte_form;` and likewise `permittivity_form`; the enums
keep their members (`None` is a declared choice). Every read site becomes `.is(...)`:
`electrolyte.hpp` (8), `born.hpp` (1), `helmholtz_model.cpp` (3), `parameter_builder.cpp` (2 reads;
its assignments are unchanged in form because `Declared` constructs from the enum). The
`TapeParameters` mirror in `helmholtz_tape.hpp` copies the two `Declared` fields.

`helmholtz_model.cpp` constructor, after the `dielectric solvent` requirement, written exactly as

```cpp
    require(has_ion == parameters_.electrolyte_form.declared(), "electrolyte_form");
    require(parameters_.electrolyte_form.declared() == parameters_.permittivity_form.declared(), "permittivity_form");
```

with `has_ion` computed as in the repair brief's §2.1 (gate `no_hidden_selector_default`
mutates the first line). `tests/native/model_support.hpp`'s two electrolyte fixtures already
declare both; the neutral fixtures and `physics_checks.cpp`'s `neutral_binary()` leave them
undeclared.

### 2.3 `core/convergence.hpp` and `core/trace.hpp`

Exactly as `CORE.md` §5.1–5.2 and §6.1. Header-only: `inline thread_local Trace* active_trace
= nullptr;`. `classify`'s first line is written exactly as

```cpp
    if (residual <= criterion.value) return Stop::converged;
```

`describe` uses `std::snprintf` into a `std::string`: `"<stop>: residual <r> against <name>
<value> after <n> iterations, contraction <c>"`, `floor`, `amplification` and `contraction`
printed as `absent` when absent. `ConvergenceFailure` derives from `std::runtime_error`, holds
`Convergence record`, `what()` is `describe(record)`. `Scope`'s constructor performs the one
thread-local read; every `probe` begins `if (trace_ == nullptr) return;`;
`steady_clock::now()` is read only when `trace_` is non-null. `Span::parent` and
`Trace::current` are `std::optional<std::uint32_t>`. A `Capture` constructed while another is
active shadows it and restores it.

### 2.4 Physics: the association Newton returns a record

`association.hpp`: `solve_linear` takes a trailing `double* pivot_ratio = nullptr` and sets
`*pivot_ratio = largest / smallest` before the solve. `solve_site_fractions`'s trailing
`std::size_t* iterations = nullptr` becomes `core::Convergence* record = nullptr`. Open
`core::Scope scope("association_newton");` after the `sites == 0` early return (which fills a
record with `residual 0`, `iterations 0`, `floor 0`, `amplification 1`, `stop converged` when
given). Keep a `ResidualHistory`; push `residual_norm` each iteration;
`scope.probe("association_residual_norm", residual_norm, iteration)`. On convergence fill the
record (`residual = residual_norm`, `criterion = &policy::ASSOCIATION_RESIDUAL`,
`floor = std::numeric_limits<double>::epsilon() * sites` — a roundoff statement, not a
threshold — `amplification` = the last Newton step's pivot ratio, absent when that step was
substitution, `contraction`, `iterations`, `evaluations` = residual evaluations made,
`stop = classify(...)`), `scope.probe("site_fractions", fractions)`, `scope.finish(*record)`
when given, and return. The domain exit throws `core::ConvergenceFailure` with
`stop = left_domain`; the budget exit throws it with `stop = classify(..., budget_hit = true)`.
`density.cpp`'s `association_failure` predicate tests the exception type and the stop instead
of the two message strings. `model_floor.cpp:169-174` passes a `Convergence` and prints
`describe(record)`.

### 2.5 Model: `evaluate`, the density root, the reference

`helmholtz_model.hpp`: `HelmholtzDerivatives` gains `double gradient_floor = 0.0;` with the
per-order comment of `CORE.md` §5.3; `DensityRoot` gains `core::Convergence convergence;`;
`State` gains `std::optional<core::Convergence> density_convergence;`.

`helmholtz_model.cpp` `evaluate`: `core::Scope scope("evaluate");` after the lock;
`scope.probe("declared_choices", values)` where `values` holds one `{field index, enum value}`
pair per *declared* selector in `Parameters` field order, so absence is absence;
`scope.probe("recording", 1.0)` beside `++counters_.recordings`; after `Forward(0)`
`scope.probe("residual_helmholtz_extensive", result.value)`; the `Gradient` and `PressureSlope`
branches probe `gradient` and, for the latter, `volume_second_derivative`; the `Hessian` branch
probes `gradient`, `hessian`, `volume_second_derivative`, `site_fractions`, and computes
`gradient_floor` as the ∞-norm over rows of `(factorization.solve(cross)).transpose() * D`,
`D = diag(total * mass_a / sites[a])`, times `policy::ASSOCIATION_RESIDUAL.value`, zero when
`sites_count == 0`; `scope.probe("gradient_floor", result.gradient_floor)`. The `PressureSlope`
branch leaves `gradient_floor` at zero — 04.repair fills it. Do not touch the dead double path
(`electro`, `born`, `direct`) or the `PressureSlope` assignment line: both are 04.repair's.
`state_at_pressure` copies `root.convergence` into `state.density_convergence`.

`density.cpp`: `core::Scope scope("solve_density");` at the top of `solve_density`; the
anchored attempt runs inside a child `core::Scope anchor("density_anchor");` in its own block
with `anchor.probe("density_newton_step", {density, r, slope, candidate}, iteration)` per step,
finishing its own record (`left_domain` when the candidate leaves `(lower, upper)`, `stalled`
when the derivative is not finite or zero, `converged` on acceptance); the scan probes
`density_scan_point` and `density_bracket`; `refine` takes the parent scope by reference and
probes `density_newton_step` per step and `density_root` per refined root; it keeps a
`ResidualHistory` and fills `DensityRoot::convergence` with `criterion =
&policy::DENSITY_ROOT_RESIDUAL`, `floor` and `amplification` **absent** (comment: *frozen-site
slope until 04.repair; DENSITY-ROOT.md §1*), `contraction`, `iterations`, `evaluations`,
`stop = classify(..., budget_hit = iterations == maximum)`. The line that stores the residual
is written exactly as

```cpp
    record.residual = std::abs(point.residual);
```

The budget exit throws `core::ConvergenceFailure{record}` in place of the bare
`runtime_error`. `solve_density` has **one exit**: every path — anchored acceptance, single
stable root, branch selection — assigns `DensityRoot result` and falls through to

```cpp
    scope.probe("density_root", result.molar_density);
    scope.finish(result.convergence);
    return result;
```

written exactly so. `evaluations` counts `at_density` calls for the whole solve.

`infinite_dilution_reference.cpp`: `core::Scope scope("infinite_dilution_reference");`; per
ladder point `scope.probe("reference_ladder_point", {molality, ln φ₁, ln φ₂}, index)`;
`scope.probe("reference_convergence_error", convergence_error)`; the witness likewise. No
record: a bias is not a stop.

### 2.6 Bindings and Python

`bindings/diagnostics_bindings.cpp` (new): `register_diagnostics(py::module_&)` binding `Stop`,
`Category`, `Tolerance` (readonly `name`, `value`, `category`, `source`, `factor`, `upstream` as
an optional reference), `Convergence` (readonly fields; `criterion` with
`py::return_value_policy::reference`; properties `residual_over_criterion`,
`residual_over_floor` (None when absent), `amplification_over_declared` (None when absent);
`__repr__` = `describe`), `Event`, `Span`, `Trace`, and `Capture` as a wrapper holding
`std::optional<core::Capture>` with `__enter__` (emplace, return self) and `__exit__` (take the
trace, reset), exposing `trace`. `module.cpp`: `ConvergenceError` registered after `EosError`
as its subclass, `DeclarationError` after `ParameterError` as its subclass; then
`register_diagnostics(module)`. `eos_bindings.cpp`: one `def_readonly("density_convergence", ...)`.

`src/epcsaft/eos.py`: `Mixture.reference(..., stoichiometry: tuple[int, int])` loses its
`= (1, 1)` default (taken from the repair brief); `tests/eos/test_activity.py` passes `(1, 1)`.
`src/epcsaft/__init__.py`: export `trace = _core.Capture`, `Stop`, `Convergence`,
`ConvergenceError`, `DeclarationError`. `_core.pyi`: the new classes and the changed signature.

### 2.7 The checker

`tools/check-engine-architecture.py`: `LAYERS = ["core", "physics", "model", "equilibrium",
"regression", "surface", "bindings"]`; `CONSTANTS_HOME` and `THRESHOLD_EXEMPT` move to
`core/`; rules 8–12 exactly as `CORE.md` §8.3, each a function in the existing style added to
the table in `main()`. Rule 9's predicate: `policy::\w+_(MAX_ITERATIONS|BUDGET)\b` over
comment-stripped text. Rule 10's forbidden include: `<(cppad|Eigen|Ipopt|coin|ceres|Highs|nlopt|boost|pybind11)`
in any `core/` file but `scalar.hpp`, and any `"..."` include from `core/` not starting
`core/`; its read-back pattern over files outside `core/` and `bindings/`:
`(?<![A-Za-z0-9_])(core::)?(Capture\b|active\s*\()`. Rule 11: collect `enum class (\w+)` names
in `core/parameters.hpp`; red on `Unspecified` anywhere in it, and on any line inside
`struct Parameters` matching `^\s+NAME\s+\w+\s*[=;{]` for a collected NAME. Rule 12: `ast.parse`
each `src/epcsaft/*.py`; every default in every `FunctionDef`/`AsyncFunctionDef` must be
`ast.Constant` with value `None`.

### 2.8 Build

`meson.build`: `surface_sources = files('native/surface/parameter_builder.cpp')` passed to
`static_library` beside the others; one executable `diagnostics_checks` in the style of
`model_checks`; `test('diagnostics', diagnostics_checks)`. `native/core/` and `native/surface/`'s
header add no source. Then `meson compile`, `meson test physics model diagnostics`,
`meson test --suite floor`, `python3 tools/check-engine-architecture.py`, `meson install` and
the Python gates — all clean, worktree free of stray files.

## 3. Checks — `tests/native/diagnostics_checks.cpp`, meson test `diagnostics`

Use the `CHECK` macro from `physics_checks.cpp`. Fixtures from `model_support.hpp`.

1. **Partition.** Construct `policy::Tolerance` values at runtime (an aggregate) and walk
   `residual ∈ {0, v/2, v, 2v, 5v, 20v, 200v}`, `floor ∈ {absent, v/10, v, 10v}`,
   `contraction ∈ {absent, 0.5, 1.0, 1.5}`, `budget_hit ∈ {false, true}`. Assert the four
   rules of `CORE.md` §8.1 and monotonicity in `residual`.
2. **Window closed form.** Push `r_0 q^k`, `q = 0.7`, 12 steps: absent for `k < STALL_WINDOW`,
   `0.7` to `1e-12` relative afterwards.
3. **Guards.** `derived()` at runtime throws for `0.5` and returns `3 × upstream` with
   `upstream` pointing at its input for `3.0`. `Declared<ElectrolyteForm>{}.require("x")` throws
   `MissingDeclaration` whose `what()` contains `"x"`; `.is(...)` is `false`;
   `aqueous_electrolyte()` with `electrolyte_form` reset to undeclared makes the model
   constructor throw `std::invalid_argument` naming `electrolyte_form`; `neutral_binary()` with
   `ElectrolyteForm::None` declared throws likewise.
4. **Record matches answer** on every fixture at `model_checks.cpp`'s `density_and_reference`
   conditions plus ethane at 200 K and 240 K: `state_at_pressure` returns `density_convergence`;
   its `residual` equals `|Z − P/(ρRT)|` from `state(T, ρ*, x)` to `1e-15` absolute; `stop ==
   converged` implies `residual <= criterion.value`; `criterion->name` is
   `"DENSITY_ROOT_RESIDUAL"`. For the association solver on the associating binary and nine
   species: `solve_site_fractions(..., &record)` then `association_residuals(...)` at the
   returned fractions has ∞-norm equal to `record.residual` to `1e-15`.
5. **Trace is a tree and changes nothing.** Per fixture: `state_at_pressure` uncaptured, then
   the same call inside a `Capture`; every `double` field `==`, every vector elementwise `==`,
   optionals agree in presence; over `trace.spans` exactly one has no `parent` and every
   other's `parent` is less than its own index. Nothing else about the trace is examined.

`tests/native/derived_inversion.cpp` (8 lines): includes `core/numerical_policy.hpp` and
declares `inline constexpr policy::Tolerance INVERTED = policy::derived(policy::ASSOCIATION_RESIDUAL,
0.5, "INVERTED", "must not compile");`. Not a Meson target; the gate compiles it with
`c++ -std=c++17 -fsyntax-only` and requires failure.

`tests/eos/test_diagnostics.py`: `with trace() as t: state = mixture.state(T, P=..., x=...,
phase="liquid")`; `state.density_convergence` is not `None`; the root span of `t` (the one with
`parent is None`) carries a `convergence` whose `residual`, `iterations` and `stop` equal the
state's; when `stop` is `Stop.converged`, `residual <= criterion.value` (two native fields
compared, which is what gate `python_trace_round_trip`'s mutation breaks);
`repr(state.density_convergence)` contains `"DENSITY_ROOT_RESIDUAL"`;
`Mixture.reference(...)` without `stoichiometry` raises `TypeError`. Python compares fields
and does no arithmetic.

## 4. Reported, not asserted — `model_floor.cpp`

For every density case: `describe(root.convergence)` on its own line. Case C: the association
record. New case G: nine-species `Liquid` at 313.15 K, 1e5 Pa cold, once uncaptured and once
inside a `Capture`: `wall_us` for both, `spans`, `events`, total `values` bytes. On this tree the
associating binary at 313.15 K, 1e5 Pa `Liquid` reports `budget_exhausted` (`DENSITY-ROOT.md`
§1); print it, do not assert it, do not "fix" it — that is 04.repair.

## 5. What to report

A JSON object per `ORCHESTRATION.md` §"The result schema": `files_written[]`, `gates[]` with
actual exit codes, `deviations[]`, `stop_detail`. Do not run the mutations. Do not claim a
number you did not measure.

Deferred deliberately, do not write: a trace renderer; a `policy` submodule; a runtime
`derived()` binding; the `IonicDispersion` selector, the `ELECTROLYTE_SERIES_THRESHOLD`
deletion, the dead double path, the density floor and amplification (all 04.repair's); the
`DensityBranch::Unspecified` rename; `profile_entry_points.py`'s deletion; the equilibrium,
regression and discovery rows of `CORE.md` §6.3.

## Appendix A. The delta 04.repair's brief takes when re-based on this stage

For the coordinator; not for this worker.

- `base`: the Stage 5 merge sha. Every path to the four moved headers and the builder changes
  (`physics/parameters.hpp` → `core/parameters.hpp`, `model/parameter_builder.*` →
  `surface/parameter_builder.*`); re-verify every marker on `density.cpp`,
  `helmholtz_model.{hpp,cpp}`, `association.hpp`, `numerical_policy.hpp` and `parameters.hpp`.
- §2.1: `IonicDispersion` has two members, no `Unspecified`; the field is
  `core::Declared<IonicDispersion> ionic_dispersion;`; the kernel line reads
  `parameters.ionic_dispersion.is(IonicDispersion::ExcludeSameSign) && charge_i * charge_j > 0`;
  the constructor line is `require(has_ion == parameters_.ionic_dispersion.declared(),
  "ionic_dispersion");` and gate `no_hidden_dispersion_default`'s sed targets it. Rule 11
  would otherwise go red on `Unspecified`.
- §2.2: the builder assigns `result.ionic_dispersion = IonicDispersion::ExcludeSameSign` etc.;
  the "missing" rejection stays as written. `Mixture.reference`'s default removal and its test
  edit are already done here; drop them from the repair's writes.
- §2.3 `PressureSlope`: with the factorisation lifted, also set `result.gradient_floor` from the
  `V` row (`|cross_v.dot(factorization.solve(D))| · ASSOCIATION_RESIDUAL`, `D` the vector
  `total * mass_a / sites[a]`) and probe it. `density.cpp`'s record then fills
  `floor = V · gradient_floor` and `amplification = 1 / |dr/dρ|` from the corrected slope; the
  *frozen-site slope until 04.repair* comment is deleted. `DENSITY_ROOT_RESIDUAL.factor` is then
  read against `floor / ASSOCIATION_RESIDUAL` in the floor harness: report it.
- §2.4: `ELECTROLYTE_SERIES_THRESHOLD` is deleted from `core/numerical_policy.hpp`.
- §2.5: the `DensityRoot` literal gains the record; the far-side walk probes
  `far_side_residual` per grid point into the `density_anchor` span, and a rejected anchor
  finishes that span's record with `stalled`.
- Ceilings: `density.cpp` +25, `helmholtz_model.cpp` +10.
- A `diagnostics` block (`CORE.md` §7) naming `solve_density` and `density_anchor`, with the
  check *`amplification` equals `1 / |dr/dρ|` recomputed from `state().pressure_density_derivative`
  at the root, and `floor` equals `V · gradient_floor` from `evaluate` at the root* — a property
  of the pair (answer, record) — and the mutation `s/record.amplification = /record.amplification = 0.0 * /`
  on `density.cpp`.

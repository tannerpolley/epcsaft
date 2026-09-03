# The core: numerical policy, the two contracts, and a debug mode that reads the tape

Status: design, 2026-09-02, against `greenfield@21b8695`. Successor to `TOLERANCES.md`, which
it absorbs (§9), and the decision document for `CATEGORIES.md`'s open questions (§10).
Specification of Stage 5 and of the `diagnostics` block every later brief carries (§7).

Evidence classes: **measured** (read from a file or a run on this tree today), **derived**
(follows from the code or the equations by an argument given here), **guessed** (no evidence
yet; each one marked). Nothing in §2–§8 is derived from the superseded tree; it appears only
in §1 as the shape of the failure.

## 1. What this replaces, and the evidence for the scope

`TOLERANCES.md` designed nested solver tolerances with a `Convergence` record and named
stops, scoped as one stage between the bindings and the minimizer. Three findings since then
make it a core with two contracts rather than a stage with one.

**The record cannot be built where it was planned.** `DENSITY-ROOT.md` §6 (*measured*): the
`amplification = 1/|dr/dρ|` the record wants is the quantity the frozen-site slope corrupts,
by 14.51× on the associating binary (`EOS-AUDIT.md` A.2), and the repair that fixes it
(`briefs/stage-04-repair.brief.md`) still ends its density loop with a bare
`throw std::runtime_error("density root refinement did not converge")` (*measured*,
`density.cpp:93`), the failure shape of `TOLERANCES.md` §1. The repair must report
convergence honestly and cannot, because the record is not defined. So the record precedes
the repair.

**Reading the engine's own numbers cost six probes and a rebuilt file.** Localising one
activity-coefficient discrepancy (`analyses/reproduction/2025-figiel-engine-crosswalk`,
*measured*) took six scripts, 203 lines, and a native probe linking the archive
(`EOS-AUDIT.md` A.2) to read the five Helmholtz terms at a fixed density, `Φ_VV` on two
derivative orders, the pressure each root actually satisfies, and the reference ladder's
convergence — every one a number in a register on the way to the answer. The requirement:
**in debug mode a caller sees every residual, iteration count, tolerance, floor and tape
intermediate along one method call, across EOS, equilibrium and regression, without
recompiling and without writing a probe.**

**The input side broke the same way.** The `ionic_dispersion` selector vanished in the
rebuild because nothing in the input surface carried it, at 49 % against the published model
and a week to find (`like-ion-dispersion-convention.md`, *measured*); `stoichiometry=(1, 1)`
is a hidden physical default (`EOS-AUDIT.md` C.5, *measured*). "No hidden defaults" and
"every solver reports its convergence" are one discipline: **a decision is declared or the
call is refused by name, and what happened is reported with numbers or the result is
refused.** Both promises had no owner. The core owns both.

So the core is

    core  =  numerical policy  +  the declared-choice contract  +  the convergence contract  +  debug mode

and `CATEGORIES.md`'s measurement says which existing files are already its tenants.

## 2. Where it lives

`engine/native/core/`, below `physics`. Seven headers, no source file; four move in, three
are written:

| file | from | owns | lines (low / est / ceiling) |
|---|---|---|---|
| `physical_constants.hpp` | `physics/` (28) | category-1 constants, defined once | 28 / 28 / 30 |
| `numerical_policy.hpp` | `physics/` (248) | every threshold, budget, step factor and discretisation number as a self-naming `constexpr` `Tolerance` with upstream and declared factor; `derived()`; `FLOOR_MARGIN`; `STALL_WINDOW` | 270 / 300 / 340 |
| `parameters.hpp` | `physics/` (86) | the parameter record; its selectors as `Declared<Enum>` | 95 / 105 / 120 |
| `scalar.hpp` | `physics/` (41) | the transcendental dispatch for `double` and the tape scalar | 41 / 41 / 45 |
| `declared.hpp` | written | `Declared<T>`, `MissingDeclaration` — the declared-choice contract | 30 / 40 / 60 |
| `convergence.hpp` | written | `Stop`, `Convergence`, `classify()`, `ResidualHistory`, `ConvergenceFailure`, `describe()` | 90 / 110 / 150 |
| `trace.hpp` | written | `Event`, `Span`, `Trace`, `Capture`, `Scope` — the debug mode | 140 / 170 / 220 |

And one directory for the input surface: `engine/native/surface/parameter_builder.{hpp,cpp}`,
moved from `model/` unchanged (30 + 395 lines, *measured*), between `regression` and
`bindings` in the include order.

**Why `core/` and why these tenants.** `CATEGORIES.md` measured it: in `physics/`, four files
carry all 45 `EqID` annotations and four carry none, and the four with none are exactly the
policy, the parameter record, the constants and the scalar dispatch. None answers "what is
true of the fluid"; all answer "how the engine computes and reports" or "what a caller
declared". `target-layout.md:86` rejected a directory for the policy header because it would
hold one file; it holds seven. A root-level file, the other alternative it named, sits outside
the include rule, which is the wrong place for the one thing everyone includes.

**The include rule becomes**

    core ← physics ← model ← equilibrium ← regression ← surface ← bindings

with two clauses the old rule did not need, both checker rule 10 (§8.3): **a `core/` file
includes only `core/` and the standard library** — except `scalar.hpp`, which exists to
include CppAD and is the one core file allowed to — and **`Capture` and `active()` are named
only under `core/` and `bindings/`**. The record and the trace are `double`, `std::size_t`,
`const char*`, `std::optional` and `std::vector`; a solver hands them numbers it already has.
*Derived:* every field in §4 and §5 is one of those types; the layer needs no matrix.

The move is mechanical and *measured*: 19 include lines under `native/` (`parameters.hpp` 5,
`physical_constants.hpp` 6, `scalar.hpp` 3, `numerical_policy.hpp` 5), one in `bindings/`
for the builder, two path strings and the `LAYERS` list in the checker, and one source path in
`meson.build`. Nothing in `tests/` includes any of the four directly.

## 3. Numerical policy, absorbed

### 3.1 The entry

```cpp
namespace epcsaft::policy {
enum class Category { ThermodynamicConstraint = 2, AlgorithmTolerance = 3 };

struct Tolerance final {
    const char* name;                      // equals the identifier; checker rule 8
    double value;
    Category category;
    const char* source;
    const Tolerance* upstream = nullptr;   // the entry this one is derived from; null for a root
    double factor = 1.0;                   // the declared amplification hypothesis; 1 for a root
};

inline constexpr Tolerance ASSOCIATION_RESIDUAL{
    "ASSOCIATION_RESIDUAL", 2.0e-13, Category::AlgorithmTolerance,
    "carried unchanged; Newton reaches it in 3-8 steps from a warm start"};
}
```

Three changes to what exists (*measured:* 27 entries, all `inline const`):

1. **`inline constexpr`.** Required by `derived()`, and by checker rule 5, whose regex matches
   `constexpr` only (*measured*, `CONSTEXPR = re.compile(...)`) and so never saw a policy
   entry — which is how `EOS-AUDIT.md` C.1's two entries with one meaning got past the rule
   written to catch them.
2. **A `name` equal to the identifier**, so a record can point at the entry it was tested
   against and a trace can print it. Rule 8 checks the string against the identifier.
3. **`upstream` and `factor` on the entry.** `TOLERANCES.md` §4 promised that a declared
   factor "beside a measured `1e4` is a contradiction the result can print", then discarded
   both in `derived()`'s return. Two defaulted fields; the 26 root entries are unchanged in form.

Categories stay 2 and 3, 4 empty. `EOS-AUDIT.md` C.9 (six entries sourced to a file name) is
carried as an open item.

**Delivered accuracy is stated on the entry.** `DENSITY-ROOT.md` §7 (*measured*): the density
tolerance is on the *reduced* residual, so relative pressure error is `|r|/z` with
`z ≈ 7e-5…7e-4` on a dense liquid root; a reduced `1e-12` buys about `1e-8` relative,
worst measured `1.7e-8`. That sentence goes in `DENSITY_ROOT_RESIDUAL`'s source string,
because "1e-12" reads tighter than what a caller receives.

### 3.2 `derived()`, kept, with two revisions

```cpp
constexpr Tolerance derived(const Tolerance& upstream, double factor,
                            const char* name, const char* why) {
    if (!(factor >= 1.0)) throw "derived tolerance tighter than its input";
    return {name, upstream.value * factor, Category::AlgorithmTolerance, why, &upstream, factor};
}
inline constexpr Tolerance DENSITY_ROOT_RESIDUAL = derived(ASSOCIATION_RESIDUAL, 5.0,
    "DENSITY_ROOT_RESIDUAL", "bracketed Newton on the reduced pressure residual ... ");
```

**Kept**, for `TOLERANCES.md` §4's reason: each dependency is one greppable edge, an
upstream change moves every downstream value by construction, and the factor is a hypothesis
the runtime record falsifies. Revised: **`constexpr`, not `consteval`** — the build is
`cpp_std=c++17` (*measured*, `meson.build`), and a `constexpr` function that throws on the
failing branch is ill-formed when the call must be a constant expression, which an
`inline constexpr` initialiser is, so the inversion still refuses to compile. **One real edge
at Stage 5:** `DENSITY_ROOT_RESIDUAL` is exactly `5 × ASSOCIATION_RESIDUAL`, the density root
is nested over association-converged evaluations, and `DENSITY-ROOT.md` says the association
contribution to the residual floor is unmeasured; after 04.repair the record's
`floor / ASSOCIATION_RESIDUAL` is that measurement, read against the declared 5. No other edge
is invented to have edges. What it cannot catch is unchanged: a wrong factor. The record is
the check on that.

### 3.3 Two new entries

`FLOOR_MARGIN` (`10.0`, *guessed*: a decade of slack across the first-order bound's norm
inequalities, to tighten once measured ratios exist) and `STALL_WINDOW` (`5.0`, *guessed*:
below Ipopt's `acceptable_iter` default of 15). From `TOLERANCES.md` §5, unchanged. Neither
changes an answer. Rule 6 is unchanged; the exempt path moves.

## 4. The declared-choice contract

### 4.1 The type

```cpp
namespace epcsaft::core {
struct MissingDeclaration final : std::invalid_argument {
    explicit MissingDeclaration(const char* name);       // "missing declaration: <name>"
};
template <typename T> class Declared final {
public:
    constexpr Declared() = default;                       // undeclared
    constexpr Declared(T value) : value_(value) {}
    constexpr bool declared() const noexcept;
    constexpr bool is(T value) const noexcept;            // declared and equal; never throws
    const T& require(const char* name) const;            // throws MissingDeclaration(name)
private:
    std::optional<T> value_;
};
}
```

Twenty lines. What it buys over `std::optional` is the same thing `derived()` buys over a
literal: **the hidden default becomes unwritable.** A kernel may ask `is(ExcludeSameSign)`
and get `false` for an undeclared set, which is the correct physics for a neutral fluid; it
cannot read a value that was never declared without naming it, and the name is what the
caller sees. An enum with an `Unspecified` member is the alternative, and it is the one the
repair brief chose for `IonicDispersion` — it works, but it puts "nothing" inside a type whose
other members are physics, and every reader must remember which member is the hole. Rule 11
(§8.3) forbids `Unspecified` in `core/parameters.hpp` and requires every enum-typed field of
`Parameters` to be wrapped.

### 4.2 Where it applies at Stage 5

The three selectors in `Parameters` (*measured*, `parameters.hpp:20-69`; `ionic_dispersion`
arrives with 04.repair):

| selector | today | after | constructor rule (`helmholtz_model.cpp`) |
|---|---|---|---|
| `electrolyte_form` | `ElectrolyteForm` defaulted to `None` | `Declared<ElectrolyteForm>` | declared iff a charged component exists; `None` stays a *declared* choice of no electrostatic term |
| `permittivity_form` | `PermittivityForm` defaulted to `None` | `Declared<PermittivityForm>` | declared iff `electrolyte_form` is declared |
| `ionic_dispersion` | (04.repair: `IonicDispersion::Unspecified`) | `Declared<IonicDispersion>` | declared iff a charged component exists — the repair's own rule, written on the type |

The struct defaults are hidden defaults in the exact sense of §1: a hand-built `Parameters`
with ions and no `electrolyte_form` computes no electrostatics and says nothing. Fourteen read
sites in physics and model (*measured:* `electrolyte.hpp` 8, `born.hpp` 1,
`helmholtz_model.cpp` 3, `parameter_builder.cpp` 2) become `.is(...)`; the builder assigns
`Declared` values from the document's declared choices as it does today; the two electrolyte
fixtures already declare both (*measured*, `model_support.hpp:103-104,141-142`). No number
changes.

### 4.3 The surface half

The rule for `src/epcsaft/`: **absence crosses the boundary as absence, and the native side
decides.** Python passes `None`, an empty string, or nothing; it never supplies a physical
value in place of one. Two mechanical consequences: checker rule 12 (canon, §8.3) walks every
function signature in `src/epcsaft/*.py` with `ast` and allows no default other than `None` —
which is red on `stoichiometry=(1, 1)` today and green when Stage 5 removes it, taking that
two-line write from the repair brief; and `parameter_builder.cpp`, now in `surface/`, remains
the one place that turns a document's declared choice into a `Declared` value or rejects by
name. That is what the coordinator's "the contract that says which model is being asked for"
means as code, and why the builder is surface: it is the only file whose job is to make an
omission loud.

**Reported as well as required.** Every declared selector is emitted once per `evaluate` span
as `declared_choices` (the enum values as doubles, in `Parameters` field order), so a trace
says which model was asked for beside how it converged. The output half of the contract
reports the input half; that is the sentence in §1 made concrete.

## 5. The convergence contract

### 5.1 The record

```cpp
namespace epcsaft::core {
enum class Stop : std::uint8_t {
    converged, converged_to_floor, stalled, budget_exhausted, left_domain, singular, infeasible
};
struct Convergence final {
    double residual = 0.0;                        // achieved, in the test's own norm
    const policy::Tolerance* criterion = nullptr; // the entry it was tested against
    std::optional<double> floor;                  // Σ amplified inner tolerances + roundoff; absent = not measurable here yet
    std::optional<double> amplification;          // the measured κ behind `floor`; absent likewise
    std::optional<double> contraction;            // per-step residual ratio over the last STALL_WINDOW steps; absent below the window
    std::size_t iterations = 0;
    std::size_t evaluations = 0;                  // inner evaluations the solve made; reported, never asserted
    Stop stop = Stop::stalled;
};
}
```

Against `TOLERANCES.md` §3: `tolerance` becomes `criterion` (the entry, with its name,
source, upstream and declared factor); `floor` and `amplification` are optional; `contraction`
and `evaluations` are added. An absent floor is the statement that none was measured; a
`0.0` there would be a claim. Same rule as 04.repair's `stable_root_count`.

**Where it goes.** `DensityRoot::convergence`; `State::density_convergence` (absent for a
density-specified state), so `Mixture.state(P=...)` reports it without debug mode; the
association solver's out-parameter `Convergence* record`, replacing the `std::size_t*
iterations` Stage 3 added (*measured*, `association.hpp:175`, one reader at
`model_floor.cpp:172`). Later: the assembly, the driver, the sensitivity, the fit, the
discovery (§7). Success returns the record on the result; failure throws
`ConvergenceFailure`, a `std::runtime_error` whose `what()` is `describe(record)`:
`"budget_exhausted: residual 3.1e-09 against DENSITY_ROOT_RESIDUAL 1.0e-12 after 60
iterations, contraction 0.735"`. The equilibrium layer will return failures rather than throw
(`SPECIFICATION-SPACE.md` §12); the record it carries is this one.

### 5.2 One classifier

```cpp
Stop classify(double residual, const policy::Tolerance& criterion,
              std::optional<double> floor, std::optional<double> contraction, bool budget_hit) noexcept;
```

| condition, in order | stop |
|---|---|
| `residual <= criterion.value` | `converged` |
| `floor` present and `residual <= FLOOR_MARGIN.value * *floor` | `converged_to_floor` |
| `budget_hit` and (`contraction` absent or `< 1.0`) | `budget_exhausted` |
| otherwise | `stalled` |

`left_domain`, `singular`, `infeasible` are set by the solver that detects them. A solver
never computes a residual stop itself; it calls `classify`. `ResidualHistory` is the shared
window statistic: a ring of `STALL_WINDOW + 1` residuals, `contraction() = (r_k /
r_{k−W})^{1/W}`, absent below the window; on `r_k = r_0 q^k` it returns exactly `q`
(*derived*; §8.1's closed form).

### 5.3 What the density root's record can honestly carry, before and after the slope repair

| | Stage 5 (before 04.repair) | after 04.repair |
|---|---|---|
| `residual` | `|r|` at the returned root, exact | same |
| `criterion` | `&DENSITY_ROOT_RESIDUAL` | same |
| `floor` | **absent** | `V · floor_grad[V]`, `floor_grad[V] = |Σ_a (Φ_XX⁻¹ Φ_XV)_a · n_tot m_a / X_a| · ASSOCIATION_RESIDUAL` from the factorisation the repair lifts into the `PressureSlope` order |
| `amplification` | **absent** — the only slope in hand is measured wrong | `1 / |dr/dρ|` with the total derivative |
| `contraction`, `iterations`, `evaluations`, `stop` | present | same |

*Derived:* `r = Z − P/(ρRT)`, `Z = 1 − V Φ_V`, so `δr = V · δΦ_V`. *Measured*
(`DENSITY-ROOT.md` §6): with a correct Newton the achieved residual is p50 `1.17e-13`, max
`9.97e-13` against `1e-12`, so `converged_to_floor` is not expected to fire on the density
root; it is exercised by the classifier's properties, never by a solve. The 102-of-630
pre-repair failures are `budget_exhausted` with `contraction ≈ 0.735`, which the floor
harness prints and no test asserts — a test that a defect exists is a test of this code.

**The association record** (leaf): `residual = ‖H‖∞`, `criterion = &ASSOCIATION_RESIDUAL`,
`floor = ε_mach · S`, `amplification` = the pivot ratio `largest / smallest` of the last
Newton step's LU (*measured:* `solve_linear` computes both, `association.hpp:155-156`),
`contraction`, `stop` from `classify` or `left_domain` where it already throws. All present
at Stage 5.

**`HelmholtzDerivatives::gradient_floor`** (`TOLERANCES.md` §3.2) lands on the `Hessian`
order as the ∞-norm over rows of `(Φ_XX⁻¹ Φ_Xy)ᵀ · diag(n_tot m_a / X_a)` times
`ASSOCIATION_RESIDUAL`; `factorization.solve(cross)` is already formed (*measured*,
`helmholtz_model.cpp:453`), so it is one `w × S` row-sum. Zero at `S = 0`. On the
`PressureSlope` order it is the `V` row only and lands with 04.repair.

## 6. Debug mode

### 6.1 The shape

```cpp
namespace epcsaft::core {
struct Event final { const char* name; std::uint32_t span; std::uint32_t iteration; std::vector<double> values; };
struct Span  final { const char* name; std::optional<std::uint32_t> parent; std::optional<Convergence> convergence; double wall_us; };
struct Trace final { std::vector<Span> spans; std::vector<Event> events; std::optional<std::uint32_t> current; };

Trace* active() noexcept;           // this thread's active trace, or nullptr; named only in core/ and bindings/

class Capture final {               // RAII: installs a fresh Trace as this thread's active trace; restores the previous
public:  Capture(); ~Capture(); const Trace& trace() const noexcept; Trace take(); };

class Scope final {                 // RAII: one span per instrumented function call
public:
    explicit Scope(const char* name) noexcept;   // one thread-local read; opens a span iff a capture is active
    ~Scope();                                     // wall time, parent restored
    bool on() const noexcept;
    void probe(const char* name, double value, std::uint32_t iteration = 0);
    void probe(const char* name, const double* values, std::size_t count, std::uint32_t iteration = 0);
    void probe(const char* name, const std::vector<double>& values, std::uint32_t iteration = 0);
    void finish(const Convergence& record);
};
}
```

**Turning it on.** Natively `core::Capture capture;` around the call; from Python

```python
with epcsaft.trace() as t:
    state = mixture.state(313.15, P=1e5, x=x, phase="liquid")
for span in t.spans: print(span.name, span.parent, span.convergence, span.wall_us)
for event in t.events: print(event.name, event.span, event.iteration, event.values)
```

No recompilation, no probe, no flag on any method. The unit is one call because the capture
is a lexical scope around it.

**Composition.** `Trace::current` is the open span; a `Scope` records `parent = current`,
becomes current, restores its parent on destruction. The density root opens `solve_density`;
each `evaluate` inside is a child; the association Newton inside each is a grandchild; the
reference ladder's fixed-pressure states nest under `infinite_dilution_reference`; Stage 7's
driver will sit above them. Nothing threads a parameter, so a later stage that opens a span
at the top of its solver is composed. The trace is a tree (§8.2 checks exactly that).

**The one rule that keeps it a debug mode.** *Nothing in the engine reads the trace.* A floor
crosses from `evaluate` to the density root and from the assembly to the driver as a return
value, never by reading an event back. `active()` and `Capture` are named only in `core/` and
`bindings/` (rule 10); engine layers see `Scope`, which is write-only. If the trace were a
data channel, debug mode would change answers and "reported, never asserted" would have
nothing to stand on. §8.2 checks that a captured call returns the same bits as an
uncaptured one.

### 6.2 What it costs when off — measured

At an instrumented site with no capture active: the `Scope` constructor reads one thread-local
pointer; each `probe` tests one member pointer. *Measured* (`scratchpad/bench/`, gcc 13.3.0,
`-O2`, one core, sources in Appendix A):

| site | cost when off | how |
|---|---|---|
| `probe(...)`: member-pointer test, not taken | **0.11 ns** | 2×10⁸ calls of a `noinline` site with and without the check |
| writing a six-field record unconditionally | 0.75 ns | same; the record is the result, written regardless |
| `Scope` constructor: thread-local read in a `-fPIC` shared object, default visibility (`__tls_get_addr` via the PLT — the pybind11 extension's case) | **2.2 ns** | a `.so` built `-fPIC -shared`, 2×10⁸ calls from a separate executable; `objdump` confirms the call |
| the same, hidden visibility | 1.9 ns | same |

Against what the engine does between sites (*measured*, this tree's archive, same harness):

| unit of engine work | wall | sites | off-cost | ratio |
|---|---|---|---|---|
| `evaluate(PressureSlope)`, neutral binary — the cheapest instrumented call | 15.8 µs | 2 `Scope` + ≤ 9 probes | ≈ 5.4 ns | **3.4 × 10⁻⁴** |
| `evaluate(Hessian)`, nine species | 663 µs | 2 `Scope` + ≤ 13 probes | ≈ 5.8 ns | 8.7 × 10⁻⁶ |
| `solve_density` cold, `Liquid`, neutral binary | 724 µs, 48 evaluations | 49 `Scope`, ≈ 250 probes | ≈ 135 ns | 1.9 × 10⁻⁴ |
| `solve_density` cold, `Liquid`, nine species | 12.1 ms | 55 `Scope`, ≈ 300 probes | ≈ 155 ns | 1.3 × 10⁻⁵ |

**Under four parts in ten thousand of the cheapest call, below the floor harness's noise.**
Synthetic-site measurement; small enough that an in-situ measurement cannot separate it from
noise, so the harness reports the *on* cost (§6.4) and this document carries the *off* cost.
The choice that makes it this small: the thread-local read happens once per function in the
`Scope`, and probes test a member. The 2.2 ns is the unrelaxed general-dynamic case and so the
ceiling.

### 6.3 What is emitted, per layer

Names are quantities, never routes: a registered `EqID`, a policy entry name, or a field name
of the result the function returns. Rule 7 already scans string literals in solver files
(*measured:* `without_comments` strips comments only), so a route word in an event name goes
red with no new rule.

| span | opened by | events | record |
|---|---|---|---|
| `association_newton` | `solve_site_fractions` (physics, Stage 5) | per iteration `association_residual_norm`; at exit `site_fractions`, `association_pivot_ratio` | §5.3 |
| `evaluate` | `HelmholtzModel::evaluate` (model, Stage 5) | `declared_choices`, `recording`, `predicate_signature`, `residual_helmholtz_extensive`, `gradient`, `volume_second_derivative`, `hessian` (Hessian order), `gradient_floor`, `site_fractions` | none: not a solve |
| `solve_density` | (model, Stage 5) | per scan point `density_scan_point` [ρ, r, dr/dρ]; per bracket `density_bracket`; per Newton step `density_newton_step` [ρ, r, dr/dρ, candidate]; per refined root `density_root` [ρ, r, dr/dρ, stable]; at exit `density_root`, `stable_root_count` | §5.3 |
| `density_anchor` | the anchored attempt, child of `solve_density` (Stage 5; 04.repair adds `far_side_residual` per grid point) | per step as above | its own record: `left_domain` or `stalled` when it falls to the cold path |
| `infinite_dilution_reference` | (model, Stage 5) | per ladder point `reference_ladder_point` [m, ln φ₁, ln φ₂]; `reference_convergence_error`; `reference_witness_error` | none: a bias, not a stop (`TOLERANCES.md` §3.6) |
| `compile` | `request_compile.cpp` (Stage 6) | `active_rows`, `degrees_of_freedom` [k, f, rank, nullity], `assembly_floor` per evaluation | none: no iteration |
| `solve_equilibrium` | the driver (Stage 7) | per Ipopt iteration from `intermediate_callback`: `primal_infeasibility`, `dual_infeasibility`, `barrier_parameter`, `primal_step`, `dual_step`, `regularization`; running `assembly_floor` | `residual = inf_pr` (k ≥ 1) or `max(inf_pr, inf_du)` (k = 0); `floor` from the assembly; `amplification = κ(F_y)` from the sensitivity; `stop` from `classify` plus Ipopt's status for `infeasible`/`singular` |
| `implicit_sensitivity` | (Stage 7) | `kkt_smallest_singular_value`, `kkt_condition` | none |
| `fit` | `least_squares.cpp` (Stage 10) | per Ceres iteration `cost`, `gradient_norm`, `trust_region_radius`; observation solves nest automatically | `residual` = gradient norm at termination; `floor = Σ_obs ‖∂obs/∂u‖ κ_sys · residual_obs` |
| `find_phases` | (Stage 11) | per start `tangent_plane_distance`; per major iteration `dual_bound`, `candidate_gap` | `budget_exhausted` here is what `TRANSFER.md` Stage 11 calls `search_budget_exhausted` |

Everything in this table is **reported and never asserted**. The only assertions in the
design are §8's four properties.

### 6.4 What it costs when on — estimated, to be reported

One `std::vector<double>` per event. A nine-species `evaluate(Hessian)` emits about 170
doubles, 1.4 KB; a cold nine-species density solve under 100 KB; a reactive solve of 300
evaluations under 1 MB. *Guessed*; the floor harness times one captured cold solve beside the
uncaptured one. Memory grows for the life of the capture and is unbounded by design.

### 6.5 Crossing the pybind11 boundary

`bindings/diagnostics_bindings.cpp`: `Stop`, `Category` as enums; `Tolerance`,
`Convergence`, `Event`, `Span`, `Trace` with `def_readonly`; `Convergence::criterion` as a
reference to a static `constexpr` object; a `Capture` wrapper with `__enter__`/`__exit__`;
`ConvergenceFailure` as `ConvergenceError`, subclass of `EosError`, registered after it so
pybind11 tries the derived type first; `MissingDeclaration` as `DeclarationError`, subclass of
`ParameterError`; `State.density_convergence`. Python does no arithmetic: the ratios
`residual / criterion.value`, `residual / floor`, `amplification / criterion.factor` are
native properties. No renderer; a `for` loop over `t.events` is twelve lines of Python, not a
probe.

One thread, one trace. The tape's mutex serialises `evaluate` anyway (*measured*,
`helmholtz_model.cpp:307`); no entry point spawns a thread.

### 6.6 What keeps it from becoming the instrumentation `TESTING.md` bans

One switch, no filters; names are quantities; nothing asserts on it; nothing reads it; and
the carry warning is taken — the superseded tree's `_core._instrumentation()` counters and its
name–value evidence list are not carried, the trace is written from §6.1, and Stage 3's
`SweepCounters` stay where they are with the trace reading nothing from them.

## 7. How each layer depends on the core, and what every brief must populate

| layer | includes from `core/` | returns | reads (as return values, never from the trace) |
|---|---|---|---|
| physics | constants, policy, `Parameters` with `Declared`, scalar; `Convergence`, `classify`, `ResidualHistory`; `Scope` | the association record | — |
| model | the same | `DensityRoot::convergence`, `State::density_convergence`, `gradient_floor` | the association record for its floor |
| equilibrium | the same | the assembly floor, the driver's record, κ_sys | `gradient_floor` per phase |
| regression | the same | the fit's record | each observation's record and κ_sys |
| surface | `Declared`, `MissingDeclaration`, `Parameters` | `Parameters` with every selector declared or a rejection by name | — |
| bindings | everything; with `core/` the only place `Capture`/`active()` are named | — | — |

**A later stage cannot skip it**, by two mechanisms. Static: **rule 9** — a file under
`native/` outside `core/` and `bindings/` whose comment-stripped text reads a budget entry
(`policy::\w+_(MAX_ITERATIONS|BUDGET)\b`) must name `core::Scope` and `Convergence`; a budget
is what makes a loop a solver, every planned solver has one, and a solver without a budget
entry is an unbounded loop that rule 6 already refuses. Procedural: **each of Stages 6–11
carries a `diagnostics` array** in its JSON, one entry per solver or assembly it writes, of
the shape in the Stage 5 brief, and the preflight refuses a brief whose `writes` include a
file matching rule 9's predicate with no entry naming it. The `check` field is always the
same shape — *the record's residual equals the residual recomputed from the returned answer
through the public equations* — a property of the pair (answer, record) that survives any
reimplementation. The preflight script is outside the repository (*measured*), so this is a
request to the coordinator.

## 8. What is tested, and what is only reported

### 8.1 Properties of the classifier, the window and the two guards (`diagnostics_checks.cpp`)

- **Partition.** Over a grid of `(residual, criterion.value, floor ∈ {absent, values},
  contraction ∈ {absent, <1, ≥1}, budget_hit)`: exactly one of the four residual stops;
  `converged` iff `residual <= value`; `converged_to_floor` only with `floor` present and
  `value < residual <= FLOOR_MARGIN · floor`; `budget_exhausted` only when `budget_hit`;
  monotone in `residual`.
- **Contraction closed form.** `ResidualHistory` fed `r_0 q^k` returns `q` to `1e-12` once
  the window is full, absent before.
- **`derived()` at runtime** throws on `factor < 1`; `derived_inversion.cpp`, a six-line
  translation unit declaring `derived(ASSOCIATION_RESIDUAL, 0.5, ...)` as an `inline
  constexpr` entry, **fails to compile**; the gate is the compiler's exit status negated.
- **`Declared`**: `require` on an undeclared value throws `MissingDeclaration` whose message
  names it; `is` on an undeclared value is `false`; the model constructor refuses an ion set
  with no `electrolyte_form` and a neutral set with one, each by name.

### 8.2 Properties of a record and a trace on real solves

- **The record matches the answer.** For every fixed-pressure solve `model_checks.cpp`
  already runs, `state.density_convergence.residual` equals `|Z − P/(ρRT)|` recomputed from
  `state(T, ρ*, x)` to roundoff, and `stop == converged` implies `residual <=
  criterion.value`. For the association solver, `record.residual` equals
  `‖association_residuals(X*)‖∞` recomputed from the returned fractions.
- **The trace is a tree that changes nothing.** A captured `state_at_pressure` on each
  fixture returns a `State` bit-identical to the uncaptured call, and its trace's `parent`
  indices each precede their span with exactly one root. Shape and non-interference; no
  check names a span, an event or a count.

### 8.3 Checker rules (`tools/check-engine-architecture.py`, +80 lines)

| rule | checks | red on |
|---|---|---|
| 3, revised | `LAYERS = ["core", "physics", "model", "equilibrium", "regression", "surface", "bindings"]` | an upward include |
| 5, 6, revised | `CONSTANTS_HOME` and the exempt paths move to `core/` | unchanged behaviour |
| 8 | every `inline constexpr Tolerance NAME{"NAME", ...}` or `= derived(..., "NAME", ...)` names itself; no `inline const Tolerance` remains | a name string differing from its identifier |
| 9 | a file reading a budget entry names `core::Scope` and `Convergence` | a solver with no record |
| 10 | a `core/` file has no `"..."` include outside `core/` and no numerical-library include (`scalar.hpp` may include `<cppad`); `Capture` and `active(` appear only under `core/` and `bindings/` | the layer reaching up; an engine layer reading the trace |
| 11 | `Unspecified` does not appear in `core/parameters.hpp`; every field of `Parameters` whose type is an `enum class` declared in that file is `Declared<...>` | a hole inside an enum; a bare selector |
| 12 | `ast` over `src/epcsaft/*.py`: every function-signature default is the constant `None` | a physical default on the surface |

### 8.4 Reported, never asserted (`model_floor.cpp`)

Each density case's record via `describe`; the association record on case C; one captured
cold nine-species solve beside the uncaptured one: wall, spans, events, bytes. On the
pre-repair tree the associating binary at 313.15 K, 1e5 Pa `Liquid` prints
`budget_exhausted` with its contraction — the number `TOLERANCES.md` §1 said the old tree
never recorded, now recorded before the fix so the fix can be seen to change it.

Test lines: `diagnostics_checks.cpp` 170–270, `derived_inversion.cpp` 8, `model_floor.cpp`
+30, `tests/eos/test_diagnostics.py` 30–50. Inside `TESTING.md`'s 1,650.

## 9. Disposition of `TOLERANCES.md`

| section | disposition | why |
|---|---|---|
| §1 the failure | kept as evidence | the measured shape the record replaces |
| §2 the principle | kept by reference | the floor derivation is right |
| §3 the record | revised → §5.1 | optionals; `criterion`; `contraction`, `evaluations` |
| §3.1–3.3 the three floors | kept; §5.3 dates each | the density floor waits for 04.repair |
| §3.4–3.5 driver and fit | kept as §6.3 and §7 commitments | the one-Ipopt-driver revision of 15482d8 taken as read |
| §3.6 bias | kept | reported on the result, never a stop |
| §3.7 cost | superseded by §6.2 | "guessed: under one percent" is now measured at under 4 × 10⁻⁴ |
| §3.8 where it lives | superseded by §2 | the policy header is not in physics |
| §3.9 one model | kept | `gradient_floor = 0` at `S = 0` |
| §4 the static guard | revised → §3.2 | `constexpr`; the entry records `upstream` and `factor`; one real edge |
| §5 named stops | revised → §5.2 | one `classify`; `floor` optional; three stops set by the solver |
| §6 couplings, §7 prior art | kept by reference | commitments and evidence |
| §8, Appendix A | discarded | argued against a framing that no longer exists; a patch to the superseded tree never applied |

## 10. The decisions `CATEGORIES.md` left open, and the sequence

### 10.1 Should the declared-choice contract live in the same core as the convergence contract?

**Yes**, and §4 is the mechanism. The argument is not filing: both contracts are a type that
makes the silent outcome unwritable (`Declared::require`, `derived()`), a report that makes
the loud outcome legible (`MissingDeclaration` by name, `describe(record)`), and a checker
rule that keeps later stages from opting out (11, 9). One core, two contracts, one discipline.
And they meet at one place: the `evaluate` span reports `declared_choices` beside every
record, so a trace answers "which model, and how did it converge" in one read.

### 10.2 `parameter_builder` — surface. Agreed.

It is the only file whose job is to make an omission loud, and filing it under `model/` is
what let a lost selector read as a mapping detail. It moves to `native/surface/`, rank
between `regression` and `bindings`, so `bindings/` may include it and nothing in the engine
can. Two files, so the no-single-file rule holds; `problem.py`'s native loader
(`SPECIFICATION-SPACE.md` §11) is its future neighbour.

### 10.3 `helmholtz_tape.hpp` — engine, and the EqID test is right once the annotation gap is closed.

The tape records `Φ(y; X) = n_tot [a_hc + a_disp + a_dh + a_born] + n_tot Q(y; X)` with the
association potential `Q = Σ m_a (ln X_a − X_a + 1) − (ρ_N/2) Σ m_a X_a B_a` and the Schur
complement `Φ_yy − Φ_yX Φ_XX⁻¹ Φ_Xy` (*measured*, Stage 3 brief §2.1). Those are equations
about the fluid — the potential *is* `a_assoc` at the site solution and the complement *is*
the total second derivative — and they carry no `EqID` because nobody annotated them, not
because there is nothing to annotate. Annotate `Q` and the complement (two lines) and the
measurement classifies the file as engine, where it belongs: it is shaped around one
Helmholtz model and cannot serve another. Core is what every layer reads; nothing above
`model/` reads the tape.

### 10.4 `physical_constants.hpp` — core.

Reference data, not equations; read by every layer; already governed by rule 5's
single-home clause. The category test is "may depend on nothing", which it satisfies; the
engine reading it is the permitted direction. Facts about the world that are *data* are what
packets are too, and packets are not engine.

### 10.5 The analyses rule, and its two violations.

Rule as proposed, with one clause: **nothing depends on an analysis at runtime or in a gate;
an analysis becomes a dependency only by promotion into a packet through the validator, and
a `source_file` citation with a hash is provenance, not a dependency.** Two instances
*measured* today:

- The LiBr digitized curve. The repair brief's §2.7 transcribes the 27 rows into
  `oracle.json` with `source_file` and its sha256; the gate reads `oracle.json`. Once landed
  that is provenance, not a read, and no exception is needed. It is not yet a packet; whether
  a digitization deserves the packet contract (provenance, lock, validator) is the promotion
  question, and the answer is yes if any second gate ever reads it — one reader lives in the
  oracle, two readers need a packet.
- `engine/tests/oracle/profile_entry_points.py:22` reads
  `analyses/2026-mea-parameter-estimation/inputs` by an absolute path into the *other*
  checkout, and imports `unit_registry` and `_core._instrumentation`, neither of which
  exists in this tree. A canon file that is dead and points outside the repository: the
  "nothing has swept the rules" finding, with a file name. Delete it or move it to the
  analysis it serves; it is not a Stage 5 write, it is a triage item.

### 10.6 `promote_mea_active_bundle.py` — canon tooling, and the rule names it.

It is the promotion boundary itself: the one sanctioned way a campaign's result becomes data
at rest. That makes it a sibling of `validate_packets.py`, not a tool of the analyses, and
the clause in 10.5 is what places it: the canon owns the crossing, so the crossing is checked
where the packet contract is checked. It resists the table only while the table has no row
for "the crossing"; the clause is that row.

### 10.7 One sequence, or per-category tracks?

**One sequence with a category column, not per-category tracks.** Three reasons, the first
decisive:

1. **A gate is one command on one tree.** Stage 4's oracle gate needs the engine, the
   surface and the canon at once; every later gate does too. Per-category tracks would have
   to synchronise at every gate, which is every stage, and the tracks collapse back into
   the sequence with more bookkeeping and the same order.
2. **One number names a stage, its brief, its branch and its worktree** (`ORCHESTRATION.md`).
   The crosswalk in `TRANSFER.md` exists because there were briefly two numberings; a track
   per category is five.
3. **What a category *can* give the sequence is a column.** Each stage row declares what it
   writes per category — engine / core / surface / canon — and the preflight refuses a write
   into a category the row did not declare. That is the thing that would have caught the
   dispersion selector: a Stage 2 row declaring engine writes only, with a registered case
   split whose selector belonged to the surface column, has a hole in the surface column that
   a reader can see. Bookkeeping that names an owner is the kind that helps.

The exception already in use stays: the lettered items P and C are canon and analyses work
whose gates do not need the engine tree, and they synchronise at named stages. That is a
per-category track for exactly the two categories where it costs nothing, and it works
because their gates are independent of the build.

### 10.8 Ordering: the EOS repair waits for Stage 5.

The repair cannot report convergence without the record (§1); the density loop should be
instrumented once, and the repair's rewrite of the anchored path should be written into a
span from the start; the `IonicDispersion` selector should be born `Declared` rather than
retrofitted from `Unspecified`; and Stage 5 changes no number, so it can be accepted on the
tree that exists. The repair is re-based on the Stage 5 merge with the delta in the brief's
Appendix A. Retrofitting the other way round would mean editing `density.cpp`,
`helmholtz_model.cpp`, `parameters.hpp`, `helmholtz.hpp` and the builder twice within a week.

### 10.9 What moves

| document | change |
|---|---|
| `TRANSFER.md` row 5 | "diagnostics core" → "the core: `native/core/` (seven files), `native/surface/`, the two contracts, `epcsaft.trace()`; design `CORE.md`"; "Stage 5 does not start until the kernel is fixed" reversed; §"Stage 5" replaced by §2 and §8 here; a category column added to the execution table |
| `target-layout.md` | `native/core/` (7 files, ~800) and `native/surface/` (2, 425) added; `physics/` becomes four files (~1,030); `model/` loses the builder; the undecided-home bullet closed; envelope +~500 |
| `README.md` | `TOLERANCES.md` row → `CORE.md`, stages "5, and the `diagnostics` block of 6–11" |
| `CATEGORIES.md` | the five decisions of §10.1–10.6 recorded; the analyses clause |
| `tools/check-engine-architecture.py` | ranks, homes, rules 8–12 |
| `briefs/stage-04-repair.brief.{md,json}` | re-based: paths for the four moved headers and the builder; `Declared<IonicDispersion>` in place of `Unspecified`; the `stoichiometry` default removal and the `ELECTROLYTE_SERIES_THRESHOLD` deletion re-pathed; the density record's floor and amplification; the `density_anchor` far-side probes; a `diagnostics` block |
| `briefs/stage-06…11` | the `diagnostics` block and the category column |
| `ORCHESTRATION.md` §"Preflight" | the rule-9 / `diagnostics`-block refusal and the category-column refusal |

Line-count honesty: the earlier Stage 5 was 260 lines; this is about 1,300 written and about
+560 net in `native/`, of which ~450 is the trace, the two contracts and the checker, and the
rest is moves.

## 11. What this design does not cover, and missing facts

- **The in-situ off-cost.** §6.2 is a synthetic-site measurement, below the floor harness's
  noise by construction.
- **Bounded memory while capturing** (§6.4). Unbounded by design; a policy entry and a named
  stop if that is ever a problem.
- **Threads.** One trace per thread; nothing today spawns one.
- **The Ipopt and Ceres rows of §6.3** are commitments; whether `intermediate_callback`
  exposes the unscaled infeasibilities as assumed is *derived* from its documentation.
- **`converged_to_floor` on a real solve.** No solver has a measured case at its floor.
- **Event-name discipline** beyond rule 7 is a brief rule; a checker rule would need the
  equation registry, which the checker does not read.
- **`DensityBranch::Unspecified`** (`helmholtz_model.hpp:55`) is a *choice* meaning "either
  branch", not a hole, so rule 11 does not reach it; the name is misleading and `Either`
  would say what it means. A rename, not a Stage 5 write.
- **Missing fact:** whether the coordinator's preflight can take the rule-9, `diagnostics`-
  block and category-column refusals (§7, §10.7).
- **Missing fact:** the association contribution to the density residual floor, which
  `DENSITY_ROOT_RESIDUAL.factor = 5` will be read against after 04.repair.
- **Missing fact:** a source for `REFERENCE_CONVERGENCE_TOLERANCE` and its siblings.
- **Missing fact:** whether the LiBr digitization is read by any gate other than the repair's
  oracle case; if a second reader appears it is promoted to a packet (§10.5).

## Appendix A. The off-cost harness

`scratchpad/bench/hook_cost.cpp` (links `engine/build/greenfield/libepcsaft_native.a` and
`tests/native/model_support.hpp`) and `scratchpad/bench/tls_lib.cpp` + `tls_main.cpp` (the
thread-local read in a `-fPIC` shared object). Output on 2026-09-02, AMD Ryzen 5 5500,
gcc 13.3.0, `OMP_NUM_THREADS=1`:

```
per-site ns: bare=1.715 record_write=2.469 record_write+null_check=2.578  (hook delta 0.109 ns)
neutral_binary:     evaluate(Hessian)=37.5 us   evaluate(PressureSlope)=15.8 us   solve_density(cold,Liquid)=723.9 us   forward_sweeps/solve=96
associating_binary: evaluate(Hessian)=82.8 us   evaluate(PressureSlope)=38.4 us   solve_density(cold,Liquid)=3949.6 us  forward_sweeps/solve=232
nine_species:       evaluate(Hessian)=663.3 us  evaluate(PressureSlope)=111.6 us  solve_density(cold,Liquid)=12089.6 us forward_sweeps/solve=216
per-call ns in a -fPIC .so: plain=2.565  tls_default_visibility=4.749 (delta 2.184)  tls_hidden=4.457 (delta 1.892)
```

Build lines:

```
g++ -O2 -std=c++17 -DNDEBUG -Iengine/native -Iengine/tests/native -Iengine/build/subprojects/CppAD-20250000.3/include \
    -Iengine/build/greenfield/build/subprojects/CppAD-20250000.3 -I/usr/include/eigen3 hook_cost.cpp \
    engine/build/greenfield/libepcsaft_native.a -lpthread -o hook_cost
g++ -O2 -std=c++17 -fPIC -shared tls_lib.cpp -o libtls.so && g++ -O2 -std=c++17 tls_main.cpp -L. -ltls -Wl,-rpath,"$PWD" -o tls_bench
```

The density solves are the pre-repair solver (frozen-site slope), which is why the associating
binary needs 232 forward sweeps; they are denominators for §6.2, not performance claims.

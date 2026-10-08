# Unified ePC-SAFT Engine

One installed Python distribution owns the Engine. Current state/equilibrium imports:

```python
import epcsaft
from epcsaft import equilibrium
```

Consumers import `epcsaft` from their environment. They do not run Meson or
load modules from this source tree.

## Canonical build and installation

From the repository root, build and install the current source non-editably:

```bash
bash .codex/environments/setup.sh
```

Setup keeps the checkout-local Ninja build warm, uses `ccache` when available,
and reuses a verified machine-local wheel keyed by source, ABI, platform, and
toolchain and dependency contents. It builds and installs before any test
starts. The root [`README.md`](../README.md) owns routine, full, native,
sanitizer, and scientific validation commands.

## Native development

`tools/native-run <program> [case|--list]` rebuilds one native program in the
checkout's Meson build directory (configured by setup with the shared HiGHS) and
runs it without building the wheel. Do not use `meson install`; the wheel is the
installation boundary.

`bash engine/tools/run-native-tests.sh` compiles and runs the fast suites
`contributions`, `properties`, `diagnostics`, `specification`, `surface_example`
and `regression`; name other tests to run them. `tools/check-affected.sh` selects
the equilibrium, boundary, Explicit and `association-comparison` tests when their
code changes (root [`README.md`](../README.md) "Checks").

## C++ usage and debugging

`tests/native/surface_example.cpp` is an ordinary Python-free caller of the
native archive. It includes only `surface/parameter_builder.hpp` and
`equilibrium/problem.hpp` and shows, in order: a `ParameterDocument` with
explicit units (species index = position in `component_ids`), one reusable
`HelmholtzModel`, a full `State` in the T-P and T-rho charts, selected
observables with first- and second-order actions, a bubble point and a
two-phase flash, exact solved-state actions with the caller-formed fugacity
`a_i R T (1 mol/m^3)`, a partial result without thermochemistry, and two
inspectable failures. Its mixture is manufactured. It returns nonzero when one
of its identities fails, and it links `epcsaft_native`; there is no installed
SDK or ABI.

```bash
tools/heavy tools/native-run surface_example   # prints every value with its unit
```

For a debugger, configure an unstripped debug tree once and run the same target:

```bash
meson setup build/engine-debug engine -Dbuildtype=debug -Dstrip=false
meson compile -C build/engine-debug surface_example
gdb --args build/engine-debug/surface_example
```

## Python counterpart

The installed package transports the same native records with the same meanings
and SI molar units; Python performs no thermodynamic arithmetic. Thermal
intervals follow `Parameters.component_ids` order and every record field is
required. A model copies its complete physical inputs at construction; changing
an input record afterwards never alters the prepared model. Before construction,
reads of primitive vectors (`coefficients`, `indices`, `composition`) are Python
copies, so `coefficients[0] = ...` or `indices.append(...)` does not change the
record: replace the whole field. Bound record elements stay live, so
`record.component_intervals[0][0].lower = ...` and `interval.correlation.source_pressure = ...`
do change that input record.

```python
from epcsaft import IdealCorrelation, IdealInterval, IdealPolynomial, Mixture, ThermochemistryRecord
from epcsaft.equilibrium import EOS_STANDARD_STATE_ID, Reaction, ReactionLogPolynomial

cp = IdealCorrelation(IdealPolynomial([80.0, 0.0], 300.0, 1000.0, 50.0), source_pressure=1.0e5)
record = ThermochemistryRecord(1.0e5, [[IdealInterval(250.0, 450.0, True, True, cp)], ...])
model = Mixture(parameters, thermochemistry=record)           # omit the record for residual-only States
state = model.state(353.15, P=1.0e5, x=(0.25, 0.75), phase="liquid")
state.total_enthalpy, state.ideal_entropy                     # None when the record does not cover T
ln_k = ReactionLogPolynomial(a, b, c, d, Tref, Tmin, Tmax, EOS_STANDARD_STATE_ID)
Reaction({"A": -1.0, "B": 1.0}, correlation=ln_k)             # exactly one of K, affinity or ln K(T)
```

Actions are order one to three contractions along caller directions and keep the
observable's unit. A state direction holds K, mol/m^3 or Pa, a zero-sum
composition increment and native family/slot parameter increments; a solved-state
direction holds K, Pa, absolute feed mol (component width) and one increment per
active parameter (active width). A free T or P is an output: its entry is zero.

```python
from epcsaft import ActiveParameter, DensityActionDirection, ParameterActionDirection, ParameterFamily
from epcsaft.equilibrium import (SolvedStateActionDirection, SolvedStateActionRequest, SolvedStateObservable,
                                 SolvedStateObservableKind as Kind, compile_problem, solved_state_actions)

pair = ActiveParameter(ParameterFamily.AssociationEnergyOverK, (1, 4))  # both slots of site pair (0, 1)
model = Mixture(parameters, thermochemistry=record, active_parameters=[pair])
move = DensityActionDirection(1.0, 0.0, (0.0, 0.0), [ParameterActionDirection(pair.family, pair.indices, 0.03)])
batch = model.density_actions(T, rho, x=x, observables=observables, directions=[move])
batch.results[0].action, batch.results[0].status, batch.diagnostics  # None beside a typed status when unavailable

central = solve_equilibrium(model, problem)
compiled = compile_problem(model, problem)                    # one compile; no solve is repeated
warm = SolvedStateActionDirection(1.0, 0.0, [0.0, 0.0], [0.0])
request = SolvedStateActionRequest([SolvedStateObservable(Kind.PhaseLogActivity, 0, 1)], [warm])
row = solved_state_actions(compiled, central, request).results[0]
fugacity = math.exp(row.value) * GAS_CONSTANT_J_PER_MOL_K * central.temperature[0]  # f_i = a_i R T (1 mol/m^3), Pa
```

`solve_equilibrium` has no time limit; `solve_equilibrium_starts` tries declared
starts in order under shared work limits and selects the first that meets the criterion:

```python
limits = equilibrium.EquilibriumWorkLimits(seconds=60.0, attempt_seconds=20.0)
tried = equilibrium.solve_equilibrium_starts(model, problem, [equilibrium.EquilibriumStart()], limits)
result = None if tried.selected is None else tried.results[tried.selected]  # tried.causes says why not
```

A default `EquilibriumStart` skips the speciated seed that `solve_equilibrium` tries
first, so it can take a different path and work count to the same state; it does not replay `solve_equilibrium`.

Invalid declarations and inputs raise `ParameterError` (a `ValueError`) before iterating. A
solve that runs reports non-convergence, structural infeasibility and rejected density
branches as `success=False` with `message` and `convergence.stop`, but a failure to
evaluate the returned point still raises `EosError` or `ParameterError`; handle both.
`solve_equilibrium_starts` records an attempt's exception in `causes` instead of raising.

For a fixed-T/P, declared two-free-liquid equilibrium, `equilibrium.asl_record`
constructs the immutable record consumed by the installed
`epcsaft_equilibrium_phase_component_amount` external function:

```python
record = equilibrium.asl_record(model, problem, continuation, feed_offset,
                               feed_matrix, anchor_controls, identity_component, minimum_gap)
# Pyomo ExternalFunction: (*controls, record, "0" or "1", str(component_index))
```

Controls and returned amounts are mol on the caller's declared batch basis.
The feed is `n = feed_offset + feed_matrix @ controls`, with component rows and
control columns. Offset and columns must be charge neutral; physical feed and
controls are nonnegative, including zero incoming reactive complexes. The anchor
must reproduce the supplied absolute feed. The selected component mole-fraction
gap (phase 0 minus phase 1) must exceed `minimum_gap`; indices and supports follow
the supplied model order. Scalar reaction K or affinity and explicit problem
standard Gibbs references are transported. Version 3 carries explicit base-EOS
and association selections; earlier records are rejected without guessed defaults.
Temperature-dependent reaction
correlations, neutral-reference transformations, active parameters and model
thermochemistry records are rejected rather than silently omitted.

For an observable `y(n)`, the callback returns `J_y B` and `B.T H_y B`, where
`B = feed_matrix`; E5 owns these implicit actions, not a finite-difference
callback. Pyomo owns any additional control products, ratios and process units.
Each new feed continues from the immutable anchor; one exact-iterate cache and
one immutable record are admitted per process. Use separate processes for
different records. Failed central/branch/action checks return nonfinite values
and derivatives with an observable-specific typed error, never a surrogate
fallback. Use an outer Ipopt executable whose C++ symbols do not interpose on
the Engine's Ipopt dependency; the retained live campaign verifies the isolated
official IDAES executable and exact option-file readback.

Diagnostics keep their native meaning. `State.thermochemistry_status` separates a
missing record from a temperature outside its intervals; `density_branch` and
`stable_root_count` exist only after a pressure-root search, so a T-rho state reports
`Unspecified` and `None`. `EquilibriumResult.multipliers` are the dimensionless
reduced-energy balance and charge coordinates of the solution;
`constraint_multipliers` are the optimizer's raw duals of its implemented, scaled
constraints (the non-chemical, non-pressure entries of `rows` in order, or every
row when a phase amount is pinned). Neither is a chemical or electric potential.
Read `requested_tolerance_met` independently of `success` (ALG-EQUILIBRIUM,
ALG-TERMINATION). A floor-limited point that Ipopt reports as successful is accepted
and keeps its state with `success = True` and `requested_tolerance_met = False`. A
tiny-step floor-limited point is rejected: `success = False`, finite row residuals, the
convergence record and the stopping cause remain, and the candidate state vectors are
cleared: in Python, `success` (equivalently non-empty state vectors) is the availability
predicate, the native `has_state()`.

Create and install wheels only through root setup. When the root environment is
already synced, `bash tools/ensure-engine-wheel.sh` owns the canonical cached
wheel build and installation.

## Native dependency graph

The package owns one Meson graph and no installed cross-module SDK. Data
packet materializations are explicit and hash-bound. Installed-wheel analyses
live in `../analyses/`. See `../AGENTS.md`, `../ARCHITECTURE.yaml`, and
`../docs/scientific/CONTEXT.md`.

Meson discovers threads natively; Ipopt 3.11.9 or newer and Eigen 3.4 or newer
through pkg-config; Ceres 2.2 through its CMake package metadata. Exact
archive wraps own CppAD 20250000.3 and HiGHS 1.15.1. CMake only resolves the
Ceres imported target; the Engine has no CMake build route.

The HiGHS tag archive replaces its same-version source release asset because
that asset embeds Git metadata. All archive URLs and SHA-256 values are pinned.
Downloads and extracted trees stay ignored under `engine/build/subprojects/`
and never enter the wheel.

## Scientific specification and availability

The [scientific owner map](../docs/scientific/CONTEXT.md) identifies the
permanent equations and methods. The
[GitHub issue roadmap](https://github.com/tannerpolley/ePC-SAFT/milestones)
records implementation gaps. Do not infer
availability from historical usage examples or a registered equation.

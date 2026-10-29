# Unified ePC-SAFT Engine

One installed Python distribution owns the Engine. Current state/equilibrium imports:

```python
import epcsaft
from epcsaft import equilibrium
```

Consumers import `epcsaft` from their environment. They do not run Meson or
load modules from this source tree.

## Build and installation

The root [`README.md`](../README.md) "Setup and installation" owns setup (an
editable install into the checkout's `.venv`) and every check.

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
native shared library. It includes only `surface/parameter_builder.hpp` and
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
tools/native-run surface_example   # prints every value with its unit
```

For a debugger, configure an unstripped debug tree once from the repository root
(the first build also compiles HiGHS, about 2 minutes) and run the same target:

```bash
export PATH="$PWD/.venv/bin:$PATH"   # Meson and the nanobind it finds live in the .venv
meson setup build/engine-debug engine -Dbuildtype=debug -Dstrip=false
meson compile -C build/engine-debug surface_example
gdb --args build/engine-debug/surface_example
```

## Python counterpart

### First calculations

Run from the repository root. Units are SI molar: T in K, P in Pa, density in
mol/m^3, energies in J/mol, molality in mol per kg of solvent. Positional
vectors (`x=`, `composition_guess=`, result lists) follow `component_ids`, which
is the `components=` order; feeds, `x=` and `composition_guess=` also take a
name mapping (unlisted names are zero, unknown names are refused).
`Parameters.molar_masses` holds kg/mol in that order (0 where the document
declares none, as for most ions). A result list with one
entry per phase and component is phase-major: entry `p*n + i` is phase `p`
(declared order) and component `i`; `result.phase_mole_fractions(p)` returns
that phase's slice. [`data/docs/parameter-index.md`](../data/docs/parameter-index.md)
lists which parameter packets hold a component and what each set is for.
`Parameters.from_packets(root, [(packet, components), ...], k_ij={(a, b): value})`
combines components from several packets and sets a constant k_ij or l_ij by
component pair, replacing any stored value of that pair;
`epcsaft.parameters.parameter_document` returns the combined document with each
packet's sources. Packets whose model families or `model_coefficients` differ
are refused by name, as are a cross-packet pair without k_ij and two associating
components from different packets without a declared cross association.

```python
import math
import epcsaft
from epcsaft import equilibrium as eq, regression as rg
pk = "data/packets/{}/parameters/parameters.json".format

water = epcsaft.Pure(epcsaft.Parameters.from_json(pk("held-cameretti-sadowski-2008/4"), components=["water"]))
print(water.state(298.15, P=1.0e5, phase="liquid").molar_density * water.parameters.molar_masses[0])  # kg/m^3
sat = eq.solve_equilibrium(water, water.saturation(373.15))        # phases (L, V) at fixed T, free P
print(sat.pressure, sat.molar_densities)                        # Pa; mol/m^3 per phase

gs = epcsaft.Parameters.from_json(pk("gross-2001-methane-ethane/3"), components=["methane", "ethane"])
mix, feed = epcsaft.Mixture(gs), eq.MoleFractions({"methane": 0.3, "ethane": 0.7})
bub = eq.solve_equilibrium(mix, eq.bubble_point(liquid=eq.Phase("L", kind="liquid"), T=200.0, feed=feed))
print(bub.pressure, bub.phase_mole_fractions(1))                # Pa; [y_CH4, y_C2H6]

brine = epcsaft.Mixture(epcsaft.Parameters.from_json(pk("held-cameretti-sadowski-2008/4"),
                                                     components=["water", "sodium-cation", "chloride-anion"]))
nacl = brine.reference(298.15, 1.0e5, solvent_fractions=(1.0, 0.0, 0.0),  # one entry per component
                       ions=("sodium-cation", "chloride-anion"), stoichiometry=(1, 1)).activity(1.0)  # mol/kg
print(nacl.mean_ionic_molality_scale, nacl.osmotic_coefficient)    # gamma+- (molality scale), phi
one_molal = eq.Molalities({"sodium-cation": 1.0, "chloride-anion": 1.0}, solvent=["water"], solvent_mass_kg=1.0)
sol = eq.solve_equilibrium(brine, eq.Problem(phases=[eq.Phase("L", kind="liquid")], T=298.15, P=1.0e5, feed=one_molal))
print(sol.mole_fractions, sol.amounts)                             # amounts in mol for 1 kg of water

kh = epcsaft.Mixture(epcsaft.Parameters.from_json(pk("khudaida-2026-figure-2-electrolyte-lle/7"),
                                                  components=["water", "isobutanol"]))
lle = eq.solve_equilibrium(kh, eq.Problem(phases=eq.Discover(), T=298.15, P=1.0e5,
                                          feed=eq.MoleFractions({"water": 0.7, "isobutanol": 0.3})))
print(lle.phase_count, lle.mole_fractions, lle.amounts)            # phases: lle.discovery.declaration.phases

what_if = epcsaft.Parameters.from_packets("data/packets", [("gross-2001-methane-ethane/3", ["methane", "ethane"])],
                                          k_ij={("methane", "ethane"): 0.02})  # a what-if k_ij; here it makes test data
truth, obs = epcsaft.Mixture(what_if), []
for T, x in [(190.0, 0.2), (200.0, 0.4), (220.0, 0.3)]:            # K, liquid methane mole fraction
    bubble = eq.bubble_point(liquid=eq.Phase("L", kind="liquid"), T=T, feed=eq.MoleFractions({"methane": x, "ethane": 1 - x}))
    measured = eq.solve_equilibrium(truth, bubble).pressure      # put a measured bubble pressure (Pa) here
    obs.append(rg.observation(gs, "pressure", bubble, observed=measured, form="log_ratio"))
kij = rg.coordinate(gs, [("k_ij", ("methane", "ethane"), 1.0)], origin=0.0, scale=0.01)
fit = rg.fit(gs, [kij], obs, controls=rg.FitControls(maximum_iterations=50, maximum_elapsed_time_seconds=60.0))
print(fit.status, fit.physical)                                    # fitted k_ij; rg.fitted_parameters(gs, fit)

def envelope_at(ln_p):                                             # bubble-side envelope point at ln(P/Pa)
    target = eq.EnvelopeCondition(); target.kind, target.target = eq.EnvelopeCoordinate.LogPressure, ln_p
    return eq.Problem(phases=[eq.Phase("L", kind="liquid"), eq.Phase("V", kind="vapor", amount=eq.Pinned(0.0))],
                      T=eq.Free(180.0), P=eq.Free(1.0e6), feed=feed, envelope=target)
steps = eq.Continuation(initial_step=0.05, minimum_step=1e-4, max_steps=2000, max_retries_per_step=5,
                        envelope_max_log_step=0.1, envelope_straddle=0.02)
trace = eq.continue_equilibrium(mix, envelope_at(math.log(1e6)), envelope_at(math.log(8e6)), steps)
print(trace.stop, len(trace.states), [(c.temperature, c.pressure) for c in trace.critical_points])
```

- LLE: `Discover()` finds neutral splits without guesses. With ions, declare
  two liquid phases with `composition_guess` vectors (one entry per component).
- Envelope mode needs `envelope_straddle > 0` and `2 * envelope_straddle <=
  envelope_max_log_step`, and T and P guesses through `Free(guess)`. It uses
  neither `initial_step` nor the target problem: a mixture trace closes at the
  start pressure on the other branch. Each ordinary step, in log units, is at
  most `envelope_max_log_step` (the predictor can take less); a rejected
  ordinary step halves, and the cap returns after an accepted step. A rejected
  critical-point crossing, closing step or terminal failure ends the trace at
  once, keeping the states accepted so far and that rejection's stop and
  message. `trace.path_values[k]` is the target of `trace.conditions[k]` (ln T,
  ln P or ln K of `component`); the traced coordinate can switch along the way.
- Ordinary continuation (no `envelope`) steps are fractions of the
  start-to-target path. A rejected step halves and never grows again, so
  `initial_step` is the largest step.
- In both modes, more than `max_retries_per_step` consecutive rejections of
  retryable steps stops with `budget_exhausted`, and a halved step below
  `minimum_step` stops with `step_too_small`; `cause_stop` and `cause_message`
  give the last rejected step's reason.

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
from epcsaft import IdealCorrelation, IdealInterval, IdealPolynomial, Mixture, Parameters, ThermochemistryRecord
from epcsaft.equilibrium import EOS_STANDARD_STATE_ID, Reaction, ReactionLogPolynomial

parameters = Parameters.from_json("data/packets/khudaida-2026-figure-2-electrolyte-lle/7/parameters/parameters.json",
                                  components=["water", "isobutanol"])
cp = IdealCorrelation(IdealPolynomial([80.0, 0.0], 300.0, 1000.0, 50.0), source_pressure=1.0e5)  # cp = 80 J/(mol K)
record = ThermochemistryRecord(1.0e5, [[IdealInterval(250.0, 450.0, True, True, cp)]] * 2)  # one list per component
model = Mixture(parameters, thermochemistry=record)           # omit the record for residual-only States
state = model.state(353.15, P=1.0e5, x=(0.25, 0.75), phase="liquid")
print(state.total_enthalpy, state.ideal_entropy)              # None when the record does not cover T
ln_k = ReactionLogPolynomial(2.0, -1000.0, 0.0, 0.0, 298.15, 273.15, 373.15, EOS_STANDARD_STATE_ID)  # a, b, c, d, T_ref, T_min, T_max
reaction = Reaction({"A": -1.0, "B": 1.0}, correlation=ln_k)  # exactly one of K, affinity or ln K(T)
```

A NASA9 record file (`nasa-glenn-nasa9-records` version 1, coefficient order a1..a7, b1, b2) loads once:
`packet = ThermochemistryPacket(path)` rescales all nine coefficients of every interval by
R_source/R_Engine in native code, so cp, h and s keep their source values, and refuses a malformed form or
version, duplicate species and an invalid gas constant, pressure, coefficient or interval chain. Intervals are
[lower, upper), the last closed. `packet.select(component_ids)` returns a `ThermochemistryRecord` in the
selection's order; it copies the prepared intervals and never rescales again.

Reaction sets come from the caller. `ReactionLogPolynomial` is
ln K(T) = a + b/T + c ln(T/T_ref) + d T + e/T² + f T², with trailing optional
`e=0.0` (K²) and `f=0.0` (K⁻²), on the Engine standard state
(`EOS_STANDARD_STATE_ID`); a K or correlation on a source basis (molality or
mole fraction at infinite dilution) also needs `Reaction(reference=ReactionReference(...))`
and the problem's `neutral_reference`. Data packets carry no complete reaction
set: `mea-co2-h2o-nine-species-estimation/5` holds candidate correlations for
R4 and R5 only. `Parameters.from_json` refuses a top-level field outside a
named set and names each such field. The set holds the fields it reads, record
metadata it admits unread (`document_id`, `document_version`, `purpose`,
`sources`, `domains`, `materialization`) and two caller-owned records it admits
unread (`reaction_correlations`, `empirical_density_correction`). Adopting reaction sets is
application work: the adopted MEA set (R1-R5, their standard-state offsets and
the Born constants `c_shell = c_dielectric = 1.0`, where the packet declares
0.0) lives in Amine-Thermodynamics, `src/MEA/common/mea_source_contracts.py`
and `analyses/mea_parameter_bundle/scripts/shared_evaluation.py`.

Actions are order one to three contractions along caller directions and keep the
observable's unit. A state direction holds K, mol/m^3 or Pa, a zero-sum
composition increment and native family/slot parameter increments; a solved-state
direction holds K, Pa, absolute feed mol (component width) and one increment per
active parameter (active width). A free T or P is an output: its entry is zero.
`regression.active_parameter(parameters, family, target)` names an active parameter
with the spelling of one `regression.coordinate` target; `ActiveParameter(family, indices)`
takes native flat slots. In `solved_state_actions`, k directions (k ≤ 3) give one mixed
derivative of order k, D^k y[d1, …, dk]: two directions give ∂²y/∂a∂b, not a gradient.
`solved_state_first_actions` returns one first-order batch per direction, sharing one
center and one factorization; that is the gradient.

```python
import math
from epcsaft import (GAS_CONSTANT_J_PER_MOL_K, DensityActionDirection, Mixture, ParameterActionDirection, Parameters,
                     PropertyObservable)
from epcsaft.equilibrium import (MoleFractions, Phase, Problem, SolvedStateActionDirection, SolvedStateActionRequest,
                                 SolvedStateObservable, SolvedStateObservableKind as Kind, compile_problem,
                                 solve_equilibrium, solved_state_actions, solved_state_first_actions)
from epcsaft.regression import active_parameter

parameters = Parameters.from_json("data/packets/khudaida-2026-figure-2-electrolyte-lle/7/parameters/parameters.json",
                                  components=["water", "isobutanol"])  # sites: isobutanol a, b, water a, b
pair = active_parameter(parameters, "association_energy_over_k", (("isobutanol", "a"), ("isobutanol", "b")))  # slots (1, 4)
model = Mixture(parameters, active_parameters=[pair])
T, rho, x = 298.15, 12000.0, (0.3, 0.7)                        # K, mol/m^3, mole fractions
move = DensityActionDirection(1.0, 0.0, (0.0, 0.0), [ParameterActionDirection(pair.family, pair.indices, 0.03)])
batch = model.density_actions(T, rho, x=x, observables=[PropertyObservable.ResidualEnthalpy], directions=[move])
result = batch.results[0]                                     # J/mol; action along move, None beside a typed status
print(result.value, result.action, result.status)

problem = Problem(phases=[Phase("L1", kind="liquid", composition_guess=(0.46, 0.54)),
                          Phase("L2", kind="liquid", composition_guess=(0.99, 0.01))],
                  T=T, P=1.0e5, feed=MoleFractions({"water": 0.7, "isobutanol": 0.3}))
central = solve_equilibrium(model, problem)
compiled = compile_problem(model, problem)                    # one compile; no solve is repeated
warm = SolvedStateActionDirection(1.0, 0.0, [0.0, 0.0], [0.0])  # +1 K; T, P, feed mol, active parameters
activity = [SolvedStateObservable(Kind.PhaseLogActivity, 0, 1)]
request = SolvedStateActionRequest(activity, [warm])
row = solved_state_actions(compiled, central, request).results[0]
print(central.mole_fractions, row.value, row.action)         # ln a of isobutanol in L1; d ln a/dT in 1/K
# Caller-side conversion of the native activity, not an Engine calculation: ln a_i is on the EOS density
# basis, ln(x_i rho / (1 mol/m^3)) + mu_i^res/RT, so f_i = a_i R T (1 mol/m^3) in Pa.
fugacity = math.exp(row.value) * GAS_CONSTANT_J_PER_MOL_K * central.temperature
print(fugacity)
bond = SolvedStateActionDirection(0.0, 0.0, [0.0, 0.0], [1.0])  # +1 K in the active epsilon_AB/k
gradient = solved_state_first_actions(compiled, central, activity, [warm, bond])  # one batch per direction
print([batch.results[0].action for batch in gradient])        # d ln a/dT, d ln a/d(epsilon_AB/k), both 1/K
```

`solve_equilibrium` has no time limit; `solve_equilibrium_starts` tries declared
starts in order under shared work limits and stops at the first that succeeds with its requested tolerance met:

```python
from epcsaft import Mixture, Parameters, equilibrium

model = Mixture(Parameters.from_json("data/packets/khudaida-2026-figure-2-electrolyte-lle/7/parameters/parameters.json",
                                     components=["water", "isobutanol"]))
problem = equilibrium.Problem(phases=[equilibrium.Phase("L1", kind="liquid", composition_guess=(0.46, 0.54)),
                                      equilibrium.Phase("L2", kind="liquid", composition_guess=(0.99, 0.01))], T=298.15, P=1.0e5,
                              feed=equilibrium.MoleFractions({"water": 0.7, "isobutanol": 0.3}))
limits = equilibrium.EquilibriumWorkLimits(seconds=60.0, attempt_seconds=20.0)
tried = equilibrium.solve_equilibrium_starts(model, problem, [equilibrium.EquilibriumStart()], limits)
result = None if tried.selected is None else tried.results[tried.selected]  # tried.causes says why not
print(tried.selected, tried.causes, result.mole_fractions, tried.iterations)
```

A default `EquilibriumStart` runs the cold start of `solve_equilibrium`. If the limits are not exhausted, it
reaches the same state with the same work (#350), except where the plain solve's speciated seeded attempt fails:
`solve_equilibrium` then retries from its own start and the bounded start does not. A partial start generates only what it omits. Seed solves,
such as a reacting bubble's liquid pre-solve, spend the same limits, and both functions include their work in
the reported counts.

Invalid declarations and inputs raise `ParameterError` (a `ValueError`) before iterating. A
solve that runs reports non-convergence, structural infeasibility and rejected density
branches as `success=False` with `message` and `convergence.stop`, but a failure to
evaluate the returned point still raises `EosError` or `ParameterError`; handle both.
`solve_equilibrium_starts` records an attempt's exception in `causes` instead of raising.

For a fixed-T/P, declared two-free-liquid equilibrium, `equilibrium.asl_record`
constructs the immutable record consumed by the installed
`epcsaft_equilibrium_phase_component_amount` external function:

```python
from epcsaft import Mixture, Parameters, equilibrium as eq

model = Mixture(Parameters.from_json("data/packets/khudaida-2026-figure-2-electrolyte-lle/7/parameters/parameters.json",
                                     components=["water", "isobutanol"]))
both = ("water", "isobutanol")                                # declared supports; an Amounts feed in mol
problem = eq.Problem(phases=[eq.Phase("L1", kind="liquid", support=both, composition_guess=(0.46, 0.54)),
                             eq.Phase("L2", kind="liquid", support=both, composition_guess=(0.99, 0.01))],
                     T=298.15, P=1.0e5, feed=eq.Amounts({"water": 0.7, "isobutanol": 0.3}))
continuation = eq.Continuation(initial_step=0.25, minimum_step=0.25, max_steps=4, max_retries_per_step=1)
feed_offset, feed_matrix = [0.0, 0.0], [[1.0, 0.0], [0.0, 1.0]]  # one control per component
anchor_controls, identity_component, minimum_gap = [0.7, 0.3], 1, 0.1  # isobutanol: x_L1 - x_L2 > 0.1
record = eq.asl_record(model, problem, continuation, feed_offset,
                       feed_matrix, anchor_controls, identity_component, minimum_gap)
# Pyomo ExternalFunction: (*controls, record, "0" or "1", str(component_index))
print(len(record))                                            # characters in the encoded record
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
Each new feed is first solved directly from the record's declared guesses and
kept only if the record's two-phase, tolerance and identity-gap checks admit it;
otherwise it continues from the immutable anchor, whose stop is the reported
failure. An admitted state need not be reachable by a straight path from the
anchor. Each distinct record, by
exact content, is prepared once per process and keeps its own model, anchor
and exact-iterate cache for the process lifetime; each preparation prints
`epcsaft_equilibrium prepared record <k>` on stderr. A failed
preparation retains no context. An evaluation failure invalidates every record's
cached response, and each record recomputes at its next call.
Failed central/branch/action checks return nonfinite values
and derivatives with an observable-specific typed error, never a surrogate
fallback. Use an outer Ipopt executable whose C++ symbols do not interpose on
the Engine's Ipopt dependency; the retained live campaign verifies the isolated
official IDAES executable and exact option-file readback.

`eos.property_asl_record(mixture)` constructs the immutable record of the installed
`epcsaft_pressure_property` external function (#73) for a mixture with document
parameters, no active parameters and a thermochemistry record in any native form
(polynomial, Shomate, NASA7, NASA9 or DIPPR107; record `epcsaft.asl.property.v2`).
`regression.fitted_parameters` writes the fitted values back into the start
document's entries, so both records export a fit; a fit the document cannot
represent (a fitted slot with no entry, a value that does not rebuild bitwise in
its entry's unit) stays usable and the record builders raise that cause. The function's real
arguments are `q = (T K, P Pa, w_1..w_C)` with `x = w / sum(w)`, then the strings
record, row (`0` molar density mol/m3, `1` total molar enthalpy J/mol, `2 + i`
`ln(phi_i)`) and phase (`0` liquid, the highest stable density root; `1` vapor,
the lowest). It returns the value, gradient and packed Hessian in physical `q`
from the native pressure actions plus the exact normalization chain. A zero
neutral weight follows the exact infinite-dilution limit when the weights sum
to a positive value. A sole stable root is returned for either phase; the
public `pressure_actions` diagnostics report `branch` and `stable_root_count`
for the caller's post-solve audit. Each distinct record is prepared once per process
and kept for its lifetime, as above (`epcsaft_property prepared record <k>`). A failure returns
nonfinite values with `epcsaft_property:<cause>:<PropertyActionStatus code>:...`,
the cause one of `invalid_input`, `support_insertion`, `no_density_root`,
`ambiguous_branch`, `nonmechanical_root`, `inner_iteration_failure` or
`unavailable_derivative`; an unavailable row leaves the other rows valid, and an
owned thermochemistry endpoint keeps its value while its temperature derivative is
unavailable. Derivatives are qualified for neutral records only: a single-ion
weight direction is not charge balanced, and native admission refuses it.

Diagnostics keep their native meaning. `State.thermochemistry_status` separates a
missing record from a temperature outside its intervals; `density_branch` and
`stable_root_count` exist only after a pressure-root search, so a T-rho state reports
`Unspecified` and `None`.

`State.requested_density_branch` retains the branch passed to pressure closure and is
`Unspecified` for a direct-density State. `returned_density_branch` reads ALG-DENSITY's
existing finite-search evidence or locally admitted anchor label and is `Unspecified`
when no side is established. `density_branch_mismatch` is true only for explicit,
opposite requested and returned branches. False means no established mismatch; it does
not prove branch agreement or physical stability. A wrongly labelled external anchor
can conceal a mismatch. These attributes are read-only in Python.

`EquilibriumResult.multipliers` are the dimensionless
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

Root setup installs this package editable; `tools/build-engine.sh` owns every
Engine build, including the wheels for application pins and A/B timing.

## Native dependency graph

The package owns one Meson graph and no installed cross-module SDK. Data
packet materializations are explicit and hash-bound. Installed-wheel analyses
live in `../analyses/`. See `../AGENTS.md`, `../ARCHITECTURE.yaml`, and
`../docs/scientific/CONTEXT.md`.

Meson discovers threads natively; Ipopt exactly 3.11.9 and Eigen >=3.4,<3.5
through pkg-config; Ceres 2.2 through its CMake package metadata. Exact
archive wraps own CppAD 20250000.3 and HiGHS 1.15.1. CMake only resolves the
Ceres imported target; the Engine has no CMake build route.

The HiGHS tag archive replaces its same-version source release asset because
that asset embeds Git metadata. All archive URLs and SHA-256 values are pinned.
Downloads and extracted trees stay ignored under `engine/build/subprojects/`
and never enter the wheel.

## Conda package

`../packaging/conda-forge/` prepares Linux x86-64 / CPython 3.13 with private
Ipopt 3.11.9, MUMPS 5.8.2/OpenBLAS, Ceres 2.2.0 CPU, Eigen 3.4.0 and
nanobind 2.12.0. It installs the Python package, extension, native and ASL
libraries, private Ipopt and notices; parameter packets, analyses, manuscript
sources and Coin-HSL are excluded. Supply parameter declarations explicitly.
The release version and public archive checksum remain publication placeholders.

Corresponding source is available from the distributor as a source-set download
alongside the matching version at
https://github.com/tannerpolley/epcsaft/releases: exact Engine and Ipopt sources,
embedded third-party sources, recipe/build instructions, patches and original
notices with SHA-256 checksums. Before binary publication, verify retrieval of
those exact bytes. See `../NOTICE` for the source-availability rule and terms.

## Scientific specification and availability

The [scientific owner map](../docs/scientific/CONTEXT.md) identifies the
permanent equations and methods. The
[GitHub issue roadmap](https://github.com/tannerpolley/ePC-SAFT/milestones)
records implementation gaps. Do not infer
availability from historical usage examples or a registered equation.

## Solver stopping tolerances

`Mixture` and `Pure` accept `solver_options={name: value}` for association,
density and equilibrium rows. `Problem` accepts the row option as a declaration
override; `row_tolerance` writes that same slot. `FitControls` accepts all six
options, and its existing `function_tolerance`, `gradient_tolerance` and
`parameter_tolerance` arguments write the corresponding slots (zero means
omitted). Declaring the same slot by both spellings is an error. Point readers
(`differentiate` with parameters, `evaluate_coordinates`, `evaluate_checks`, and
the native point readers) accept the association, density and row names.

| Name | Default | Unit and stopping meaning |
|---|---:|---|
| `ASSOCIATION_RESIDUAL` | `2e-13` | Dimensionless mass-action infinity norm. |
| `DENSITY_ROOT_RESIDUAL` | `1e-12` | Dimensionless reduced pressure residual `(P_calc-P)/(rho R T)`, with a bracket-width check. |
| `EQUILIBRIUM_ROW_RESIDUAL` | `1e-10` | Dimensionless classified equilibrium-row infinity norm. |
| `FIT_FUNCTION_TOLERANCE` | `1e-12` | Relative least-squares cost decrease. |
| `FIT_GRADIENT_TOLERANCE` | `1e-12` | Projected gradient infinity norm, scaled objective per dimensionless optimizer coordinate. |
| `FIT_PARAMETER_TOLERANCE` | `1e-10` | Relative optimizer-coordinate step. |

For example, `Pure(parameters, solver_options={"DENSITY_ROOT_RESIDUAL": 1e-8})`
requests less density refinement. Model settings are immutable. A declaration
inherits its model's row setting; an explicit default value resets an inherited
custom value. Fit settings reach every constructed model and worker process;
individual observation declarations can still override their own row target.
Changing association stopping does not change the independent density default.

Values must be finite positive numbers, at most the largest representable
`double`. Unknown or fixed policies, names outside the owning input's scope,
booleans and nonnumeric values are rejected. A derived Ipopt target that
underflows or overflows is also rejected. There is no universal observable-error
bound or guaranteed speedup: reduced pressure error can amplify near a spinodal,
and loose central equilibrium states can fail the fixed derivative-admission
check. Search geometry, physical admission and derivative checks stay fixed.

`epcsaft.tolerances()` reports immutable source defaults. Each owning result's
`solver_options` lists effective names, values, units and defaults. Equilibrium
results also report consumed Ipopt `solver_settings` (empty for Newton-only
solves) and the diagnostic options-file path when supplied. Fit results' common
options do not replace heterogeneous row overrides retained by observations.
Achieved residuals, numerical floors and conditioning remain the evidence for
assessing a returned value; requested tolerances alone are not error estimates.

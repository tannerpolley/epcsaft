# Stage 4 — `bindings/` and the Python EOS surface

You are implementing Stage 4 of a greenfield rebuild in this worktree (branch `stage-3`).
You have no one to ask. Read, in order, before writing anything:

1. `engine/docs/greenfield/TESTING.md` — a hard restriction on what may be tested
2. `engine/docs/greenfield/TRANSFER.md` §"Stage 4", §"Units: C++ resolves them, Python only
   reads the file", §"The MEA parameter set is the living bundle", §"Naming"
3. `engine/docs/greenfield/STAGE-LOG.md` — Stages 0–2 decisions; you build on them
4. `engine/native/model/helmholtz_model.hpp` completely — it is the surface you bind
5. `engine/native/physics/parameters.hpp` — the arrays your builder fills
6. `engine/tests/native/model_checks.cpp` §`density_and_reference` — the native reference call
7. `tools/check-engine-architecture.py` — seven rules; rule 6 (no threshold literals) applies to
   every file you write under `engine/native/`
8. `engine/docs/greenfield/briefs/stage-3.brief.json` — sources, budgets, gates, mutations

Verify every source marker and every precondition in the JSON against the tree **before**
writing. If a claim is false, stop and say which.

## 0. Non-negotiable

**Layers.** `physics <- model <- equilibrium <- regression <- bindings`. `native/bindings/` is
the only directory that may include pybind11; nothing includes `bindings/`. The parameter
builder lives in `native/model/` and includes `"physics/..."` and `"model/..."` only, rooted at
`engine/native`, never relative.

**Python does no arithmetic.** Any Python method takes the inputs it is given, hands them to
`_core`, and returns what comes back. No `math`, `numpy`, `pint`, `scipy`, `decimal`; no unit
factor, no `exp`, no normalisation, no combining rule, no index arithmetic on physical
quantities. Python walks JSON and resolves *names* to *indices*; every number crosses the
boundary untouched, as `(identity, field, index, magnitude, unit_string)`. The gate
`python_does_no_arithmetic` greps for this and its mutation is `import math`.

**Units resolve natively, once.** The builder in `model/parameter_builder.cpp` owns the unit
table (§2.1). A unit it does not accept for a field is rejected by the record's identity and
the unit string. Pint leaves `engine/pyproject.toml`. `units.py` is deleted, not rewritten:
with Pint gone there is nothing in it.

**Bindings are `py::class_`**, never `py::dict` round trips. Native result structs are bound
with `def_readonly` and returned as they are; Python does not copy them into dataclasses.

**Constants and tolerances.** No numeric literal acting as a threshold in `engine/native/`.
Rule 6 catches a literal adjacent to a comparison operator; conversion factors in the unit
table (`1.0e-3`) are not thresholds and are permitted. Nothing in this stage should need a new
`numerical_policy.hpp` entry; if you believe one is needed, add it in the existing style and
report it.

**Naming.** No proper nouns in identifiers. Sources are cited in comments and docstrings.

**Testing.** Before writing a check ask: would it survive a legitimate reimplementation of the
module from its equations? Only real measured data, exact symmetries, closed-form identities,
and fail-closed validation qualify. No golden files, no assertions on internals. Counts and
timings are reported, never pass/fail.

**Environment is not yours to repair.** You may edit `engine/meson.build` only to add
`native/model/parameter_builder.cpp` to `model_sources` and to replace the empty `files(),`
of the `_core` extension with your two binding sources. You may **not** touch its dependency
block, `default_options`, the virtual environment, `uv.lock`, `engine/build/subprojects`, or
how you were invoked, in order to make a build succeed. `pip`, `uv pip`, `uv sync`, and
`uv lock` are not yours to run. Stop with `environment_broken`.

**Three environment facts you must meet before the first build** (§6): `PATH` leads with the
worktree's `.venv/bin`; `$HOME/.cache` is writable; the build directory's `prefix` is the
worktree's `.venv`.

**Model and physics are closed.** Five edits permitted in `native/model/` and none in
`native/physics/`: (a) the new `parameter_builder.{hpp,cpp}`; (b) one field
`double mean_ionic_molality_scale = 0.0;` added to `ElectrolyteActivity` in
`helmholtz_model.hpp`; (c) one assignment `result.mean_ionic_molality_scale =
std::exp(mean_log_m);` beside the existing `result.mean_ionic_log_molality_scale = mean_log_m;`
in `infinite_dilution_reference.cpp`; (d) in `helmholtz_model.cpp`, the line
`require(positive_vector(parameters_.molar_mass_kg_per_mol), "molar_mass");` becomes

```cpp
    for (std::size_t component = 0; component < components; ++component) {
        const double mass = parameters_.molar_mass_kg_per_mol[component];
        require(std::isfinite(mass) && mass >= 0.0, "molar_mass");
        require(parameters_.charge_number[component] != 0 || mass > 0.0, "molar_mass");
    }
```

(a charged component may carry 0.0, meaning not provided; a neutral one may not); (e) in
`formula_composition`, the condition `|| molar_masses[component] <= 0.0` becomes
`|| (solvent_fractions[component] > 0.0 && molar_masses[component] <= 0.0)` — a mass is
required where a solvent fraction reads it, and nowhere else. Python must not `exp()` a
log; the native side hands over the number the oracle is stated in.

**Why (d) and (e).** Molar mass never enters the Helmholtz energy. Every native consumer
reads solvent masses only: the three permittivity-mixing sites in `physics/electrolyte.hpp`
are all inside `dielectric_solvent_component()`, which returns false for any charged
component, and both `infinite_dilution_reference.cpp` sites multiply by a solvent fraction
the reference composition forces to zero for the ions. Requiring a charged component's mass
is a false rejection of sufficient data, and it is what stopped this stage's first attempt.
These are relaxations of Stage 3's validation, not new arithmetic; every existing fixture
carries positive masses, so the Stage 3 checks stay green under them.

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

### 2.1 The parameter builder — `model/parameter_builder.{hpp,cpp}`

```cpp
namespace epcsaft {
struct ParameterEntry final {
    std::string identity;            // the JSON identity, for error messages only
    std::string field;               // a name from the table below
    std::vector<std::size_t> index;  // positions, resolved by Python from ids
    double magnitude;
    std::string unit;
};
struct ParameterDocument final {
    std::vector<std::string> component_ids;   // selection order = model order
    std::vector<int> site_component;          // S entries
    std::vector<int> site_multiplicity;       // S entries
    std::string electrolyte_choice;           // "none" | "fully-dissociated-debye-huckel" | "born"
    std::string permittivity_choice;          // "none" | "solvent-only" | "ion-fraction-suppression"
                                              // | "ion-specific-suppression" | "component-permittivity-mixing"
    std::vector<ParameterEntry> entries;
};
Parameters build_parameters(const ParameterDocument&);   // throws std::invalid_argument
}
```

The builder allocates every `Parameters` array at its size (`C`, `C*C`, `S`, `S*S`, terms
appended), fills from entries, then checks completeness. Every rejection is
`std::invalid_argument` whose message names the identity (or `component_id` / pair) and what
was wrong: `"unit gram for component/water/molar_mass: expected kilogram / mole"`,
`"missing pair/ethanol/methanol/k_ij"`, `"duplicate component/water/segment_diameter"`.
`HelmholtzModel`'s constructor validates the rest; do not duplicate its checks.

**The unit table.** Nine strings are the whole parameter vocabulary in `data/packets/`;
pressure units appear only in domains and observations and never reach a parameter.

| dimension | accepted strings → factor |
|---|---|
| dimensionless | `dimensionless` → 1 |
| temperature / energy over k | `kelvin` → 1 |
| length | `angstrom` → 1 (native arrays are in Å) |
| charge | `elementary-charge` → 1, magnitude must be integral |
| molar mass | `kilogram / mole` → 1, `kg / mol` → 1, `gram / mole` → 1e-3 |
| reciprocal temperature | `1 / kelvin` → 1 |
| reciprocal temperature to power p | `1 / kelvin ** p` for the recorded `power` p (p = 1 is spelled `1 / kelvin`) → 1 |

The molar-mass row must be written literally as `{"gram / mole", 1.0e-3}` — a gate mutation
targets that text.

**Fields.** `i, j` are component positions, `s, t` site positions, `n` a term position
assigned by Python in document order.

| `field` | `index` | unit | fills | required |
|---|---|---|---|---|
| `molar_mass` | `[i]` | molar mass | `molar_mass_kg_per_mol[i]` | every component with `charge_number == 0`; a charged component may omit it and the builder stores 0.0 — molar mass never enters the Helmholtz energy and every native consumer (solvent permittivity mixing, molality) reads solvent masses only; an operation that needs a charged component's mass rejects by name at its own boundary |
| `charge_number` | `[i]` | charge | `charge_number[i]` | every component |
| `segment_count` | `[i]` | dimensionless | `segment_count[i]` | every component |
| `sigma_constant` | `[i]` | length | `sigma_constant[i]` | every component, exactly once (a `segment_diameter` coefficient, or the `constant` of a `segment_diameter` correlation) |
| `sigma_term_amplitude` | `[n, i]` | length | `sigma_term_amplitude[n]`, `sigma_term_component[n] = i` | per exponential term |
| `sigma_term_exponent` | `[n, i]` | reciprocal temperature | `sigma_term_exponent[n]` | per exponential term |
| `dispersion_energy_over_k` | `[i]` | kelvin | `dispersion_energy_over_k[i]` | every component |
| `solvation_factor` | `[i]` | dimensionless | `solvation_factor[i]` | every component under `born`; else absent → 0 |
| `relative_permittivity` | `[i]` | dimensionless | `relative_permittivity[i]` | absent → 0 (non-dielectric) |
| `permittivity_reference_temperature` | `[i]` | kelvin | `permittivity_reference_temperature_k[i]` | with a polynomial correlation |
| `permittivity_term_coefficient` | `[n, i, p]` | reciprocal temperature to power p | `permittivity_term_coefficient[n]`, `permittivity_term_power[n] = p`, `permittivity_term_component[n] = i` | per polynomial term |
| `packing_diameter`, `debye_huckel_diameter` | `[i]` | length | the two diameter arrays | every charged component; absent → 0 |
| `born_diameter` | `[i]` | length | `born_diameter_angstrom[i]` | every charged component under `born` |
| `ion_suppression_coefficient` | `[i]` | dimensionless | `ion_suppression_coefficient[i]` | see permittivity choices |
| `k_ij` | `[i, j]` | dimensionless | `k_ij[i*C+j]` **and** `k_ij[j*C+i]` | every unordered pair of selected components |
| `l_ij` | `[i, j]` | dimensionless | both orders | absent → 0 |
| `k_ij_reciprocal_temperature_slope` | `[i, j]` | kelvin | both orders | absent → 0 |
| `k_ij_reference_temperature` | `[i, j]` | kelvin | `k_ij_reference_temperature_k`, both orders | with a slope; absent → 0 |
| `association_energy_over_k` | `[s, t]` | kelvin | `[s*S+t]` and `[t*S+s]`; sets `association_active` to 1 at both | per explicit edge |
| `association_volume` | `[s, t]` | dimensionless | both orders | per explicit edge |
| `association_combining_rule` | `[s, t, s1, t1, s2, t2]` | (magnitude and unit ignored) | see below | per combining-rule edge |
| `ionic_region_relative_permittivity` | `[]` | dimensionless | the scalar | under `born` |
| `born_c_shell`, `born_c_dielectric` | `[]` | dimensionless | the two scalars | under `born` |

**Unconsumed records are rejected.** A record whose field the declared choices do not read
is a probable transcription loss and is rejected by identity: `ion_suppression_coefficient`
(either origin) under `none` or `solvent-only`; `ionic_region_relative_permittivity`,
`born_c_shell`, `born_c_dielectric`, `born_diameter`, `solvation_factor` under an
electrolyte choice other than `born`. A present coefficient whose value is zero is data —
packets retain printed zeros — and is never rejected for being zero.

Symmetric fills must be written as two assignments, e.g. `k_ij[i * components + j] = value;`
and `k_ij[j * components + i] = value;`. A field given twice at the same index is a
`duplicate`.

**The combining rule.** `association_combining_rule` at `[s, t, s1, t1, s2, t2]` names the
cross edge `(s, t)` and the two explicit self edges `(s1, t1)` and `(s2, t2)` it derives from.
After all explicit edges are filled:

    energy(s, t) = 0.5 * (energy(s1, t1) + energy(s2, t2))
    volume(s, t) = sqrt(volume(s1, t1) * volume(s2, t2))

both orders, `association_active` 1. This is the equation the superseded engine evaluated at
run time (`_dictionary_io.py`: `volume = (volume_values[0] * volume_values[1]).sqrt()`) and
the one the Stage 2 nine-species fixture encodes: `0.0411037990944876 = sqrt(0.03747 *
0.04509)`. It has **no** diameter factor. (The `figiel-yu-held-2025/1` CSV materialised its
solvent–solvent cross volumes with Gross–Sadowski 2002 Eq. 3's diameter factor instead; that is
a packet-revision divergence outside this stage — no gate here selects a solvent–solvent cross
edge. Do not add the factor.) Reject a rule whose self edges are not both active by the time it
is evaluated, naming the identity.

**Model choices.** `electrolyte_choice`: `none` → `ElectrolyteForm::None`;
`fully-dissociated-debye-huckel` → `DebyeHuckel`; `born` → `DebyeHuckelBorn`, which requires
`ionic_region_relative_permittivity`, `born_c_shell`, `born_c_dielectric`, every component's
`solvation_factor`, and every charged component's `born_diameter`. `permittivity_choice`:
`none` → `PermittivityForm::None`, all `ion_suppression_coefficient` 0; `solvent-only` →
`IonSuppression` with all 0; `ion-fraction-suppression` → `IonSuppression` with the one
`ion_suppression_coefficient` entry at `index []` — Python emits it once with an empty index
from `model_coefficients.ion_fraction_suppression_coefficient`, and the builder copies it to
every charged component and 0 to every neutral one; `ion-specific-suppression` →
`IonSuppression` with a per-component entry required for every charged component;
`component-permittivity-mixing` → `ComponentMixing`. Any other string → reject by name. The
five arrays the two unmapped forms would read (`dielectric_volume_m3_per_mol`,
`dielectric_temperature_k`, `dipole_scaling_debye2`, `polarizability_scaling_angstrom3`,
`correlation_integral_parameter`) are zeros of length `C`: no packet in the tree carries their
families, so nothing maps them.

### 2.2 The bindings — `bindings/module.cpp`, `bindings/eos_bindings.cpp`

`module.cpp` holds `PYBIND11_MODULE(_core, m)`, registers the two exceptions, and calls
`register_eos(m)` from `eos_bindings.cpp`. Use `<pybind11/stl.h>` for vectors, arrays,
optionals, and tuples. Exact names, because the stub and two mutations depend on them:

| Python name | binds | notes |
|---|---|---|
| `_core.ParameterError` | `std::invalid_argument` via `py::register_exception<std::invalid_argument>(m, "ParameterError", PyExc_ValueError)` | invalid document, unit, field, or call argument |
| `_core.EosError` | `std::runtime_error` via `py::register_exception<std::runtime_error>(m, "EosError", PyExc_RuntimeError)` | a solve failed: no or ambiguous density root, association, reference limit |
| `_core.ParameterEntry`, `_core.ParameterDocument` | the two structs, `py::init<>()` plus `def_readwrite` on every field | Python fills them |
| `_core.build_parameters(document)` | `build_parameters` | returns `_core.Parameters` |
| `_core.Parameters` | `epcsaft::Parameters` | opaque; `def_property_readonly("component_count", ...)`, `("site_count", ...)` |
| `_core.DensityBranch` | `py::enum_<DensityBranch>` with `Unspecified`, `Vapor`, `Liquid` | |
| `_core.State` | `epcsaft::State`, `def_readonly` on every field, same names | the line `.def_readonly("compressibility_factor", &epcsaft::State::compressibility_factor)` must appear verbatim — a mutation targets it |
| `_core.SweepCounters` | `def_readonly` on the four fields | |
| `_core.HelmholtzModel(parameters)` | `py::init<Parameters>()`, methods `state(T, molar_density, x)`, `state_at_pressure(T, P, x, branch, anchor_density: float | None)`, `counters()`, `reset_counters()` | `anchor_density` is `std::optional<double>`, default `None` |
| `_core.InfiniteDilutionReference`, `_core.ElectrolyteActivity` | `def_readonly` on every field, including the new `mean_ionic_molality_scale` | |
| `_core.infinite_dilution_reference(model, T, P, solvent_fractions, ion_indices, stoichiometry)` | the free function; `ion_indices` is `std::array<std::size_t, 2>`, `stoichiometry` `std::array<int, 2>` | |
| `_core.electrolyte_activity(model, T, P, solvent_fractions, ion_indices, stoichiometry, molality, reference)` | the free function | |

Do not bind `evaluate`, `solve_density`, `formula_composition`, active parameters, or the
tape. Do not release the GIL; the model's mutex is enough and nothing here threads.

### 2.3 Python — `src/epcsaft/`

Eight exports in `__init__.py`: `Parameters`, `ParameterError`, `EosError`, `Pure`,
`Mixture`, `State`, `Reference`, `Activity`. `State = _core.State`,
`Activity = _core.ElectrolyteActivity`, `ParameterError = _core.ParameterError`,
`EosError = _core.EosError` — re-exports, not wrappers. TRANSFER's "eleven" is a cap:
`Saturation` (Stage 5) and `ActiveParameter` (Stage 6) arrive with the stages that consume
them; `Derivatives` has no Python consumer and is not exposed.

**`parameters.py`.**

```python
class Parameters:
    component_ids: tuple[str, ...]           # selection order
    site_ids: tuple[tuple[str, str], ...]    # (component_id, site_id), document order
    @classmethod
    def from_json(cls, path, *, components: Sequence[str] | None = None) -> Parameters
    @classmethod
    def from_mapping(cls, mapping, *, components: Sequence[str] | None = None) -> Parameters
```

`from_mapping` walks the `epcsaft.parameters` document (`schema == "epcsaft.parameters"`,
`schema_version == 1`, else `ParameterError`). `components` selects and orders; `None`
means document order. It emits one `ParameterEntry` per numeric value in §2.1's table —
`identity` from the JSON, `field` from the table, `index` from `component_ids.index(...)` and
`site_ids.index(...)`, `magnitude` and `unit` copied — plus the structural lists, and calls
`_core.build_parameters`. What it reads: `components[].fixed.{charge_number, and molar_mass when present}`,
`components[].coefficients[]`, `pairs[]` (a pair is emitted only when both components are
selected; `component_id_a < component_id_b` in the document and the pair index is `[i, j]`
in selection positions; `k_ij_reciprocal_temperature_slope` carries a sibling
`reference_temperature` emitted as `k_ij_reference_temperature`), `correlations[]`
(`segment_diameter` in forms `constant` and `constant-plus-sum-of-exponentials`;
`relative_permittivity` in forms `constant` and `constant-plus-polynomial-in-temperature`;
any other `family`/`form` → `ParameterError` by `correlation_id`), `model_families[]`
(`base` must be `pc-saft`; `association` `none` or `general-site`; the `electrolyte` record's
`c_shell` and `c_dielectric` are emitted as `born_c_shell` / `born_c_dielectric` when the
choice is `born`), `model_coefficients[]`, `topology.sites[]` and `topology.edges[]` for
selected components (`source.kind` `explicit` → the two value entries; `combining-rule` with
`rule_id` `arithmetic-energy-geometric-volume` → one `association_combining_rule` entry
whose `inputs.energy_over_k` and `inputs.volume` identities are resolved to the two self
edges; any other `rule_id` → `ParameterError`). What it rejects by name: non-empty
`topology.presets`, a `k_hb_ij` pair coefficient, an unknown family, form, or model choice.
What it does not read: `domains`, `sources`, `purpose`, `reaction_correlations` (Stage 4),
`site_role`, and every provenance field.

**`eos.py`.**

```python
class Mixture:
    def __init__(self, parameters: Parameters)          # >= 2 components else ParameterError
    component_ids, site_ids, parameters
    def state(self, T: float, *, rho: float | None = None, P: float | None = None,
              x: Sequence[float], phase: str | None = None, anchor: State | None = None) -> State
    def reference(self, T: float, P: float, *, solvent_fractions: Sequence[float],
                  ions: tuple[str, str], stoichiometry: tuple[int, int] = (1, 1)) -> Reference
    def counters(self) -> _core.SweepCounters
    def reset_counters(self) -> None
class Pure:
    def __init__(self, parameters: Parameters)          # exactly 1 component
    def state(self, T, *, rho=None, P=None, phase=None) -> State     # x = (1.0,)
```

Exactly one of `rho`, `P`; `phase` in `{"liquid", "vapor", None}` maps to
`_core.DensityBranch`; anything else is `ParameterError`. With `P`, `state` calls
`self._model.state_at_pressure(T, P, x, branch, anchor.molar_density if anchor else None)`
— that method name must appear in `eos.py`, a mutation targets the stub. SI floats in, SI
floats out: K, Pa, mol/m³, mole fractions, mol/kg.

**`reference.py`.** `Reference` holds the model, the call inputs, the resolved
`ion_indices`, and the native `InfiniteDilutionReference`; exposes `log_fugacity`,
`molality`, `convergence_error`, `pure_solvent_log_fugacity` as attributes of the native
object; and

```python
    def activity(self, molality: float) -> Activity:
        return _core.electrolyte_activity(self._model, self.T, self.P, self.solvent_fractions,
                                          self.ion_indices, self.stoichiometry, molality, self._native)
```

**`_core.pyi`.** One stub for every name in §2.2, complete enough for `mypy --strict` on
`src/epcsaft`. It must contain `def state_at_pressure(`.

Not written: `units.py` (deleted), `Parameters.select`, `to_mapping`, `with_values`,
`ActiveParameter`, `parameter_hash` (Stage 6 is its first consumer), `Pure.saturation`
(Stage 5), domain validation, `transfer_to`, solvation free energy, `Activity.mole_fractions`,
any `__eq__`/`__hash__`/fingerprint on `Parameters`.

## 3. The tolerance: `derived(activity_reproduction_tolerance)`

The 24 rows of `aqueous_alkali_halide_activity` were captured with
`figiel-2025-reference-electrolytes/2` (`capture2.py`: `pk("figiel-2025-reference-electrolytes")`,
`solvent_x=(1,0,0)`, `formula_stoichiometry=(1,1)`, `T = 298.15 K`, `P = 1 bar`), so that
packet, selected as `("water", cation, anion)`, is what the gate evaluates them with. It is
the published model whose reproduction is a fact about the world in the sense of
`TESTING.md`.

**What the source reports.** Figiel, Yu & Held 2025 (transcription sha256
`ce80533925a91bc59d8d0d8056113c40611ca26c2edf04aced76986d50bd4bae`, re-read for this brief)
reports **no deviation statistic** for the aqueous MIACs of its Figure 5 — the text says the
approach describes them "very well" and reaches "quantitative agreement", and cites Hamer–Wu
1972 as the data. The one number the source publishes for MIAC adequacy is in its supporting
information, section S2: candidate Born-term approaches were screened against NaBr MIACs in
water and retained when **AARD < 10 %** (SI transcription sha256
`85bd39f727158d5a9d6eea6828c1673f73850e783a655b09660cc9b66d84321a`, lines 55–80, as recorded
in `engine/docs/regression/research/figiel-table-4-5-fitting-method.md`; the SI file itself
is not present on this machine, so the quotation is taken from that document). Missing fact,
stated plainly: a per-salt or pooled deviation for LiCl, LiBr, KCl, KBr reported by the
parameter source does not exist in the paper text.

**The constructor.** The gate applies the source's own adequacy criterion, as the same
statistic over the same kind of scope: for each salt, the average absolute relative deviation
of `mean_ionic_molality_scale` from `measured_mean_ionic_activity_coefficient` over that
salt's six rows must not exceed the source's criterion.

    aard(salt) = (1/6) * sum_rows |gamma_calc / gamma_measured - 1|  <=  SOURCE_ADEQUACY_AARD = 0.10

There is no factor: the source's statistic is an average, the gate's is the same average,
and the source applied it per salt (NaBr). The number is defined exactly once, as data with
its citation and hashes, in `engine/tests/eos/test_activity.py`; the check reads it from
there. Rule 6 governs `engine/native/` and this is a test constant in the style of Stage 2's
checks; no native source carries it.

**Why the average and not the maximum.** The superseded engine, evaluated on these 24 rows
with this packet (`analyses/reproduction/2025-figiel/results/data/figiel-2025-miac.csv`,
evidence only, never a gate value), gives per-salt AARD 3.36 % (LiCl), 4.13 % (LiBr),
3.34 % (KCl), 2.44 % (KBr) — and a single-row maximum of 10.11 % (LiBr at 2 mol/kg). A
per-row gate at the source's number would fail the parameter set itself on that row; the
average is the statistic the source used and the one the parameter set meets with at least
2.4× headroom.

**The floor argument.** Every solve under `Mixture.reference` carries its own tolerance, and
the gate sits above all of them:

| floor | value | effect on gamma |
|---|---|---|
| reference limit truncation at the last ladder point `m4 = 1e-12` | `A_gamma * sqrt(m4)` ≈ 1.17 × 1e-6 | 1.2e-6 relative |
| accepted reference convergence `REFERENCE_CONVERGENCE_TOLERANCE` | 5e-5 in ln phi | ≤ 5e-5 relative |
| density root `DENSITY_ROOT_RESIDUAL` (reduced pressure) | 1e-12 | negligible |
| association residual `ASSOCIATION_RESIDUAL` | 2e-13 | negligible |
| Hamer–Wu tabulation, three decimals | ±0.0005 absolute | ≤ 8.7e-4 relative at gamma = 0.573 |

The numerical floor is ≤ 5e-5 relative; the data floor 8.7e-4. The gate at 0.10 is 2,000×
the first and 115× the second: it tests the equations and the parameter set, not the solves
under them. The check also asserts every gamma finite and positive — a solve that fails
raises `EosError` before any statistic is formed.

## 4. Checks

### 4.1 `tests/eos/test_activity.py` — the oracle gate

Load `engine/tests/oracle/oracle.json`; group `aqueous_alkali_halide_activity` by `salt`;
for each salt build `Mixture(Parameters.from_json(<figiel packet>, components=("water",
cation, anion)))`, take `T` and `P` from the rows (assert all six agree), call
`reference(T, P, solvent_fractions=(1.0, 0.0, 0.0), ions=(cation, anion),
stoichiometry=(1, 1))`, then `reference.activity(row["molality_mol_kg"]).mean_ionic_molality_scale`
per row. Assert finite and positive; assert `aard(salt) <= SOURCE_ADEQUACY_AARD` with the
per-row deviations in the failure message. No literal other than `SOURCE_ADEQUACY_AARD`
and its citation comment. Print the 24 deviations under `-s` so acceptance can record them.

### 4.2 `tests/eos/test_state.py` — the mapping layer

1. **Relabelling is a symmetry through the mapping.** The MEA bundle
   (`mea-co2-h2o-nine-species-estimation/2`) loaded in the Stage 2 fixture order
   `[water, monoethanolamine, carbon-dioxide, protonated-monoethanolamine, carbamate-anion,
   bicarbonate-anion, carbonate-anion, hydronium-cation, hydroxide-anion]` and in the
   reversed order; `state(313.15, rho=3.0e4, x=<floor composition>)` with `x` permuted to
   match: `pressure`, `residual_helmholtz`, `compressibility_factor`, `bulk_relative_permittivity`
   equal to 1e-12 relative; `residual_chemical_potential_over_rt` and
   `log_fugacity_coefficient` equal after permutation to 1e-12. Floor composition
   `(0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0)`.
2. **A state satisfies its own identities.** On the same state: `sum(x_i * mu_i) ==
   residual_helmholtz + compressibility_factor - 1` to 1e-12; `pressure ==
   compressibility_factor * molar_density * R * T` to 1e-12 relative with `R = 8.314462618`
   written in the test with its source (CODATA 2018, exact); `state(313.15, P=1e5, x,
   phase="liquid").pressure == 1e5` to **1e-9** relative — `DENSITY_ROOT_RESIDUAL` is 1e-12 in
   reduced units, and a gate written at the floor would test the floor; then the anchored solve
   `state(313.15, P=1e5, x, phase="liquid", anchor=<that state>)` returns the same
   `molar_density` to 1e-10.
3. **Pure is the one-component mixture.** `Pure(Parameters.from_json(<gross-2001-methane-ethane>,
   components=("ethane",))).state(300.0, rho=1.0e4).pressure` equals
   `Mixture(Parameters.from_json(<same>, components=("methane", "ethane"))).state(300.0,
   rho=1.0e4, x=(0.0, 1.0)).pressure` to 1e-12 relative.
4. **Fails closed, by name.** (a) `Parameters.from_json(<figiel>, components=("methanol",
   "ethanol"))` raises `ParameterError` whose message contains `pair/ethanol/methanol/k_ij` —
   the source does not report that cell and the packet has no record for it; (b) the
   methane–ethane mapping with ethane's `fixed.molar_mass.unit` rewritten to `"gram"` raises
   `ParameterError` naming `component/ethane/molar_mass` and `gram`; (c)
   `Mixture.state(..., P=1e5, x=..., phase="gas")` raises `ParameterError`; (d) the
   methane–ethane mapping with ethane's `fixed.molar_mass` deleted raises `ParameterError`
   naming `component/ethane/molar_mass` — a neutral must carry its mass; the oracle gate is
   the evidence that a charged component need not.

`tests/packets.py` (≤ 40 lines): `PACKETS = Path(__file__).resolve().parents[2] / "data" /
"packets"`, `ORACLE = Path(__file__).resolve().parent / "oracle" / "oracle.json"`, and
`def parameters_json(packet_id: str, version: int) -> Path`. `tests/eos/__init__.py` is empty
so that `tests.packets` imports with `engine/` as the root.

## 5. Performance — `tests/oracle/measure_floor.py`, rewritten

The script in the tree opens the superseded API (`_core._instrumentation()`, Pint, old
fixture paths) and cannot run; rewrite it in place against §2.3. Reported, never pass/fail.
Assert `OMP_NUM_THREADS == "1"` at import. Median of 20 after 3 warm-ups,
`time.perf_counter_ns`. One row per case: `wall_us`, and the per-call delta of the four
counters. Print `git rev-parse HEAD`, the `/proc/cpuinfo` model name, `sys.version`,
`OMP_NUM_THREADS`; write the table as JSON to `argv[1]` if given. Compositions are passed
literally, never normalised; the nine-species tuple must appear in the script as
`0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0` — a mutation targets it.

| case | what | today | floor |
|---|---|---|---|
| A | methane–ethane (`gross-2001-methane-ethane/2`, order `("methane", "ethane")`), 300 K, rho = 1.0e4, x = (0.4, 0.6), `state` | 458 µs | ≤ 110 µs |
| B | nine species (MEA bundle, fixture order), 313.15 K, rho = 3.0e4, floor composition, `state` | 5,250 µs | ≤ native case A wall + 50 µs (449 µs on this machine today) |
| C | nine species, 313.15 K, 1e5 Pa, `phase="liquid"`, cold | 382 ms | ≤ native case D wall + 50 µs (7,298 µs today) |
| D | case C anchored on its own result | — | ≤ native case D anchored + 50 µs |
| E | `reference` for water / sodium-cation / chloride-anion from the figiel packet, 298.15 K, 1e5 Pa | 108 ms | ≤ native case F wall + 50 µs (5,687 µs today) |
| F | one `activity(1.0)` on E's reference | — | reported |

TRANSFER's "≤ 5 ms" for case C is below the native harness's own measured number, so the
Python figure is judged against the native harness run in the same session, not against a
fixed literal: Python cannot beat the native it wraps, and the difference is the marshalling
this stage owns. `recordings` per call must be 0 in steady state for A and B; report it.

## 6. Build, install, run

Before anything else, in the worktree root:

```bash
test "$(command -v python3)" = "$PWD/.venv/bin/python3"      # PATH leads with .venv/bin
touch "$HOME/.cache/probe" && rm "$HOME/.cache/probe"        # --add-dir "$HOME/.cache" is in effect
.venv/bin/meson configure engine/build/greenfield -Dprefix="$PWD/.venv"
```

The first is the regeneration trap: any edit to `meson.build` makes Meson re-resolve
`python3`, and `/usr/bin/python3` has no pybind11. The third routes `meson install` into
`.venv/lib/python3.13/site-packages/epcsaft/` (`python.install_env` stays `prefix`; the
`venv` setting needs `VIRTUAL_ENV` at every regeneration and is not used). `meson install`
rebuilds before it copies, so every Python gate is
`.venv/bin/meson install -C engine/build/greenfield --quiet && <python command>`; a stale
installed copy is the failure mode this avoids. Never `pip install`, never `PYTHONPATH`.

`meson.build`: add `'native/model/parameter_builder.cpp'` to `model_sources`; replace the
`_core` extension's `files(),` with `files('native/bindings/module.cpp',
'native/bindings/eos_bindings.cpp'),`. Nothing else.

`model_support.hpp`: refresh `nine_species()` to the living bundle, four edits —
`dispersion_energy_over_k[2]` `169.21` → `173.44025`; `relative_permittivity[2]` `1.4122` →
`0` (revision 2 has no permittivity record for carbon dioxide; it is no longer a dielectric
solvent in the mixing rule); `k_ij[2 * 9 + 0] = k_ij[0 * 9 + 2] = 0.013262879176919628`;
`ion_suppression_coefficient` all zero (`solvent-only`), and the fixture carries the comment `// revision 2 declares solvent-only; the physics layer spells that IonSuppression with alpha = 0` on that line. `born_c_shell = born_c_dielectric =
0` already matches revision 2. Then `meson test physics model` must still pass — every Stage
2 check is parameter-independent.

Then: `meson compile`, `meson test physics model`, `meson install --quiet`, `python -m pytest
-q engine/tests/eos`, `(cd engine && ../.venv/bin/mypy)`, the no-arithmetic grep,
`python3 tools/check-engine-architecture.py`, `OMP_NUM_THREADS=1 python
engine/tests/oracle/measure_floor.py` — all clean, worktree free of stray files
(`__pycache__` under `src/` and `tests/` is ignored; nothing else may appear in `git status`).

## 7. What to report

A JSON object: `files[]` {path, lines from `wc -l`, or `deleted`}; `gates[]` {id, command,
passed} for every gate in the JSON, run by you; `activity_rows[]` — the 24 `(salt, molality,
measured, calculated, relative_deviation)` tuples and the four per-salt AARDs; `floor` — the
§5 table verbatim plus the environment lines, and the native `meson test --suite floor`
lines from the same session for the comparison; `deviations[]` — every departure from this
brief and why; `policy_entries_added[]` (expected empty); `stopped_because` in {none,
budget_ceiling, gate_impossible, source_contradicted, environment_broken} with `detail`. Do
not run the mutations; acceptance runs them. Do not claim a number you did not measure.

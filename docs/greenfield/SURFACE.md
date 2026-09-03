# Surface — what a caller declares, and what a caller reads

Status: design, 2026-09-03, against `greenfield@152476b` (read at `41eb25b`, one brief-only
commit later; no source under `engine/` differs). Stage 3's design document, written after
the stage: `TRANSFER.md` named it and Stage 3 was built from `briefs/stage-03.brief.md`
alone, closed at `b14851f`, and reopened on the three findings in
`briefs/stage-03-repair.brief.md`. This document records the design those two briefs and
`STAGE-LOG.md` §"Stage 3" decided, states the tree as it stands, and names where the repair
changes it. It depends on `CORE.md` §4 (the declared-choice contract), §5 (the convergence
record), §6 (debug mode) and §8.3 (checker rules 10–12), and restates none of them.

Evidence classes as in `CORE.md`: **measured** (read from a file or a run on this tree
today), **derived** (follows from the code by an argument given here), **guessed** (marked).

## 1. What the surface is, in one rule

The surface is the three directories a caller touches and nothing in the engine touches:

| directory | rank | owns | lines (*measured*) |
|---|---|---|---|
| `native/surface/` | between `regression` and `bindings` | `parameter_builder.{hpp,cpp}`: turns a document into `Parameters`, or rejects it by name | 31 + 426 |
| `native/bindings/` | top | `module.cpp`, `eos_bindings.cpp`, `diagnostics_bindings.cpp`: the only three files that include pybind11 | 20 + 84 + 112 |
| `src/epcsaft/` | outside `native/` | `parameters.py`, `eos.py`, `reference.py`, `__init__.py`, `_core.pyi`: a mapping layer | 323 + 113 + 45 + 26 + 149 |

**Python does structure, C++ does numbers, and the surface is where a caller's decision
becomes loud.** Every entry point does one of three things: it resolves a *name* to an
*index* (a component id to a position, a site pair to two positions, a phase word to an enum
member), it hands what it was given to `_core` untouched, or it returns what `_core` returned
untouched. There is no fourth thing. The rule the owner set for the old tree's #190 is the
whole design: *any Python method is agnostic to what it is solving; it takes the inputs it is
given, gives them to C++, gets the answer back, and returns the output.*

Three consequences that are checked rather than hoped for:

- **No arithmetic in Python.** No `math`, `numpy`, `pint`, `scipy`, `decimal`, `fractions`,
  `statistics`, `cmath` import anywhere under `src/epcsaft/` — the `python_does_no_arithmetic`
  gate is a grep and its mutation is `import math`. No unit factor, no `exp`, no
  normalisation, no combining rule, no index arithmetic on a physical quantity. *Measured:*
  the only `import`s in the five files are `json`, `pathlib`, `collections.abc`, `typing`,
  and `_core`.
- **No default but `None`.** Checker rule 12 walks every function signature under
  `src/epcsaft/*.py` with `ast` and refuses any default other than the constant `None`
  (`CORE.md` §4.3). The surface can therefore not supply a physical value in place of a
  missing one; it can only pass the absence on.
- **Validation happens once, natively.** The builder rejects what only the document can show
  (§3.2); the model constructor rejects what only the arrays can show; Python rejects only
  what is not a document at all. The three-pass, two-language validation of the old tree is
  not written.

The include rank puts `surface/` above every engine layer so that `bindings/` may include it
and nothing in the engine can (`CORE.md` §10.2). The reason it is not `model/` is the reason
the dispersion selector was lost: a builder filed under the model reads as a mapping detail,
and a mapping detail's omissions are quiet.

## 2. The declared-choice contract at the surface

`CORE.md` §4 defines `core::Declared<T>` and the rule *a decision is declared or the call is
refused by name*. This section is what that rule means for the three things a packet
declares and for the strings that carry them across the boundary.

### 2.1 What a packet declares

An `epcsaft.parameters` document (`schema_version 1`) carries its model choices in
`model_families[]`, one record per `kind`. *Measured* over the 21 packet revisions under
`data/packets/`:

| `kind` | `choice` strings in use | carried how |
|---|---|---|
| `base` | `pc-saft` | must be `pc-saft`; anything else is rejected by `parameters.py` |
| `association` | `none`, `general-site` | read for the reject; the topology block is the actual declaration |
| `electrolyte` | `none` (4 revisions, all with zero ions), `fully-dissociated-debye-huckel` (2), `born` (15) | `ParameterDocument::electrolyte_choice` |
| `electrolyte` → `ionic_dispersion` key | `exclude-same-sign-ion-pairs` on every current revision with ions (8); absent on every neutral packet and on every superseded ion revision | `ParameterDocument::ionic_dispersion_choice`, the empty string when the key is absent |
| `permittivity` | `none` (4), `solvent-only` (3), `ion-fraction-suppression` (10), `ion-specific-suppression` (1), `component-permittivity-mixing` (2) | `ParameterDocument::permittivity_choice` |

Three of those are *selectors* in `core/parameters.hpp` — `ElectrolyteForm`,
`IonicDispersion`, `PermittivityForm` — and each is a `core::Declared<...>` field of
`Parameters` (*measured*, `parameters.hpp:76-80`). The builder is the one place that turns
the string into the `Declared` value.

### 2.2 Absence crosses as absence

`parameters.py` reads the `ionic_dispersion` key with
`str(electrolyte_record.get("ionic_dispersion", ""))` (*measured*, `parameters.py:320`) and
passes the empty string. It does not know that `exclude-same-sign-ion-pairs` is the
literature's convention, and it must not: the day it does, a packet that forgot the key gets
the right physics by accident and the next packet that meant `all-ion-pairs` gets the wrong
physics silently. The native side decides (*derived*, `parameter_builder.cpp:163-175`):

| `electrolyte_choice` | `ionic_dispersion_choice` | result |
|---|---|---|
| `none` | empty | `ionic_dispersion` undeclared |
| `none` | anything | rejected: `unconsumed model/electrolyte/ionic_dispersion` |
| not `none` | empty | rejected: `missing model/electrolyte/ionic_dispersion` |
| not `none` | `all-ion-pairs` / `exclude-same-sign-ion-pairs` | declared |
| not `none` | anything else | rejected by the string |

That table is the contract in miniature. An absent key is never an error in itself and never
a value; it is an error exactly when the other declarations say the key was needed, and the
message names the key. `tests/eos/test_state.py` holds this red-first against
`figiel-2025-reference-electrolytes/2`, the superseded revision that has ions and no key.

### 2.3 Why no selector carries an `Unspecified` or a `None` member

Two different mistakes, one rule against both (checker rule 11, `CORE.md` §8.3).

**`Unspecified` is a hole inside a type whose other members are physics.** The EOS repair
first wrote `IonicDispersion::Unspecified` and it worked; `CORE.md` §4.1 rejected it because
every reader must then remember which member is the hole, and a kernel can compare against
it without anyone noticing. `Declared<T>` puts the hole *outside* the type: `declared()` is
the only way to ask, and `require(name)` is the only way to read an undeclared value, which
throws `MissingDeclaration` naming it. `IonicDispersion` was therefore born with two members.
`DensityBranch::Unspecified` (`helmholtz_model.hpp:59`) is not a violation of this rule: it is
a *choice* meaning "either branch", and the caller who passes it gets the branch the root
finder finds. The name is wrong and `Either` would say what it means; rule 11 does not reach
`model/`, so this is a rename, not a gate.

**`None` is a member the surface cannot construct.** `PermittivityForm::None = 0` and
`ElectrolyteForm::None = 0` still exist at this commit (*measured*,
`core/parameters.hpp:23,31`). The builder never assigns either: the `else if (... != "none")`
guard falls through on `"none"` and leaves the field undeclared (*measured*,
`parameter_builder.cpp:156-185`). So `None` is dead inside a live selector, and the one
kernel read of it is a guard that never guards:

```cpp
if (parameters_.permittivity_form.is(PermittivityForm::None) == false) {   // helmholtz_model.cpp:242
```

`is()` on an undeclared value is `false`, so the condition is `true` for every parameter
set that exists, declared or not. That is the repair brief's finding 1, and its fix is two
deletions plus `declared()` in that line. After it, "nothing" has exactly one spelling on the
native side — an undeclared `Declared` — and one spelling on the packet side, the string
`none`, and the builder is the only translation between them.

### 2.4 What a packet saying `none` means

`none` is not "unspecified". It is a positive statement about the model, made by the packet
author, that the physics it names is absent — and it crosses the boundary as *undeclared*
because an undeclared selector is what the kernels already read as "no term". The model
constructor then holds the packet to its statement:

| selector | packet says | native value | constructor rule (`helmholtz_model.cpp`) | what a violation says |
|---|---|---|---|---|
| `electrolyte_form` | `none` | undeclared | declared **iff** a charged component exists | an ion under `none` → `missing declaration: electrolyte_form`; `born` with no ion → rejected by name |
| `ionic_dispersion` | (no key) | undeclared | declared **iff** a charged component exists (`require(has_ion == declared(), "ionic_dispersion")`) | same shape |
| `permittivity_form` | `none` | undeclared | today: none; after the repair, `declared()` gates the bulk-permittivity mixing | a neutral set is no longer required to carry a dielectric solvent — which is correct, and the repair records it as a behaviour change |

So `none` on a neutral packet is consistent by construction: there is nothing to declare, and
the *measured* fact that all four `none` packets have zero ions is what makes the `None`
members deletable rather than merely unused. A packet that says `none` and carries an ion is
lying about its model, and the constructor says so by the selector's name.

**Reported as well as required.** Every `evaluate` span emits `declared_choices`
(`CORE.md` §4.3), so a trace prints which model was asked for beside how it converged. The
surface makes an omission loud on the way in; the trace makes the declaration visible on the
way out.

## 3. The parameter builder

### 3.1 Its input is a list of facts, not a document

```cpp
struct ParameterEntry final {                 // surface/parameter_builder.hpp
    std::string identity;                     // the JSON identity — for messages only
    std::string field;                        // one of the names in §3.4
    std::vector<std::size_t> index;           // positions, resolved by Python from ids
    double magnitude = 0.0;
    std::string unit;                         // the packet's string, untouched
};
struct ParameterDocument final {
    std::vector<std::string> component_ids;   // selection order = model order
    std::vector<int> site_component, site_multiplicity;
    std::string electrolyte_choice, ionic_dispersion_choice, permittivity_choice;
    std::vector<ParameterEntry> entries;
};
Parameters build_parameters(const ParameterDocument&);   // throws std::invalid_argument
```

Python walks the JSON and emits one `ParameterEntry` per numeric value it finds, with the
`index` already resolved from ids to positions in the caller's selection order. The
`identity` travels only so that a rejection can name the record the packet author wrote. The
`unit` travels as the packet's string so the native side can refuse it. Nothing else about the
document — provenance, domains, sources, purpose, `site_role`, `reaction_correlations` —
reaches C++, because nothing in `Parameters` reads it.

This shape is why the builder is testable without Python and why the Python side has no
number to get wrong: a `ParameterEntry` with `magnitude 18.015` and `unit "gram / mole"` is
the same fact whichever language built it.

### 3.2 What the builder may reject, and why each rejection exists

Every rejection is a `std::invalid_argument` — `ParameterError` in Python — whose message
names the identity, the component, or the pair, and what was wrong. The categories, with the
reason each one is the builder's and not someone else's (*measured* line ranges in
`parameter_builder.cpp`):

| category | example message | why it is the builder's |
|---|---|---|
| **unknown choice string** (`:156-185`) | `unknown electrolyte choice X` | the packet vocabulary is closed; a string outside it is a transcription error, not a request for a default |
| **wrong unit for a field's dimension** (`:44-79`) | `unit gram for component/ethane/molar_mass: expected molar mass` | only the field knows its dimension; a document in g/mol read as kg/mol is wrong by 1000 and nothing downstream can see it |
| **index shape or range** (`:81-114`) | `invalid component pair for ...`, `site index out of range for ...` | Python resolved names to positions; a malformed position is a mapping bug and must not become a silent write past the end |
| **duplicate** (`:206`) | `duplicate component/water/segment_diameter` | two records for one cell means the packet disagrees with itself; last-writer-wins would hide it |
| **non-integral charge** (`:217-219`) | `charge ... is not integral` | `charge_number` is an `int` array and a fractional charge has no meaning in the electrostatic terms |
| **missing required record** (`:346-359, 389-394, 405`) | `missing component/lithium-cation/packing_diameter`, `missing pair/ethanol/methanol/k_ij`, `missing born_c_shell` | what is *required* depends on the declared choices and the charge — the diameters for every ion, the Born scalars and every `solvation_factor` under `born`, a per-ion `ion_suppression_coefficient` under `ion-specific-suppression`, a `k_ij` for every unlike pair that the dispersion convention does not exclude; only the builder holds both the choices and the coverage |
| **unconsumed record** (`:194-204, 360-376`) | `unconsumed model/relative_permittivity/...`, `unconsumed pair/lithium-cation/sodium-cation/k_ij under exclude-same-sign-ion-pairs` | a record the declared choices do not read is a probable transcription loss: either the value or the choice is wrong, and the builder cannot know which, so it refuses both |
| **combining rule without its inputs** (`:407-411`) | `combining rule self edge inactive: ...` | the rule derives a cross edge from two explicit self edges; evaluating it before both exist would fabricate a number |

**Unconsumed is the load-bearing category.** Rejecting `ion_suppression_coefficient` under
`solvent-only`, the Born fields under a non-Born electrolyte, and a same-sign `k_ij` under
the exclusion is the surface refusing to let two declarations disagree quietly. The same-sign
case is the one that matters historically: the EOS audit (`EOS.md` C.7) found the builder
*requiring* a `k_ij` record for same-sign ion pairs, which no source reports; with the
convention declared, those pairs are excluded from the term, and a record for one is
unconsumed rather than required. The builder cannot synthesise a self-pair, and for a two-ion
salt the same-sign pairs *are* the diagonal — which is why the convention had to be a
`Declared` selector and could not be packet data (`TRANSFER.md` §"Execution order").

**What the builder must not reject.** Three decisions, each recorded because each was once
the other way:

- **A zero value is data.** Packets retain printed zeros (an `ion_suppression_coefficient`
  of 0, a `k_ij` of 0). A present coefficient is never rejected for being zero; only an
  *absent* required one is `missing`. The `k_ij` completeness check therefore looks for the
  record's presence, not for a non-zero cell (`:378-386`).
- **A charged component may omit its molar mass** and the builder stores `0.0`
  (`:347-350`). Stage 3's first attempt stopped on this: no ion in any packet carries a mass,
  and molar mass never enters the Helmholtz energy — every native consumer (the permittivity
  mixing, the molality basis) reads solvent masses only. The constructor's rule is
  `charge_number != 0 || mass > 0` (`helmholtz_model.cpp:162-163`); a neutral must carry
  its mass, and the fail-closed test holds that. Whether molar mass is a *parameter* at all
  or a species property the engine should carry once was referred to the planner and is
  still open (§8).
- **A neutral component carrying `packing_diameter` or `born_diameter` is accepted** and
  the value is stored and never read. The old loader rejected it (`ionic-component-required`);
  `EOS.md` C.8 lists the drop so the decision is on record. It stays dropped: the field is
  harmless, and a rule that rejects a harmless value for tidiness is one more thing a
  packet author has to know.

**What the builder leaves to the constructor.** Array sizes, finiteness, positivity of
`segment_count` and `sigma_constant`, symmetry of the pair arrays, site-component range, and
the three declared-iff-ions rules of §2.4 (*measured*, `helmholtz_model.cpp:105-199`). Those
are facts about the arrays, and a hand-built `Parameters` that never passed through a
document must be held to them too. The builder does not duplicate them.

### 3.3 The arithmetic the builder is allowed

Three operations, all of them decisions and none of them physics a caller could vary:

1. **The unit factor** (§3.4). `magnitude * factor`, once, at the boundary.
2. **The association combining rule.** For `association_combining_rule` at
   `[s, t, s1, t1, s2, t2]`, after every explicit edge is filled:
   `energy(s,t) = ½(energy(s1,t1) + energy(s2,t2))`, `volume(s,t) = √(volume(s1,t1)·volume(s2,t2))`,
   both orders, `association_active = 1` (`:412-422`). No diameter factor: the Stage 4
   nine-species fixture encodes `0.0411037990944876 = √(0.03747·0.04509)` and the old engine
   evaluated exactly this at load. A packet revision that materialised Gross–Sadowski's
   diameter factor is a packet divergence, not a builder option.
3. **The ion-fraction fan-out.** Under `ion-fraction-suppression` Python emits the single
   model-level coefficient once, with an empty index, and the builder copies it to every
   charged component and `0` to every neutral one (`:395-406`). `solvent-only` and `none`
   zero the array; `ion-specific-suppression` requires one record per ion. `EOS.md` §B calls
   this the healthy counter-example to the dispersion defect: three published forms became
   one kernel equation with α as per-component data, and the *resolution* of the choice into
   data lives here, at the surface, where a missing coefficient is loud.

All three live in C++ for one reason: the moment one of them is written in
`parameters.py`, the rule "Python does no arithmetic" has an exception clause, and the value
of that rule is that it has none.

### 3.4 The unit table

The whole parameter vocabulary across every packet is nine strings; pressure units appear
only in domains and observations and never reach a parameter (`TRANSFER.md` §"Units").
*Measured*, `parameter_builder.cpp:44-79`:

| dimension | accepted → factor | fields |
|---|---|---|
| dimensionless | `dimensionless` → 1 | `segment_count`, `solvation_factor`, `relative_permittivity`, `ion_suppression_coefficient`, `k_ij`, `l_ij`, `association_volume`, the three Born/ionic-region scalars |
| kelvin | `kelvin` → 1 | `dispersion_energy_over_k`, `association_energy_over_k`, `k_ij_reciprocal_temperature_slope`, `k_ij_reference_temperature`, `permittivity_reference_temperature` |
| length | `angstrom` → 1 (the arrays are in Å) | `sigma_constant`, `sigma_term_amplitude`, the three diameters |
| charge | `elementary-charge` → 1, integral | `charge_number` |
| molar mass | `kilogram / mole` → 1, `kg / mol` → 1, `gram / mole` → 1e-3 | `molar_mass` |
| reciprocal temperature | `1 / kelvin` → 1 | `sigma_term_exponent` |
| reciprocal temperature to power *p* | `1 / kelvin ** p` for the term's recorded `power` (p = 1 spelled `1 / kelvin`) → 1 | `permittivity_term_coefficient` |

The `1e-3` is a conversion factor, not a threshold; rule 6 permits it. `kelvin` and
`angstrom` are accepted only for their own dimension: a length in kelvin is rejected even
though both factors are 1, because the check is on meaning, not on the number.

### 3.5 What is deliberately not mapped

- `PermittivityForm::DensityTemperature` and `::DipolarPerturbation`, and the five arrays
  only they read (`dielectric_volume_m3_per_mol` … `correlation_integral_parameter`). No
  packet carries their families; the builder zero-fills the arrays and no choice string
  reaches the forms. `EOS.md` C.4 names the two options — the builder learns them, or the
  forms and arrays go — and the answer belongs to whichever stage first needs a packet that
  declares one. The surface does not invent a choice string for a form nothing asks for.
- `domains`: the source's stated `T`/`P` range and ion-fraction limits. The old loader
  enforced them; the greenfield surface has no domain concept (`EOS.md` C.8). A domain is a
  fact about where the parameters were fitted, which is provenance, and the place to act on
  provenance is the analysis that chooses a packet, not the constructor of a model.
- `k_hb_ij` on a pair: rejected by name. The old `cross_association_energy_correction` has
  no successor and no packet carries the family as a coefficient (`EOS.md` C.12).
- `topology.presets`: rejected when non-empty. A preset is a name for a site set, and the
  surface maps sites, not names for sets of sites.

## 4. The bindings

### 4.1 Shape

`native/bindings/` is the only directory that includes pybind11 (checker rule 3 on the
include direction; rule 10 forbids the include under `core/`). Three files, one job each:
`module.cpp` registers the module and its four exceptions; `eos_bindings.cpp` binds the
model surface; `diagnostics_bindings.cpp` binds the record, the trace and the capture
(`CORE.md` §6.5). *Measured* against `eos_bindings.cpp`:

| bound | as | rule |
|---|---|---|
| `ParameterEntry`, `ParameterDocument` | `py::class_`, `py::init<>()`, `def_readwrite` on every field | the only structs Python *fills* |
| `build_parameters` | free function | returns an opaque `Parameters` |
| `Parameters` | opaque; `component_count`, `site_count` read-only | Python never reads an array back |
| `DensityBranch` | `py::enum_` with `Unspecified`, `Vapor`, `Liquid` | the one enum Python constructs, from a phase word |
| `State`, `SweepCounters`, `InfiniteDilutionReference`, `ElectrolyteActivity` | `def_readonly` on every field, same names, returned as they are | results are never copied into dataclasses or dicts |
| `HelmholtzModel` | `py::init<Parameters>()`; `state`, `state_at_pressure` (with `anchor_density = None`), `counters`, `reset_counters` | |
| `infinite_dilution_reference`, `electrolyte_activity` | free functions; `ion_indices` as `std::array<size_t,2>`, `stoichiometry` as `std::array<int,2>` | |

Exceptions, in registration order so pybind11 tries the derived type first:
`ParameterError` ← `std::invalid_argument` (a `ValueError`); `EosError` ←
`std::runtime_error` (a `RuntimeError`); `ConvergenceError` ← `core::ConvergenceFailure`,
under `EosError`; `DeclarationError` ← `core::MissingDeclaration`, under `ParameterError`.
A caller who catches `ParameterError` therefore also catches a missing declaration, and one
who wants to distinguish it can.

**Not bound, by decision:** `evaluate`, `solve_density`, `formula_composition`, the tape,
active parameters, `packing_fraction_coefficient` and its siblings. Each is either an engine
internal or a Stage 5/6 consumer's concern. `<pybind11/stl.h>` converts vectors, arrays,
optionals and tuples; no custom caster is written. The GIL is not released: the model's mutex
serialises `evaluate` and nothing on the surface threads.

### 4.2 The Scope rule: every bound function opens exactly one root

`CORE.md` §6 makes debug mode a lexical capture around one call, and a captured call's
trace is a tree with one root. That holds only if **every function reachable from
`bindings/` opens a `core::Scope` at its top**. If one does not, its children — `evaluate`,
`solve_density`, `association_newton` — become roots of their own, and a capture around one
Python call has several roots with nothing tying them to the call that made them.

*Measured* at this commit, over the five functions `eos_bindings.cpp` binds to a definition
under `native/`:

| bound function | opens | state |
|---|---|---|
| `HelmholtzModel::state_at_pressure` | `state_at_pressure` (`helmholtz_model.cpp:623`) | one root |
| `infinite_dilution_reference` | `infinite_dilution_reference` (`infinite_dilution_reference.cpp:69`) | one root |
| `HelmholtzModel::state` | nothing | **its `evaluate` is the root** — a capture around `Mixture.state(T, rho=…)` has as many roots as evaluations |
| `electrolyte_activity` | nothing | **its states are the roots** |
| `build_parameters` | nothing | reaches no instrumented site; a capture around it is empty |

The two missing spans survived nineteen gates because `trace_case` in
`diagnostics_checks.cpp` drives only `state_at_pressure` (*measured*, `:206-236`): the check
and the defect were written against the same blind spot. That is the repair brief's finding
2, and its fix has three parts, of which the third is the one that matters: two `Scope`
lines, two more `trace_invariance` cases, and **checker rule
`rule_bound_functions_open_a_scope`** — for each symbol bound in `native/bindings/`, its
definition opens a `core::Scope`. A rule in prose with nothing reading it is the defect, not
the fix.

The rule is written without exceptions. `build_parameters` reaches no span and so cannot
produce a second root; giving it one anyway costs one thread-local read (2.2 ns,
`CORE.md` §6.2) and gives the trace a natural place for a surface-side event if one is ever
wanted. The repair's `writes` list names only the two solver functions, so whether the
checker's predicate counts accessors (`build_parameters`, `counters`, `reset_counters`) is
decided by the rule as written at acceptance; §8 carries it as an open item rather than
pre-empting the coordinator.

**What the bindings themselves do not do.** They never open a span, never read the trace,
and never name `active()`. `Capture` is named in `diagnostics_bindings.cpp` and nowhere else
outside `core/` (rule 10). The `Capture` wrapper's `trace`, `spans` and `events` are
`reference_internal` to the wrapper so a trace cannot outlive its capture — the dangling
reference the Stage 2 review found and `53b9007` fixed.

## 5. The Python layer

### 5.1 The API, and what each call maps

```python
Parameters.from_json(path, *, components=None) -> Parameters
Parameters.from_mapping(mapping, *, components=None) -> Parameters
    .component_ids: tuple[str, ...]              # selection order = model order
    .site_ids: tuple[tuple[str, str], ...]       # (component_id, site_id), document order

Mixture(parameters)                              # >= 2 components, else ParameterError
Pure(parameters)                                 # exactly 1
    .state(T, *, rho=None, P=None, x, phase=None, anchor=None) -> State      # Pure: x = (1.0,)
    .reference(T, P, *, solvent_fractions, ions, stoichiometry) -> Reference # Mixture only
    .counters() -> SweepCounters; .reset_counters()

Reference.activity(molality) -> Activity
    .log_fugacity, .molality, .convergence_error, .pure_solvent_log_fugacity   # the native fields
```

SI floats in, SI floats out: K, Pa, mol/m³, mole fractions, mol/kg. Each method is one
mapping (*measured*, `eos.py`, `reference.py`):

| call | what Python does | what it does not do |
|---|---|---|
| `from_mapping` | checks `schema`/`schema_version`; selects and orders components; walks `components[].fixed`, `coefficients[]`, `correlations[]`, `pairs[]`, `model_families[]`, `model_coefficients[]`, `topology.{sites,edges}`; emits entries with positions; calls `build_parameters` | read a magnitude for any purpose other than copying it; supply a value for an absent record — the one `charge_number = 0` it emits when `fixed` has no charge is a schema convention exercised by 82 of 104 neutral components (`EOS.md` C.5), and it is emitted as a record so the builder still sees it |
| `state` | refuses both or neither of `rho`/`P`; maps `phase` ∈ {`None`, `"vapor"`, `"liquid"`} to `DensityBranch`, anything else `ParameterError`; forwards `anchor.molar_density` or `None` | normalise `x`; check its length; check `T > 0` — the native `validate_point` does all three and names the field |
| `reference` | resolves the two ion ids to positions, `ParameterError` if unknown; forwards | default the stoichiometry — `(1, 1)` was the hidden physical default `EOS.md` C.5 found, and its removal is held by `test_reference_requires_stoichiometry` (a missing keyword is a `TypeError` at the call, before any native code runs) |
| `activity` | forwards the stored inputs and the native reference | `exp()` the log-scale value — the native side returns `mean_ionic_molality_scale` because the oracle is stated in it |

`Parameters.__init__` takes the native object and the two id tuples; it is the one
constructor a caller never calls directly. `Pure` is the one-component `Mixture` with `x`
fixed at `(1.0,)`, and `test_pure_is_one_component_mixture` holds that they are the same
model to 1e-12.

### 5.2 Exports

*Measured*, `__init__.py`: `Parameters`, `ParameterError`, `EosError`, `Pure`, `Mixture`,
`State`, `Reference`, `Activity` — the eight the Stage 3 brief named, of which `State`,
`Activity`, `ParameterError` and `EosError` are re-exports of `_core` names, not wrappers —
plus `trace`, `Stop`, `Convergence`, `ConvergenceError`, `DeclarationError` from Stage 2.
Thirteen names. `TRANSFER.md` §"Stage 3" says *eleven*, as a cap on the EOS surface with
`Saturation` and `ActiveParameter` still to come; the five diagnostics names are the debug
mode's surface and were not in that count. §8 asks for the ruling rather than assuming one.

Not written, and each was a decision in the brief: `units.py` (Pint is gone, so there was
nothing in it), `Parameters.select`, `to_mapping`, `with_values`, `parameter_hash`,
`ActiveParameter`, `Pure.saturation`, `transfer_to`, any `__eq__`/`__hash__` on
`Parameters`, any domain validation. A name arrives with the stage that consumes it.

### 5.3 The stub is the contract

`_core.pyi` carries one declaration for every bound name, and `mypy --strict` runs over
`src/epcsaft` **and** `tests` (`STAGE-LOG.md` §"Stage 2": a stub nothing consumes is a stub
nothing checks; the config was widened after a gate that could not fail was found). A field
renamed in a binding and not in the stub is red in the tests that read it. The stub is
therefore as much a part of the surface as the binding, and a stage that binds a new name
edits both or its `types_consistent` gate goes red.

## 6. Marshalling, measured

The Stage 3 brief judged Python's overhead against the native floor harness run in the same
session, because Python cannot beat the native it wraps and the difference is the marshalling
this layer owns. *Measured* on this tree (AMD Ryzen 5 5500, gcc 13.3.0, Python 3.13.14,
`OMP_NUM_THREADS=1`, median of 20 after 3 warm-ups):

| Python case (`measure_floor.py`) | wall | native case (`model_floor`) | wall | difference |
|---|---|---|---|---|
| B — nine species, 313.15 K, ρ = 3.0e4, `state` | 425.0 µs, 12 forward sweeps, 0 recordings | A — the same evaluation on the same fixture | 422.0 µs, 25 sweeps | **≈ 3 µs** |
| A — methane–ethane, 300 K, ρ = 1.0e4, `state` | 34.5 µs | — | — | against the brief's ≤ 110 µs floor; the old tree's 458 µs |
| C — nine species, 1e5 Pa, `liquid`, cold | 3,349 µs, 85 sweeps, 4 recordings | D — cold density on the same fixture | 2,832 µs, 152 sweeps, 1 density iteration | not the same solve |
| E — reference, water / sodium / chloride | 2,645 µs, 199 sweeps | F — the native reference | 2,304 µs, 393 sweeps | not the same solve |

The state call is the one the brief singled out — "Python marshalling, `Mixture.state`
only" (`TRANSFER.md` §"Order of the performance work") — and its overhead is three
microseconds on a 422 µs call, below the harness's own run-to-run noise. The old tree's
590 µs of Python per call is gone because the call is a positional forward of three
arguments.

The other two pairs cannot be read as marshalling: the sweep and iteration counts show the
two harnesses solving different things on the same fixture (the native cold case takes one
density iteration; the Python one takes a scan). The brief's "native case wall + 50 µs"
floors for C, D and E are therefore not evaluable as written, and the fix is on the harness
side — one case list shared by both — not on the surface. Reported, never asserted, per
`CANON.md`.

## 7. What is tested, and what is only reported

The surface's checks are the ones `CANON.md` permits: measured rows, exact symmetries,
identities that hold for any parameter set, and fail-closed validation by name. *Measured*,
`tests/eos/`:

| check | kind | what would break it |
|---|---|---|
| `aqueous_alkali_halide_activity`: per-salt AARD of `mean_ionic_molality_scale` ≤ the source's own 10 % adequacy criterion, 24 rows, `figiel-2025-reference-electrolytes/3` | oracle | a wrong unit factor, a lost selector, a wrong combining rule, a wrong ion index |
| `aqueous_lithium_bromide_published_model`: 27 digitised points within their own digitisation uncertainty on average (added by the EOS repair; 0.2723 % AARD) | published-model reproduction | the dispersion convention not reaching the kernel |
| relabelling: the nine-species bundle loaded in fixture order and reversed, `x` permuted, every scalar equal and every vector equal after permutation to 1e-12 | symmetry through the mapping | any index arithmetic error in `parameters.py` or the builder |
| identities: `Σ x_i μ_i = a_res + Z − 1`; `P = Z ρ R T` with the CODATA gas constant; `state(P=…).pressure = P` to 1e-9; the anchored solve returns the same root to 1e-10 | closed form | a field bound to the wrong member |
| `Pure` = one-component `Mixture` to 1e-12 | identity | |
| fails closed by name: the superseded ion revision with no dispersion key; a selected pair with no `k_ij`; a same-sign `k_ij` under the exclusion; `gram` for a molar mass; `phase="gas"`; `reference` without `stoichiometry` | validation | a default creeping in anywhere on the path |
| trace round trip: a captured `state(P=…)` has one root whose record equals `state.density_convergence` | shape | the Scope rule (for this one function) |
| `python_does_no_arithmetic` (grep), `types_consistent` (mypy), checker rules 11 and 12 | static | |

Each gate has a named mutation that turns it red; acceptance re-runs both
(`ORCHESTRATION.md`). The activity criterion is the one number the design owed: the source
reports no per-salt deviation for these salts, and its SI's `AARD < 10 %` screening
criterion is the only adequacy statement it publishes, applied here as the same statistic
over the same scope (`briefs/stage-03.brief.md` §3). The floor argument under it — reference
truncation 1.2e-6, reference convergence 5e-5, tabulation 8.7e-4 relative — puts the gate
115× above the data floor, so it tests the equations and the parameter set, not the solves.

**Reported, never asserted:** every timing in §6, the counters per call, the per-row activity
deviations printed under `-s`. The open finding in `STAGE-LOG.md` §"Stage 3" — the activity
drift with molality that the EOS repair then explained — is the case for why: the gate
passed at 5.1–7.1 % while the physics was wrong, and the number that found the defect was a
comparison the gate could not make. An average over six rows cannot see a missing upturn in
the last one.

## 8. Open items and missing facts

- ~~**The repair is not landed at this commit.**~~ **Landed at `3c4e2be`**, after this
  document was written. Every item it listed is now true of the tree: both `None = 0`
  members deleted, the guard reading `.declared()`, a `Scope` in `HelmholtzModel::state`
  and in `electrolyte_activity`, and `rule_bound_functions_open_a_scope` as rule 13 with
  `rule_no_hand_named_residual_stop` as rule 14. §2.3 and §4.2 now describe the tree, not
  only the design. *Coordinator, on landing.*
- **Accessors under the Scope rule: ruled 2026-09-03, no exception.** The repair had
  written rule 13 with an allowlist of three names; the owner removed it. `build_parameters`,
  `counters` and `reset_counters` open a span, the allowlist is an empty set, and `CORE.md`
  §6.1 now says the same as §4.2 here.
- **Export count.** Thirteen names against `TRANSFER.md`'s "eleven"; the five diagnostics
  names were not in that count (§5.2). A ruling, not a change.
- **Molar mass as a species property.** Referred to the planner by the Stage 3 first attempt
  and still open. Today: a parameter, optional for ions, required for neutrals (§3.2).
- **The two unreachable permittivity forms** (`EOS.md` C.4): the surface takes no position
  until a packet declares one.
- **Marshalling for the cold-density and reference calls** is unmeasured because the two
  harnesses do not run the same case (§6). Owed by whichever brief next touches
  `model_floor.cpp` or `measure_floor.py`: one shared case list.
- **`DensityBranch::Unspecified` → `Either`.** A rename; rule 11 does not reach it (§2.3).
- **The `permittivity: none` behaviour change.** After the repair, a neutral set is no longer
  required to carry a dielectric solvent. Correct, and recorded in the repair brief; no
  packet in the tree is affected (*measured:* the four `none` revisions have no ions and
  every `solvent-only`/suppression revision has a solvent).

## 9. Where the tree departs from the Stage 3 brief, and why each is right

| brief said | tree has | why |
|---|---|---|
| builder in `native/model/` | `native/surface/` | `CORE.md` §10.2: the only file whose job is to make an omission loud is not a model detail |
| `ParameterDocument` with two choice strings | three: `ionic_dispersion_choice` added | the EOS repair; the convention had to be declared, and the surface is where declarations enter |
| `ElectrolyteForm::None` / `PermittivityForm::None` assigned on `none` | never assigned; undeclared | §2.3; the members are deleted by the repair |
| `stoichiometry=(1, 1)` default | no default; rule 12 | `EOS.md` C.5; a 2:1 salt with the default computed a 1:1 reference silently |
| every component requires `molar_mass` | neutrals only | no ion in any packet carries one; the constructor's `charge != 0 || mass > 0` is the rule |
| eight exports | thirteen | Stage 2 added the debug mode's five (§5.2, §8) |
| `units.py` rewritten | deleted | Pint left; nothing remained |
| `fixture_reads_living_bundle` PASS on first report | three of four refreshes had landed; caught by re-running the gate | `STAGE-LOG.md`: a gate summary is never acceptance |
| one oracle gate | two: the LiBr published curve added | the repair's evidence that the convention reaches the kernel, now a permanent row |

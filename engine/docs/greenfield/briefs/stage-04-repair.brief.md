# Stage 4 repair — the EOS layers re-accepted together

You are repairing Stages 3 and 4 of a greenfield rebuild in this worktree (branch
`stage-04-repair`, cut from `greenfield@eb7c728`). You have no one to ask. Read, in order,
before writing anything:

1. `engine/docs/greenfield/CANON.md` — a hard restriction on what may be tested
2. `engine/docs/greenfield/EOS.md` §1, §3 and §5 — the slope and anchor defects, measured
3. `analyses/reproduction/2025-figiel-engine-crosswalk/like-ion-dispersion-convention.md`,
   Addenda 1–4 — the dispersion defect, its registry entry, and why the fix is not packet data
4. `engine/native/physics/helmholtz.hpp`, `parameters.hpp`; `engine/native/model/*.cpp,*.hpp`;
   `engine/tests/native/model_support.hpp`, `model_checks.cpp`, `physics_checks.cpp`
5. `engine/src/epcsaft/parameters.py`, `engine/tests/eos/*.py`, `engine/tests/packets.py`
6. `data/docs/packet-contract.md` §"Manifest" and §"Consumer lock"
7. `engine/docs/greenfield/briefs/stage-04.repair.brief.json` — sources, budgets, gates, mutations

Verify every source marker and every precondition in the JSON against the tree **before**
writing. If a claim is false, stop and say which.

## 0. Non-negotiable

**Layers.** `physics <- model <- equilibrium <- regression <- bindings`. Includes rooted at
`engine/native`, never relative. Only `native/bindings/` includes pybind11.

**Python does no arithmetic.** It reads one more string from the packet and hands it to
`_core` untouched. Every check on it happens natively, once.

**No hidden physical default.** The dispersion convention is declared or the set is rejected
by name — natively for a `Parameters` struct, natively for a packet. No enum member, builder
branch, Python `.get(..., default)`, or fixture may supply one. This brief makes the missing
declaration a constructor error (§2.1), which is what makes the "no default" gate red-first.

**Constants and tolerances.** No new numeric literal acting as a threshold in `engine/native/`.
This repair adds no policy entry and deletes one (§2.4). Rule 6 is enforced by the checker.

**Naming.** No proper nouns in identifiers. `IonicDispersion::ExcludeSameSign`, not `Figiel`;
`AllIonPairs`, not `Ji`. Sources are cited in comments and `// EqID:` annotations only.

**Testing.** Every check below is a closed form, an exact limit, a symmetry, derivative
agreement, fail-closed validation, or a published-model reproduction with its own stated
uncertainty. No golden files, no assertions on internals, no iteration counts as correctness.
Counts and timings are reported, never pass/fail.

**Physics is open for exactly two edits**, both in §2.1: one enum and one field in
`parameters.hpp`, and the branch on the annotated line of `helmholtz.hpp`. No other
kernel arithmetic, tolerance value or signature changes.

**Data is in scope.** `engine/AGENTS.md`: an engine task owns its full scientific slice across
`data/`, `engine/` and `analyses/`. Eight packets gain a revision (§2.6). Existing revisions
are immutable: never edit a file under an existing `<packet>/<version>/`.

**Environment is not yours to repair.** `meson.build` is not edited at all in this stage.
You may not touch the dependency block, `default_options`, the virtual environment,
`uv.lock`, `engine/build/subprojects`, or how you were invoked. `PATH` leads with the
worktree's `.venv/bin`; `$HOME/.cache` is writable. Stop with `environment_broken`.

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

Three defects, all verified on `c860395` before this brief was written (measured):

| defect | where | measurement |
|---|---|---|
| like-ion dispersion not excluded | `physics/helmholtz.hpp:243-251` | greenfield dispersion equals the superseded engine's `all-ion-pairs` branch to 12 digits (−8.852592255706 at KCl 1 mol/kg, fixed ρ); 49.19 % AARD against Figiel's published LiBr curve, 0.27 % under the exclusion |
| `PressureSlope` is a frozen-site derivative | `model/helmholtz_model.cpp:385` | Φ_VV ratio `PressureSlope/Hessian`: 1.0000 ethane, **14.51** associating binary, 0.679 aqueous electrolyte, 0.807 nine species; the centred difference agrees with `Hessian` in every case; NaCl 280 K 5e5 Pa throws through the public API |
| anchor bypasses branch selection | `model/density.cpp:318-340` | ethane 280 K 1e6 Pa and 250 K 5e5 Pa: `liquid` with a vapour anchor returns the vapour root (475.5, 257.3) labelled `Liquid`, count 1 |

### 2.1 The selector — the registered case split gets its selector (Stage 4)

`equations_registry.yaml` registers `same_sign_ionic_dispersion` as a case split whose caller
"selects all ion pairs or same-sign exclusion as part of the electrolyte model family". The
annotated line implements the first branch unconditionally. The selector lived in the deleted
Python layer (`main@9eccbfd:_resolve.py:735-753`), which forced `k_ij = 1` on every same-sign
pair *including the diagonal*. That mechanism is not carried: it needs the temperature slope
zeroed with it, and it leaves a same-sign `k_ij` activatable as a fit parameter with a
non-zero derivative. The branch goes on the equation, where the registry says it is.

`physics/parameters.hpp` — add, above `struct Parameters`:

```cpp
// The registered case split same_sign_ionic_dispersion. Declared by the parameter
// source, never defaulted: a set with a charged component and Unspecified is rejected
// by the model constructor, and a set with no charged component may not declare one.
enum class IonicDispersion {
    Unspecified = 0,
    AllIonPairs = 1,      // every pair takes sqrt(eps_i eps_j)(1 - k_ij)
    ExcludeSameSign = 2,  // pairs with z_i z_j > 0, the diagonal included, take 0
};
```

and in `struct Parameters`, under `// Electrolyte`:

```cpp
    IonicDispersion ionic_dispersion = IonicDispersion::Unspecified;
```

`physics/helmholtz.hpp` — the pair loop, replacing the two lines that compute `epsilon_ij`:

```cpp
            // EqID: cross_dispersion_energy
            // EqID: same_sign_ionic_dispersion
            const bool same_sign_excluded =
                parameters.ionic_dispersion == IonicDispersion::ExcludeSameSign
                && parameters.charge_number[i] * parameters.charge_number[j] > 0;
            const Scalar epsilon_ij = same_sign_excluded
                ? Scalar(0.0)
                : scalar_sqrt(epsilon_i * epsilon_j) * (Scalar(1.0) - k_ij);
```

The loop runs `i` and `j` over every component including `i == j`, so the diagonal is covered.
The predicate is built from `int` parameters, so it is parameter-constant and needs no tape
cache key. The `k_ij(T)` slope is unreachable for an excluded pair by construction.

`model/helmholtz_tape.hpp` — `TapeParameters` gains `IonicDispersion ionic_dispersion;` and
`make_tape_parameters` copies it. Nothing else in the tape changes.

`model/helmholtz_model.cpp` constructor — after the `dielectric solvent` requirement:

```cpp
    const bool has_ion = std::any_of(
        parameters_.charge_number.begin(), parameters_.charge_number.end(),
        [](int charge) { return charge != 0; }
    );
    const bool declared = parameters_.ionic_dispersion != IonicDispersion::Unspecified;
    require(has_ion == declared, "ionic_dispersion");
```

Exactly that line: the mutation for gate `no_hidden_dispersion_default` is a `sed` on it.

`tests/native/model_support.hpp` — `aqueous_electrolyte()` and `nine_species()` each set
`p.ionic_dispersion = epcsaft::IonicDispersion::ExcludeSameSign;` (both fixtures are
transcribed from Held-lineage sources; `like-ion-dispersion-convention.md`). The neutral
fixtures and `physics_checks.cpp`'s `neutral_binary()` leave it `Unspecified`.

### 2.2 The declaration reaches the builder and the packet (Stage 3)

`model/parameter_builder.hpp` — `ParameterDocument` gains
`std::string ionic_dispersion_choice;` (empty when the packet carries no key).

`model/parameter_builder.cpp`:

- after the electrolyte/permittivity choice mapping: if `electrolyte_choice == "none"` and
  the string is non-empty, `reject("unconsumed model/electrolyte/ionic_dispersion")`; if
  `electrolyte_choice != "none"` and the string is empty,
  `reject("missing model/electrolyte/ionic_dispersion")`; `"all-ion-pairs"` →
  `AllIonPairs`, `"exclude-same-sign-ion-pairs"` → `ExcludeSameSign`, anything else →
  `reject("unknown ionic_dispersion choice " + value)`.
- in the post-pass over pairs (the loop that today rejects `missing pair/.../k_ij`): compute
  `same_sign = result.charge_number[i] * result.charge_number[j] > 0`. Under
  `ExcludeSameSign` a same-sign pair is **not required**, and a `k_ij`,
  `k_ij_reciprocal_temperature_slope` or `k_ij_reference_temperature` entry on it is rejected:

  ```cpp
  reject("unconsumed " + identity + " under exclude-same-sign-ion-pairs");
  ```

  where `identity` is the entry's `where(entry)`. Exactly that line — gate
  `unconsumed_same_sign_record_rejected` mutates it. Under `AllIonPairs` same-sign pairs are
  required like any other pair. This is the old loader's `not_applicable` set, fail-closed
  instead of silently overwritten.

`bindings/eos_bindings.cpp` — one `def_readwrite("ionic_dispersion_choice", ...)`.
`src/epcsaft/_core.pyi` — one line.

`src/epcsaft/parameters.py` — the electrolyte family record is located for every electrolyte
choice (today only for `"born"`); the document gets
`native_document.ionic_dispersion_choice = str(electrolyte_record.get("ionic_dispersion", ""))`
— the empty string means "absent" and the native side decides. No default value may appear
here; gate `python_eos_tests` mutates exactly this `.get`.

`src/epcsaft/eos.py` — `Mixture.reference(..., stoichiometry: tuple[int, int])` loses its
`= (1, 1)` default. A formula stoichiometry is a physical fact about the salt, not an API
convenience. Update the one call in `tests/eos/test_activity.py` to pass `(1, 1)`.

### 2.3 The pressure slope is the total derivative (Stage 4)

`model/helmholtz_model.cpp` `evaluate`: the `Hessian` branch already forms the S×S site
Hessian `Φ_XX` and its `FullPivLU`. Lift that block into an anonymous-namespace helper
returning the factorisation (its rank check included), called by both orders. In the
`PressureSlope` branch, `first[1 + site]` from `Forward(1, e_V)` is `Φ_{X_a V}` for every
site — it is already computed and discarded today. The branch ends:

```cpp
        double site_correction = 0.0;
        if (!sites.empty()) {
            Eigen::VectorXd cross_v(static_cast<Eigen::Index>(sites.size()));
            for (std::size_t site = 0; site < sites.size(); ++site) {
                cross_v(static_cast<Eigen::Index>(site)) = first[1 + site];
            }
            site_correction = cross_v.dot(factorization.solve(cross_v));
        }
        result.volume_second_derivative = second[3] - site_correction;
        return result;
```

`second[3]` is `second[2 * 1 + 1]`, the raw tape value, and the assignment line is written
exactly as shown: gate `pressure_slope_is_total` mutates it back to `second[3];`. The sweep
count for `PressureSlope` stays at 4; one S×S solve is added.

**Delete the dead double path** in `evaluate`: the `electro`, `born`, `direct` and
`direct_value` locals (today lines 324–338 and the `(void)direct_value;`) are computed on
every evaluation and never read. `sites` stays — it is the tape's dynamic parameter. This is
a deletion with no gate; the floor harness reports the change.

`density.cpp` `at_density` and the far-side check (§2.5) evaluate at
`DerivativeOrder::Gradient` when only the residual is needed, `PressureSlope` when the Newton
step needs the slope (`EOS.md` §5 change 2; answers bit-identical on 20 of 20 rows).

### 2.4 One threshold, one name

`predicate_signature` reads `policy::ELECTROLYTE_SERIES_THRESHOLD` while the kernel it
mirrors reads `policy::SCREENING_SERIES_SWITCH` — two entries, one meaning, both `1e-3`. If
either moves, the tape cache key and the kernel branch disagree and a cached tape is reused
on the wrong side of the switch with no error. `predicate_signature` reads
`SCREENING_SERIES_SWITCH`; `ELECTROLYTE_SERIES_THRESHOLD` is deleted from
`numerical_policy.hpp`. No gate: the checker's rule 5 matches `constexpr` only, and a
`Tolerance` is `inline const`. Report this gap; do not widen the checker in this stage.

### 2.5 The anchor is honoured only on its own branch (Stage 4)

`model/density.cpp` `solve_density`:

- An anchor is used only with a declared branch. `Unspecified` with an anchor ignores the
  anchor and runs the cold path (the cold path is the only one that can enumerate roots and
  therefore the only one that can say "ambiguous").
- The anchored Newton runs as today. On convergence at `ρ*` with `slope > 0`, the far side
  of the grid is checked before the root is returned:

  ```cpp
  const bool far_side_clear = far_side_has_no_crossing(
      model, temperature, pressure, mole_fractions, density, branch,
      log_lower, log_upper, intervals);
  if (far_side_clear) {
      return {density, point.residual, point.slope, iteration + 1, std::nullopt, branch};
  }
  // otherwise fall through to the cold scan below
  ```

  written exactly so — gate `anchor_respects_branch` mutates `if (far_side_clear)`.
- `far_side_has_no_crossing` walks the same logarithmic grid the cold scan uses (`intervals`
  points between `log_lower` and `log_upper`), restricted to grid densities **above** `ρ*`
  for `Liquid` and **below** `ρ*` for `Vapor`, at `DerivativeOrder::Gradient`. For `Liquid`
  every far-side residual must be `> 0.0`; for `Vapor` every far-side residual must be
  `< 0.0`. Any other sign, or an association failure at a grid point, returns `false`.
  Why this is exact for the convention the cold path already uses (derived): the residual
  is negative at the dilute end and positive at the dense end, so the smallest root and the
  largest root both cross upward and are stable; a crossing beyond `ρ*` on the requested side
  means a root exists that the cold path would have chosen instead. The anchor cannot be
  trusted to enumerate roots; it is trusted to land on one, and the far side is checked.
  The scan-granularity caveat of `EOS.md` §2 applies equally to both paths and is
  not addressed here.
- `DensityRoot::stable_root_count` becomes `std::optional<std::size_t>`. The cold path sets
  it; the anchored path leaves it empty. Nothing outside `density.cpp` reads it (verified by
  grep). No sentinel value.

Expected cost (reported, not asserted): a liquid anchor sits in the top decade of a grid that
spans about ten decades at four points per decade, so the far-side check is two to five
residual-only evaluations; a vapour anchor is symmetric at the dilute end.

### 2.6 The packets declare the convention

Eight new packet revisions, each a byte-for-byte copy of its predecessor with one edit to
`parameters/parameters.json`, rehashed in `packet.toml` and `provenance/*.toml`, and one new
entry in `data/packet-locks.json`. The key goes on the electrolyte record of `model_families`:

```json
"kind": "electrolyte", "choice": "born", "ionic_dispersion": "exclude-same-sign-ion-pairs", ...
```

with the record's existing `provenance.locator` extended by the sentence that states the
convention. Every source in the tree is the Held lineage; none is Ji 2012:

| packet | new version | declares | also dropped from `pairs` | locator for the declaration |
|---|---|---|---|---|
| `figiel-2025-reference-electrolytes/2` | `3` | `exclude-same-sign-ion-pairs` | nothing (four same-sign pairs absent) | Figiel, Yu & Held 2025, §2: "the dispersion energy between cation-cation pairs and anion-anion pairs is set to zero (u_ii = 0)" |
| `figiel-2025-validation-parameters/2` | `3` | same | 6 same-sign `k_ij = 0.0` records | same |
| `figiel-2025-aqueous-iodides/2` | `3` | same | nothing | same |
| `figiel-yu-held-2025/2` | `3` | same | nothing | same |
| `held-cameretti-sadowski-2008/2` | `3` | same | 198 same-sign `k_ij = 1.0` records; the 208 cation–anion `1.0` records **stay** (Held 2008 neglects all ion–ion dispersion; on cross pairs that is data the kernel reads) | Held, Cameretti & Sadowski 2008: "Dispersion between two ions is neglected in this work" |
| `ascani-2022-case-study-2-electrolyte-lle/2` | `3` | same | nothing | ePC-SAFT advanced (Bülow, Ascani & Held 2021) retains the 2014 combining rules; the like-ion rule is unchanged since Held 2014 |
| `khudaida-2026-figure-2-electrolyte-lle/5` | `6` | same | nothing (two ions; the diagonal is the whole exclusion) | Khudaida et al. 2026 use ePC-SAFT with cross-pair `k_ij` only |
| `mea-co2-h2o-nine-species-estimation/2` | `3` | same | 7 same-sign `k_ij = 0.0` records | estimated under the superseded engine, whose loader defaulted to the exclusion (`main@9eccbfd:_dictionary_io.py:750`); the bundle's values are only valid under it |

`mea-co2-h2o-nine-species-estimation/1` is historical and gets no revision; it remains
loadable nowhere in the greenfield tree, which is the fail-closed outcome for a set that does
not declare its convention. Provenance kind stays `repository-migration-v1` with
`copy_mode = "strict-schema-translation"`: the values are untouched and the added key
transcribes a sentence of the same source. If Stage P lands first, its lighter form carries the
same key and this table is its migration input.

Update every test and script reference: `figiel-2025-reference-electrolytes` 2 → 3 and
`mea-co2-h2o-nine-species-estimation` 2 → 3 in `tests/eos/*.py`, `tests/oracle/measure_floor.py`
and `tests/oracle/capture2.py`. `tests/packets.py` is unchanged.

### 2.7 The published-model row set

Add to `engine/tests/oracle/oracle.json` a case `aqueous_lithium_bromide_published_model`:
the 27 rows of `analyses/reproduction/2025-figiel/figures/figure-05/source/libr-model-digitized.csv`
(`molality_mol_kg`, `published_mean_ionic_activity_coefficient`, `digitization_uncertainty`),
`source_file` and its sha256 `a4fc73b54ed77598cc5ff2f30d5eeff93daaf8cf274d0847157cb5f21a0e5ed9`,
temperature 298.15 K, pressure 1e5 Pa, `kind: "published-model-digitized"`, parameters packet
`figiel-2025-reference-electrolytes/3`. This is a published-model reproduction in the sense
`TRANSFER.md` §"Oracle" already admits for Baygi 2015; it is not measured data and the row
says so.

## 3. Checks

`tests/native/physics_checks.cpp` — one new function, `same_sign_exclusion_is_a_limit`:
a two-component set (a neutral, one cation) with `ExcludeSameSign`; at `x = (0, 1)` the
dispersion term of `residual_helmholtz` is exactly `0.0` (a fluid of one ion species has no
dispersion under the rule); with `AllIonPairs` on the same set it is `< 0.0`; and at
`x = (0.5, 0.5)` the two conventions differ. Debye–Hückel and Born are passed as `0.0` as the
existing checks do.

`tests/native/model_checks.cpp` — three additions:

- in `hessian_fd_and_pressure_slope`, for every case, `evaluate(..., PressureSlope)
  .volume_second_derivative` against the centred difference of
  `evaluate(..., Gradient).gradient[1]` in `V`, to `1e-6` — the same construction the Hessian
  check already uses, on the order the density solve consumes.
- `same_sign_pair_is_inert`: on `nine_species()` (components 3 and 7 are both `+1`), setting
  `k_ij[3 * 9 + 7] = k_ij[7 * 9 + 3] = 0.5` leaves `state(...).dispersion` unchanged to
  `1e-14` under `ExcludeSameSign` and changes it under `AllIonPairs`.
- in `parameter_derivative_and_helpers`: `aqueous_electrolyte()` with
  `ionic_dispersion = Unspecified` throws `std::invalid_argument`; `neutral_binary()` with
  `ExcludeSameSign` throws.
- in `density_and_reference`: on ethane at 280 K, 1e6 Pa, `state_at_pressure(..., Liquid,
  anchor = vapour root)` returns the cold liquid root to `1e-10` relative, and
  `(..., Vapor, anchor = liquid root)` returns the cold vapour root; on `associating_binary()`
  at 313.15 K, 1e5 Pa, `Liquid` from cold converges (this is a case that throws today).

`tests/eos/test_activity.py` — `test_published_model_curve_reproduces_source`: through
`Mixture.reference(..., stoichiometry=(1, 1))` and `.activity(m)` on the 27 rows, AARD of
`mean_ionic_molality_scale` against the published value is at most the mean of
`digitization_uncertainty / published` over the same rows. No literal: the tolerance is the
source's own. Measured on the superseded engine, which carries the rule: AARD 0.27 %, max
0.56 %, against a mean digitization uncertainty of 0.6 %; the same engine with
`all-ion-pairs` gives 49.19 %. Print each row and the per-point ratio to its uncertainty.

`tests/eos/test_state.py` — extend `test_mapping_fails_closed_by_name`: loading
`figiel-2025-reference-electrolytes/2` raises `ParameterError` matching
`model/electrolyte/ionic_dispersion`; a `from_mapping` of a deep copy of revision 3 with an
injected `pair/lithium-cation/sodium-cation/k_ij = 0.0` record and components
`("water", "lithium-cation", "sodium-cation", "chloride-anion")` raises `ParameterError`
matching `lithium-cation/sodium-cation` and `exclude-same-sign-ion-pairs`.

## 4. Reported, not asserted

- `tests/native/model_floor.cpp` and `measure_floor.py`: the nine-species fixed-pressure
  cold solve's evaluation count and wall, before and after (`EOS.md` §4 gives
  86 evaluations / 6.4 ms before; §5 predicts 4–8 polish evaluations after).
- The `aqueous_alkali_halide_activity` per-salt AARDs after the repair, beside the
  superseded engine's 2.4–4.1 %.
- The count of far-side evaluations spent by anchored solves in `measure_floor.py`'s
  anchored case.

## 5. What this brief does not do

- `EOS.md` §2 (two roots inside one scan interval at ethane 280 K 1e4 Pa) is a
  separate policy decision on `DENSITY_BRACKET_INTERVALS` with its own measurement.
- The `k_hb_ij` cross-association correction (`cross_association_energy_correction`, registered,
  no successor) stays unsupported; no packet carries the family as a coefficient.
- Per-ion Zuber suppression coefficients are Stage 6's re-regression
  (`EOS.md` §3).

---

## Re-base delta — read this before the body, and prefer it where they differ

This brief was written against the tree before Stage 2. Stage 2 has since landed the core:
four headers moved from `physics/` to `core/`, the parameter builder moved from `model/` to
`surface/`, and every solver now returns a `core::Convergence`. The paths, sources,
preconditions, gates and ceilings above are already re-based and preflight clean. **What
follows is the design delta**, written by the core planner as part of Stage 2 and applied
here verbatim. Where it contradicts the body above, it wins.

- `base`: the Stage 2 merge sha. Every path to the four moved headers and the builder changes
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

Two further facts about the tree you are working on, from the Stage 2 review:

- `density.cpp`'s `association_failure` now asks *who threw* rather than *how it stopped* —
  `failure->record.criterion == &policy::ASSOCIATION_RESIDUAL`. Do not narrow it back to an
  enumeration of stops; the previous enumeration silently stopped recognising `stalled` when
  `classify` gained it, and the scan propagated instead of skipping a point.
- The anchored loop's budget exit currently names `stalled` by hand and keeps no
  `ResidualHistory`, against `CORE.md` §5.2's rule that a solver never names a residual stop
  itself. You are touching that block; fix it while you are there — build the history and let
  `classify` name the stop.


# Generic Relative-Permittivity Formulation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace provenance-bearing relative-permittivity runtime names with one explicit immutable generic formulation record per bundle while preserving every admitted equation and numerical result.

**Architecture:** `ModelParameterRecord(family="relative_permittivity_formulation", value=...)` becomes the sole formulation selector. Python resolves exact generic family coverage into a generic `_ResolvedModel`; the Python/native boundary maps the seven generic values to unchanged integer codes; native C++ retains the equations under generic symbols. Provenance metadata remains data-only under `src/epcsaft/data/bundles/`, enforced by a static repository check.

**Tech Stack:** Python 3.11+, frozen dataclasses, CSV/TOML bundle resources, Pint, C++17, pybind11, CppAD, pytest, Ruff, mypy, CMake/scikit-build-core.

## Global Constraints

- Native C++ remains the sole production owner of nonlinear relative-permittivity equations.
- CppAD remains the sole nonlinear production derivative authority.
- Preserve native formulation integer codes `0` through `6` and every scientific equation, parameter value, unit, domain check, and numerical tolerance.
- Require exactly one immutable `relative_permittivity_formulation` record in every installed and caller-provided bundle; neutral and associating bundles use `none`.
- Route exclusively on generic formulation and parameter-family names; never inspect bundle IDs, source IDs, locators, citations, authors, publication years, or domains to select equations.
- Remove `EPCSAFT.permittivity_model` with no compatibility alias and expose `EPCSAFT.relative_permittivity_formulation`.
- Do not accept or translate legacy formulation, family, topology, serialized-key, native-symbol, or EqID names.
- Do not add a mutable scientific dictionary, backend selector, hidden default, generic parameter overlay, duplicate equation owner, or new runtime dependency.
- Observe every focused test failing for the intended reason before changing its production owner.

---

## File Structure

- `src/epcsaft/records.py`: authoritative literal families, formulation values, units, and immutable record validation.
- `src/epcsaft/_resolve.py`: sole bundle-to-generic-resolved-model topology and family-coverage validation.
- `src/epcsaft/eos.py`: public property, generic-to-native integer mapping, native argument/export keys, and structural staged admission.
- `src/epcsaft/_native/model.hpp`: stable integer constants and generic native parameter storage.
- `src/epcsaft/_native/electrolyte.hpp`: unchanged seven relative-permittivity equations under generic symbols and EqIDs.
- `src/epcsaft/_native/eos.cpp`, `module.cpp`, `native_model.cpp`, `born.hpp`, `helmholtz.hpp`, `include/epcsaft/native_sdk_v1.h`: generic native validation, bindings, serialization, comments, and ionic-diameter EqIDs.
- `src/epcsaft/data/bundles/*/1/model.csv`: one explicit formulation record per installed bundle.
- `src/epcsaft/data/bundles/provenance_names.toml`: data-only prohibited runtime-name tokens.
- `scripts/check_runtime_provenance_names.py`: deterministic runtime-source policy checker.
- `docs/latex/equations.tex`: authored generic EqIDs; generated views remain derived.
- `docs/equations.md`, `docs/equations_registry.yaml`: regenerated equation views.
- `docs/python-api.md`, `docs/science/full-provider-eos.md`: generic public API and scientific terminology.

---

### Task 1: Make the formulation an explicit immutable bundle record

**Files:**
- Modify: `tests/test_records.py`
- Modify: `tests/test_parameters.py`
- Modify: `src/epcsaft/records.py`
- Create or modify: `src/epcsaft/data/bundles/*/1/model.csv`
- Create or modify: `tests/data/bundles/figiel-2025-contract/1/model.csv`

**Interfaces:**
- Produces: `RelativePermittivityFormulation` literal values and `ModelParameterRecord(family="relative_permittivity_formulation", value=<generic value>)`.
- Consumes: existing immutable `ModelParameterRecord` and CSV bundle loader.

- [ ] **Step 1: Write failing record tests**

Add a parametrized test containing exactly:

```python
FORMULATIONS = (
    "none",
    "ion-fraction-suppression",
    "component-permittivity-mixing",
    "density-volume-temperature",
    "ion-specific-suppression",
    "dipolar-density-perturbation",
    "solvent-only",
)


@pytest.mark.parametrize("value", FORMULATIONS)
def test_relative_permittivity_formulation_record_accepts_generic_values(value: str) -> None:
    record = ModelParameterRecord(
        "permittivity-formulation",
        "relative_permittivity_formulation",
        value,
        **PROVENANCE,
    )
    assert record.value == value


@pytest.mark.parametrize(
    "value",
    (
        "figiel-composition",
        "ascani-mole-fraction",
        "schreckenberg-density",
        "zuber-ion-specific",
        "rueben-perturbation",
        "held-water",
    ),
)
def test_relative_permittivity_formulation_rejects_legacy_values(value: str) -> None:
    with pytest.raises(ValueError, match="invalid-value"):
        ModelParameterRecord(
            "permittivity-formulation",
            "relative_permittivity_formulation",
            value,
            **PROVENANCE,
        )
```

Add catalog tests asserting every discovered installed bundle has exactly one
model record with family `relative_permittivity_formulation`, including
`none` for both neutral bundles.

- [ ] **Step 2: Run the focused tests and observe the schema failure**

Run:

```bash
python -m pytest tests/test_records.py tests/test_parameters.py -k 'relative_permittivity_formulation or installed_bundle' -q
```

Expected: FAIL because the family/value literals and required bundle records do not exist.

- [ ] **Step 3: Implement the generic record and family literals**

In `records.py`, define:

```python
RelativePermittivityFormulation: TypeAlias = Literal[
    "none",
    "ion-fraction-suppression",
    "component-permittivity-mixing",
    "density-volume-temperature",
    "ion-specific-suppression",
    "dipolar-density-perturbation",
    "solvent-only",
]
```

Add `relative_permittivity_formulation` to `ModelParameterFamily`, validate it
against `frozenset(get_args(RelativePermittivityFormulation))`, and replace the
seven provenance-bearing parameter families with:

```python
"ion_fraction_suppression_coefficient"
"dielectric_volume"
"dielectric_temperature"
"ion_specific_suppression_coefficient"
"dipole_scaling"
"polarizability_scaling"
"correlation_integral_parameter"
```

Rename the existing unit and dimensionless-family entries without changing
their units or validation.

- [ ] **Step 4: Migrate every bundle record**

Create `model.csv` where absent with the existing model CSV columns. Add
exactly one formulation record to every installed and test-data bundle.
Use `none` for neutral/associating bundles,
`component-permittivity-mixing` for the original-Born LLE bundle,
`solvent-only` for the fully dissociated water bundle, and
`ion-fraction-suppression` for the currently installed corrected-Born aqueous
bundles. Preserve source, locator, and domain fields verbatim. Rename all
formulation-specific family values in `single.csv` and `model.csv`.

- [ ] **Step 5: Run schema and catalog tests**

Run:

```bash
python -m pytest tests/test_records.py tests/test_parameters.py -q
```

Expected: PASS.

- [ ] **Step 6: Commit the schema slice**

```bash
git add src/epcsaft/records.py src/epcsaft/data/bundles tests/data/bundles tests/test_records.py tests/test_parameters.py
git commit -m "refactor: add explicit permittivity formulation records"
```

### Task 2: Resolve only generic formulations and expose the new public property

**Files:**
- Modify: `tests/test_model_resolution.py`
- Modify: `tests/test_eos.py`
- Modify: `tests/test_permittivity_models.py`
- Modify: `tests/typecheck/public_api.py`
- Modify: `tests/data/native_models/gross-methane-ethane.json`
- Modify: `src/epcsaft/_resolve.py`
- Modify: `src/epcsaft/eos.py`

**Interfaces:**
- Consumes: Task 1 formulation and generic family literals.
- Produces: `_ResolvedModel.relative_permittivity_formulation`, generic resolved arrays, `EPCSAFT.relative_permittivity_formulation`, and generic native-export keys.

- [ ] **Step 1: Rewrite fixtures and add failing resolver/public-API tests**

Change the independent oracle fixture to always append:

```python
ModelParameterRecord(
    "permittivity-formulation",
    "relative_permittivity_formulation",
    kind,
    **PROVENANCE,
)
```

Use only the seven generic values and generic family names. Assert:

```python
assert model.relative_permittivity_formulation == kind
assert not hasattr(model, "permittivity_model")
```

Add negative tests for a missing formulation record, two formulation records,
an extraneous formulation-specific family, and source/locator changes that
leave resolution unchanged.

- [ ] **Step 2: Run resolver tests and observe legacy inference/API failures**

Run:

```bash
python -m pytest tests/test_model_resolution.py tests/test_eos.py tests/test_permittivity_models.py -q
```

Expected: FAIL because `_resolve_model` still infers author-bearing values and
the public property/export keys are legacy names.

- [ ] **Step 3: Replace inference with exact formulation dispatch**

In `_ResolvedModel`, use:

```python
relative_permittivity_formulation: RelativePermittivityFormulation
dielectric_volume_m3_per_mol: tuple[float, ...] = ()
dielectric_temperature_k: tuple[float, ...] = ()
ion_specific_suppression_coefficient: tuple[float, ...] = ()
dipole_scaling_debye2: tuple[float, ...] = ()
polarizability_scaling_angstrom3: tuple[float, ...] = ()
correlation_integral_parameter: tuple[float, ...] = ()
```

Read exactly one `relative_permittivity_formulation` model record before the
neutral/ionic split. Dispatch on that value and validate required/extraneous
generic families explicitly. Use generic topology values:

```python
"ionic-original-born"
"ionic-fully-dissociated-debye-huckel"
"ionic-ssm-ds"
```

Rename `electrolyte_formulation` value to
`fully-dissociated-debye-huckel`, including its exact component/charge and
`k_ij` validation. Remove all family-presence formulation inference.

- [ ] **Step 4: Replace the Python public and serialization boundary**

Define:

```python
@property
def relative_permittivity_formulation(self) -> str:
    return self._model.relative_permittivity_formulation
```

Delete `permittivity_model`. Rename Python-side argument variables and
`export_native_model` keys to the generic resolved names. Map formulations to
the unchanged codes:

```python
{
    "none": 0,
    "ion-fraction-suppression": 1,
    "component-permittivity-mixing": 2,
    "density-volume-temperature": 3,
    "ion-specific-suppression": 4,
    "dipolar-density-perturbation": 5,
    "solvent-only": 6,
}
```

- [ ] **Step 5: Run Python-side focused tests**

Run:

```bash
python -m pytest tests/test_records.py tests/test_parameters.py tests/test_model_resolution.py tests/test_eos.py tests/test_permittivity_models.py -q
python -m mypy --no-incremental
```

Expected: resolver tests pass; tests constructing the old native signature may
still fail until Task 3 and must fail only at that boundary.

- [ ] **Step 6: Commit the resolver/API slice**

```bash
git add src/epcsaft/_resolve.py src/epcsaft/eos.py tests/test_model_resolution.py tests/test_eos.py tests/test_permittivity_models.py tests/typecheck/public_api.py tests/data/native_models/gross-methane-ethane.json
git commit -m "refactor: resolve generic permittivity formulations"
```

### Task 3: Rename the native boundary without changing equations

**Files:**
- Modify: `src/epcsaft/_native/model.hpp`
- Modify: `src/epcsaft/_native/electrolyte.hpp`
- Modify: `src/epcsaft/_native/eos.cpp`
- Modify: `src/epcsaft/_native/module.cpp`
- Modify: `src/epcsaft/_native/native_model.cpp`
- Modify: `src/epcsaft/_native/born.hpp`
- Modify: `src/epcsaft/_native/helmholtz.hpp`
- Modify: `src/epcsaft/include/epcsaft/native_sdk_v1.h`
- Modify: `tests/test_native_sdk.py`
- Modify: `tests/test_permittivity_models.py`

**Interfaces:**
- Consumes: Task 2 generic Python argument order and export keys.
- Produces: stable native codes under `RELATIVE_PERMITTIVITY_FORMULATION_*` constants and unchanged scalar-generic numerical behavior.

- [ ] **Step 1: Add failing native-code and oracle assertions**

Assert the exported native JSON uses
`relative_permittivity_formulation` and the six generic parameter-array keys,
and contains integer codes `0..6` for the seven generic fixtures. Retain the
independent equation expected values and pressure/CppAD derivative checks.

- [ ] **Step 2: Run tests and observe native signature/key failures**

Run:

```bash
python -m pytest tests/test_native_sdk.py tests/test_permittivity_models.py -q
```

Expected: FAIL at legacy native constructor fields, constants, and JSON keys.

- [ ] **Step 3: Rename native storage and stable constants**

Rename code constants while retaining values:

```cpp
inline constexpr int RELATIVE_PERMITTIVITY_FORMULATION_NONE = 0;
inline constexpr int RELATIVE_PERMITTIVITY_FORMULATION_ION_FRACTION_SUPPRESSION = 1;
inline constexpr int RELATIVE_PERMITTIVITY_FORMULATION_COMPONENT_MIXING = 2;
inline constexpr int RELATIVE_PERMITTIVITY_FORMULATION_DENSITY_VOLUME_TEMPERATURE = 3;
inline constexpr int RELATIVE_PERMITTIVITY_FORMULATION_ION_SPECIFIC_SUPPRESSION = 4;
inline constexpr int RELATIVE_PERMITTIVITY_FORMULATION_DIPOLAR_DENSITY_PERTURBATION = 5;
inline constexpr int RELATIVE_PERMITTIVITY_FORMULATION_SOLVENT_ONLY = 6;
```

Rename every corresponding struct field, pybind argument, JSON key,
diagnostic, local variable, comment, and C API comment generically. Preserve
constructor order, vector dimensions, inequalities, and arithmetic.

- [ ] **Step 4: Rename equation-owner locals and EqIDs**

Apply the ten EqID mappings from the approved design. In
`electrolyte.hpp`, rename only identifiers and comments; compare the
pre/post arithmetic expressions to ensure each coefficient, operation, and
branch is unchanged.

- [ ] **Step 5: Rebuild the editable native module**

Run:

```bash
python -m pip install -e .
```

Expected: successful C++17 build and editable installation.

- [ ] **Step 6: Run native and numerical gates**

Run:

```bash
python -m pytest tests/test_permittivity_models.py tests/test_native_sdk.py tests/test_born.py tests/test_debye_huckel.py tests/test_density.py -q
```

Expected: PASS at existing numerical tolerances.

- [ ] **Step 7: Commit the native rename**

```bash
git add src/epcsaft/_native src/epcsaft/include/epcsaft/native_sdk_v1.h tests/test_native_sdk.py tests/test_permittivity_models.py
git commit -m "refactor: make native permittivity names generic"
```

### Task 4: Make staged aqueous admission provenance-independent

**Files:**
- Modify: `src/epcsaft/eos.py`
- Modify: `src/epcsaft/reference.py`
- Modify: `tests/test_eos.py`
- Modify: `tests/test_native_sdk.py`
- Modify: the existing staged callback test modules selected by `rg -l '_figiel_staged|staged_trial_admitted|active_born' tests`

**Interfaces:**
- Consumes: Task 2 generic resolved formulation/topology/families.
- Produces: `_staged_aqueous_trial_bundle_admitted(parameters, model) -> bool`, based only on generic structural facts.

- [ ] **Step 1: Add failing provenance-invariance tests**

Clone a valid staged bundle through `ParameterBundle.from_records`, replacing
every `source_id`, `locator`, and `domain_id` consistently while retaining
components, topology, formulation, families, values, and units. Assert the
same callback is admitted and returns the same value/derivative. Add a
structurally incomplete clone and assert rejection.

- [ ] **Step 2: Run staged callback tests and observe source-ID dependence**

Run:

```bash
python -m pytest tests -k 'staged and (aqueous or born or kij or solvation)' -q
```

Expected: the renamed-provenance valid bundle FAILS admission.

- [ ] **Step 3: Replace catalog/source comparisons with structural checks**

Rename the admission helper and stored boolean generically. Remove catalog
loading and all source/domain string comparisons. Check only exact component
IDs and charges required by the callback, generic topology, generic
formulation, complete generic family coverage, active parameter ordering, and
the existing finite/unit/domain-value constraints.

Rename author-bearing private constants, diagnostics, and comments in
`eos.py` and `reference.py`.

- [ ] **Step 4: Run callback and integration tests**

Run:

```bash
python -m pytest tests -k 'staged or active_born or active_kij or solvation_factor' -q
```

Expected: PASS, including provenance-invariance and incomplete-structure rejection.

- [ ] **Step 5: Commit the admission slice**

```bash
git add src/epcsaft/eos.py src/epcsaft/reference.py tests
git commit -m "refactor: admit staged aqueous bundles structurally"
```

### Task 5: Enforce the runtime provenance-name boundary

**Files:**
- Create: `src/epcsaft/data/bundles/provenance_names.toml`
- Create: `scripts/check_runtime_provenance_names.py`
- Create: `tests/test_runtime_provenance_names.py`
- Modify: any remaining file under `src/epcsaft/` reported by the new check

**Interfaces:**
- Produces: `scan_runtime_source(repo_root: Path) -> tuple[Violation, ...]` and CLI exit `0` for clean, `1` for violations.
- Consumes: TOML `tokens = [...]` data inside the bundle provenance boundary.

- [ ] **Step 1: Write failing checker tests**

Use a temporary repository tree and assert:

```python
def test_scan_rejects_provenance_name_in_runtime_source(tmp_path: Path) -> None:
    # token policy contains "exampleauthor"
    # src/epcsaft/runtime.py contains "ExampleAuthor"
    assert scan_runtime_source(tmp_path)[0].line == 1


def test_scan_excludes_bundle_provenance(tmp_path: Path) -> None:
    # the same token under src/epcsaft/data/bundles/example/1/bundle.toml
    assert scan_runtime_source(tmp_path) == ()
```

Also assert CLI diagnostics contain the relative file, line, and token.

- [ ] **Step 2: Run checker tests and observe missing module failure**

Run:

```bash
python -m pytest tests/test_runtime_provenance_names.py -q
```

Expected: FAIL because the checker does not exist.

- [ ] **Step 3: Implement the deterministic static checker**

Load `provenance_names.toml` with `tomllib`, recursively scan regular text
files below `src/epcsaft`, skip the entire `data/bundles` subtree, and match
case-insensitive identifier tokens with boundaries equivalent to
`(?<![A-Za-z0-9])TOKEN(?![A-Za-z0-9])`. Return sorted immutable violations
and print `path:line: prohibited provenance name: token`.

Keep every real author token only in `provenance_names.toml`; production and
checker source use generic fixture tokens.

- [ ] **Step 4: Remove every reported runtime occurrence**

Run:

```bash
python scripts/check_runtime_provenance_names.py
```

Rename remaining source identifiers, comments, diagnostics, locals, and
strings generically until exit status is zero. Do not edit bundle provenance
to make the check pass.

- [ ] **Step 5: Run checker tests, checker, lint, and type checking**

Run:

```bash
python -m pytest tests/test_runtime_provenance_names.py -q
python scripts/check_runtime_provenance_names.py
python -m ruff check src scripts tests
python -m mypy --no-incremental
```

Expected: PASS.

- [ ] **Step 6: Commit the policy gate**

```bash
git add src/epcsaft scripts/check_runtime_provenance_names.py tests/test_runtime_provenance_names.py
git commit -m "test: enforce generic runtime provenance names"
```

### Task 6: Synchronize equation/docs surfaces and run the complete gate

**Files:**
- Modify: `docs/latex/equations.tex`
- Regenerate: `docs/equations.md`
- Regenerate: `docs/equations_registry.yaml`
- Modify: `docs/python-api.md`
- Modify: `docs/science/full-provider-eos.md`
- Modify: `tests/test_equation_registry.py`
- Modify: any tests still containing legacy public/runtime identifiers

**Interfaces:**
- Consumes: Tasks 1-5 completed generic runtime surface.
- Produces: current generic documentation, equation registry, and complete verified repository.

- [ ] **Step 1: Add failing registry and negative-space assertions**

Update registry expectations to the ten generic EqIDs. Add a repository
negative-space assertion that the legacy public property, formulation values,
family values, topology values, native keys, and EqIDs do not occur in
runtime source or generated registry.

- [ ] **Step 2: Run registry tests and observe stale authored/generated views**

Run:

```bash
python -m pytest tests/test_equation_registry.py -q
python scripts/docs/sync_equation_registry.py --check
```

Expected: FAIL because LaTeX and generated views still use legacy EqIDs.

- [ ] **Step 3: Rename authored EqIDs and update public/science docs**

Apply the approved ten EqID mappings in `docs/latex/equations.tex`; keep
source citations, equations, units, and locators unchanged. Document the
seven generic public formulation values and
`EPCSAFT.relative_permittivity_formulation` in `docs/python-api.md` and
`docs/science/full-provider-eos.md`.

- [ ] **Step 4: Regenerate and verify equation views**

Run:

```bash
python scripts/docs/sync_equation_registry.py
python scripts/docs/sync_equation_registry.py --check
python -m pytest tests/test_equation_registry.py -q
```

Expected: PASS.

- [ ] **Step 5: Run complete verification**

Run:

```bash
python -m pip install -e .
python scripts/check_runtime_provenance_names.py
python scripts/docs/sync_equation_registry.py --check
python -m pytest -q
python -m ruff check .
python -m mypy --no-incremental
git diff --check
bash ~/.codex/hooks/codex-cleanup.sh --repo-root .
```

Expected: all commands pass. Cleanup candidates not created by this task are
reported but not removed.

- [ ] **Step 6: Audit exact negative space and scientific preservation**

Run:

```bash
rg -n -i 'figiel|ascani|schreckenberg|zuber|rueben|held' src/epcsaft --glob '!data/bundles/**'
rg -n 'permittivity_model|figiel-composition|ascani-mole-fraction|schreckenberg-density|zuber-ion-specific|rueben-perturbation|held-water' src tests docs/equations.md docs/equations_registry.yaml
git diff fea1dc4 -- src/epcsaft/_native/electrolyte.hpp
```

Expected: the first two searches return no runtime/generated identifier
matches; the native diff shows naming and EqID changes only, with arithmetic
expressions preserved.

- [ ] **Step 7: Commit documentation and final verification state**

```bash
git add docs tests src/epcsaft/data/bundles scripts
git commit -m "docs: publish generic permittivity formulation surface"
```

- [ ] **Step 8: Verify final branch state**

Run:

```bash
git status --short --branch
git log --oneline --decorate -8
```

Expected: clean `codex/generic-relative-permittivity-formulation` with the
design, plan, and implementation commits ahead of its starting point.

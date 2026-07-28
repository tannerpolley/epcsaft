# Provider Equation Documentation Mirror Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a complete, source-traceable, mechanically synchronized LaTeX mirror of every meaningful mathematical relation implemented by the current Provider.

**Architecture:** Current C++ and Python implementation files remain the equation owners. Stable `EqID` comments link those owners to `docs/latex/equations.tex`; a dependency-free Python synchronizer validates complete relation/exclusion classification and generates searchable Markdown and YAML indexes. CppAD remains the sole nonlinear derivative authority, so the mirror documents defining equations and exact property transformations but contains no expanded handwritten nonlinear derivative catalogue.

**Tech Stack:** C++17, Python 3.11+, CppAD, pytest, standard-library Python parsing/rendering, LaTeX/biblatex for optional local rendering.

## Global Constraints

- Work only on `codex/provider-equation-documentation-mirror` in the attached Provider checkout.
- Preserve every current runtime equation, parameter, numerical algorithm, public API, native ABI, and package dependency.
- Treat Migration-pinned legacy `docs/latex/equations.tex` at commit `7920e56dbb0fc1fb418bae44e67ef50f334b75cd`, Git blob `2e970476258558de05391d9d77f31d1736d78de8`, as evidence and structural guidance only.
- Reconstruct the document from current Provider owners; do not bulk-copy the 172 legacy equations or generated registries.
- Track every production arithmetic relation that defines a thermodynamic, model, property, coordinate, conversion, limiting, polynomial, mixing, or source-backed helper quantity.
- Group multiple implementation statements under one `EqID` only when they evaluate one displayed mathematical relation.
- Classify every function in declared mathematical-owner files with at least one `EqID` or an explicit `EqExclude` reason.
- Allow `EqExclude` only for scalar dispatch, numerical iteration mechanics, validation, transport/buffer mechanics, or serialization.
- Add stable implementation comments only; do not modify C++ or Python function bodies.
- Do not include expanded analytical nonlinear gradients, Hessians, or higher derivatives. Document that CppAD differentiates the defining equations.
- Track exact thermodynamic identities and exact algebraic/implicit transformations that Provider actually consumes.
- Track `docs/equations.md` and `docs/equations_registry.yaml`; do not track a generated PDF.
- Add no runtime or wheel dependency on LaTeX, YAML libraries, the legacy repository, Migration, or the lab.
- Keep `docs/science/full-provider-eos.md` as the normative scientific contract and `docs/latex/equations.tex` as its exhaustive checked mirror.

---

## File Structure

**Create**

- `scripts/docs/sync_equation_registry.py`: parse equation metadata and implementation annotations, enforce traceability, and render generated views.
- `tests/test_equation_registry.py`: focused parser, linkage, exclusion, metadata, citation, coverage, and freshness tests.
- `docs/latex/equations.tex`: current Provider equation catalogue.
- `docs/latex/references.bib`: bibliography entries used by the catalogue.
- `docs/latex/README.md`: provenance, authority boundary, synchronization, and optional render commands.
- `docs/equations_registry.yaml`: generated machine-readable equation and exclusion index.
- `docs/equations.md`: generated searchable Markdown equation index.

**Modify with comments only**

- `src/epcsaft/_native/helmholtz.hpp`
- `src/epcsaft/_native/association.hpp`
- `src/epcsaft/_native/electrolyte.hpp`
- `src/epcsaft/_native/born.hpp`
- `src/epcsaft/_native/derivatives.hpp`
- `src/epcsaft/_native/density.hpp`
- `src/epcsaft/_native/electrolyte_reference.hpp`
- `src/epcsaft/_native/eos.cpp`
- `src/epcsaft/_native/phase_block.hpp`
- `src/epcsaft/_native/native_sdk.cpp`
- `src/epcsaft/reference.py`

**Modify as documentation**

- `docs/science/full-provider-eos.md`
- `README.md`

---

### Task 1: Build the synchronization and classification engine

**Files:**

- Create: `scripts/docs/sync_equation_registry.py`
- Create: `tests/test_equation_registry.py`

**Interfaces:**

- Consumes: LaTeX metadata comments, `// EqID:` and `// EqExclude:` C++ comments, and `# EqID:` / `# EqExclude:` Python comments.
- Produces:
  - `EquationEntry`
  - `SourceReference`
  - `SourceExclusion`
  - `parse_equations(tex_path: Path, repo_root: Path) -> list[EquationEntry]`
  - `parse_source_annotations(source_root: Path, repo_root: Path) -> tuple[dict[str, list[SourceReference]], list[SourceExclusion]]`
  - `load_repository_inputs(repo_root: Path) -> tuple[list[EquationEntry], dict[str, list[SourceReference]], list[SourceExclusion]]`
  - `load_bibliography_keys(repo_root: Path) -> set[str]`
  - `unclassified_owner_functions(repo_root: Path) -> list[str]`
  - `validate_registry(entries, references, exclusions, bibliography_keys) -> None`
  - `render_registry(entries, exclusions) -> str`
  - `render_markdown(entries, exclusions) -> str`
  - `synchronize(repo_root: Path, *, check: bool) -> None`

- [ ] **Step 1: Write the failing metadata and linkage tests**

Add a loader and minimal fixtures:

```python
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "docs" / "sync_equation_registry.py"


def _load_sync_module():
    spec = importlib.util.spec_from_file_location("sync_equation_registry", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_equation_parser_requires_unique_complete_metadata(tmp_path: Path) -> None:
    module = _load_sync_module()
    tex = tmp_path / "equations.tex"
    tex.write_text(
        r"""
\section{Residual Helmholtz}
% EqID: residual_helmholtz
% Source: \cite{Figiel2025}, Eq.~(1)
% Status: Source equation
% Description: Sum of active residual contributions.
% Change note: Reconciled to current Provider.
\begin{equation}
a^r=a^{hc}+a^{disp}+a^{assoc}+a^{DH}+a^{Born}
\label{eq:residual_helmholtz}
\end{equation}
""",
        encoding="utf-8",
    )
    entries = module.parse_equations(tex, tmp_path)
    assert [entry.eqid for entry in entries] == ["residual_helmholtz"]
    assert entries[0].status == "Source equation"
```

Add a second fixture with a duplicate `EqID` and assert `ValueError` contains
`duplicate EqID`.

- [ ] **Step 2: Run the test and verify the missing-script failure**

Run:

```bash
python -m pytest tests/test_equation_registry.py::test_equation_parser_requires_unique_complete_metadata -q
```

Expected: FAIL because `scripts/docs/sync_equation_registry.py` does not exist.

- [ ] **Step 3: Implement immutable entries and LaTeX parsing**

Use standard-library dataclasses and fixed statuses:

```python
@dataclass(frozen=True)
class EquationEntry:
    eqid: str
    section: str
    subsection: str
    label: str
    source: str
    status: str
    description: str
    change_note: str
    latex: str
    tex_file: str
    tex_line: int
    implementation_refs: tuple["SourceReference", ...] = ()


ALLOWED_STATUSES = frozenset(
    {
        "Source equation",
        "Exact derivation",
        "Implementation convention",
        "Bounded empirical relation",
        "Documentation-only",
    }
)
```

Distill only the equation-block parser from the pinned legacy
`scripts/docs/sync_equation_registry.py`. Require all five metadata fields,
an equation label, a supported status, and an `EqID` matching
`[a-z][a-z0-9_]*`.

- [ ] **Step 4: Run the parser tests**

Run:

```bash
python -m pytest tests/test_equation_registry.py -q
```

Expected: parser tests PASS; later linkage tests are not yet present.

- [ ] **Step 5: Write failing implementation-link and exclusion tests**

Add fixtures proving:

```python
def test_source_annotations_link_real_owners_and_require_function_classification(
    tmp_path: Path,
) -> None:
    module = _load_sync_module()
    native = tmp_path / "src" / "epcsaft" / "_native"
    native.mkdir(parents=True)
    (native / "owner.hpp").write_text(
        """
// EqID: residual_helmholtz
double residual_helmholtz(double a, double b) {
    return a + b;
}

// EqExclude: scalar utility dispatch
double scalar_log(double value) {
    return std::log(value);
}
""",
        encoding="utf-8",
    )
    refs, exclusions = module.parse_source_annotations(tmp_path / "src", tmp_path)
    assert refs["residual_helmholtz"][0].symbol == "residual_helmholtz"
    assert exclusions[0].reason == "scalar utility dispatch"
```

Also assert that an unannotated function in a declared mathematical-owner file
fails with `unclassified mathematical-owner function`.

- [ ] **Step 6: Implement C++/Python annotation parsing and classification**

Support `.cpp`, `.hpp`, `.h`, and `.py`. Use Python `ast` for Python function
names and a bounded C++ function-signature recognizer matching the repository's
current free functions and `NativeModel::method` definitions. Associate
pending comments with the next function or arithmetic statement. Record:

```python
@dataclass(frozen=True)
class SourceReference:
    file: str
    line: int
    comment_line: int
    symbol: str
    context: str


@dataclass(frozen=True)
class SourceExclusion:
    file: str
    line: int
    symbol: str
    reason: str
```

Use this exact exclusion vocabulary:

```python
ALLOWED_EXCLUSIONS = frozenset(
    {
        "scalar utility dispatch",
        "numerical iteration mechanics",
        "input and domain validation",
        "transport and buffer mechanics",
        "serialization",
        "CppAD-generated derivative machinery",
    }
)
```

Implement `load_repository_inputs()` as the single owner of the declared
mathematical-owner file set, `load_bibliography_keys()` as the biblatex-key
reader, and `unclassified_owner_functions()` as the complete function-level
coverage audit used by the final gate. Do not duplicate the owner-file list in
tests or rendering code.

- [ ] **Step 7: Write failing rendering, citation, and freshness tests**

Add tests that:

- reject a source `EqID` absent from LaTeX;
- reject an implemented LaTeX entry without an implementation reference;
- exempt only `Status: Documentation-only`;
- reject a `\cite{MissingKey}` absent from `references.bib`;
- render deterministic YAML and Markdown;
- make `synchronize(..., check=True)` fail when either generated file differs.

- [ ] **Step 8: Implement validation and deterministic rendering**

Render YAML without adding PyYAML. Use JSON string quoting so the output remains
valid YAML. Include `implementation_refs` and a top-level
`source_exclusions` section. The Markdown view must include the formula,
metadata, LaTeX location, implementation owners, and audited exclusions.

Add CLI behavior:

```bash
python scripts/docs/sync_equation_registry.py
python scripts/docs/sync_equation_registry.py --check
```

The first writes generated views; the second performs no writes.

- [ ] **Step 9: Run focused tests and static checks**

Run:

```bash
python -m pytest tests/test_equation_registry.py -q
python -m ruff check scripts/docs/sync_equation_registry.py tests/test_equation_registry.py
git diff --check
```

Expected: all PASS.

- [ ] **Step 10: Commit the synchronization engine**

```bash
git add scripts/docs/sync_equation_registry.py tests/test_equation_registry.py
git commit -m "test: add equation traceability engine"
```

---

### Task 2: Document neutral PC-SAFT and association mathematics

**Files:**

- Create: `docs/latex/equations.tex`
- Create: `docs/latex/references.bib`
- Create: `docs/latex/README.md`
- Create: `docs/equations_registry.yaml`
- Create: `docs/equations.md`
- Modify with comments only: `src/epcsaft/_native/helmholtz.hpp`
- Modify with comments only: `src/epcsaft/_native/association.hpp`
- Test: `tests/test_equation_registry.py`

**Interfaces:**

- Consumes: Task 1 synchronizer.
- Produces: the initial compilable LaTeX document and linked neutral/association equation families.

- [ ] **Step 1: Add the failing real-repository coverage test**

```python
def test_repository_equation_registry_is_current_and_complete() -> None:
    module = _load_sync_module()
    module.synchronize(REPO_ROOT, check=True)
```

Run:

```bash
python -m pytest tests/test_equation_registry.py::test_repository_equation_registry_is_current_and_complete -q
```

Expected: FAIL because `docs/latex/equations.tex` and generated outputs are
absent.

- [ ] **Step 2: Create the LaTeX shell and provenance README**

Create a standalone article using `amsmath`, `amssymb`, `mathtools`,
`siunitx`, `biblatex`, and `hyperref`. Its introduction must state:

```text
This document mirrors the current Provider implementation. It is not a runtime
input and does not override docs/science/full-provider-eos.md. The legacy
catalogue at commit 7920e56d, blob 2e970476..., is provenance only.
CppAD is the sole nonlinear derivative authority.
```

Add bibliography entries for `Gross2001`, `Gross2002`, and `Chapman1990`.
Document these commands in `docs/latex/README.md`:

```bash
python scripts/docs/sync_equation_registry.py --check
cd docs/latex
latexmk -pdf -outdir=out equations.tex
```

State that `out/` and generated PDFs are untracked.

- [ ] **Step 3: Add neutral and association EqIDs beside the defining expressions**

Use these stable identifiers:

```text
angstrom_number_density
segment_sigma_correlation
neutral_segment_diameter
figiel_ionic_segment_diameter
ascani_ionic_segment_diameter
mean_segment_count
packing_moments
packing_fraction
close_packing_fraction
packing_fraction_coefficient
phase_packing_fraction
hard_sphere_helmholtz
hard_sphere_contact
hard_chain_helmholtz
dispersion_chain_coefficients
dispersion_polynomial
dispersion_integral_one
dispersion_integral_two
dispersion_compressibility_factor
cross_segment_diameter
cross_dispersion_energy
same_sign_ionic_dispersion
dispersion_first_moment
dispersion_second_moment
dispersion_helmholtz
association_contact
association_strength
association_mass_action
association_helmholtz
residual_helmholtz
```

Place comments directly before the corresponding current expressions.
`polynomial` owns `dispersion_polynomial`; `close_packing_fraction`,
`packing_fraction_coefficient`, and `phase_packing_fraction` likewise receive
their listed EqIDs. Mark only `scalar_log`, `scalar_exp`, `scalar_sqrt`,
`association_exp`, `solve_site_fractions`, and `solve_linear` with the
applicable exact `EqExclude` reason rather than inventing thermodynamic
equations for utility or iteration mechanics.

- [ ] **Step 4: Write the matching LaTeX sections**

Follow the legacy structure:

```text
Symbols, Bases, and Units
Parameter Setup
Density and Packing
Contribution Intermediates
  Hard Chain
  Dispersion
  Association
Residual Helmholtz Energy
```

Use the current formulas in `helmholtz.hpp` and `association.hpp`, including
the current Å-based number-density factor, current `l_ij` diameter rule,
same-sign ionic dispersion suppression, and current association-site
multiplicity. Do not import any derivative-suffixed legacy EqIDs.

- [ ] **Step 5: Generate the registry and Markdown view**

Run:

```bash
python scripts/docs/sync_equation_registry.py
python scripts/docs/sync_equation_registry.py --check
```

Expected: both commands succeed and every implemented entry has a source
reference.

- [ ] **Step 6: Run focused equation and neutral scientific tests**

```bash
python -m pytest tests/test_equation_registry.py tests/test_eos.py tests/test_association.py -q
python -m ruff check scripts/docs/sync_equation_registry.py tests/test_equation_registry.py
git diff --check
```

Expected: all PASS with unchanged numerical tests.

- [ ] **Step 7: Commit the neutral/association mirror**

```bash
git add docs/latex docs/equations.md docs/equations_registry.yaml \
  src/epcsaft/_native/helmholtz.hpp src/epcsaft/_native/association.hpp \
  tests/test_equation_registry.py
git commit -m "docs: mirror neutral and association equations"
```

---

### Task 3: Document electrolyte and Born mathematics

**Files:**

- Modify: `docs/latex/equations.tex`
- Modify: `docs/latex/references.bib`
- Modify with comments only: `src/epcsaft/_native/electrolyte.hpp`
- Modify with comments only: `src/epcsaft/_native/born.hpp`
- Regenerate: `docs/equations_registry.yaml`
- Regenerate: `docs/equations.md`
- Test: `tests/test_equation_registry.py`

**Interfaces:**

- Consumes: Task 2 equation document and synchronizer.
- Produces: complete admitted dielectric, Debye--Hückel, and Born equation coverage.

- [ ] **Step 1: Add a failing required-family test**

Add:

```python
def test_required_provider_equation_families_are_present() -> None:
    module = _load_sync_module()
    entries = module.parse_equations(
        REPO_ROOT / "docs" / "latex" / "equations.tex",
        REPO_ROOT,
    )
    eqids = {entry.eqid for entry in entries}
    assert {
        "figiel_bulk_permittivity",
        "ascani_bulk_permittivity",
        "debye_screening_parameter",
        "debye_huckel_helmholtz",
        "original_born_helmholtz",
        "ssm_ds_born_helmholtz",
    } <= eqids
```

Run the test and expect failure on the missing identifiers.

- [ ] **Step 2: Add current electrolyte EqIDs**

Use:

```text
salt_free_solvent_permittivity
figiel_bulk_permittivity
ascani_bulk_permittivity
ionic_strength_charge_sum
debye_screening_parameter
debye_finite_size_chi
debye_huckel_helmholtz
```

Track both branches of `finite_size_chi`: the exact logarithmic expression and
its small-argument series are one piecewise displayed relation. Classify
`electrolyte_sqrt` and `electrolyte_log` as scalar utility dispatch.

- [ ] **Step 3: Add current Born EqIDs**

Use:

```text
mixture_solvation_factor
born_shell_diameter
ssm_ds_radial_work
ssm_ds_born_helmholtz
original_born_radial_work
original_born_helmholtz
```

The active-Born and active-solvation-factor overloads reuse the same defining
EqIDs and receive additional implementation references; do not duplicate the
equations for trial coordinates.

- [ ] **Step 4: Add source-faithful LaTeX and bibliography entries**

Add only the current sources:

```text
Bulow2019
Figiel2025
Figiel2026Correction
Ascani2022
AscaniHeld2021
```

Represent the corrected SSM+DS bracket from Figiel 2026 and distinguish
Figiel's ionic-region permittivity from Ascani's ion-mixture permittivity.
Preserve diameter and SI conversion conventions exactly as implemented.

- [ ] **Step 5: Regenerate and run focused scientific checks**

```bash
python scripts/docs/sync_equation_registry.py
python scripts/docs/sync_equation_registry.py --check
python -m pytest tests/test_equation_registry.py tests/test_debye_huckel.py tests/test_born.py -q
git diff --check
```

Expected: all PASS.

- [ ] **Step 6: Commit the electrolyte/Born mirror**

```bash
git add docs/latex/equations.tex docs/latex/references.bib \
  docs/equations.md docs/equations_registry.yaml \
  src/epcsaft/_native/electrolyte.hpp src/epcsaft/_native/born.hpp \
  tests/test_equation_registry.py
git commit -m "docs: mirror electrolyte and Born equations"
```

---

### Task 4: Document properties, density closure, and CppAD authority

**Files:**

- Modify: `docs/latex/equations.tex`
- Modify with comments only: `src/epcsaft/_native/eos.cpp`
- Modify with comments only: `src/epcsaft/_native/density.hpp`
- Modify with exclusions only: `src/epcsaft/_native/derivatives.hpp`
- Regenerate: `docs/equations_registry.yaml`
- Regenerate: `docs/equations.md`
- Test: `tests/test_equation_registry.py`

**Interfaces:**

- Consumes: current `NativeModel::evaluate`, density, and CppAD owners.
- Produces: complete state/property identities and an explicit negative boundary around generated derivatives.

- [ ] **Step 1: Add failing property and derivative-boundary assertions**

```python
def test_property_identities_and_cppad_boundary_are_documented() -> None:
    tex = (REPO_ROOT / "docs" / "latex" / "equations.tex").read_text()
    for eqid in (
        "compressibility_factor",
        "pressure",
        "residual_chemical_potential",
        "log_fugacity_coefficient",
        "residual_entropy",
        "residual_enthalpy",
        "residual_heat_capacity_cp",
        "density_residual",
        "mechanical_stability",
    ):
        assert f"% EqID: {eqid}" in tex
    assert "CppAD is the sole nonlinear derivative authority" in tex
    assert "% Status: Analytical derivative" not in tex
    assert "% Status: Manual derivative" not in tex
```

Run and expect failure on the absent entries.

- [ ] **Step 2: Link state and property identities in `eos.cpp`**

Use:

```text
compressibility_factor
pressure
pressure_density_derivative
residual_chemical_potential
log_fugacity_coefficient
residual_entropy
residual_internal_energy
residual_enthalpy
residual_heat_capacity_cv
residual_heat_capacity_cp
fixed_pressure_total_derivative
```

These are exact thermodynamic identities or exact transformations of CppAD
partials. Do not expand the underlying nonlinear derivatives.

- [ ] **Step 3: Link density-domain relations and classify solver mechanics**

Use:

```text
density_residual
ideal_gas_density
density_upper_bound
mechanical_stability
```

Classify logarithmic scan construction, bracket updates, Newton candidates,
iteration limits, and root deduplication as `numerical iteration mechanics`.
Do not present the scan algorithm as thermodynamics.

- [ ] **Step 4: Classify `derivatives.hpp` without documenting generated formulas**

Add `EqExclude: CppAD-generated derivative machinery` to
`directional_second`, `thermodynamic_derivatives_impl`, and its thin active
coordinate wrappers. Add scalar/parallel utility exclusions to execution
helpers. The LaTeX derivative section must state the taped independent
coordinates and dependent defining equations in prose, without derivative
EqIDs or expanded Jacobian/Hessian formulas.

- [ ] **Step 5: Add property and density LaTeX sections**

Add:

```text
Thermodynamic Properties
Density Closure and Mechanical Domain
Derivative Authority
```

Use Gross2001 equation locators already recorded in
`docs/science/full-provider-eos.md`. Clearly distinguish residual from total
caloric properties.

- [ ] **Step 6: Regenerate and run focused property tests**

```bash
python scripts/docs/sync_equation_registry.py
python scripts/docs/sync_equation_registry.py --check
python -m pytest tests/test_equation_registry.py tests/test_properties.py tests/test_density.py -q
git diff --check
```

Expected: all PASS.

- [ ] **Step 7: Commit property/density documentation**

```bash
git add docs/latex/equations.tex docs/equations.md \
  docs/equations_registry.yaml src/epcsaft/_native/eos.cpp \
  src/epcsaft/_native/density.hpp src/epcsaft/_native/derivatives.hpp \
  tests/test_equation_registry.py
git commit -m "docs: mirror Provider property equations"
```

---

### Task 5: Document reference-state and native phase-coordinate relations

**Files:**

- Modify: `docs/latex/equations.tex`
- Modify: `docs/latex/references.bib`
- Modify with comments only: `src/epcsaft/_native/electrolyte_reference.hpp`
- Modify with comments only: `src/epcsaft/_native/phase_block.hpp`
- Modify with comments only: `src/epcsaft/_native/native_sdk.cpp`
- Modify with comments only: `src/epcsaft/reference.py`
- Regenerate: `docs/equations_registry.yaml`
- Regenerate: `docs/equations.md`
- Test: `tests/test_equation_registry.py`

**Interfaces:**

- Consumes: Task 1 support for both C++ and Python owner comments.
- Produces: complete reference/property transformation and native coordinate coverage.

- [ ] **Step 1: Add failing reference/phase coverage assertions**

Require:

```text
formula_unit_composition
reference_molality_sequence
single_ion_activity_coefficient
mean_ionic_activity_coefficient
molality_activity_conversion
practical_osmotic_coefficient
infinite_dilution_solvation_energy
transfer_free_energy
phase_total_amount
phase_mole_fraction
phase_molar_density
phase_dimensionless_helmholtz
phase_pressure
packing_volume_bounds
neutral_reference_stoichiometric_basis
```

Run `test_required_provider_equation_families_are_present` and expect failure
until the entries are added.

- [ ] **Step 2: Link reference-sequence construction**

In `electrolyte_reference.hpp`, link the formula-unit composition relation and
the fixed four-molality sequence. Classify convergence comparisons and
sequence iteration as numerical iteration mechanics.

- [ ] **Step 3: Link Python-owned exact transformations**

In `reference.py`, use Python comments for:

```text
single_ion_activity_coefficient
mean_ionic_activity_coefficient
molality_activity_conversion
practical_osmotic_coefficient
infinite_dilution_solvation_energy
transfer_free_energy
```

Place each comment immediately beside its defining arithmetic. Do not relocate
or rewrite function bodies.

- [ ] **Step 4: Link native phase-coordinate relations**

In `phase_block.hpp`, track amount summation, mole-fraction conversion,
`rho_m=n/V`, `Phi=A/(RT)=n a^r`, pressure conversion, and packing-volume-bound
relations. Classify CppAD gradient/Hessian extraction as generated derivative
machinery and classify buffer/index bookkeeping as transport mechanics.

- [ ] **Step 5: Classify native SDK transport and link its one mathematical basis**

In `native_sdk.cpp`, link the charge-balanced neutral-reference basis
construction as `neutral_reference_stoichiometric_basis`. Classify callback
buffer initialization, capacity checks, status/error propagation, table
layout, result copying, and batch indexing as transport and buffer mechanics.
Do not duplicate equations already owned by `eos.cpp` or `phase_block.hpp`.

- [ ] **Step 6: Add reference and phase LaTeX sections**

Add bibliography entry `HamerWu1972`. Use Figiel2025 for MIAC, solvation, and
transfer relations. State formula-unit molality, salt-free-solvent mass basis,
stoichiometric coefficient sum, pure-solvent reference, and units explicitly.

- [ ] **Step 7: Regenerate and run focused reference/native tests**

```bash
python scripts/docs/sync_equation_registry.py
python scripts/docs/sync_equation_registry.py --check
python -m pytest tests/test_equation_registry.py tests/test_reference_states.py tests/test_native_sdk.py -q
git diff --check
```

Expected: all PASS.

- [ ] **Step 8: Commit reference and native-coordinate documentation**

```bash
git add docs/latex/equations.tex docs/latex/references.bib \
  docs/equations.md docs/equations_registry.yaml \
  src/epcsaft/_native/electrolyte_reference.hpp \
  src/epcsaft/_native/phase_block.hpp src/epcsaft/_native/native_sdk.cpp \
  src/epcsaft/reference.py \
  tests/test_equation_registry.py
git commit -m "docs: mirror reference and phase equations"
```

---

### Task 6: Wire the mirror into canonical documentation and verify the full slice

**Files:**

- Modify: `docs/science/full-provider-eos.md`
- Modify: `README.md`
- Modify: `tests/test_equation_registry.py`
- Regenerate: `docs/equations_registry.yaml`
- Regenerate: `docs/equations.md`

**Interfaces:**

- Consumes: complete Tasks 1–5 mirror.
- Produces: the final discoverable, checked documentation feature.

- [ ] **Step 1: Add failing discoverability and complete-classification tests**

Add:

```python
def test_science_contract_links_the_checked_equation_mirror() -> None:
    science = (REPO_ROOT / "docs" / "science" / "full-provider-eos.md").read_text()
    readme = (REPO_ROOT / "README.md").read_text()
    assert "../equations.md" in science
    assert "docs/equations.md" in readme


def test_every_declared_math_owner_function_is_classified() -> None:
    module = _load_sync_module()
    entries, refs, exclusions = module.load_repository_inputs(REPO_ROOT)
    module.validate_registry(entries, refs, exclusions, module.load_bibliography_keys(REPO_ROOT))
    assert module.unclassified_owner_functions(REPO_ROOT) == []
```

Run and expect failure before the documentation links are added.

- [ ] **Step 2: Replace the stale legacy-relative citation**

In `docs/science/full-provider-eos.md`, replace the ambiguous
`docs/latex/equations.tex` archive reference with:

- the local checked Provider mirror, linked as `../equations.md`;
- the pinned legacy commit/blob as provenance; and
- a statement that C++/Python owners execute while LaTeX mirrors.

Add a compact README link under scientific documentation. Do not add a new
documentation hierarchy or duplicate the scientific contract.

- [ ] **Step 3: Perform the complete relation audit**

Review every function and arithmetic block in the declared owner files.
For each, verify one of:

- its defining mathematical relation has an `EqID`;
- it is part of a clearly bounded multi-statement relation carrying an `EqID`;
- it has an allowed `EqExclude` reason.

Use the generated exclusion report to reject vague or overbroad exclusions.
In particular, do not classify thermodynamic transformations as transport
mechanics.

- [ ] **Step 4: Regenerate and run the strict documentation gate**

```bash
python scripts/docs/sync_equation_registry.py
python scripts/docs/sync_equation_registry.py --check
python -m pytest tests/test_equation_registry.py -q
python -m ruff check scripts/docs/sync_equation_registry.py tests/test_equation_registry.py
git diff --check
```

Expected: all PASS with zero unclassified owner functions and zero implemented
equations lacking owners.

- [ ] **Step 5: Build the optional local LaTeX rendering when available**

```bash
if command -v latexmk >/dev/null 2>&1; then
  (cd docs/latex && latexmk -pdf -outdir=out equations.tex)
fi
```

Expected: success when `latexmk` exists. Remove the task-created `docs/latex/out`
directory afterward; do not stage it.

- [ ] **Step 6: Run full Provider verification**

```bash
python -m pytest -q
python -m ruff check .
python -m mypy --no-incremental
git diff --check
bash "$HOME/.codex/hooks/codex-cleanup.sh" --repo-root .
```

Expected: all repository checks PASS. The cleanup audit may list unrelated
pre-existing caches; remove only ignored artifacts created by this task.

- [ ] **Step 7: Review the diff for runtime neutrality**

Run:

```bash
git diff --stat e80c14f..HEAD
git diff e80c14f..HEAD -- src/epcsaft/_native src/epcsaft/reference.py
```

Verify source changes contain comments only. Confirm no public header,
function body, CMake target, Python dependency, package-data rule, or generated
PDF changed.

- [ ] **Step 8: Commit the completed documentation mirror**

```bash
git add README.md docs/science/full-provider-eos.md docs/latex \
  docs/equations.md docs/equations_registry.yaml \
  scripts/docs/sync_equation_registry.py tests/test_equation_registry.py \
  src/epcsaft/_native src/epcsaft/reference.py
git commit -m "docs: complete Provider equation mirror"
```

- [ ] **Step 9: Record final identities**

```bash
git status --short --branch
git rev-parse HEAD
git rev-parse HEAD^{tree}
sha256sum docs/latex/equations.tex docs/equations_registry.yaml docs/equations.md
```

Expected: clean feature branch and reproducible document hashes.

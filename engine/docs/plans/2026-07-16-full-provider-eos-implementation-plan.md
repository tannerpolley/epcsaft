# Full Provider EOS Implementation Plan

> **Execution:** Use `superpowers:subagent-driven-development` in this session
> or `superpowers:executing-plans` in a fresh session. Apply
> `superpowers:test-driven-development` to each runtime leaf and
> `superpowers:verification-before-completion` before completion claims.

**Goal:** Implement the approved full provider EOS and property boundary in
the clean `epcsaft` package, including association, Debye--Huckel, Born SSM+DS,
residual properties, and branch-safe density closure.

**Design:**
`docs/designs/2026-07-16-full-provider-eos-design.md`

**Repository:** `~/Workspaces/Engineering/ePC-SAFT-project/ePC-SAFT`

**Evidence source:** `~/Workspaces/Engineering/ePC-SAFT`

**Working rule:** Work on provider `main` through small local commits. Do not
push, release, or change runtime authority without explicit user approval.
Keep one Python package, one `_core` module, one native target, one parameter
lifecycle, and one Helmholtz owner.

## Execution status

Tasks 1--9 are complete in local commits through branch-safe specified-pressure
closure. Task 10 is deferred because reference-state orchestration is not part
of the current EOS state boundary and requires its own admitted dilution paths
and literature campaign. Tasks 11--12 remain the independent-validation and
promotion gates; neither is implied by a green package test suite.

## Program verification commands

Use no-bytecode/no-cache tests during development:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src \
  uv run --with pytest --with Pint python -m pytest -p no:cacheprovider -q
```

Native changes require a wheel rebuild and isolated test:

```bash
ROOT=$PWD
TMP=$(mktemp -d)
uv build --wheel --out-dir "$TMP/dist"
uv venv "$TMP/venv" --python 3.13
uv pip install --python "$TMP/venv/bin/python" "$TMP"/dist/*.whl pytest
(cd "$TMP" && PYTHONDONTWRITEBYTECODE=1 \
  "$TMP/venv/bin/python" -I -m pytest "$ROOT/tests" -p no:cacheprovider -q)
rm -rf "$TMP"
```

Run after each commit:

```bash
git diff --check HEAD^
bash "$HOME/.codex/hooks/codex-cleanup.sh" --repo-root .
git status --short --branch
```

## Task 1: Freeze the slice-local scientific specification

**Files**

- Create: `docs/science/full-provider-eos.md`
- Evidence: lab `docs/latex/equations.tex`, retained primary papers, and the
  archived native contribution files

**Steps**

1. Record every equation consumed by hard-chain, dispersion, association,
   Debye--Huckel, Born SSM+DS, pressure, fugacity, and residual caloric
   properties. Each entry must name its primary source, locator, status, units,
   variables, implementation owner, and bounded domain.
2. Inspect the primary Figiel--Yu--Held article or source PDF, not the Markdown
   extraction alone, and resolve equation 7 against its ionic-region
   permittivity prose. Record the exact SSM+DS expression.
3. Exclude handwritten Debye--Huckel composition derivatives from production;
   CppAD will derive them from the defining screening and Helmholtz equations.
4. Record the selected association mixing and mass-action equations, including
   site-pair indexing and cross-association corrections.
5. Record the residual property identities and the exact derivatives they
   consume. Mark total ideal-gas properties out of scope.
6. Review every entry against the primary source and approved design, then
   commit:

```bash
git add docs/science/full-provider-eos.md
git commit -m "docs: specify full provider equations"
```

**Stop condition:** The primary paper cannot resolve the SSM+DS expression or
an equation required by the chosen model.

## Task 2: Add bundle-level scientific model records

**Files**

- Modify: `src/epcsaft/records.py`
- Modify: `src/epcsaft/parameters.py`
- Modify: `src/epcsaft/_bundle_io.py`
- Modify: `src/epcsaft/_canonical.py`
- Modify: `tests/test_records.py`
- Modify: `tests/test_parameters.py`

**Contract**

Add one typed `ModelParameterRecord` in `epcsaft.records` with closed families:

- `ionic_region_relative_permittivity`, dimensionless;
- `dielectric_ion_suppression_coefficient`, dimensionless.

Bundles store these in optional `model.csv` with columns
`record_id,family,value,unit,source_id,locator,domain_id`.
`ParameterBundle.from_records` accepts `models=()`. A selected `ParameterSet`
retains all bundle-level model records because they do not belong to one
component. Fingerprints include them. No dictionary API or default appears.

**Steps**

1. Add red tests for type, closed families, finite values, source/domain
   validation, strict `model.csv` inventory, row-order independence, and
   fingerprint changes.
2. Implement the record and parser with the existing normalization helpers.
3. Update canonical projections and bundle selection without adding a second
   parameter lifecycle.
4. Run `tests/test_records.py tests/test_parameters.py` and commit:

```bash
git add src/epcsaft/records.py src/epcsaft/parameters.py \
  src/epcsaft/_bundle_io.py src/epcsaft/_canonical.py \
  tests/test_records.py tests/test_parameters.py
git commit -m "feat: add source-backed model parameters"
```

## Task 3: Resolve one immutable full-model input

**Files**

- Create: `src/epcsaft/_resolve.py`
- Modify: `src/epcsaft/eos.py`
- Modify: `src/epcsaft/__init__.py`
- Create: `tests/test_model_resolution.py`
- Modify: `tests/test_eos.py`

**Contract**

Move model-specific record resolution out of `eos.py` into one private
resolver. It produces one immutable internal model containing ordered
component values, pair matrices, temperature correlations, site topology,
site-pair parameters, and the two model parameters.

The resolver recognizes three complete topologies: neutral, associating, and
ionic SSM+DS. It rejects partial parameter sets. It never supplies a missing
zero, one, diagonal, combining rule, charge, dielectric value, or association
pair. Source equations may generate diagonal/cross quantities only after all
their operands exist.

Add public `EosError(ValueError)` with stable `code` and immutable `context`.
Keep `ParameterError` for construction/resolution failures.

**Steps**

1. Add red tests for all three complete topologies and every incomplete family.
2. Add tests proving current methane/ethane fingerprints and results do not
   change.
3. Implement `_resolve.py`, keep `EPCSAFT(parameters)` as the only constructor,
   and update root exports only for `EosError`.
4. Run the complete Python suite and commit:

```bash
git add src/epcsaft/_resolve.py src/epcsaft/eos.py src/epcsaft/__init__.py \
  tests/test_model_resolution.py tests/test_eos.py
git commit -m "feat: resolve full provider model input"
```

## Task 4: Add one complete source-backed electrolyte bundle

**Files**

- Create: `src/epcsaft/data/bundles/figiel-2025-reference-electrolytes/1/*`
- Modify: `tests/test_parameters.py`
- Modify: `tests/test_model_resolution.py`

**Steps**

1. Build the compact water/sodium/chloride bundle from Figiel--Yu--Held 2025
   Tables 1--5 and cited primary sources. Include complete PC-SAFT, molar-mass,
   charge, dielectric, Born, solvation, association, pair, model, and domain
   records. Preserve every source locator and do not fill blank cells.
2. Keep the existing disconnected parser fixture. Do not reinterpret it as a
   physical system.
3. Add exact parameter, fingerprint, complete-topology, and component-order
   tests. Do not add model-output goldens before an equation implementation
   exists.
4. Run parameter and model-resolution tests and commit:

```bash
git add src/epcsaft/data tests/test_parameters.py tests/test_model_resolution.py
git commit -m "test: add source-backed electrolyte bundle"
```

**Stop condition:** A required value lacks a defensible primary source,
locator, unit, domain, or use basis.

## Task 5: Refactor the shared kernel and implement association

**Files**

- Create: `src/epcsaft/_native/association.hpp`
- Create: `src/epcsaft/_native/model.hpp`
- Create: `src/epcsaft/_native/helmholtz.hpp`
- Modify: `src/epcsaft/_native/eos.hpp`
- Modify: `src/epcsaft/_native/eos.cpp`
- Modify: `src/epcsaft/_native/module.cpp`
- Modify: `tests/test_eos.py`
- Create: `tests/test_association.py`

**Steps**

1. Confirm the current neutral golden and permutation tests pass, then add red
   tests for one-site and 2B water topologies, mass-action residuals,
   nonzero association energy, low-density reduction, and invalid topology.
2. Add independent finite-difference checks of density, temperature, and
   composition directions. Keep these as test oracles, never runtime modes.
3. Move immutable native parameters to `model.hpp` and the promoted scalar
   neutral equations to `helmholtz.hpp` only as needed to add association.
   Preserve one `residual_helmholtz` entry, one CppAD owner, one native model,
   and the exact Slice 1 neutral results.
4. Implement the site-fraction value solve outside CppAD recording. Return
   convergence residual, iteration count, and Jacobian conditioning to the
   native orchestrator.
5. Implement CppAD partials of the defining mass-action and Helmholtz equations
   and exact implicit linear solves for consumed sensitivities.
6. Reject nonphysical fractions, nonconvergence, and ill-conditioned
   Jacobians through `EosError` translation.
7. Build, test in isolation, and commit:

```bash
git add src/epcsaft/_native tests/test_association.py
git commit -m "feat: evaluate associating PC-SAFT"
```

## Task 6: Implement the dielectric model and Debye--Huckel term

**Files**

- Create: `src/epcsaft/_native/electrolyte.hpp`
- Modify: `src/epcsaft/_native/model.hpp`
- Modify: `src/epcsaft/_native/helmholtz.hpp`
- Modify: `src/epcsaft/_native/eos.cpp`
- Create: `tests/test_debye_huckel.py`

**Steps**

1. Add red tests for salt-free solvent mass-fraction mixing, ion-fraction
   suppression, SI number-density conversion, screening length, zero-charge
   reduction, charge-sign symmetry, and source-domain rejection.
2. Add composition-direction finite differences against CppAD output. Do not
   copy the archived handwritten `dkappa/dx` formula.
3. Implement the bulk dielectric and Debye--Huckel equations from the frozen
   scientific specification.
4. Reject nonpositive dielectric values, non-electroneutral ordinary states,
   absent molar masses, and incomplete charged topology.
5. Build, test in isolation, and commit:

```bash
git add src/epcsaft/_native tests/test_debye_huckel.py
git commit -m "feat: add Debye-Huckel electrostatics"
```

## Task 7: Implement Born SSM+DS

**Files**

- Create: `src/epcsaft/_native/born.hpp`
- Modify: `src/epcsaft/_native/model.hpp`
- Modify: `src/epcsaft/_native/helmholtz.hpp`
- Modify: `src/epcsaft/_native/eos.cpp`
- Create: `tests/test_born.py`

**Steps**

1. Add red equation-level tests for fitted Born diameter, charge-normalized
   shell thickness, component solvation-factor mixing, separate bulk and ionic
   dielectric regions, and the combined SSM+DS energy.
2. Prove reductions for zero charge and unit solvent solvation factors. Check
   charge-squared and diameter scaling.
3. Implement the frozen Figiel equations. Read both model parameters from the
   resolved input; a lexical gate must reject `eps_ion = 8.0` and an embedded
   `7.01` in native production source.
4. Add CppAD density, temperature, and composition directional checks.
5. Build, test in isolation, and commit:

```bash
git add src/epcsaft/_native tests/test_born.py
git commit -m "feat: add Born SSM DS contribution"
```

## Task 8: Add residual chemical, fugacity, and caloric properties

**Files**

- Create: `src/epcsaft/_native/derivatives.hpp`
- Modify: `src/epcsaft/_native/eos.hpp`
- Modify: `src/epcsaft/_native/eos.cpp`
- Modify: `src/epcsaft/_native/module.cpp`
- Modify: `src/epcsaft/eos.py`
- Create: `tests/test_properties.py`

**Contract**

Expand `EosResult` rather than adding a parallel result. Add component-ordered
dimensionless residual chemical potentials and log fugacity coefficients plus
Pint residual entropy, internal energy, enthalpy, `C_V`, and `C_P` quantities.
Do not claim total ideal-gas properties.

**Steps**

1. Add red thermodynamic identity tests and independent directional checks for
   neutral, associating, and ionic states.
2. Record one CppAD function over temperature, molar density, and independent
   composition coordinates. Reconstruct the constrained composition
   derivatives used by chemical potentials from that owner.
3. Include association implicit sensitivities at each consumed order. Reject
   singular association and density states before returning heat capacities.
4. Keep contribution fields as views of the same derivative evaluation.
5. Build, run all tests in isolation, and commit:

```bash
git add src/epcsaft/_native src/epcsaft/eos.py tests/test_properties.py
git commit -m "feat: add residual EOS properties"
```

## Task 9: Add branch-safe density closure and T/P evaluation

**Files**

- Create: `src/epcsaft/_native/density.hpp`
- Modify: `src/epcsaft/_native/eos.hpp`
- Modify: `src/epcsaft/_native/eos.cpp`
- Modify: `src/epcsaft/_native/module.cpp`
- Modify: `src/epcsaft/eos.py`
- Create: `tests/test_density.py`

**Contract**

Add immutable public `DensityDiagnostics` and extend
`EPCSAFT.evaluate` to accept exactly one of `molar_density` or `pressure`.
`phase=None` accepts one stable root; `"vapor"` and `"liquid"` select the
lowest- and highest-density stable roots when a multi-root state requires
disambiguation. A lone root remains labelled `"single"` even when a phase hint
was supplied.

**Steps**

1. Add red tests for one-root, three-root, root-slope, and near-packing
   behavior. Prove that phase hints cannot relabel a lone root.
2. Implement a deterministic log-density scan bounded by positive density and
   the mixture packing limit, then refine candidate brackets.
3. Report bracket, residual, `dR/d(rho)`, stability, iterations, and selected
   branch for successful states; failed solves report through `EosError`. Do
   not expose warm-start caches until a demonstrated workload needs them.
4. Run the same density path for neutral, associating, and ionic kernels.
5. Build, run all tests in isolation, and commit:

```bash
git add src/epcsaft/_native src/epcsaft/eos.py tests/test_density.py
git commit -m "feat: solve branch-safe density states"
```

## Task 10: Add explicit electrolyte reference-state quantities

**Files**

- Create: `src/epcsaft/reference.py`
- Modify: `src/epcsaft/eos.py`
- Modify: `src/epcsaft/__init__.py`
- Create: `tests/test_reference_states.py`

**Contract**

Expose only source-defined reference calculations needed by the bounded
electrolyte campaign: x-process infinite-dilution solvation, transfer between
two explicit solvent states, single-ion activity relative to its declared
infinite-dilution convention, and stoichiometric mean ionic activity. Each
method requires an electroneutral dilution path and names its output scale.

**Steps**

1. Add red tests for reference-state identity, standard-state conversion,
   electroneutral limiting paths, stoichiometric means, and rejection of a
   lone non-electroneutral ion path.
2. Implement these calculations as orchestration over `EPCSAFT.evaluate` and
   its fugacity output. Do not copy Helmholtz equations into Python.
3. Keep concentration-scale conversion formulas separate from EOS state
   evaluation and preserve units in results.
4. Run all Python and isolated-wheel tests and commit:

```bash
git add src/epcsaft/reference.py src/epcsaft/eos.py \
  src/epcsaft/__init__.py tests/test_reference_states.py
git commit -m "feat: add electrolyte reference states"
```

## Task 11: Run installed-artifact literature validation

**Repositories**

- Provider: current repository
- Validation: `../ePC-SAFT-validation`

**Provider files**

- Modify: `README.md`
- Modify: `CONTEXT.md`
- Modify: `AGENTS.md`
- Modify: `ARCHITECTURE.yaml`
- Add only compact package anchors required by Tasks 4--10

**Validation files**

- Add a bounded associating campaign
- Add a bounded water/sodium-chloride electrolyte campaign
- Retain source data, transformation record, tolerances, CSV output, and plot

**Steps**

1. Build an immutable provider wheel from the candidate commit and record its
   hash.
2. Install it into an isolated validation environment with no provider source,
   lab, or migration checkout on `sys.path`.
3. Compare retained literature data with model predictions for at least one
   associating property and one electrolyte property. Produce the required
   retained plots and numerical tables.
4. Keep solver, numerical, and physical acceptance separate. Investigate any
   disagreement before changing tolerance or scope.
5. Update provider documentation and the architecture cost ledger. Confirm no
   alternate derivative modes, hidden model constants, compatibility seams,
   duplicate kernels, sibling dependencies, or unreachable files.
6. Commit provider and validation changes in their owning repositories. Do not
   push or promote yet.

## Task 12: Final verification and promotion preparation

**Files**

- Provider: update `ARCHITECTURE.yaml` with final measured surface
- Migration: prepare candidate validation and promotion receipts after the
  immutable validation passes

**Steps**

1. Run the complete provider source suite, wheel/sdist build, isolated wheel
   suite, Git diff check, negative-space scan, and cleanup audit.
2. Run the complete validation campaigns against the immutable wheel.
3. Review public symbols, schemas, dependencies, native targets, source files,
   tests, and documentation against the approved design. Remove obsolete
   owners in the same checkpoint.
4. Prepare receipts with exact commits, artifact hashes, commands,
   environments, datasets, tolerances, and reviewer identity.
5. Stop for explicit user approval before push, publication, authority
   promotion, or release.

## Completion evidence

Implementation is ready for promotion only when:

- all five residual Helmholtz contributions share one source-backed kernel;
- explicit-density and T/P routes return the same typed result;
- association derivatives use implicit differentiation of the defining
  residual;
- density requests preserve branch identity;
- the full admitted residual property set passes independent derivative and
  thermodynamic-identity checks;
- source-backed associating and electrolyte installed-artifact campaigns pass;
- the package contains no hidden physical constants, alternate derivative
  backends, compatibility shims, lab dependencies, or duplicate owners; and
- the final architecture cost ledger justifies every increase from Slice 1.

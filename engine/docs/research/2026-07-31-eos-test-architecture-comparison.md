# Thermodynamic EOS test-architecture comparison

Date: 2026-07-31

## Scope and method

This note compares the current test source and CI configuration of five official
repositories at pinned commits. All project observations below are **verified**
against those primary sources; recommendations are **inferences** for
`ePC-SAFT/ePC-SAFT-eos`.

| Project | Snapshot | Why it is useful here |
|---|---|---|
| teqp | [`5f62a6f`](https://github.com/usnistgov/teqp/tree/5f62a6f515d517e39c3fb035c11a03524ffa3ad6) | C++ Helmholtz-EOS core, automatic differentiation, JSON model construction, C and Python build surfaces |
| feos | [`c658aea`](https://github.com/feos-org/feos/tree/c658aeab484f7a7096bfbf5425e40effd60da167) | Rust EOS/DFT core, PC-SAFT and ePC-SAFT, exact derivative propagation, Python bindings |
| Clapeyron.jl | [`fcadeff`](https://github.com/ClapeyronThermo/Clapeyron.jl/tree/fcadeff7289bdfe1492bd0b3f4f69cc98103ec3b) | Broad Julia thermodynamics API, many EOS families, reusable invariant helpers, database and solver tests |
| zmeri/PC-SAFT | [`db1a229`](https://github.com/zmeri/PC-SAFT/tree/db1a2293499236aeafe736253f0be008a455888e) | The user-requested simplicity reference: Cython/C++ PC-SAFT with electrolyte and polar terms and raw NumPy parameter dictionaries |
| CoolProp | [`e8fff56`](https://github.com/CoolProp/CoolProp/tree/e8fff5694e8fba11d7cbd9a80789a3b194ebbe59) | Additional comparator for a mature native core, installed Python bindings, wheel matrices, and opt-in downstream compile contracts |

Raw test counts are not directly comparable across Catch2 sections, Rust test
functions, Julia `@testset`s, and parametrized pytest cases. The useful unit of
comparison is the failure claim each architecture assigns to unit, integration,
installed-artifact, or broad validation evidence.

The source-level conclusions were also checked against the local forks on the
same date. `~/Workspaces/Forks/teqp` is a clean clone of the
pinned teqp commit above; the current FeOs and Clapeyron working trees under
`~/Workspaces/Forks` retain the same test-architecture patterns
described here. The existing FeOs and Clapeyron working trees were inspected
read-only because they contain user work.

## Project findings

### teqp

**Organization and naming.** teqp places its C++ Catch2 suite under
`src/tests/`, mostly as `catch_test_<model-or-capability>.cxx`: PC-SAFT,
generic SAFT, polar SAFT, composition derivatives, cubics, multifluid models,
phase equilibrium, and other model families each have a file. The PC-SAFT file
uses short behavior-oriented case names such as “Single alphar check value,”
“Check 0n derivatives,” and “Check PCSAFT with kij”; the naming is not perfectly
uniform, but the file boundary supplies the model context
([PC-SAFT tests](https://github.com/usnistgov/teqp/blob/5f62a6f515d517e39c3fb035c11a03524ffa3ad6/src/tests/catch_test_PCSAFT.cxx#L20-L149)).

**Evidence emphasis and parameters.** The suite combines a small number of
fixed residual-Helmholtz values with generic derivative agreement across
autodiff, multicomplex, complex-step, and multiprecision routes. It also checks
virial limits and a model matrix that applies the same virial fixture to several
EOS families
([derivative checks](https://github.com/usnistgov/teqp/blob/5f62a6f515d517e39c3fb035c11a03524ffa3ad6/src/tests/catch_test_PCSAFT.cxx#L64-L158),
[generic model matrix](https://github.com/usnistgov/teqp/blob/5f62a6f515d517e39c3fb035c11a03524ffa3ad6/src/tests/catch_test_modelset.cxx#L110-L145)).
Parameters are usually inline C++ arrays or inline JSON model specifications;
multifluid tests additionally use repository JSON databases. Failure tests are
present for malformed matrix shape and unsupported operations, but they are
secondary to values, derivatives, and thermodynamic limits
([inline `k_ij` contracts](https://github.com/usnistgov/teqp/blob/5f62a6f515d517e39c3fb035c11a03524ffa3ad6/src/tests/catch_test_PCSAFT.cxx#L199-L246)).

**Execution.** One CMake target glob-compiles the Catch sources and registers a
single CTest entry
([CMake](https://github.com/usnistgov/teqp/blob/5f62a6f515d517e39c3fb035c11a03524ffa3ad6/CMakeLists.txt#L295-L313)).
The visible CI builds Release mode with Clang on macOS, runs CTest, builds the C
shared library, and verifies that a Python wheel builds on Python 3.9; it does
not run a broad installed-wheel test matrix in that workflow
([workflow](https://github.com/usnistgov/teqp/blob/5f62a6f515d517e39c3fb035c11a03524ffa3ad6/.github/workflows/runcatch.yml#L5-L34)).

### feos

**Organization and naming.** feos uses two layers. Small equation-contribution
and parameter tests live beside Rust implementation modules under `#[cfg(test)]`;
public workflows live under feature-gated integration modules such as
`tests/pcsaft/{properties,state_creation,critical_point,tp_flash,...}.rs`
([PC-SAFT module index](https://github.com/feos-org/feos/blob/c658aeab484f7a7096bfbf5425e40effd60da167/crates/feos/tests/pcsaft/mod.rs#L1-L9)).
Names are concise claims—`ideal_gas_pressure`, `hard_sphere_mix`,
`test_dln_phi_dp`, `test_critical_point_limits`—because module paths carry the
rest of the context.

**Evidence emphasis and parameters.** feos strongly exercises generic
thermodynamic identities: ideal plus residual pressure equals total pressure,
a pure model equals the zero-fraction limit of a mixture, solved pressure
reconstructs the requested pressure, and critical-mixture endpoints equal the
pure critical states
([core PC-SAFT invariants](https://github.com/feos-org/feos/blob/c658aeab484f7a7096bfbf5425e40effd60da167/crates/feos/src/pcsaft/eos/mod.rs#L396-L525),
[critical limits](https://github.com/feos-org/feos/blob/c658aeab484f7a7096bfbf5425e40effd60da167/crates/feos/tests/pcsaft/critical_point.rs#L51-L76)).
It pairs those invariants with contribution anchors for hard-sphere, ionic,
Born, dipole, and quadrupole terms. Derivative-sensitive flash tests construct
PC-SAFT directly from raw nested arrays, seed an active parameter, and compare
the propagated derivative with a finite difference
([raw flash fixture](https://github.com/feos-org/feos/blob/c658aeab484f7a7096bfbf5425e40effd60da167/crates/feos/tests/pcsaft/px_flashes.rs#L12-L71)).
Other integration tests intentionally load JSON parameter files, while generic
parameter-parser tests use inline JSON records. Thus scientific numerics and
database integration are separate claims rather than every numerical test
depending on the database
([property derivative using JSON](https://github.com/feos-org/feos/blob/c658aeab484f7a7096bfbf5425e40effd60da167/crates/feos/tests/pcsaft/properties.rs#L9-L43),
[parameter records](https://github.com/feos-org/feos/blob/c658aeab484f7a7096bfbf5425e40effd60da167/crates/feos-core/tests/parameters.rs#L17-L63)).
Malformed-input testing is concentrated in the parameter layer; it does not
dominate every numerical capability.

**Execution.** CI runs release-mode `cargo test` separately for core crates,
Python bindings, each EOS feature, and PC-SAFT/GC-PC-SAFT DFT features on Linux.
The Python-binding job uses Python 3.14. Criterion performance targets are
declared as benches rather than ordinary tests
([test workflow](https://github.com/feos-org/feos/blob/c658aeab484f7a7096bfbf5425e40effd60da167/.github/workflows/test.yml#L12-L64),
[bench targets](https://github.com/feos-org/feos/blob/c658aeab484f7a7096bfbf5425e40effd60da167/crates/feos/Cargo.toml#L64-L92)).

### Clapeyron.jl

**Organization and naming.** Clapeyron separates tests by responsibility:
database parsing, differentials, solvers, model families, EOS methods, flash
methods, electrolytes, estimation, and reported issues. `runtests.jl` is an
explicit suite index rather than one giant capability file
([suite index](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/fcadeff7289bdfe1492bd0b3f4f69cc98103ec3b/test/runtests.jl#L216-L234)).
Nested `@testset` names identify the model or method; individual assertions are
usually unnamed.

**Evidence emphasis and parameters.** Reusable helpers test extensivity/Gibbs–
Duhem behavior, pure-component limits, volume-pressure closure, scale
invariance, recombination, and binary-parameter round trips
([generic helpers](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/fcadeff7289bdfe1492bd0b3f4f69cc98103ec3b/test/runtests.jl#L39-L153)).
Model files then combine those invariants with a few contribution anchors; the
PC-SAFT family is a representative example
([PC-SAFT model tests](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/fcadeff7289bdfe1492bd0b3f4f69cc98103ec3b/test/test_models_saft_pc.jl#L1-L89)).
A manufactured analytic model checks first-, second-, and third-order
thermodynamic differential machinery generically
([differential tests](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/fcadeff7289bdfe1492bd0b3f4f69cc98103ec3b/test/test_differentials.jl#L5-L85)).
Most model tests use named database constructors, but the database suite also
proves direct in-memory `NamedTuple` arrays and association `Dict`s, which is
close to the desired raw-fixture style
([in-memory PC-SAFT parameters](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/fcadeff7289bdfe1492bd0b3f4f69cc98103ec3b/test/test_database.jl#L446-L475)).

**Execution and broad validation.** CI spans Julia LTS, current, and nightly on
Linux and Windows plus Apple Silicon macOS. The suite distributes test files
across that matrix while local runs execute all files
([CI matrix](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/fcadeff7289bdfe1492bd0b3f4f69cc98103ec3b/.github/workflows/ci.yml#L23-L67),
[distribution logic](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/fcadeff7289bdfe1492bd0b3f4f69cc98103ec3b/test/utils.jl#L1-L61)).
The expensive multi-seed TP-flash tuning and validation workflow has a separate
manual runner and result aggregation; it is not included by `runtests.jl`
([benchmark workflow](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/fcadeff7289bdfe1492bd0b3f4f69cc98103ec3b/test/tpflash/README.md#L1-L38)).

### zmeri/PC-SAFT

**Organization and naming.** This repository has one 1,491-line pytest module
with twelve property-oriented tests: residual enthalpy, entropy and Gibbs
energy, density, limiting-index behavior, two flash routes, osmotic
coefficient, heat of vaporization, temperature derivative, heat capacity, and
pressure
([test module](https://github.com/zmeri/PC-SAFT/blob/db1a2293499236aeafe736253f0be008a455888e/tests/test_cython.py#L13-L166)).

**Evidence emphasis and parameters.** Nearly every case declares `m`, `sigma`,
`epsilon/k`, association, dipole, charge, dielectric, and binary matrices as
local NumPy arrays, puts them into a plain dictionary, and compares one or more
results with Aspen Plus, CoolProp, IAPWS, DIPPR, or literature values. The
electrolyte density case shows this style in one place
([raw electrolyte case](https://github.com/zmeri/PC-SAFT/blob/db1a2293499236aeafe736253f0be008a455888e/tests/test_cython.py#L292-L319)).
The implementation also combines the dipolar contribution into one
`ares_polar` result rather than exposing a test abstraction per internal polar
subterm
([polar implementation](https://github.com/zmeri/PC-SAFT/blob/db1a2293499236aeafe736253f0be008a455888e/pcsaft_electrolyte.cpp#L680-L788)).
The simplicity is real, but the suite has little generic derivative,
permutation, ABI/FFI, malformed-buffer, or installed-consumer evidence. Several
large property functions repeat the same raw parameters rather than sharing a
small case matrix.

**Execution.** The official workflow builds wheels on Linux, Windows, Intel
macOS, and Apple Silicon macOS and builds an sdist, but it does not configure a
wheel test command or invoke pytest
([wheel workflow](https://github.com/zmeri/PC-SAFT/blob/db1a2293499236aeafe736253f0be008a455888e/.github/workflows/build_wheels.yml#L10-L63)).
It is therefore a useful fixture-style reference, not a sufficient evidence
model for this repository.

### CoolProp

**Organization and naming.** CoolProp separates a broad C++ Catch2 suite under
`src/Tests/`, an installed Python-interface pytest suite, wheel smoke tests, and
developer/consumer contracts. Its modern binding tests are explicitly scoped
to high-level property calls, `AbstractState`, capsule-shim unit conventions,
and approved interface removals
([Python suite scope](https://github.com/CoolProp/CoolProp/blob/e8fff5694e8fba11d7cbd9a80789a3b194ebbe59/wrappers/Python/pytest/README.md#L1-L33),
[representative API tests](https://github.com/CoolProp/CoolProp/blob/e8fff5694e8fba11d7cbd9a80789a3b194ebbe59/wrappers/Python/pytest/test_nanobind_interface.py#L31-L108)).

**Execution and separation.** Native Catch tests build and run as their own CI
job. Wheel CI runs a four-platform representative matrix for ordinary pull
requests and expands to the full matrix on protected pushes or an explicit
label; each wheel gets a compact import/property/capsule smoke test
([Catch workflow](https://github.com/CoolProp/CoolProp/blob/e8fff5694e8fba11d7cbd9a80789a3b194ebbe59/.github/workflows/test_catch2.yml#L24-L57),
[matrix policy](https://github.com/CoolProp/CoolProp/blob/e8fff5694e8fba11d7cbd9a80789a3b194ebbe59/.github/workflows/python_buildwheels.yml#L1-L12),
[wheel smoke](https://github.com/CoolProp/CoolProp/blob/e8fff5694e8fba11d7cbd9a80789a3b194ebbe59/.github/workflows/python_cibuildwheel.yml#L93-L142)).
A downstream Cython compile contract is deliberately opt-in because it builds
against an installed artifact and should not enter default pytest collection
([opt-in contract](https://github.com/CoolProp/CoolProp/blob/e8fff5694e8fba11d7cbd9a80789a3b194ebbe59/dev/pdsim_cimport_contract/test_pdsim_contract.py#L1-L38)).
That separation is directly relevant to the ePC-SAFT installed C++ consumer.

## Comparison with ePC-SAFT EOS

The current EOS repository owns more than the simple PC-SAFT comparator. It has
a Python API and immutable parameter lifecycle, a C++/CppAD numerical owner,
density and electrolyte-reference roots, and an append-only C capsule plus
installed declaration header and C++ consumer. The ABI explicitly requires
prefix-size negotiation, caller-owned result buffers, exact coordinate order,
and a new major version for semantic layout changes
([native transport ADR](https://github.com/ePC-SAFT/ePC-SAFT-eos/blob/86161c963cc7714a45af0eccc550c36cc3038808/docs/adr/0001-native-consumer-transport.md#L242-L290)).
Its admission gate also requires installed-wheel evidence, public-EOS parity,
exact derivative payloads, and fingerprint/domain rejection
([admission gate](https://github.com/ePC-SAFT/ePC-SAFT-eos/blob/86161c963cc7714a45af0eccc550c36cc3038808/docs/adr/0001-native-consumer-transport.md#L292-L319)).
Those are legitimate unique responsibilities, not test-suite ceremony.

The imbalance is architectural rather than simply numerical. At the pinned
snapshot, `tests/test_native_sdk.py` is 8,370 lines—more than half of the
14,413 lines across `tests/test_*.py`—and mixes ctypes declarations, model
fixtures, installed-header compilation, descriptor metadata, numerical parity,
derivatives, batching, roots, and many callback-specific failure matrices in
one module
([current SDK tests](https://github.com/ePC-SAFT/ePC-SAFT-eos/blob/86161c963cc7714a45af0eccc550c36cc3038808/tests/test_native_sdk.py)).
The repository's tracked local actions run the entire pytest suite, Ruff plus
mypy, and a wheel build, but the pinned tree has no repository-owned GitHub
Actions workflow that visibly assigns fast, installed-consumer, and broad
validation tiers
([local actions](https://github.com/ePC-SAFT/ePC-SAFT-eos/blob/86161c963cc7714a45af0eccc550c36cc3038808/.codex/environments/environment.toml#L8-L23)).

The external pattern is consistent:

- teqp and feos keep model-specific numerical anchors, but reusable derivative
  or invariant machinery carries most of the general claim;
- Clapeyron makes generic thermodynamic properties explicit helpers and applies
  them to many models;
- feos and Clapeyron test both raw/in-memory parameters and database loading,
  but do not make every equation test a database test;
- zmeri/PC-SAFT demonstrates the desired low-friction raw parameter style, but
  not the evidence needed for an exact native derivative ABI; and
- CoolProp separates native science, installed binding smoke tests, and costly
  downstream compile contracts.

## Recommendations

All recommendations in this section are **inferences** from the verified source
comparison and this repository's stated ownership.

### 1. Organize by evidence owner

Use directories to remove repetitive prefixes and keep each module below a
reviewable size:

```text
tests/
  fixtures.py                    # one raw-parameter adapter and case matrices
  unit/
    test_records.py              # typed records and canonicalization
    test_parameter_io.py         # bundle/schema/data integration only
    test_helmholtz.py             # contribution values and limiting behavior
  numerical/
    test_invariants.py            # extensivity, permutation, pure limits, sums
    test_derivatives.py           # generic CppAD vs independent oracle matrices
    test_density_roots.py         # branch identity and pressure closure
    test_reference_states.py      # electrolyte-reference roots and identities
  api/
    test_frontend.py
    test_properties.py
  sdk/
    test_abi.py                   # capsule, prefix/layout, no-write semantics
    test_callbacks.py             # table-driven public parity and derivatives
    test_batches.py               # scalar parity, cancellation, serialization
    test_installed_consumer.py    # one installed-header compile/link/run gate
```

Keep full literature grids, all-species catalog sweeps, plots, and cross-package
acceptance in `ePC-SAFT-validation`. Keep reusable canonical parameter packets
and their exhaustive schema/provenance checks in `ePC-SAFT-data`. EOS should
retain only the compact runtime-copy result checks, generic numerical tests,
and a few scientific anchors required to falsify its equations.

### 2. Adopt claim-oriented names

Within a category, use `test_<subject>__<claim>` and let the directory provide
the repeated context. Examples:

- `sdk/test_abi.py::test_table_prefix__remains_append_only`
- `sdk/test_callbacks.py::test_phase_block__matches_python_owner`
- `numerical/test_derivatives.py::test_helmholtz__gradient_and_hessian_match_oracle`
- `numerical/test_invariants.py::test_residual_helmholtz__is_permutation_invariant`
- `numerical/test_density_roots.py::test_pressure_closure__preserves_requested_branch`

Use parametrization IDs for the varying regime—`neutral-pure`,
`neutral-binary`, `associating-2b`, `aqueous-1-1`, `multivalent`—instead of
embedding paper, species, callback, and failure details into every function
name. Use “exact” only for an analytic/AD identity, not ordinary floating-point
agreement.

### 3. Use the public dictionary schema for raw numerical fixtures

Numerical equation tests should not load catalog bundles. Store explicit test
numbers once in immutable tuples inside a module-level mapping and construct
the public parameter object through `Parameters.from_dictionary(...)`:

```python
RAW_MODELS = {
    "neutral_binary": {
        "schema": "epcsaft.parameters",
        "schema_version": 1,
        "components": ("a", "b"),
        "parameters": {
            "m": (1.5, 2.2),
            "s": (3.4, 3.8),
            "e": (180.0, 250.0),
            "z": (0, 0),
            "k_ij": ((0.0, 0.12), (0.12, 0.0)),
            "l_ij": ((0.0, 0.0), (0.0, 0.0)),
        },
        "options": {"relative_permittivity_formulation": "none"},
        "validity": {"kind": "unknown"},
    },
}
```

Component-indexed fields are tuples and pair fields are square tuples. Each
field has one documented fixed unit. The parser accepts both the established
abbreviation and one descriptive name for each scientific field, normalizes
either spelling immediately, and rejects a payload that supplies both names.
Association topology, correlations, electrolyte inputs, and all dipole,
quadrupole, and dipole-quadrupole inputs remain explicit. The public result is
immutable, and the exposed Helmholtz projection combines the active polar
physics into one `ares_polar` contribution.

The schema is strict and versioned. It rejects unknown fields, missing required
fields, inconsistent vector or matrix sizes, nonfinite values, and ambiguous
alias pairs. It carries no citations or source provenance; those belong to
ePC-SAFT-data. It does require a validity mapping that either declares explicit
bounds or explicitly declares `unknown`.

Use this schema for manufactured invariants, contribution anchors, derivative
checks, density closure, and most SDK numerical parity. Retain actual files only
for `Parameters.from_bundle(path)` and bundle-template I/O tests. The EOS
package has no installed parameter catalog; callers obtain reusable scientific
packets elsewhere and materialize them as a dictionary or bundle path.

### 4. Consolidate the SDK suite around unique contracts

The necessarily unique SDK evidence is smaller than the current list of
callback-specific tests, but it is not zero:

1. one installed-wheel/header C++ consumer compile-link-run test;
2. one table-driven ABI test for capsule name/version, result sizes, field
   offsets, and append-only prefix compatibility;
3. one common fail-closed matrix for null pointers, undersized buffers,
   nonfinite values, and “no writes on failure” across callback shapes;
4. one capability-admission matrix for genuinely distinct model topology,
   required parameter fields, component order, and applicability-domain rules;
5. one numerical matrix per result shape that compares the C callback with the
   Python owner and an independent directional derivative oracle;
6. batch-only scalar/row parity, bounded cancellation, and serialized CppAD
   ownership; and
7. root-tail identity/continuity tests where density or electrolyte-reference
   branch selection is part of the returned contract.

Descriptor presence, header layout, and malformed-buffer tests should be
parameterized over callback descriptors instead of copied for each appended
entry. Scientific cases should vary only when they exercise a different
equation family, topology, derivative order, root structure, or admission rule.

### 5. Merge, move, and retire current families

**Merge in EOS:**

- merge all `test_native_sdk_describes_*` functions into one descriptor matrix;
- merge the separate installed-header-layout tests into one ABI layout matrix;
- merge shared `rejects_invalid`, `fails_closed`, `malformed_buffers`, and
  `without_writes` cases into the common transport-failure matrix;
- merge phase/public parity and derivative checks into callback/result-shape
  matrices with explicit regime IDs; and
- merge hard-chain, dispersion, association, Debye–Hückel, Born, permittivity,
  and future polar cases around generic contribution-sum, zero-reduction,
  permutation, limit, and derivative claims while retaining distinct equations.

**Move after Data/Validation gates exist:**

- move all-species and all-molality literature matrices for Figiel, Held,
  Ascani, and Gross to installed-artifact Validation campaigns;
- move exhaustive canonical bundle contents, provenance, and source-fingerprint
  checks to `ePC-SAFT-data`;
- keep in EOS the smallest independent equation anchor required by each
  distinct contribution or specialized callback; and
- update the ADR admission gates before deleting any matrix they currently
  require. A test is not obsolete merely because the desired future ownership
  differs from the accepted contract today.

**Retire once replaced:**

- retire `test_native_sdk_header_is_installed_package_data` when the installed
  C++ consumer proves header discovery and compilation;
- retire callback-by-callback repetitions of common pointer, size, nonfinite,
  and no-write failures after the generic matrix covers every callback shape;
- retire checks of exact private wording, helper layout, or source path where a
  public behavior or generated consistency command already fails closed; and
- retire duplicated source-number copies that are neither independent oracles
  nor runtime-data result checks.

### 6. Preserve the risks that raw parameters cannot cover

Replacing every source-bound case with inline numbers would remove important
evidence:

- inline tuples do not prove external bundle schema, provenance, Data-release
  pinning, or reusable-packet reproducibility;
- copied canonical numbers can drift from `ePC-SAFT-data` while both the test
  and implementation continue to agree;
- a value copied from production data is not an independent numerical oracle;
  raw cases must be labeled manufactured, analytic, or externally anchored;
  and
- bundle loading and template round trips inherently require real temporary
  bundle artifacts.

The target architecture should therefore use raw tuples for equation and
derivative tests, narrowly scoped bundle-I/O tests, Data-owned canonical packet
tests, and Validation-owned broad literature campaigns. SDK admission becomes
generic over topology, required fields, component order, and declared validity
rather than membership in a paper-specific EOS catalog. That preserves the
simplicity of zmeri/PC-SAFT where it helps while retaining the native,
derivative, and root contracts that make this package materially different.

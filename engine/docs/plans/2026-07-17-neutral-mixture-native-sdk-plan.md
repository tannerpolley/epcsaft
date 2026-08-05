# Neutral Mixture Native SDK Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Append one exact neutral-mixture Helmholtz value/gradient/Hessian entry to `epcsaft.native_sdk.v1` without changing its accepted prefix or adding a second EOS owner.

**Architecture:** The provider remains the sole equation and CppAD owner. A model-bound capsule accepts caller-owned amount and output buffers, checks their sizes against the immutable component count, and evaluates one extensive fixed-temperature phase function in coordinates `(n_1, ..., n_C, V)`.

**Tech Stack:** Python 3.13, C++17, CppAD 20250000.3, pybind11, CMake, pytest, scikit-build-core.

## Global Constraints

- Read `AGENTS.md`, `CONTEXT.md`, `ARCHITECTURE.yaml`, `docs/adr/0001-native-consumer-transport.md`, and the canonical Slice 2 design before editing.
- Run the chemical-engineer preflight with registry checking before equation-sensitive work.
- Preserve the accepted C table prefix byte-for-byte and append fields only.
- Keep one Python package, one extension module, one native target, and Pint as the sole runtime dependency.
- Admit only neutral, nonassociating models; reject association, ions, empty models, and malformed buffers.
- Add no Python API beyond the existing `native_sdk(model)` capsule.
- Use CppAD as the sole nonlinear production derivative authority.
- Do not copy consumer residuals, Ipopt, Ceres, phase identity, density closure, parameter fitting, or third derivatives into this entry.
- Build and test from an installed wheel without sibling-source linkage.

---

### Task 1: Freeze the append-only C contract

**Files:**
- Modify: `src/epcsaft/include/epcsaft/native_sdk_v1.h`
- Modify: `tests/test_native_sdk.py`

**Interfaces:**
- Consumes: accepted `epcsaft_native_sdk_v1` prefix through `evaluate_pure_phase_parameters`.
- Produces: `component_count`, `mixture_result_size`, and `evaluate_mixture_phase` appended to the table.

- [ ] **Step 1: Add a failing ABI test**

Extend the ctypes table with these trailing fields and assert that the original prefix offsets are unchanged:

```python
class MixturePhaseBlockResult(ctypes.Structure):
    _fields_ = (
        ("struct_size", ctypes.c_uint32),
        ("status", ctypes.c_int32),
        ("coordinate_count", ctypes.c_size_t),
        ("gradient_capacity", ctypes.c_size_t),
        ("hessian_capacity", ctypes.c_size_t),
        ("gradient", ctypes.POINTER(ctypes.c_double)),
        ("hessian", ctypes.POINTER(ctypes.c_double)),
        ("helmholtz_over_rt_reference_amount", ctypes.c_double),
        ("pressure_pa", ctypes.c_double),
        ("parameter_fingerprint", ctypes.c_char * 72),
        ("error", ctypes.c_char * 160),
    )

class NativeSdkTable(ctypes.Structure):
    _fields_ = EXISTING_FIELDS + (
        ("component_count", ctypes.c_size_t),
        ("mixture_result_size", ctypes.c_size_t),
        ("evaluate_mixture_phase", ctypes.c_void_p),
    )
```

Create a methane/ethane model and require `component_count == 2`, an exact result size, and a non-null callback. Preserve the existing pure-prefix test.

- [ ] **Step 2: Run the focused test and verify RED**

Run: `python -m pytest tests/test_native_sdk.py -k 'mixture or prefix' -q`

Expected: failure because binary models are rejected and the table tail is absent.

- [ ] **Step 3: Add the exact C declarations**

Append this shape to the public header; keep all existing declarations unchanged:

```c
typedef struct epcsaft_mixture_phase_block_result_v1 {
    uint32_t struct_size;
    int32_t status;
    size_t coordinate_count;
    size_t gradient_capacity;
    size_t hessian_capacity;
    double* gradient;
    double* hessian;
    double helmholtz_over_rt_reference_amount;
    double pressure_pa;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_mixture_phase_block_result_v1;

typedef int (*epcsaft_evaluate_mixture_phase_v1)(
    const void* model_context,
    double temperature_k,
    const double* amounts_mol,
    size_t amount_count,
    double volume_m3,
    epcsaft_mixture_phase_block_result_v1* result
);
```

Append `size_t component_count`, `size_t mixture_result_size`, and the callback to `epcsaft_native_sdk_v1`.

- [ ] **Step 4: Re-run the focused test**

Run: `python -m pytest tests/test_native_sdk.py -k prefix -q`

Expected: the pure prefix test passes; the binary callback test remains RED until Task 3.

- [ ] **Step 5: Commit the frozen ABI contract**

```bash
git add src/epcsaft/include/epcsaft/native_sdk_v1.h tests/test_native_sdk.py
git commit -m "test: freeze neutral mixture SDK tail"
```

### Task 2: Implement the scalar-generic mixture phase block

**Files:**
- Modify: `src/epcsaft/_native/phase_block.hpp`
- Modify: `src/epcsaft/_native/eos.hpp`
- Modify: `src/epcsaft/_native/eos.cpp`
- Modify: `tests/test_native_sdk.py`

**Interfaces:**
- Consumes: `NativeParameters` and the canonical `residual_helmholtz` template.
- Produces: `NativeMixturePhaseBlock { value, gradient, hessian, pressure }` and `NativeModel::evaluate_mixture_phase`.

- [ ] **Step 1: Add failing value and derivative tests**

At the methane/ethane state `T=220 K`, `n=(0.4,0.6) mol`, `V=1.25e-4 m3`, compare the future block with public explicit-density evaluation:

```python
rho = sum(amounts) / volume
x = tuple(value / sum(amounts) for value in amounts)
state = model.evaluate(
    temperature=220.0 * unit_registry.kelvin,
    molar_density=rho * unit_registry.mole / unit_registry.meter**3,
    mole_fractions=x,
)
expected = sum(n * (math.log(n / volume) - 1.0) for n in amounts)
expected += sum(amounts) * state.residual_helmholtz
```

Check value and pressure, then contract the gradient and Hessian with direction `(0.02,-0.01,-2e-6)` against centered differences. Add a reversed-component model test that permutes amounts and derivative rows/columns.

- [ ] **Step 2: Run tests and verify RED**

Run: `python -m pytest tests/test_native_sdk.py -k mixture -q`

Expected: failure because `evaluate_mixture_phase` is not implemented.

- [ ] **Step 3: Add the narrow model support predicate and extensive function**

Implement a predicate requiring at least two components, no association sites, and zero charge for every component. Implement the scalar-generic function:

```cpp
const Scalar total = std::accumulate(amounts.begin(), amounts.end(), Scalar(0.0));
const Scalar density = total / volume;
std::vector<Scalar> composition;
for (const Scalar& amount : amounts) composition.push_back(amount / total);
const auto sigma = segment_sigma(Scalar(temperature), parameters);
const auto diameter = segment_diameter(Scalar(temperature), parameters, sigma);
const auto contributions = residual_helmholtz(
    Scalar(temperature), density, composition, parameters, sigma, diameter,
    {}, Scalar(0.0), Scalar(0.0));
Scalar ideal = Scalar(0.0);
for (const Scalar& amount : amounts) {
    ideal += amount * (scalar_log(amount / volume) - Scalar(1.0));
}
return ideal + total * (contributions.hard_chain + contributions.dispersion);
```

Record one CppAD tape in `C+1` coordinates, then return its value, Jacobian, and Hessian. Set `pressure_pa = -R*T*gradient[V]`. Reject nonpositive amounts/volume and nonfinite tensors.

- [ ] **Step 4: Expose the model method**

Declare and define:

```cpp
bool supports_mixture_native_sdk_v1() const;
NativeMixturePhaseBlock evaluate_mixture_phase(
    double temperature,
    const std::vector<double>& amounts,
    double volume
) const;
```

Do not weaken `supports_native_sdk_v1()` for the accepted pure entries; `native_sdk(model)` admission is changed in Task 3.

- [ ] **Step 5: Run the model-level tests**

Run: `python -m pytest tests/test_native_sdk.py -k mixture -q`

Expected: callback tests remain RED, while direct model implementation compiles under the full extension build.

- [ ] **Step 6: Commit the equation-owner change**

```bash
git add src/epcsaft/_native/phase_block.hpp src/epcsaft/_native/eos.hpp src/epcsaft/_native/eos.cpp tests/test_native_sdk.py
git commit -m "feat: evaluate neutral mixture phase Hessians"
```

### Task 3: Wire caller-owned buffers through the capsule

**Files:**
- Modify: `src/epcsaft/_native/native_sdk.cpp`
- Modify: `tests/test_native_sdk.py`

**Interfaces:**
- Consumes: `NativeModel::evaluate_mixture_phase`.
- Produces: the appended C callback over the existing model-bound capsule.

- [ ] **Step 1: Add malformed-contract tests**

Require explicit failures for null amount pointer, wrong amount count, null gradient/Hessian pointer, short capacities, wrong result size, temperature outside source domain, nonpositive amounts, and a pure or associating model passed to the mixture callback.

- [ ] **Step 2: Run malformed tests and verify RED**

Run: `python -m pytest tests/test_native_sdk.py -k 'mixture and (reject or invalid or malformed)' -q`

Expected: failure because the callback is absent.

- [ ] **Step 3: Append endpoint metadata and callback**

Populate the table tail with the immutable component count, exact result size, and callback. The callback must validate all pointers/capacities before zeroing scalar output fields, preserve caller pointers and capacities, copy exactly `C+1` gradient values and `(C+1)^2` Hessian values, and translate exceptions to the existing status enum.

Change `make_native_sdk_v1` admission to accept either the existing pure predicate or the new neutral-mixture predicate. Pure callbacks invoked on a mixture must return `UNSUPPORTED_MODEL`; the mixture callback invoked on a pure model must do the same.

- [ ] **Step 4: Run the complete provider SDK test**

Run: `python -m pytest tests/test_native_sdk.py -q`

Expected: all SDK tests pass, including old prefix and pure consumer behavior.

- [ ] **Step 5: Commit capsule wiring**

```bash
git add src/epcsaft/_native/native_sdk.cpp tests/test_native_sdk.py
git commit -m "feat: append neutral mixture SDK entry"
```

### Task 4: Installed-artifact proof and capability records

**Files:**
- Modify: `docs/adr/0001-native-consumer-transport.md`
- Modify: `AGENTS.md`
- Modify: `CONTEXT.md`
- Modify: `ARCHITECTURE.yaml`
- Modify: `README.md`
- Create: `docs/science/neutral-mixture-native-sdk.md`

**Interfaces:**
- Consumes: final installed provider wheel and SDK tests.
- Produces: one candidate artifact and an exact documented capability delta for migration review.

- [ ] **Step 1: Update authority-neutral candidate documentation**

Document coordinate order, extensive ideal-plus-residual function, units, buffer ownership, derivative order two, admitted fingerprint, and all exclusions. State that implementation does not become authoritative until its own promotion receipt is accepted.

- [ ] **Step 2: Build and test an isolated wheel**

Run the repository's documented non-editable wheel build, install only that wheel plus test dependencies in an isolated Python 3.13 environment, and run:

```text
python -m pytest tests/test_native_sdk.py -q
python -m pytest -q
```

Expected: all provider tests pass from installed artifacts; import origins resolve under site-packages.

- [ ] **Step 3: Run negative-space and cleanup checks**

Confirm no new Python export, target, runtime dependency, provider source in a consumer, backend selector, density route, third tensor, Ipopt, or Ceres symbol was added. Run `git diff --check`, exact tracked-file review, and the cleanup hook.

- [ ] **Step 4: Commit the candidate records**

```bash
git add AGENTS.md CONTEXT.md ARCHITECTURE.yaml README.md docs/adr/0001-native-consumer-transport.md docs/science/neutral-mixture-native-sdk.md
git commit -m "docs: record neutral mixture SDK candidate"
```

- [ ] **Step 5: Stop for independent review**

Report the final commit, tree, wheel hash, test receipt hash, exact commands, and remaining exclusions. Do not push, promote, or begin equilibrium work from a source checkout.

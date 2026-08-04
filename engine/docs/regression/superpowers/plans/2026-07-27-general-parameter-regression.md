# General Parameter Regression Implementation Plan

> Execute this plan test-first. Keep Provider equations and derivatives in
> Provider, Ceres and regression diagnostics in Regression, and source data in
> caller-owned source-bound contracts.

**Goal:** Deliver the first reusable parameter-family slice: caller-supplied,
source-bound, fixed-composition neutral-binary VLE fitting of one shared
`k_ij`, backed by an installed Provider capability descriptor and exact native
derivatives.

**Architecture:** Provider appends an immutable, model-bound capability
descriptor to its v1 native SDK and remains the only EOS owner. Regression
adds closed typed contracts and one generic public entry point. The Ceres loop
stays in the existing `_native` module and target. The first native contract
has one active `k_ij`, two lifted log-volumes per training row, four residuals
per row, and a `4R x (1 + 2R)` exact Jacobian.

**Tooling:** C++17, Ceres 2.2, CPython 3.13 extension API, Python frozen
dataclasses/enums, pytest, scikit-build-core/CMake.

## Task 1: Provider model-bound capability descriptor

**Provider files**

- Modify: `src/epcsaft/include/epcsaft/native_sdk_v1.h`
- Modify: `src/epcsaft/_native/native_sdk.cpp`
- Modify: `src/epcsaft/_native/native_model.hpp`
- Modify: `src/epcsaft/_native/native_model.cpp`
- Modify: `src/epcsaft/eos.py`
- Test: `tests/test_native_sdk.py`
- Document: `docs/science/neutral-mixture-native-sdk.md`

1. Merge Provider local `main` into its existing
   `codex/regression-provider-integration` worktree and confirm a clean
   baseline.
2. Add failing tests for an append-only v1 descriptor tail, a separate
   nonempty topology fingerprint, exact neutral-binary coordinate order
   `(n1,n2,V,k_ij)`, row-major value/gradient/Hessian layout, canonical
   component mapping, derivative order 2, fixed neutral/nonassociating domain,
   and fail-closed unsupported models.
3. Append closed POD enums and `epcsaft_native_capability_descriptor_v1`, plus
   descriptor count/pointer fields, to `epcsaft_native_sdk_v1`. Preserve every
   existing field offset.
4. Compute a deterministic Provider-owned topology fingerprint from resolved
   component order, charge/site topology, and formulation. Do not reuse the
   parameter fingerprint.
5. Emit only truthful model-bound descriptors. The current neutral-binary
   descriptor is authority-neutral and derivative-ready; it does not claim
   regression validation or admission.
6. Run the focused Provider tests and its native SDK suite. Commit the bounded
   Provider change locally.

## Task 2: Regression closed contracts and capability reader

**Regression files**

- Create: `src/epcsaft_regression/contracts.py`
- Create: `src/epcsaft_regression/parameter_regression.py`
- Modify: `src/epcsaft_regression/__init__.py`
- Test: `tests/test_parameter_regression_contracts.py`

1. Add failing tests for closed parameter-family and observation enums,
   canonical unordered pair identity, duplicate/reversed duplicate rejection,
   SHA-256/source requirements, canonical transformed-dataset identity,
   explicit units/scales/partitions, affine transforms, bounds, and
   primary/confirmation starts.
2. Implement frozen slotted contracts for `SourceDescriptor`,
   `PairParameterIdentity`, `AffineParameterTransform`,
   `ParameterCoordinate`, `FixedCompositionVleObservation`,
   `RegressionProblem`, capability metadata, row diagnostics, parameter
   diagnostics, and the canonical `RegressionResult`.
3. Implement `parameter_capabilities(model)` by reading only the installed
   model capsule descriptor. Reject unknown schema members and incomplete,
   mismatched, or unrecognized metadata; never infer support from callback
   presence.
4. Export only the two intended public functions and the contract types needed
   to construct a problem. Keep legacy public entry points unchanged.
5. Run the focused contract tests and commit the Regression contract layer.

## Task 3: Exact native neutral-binary evaluator

**Regression files**

- Create: `src/epcsaft_regression/native/general_fit.hpp`
- Create: `src/epcsaft_regression/native/general_fit.cpp`
- Modify: `src/epcsaft_regression/native/module.cpp`
- Modify: `CMakeLists.txt`
- Test: `tests/test_parameter_regression_native.py`

1. Add failing native tests for arbitrary row count and order, Provider
   descriptor/fingerprint/component mismatch, exact residual shape, exact
   Jacobian shape, and one directional derivative check.
2. Parse one validated payload into contiguous C++ storage. Use observed
   `x/y`, one mole per phase, caller-declared volume transforms/bounds, and no
   chemistry defaults.
3. For each row form:

   ```text
   r_PL = (P_L - P_obs) / s_P
   r_PV = (P_V - P_obs) / s_P
   r_mu1 = (mu_1,L/RT - mu_1,V/RT) / s_mu1
   r_mu2 = (mu_2,L/RT - mu_2,V/RT) / s_mu2
   ```

   with `P = -RT*Phi_V`, `mu_i/RT = Phi_ni`, lifted
   `V = V_origin*exp(u)`, and `k_ij = origin + scale*z`.
4. Assemble the exact Jacobian from Provider's row-major Hessian in
   `(n1,n2,V,k_ij)`. Use `Phi_VV`, `Phi_Vk`, `Phi_ni,V`, and `Phi_ni,k`;
   no third derivatives or numerical differencing.
5. Add private `_native` evaluate/solve entry points to the existing module
   and target. Do not add a second extension, registry, or backend selector.
6. Run focused native tests and commit.

## Task 4: Canonical Ceres solve and diagnostics

**Regression files**

- Modify: `src/epcsaft_regression/native/general_fit.cpp`
- Modify: `src/epcsaft_regression/parameter_regression.py`
- Test: `tests/test_parameter_regression.py`

1. Add failing tests for deterministic primary and confirmation starts,
   full/projected rank, conditioning, active bounds, incomplete-row failure,
   and exact training/held-out/stress accounting.
2. Solve one parameter block plus two nuisance log-volume coordinates per
   training row using Ceres. Reporting rows are evaluated only after fitting.
3. Compute singular values and ranks using
   `tau = 100*epsilon*max(M,N)*sigma_max`; project parameter columns against
   the nuisance-column span before computing parameter rank.
4. Return immutable diagnostics with separate solver, numerical, workflow,
   physical/scientific, and predictive statuses. A local successful fit is
   authority-neutral and predictive status remains unadjudicated without an
   approved held-out cutoff.
5. Verify that mutating held-out data cannot change fitted values or training
   cost. Fail closed on partition overlap or missing exact derivative paths.
6. Run focused tests and commit.

## Task 5: Installed-artifact reference and performance evidence

**Regression files**

- Modify only if needed: `tests/test_parameter_regression.py`
- Create: `docs/reviews/general-parameter-regression-first-slice.md`
- Modify: `CONTEXT.md`
- Modify: `ARCHITECTURE.yaml`
- Modify: `docs/science/general-parameter-regression.md`

1. Build one exact Provider wheel from Task 1 and install it into an isolated
   environment; build Regression against that installed public header and
   artifact only.
2. Run a bounded source-backed methane/ethane fixed-composition campaign as
   reference evidence, not as a mixture restriction or predictive claim.
   Preserve retained historical blocked evidence and do not tune its rows,
   scales, bounds, or tolerances.
3. Record exact derivative checks, row accounting, `4R x (1+2R)` dimensions,
   ranks, conditioning, active-bound status, start agreement, artifacts, and
   commands. A converged reference campaign does not admit a capability.
4. Demonstrate that Python Provider methods are not invoked inside native
   residual evaluation after payload preparation. Run a separate bounded
   benchmark against the legacy native path; do not add wall-clock assertions
   to unit tests.
5. Reconcile the status owners to distinguish represented,
   derivative-ready, fit-ready, reference-validated, and admitted.

## Task 6: Whole-repository verification and review

1. Run focused Provider tests, then its full isolated installed-wheel suite.
2. Run focused Regression tests, then its full isolated installed-wheel suite.
3. Audit wheel contents, imports/linkage, one native target/module, descriptor
   offsets, and negative space.
4. Run independent scientific/correctness review and address actionable
   findings.
5. Apply cutthroat cleanup and minimum-surface review; remove obsolete
   planning or compatibility surface introduced by the implementation.
6. Run repository cleanup hooks in both repositories, verify clean status,
   commit final documentation/evidence, and prepare a non-draft PR only after
   the local subjects are green. Do not push or publish without explicit
   authorization.

## Explicit exclusions

- no `l_ij`, association, `k_hb_ij`, polar, dielectric, reactive, uncertainty,
  global-identifiability, or catalog-persistence implementation in this slice;
- no copied EOS or equilibrium equations;
- no numerical production derivative backend;
- no arbitrary residual plugins or mutable registry;
- no hidden data, bounds, scales, starts, volume policy, or chemistry defaults;
- no paper-named general API and no predictive or downstream-readiness claim.

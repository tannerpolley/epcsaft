# Model-Bound Aqueous Active-`k_ij` MIAC Design

Status: local implementation candidate pending bounded independent subagent
review; no retained installed artifact, Regression fit, Validation campaign,
or authority exists.

Project Truss outcome: `ePC-SAFT/ePC-SAFT-regression#1`; Provider leaf:
`ePC-SAFT/ePC-SAFT#1`. The Provider implementation base is
`43ebe47366bffff9f595c4956a6a024967a5c6ee`.

## Outcome and ownership boundary

The minimum Provider seam for the approved aqueous Figiel interaction recovery
is one model-bound native callback that evaluates the existing molality-scale
stoichiometric mean ionic activity coefficient while replacing exactly three
binary dispersion interactions for one water/cation/anion model. It returns
`ln(gamma_pm^m)` and its exact total first derivatives with respect to those
three trial interactions.

Provider remains the sole owner of the ePC-SAFT Helmholtz equation, association
state, pressure density closure, reference-limit sequence, component ordering,
parameter resolution, and CppAD derivative propagation. Regression owns the
164 residuals, global 11-parameter mapping, bounds, scaling, starts, Ceres,
rank/conditioning adjudication, fitted values, and result. Validation owns the
installed-artifact campaign. This callback does not admit a parameter fit or
write a trial value into the Provider catalog.

The design appends one result structure, one coordinate enum, one result-size
field, and one function pointer to `epcsaft.native_sdk.v1`. It adds no Python
public API, native module or target, generic active-parameter descriptor,
parameter overlay, density solver, reference-state implementation, numerical
derivative backend, or compatibility route.

## Source and admitted domain

The parameter classification comes from Figiel, Yu, and Held (2025), DOI
`10.1021/acs.iecr.5c00475`, Tables 4 and 5 and the **Model Parameters**
discussion. The verified 71,826-byte Markdown artifact has SHA-256
`ce80533925a91bc59d8d0d8056113c40611ca26c2edf04aced76986d50bd4bae`.
The experimental targets come from Hamer and Wu (1972), DOI
`10.1063/1.3253108`, through approved Validation packet
`8944d34f7002cda1bb8760e606cc1f11696f58cd`. Its retained 164-row CSV has
SHA-256 `2f63e13f06a5b0f4e8bca2980b6a8d9d7fb0f839153c43e3a71952daf9796595`.

The callback admits only:

- `T = 298.15 K` and `P = 100000 Pa`;
- liquid water as component 0, one monovalent cation as component 1, and one
  monovalent anion as component 2;
- formula-unit molality from `0.001` through `6 mol kg^-1` of salt-free water;
- the six 1:1 salts LiCl, NaCl, KCl, LiBr, NaBr, and KBr; and
- the exact current `figiel-2025-reference-electrolytes/1` inactive catalog,
  including water `f_solv = 1.5`, current Table 3 ion and Born values, and the
  current dielectric treatment.

The six admitted ordered models and immutable source fingerprints are:

| Salt | Ordered components | Source fingerprint |
| --- | --- | --- |
| LiCl | `(water, lithium-cation, chloride-anion)` | `sha256:1bb528ebe8f5612757e148608fc55821f9fb03737dbcec6d0bc4fffd0f4cbc4c` |
| NaCl | `(water, sodium-cation, chloride-anion)` | `sha256:7c637771bc9f717b8f47b44bb2a61044c3fe83084dca7c3c16102fba0989912d` |
| KCl | `(water, potassium-cation, chloride-anion)` | `sha256:d29cef0c0f63034436d547d0aafa57934effe06783c8dffd89c94fa85e6940f4` |
| LiBr | `(water, lithium-cation, bromide-anion)` | `sha256:1258d8f76ddb4300d4983e70a2efa107b5cdc468748fe8760a04f82c4de38668` |
| NaBr | `(water, sodium-cation, bromide-anion)` | `sha256:b46a48ee5025e15a2d063b6f3ca4de54e29c33d429b94940913a9d225c78456b` |
| KBr | `(water, potassium-cation, bromide-anion)` | `sha256:1d498278767a568e2120286f27f98a0b79268c92a68805d85ca30a21f14f3284` |

Admission requires the fingerprint and the structural invariants: three
components in the exact order above, charges `(0,+1,-1)`, `ionic-ssm-ds`
topology, water as the sole solvent, and the Figiel source domain. A trial
interaction does not alter or replace the returned source fingerprint.

The callback excludes I-, H+, sulfate, expanded ions, multiple salts,
methanol, ethanol, organic ion-solvent interactions, solvent-solvent
interactions, other temperatures or pressures, molalities outside the retained
source overlap, osmotic targets, Born diameters, solvation factors, dielectric
parameters, neutral parameters, association parameters, and simultaneous
Tables 2--5 recovery.

## Exact observable and coordinates

For formula-unit molality `m`, let the existing Provider reference owner return
the infinite-dilution ion log fugacity coefficients and let the finite liquid
state return the corresponding values at `m`. For ion `i`:

```text
ln(gamma_i^*) = ln(phi_i(T,P,m;k)) - ln(phi_i^infinity(T,P;k)).
```

All admitted salts are 1:1. Figiel Eqs. (20)--(22) therefore give

```text
ln(gamma_pm^m)
  = 1/2 [ln(gamma_cation^*) + ln(gamma_anion^*)]
    - ln(1 + 2 M_water m),
```

where `M_water` is the resolved salt-free water molar mass in `kg mol^-1` and
`m` is in `mol kg^-1`. The callback returns the logarithm, which is
dimensionless. It does not exponentiate inside the Ceres loop.

The three dimensionless active coordinates are immutable and symmetric in the
Provider parameter matrix:

```text
0: k_water_cation = k_cation_water
1: k_water_anion   = k_anion_water
2: k_cation_anion  = k_anion_cation.
```

Every coordinate must be finite. Provider does not own Regression's `[-1,1]`
bounds, parameter scale, starts, or acceptance gates. It evaluates the trial
coordinates or returns a structured domain/evaluation failure.

The composition is generated only by
`epcsaft::detail::formula_composition` for `(water,cation,anion)` and
stoichiometry `(1,1)`. The infinite-dilution reference is the existing
`(1e-6, 1e-8, 1e-10, 1e-12) mol kg^-1` sequence. The terminal value is accepted
only under the existing `5e-5` maximum ion-log-fugacity change criterion.

## Exact derivative closure

Let `a^r(T,rho,x,k)` be the one canonical dimensionless residual molar
Helmholtz function with the three trial interactions substituted only into the
six symmetric unlike-pair matrix cells. With

```text
Z = 1 + rho a^r_rho
P = rho R T Z
ln(phi_i) = mu_i^r/(R T) - ln(Z),
mu_i^r/(R T)
  = a^r + rho a^r_rho + a^r_xi - sum_j x_j a^r_xj,
```

the density root is solved outside the tape by the existing Provider density
owner. For each active coordinate `k_alpha`, its exact fixed-`T,P,x`
sensitivity is

```text
rho_kalpha = -P_kalpha / P_rho,

d ln(phi_i)/d kalpha |T,P,x
  = partial ln(phi_i)/partial kalpha |rho
    + partial ln(phi_i)/partial rho |k * rho_kalpha.
```

The required pressure and fugacity partials use the gradient and reduced
Hessian of `a^r` in `(T,rho,x_water,x_cation,x_anion,k0,k1,k2)`. The current
association site fractions remain eliminated states in the same derivative
owner: their first and second implicit sensitivities are obtained from the
mass-action Jacobian exactly as they are for existing thermodynamic
derivatives. The three `k_ij` coordinates affect the canonical dispersion
term; the implementation must parameterize the existing templated parameter
view rather than copy the dispersion equation.

At the finite and terminal reference states the callback forms

```text
d ln(gamma_pm^m)/d kalpha
  = 1/2 sum_i in {cation,anion}
      [d ln(phi_i(T,P,m;k))/d kalpha
       - d ln(phi_i^infinity(T,P;k))/d kalpha].
```

The molality conversion term is independent of `k_alpha`. CppAD value,
gradient, and Hessian plus the first implicit density solve therefore produce
the exact three-column callback Jacobian. No third derivative is required.

The reference derivative sequence is accepted only when all six terminal
ion/coordinate derivatives are finite and

```text
max_i,alpha |dlnphi_i/dk_alpha(1e-12) - dlnphi_i/dk_alpha(1e-10)|
  <= 5e-5 * max(1, max_i,alpha |dlnphi_i/dk_alpha(1e-12)|).
```

The callback reports the unscaled maximum difference as a diagnostic.

## Native ABI candidate

Version 1 appends the following coordinate enum, result, result-size field,
and callback after the current active-Born tail:

```c
#define EPCSAFT_AQUEOUS_ACTIVE_KIJ_COORDINATE_COUNT_V1 3u

enum epcsaft_aqueous_active_kij_coordinate_v1 {
    EPCSAFT_AQUEOUS_KIJ_WATER_CATION_V1 = 0,
    EPCSAFT_AQUEOUS_KIJ_WATER_ANION_V1 = 1,
    EPCSAFT_AQUEOUS_KIJ_CATION_ANION_V1 = 2
};

typedef struct epcsaft_aqueous_miac_kij_result_v1 {
    uint32_t struct_size;
    int32_t status;
    double log_mean_ionic_activity_coefficient_molality;
    double derivative[EPCSAFT_AQUEOUS_ACTIVE_KIJ_COORDINATE_COUNT_V1];
    double reference_molality_mol_per_kg;
    double reference_convergence_error;
    double reference_derivative_convergence_error;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_aqueous_miac_kij_result_v1;

typedef int (*epcsaft_evaluate_aqueous_miac_kij_v1)(
    const void* model_context,
    const char* expected_parameter_fingerprint,
    double temperature_k,
    double pressure_pa,
    double formula_unit_molality_mol_per_kg,
    const double* k_ij,
    size_t k_ij_count,
    epcsaft_aqueous_miac_kij_result_v1* result
);

/* appended to epcsaft_native_sdk_v1 */
size_t aqueous_miac_kij_result_size;
epcsaft_evaluate_aqueous_miac_kij_v1 evaluate_aqueous_miac_kij;
```

On the current 64-bit layout, the existing 168-byte active-Born table prefix
is unchanged. The appended result-size field is at offset 168, the callback is
at offset 176, and the full candidate table is 184 bytes. The result is 296
bytes: scalar fields begin at offsets 8, 16, 40, 48, and 56; the fingerprint
begins at 64 and the error buffer at 136. These are design inspection values;
implementation tests must prove them with installed-header `offsetof` and
`sizeof` checks.

The expected fingerprint input is mandatory and must equal the bound source
fingerprint byte-for-byte. The existing table-owned component IDs and charges
declare component order. The fixed enum declares active-coordinate order. No
second parameter-order representation is added.

## Failure and status semantics

The callback follows the current fail-closed SDK rule:

- null context, fingerprint, coordinate pointer, or result is
  `INVALID_ARGUMENT`;
- result-size mismatch, coordinate count other than three, or an undersized
  table is `ABI_MISMATCH`; a result-size mismatch does not write the untrusted
  result layout;
- a nonfinite state, molality, or trial coordinate is `INVALID_ARGUMENT`;
- a fingerprint, component order, charge order, topology, salt, or source
  model outside the six-model allowlist is `UNSUPPORTED_MODEL`;
- a finite state other than exactly 298.15 K and 100000 Pa, molality outside
  `[0.001,6] mol kg^-1`, source ion-fraction violation, unavailable liquid
  branch, unstable/ill-conditioned pressure root, or close-packing violation
  is `DOMAIN_ERROR`; and
- association, density, reference-limit, derivative-limit, fugacity, or other
  numerical failure is `EVALUATION_ERROR`.

After exact-size acceptance, every failure returns the same code in the
function return and `result.status`, writes a bounded deterministic error,
sets all scalar values and derivatives to quiet NaN, and clears the
fingerprint. Success returns/statuses `OK`, clears the error, reports the exact
source fingerprint and reference diagnostics, and returns finite value and
derivatives. No failure silently drops a row or substitutes catalog values.

Solver convergence belongs to the density/reference owners; callback numerical
validity requires finite exact derivatives and the reference derivative gate;
physical validity requires the admitted model/domain and certified liquid
roots. Regression solver and scientific acceptance statuses remain separate.

## Evidence gate before implementation handoff

The implementation may begin only after a bounded independent subagent review
approves this exact stable design subject. That approval opens only the
separate implementation leaf and cannot promote runtime, catalog, regression,
scientific, or publication authority. The minimum retained Provider proof is
then:

1. **Source and identity.** Installed-header tests bind the six fingerprints,
   exact component/charge order, coordinate enum, state, molality domain,
   callback/result sizes, offsets, and prefix preservation. Unsupported
   fingerprints, reordered or synthetic models, other salts/solvents, wrong
   states, bad counts, and nonfinite inputs reject.
2. **Value parity.** At the source three-parameter tuple, representative low,
   middle, and high molalities for all six salts reproduce the current public
   `EPCSAFT.reference_state(...).at_molality(...)` logarithmic MIAC within the
   existing owner path's binary64 resolution. This is owner parity, not an
   independent physical oracle.
3. **Exact derivatives.** At one nonpublished interior triple, at least LiCl
   and KBr cover low and high molality. Each returned column is compared with
   complete-callback centered differences at `h=1e-4` and `h/2=5e-5`:

   ```text
   |J_exact - J_h/2|
     <= max(1e-8, 20 |J_h - J_h/2|, 2e-8 |J_exact|).
   ```

   Finite differences are test-only oracles. They are never a production
   backend.
4. **Derivative ownership.** A white-box fixed-state check proves the active
   dispersion coordinates enter the one Helmholtz tape symmetrically and that
   the returned total derivative includes pressure-density and reference-limit
   terms. Association implicit solves and density conditioning must fail
   closed. No copied dispersion, density, or reference equation may appear in
   the SDK transport owner.
5. **Scientific anchors.** The existing 164-row published-tuple campaign
   remains Validation evidence. Provider package tests retain only compact
   value-parity and derivative cases; they do not duplicate the Hamer--Wu
   training campaign or Regression rank/fit acceptance.
6. **Installed artifact.** A normal wheel and public header are retained with
   commit, tree, SHA-256, import/linkage, prefix, wheel-content, negative-space,
   full-suite, cleanup, and clean-status evidence. No Provider source checkout
   is consumed by Regression or Validation.

The later Regression preflight owns the complete 164x11 residual/Jacobian,
rank 11, singular values, conditioning, sensitivity correlations, and
three-start fit. A Provider pass cannot establish those claims.

## Canonical implementation owners and surface budget

An approved implementation changes only the existing owners needed by the
callback:

- `src/epcsaft/_native/helmholtz.hpp`: unchanged equation owner, parameterized
  through its existing templated parameter access;
- `src/epcsaft/_native/derivatives.hpp`: three active interaction coordinates,
  reduced association derivatives, and state partials;
- `src/epcsaft/_native/electrolyte_reference.hpp`: the one existing
  composition/reference schedule, generalized privately to carry the six
  ion/coordinate derivative values without a second sequence;
- `src/epcsaft/_native/eos.hpp` and `eos.cpp`: one closed native MIAC result and
  fixed-pressure total-derivative assembly;
- `src/epcsaft/_native/native_sdk.cpp` and
  `src/epcsaft/include/epcsaft/native_sdk_v1.h`: one append-only transport
  result and callback;
- existing `eos.py`/native SDK admission sets: the exact six fingerprints only;
  and
- existing Provider native SDK/package tests: at most three new persistent
  parameterized claim owners plus installed-artifact verification.

No new production source file, Python public result, dependency, native target,
module, registry, overlay, compatibility shim, or catalog record is justified.
If implementation cannot remain within those owners, it returns for design
review rather than widening surface.

## Batch throughput amendment

Project Truss Provider issue `ePC-SAFT/ePC-SAFT#29` adds two typed append-only
tails to the same 64-bit SDK table:

- offset 224: active aqueous `k_ij` over an ordered molality array;
- offset 232: active water-solvation-factor over an ordered molality array;
- complete table size: 240 bytes.

Both callbacks reuse the existing scalar result records. They require equal
positive molality/result counts, preserve input order, initialize every
correctly size-declared result fail-closed, and apply the existing model,
fingerprint, state, and `0.001-6 mol/kg` gates. If any result size is unknown,
the entire result array remains untouched because Provider cannot safely write
an unknown ABI record. The active model and derivative reference sequence are
constructed once per salt and parameter vector. At the exact catalog
solvation-factor coordinate, the pre-existing public-value parity branch also
runs its value reference once per batch, never once per row. Finite-row values
and exact total fixed-pressure derivatives then execute through bounded
Provider-owned CppAD concurrency; scalar callbacks delegate to the same
calculation owner with a one-row array.

The consumer must require the 240-byte tail and may not fall back to scalar
row calls. Installed evidence must establish bitwise scalar/batch parity over
the frozen 21-row and 164-row inputs and a bounded Stage-C throughput result.
This amendment changes transport cost only: equations, rows, derivatives,
diagnostics, bounds, targets, fit policy, and authority remain unchanged.

## Cooperative evaluation-budget tail

The Provider appends one 24-byte evaluation-budget descriptor and two typed
bounded-batch callbacks after the 240-byte throughput prefix, making the
current complete table 264 bytes on the supported 64-bit ABI. Existing
consumers that require only the 240-byte prefix remain valid and retain the
original batch behavior. The descriptor carries an atomic work budget counted
at safe reference-sequence and row checkpoints plus an absolute monotonic
deadline in nanoseconds; `UINT64_MAX` and zero disable the respective limits.

The bounded callbacks return the existing result records with a uniform
`CANCELLED` or `DEADLINE_EXCEEDED` status when a checkpoint rejects further
work. Provider-owned CppAD workers stop claiming new rows after the first
checkpoint exception and are joined before the callback returns. No equation,
derivative, row ordering, fingerprint, model admission, scalar callback, or
unbounded batch contract changes.

## Alternatives and decision

| Alternative | Consequence | Decision |
| --- | --- | --- |
| Generic active-parameter descriptors or a mutable overlay | Adds parameter identity, catalog mutation, validation, and lifetime policy far beyond one source-backed family | Reject |
| Expose active electrolyte phase tensors and make Regression repeat density closure and the reference sequence | Smaller-looking Provider callback, but duplicates Provider ownership and requires implicit branch logic in Regression | Reject |
| Append one model-bound aqueous MIAC value/derivative callback | Exact three-column Ceres input, one existing reference/density owner, six fingerprints, no generic surface | Approve for review |

The current handoff status is
`PROVIDER_AQUEOUS_ACTIVE_KIJ_IMPLEMENTATION_CANDIDATE_PASSED_LOCAL_GATES_INDEPENDENT_REVIEW_PENDING`.
The callback remains a local candidate until bounded independent subagent
review approves an exact stable commit; the approval changes sequencing only
and transfers no scientific or production authority.

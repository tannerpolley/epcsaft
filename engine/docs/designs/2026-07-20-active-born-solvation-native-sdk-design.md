# Model-Bound Active Born-Diameter Solvation Design

Status: provider design candidate; implementation and authority not approved

Migration input: D-023 at migration commit
`2ab33646c9ef6069670fc37279dc19e0e7bc69ea`, tree
`dcb96f150e59c32e7453fa5eadf2b6102f39e715`

## Outcome and boundary

The minimum Provider seam for the first Figiel regression tracer is one
model-bound native callback that evaluates the existing single-ion
mole-fraction-reference solvation Gibbs energy in water while replacing only
the selected ion's Born diameter with one positive finite trial value. It
returns the scalar value and its first total derivative with respect to that
trial diameter. No second derivative is consumed by the proposed least-squares
residual and none belongs in the interface.

This is a design-only, authority-neutral candidate. It does not admit a target,
residual definition, scale, weight, bound, Ceres object, fit result, catalog
write, density solver, generic active-parameter overlay, or additional ion.
Implementation must wait for permanent-lab review and Validation's immutable
Figiel main/SI target packet.

The read-only primary-source cross-check is Figiel, Yu, and Held (2025), DOI
`10.1021/acs.iecr.5c00475`, Eq. 17 and the **Model Parameters** discussion. The
verified Markdown artifact has SHA-256
`ce80533925a91bc59d8d0d8056113c40611ca26c2edf04aced76986d50bd4bae`.
That discussion classifies the ion Born diameters as adjusted to
infinite-dilution water solvation Gibbs energies, while solvent `f_solv` and
ion/water or ion/ion `k_ij` use separate MIAC fits. Those parameters, neutral
solvent parameters, ordinary ion segment parameters, dielectric parameters,
and V species therefore do not enter this seam. Local source paths are
execution metadata and are not part of the contract.

## Exact observable

For active ion `i` in the water reference state `S`, the callback evaluates the
same x-process quantity already owned by `EPCSAFT.reference_state`:

```text
Delta_solv,x G_i^infinity(d_i^Born)
    = R T ln(phi_i^infinity,S(d_i^Born))                    [J mol^-1]
```

This is Figiel--Yu--Held Eq. 17 under the Provider's existing path-defined
single-ion convention. It is not a gas-phase absolute standard chemical
potential and it is not a stoichiometric mean-ion observable. The reference
path is exactly:

- `T = 298.15 K` and `P = 100000 Pa` (`1 bar`);
- liquid water with salt-free solvent mole fractions `(1, 0, 0)` in the bound
  component order;
- one electroneutral 1:1 formula unit with the active monovalent ion and its
  fixed catalog counterion;
- formula-unit molalities `(1e-6, 1e-8, 1e-10, 1e-12) mol kg^-1`;
- the terminal `1e-12 mol kg^-1` value, accepted only when the maximum change
  in either ion's `ln(phi)` from `1e-10 mol kg^-1` is at most `5e-5`; and
- `R = 8.31446261815324 J mol^-1 K^-1`, matching the existing reference owner.

The callback takes `temperature_k`, `pressure_pa`, and
`born_diameter_angstrom` explicitly. The first two must equal the source-domain
state above; they are not active coordinates and accepting other values would
silently broaden this tracer. The sole active coordinate is
`d_i^Born [angstrom]`. The reported derivative is therefore

```text
d Delta_solv,x G_i^infinity / d d_i^Born
    at fixed T, P, reference composition path, and all other source records
                                                        [J mol^-1 angstrom^-1].
```

At each reference molality, the composition is fixed by the salt-free water
molar mass, the two integer charges, and electroneutral 1:1 stoichiometry.
Every other input is immutable model state: ordered component identity; molar
mass; segment count; temperature-dependent segment diameter; dispersion
energy; all `k_ij`; water association sites, energy, and volume; charge;
solvent relative-permittivity correlation; dielectric ion-suppression
coefficient; ionic-region relative permittivity; every solvation factor; the
counterion Born diameter; source-domain limits; and the Provider physical
constants. The `0.38` total-ion mole-fraction source cap remains enforced.
Only the active ion's source Born diameter is substituted before the canonical
Helmholtz owner is recorded.

## Single reference-sequence owner

A future implementation must first move the existing reference-sequence
algorithm behind one private native owner:
`epcsaft::detail::electrolyte_reference_sequence` in
`src/epcsaft/_native/electrolyte_reference.hpp`. That owner alone defines the
electroneutral formula-unit composition construction, the exact four-molality
schedule, the `5e-5` convergence adjudication, and selection of the terminal
ion fugacity values and diagnostics. Its internal scalar/evaluator parameter
may carry either the immutable source model or the one active Born diameter;
it is not a generic active-parameter API.

Delegation is one-way:

1. The existing public `EPCSAFT.reference_state` path validates its current
   inputs and delegates through the existing private `epcsaft._core` module to
   this native owner. `src/epcsaft/reference.py` retains public units, immutable
   result construction, activity/transfer behavior, and error translation, but
   no longer owns the schedule, composition equation, convergence decision, or
   terminal selection. Its finite-molality activity path uses the same owner's
   internal composition operation rather than retaining the Python formula.
2. The future `evaluate_ion_solvation_born` callback calls that owner directly
   with the active-diameter evaluator and reads the same terminal selection and
   diagnostics. It does not reproduce the loop in `native_sdk.cpp`.

The refactor changes no public Python signature, result type, units, schedule,
tolerance, numerical convention, or failure behavior. It adds no public
reference abstraction and no second implementation; the private `_core`
binding is only the delegation route to the one Provider-internal owner.

## Fixed models, identity, and order

The endpoint admits exactly these selections from
`figiel-2025-reference-electrolytes/1`. The active ion is always component 1;
component 0 is water; component 2 is the fixed counterion. The returned
fingerprint remains the immutable source-model fingerprint. A trial coordinate
never masquerades as a fitted catalog record.

| Active ion | Charge | Source `d_i^Born` [angstrom] | Bound component order | Source fingerprint |
|---|---:|---:|---|---|
| `lithium-cation` | +1 | 2.784 | `(water, lithium-cation, chloride-anion)` | `sha256:1bb528ebe8f5612757e148608fc55821f9fb03737dbcec6d0bc4fffd0f4cbc4c` |
| `sodium-cation` | +1 | 3.445 | `(water, sodium-cation, chloride-anion)` | `sha256:7c637771bc9f717b8f47b44bb2a61044c3fe83084dca7c3c16102fba0989912d` |
| `potassium-cation` | +1 | 4.150 | `(water, potassium-cation, chloride-anion)` | `sha256:d29cef0c0f63034436d547d0aafa57934effe06783c8dffd89c94fa85e6940f4` |
| `chloride-anion` | -1 | 4.100 | `(water, chloride-anion, sodium-cation)` | `sha256:7551f1eee5903b66061cf7520f3bb59b169896ce372f3df3d48aa7ec778c39d4` |
| `bromide-anion` | -1 | 4.480 | `(water, bromide-anion, sodium-cation)` | `sha256:70ae04599dfa8338175e793bac6b9e4dfad37a9b96a568b5484dc87f104ef1a9` |

Admission requires both the exact fingerprint and the structural invariants:
three components in the listed order, `ionic-ssm-ds` topology, water as the
only neutral solvent, charges `(0,+1,-1)` or `(0,-1,+1)`, and the existing
Figiel source domain. This defense in depth rejects reordered selections,
synthetic record copies, mixtures with another solvent or counterion, and
models containing more ions. In particular H+, I-, sulfate, transition-metal
ions, and vanadium species return `UNSUPPORTED_MODEL`; absence from the
current catalog is not permission to accept a structurally similar model.

## Derivative closure

Let `a^r(T,rho,x,d)` be the dimensionless residual molar Helmholtz function,
with `d = d_i^Born`, and define

```text
Z = 1 + rho a^r_rho
P = rho R T Z
mu_i^r/(R T) = a^r + Z - 1 + a^r_xi - sum_j x_j a^r_xj
ln(phi_i) = mu_i^r/(R T) - ln(Z).
```

The pressure root remains Provider-owned and outside the CppAD tape. At fixed
`T`, `P`, and reference composition, its exact implicit sensitivity is

```text
rho_d = -P_d / P_rho

d Delta_solv,x G_i^infinity/dd
    = R T [(partial ln(phi_i)/partial d)_rho
           + (partial ln(phi_i)/partial rho)_d rho_d].
```

`P_d`, `P_rho`, and the two fugacity partials require only the Hessian of the
single canonical Helmholtz function in density, composition, and `d`. No third
derivative is involved. The existing association solution remains part of
that EOS evaluation; it is not a second derivative owner. The callback's
terminal derivative is the derivative of the accepted finite
`1e-12 mol kg^-1` approximation. Admission evidence must also show that the
derivative is stable between the final two reference molalities:

```text
abs(g'_1e-12 - g'_1e-10)
    <= 5e-5 * max(1 J mol^-1 angstrom^-1, abs(g'_1e-12)).
```

For a Regression-owned residual such as
`r(d) = [G(d) - G_target] / scale`, Ceres consumes `r` and
`dr/dd = G'(d)/scale`. Targets and scaling do not change the Provider
derivative order. A scalar second derivative would be unused and is excluded.

## Native ABI candidate

Version 1 can preserve every existing prefix by appending a result-size field
and one callback after `evaluate_mixture_phase_kij`:

```c
typedef struct epcsaft_ion_solvation_born_result_v1 {
    uint32_t struct_size;
    int32_t status;
    double solvation_gibbs_j_per_mol;
    double derivative_j_per_mol_per_angstrom;
    double reference_molality_mol_per_kg;
    double reference_convergence_error;
    char parameter_fingerprint[EPCSAFT_NATIVE_SDK_V1_FINGERPRINT_SIZE];
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE];
} epcsaft_ion_solvation_born_result_v1;

typedef int (*epcsaft_evaluate_ion_solvation_born_v1)(
    const void* model_context,
    double temperature_k,
    double pressure_pa,
    double born_diameter_angstrom,
    epcsaft_ion_solvation_born_result_v1* result
);

/* appended to epcsaft_native_sdk_v1 */
size_t ion_solvation_born_result_size;
epcsaft_evaluate_ion_solvation_born_v1 evaluate_ion_solvation_born;
```

On the current 64-bit ABI, the accepted pure and local mixture prefixes remain
unchanged, the current full table remains the 88-byte prefix, the appended
result-size field is at offset 88, the callback is at offset 96, and the new
full table is 104 bytes. The proposed result is 272 bytes. These are candidate
inspection values, not authority; implementation tests must prove them using
`offsetof` and `sizeof` against the installed header.

One semantically exact result type is necessary because every existing SDK
result represents an extensive phase Helmholtz block and its coordinate
tensors. Reusing one would falsely label a single-ion reference-state Gibbs
energy as `A/(RT)`. There is no second representation of this result.

The callback follows the existing fail-closed ABI rule:

- a null result pointer returns `INVALID_ARGUMENT` and cannot write a result;
- a result `struct_size` mismatch returns `ABI_MISMATCH` without writing the
  untrusted result layout;
- after exact-size acceptance, every failure returns the same code in both the
  function return and `result.status`, writes a bounded deterministic error,
  leaves every scalar observable as a quiet NaN, clears the fingerprint, and
  does not change table/context ownership;
- success returns/statuses `OK`, clears `error`, reports the terminal reference
  molality and convergence error, and copies the exact source fingerprint; and
- finite positive trial diameter is Provider domain validation, while any fit
  bounds or transform remain Regression policy.

Null context and nonfinite or nonpositive `T`, `P`, or diameter are
`INVALID_ARGUMENT`. An unlisted bound model is `UNSUPPORTED_MODEL`. A finite
state other than exactly 298.15 K and 100000 Pa, violation of the source ion
fraction cap, or absence of a valid mechanically stable density root is
`DOMAIN_ERROR`. Unavailable/nonfinite fugacity, reference-limit
nonconvergence, or an unexpected numerical failure is `EVALUATION_ERROR`.
The capsule and opaque context retain the existing Python lifetime rule.
Existing callbacks remain unsupported for these ionic models; admitting this
exact tail does not broaden their contracts.

## Alternatives

| Shape | Surface | Decision |
|---|---|---|
| Add a trial diameter to the Python `reference_state` API and call Python from Regression | No C ABI change, but Python remains in the Ceres loop and fitted-candidate semantics leak into the public reference object | Reject |
| Add a generic native active-parameter descriptor/overlay | One reusable mechanism, but creates parameter registry, identity, validation, and result semantics far beyond one tracer | Reject |
| Append one model-bound solvation callback and one exact scalar result | One active scalar, one observable, one derivative, immutable five-model allowlist, existing EOS/density/reference owners | Recommend |

Direct output pointers were also considered. They save one small structure but
lose size-gated failure semantics, fingerprint transport, and the reference
limit diagnostics needed to audit the single-ion convention. That reduction
is not maintainable minimality.

## Implementation evidence gate

Before runtime implementation can be reviewed, the smallest sufficient proof
set is:

1. At each of the five source Born diameters, callback value and reference
   diagnostics agree with the public `EPCSAFT.reference_state` owner at the
   exact 298.15 K, 1 bar water point. This is owner parity, not an independent
   physical oracle.
2. Once Validation freezes SI Table S5, at least one cation and one anion
   source point agree with the independently retained value, units, sign, and
   x-process convention. The broad target ledger remains Validation-owned.
3. A test-only high-precision implementation of the published corrected
   two-region Born expression checks the active target-ion contribution at a
   real dilute water state. It may share constants and test fixtures, but not
   production helper code or runtime ownership.
4. A five-point centered difference of the complete callback value checks the
   returned derivative for one cation and one anion at non-source trial
   diameters. Step halving must exhibit a stable truncation/roundoff window;
   acceptance is
   `|g'_CppAD-g'_FD| <= max(1e-3 J mol^-1 angstrom^-1,
   20|g'_h-g'_(h/2)|, 2e-8|g'_CppAD|)`.
5. The public reference sequence passes its existing `5e-5` log-fugacity
   limit, the derivative passes the terminal relative test above, and the
   direct Born partial tends toward zero as positive diameter is increased.
   No large-diameter point is promoted as a source-domain physical claim.
6. Source and isolated installed-wheel tests prove the 88-byte old prefix,
   104-byte new table, 272-byte result, exact offsets, no-write undersized
   result behavior, capsule lifetime, and success/error status parity.
7. Installed-artifact negative space covers every unlisted fingerprint,
   reordered components, synthetic copies, methanol/ethanol, extra salts,
   nonfinite and nonpositive inputs, wrong state, H+, I-, sulfate,
   transition-metal/V species, missing source files, and absence of Provider
   implementation sources from the wheel.
8. A production-source ownership check finds the molality schedule and
   convergence-tolerance literals, electroneutral composition arithmetic, sequence loop,
   convergence comparison, and terminal-fugacity selection only in
   `electrolyte_reference.hpp`. `reference.py`, `_native/module.cpp`, and
   `_native/native_sdk.cpp` may contain delegation and result translation but
   none of that algorithm. The same check proves the public fixed-model path
   and active-Born callback both reach this owner and that no copied schedule
   or composition helper remains.

No frozen callback output generated by the production path qualifies as the
independent oracle. No new native target, density-root owner, third derivative,
solver, persistence layer, or compatibility layout is required.

## Canonical owners affected by a later implementation

- equation and active derivative: `src/epcsaft/_native/born.hpp`,
  `src/epcsaft/_native/derivatives.hpp`, and the existing canonical Helmholtz
  call path;
- reference sequence and composition:
  `epcsaft::detail::electrolyte_reference_sequence` in
  `src/epcsaft/_native/electrolyte_reference.hpp` as the sole future internal
  owner, with `src/epcsaft/reference.py` delegating through the existing private
  `_core` binding and retaining only its public behavior;
- density closure: `src/epcsaft/_native/density.hpp`, without copying its
  algorithm;
- model-bound transport and declaration:
  `src/epcsaft/_native/native_sdk.cpp` and
  `src/epcsaft/include/epcsaft/native_sdk_v1.h`;
- scientific contract: `docs/science/full-provider-eos.md`;
- SDK compatibility/status: `docs/adr/0001-native-consumer-transport.md`,
  `CONTEXT.md`, and `ARCHITECTURE.yaml`; and
- package evidence: the existing native SDK, Born, reference, and distribution
  test owners.

Regression owns the eventual target residual and Ceres integration. Validation
owns the immutable primary-source target packet and installed-artifact
cross-check. Migration and permanent lab alone can adjudicate authority.

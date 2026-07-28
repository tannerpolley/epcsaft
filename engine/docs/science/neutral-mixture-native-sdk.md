# Neutral Mixture Native SDK Candidates

Status: local candidate pending independent review and promotion receipt

Authority effect: none

## Scope

The append-only `evaluate_mixture_phase` entry in `epcsaft.native_sdk.v1`
evaluates one fixed-temperature, neutral, nonassociating phase for the
source-backed Gross--Sadowski methane/ethane parameter set. The capsule is
bound to an immutable model, parameter fingerprint, source temperature domain,
component count, and component order. One further append-only
`evaluate_mixture_phase_kij` pointer evaluates the same phase with the single
unlike interaction parameter active.

The admitted candidate fingerprints are:

- methane, ethane: `sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286`;
- ethane, methane: `sha256:cb0c63530f96d51d5a24f48a2d02ae6567d163abdc373099ce096f7ff0a896f5`.

Both fingerprints describe an ordered selection from the same installed
source bundle. They do not admit synthetic copies or fitted candidates.

The SDK table also appends one immutable capability descriptor for either
ordered selection. It identifies the callback as derivative-ready for one
unordered binary `k_ij` on a neutral, nonassociating binary model, while its
authority effect remains `none`. Pure, associating, ionic, or otherwise
unsupported model shapes advertise no such descriptor.

The descriptor binds a topology fingerprint separately from the parameter
fingerprint. Its canonical payload contains the ordered component IDs,
normalized component charges, resolved topology, site IDs and
multiplicities, association-pair count, and the fixed ideal, hard-chain, and
dispersion formulation identifiers. Parameter values and source records are
not part of this identity. Reordering the components therefore changes the
topology fingerprint; changing a value without changing topology does not.

## Thermodynamic function

For ordered positive component amounts `n_i [mol]`, positive phase volume
`V [m^3]`, and fixed source-domain temperature `T [K]`, define

```text
n       = sum_i n_i
rho_m   = n / V
x_i     = n_i / n
n_ref   = 1 mol
rho_ref = 1 mol m^-3

Phi(T,n_1,...,n_C,V) = A/(R T n_ref)
  = sum_i (n_i/n_ref) * [ln((n_i/V)/rho_ref) - 1]
    + (n/n_ref) * [a_hc^r(T,rho_m,x) + a_disp^r(T,rho_m,x)]
```

The ideal term contains ideal mixing through the component densities. The
residual term calls the canonical scalar-generic hard-chain and dispersion
owner; there is no copied EOS equation. One CppAD tape over numerical
coordinates `(n_1_mol,...,n_C_mol,V_m3)` supplies the value, gradient, and
symmetric row-major Hessian. The pressure diagnostic is

```text
P [Pa] = -R T n_ref * partial(Phi)/partial(V_m3).
```

The numerical amount coordinates are magnitudes in mol and the volume
coordinate is the magnitude in m3. The returned derivative numbers therefore
carry the corresponding implied inverse-coordinate units. The fixed entry
returns no third tensor and no parameter derivatives.

The active entry uses coordinates `(n_1_mol,n_2_mol,V_m3,k_ij)`. The trial
value replaces only the unlike pair and is assigned symmetrically to matrix
positions `(1,2)` and `(2,1)`; diagonal values remain those of the bound model.
The same defining function and CppAD owner supply its value, four-element
gradient, and 4-by-4 row-major Hessian. The source fingerprint continues to
identify the immutable catalog model, not the trial coordinate.

## Native buffer contract

The function table appends `component_count`, `mixture_result_size`, and the
fixed callback after the accepted pure prefix, then appends exactly the active-
`k_ij` callback pointer. Both callbacks use the same result structure. A caller
must:

- keep the capsule alive for every call;
- pass exactly `component_count` positive finite amounts;
- set `coordinate_count = component_count + 1`;
- own non-null gradient and Hessian buffers with capacities of at least
  `coordinate_count` and `coordinate_count * coordinate_count`;
- set the exact result `struct_size`; and
- check both the returned status and result status before consuming outputs.

A result `struct_size` mismatch is the sole return-code-authoritative failure:
the provider returns `EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1` without writing
the untrusted result layout. After exact-size acceptance, every failure also
sets the matching result status and error text.

For the active callback, `coordinate_count` is exactly four and the caller
provides gradient and Hessian capacities of at least 4 and 16. The provider
accepts any finite `k_ij`; regression owns fit bounds, transforms, scaling,
targets, and acceptance.

The provider preserves the caller's pointers and capacities, copies exactly
the advertised number of gradient and Hessian values, and propagates neither
exceptions nor ownership across the C boundary. The Hessian layout is row
major in the documented coordinate order.

The descriptor repeats that exact coordinate contract as four typed records:
amounts for component indices zero and one in `mol`, volume in `m3`, and the
unordered pair `(0,1)` `k_ij` in `dimensionless` units. It reports three state
coordinates, one active parameter, derivative order two, the source
temperature limits, row-major tensor layout, component mapping, Helmholtz
basis, and the Provider status codes used for unsupported models and domain
failures. It defines no data target, fit variable, scale, bound, solver, or
acceptance policy.

## Evidence and boundary

Package evidence checks the candidate at `T = 220 K`,
`n = (0.4, 0.6) mol`, and `V = 1.25e-4 m3`. It compares the value and pressure
with the public explicit-density provider path, contracts the gradient and
Hessian against centered directional differences, and verifies component
permutation. A compact malformed-contract test covers null pointers, count and
capacity mismatches, source-domain rejection, nonpositive amounts, result-size
mismatch, and pure/mixture scope separation. Accepted-prefix field offsets are
unchanged.

Active-`k_ij` evidence contracts all four gradient and Hessian coordinates
against centered differences at 220 K and 280 K. At the catalog value
`k_ij = 0`, value, pressure, state gradient, and state Hessian agree with the
fixed callback. Negative evidence covers malformed four-coordinate buffers,
nonfinite trial values and state inputs, pure models, association, and ions;
a finite negative trial confirms that the provider applies no fit bound.
The descriptor evidence additionally freezes the append-only table offsets,
the exact Gross methane/ethane component mapping and topology fingerprint,
and zero advertised capabilities for an unsupported pure-model shape.

These candidates do not own or admit association, ions, additional mixture
families, density closure, phase identity, phase discovery, third derivatives,
Ipopt, Ceres, equilibrium residuals, regression targets, fit bounds, scaling,
or fit acceptance. Implementation and a passing wheel do not transfer
authority. Separate independent review, a promotion receipt, and each
consumer's installed-wheel gate remain required.

# Neutral Parameter Native SDK Candidates

Status: historical internal-transport design evidence; superseded by the
Unified Engine direct cutover

Authority effect: none

The names and capsule layout below remain internal implementation details.
They are not a public Python API, installed header set, external ABI, or
separately versioned SDK.

## Scope

The append-only `evaluate_mixture_phase` entry in `epcsaft.native_sdk.v1`
evaluates one fixed-temperature, neutral, nonassociating phase for the
source-bound model supplied by the caller. The capsule is bound to an
immutable model, parameter fingerprint, topology fingerprint, source
temperature domain, component count, and component order. One further
append-only `evaluate_mixture_phase_kij` and `evaluate_mixture_phase_lij` pointers evaluate
the same phase with one distinct unlike interaction parameter active. A final
append-only `evaluate_pure_phase_parameter` pointer evaluates a supported pure
nonassociating model with exactly one of `m`, `sigma`, or `epsilon/k` active.
Associating models use one generic fixed-topology descriptor. Pure models use
the original variable-width phase callback; multicomponent models use the
append-only associating-mixture callback with the same descriptor-owned slots.

The first reference-evidence fingerprints are:

- methane, ethane: `sha256:3a840001adcb8b82f44e48307ad61e566f6a65d9b82d8312299a439dbce09195`;
- ethane, methane: `sha256:99ae0ffd9f61aee8c2e0fb718960534a2d0700c1adf5b72c1984517eded21f5b`.

Both fingerprints describe an ordered selection from the same source packet.
Other source-bound caller models carry their own parameter and topology
fingerprints; transport admission does not make their values predictive or
authoritative.

The SDK table advertises two immutable capability descriptors for each
supported ordered model. They identify the callbacks as derivative-ready for
one unordered binary `k_ij` and one unordered binary `l_ij` on a neutral,
nonassociating binary model, while their authority effect remains `none`.
Pure, associating, ionic, or otherwise unsupported model shapes advertise
neither binary-interaction descriptor. A neutral, nonassociating,
constant-diameter pure
model advertises three component-identity descriptors, one for each scalar
pure family; it advertises no binary descriptor. A neutral constant-diameter associating model
instead advertises one structural fixed-topology descriptor. That descriptor
contains three ordinary slots for every component plus one energy and one
volume slot for every canonical allowed site pair, and the caller may request
any nonempty ordered subset. The model fixes component order, site order,
multiplicity, and the complete allowed-pair graph. No scalar association
capability, species-specific callback, dynamic topology, or compatibility
adapter is introduced.

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

Each active entry uses coordinates
`(n_1_mol,n_2_mol,V_m3,pair_parameter)`. The trial value replaces only the
requested unlike `k_ij` or `l_ij` pair and is assigned symmetrically to matrix
positions `(1,2)` and `(2,1)`; diagonal values remain those of the bound model.
The same defining function and CppAD owner supply its value, four-element
gradient, and 4-by-4 row-major Hessian. The source fingerprint continues to
identify the immutable resolved model, not the trial coordinate.

For a scalar pure trial `theta`, the callback holds the other two pure
parameters at their immutable model values and projects the existing exact
five-coordinate pure result onto

```text
(n_mol, V_m3, theta)
```

It returns the same `Phi` and pressure with a three-element gradient and
3-by-3 row-major Hessian. The family selects `m [dimensionless]`,
`sigma [angstrom]`, or `epsilon/k [kelvin]`. This path supports independent
one-family fits only; the existing joint callback remains the owner when all
three parameters vary in one solve.

## Native buffer contract

The function table appends `component_count`, `mixture_result_size`, and the
fixed callback after the accepted pure prefix, then appends the active-`k_ij`
and active-`l_ij` callback pointers, capability descriptors, and the scalar
pure callback. The fixed mixture, active binary, and scalar pure callbacks use
the same result structure. Every caller must:

- keep the capsule alive for every call;
- own non-null gradient and Hessian buffers with capacities of at least
  `coordinate_count` and `coordinate_count * coordinate_count`;
- set the exact result `struct_size`; and
- check both the returned status and result status before consuming outputs.

The fixed mixture callback takes exactly `component_count` positive finite
amounts and uses `coordinate_count = component_count + 1`.

A result `struct_size` mismatch is the sole return-code-authoritative failure:
the provider returns `EPCSAFT_NATIVE_STATUS_ABI_MISMATCH_V1` without writing
the untrusted result layout. After exact-size acceptance, every failure also
sets the matching result status and error text.

For either active binary callback, the caller passes exactly two positive
finite amounts, `coordinate_count` is exactly four, and the gradient and
Hessian capacities are at least 4 and 16. The provider
accepts any finite requested `k_ij` or `l_ij`; regression owns fit bounds,
transforms, scaling, targets, and acceptance.

For the scalar pure callback, `coordinate_count` is exactly three and the
caller passes one positive finite amount and provides capacities of at least
3 and 9. Unknown family codes and
nonfinite trial values fail as invalid arguments. Zero or negative `m`,
`sigma`, or `epsilon/k` fails the Provider's positive EOS domain; regression
still owns the narrower scientific fit bounds.

For the fixed-topology association callbacks, the descriptor owns the immutable
component/site order, multiplicities, allowed-pair graph, slot identity, unit,
positive physical domain, parameter fingerprint, topology fingerprint, and
installed-artifact fingerprint. Pure requests must be a nonempty subset in
advertised slot order. Mixture requests may use any nonempty unique subset of
advertised slots, and the caller's request order defines the parameter
coordinates. A pure caller supplies exactly `2 + request_count` gradient
entries for `(n,V,requested slots...)`. A mixture caller supplies exactly
`component_count + 1 + request_count` entries for
`(n_1,...,n_C,V,requested slots...)` and repeats the exact ordered component
IDs at the call boundary. Its result binds component-order and source-domain
fingerprints in addition to the descriptor identities. The square row-major
Hessian uses the same order. Unknown, duplicated, nonfinite, out-of-domain, or
identity-mismatched requests fail closed; reordered pure requests also fail.

The provider preserves the caller's pointers and capacities, copies exactly
the advertised number of gradient and Hessian values, and propagates neither
exceptions nor ownership across the C boundary. The Hessian layout is row
major in the documented coordinate order.

Each binary descriptor repeats its exact coordinate contract as four typed records:
amounts for component indices zero and one in `mol`, volume in `m3`, and the
unordered pair `(0,1)` `k_ij` or `l_ij` in `dimensionless` units. Each pure
descriptor has three records: amount in `mol`, volume in `m3`, and the active
component parameter in its declared unit. Binary descriptors report three
state coordinates; pure descriptors report two. All report one active
parameter, derivative order two, the source
temperature limits, row-major tensor layout, component mapping, Helmholtz
basis, and the Provider status codes used for unsupported models and domain
failures. It defines no data target, fit variable, scale, bound, solver, or
acceptance policy.

## Evidence and boundary

Package evidence checks the candidate at `T = 220 K`,
`n = (0.4, 0.6) mol`, and `V = 1.25e-4 m3`. It compares the value and pressure
with an independent explicit-density calculation, compares the gradient and
Hessian against centered directional differences, and verifies component
permutation. Production validation rejects malformed buffers, invalid domains,
and nonpositive amounts, but NVDD does not freeze those rejection paths as a
compatibility test surface.

Active-pair evidence contracts all four gradient and Hessian coordinates
against centered differences at 220 K and 280 K. At the resolved baseline value
`k_ij = 0` or `l_ij = 0`, value, pressure, state gradient, and state Hessian
agree with the fixed callback. Negative evidence covers malformed four-coordinate buffers,
nonfinite trial values and state inputs, pure models, association, and ions;
a finite negative trial confirms that the provider applies no fit bound.
Descriptor metadata records the current coordinate meaning, reference mapping,
and provenance fingerprints. It is implementation documentation, not an
append-only ABI or table-layout test contract.

Scalar pure evidence checks each family against the corresponding slice of the
joint five-coordinate result, contracts its value gradient and Hessian against
centered directional differences, freezes descriptor units and component
identity, compiles the installed C header layout, and rejects unknown families
and nonfinite trials. It does not claim that three independent scalar fits are
equivalent to the accepted joint pure-parameter solve.

For any admitted neutral constant-diameter fixed topology, its phase callback
accepts a nonempty ordered subset of descriptor slots. Slots identify ordinary component
parameters or association energy/volume by canonical component/site endpoints.
The returned coordinate order is `(n,V,requested slots...)` for a pure model
or `(n_1,...,n_C,V,requested slots...)` for a mixture, with a complete exact
row-major Hessian on the ordinary `sigma_ij^3 kappa` basis. Compact pure 2B and
three-site cases check every derivative against an independent equation-level
oracle. Anonymous binary and ternary cross-associating cases finite-difference
every mixture gradient and Hessian column and verify component permutation.
No scalar association callback, species name, alternate basis adapter,
topology name, or compatibility route remains.

Each EOS slot remains a distinct derivative coordinate. A Regression sharing
map may assign one fitted value to several requested slots and accumulate the
returned columns; EOS accepts those equal-valued slots without owning the
consumer's sharing policy or collapsing endpoint identities.

These candidates do not own topology discovery, induced-association inference,
ions, multicomponent interaction tensors, density closure,
phase identity, phase discovery, third derivatives, Ipopt, Ceres, equilibrium
residuals, regression targets, fit bounds, scaling, or fit acceptance.
Implementation and a passing wheel do not transfer authority. Independent
review, a promotion receipt, and each consumer's installed-wheel gate remain
required.

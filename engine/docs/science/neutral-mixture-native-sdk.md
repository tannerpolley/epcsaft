# Neutral Mixture Native SDK Candidate

Status: local candidate pending independent review and promotion receipt

Authority effect: none

## Scope

The append-only `evaluate_mixture_phase` entry in `epcsaft.native_sdk.v1`
evaluates one fixed-temperature, neutral, nonassociating phase for the
source-backed Gross--Sadowski methane/ethane parameter set. The capsule is
bound to an immutable model, parameter fingerprint, source temperature domain,
component count, and component order.

The admitted candidate fingerprints are:

- methane, ethane: `sha256:307fcb28d535b94782f3e3caf4012c0c8c0dc87ee4239d6c316de56553543286`;
- ethane, methane: `sha256:cb0c63530f96d51d5a24f48a2d02ae6567d163abdc373099ce096f7ff0a896f5`.

Both fingerprints describe an ordered selection from the same installed
source bundle. They do not admit synthetic copies or fitted candidates.

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
carry the corresponding implied inverse-coordinate units. The entry returns
no third tensor and no parameter derivatives.

## Native buffer contract

The function table appends `component_count`, `mixture_result_size`, and the
callback after the accepted pure prefix. A caller must:

- keep the capsule alive for every call;
- pass exactly `component_count` positive finite amounts;
- set `coordinate_count = component_count + 1`;
- own non-null gradient and Hessian buffers with capacities of at least
  `coordinate_count` and `coordinate_count * coordinate_count`;
- set the exact result `struct_size`; and
- check both the returned status and result status before consuming outputs.

The provider preserves the caller's pointers and capacities, copies exactly
the advertised number of gradient and Hessian values, and propagates neither
exceptions nor ownership across the C boundary. The Hessian layout is row
major in the documented coordinate order.

## Evidence and boundary

Package evidence checks the candidate at `T = 220 K`,
`n = (0.4, 0.6) mol`, and `V = 1.25e-4 m3`. It compares the value and pressure
with the public explicit-density provider path, contracts the gradient and
Hessian against centered directional differences, and verifies component
permutation. A compact malformed-contract test covers null pointers, count and
capacity mismatches, source-domain rejection, nonpositive amounts, result-size
mismatch, and pure/mixture scope separation. Accepted-prefix field offsets are
unchanged.

This candidate does not own or admit association, ions, additional mixture
families, binary parameter fitting, density closure, phase identity, phase
discovery, third derivatives, parameter derivatives, Ipopt, Ceres,
equilibrium residuals, regression targets, or fit acceptance. Implementation
and a passing wheel do not transfer authority. A separate independently
reviewed promotion receipt and each consumer's installed-wheel gate remain
required.

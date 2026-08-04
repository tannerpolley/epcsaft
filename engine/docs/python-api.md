# Python API guide

The public parameter routes in `0.2.0.dev0` are:

| Entry point | Purpose |
| --- | --- |
| `Parameters.from_dictionary(mapping)` | Resolve compact source-neutral values already held in memory. |
| `Parameters.from_bundle(path, components=None)` | Load a strict external bundle from an explicit path. |
| `Parameters.create_bundle_template(...)` | Create a non-overwriting external-bundle authoring skeleton. |
| `Mixture.from_dictionary(mapping)` | Construct a mixture through the same strict dictionary loader. |

There is no installed catalog, record-construction facade, Data import, network
lookup, or sibling-checkout discovery. Both loaders copy caller containers into
immutable records and produce the same scientific fingerprint from normalized
values, component order, options, and validity. Bundle citation metadata does
not alter that equation identity.

Immutable Data migration packets retain a historical `purpose` value. The path
loader accepts that exact legacy value as inert external-file metadata only; it
does not expose catalog discovery or use purpose in scientific identity or SDK
admission.

## Dictionary schema

The top-level mapping is exact; every field is required and unknown fields fail
closed:

```python
{
    "schema": "epcsaft.parameters",
    "schema_version": 1,
    "components": ("component-a", "component-b"),
    "parameters": {...},
    "options": {...},
    "validity": {"kind": "unknown"},
}
```

`components` establishes all vector, matrix, composition, and result ordering.
Pure fields are parallel vectors of that length. Pair fields are symmetric
square matrices with an exact zero diagonal. Use `None` only in optional vector
positions where the selected topology does not require that component value.

Each scientific field accepts exactly one abbreviation or descriptive spelling.
Supplying both is `duplicate-alias`.

| Abbreviation | Descriptive spelling | Fixed unit |
| --- | --- | --- |
| `mw` | `molar_mass` | kg/mol |
| `m` | `segment_count` | dimensionless |
| `s` | `segment_diameter` | Å |
| `e` | `dispersion_energy_over_k` | K |
| `z` | `charge_number` | exact integer |
| `epsilon_r` | `relative_permittivity` | dimensionless |
| `d_born` | `born_diameter` | Å |
| `f_solv` | `solvation_factor` | dimensionless |
| `v_dielectric` | `dielectric_volume` | m³/mol |
| `t_dielectric` | `dielectric_temperature` | K |
| `a_ion` | `ion_specific_suppression_coefficient` | dimensionless |
| `mu2` | `dipole_scaling` | debye² |
| `alpha` | `polarizability_scaling` | Å³ |
| `a2` | `correlation_integral_parameter` | dimensionless |
| `mu` | `dipole_moment` | debye |
| `Q` | `quadrupole_moment` | debye Å |
| `k_ij` | `dispersion_interaction` | dimensionless matrix |
| `l_ij` | `diameter_interaction` | dimensionless matrix |

`parameters.sites` contains explicit component/site/class/multiplicity rows.
`parameters.association` contains explicit endpoint pairs with
`association_energy_over_k` in K and dimensionless `association_volume`.
`parameters.correlations` accepts `constant` correlations for segment diameter,
relative permittivity, or solvation factor, plus the established
constant-plus-exponentials segment-diameter form.

Options use the same one-to-one alias rule:

| Abbreviation | Descriptive spelling |
| --- | --- |
| `electrolyte` | `electrolyte_formulation` |
| `permittivity_model` | `relative_permittivity_formulation` |
| `epsilon_r_ion` | `ionic_region_relative_permittivity` |
| `a_dh` | `ion_fraction_suppression_coefficient` |
| `polar_model` | `polar_formulation` |

`polar_formulation="gross-vrabec-point-multipole"` admits optional physical
`dipole_moment` and `quadrupole_moment` values for neutral components. Missing
moments are exactly inactive. Negative, nonfinite, charged-component, reduced-
moment, and multiplicity inputs fail closed. The dielectric `dipole_scaling`
field is separate and cannot activate the polar Helmholtz contribution.

Polar mixture states require strictly positive mole fractions for every
model component. Composition-boundary states are rejected because the
family-wise Padé form does not provide a single CppAD derivative graph at an
incipient polar component; interior compositions remain fully differentiable.

Validity is either exactly `{"kind": "unknown"}` or a `fit-range` /
`reported-conditions` mapping with at least one explicit bound. Bounds are
`temperature_min_k`, `temperature_max_k`, `pressure_min_pa`,
`pressure_max_pa`, and `ion_mole_fraction_max`.

Wrong versions, unknown keys, nonfinite values, malformed dimensions,
asymmetric matrices, inconsistent topology, and extraneous formulation fields
raise `ParameterError` with a stable `code` and context.

## External bundles

`Parameters.from_bundle(path)` validates the manifest and table inventory at
that path. `components=` may select and order a complete subset; omitting it
uses the declared component order. The loader does not require the Data
repository. A caller may materialize a packet from ePC-SAFT-data, then pass its
`parameters/` directory directly.

`parameters.to_bundle(path)` writes the source bundle to a new path, and
`Parameters.create_bundle_template(path, components=..., profile=...)` creates
one of the `neutral`, `associating`, `electrolyte`, or
`associating-electrolyte` skeletons. Neither operation overwrites an existing
target.

## State and properties

Construct `Mixture(parameters)`, then provide exactly one of `rho` or `P`:

```python
state = mixture.state(T=298.15 * u.K, rho=1000 * u.mol / u.m**3, x=(0.4, 0.6))
state = mixture.state(T=298.15 * u.K, P=1 * u.bar, x=(0.4, 0.6), phase="vapor")
```

`phase` is valid only for pressure closure. `State.ares`, `ures`, `hres`,
`gres`, `sres`, `cvres`, and `cpres` return Pint molar quantities by default
and reduced scalars with `basis="reduced"`. `State.fugacity` preserves component
order. Admitted electrolyte models can create a reference with
`Mixture.reference(...)` and project activities with `reference.state(...)`.
`State.ares_polar` is the dimensionless aggregate
`a_DD + a_QQ + a_DQ35`; no family-specific public contribution is exposed.

## Native SDK

`native_sdk(mixture)` binds the resolved model to the append-only v1 capsule.
Availability is determined by structural capability: topology, active fields,
exact component order, and declared applicability. Parameter and topology
fingerprints remain deterministic handshake values, but paper, packet, and
bundle membership never admit or reject an endpoint.

For an admitted neutral constant-diameter associating model, the installed header exposes one
`fixed_topology_association` descriptor. The descriptor lists the immutable
component/site graph and typed active slots. Pure callers use
`evaluate_fixed_topology_association` and receive exact derivatives in
`(amount, volume, requested slots...)`; multicomponent callers use
`evaluate_fixed_topology_associating_mixture` and receive exact derivatives in
`(component amounts in model order, volume, requested slots...)`. Exact
component IDs are passed back to the mixture callback; its result binds
parameter, topology, component-order, domain, and installed-artifact
fingerprints. There are no species-specific callbacks, topology inference,
2B-specific entry points, or alternate association-basis selector.

Native phase blocks evaluate the same fingerprint-distinct aggregate for polar
models; no polar callback or family-specific result type is added.

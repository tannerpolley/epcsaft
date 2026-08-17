# Python API

`epcsaft` owns one strict parameter document, one EOS, and one installed wheel.
It does not discover catalogs, sibling checkouts, callbacks, or alternate
providers.

## Parameters

```python
from epcsaft import Parameters

parameters = Parameters.from_json("parameters.json")
same_model = Parameters.from_mapping(document)
subset = parameters.select(("water", "sodium-cation", "chloride-anion"))
```

Both loaders validate the same versioned `epcsaft.parameters` document. The
document separates fixed component/model specifications from continuous
coefficients and records exact canonical identities for selectable
coefficients. Component order is explicit and survives selection.

Continuous identity examples are:

```text
component/propane/segment_count
component/water/segment_diameter/constant-plus-sum-of-exponentials/constant
pair/ethane/methane/k_ij
association/ethanol/a/ethanol/b/energy_over_k
model/relative_permittivity/ionic_region_relative_permittivity
```

Molar mass, charge, component identity, molecular moments, association sites,
multiplicity, topology, model-family choices, citations, and provenance are
fixed specifications. Unknown fields, incomplete provenance, ambiguous edge
sources, invalid units, nonfinite values, and inconsistent topology fail
closed with `ParameterError`.

## EOS state and saturation

```python
from epcsaft import Mixture, Pure, unit_registry as u

pure = Pure(parameters.select(("propane",)))
state = pure.state(T=300 * u.kelvin, rho=4_000 * u.mole / u.meter**3)
saturation = pure.saturation(T=300 * u.kelvin)

mixture = Mixture(parameters)
state = mixture.state(
    T=298.15 * u.kelvin,
    P=1 * u.bar,
    x=(0.98, 0.01, 0.01),
    phase="liquid",
)
```

`Pure` requires one component. `Mixture` requires at least two and always
receives an explicit composition. Both return `State`. A state accepts exactly
one of molar density or pressure; `phase` selects a pressure-closure branch.

`State.ares`, `ures`, `hres`, `gres`, `sres`, `cvres`, and `cpres` return Pint
molar quantities by default and reduced scalars with `basis="reduced"`.
`State.fugacity` preserves component order. `Pure.saturation` returns the
common pressure, vapor and liquid densities, convergence diagnostics, and
physical acceptance evidence. EOS owns this solve and its packing domain.

## Exact selected-parameter derivatives

```python
from epcsaft import ActiveParameterSet

pure_parameters = parameters.select(("propane",))
pure = Pure(pure_parameters)
active = ActiveParameterSet(
    pure_parameters,
    (
        "component/propane/segment_count",
        "component/propane/segment_diameter",
        "component/propane/dispersion_energy_over_k",
    ),
)

pair_parameters = parameters.select(("dipolar", "quadrupolar"))
pair_active = ActiveParameterSet(
    pair_parameters,
    (
        "pair/dipolar/quadrupolar/k_ij",
        "pair/dipolar/quadrupolar/l_ij",
    ),
    values=(0.12, 0.01),
)

fixed_state = pure.observable_block(
    active,
    identity="eos.fixed-state",
    T=300 * u.kelvin,
    rho=4_000 * u.mole / u.meter**3,
)
saturation = pure.observable_block(
    active,
    identity="eos.pure-saturation",
    T=300 * u.kelvin,
)
```

`ActiveParameterSet` freezes an ordered, nonempty selection of exact continuous
identities from one `Parameters` value. An explicit `values` sequence is
interpreted in descriptor units and does not mutate the base `Parameters`.
Each EOS-owned block retains its versioned descriptor, exact active identity
and order, row units and bases, state metadata, values, Jacobian layout, and
structured row/shared-solve status. Saturation blocks additionally carry the
two density-branch brackets, root residuals, and residual-density derivatives.
Any number of constant unordered `k_ij` and `l_ij` identities may be selected
for a fixed-topology mixture; each coordinate records its canonical component
endpoints, resolved reciprocal matrix application, source validity domain, and
intrinsic domain (`k_ij <= 1`, `l_ij < 1`).

A cross-association energy correction has the canonical identity
`pair/<component-a>/<component-b>/association/<role-a>/<role-b>/k_hb_ij`,
with both pairs unordered. It applies `epsilon_AB = epsilon_AB,0 * (1 -
k_hb_ij)` to every combining-rule edge with the declared component and role
pairs without changing association volume. The shorthand
`pair/<component-a>/<component-b>/association/k_hb_ij` resolves only when the
parameter topology contains exactly one eligible declared role pair; missing,
ambiguous, directed, and unmatched selections fail closed. Its intrinsic
domain is `k_hb_ij < 1`.

Every derivative result contains:

- values and row-major selected-parameter derivatives;
- ordered active identities, trial values and units, plus observation row,
  value-unit, and shape metadata;
- exact parameter and active-set fingerprints; and
- the active EOS topology fingerprint.

`Pure` and `Mixture` provide fixed-state Helmholtz derivatives. Admitted base
coefficients also provide fixed-pressure fugacity derivatives, including the
implicit density closure. `Pure` provides exact implicit saturation
sensitivities for pressure and both coexisting densities.

Association sensitivities include the implicit site-fraction solve for pure,
mixture, neutral, electrolyte, constant-diameter, and temperature-dependent
diameter contexts. Auxiliary electrolyte methods expose EOS derivatives for
Born, solvation-factor, and relative-permittivity quantities, but their
presence does not make those quantities admitted Regression observations.
Unsupported identity, observation, and model combinations fail before an
optimizer is built.

CppAD is the production derivative authority. Finite differences exist only as
independent numerical evidence in Engine tests; they are never runtime
fallbacks.

## Regression

```python
from epcsaft.regression import FitResult, ParameterFit, fit

specification = ParameterFit.from_json("parameter-fit.json")
result: FitResult = fit(specification, base_path="campaign")
candidate_parameters = result.to_parameters()
```

`ParameterFit` declares one parameter document, model kind, active selectors,
bounds, scales, ties, observations, source hashes, residuals, multistart policy,
and optional optimizer controls. Defaults are 25 maximum iterations and `1e-6`
function, gradient, and parameter stopping tolerances. `fit` resolves each
observation against an owner descriptor, sends
the owner's values and exact ordered Jacobian blocks through the private native
Ceres core, and returns every start without mutating the input parameters.
Each start retains predictions, scaled residuals and Jacobians, rank and
conditioning, active bounds, and typed owner-evaluation failures.
`FitResult.to_parameters()` creates a new immutable candidate parameter
document with fit provenance. It does not accept, promote, or adopt the values.

Selectable coordinates include continuous component and association parameters
and supported pair or model coefficients when the chosen owner output has the
required exact derivative. Generic configured-equilibrium outputs use the same
fit path; Regression does not select a phase topology or create a case-specific
public route.

## Native boundary

The compiled EOS and solver code is private to the wheel. Python consumers use
`Parameters`, `Pure`, `Mixture`, `State`, and the public Equilibrium and
Regression modules. Raw native slots, capability callbacks, C headers, and
provider capsules are not public interfaces.

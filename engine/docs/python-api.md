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

fixed_state = pure.phase_parameter_derivatives(
    active,
    T=300 * u.kelvin,
    rho=4_000 * u.mole / u.meter**3,
)
saturation = pure.saturation_parameter_derivatives(
    active,
    T=300 * u.kelvin,
)
```

`ActiveParameterSet` freezes an ordered, nonempty selection of exact continuous
identities from one `Parameters` value. Every derivative result contains:

- values and row-major selected-parameter derivatives;
- observation, value-label, unit, and shape metadata;
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
fitted_parameters = result.to_parameters()
```

`ParameterFit` declares one parameter document, model kind, active selectors,
bounds, scales, ties, observations, source hashes, residuals, and multistart
policy. `fit` resolves that declaration, executes the private native Ceres
owner against exact EOS derivatives, and returns every start without mutating
the input parameters. `FitResult.to_parameters()` creates a new immutable
parameter document with fit provenance.

This cutover admits `pure_saturation` and `fixed_state_pressure` observation
families. Selectable coordinates include continuous component and association
parameters and constant neutral `k_ij` or `l_ij` pairs when the chosen
observation has the required exact EOS derivative. Mixture phase-equilibrium,
charged/electrolyte pair fitting, direct permittivity or Born calibration, and
reaction observations are not admitted. Adding one requires a source-neutral
Observation Family and scientific evidence; it never creates another public
fit route.

## Native boundary

The compiled EOS and solver code is private to the wheel. Python consumers use
`Parameters`, `Pure`, `Mixture`, `State`, and the public Equilibrium and
Regression modules. Raw native slots, capability callbacks, C headers, and
provider capsules are not public interfaces.

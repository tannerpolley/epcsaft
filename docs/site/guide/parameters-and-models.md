# Parameters and models

A parameter set declares the physical model as well as its numbers. Select
species, parameter source and version, model families, association sites and
unlike pairs before constructing a model. The Engine does not select physics
from an application name. The
[water notebook](../examples/parameters-and-water-state.ipynb) demonstrates a
published source selection; the [parameter algorithm](../theory/algorithms.md#alg-parameters)
and [EOS conventions](../theory/eos.md) define the units and model terms.

## Select source records

`Parameters.from_json` reads a parameter document with an optional component
selection. `Parameters.from_packets` combines named component selections from
versioned packets below a supplied directory. `parameter_document` exposes the
combined document before native construction. These calls are documented in
[parameter inputs](../api/parameters.rst).

The packet root in the quickstart is `data/packets`; it is a checkout path,
not a resource bundled with the wheel. Keep the selected packet versions,
component order, source locators and caller overrides with the calculation.
A combined document retains materialization source information. A successful
load establishes mathematical admission, not the validity of extrapolation.

| Quantity | Input unit |
|---|---|
| Segment count, unlike `k_ij` and `l_ij` | dimensionless |
| Segment and ionic diameters | Å |
| Dispersion and association energy divided by Boltzmann's constant | K |
| Association volume | dimensionless |
| Molar mass | kg/mol |
| State pressure and molar density | Pa and mol/m³ |

Source conversions must be explicit. The parameter reader and native builder
check the declared units; do not silently treat an Å coefficient as metres.

## Combine model choices deliberately

Packet combination refuses conflicting base, electrolyte, permittivity or
polar choices. A component without association sites can join a general-site
association model, but arbitrary model conflicts are not resolved by taking
the first source. Selecting two copies of the same component is also refused.

Sites and their interaction edges are explicit. A missing cross-association
edge is not permission to invent one. Declare the intended edge or the
source's deliberate omission. Unlike-component `k_ij` and `l_ij` overrides
refer to two selected, distinct components. A constant override replaces that
pair's corresponding temperature law as well as its constant coefficient.
Keep a caller override labelled as an assumption or a fitted value.

See the [ammonia parameter combination](../examples/index.md#ammonia-parameter-combination)
for an explicit combined declaration.

## Construct and reuse the model

Use `Pure(parameters)` for exactly one component and `Mixture(parameters)` for
two or more. `Pure` fixes the composition at one and offers a saturation
problem declaration. Both copy their complete physical inputs into the native
model. Editing a Python input afterwards does not update a prepared model.
Construct a new model for a different parameter set.

Pass thermal records in component order when the requested outputs consume
them. Declare active parameters when requesting parameter actions; selecting
an active coordinate changes what is differentiated, not the physical source
record. See [states](states-and-properties.md),
[parameter estimation](parameter-estimation.md), and the
[model reference](../api/eos.rst).

## Choose electrolyte and polar terms by their equations

The EOS includes the declared hard-chain, dispersion and association terms,
and optional electrostatic and polar contributions. A dipole-dipole residual
term and a bulk-permittivity choice are distinct declarations. The
[EOS theory](../theory/eos.md) defines the admitted families and exclusions.
See [aqueous NaCl inputs](../examples/index.md#aqueous-nacl-inputs)
for the Held source selection.

Charged ordinary states must be electroneutral. Mathematical admission does
not establish parameter accuracy at a temperature, concentration or solvent
composition. Read the source's fit domain and retain unknown applicability
as unknown; the [limits](../limits.md) explain this distinction.

# Examples

Every notebook on this site was executed when the site was built. To run one
yourself, open it from `docs/site/examples/` in the unpacked source archive:
each notebook locates the archive root from its own directory and reads
parameter packets and observations under `data/` and `analyses/`, which the
conda package does not contain. See the
[installation guide](../installation.md) for setup and the
[optional process tools](../contributing.md#optional-process-example-tools)
for the IDAES settler.

Each notebook identifies its sources, conditions, units, numerical comparisons
and claim limits. Numerical verification checks equations, balances or retained
reference calculations. The halite example includes a physical comparison;
the propane fit is calibration against the observations used in that fit.

```{toctree}
:maxdepth: 1

parameters-and-water-state
methane-caloric-properties
methane-ethane-phase-equilibrium
aqueous-nacl-activity
halite-solubility
ammonia-synthesis-equilibrium
water-isobutanol-discovery
water-isobutanol-solved-state-actions
propane-parameter-estimation
methane-asl-properties
idaes-two-liquid-settler
```

## Download the original notebooks

- {download}`Parameters and a water state <parameters-and-water-state.ipynb>`
- {download}`Methane caloric properties <methane-caloric-properties.ipynb>`
- {download}`Methane–ethane phase equilibrium <methane-ethane-phase-equilibrium.ipynb>`
- {download}`Aqueous NaCl activity <aqueous-nacl-activity.ipynb>`
- {download}`Halite solubility <halite-solubility.ipynb>`
- {download}`Ammonia synthesis equilibrium <ammonia-synthesis-equilibrium.ipynb>`
- {download}`Water–isobutanol phase discovery <water-isobutanol-discovery.ipynb>`
- {download}`Water–isobutanol solved-state actions <water-isobutanol-solved-state-actions.ipynb>`
- {download}`Propane parameter estimation <propane-parameter-estimation.ipynb>`
- {download}`Methane ASL properties <methane-asl-properties.ipynb>`
- {download}`IDAES two-liquid settler <idaes-two-liquid-settler.ipynb>`

## Notes on individual examples

Case-specific declarations and interpretation for the notebooks above. General
methods are in the [user guide](../guide/index.md) and the
[theory pages](../theory/index.md).

## Ammonia parameter combination

The [ammonia notebook](ammonia-synthesis-equilibrium.ipynb) shows
why compatible physical contributions may still require an explicit combined
declaration when source packets declare different polar choices. It preserves
the published pure records and explicitly assumes zero unlike interactions;
it is not a published ternary parameter fit.

## Aqueous NaCl inputs

The [aqueous NaCl notebook](aqueous-nacl-activity.ipynb) keeps the
Held source's electrolyte choices with its water/ion parameters.

## Formula and explicit-ion feeds

A salt formula feed and an explicit-ion feed have different total species
amounts. One mol NaCl supplies one mol Na⁺ and one mol Cl⁻. Keep formula-unit
molality distinct from the sum of ionic molalities; the
[NaCl notebook](aqueous-nacl-activity.ipynb) demonstrates the basis.

## Halite complementarity

The [halite notebook](halite-solubility.ipynb) checks conservation,
saturation and precipitate amount for one salt at 298.15 K; it does not
represent a solid solution.

## Ammonia fixed-temperature standard

The ammonia example is isothermal at 573.15 K. It checks the ideal-gas
thermochemical scalar and converts the pressure standard to the native density
standard explicitly. It then solves with the scalar density-basis constant.

## Ammonia and halite interpretation

For ammonia, check atom balances, positive amounts, pressure closure and the
independently reconstructed mass-action residual. Agreement with the
low-pressure analytic ideal limit verifies the limiting mathematics. The
dense-gas result uses assembled pure records, assumed zero unlike interactions
and stated extrapolations; it is a model prediction with unquantified physical
uncertainty. It is not a validated industrial reactor calculation.

For halite, check both dissolved formula amount and precipitated solid amount.
An absent solid must satisfy its complementary saturation condition; a
positive solid amount must satisfy saturation. The one-temperature source
comparison does not establish general salt accuracy, homogeneous reaction
accuracy or solid-solution behavior.

## Propane calibration

The propane notebook starts from a published parameter packet, fits pressure
and mass density at three temperatures, checks a reduced objective, then
reconstructs fitted parameters and re-evaluates predictions. These checks
verify the estimate's numerical consistency. They do not adopt a new canonical
parameter set or establish predictive accuracy beyond those fitted rows.

Unit log-residual
weights in the propane example are a calibration choice, not a claim about
measurement uncertainty.

## Water/isobutanol discovery

The water/isobutanol notebook uses the neutral binary subset at 293.15 K and
100,000 Pa. It compares the returned phases without relying on their output
order and checks feed conservation, then repeats the discovery with
`Discover(mode="fast")` and requires the same split. This is numerical
verification for that case, not a measured tie-line validation or a universal
search guarantee.

## Methane property transport

The [methane ASL notebook](methane-asl-properties.ipynb) loads
`libepcsaft_asl.so` beside the installed native extension and prepares
`property_asl_record`. It evaluates molar density (mol/m³), total molar
enthalpy (J/mol) and component ln φ (dimensionless) at 298.15 K and 100,000 Pa.
The callback uses the declared branch and the pressure-property coordinate
chart, including T, P and its supplied composition weights.

Pyomo `ExternalFunction` transports the values, gradients and Hessians. The
notebook compares them with public native property actions, retaining units
and coordinate order. Finding the installed library avoids accidentally
loading a different checkout's build. These are derivative-transport checks
for the neutral methane case; they do not qualify an ionic or reacting
flowsheet. Calling a property callback alone is not a solved process model.

## Two-liquid settler

The [two-liquid settler notebook](idaes-two-liquid-settler.ipynb)
uses a separate `equilibrium.asl_record` and the phase-component-amount
callback. It declares an affine feed map, phase supports, an anchor and an
identity gap to distinguish the water-rich and organic phases. The native
callback performs the declared equilibrium calculation and supplies exact
first and second actions for the qualified feed controls.

The model uses water/ethanol/isobutanol/Na⁺/Cl⁻ at 293.15 K and 100,000 Pa.
Formula-feed flows and explicit-ion flows have different totals: a NaCl
formula flow supplies the same flow to each ion. A one-mol/s formula feed
therefore contains more than one mol/s of explicit species. Preserve this
mapping in every inlet and outlet balance.

IDAES owns the flowsheet, ports, outlet equations, recovery constraint and
feed-scale optimization. The notebook uses one start with fixed T/P and
checks component/charge closure, derivative transport and independent local
optimality quantities. The Engine callback's phase amounts are mapped onto
the declared mol/s feed basis by the process equations.

## Short IDAES VLE flowsheet

The {download}`standard IDAES Flash script <idaes-vle-flash.py>` separates
process declaration from the longer numerical verification notebooks.
Run it from the checkout root after the
[process-model setup](../guide/process-models.md#library-and-external-solver-setup).
It prints the neutral binary phase split and native comparisons; numerical
assertions remain in the separate qualification script.

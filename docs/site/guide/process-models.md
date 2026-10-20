# Process models with IDAES and ASL

The Engine can supply native properties or a qualified solved-equilibrium
mapping to a Pyomo/IDAES model through an AMPL Solver Library (ASL) callback.
The process model owns its material/energy links, unit equations, objective
and outer solve. The Engine owns the thermodynamics and declared native
implicit calculations. See the [action algorithm](../theory/algorithms.md#alg-actions),
[equilibrium formulation](../theory/equilibrium.md) and
[model reference](../api/eos.rst).

## Pressure-property callbacks

The [methane ASL notebook](../examples/methane-asl-properties.ipynb) loads
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

## Equilibrium callbacks

The [two-liquid settler notebook](../examples/idaes-two-liquid-settler.ipynb)
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

The [equilibrium reference](../api/equilibrium.rst) and
[native equilibrium callback fields](../api/native/equilibrium.rst) define
the record and its controls. Keep the qualified phase-identity conditions;
coalescence or a changed support can make a derivative unavailable rather
than returning a meaningful zero.

## Preserve the calculation ownership

Python may declare the external process equations in Pyomo/IDAES and invoke
the outer process solver. This is the explicit process-model exception to
the native calculation rule. It does not authorize Python implementations of
Engine density, equilibrium, fitting or continuation calculations.

A property-based external phase-equilibrium model owns its own process rows
and their claims. Do not describe its result as an Engine equilibrium solve
unless it actually consumes the Engine's declared-equilibrium mapping.

## Optional tools and numerical limits

The property callback example needs the installed Engine library and the
Python packages; the settler also needs an external Ipopt executable and
PyNumero support. The repository locks the Python packages used by the full
docs build. Obtain IDAES extension tools directly from upstream and preserve
their notices, as described in
[optional process-example tools](../contributing.md#optional-process-example-tools).
Both the Engine and outer Ipopt select MUMPS for these examples.

The upstream executable contains HSL even when MUMPS is selected. It is an
external build-time tool: do not put it, its libraries, archives or environment
in Engine wheels, documentation downloads or uploaded site outputs.

The methane pressure-property comparison establishes neutral derivative
transport; the settler establishes the stated fixed-T/P ionic two-liquid
feed mapping. Neither implies arbitrary ionic/reacting thermal process
qualification, global optimality, or physical validation. Preserve the
case-specific supports, references, controls and checks before extending a
process calculation.

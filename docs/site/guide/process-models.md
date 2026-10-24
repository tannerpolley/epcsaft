# Process models with IDAES and ASL

The Engine can supply native properties or a qualified solved-equilibrium
mapping to a Pyomo/IDAES model through an AMPL Solver Library (ASL) callback.
The process model owns its material/energy links, unit equations, objective
and outer solve. The Engine owns the thermodynamics and declared native
implicit calculations. See the [action algorithm](../theory/algorithms.md#alg-actions),
[equilibrium formulation](../theory/equilibrium.md) and
[model reference](../api/eos.rst).

## Pressure-property callbacks

Load `libepcsaft_asl.so` beside the installed native extension and prepare
`property_asl_record`. The callback uses the declared branch and the
pressure-property coordinate chart, including T, P and supplied composition
weights. Pyomo `ExternalFunction` transports values, gradients and Hessians.
Finding the installed library avoids loading a different checkout's build.
Calling a property callback alone is not a solved process model. See
[methane property transport](../examples/index.md#methane-property-transport)
for the notebook's conditions, derivative checks and limits.

## Equilibrium callbacks

Prepare `equilibrium.asl_record` for the phase-component-amount callback.
Declare an affine feed map, phase supports, an anchor and an identity gap.
The native callback performs the declared equilibrium calculation and supplies
exact first and second actions for the qualified feed controls. IDAES owns
the flowsheet, ports, outlet equations, recovery constraint and optimization.
Process equations map the Engine's phase amounts onto the declared mol/s feed
basis. See the [two-liquid settler](../examples/index.md#two-liquid-settler)
for the formula-flow mapping, numerical checks and qualification limits.

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

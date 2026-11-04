# Process models with IDAES and ASL

The Engine can supply native properties or a qualified solved-equilibrium
mapping to a Pyomo/IDAES model through an AMPL Solver Library (ASL) callback.
The process model owns its material/energy links, unit equations, objective
and outer solve. The Engine owns the thermodynamics and declared native
implicit calculations. See the [action algorithm](../theory/algorithms.md#alg-actions),
[equilibrium formulation](../theory/equilibrium.md) and
[model reference](../api/eos.rst).

## Generic IDAES property package and VLE

Install the optional Python tools with the selected Engine wheel:

```bash
python -m pip install '/absolute/path/epcsaft-<version>-<tags>.whl[idaes]'
```

In a source checkout, `bash tools/setup.sh` creates that checkout's locked
`.venv` and editable Engine install. Run examples from the checkout root;
parameter packets under `data/` are not included in the wheel. The optional
`epcsaft.idaes` module uses IDAES 2.12.0 and Pyomo 6.10.1. It constructs a
`GenericParameterBlock` configuration with the existing FTPx state definition,
LogBubbleDew fugacity equalities and SmoothVLE phase-equilibrium formulation.
The Engine supplies density, total molar enthalpy and component ln(phi),
including exact first/second derivatives in T, P and composition weights.

Run the short standard-unit example:

```bash
.venv/bin/python docs/site/examples/idaes-vle-flash.py
```

It solves a neutral methane/ethane TP flash at 180 K and 1 MPa with a
1 mol/s equimolar feed, using a standard IDAES `Flash`. It prints phase
flows, compositions, densities, native comparisons and balance residuals.
The separate qualification script retains numerical assertions:
`analyses/2026-idaes-asl-adapter/scripts/qualify_properties.py`.

State variables are molar flow (mol/s), T (K), P (Pa) and component mole
fractions. Choose `phases=("liquid",)`, `("vapor",)` or both in
`property_configuration`; supply state bounds, initial values and reference
scales explicitly. Density has units mol/m3 and total molar enthalpy J/mol.
Thermochemistry is optional for density, fugacity and material-only VLE;
requesting enthalpy without it raises `MissingThermochemistry`. Caloric
quantities retain the supplied thermal record's anchors and standard pressure.
Enthalpy flows are not comparable with another IDAES package's
`enth_mol_phase`; mixing or translating between packages is unqualified.

The generic initializer requires additional pure saturation correlations and
is unavailable here. Supply state starts, as in the example; native Engine
bubble/dew solves may supply boundary guesses without a Python root loop.
The example uses distinct phase starts rather than starting from the final
native TP solution. Native calculations remain the Engine's responsibility;
IDAES owns the process balances, fugacity equalities and outer solver.

### Library and external solver setup

`epcsaft.idaes.asl_library()` resolves the adapter beside the loaded `_core`
extension: `libepcsaft_asl.so` on Linux, `libepcsaft_asl.dylib` on macOS or
`epcsaft_asl.dll` on Windows. It raises if the installed library is missing.
Set `AMPLFUNC` to its absolute path for the outer ASL-enabled solver;
multiple existing entries are separated by newlines (Pyomo's `amplfunc_merge`
uses this convention). In-process Python evaluation alone does not configure
an external executable.

Obtain IDAES extension tools directly from upstream using `idaes get-extensions`
and its supported platform selection, outside the checkout. The normal Ipopt
path is `Path(idaes.bin_directory) / "ipopt"`; the example uses the upstream
IDAES Ipopt 3.13.2 with MUMPS. The wheel and Python extra do not install this
executable. See [optional tools](../contributing.md#optional-process-example-tools)
for external installation and notices. Keep `EPCSAFT_HSL_LIBRARY` and
`EPCSAFT_IPOPT_OPTIONS_FILE` unset for this example. Upstream binaries may
contain Coin-HSL even with MUMPS selected; never redistribute them or copy them
into Engine wheels, documentation outputs or shared caches.

Put Ipopt options in an option file and pass `option_file_name`; this executable's
AMPL command-line option interface does not recognize every native Ipopt option.
Select `hessian_approximation exact`, `linear_solver mumps`, `tol 1e-11`,
`constr_viol_tol 1e-11` and `acceptable_iter 0`. Read the runtime user-options
log with `print_user_options yes` and `print_level 5`, rather than `ipopt -=`.
The derivative verification additionally selects `derivative_test second-order`,
`derivative_test_tol 1e-4`, `derivative_test_perturbation 1e-6` and
`point_perturbation_radius 0` to check at the declared in-domain start.

Use `solve_with_diagnostics(solver, flowsheet, logfile=path)` to preserve the
combined solver output and receive an exception with native callback messages
and the log path if the solve fails. A trial refusal recovered by Ipopt remains
in the log and is not a failed final solve. For direct Python evaluation,
`PropertyExternalFunction` reads ASL `Errmsg` and raises the native cause,
instead of silently returning NaNs. `check_state(state, phase)` reads native
endpoint diagnostics and warns if a sole root has no classified phase side.

### Phase identity and limits

The initial VLE domain is neutral, inactive-parameter, positive-composition,
steady molar-flow liquid/vapor equilibrium with classified local density roots.
Each phase's own mixture density root must classify as its requested liquid or
vapor side on its own (T, P, x) isotherm, or the strict equilibrium call is
refused. This is mixture-root classification, not a pure-component critical
test. A component above its own critical temperature is admitted when both
phase roots classify, including methane at the retained 219.459 K near-dew
point. Unclassified VLE roots and resolved opposite branches are refused.
Single-phase packages admit a sole unclassified root with an explicit endpoint
warning. No arbitrary density cutoff selects a phase.

SmoothVLE evaluates fugacity equalities at a smoothed temperature between the
mixture bubble and dew temperatures. Incipient-phase flow remains small rather
than exactly zero at phase disappearance. The qualification uses the finest
smoothing pair (1e-6 K, 5e-8 K) for its absent-fraction criterion; defaults
(0.01 K, 0.0005 K) and the intermediate pair expose smoothing bias. A small
phase flow is not evidence of exact coexistence. Critical passage, structural
zero-component supports, ionic/reacting VLE and global stability are unqualified.

Inactive-phase calorics at actual stream T can be unavailable even when its
incipient fugacity calculation exists at the smoothed temperature. Material-only
flash does not request them. Energy balances require valid total enthalpy in
every phase they consume; unavailable enthalpy is never replaced by zero or
residual enthalpy. Strict and ordinary records retain two native models/caches;
coincident points can incur an extra density search. The longer
[methane property notebook](../examples/index.md#methane-property-transport)
retains the local property-transport verification.

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

## What runs in Pyomo and what runs in the Engine

Python may declare the external process equations in Pyomo/IDAES and invoke
the outer process solver. This is the only place where equations outside the
Engine are solved; density, equilibrium, fitting and continuation
calculations still run in the Engine.

A property-based external phase-equilibrium model is responsible for its own
process equations and their claims. Do not describe its result as an Engine equilibrium solve
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

# Quickstart

After [installation](installation.md), the first calculation below runs from any
directory with the conda package. Later examples use the source archive's packets.

## First saturation and flash

The conda package ships `feos-gross-sadowski-2001/1` and
`feos-gross-sadowski-2002/1`. This example uses the first packet's pure propane,
methane and ethane records. The second packet supplies associating fluids for
other calculations. The illustrative binary interaction assumption is explicitly
$k_{ij}=0$; it is not a fitted methane/ethane interaction.

```python
from pathlib import Path
import epcsaft
from epcsaft import Parameters, Pure, Mixture, equilibrium as eq

packets = Path(epcsaft.__file__).parent / "data" / "packets"
packet = "feos-gross-sadowski-2001/1"
propane = Pure(Parameters.from_packets(packets, [(packet, ["propane"])]))
saturation = eq.solve_equilibrium(propane, propane.saturation(300.0))  # K
mixture = Mixture(Parameters.from_packets(
    packets, [(packet, ["methane", "ethane"])],
    k_ij={("methane", "ethane"): 0.0}))
problem = eq.Problem(T=250.0, P=3_750_000.0,  # K, Pa
    feed=eq.Amounts({"methane": 0.3, "ethane": 0.7}),  # mol
    phases=[eq.Phase("L", kind="liquid"), eq.Phase("V", kind="vapor")])
flash = eq.solve_equilibrium(mixture, problem)
assert saturation.success and flash.success
assert abs(saturation.pressure / 998660.89574 - 1) < 1e-8
assert abs(flash.amounts[1] - 0.1646593108) < 1e-8
assert max(map(abs, saturation.classified_residuals)) <= 1e-10
assert max(map(abs, flash.classified_residuals)) <= 1e-10
print(saturation.pressure, flash.amounts)  # Pa, mol in L/V order
```

These are model predictions and numerical checks, with unquantified physical
uncertainty. A declared liquid/vapor split performs no global stability search.
On source installations without bundled packets, use the source archive inputs
below. The equivalent shorter call for the flash is
`eq.tp_flash(mixture, T=250.0, P=3_750_000.0, feed=problem.feed)`.

## Parameters, model and state

```python
from epcsaft import Parameters, Pure

parameters = Parameters.from_packets(
    "data/packets", [("held-cameretti-sadowski-2008/4", ["water"])])
water = Pure(parameters)
state = water.state(298.15, P=100_000.0, phase="liquid")
assert abs(state.pressure / 100_000.0 - 1.0) < 1e-10
mass_density = state.molar_density * parameters.molar_masses[0]
print(mass_density)  # kg/m³; mol/m³ times kg/mol
```

Temperature is in K, pressure in Pa and molar density in mol/m³. Specify
exactly one of pressure `P` and molar density `rho`. A liquid branch declaration
selects a density branch; it does not establish global phase stability.
The [state guide](guide/states-and-properties.md) explains Ideal, Residual
and Total quantities. This model has no ideal thermal record, so its total
enthalpy is unavailable; residual enthalpy remains a separate quantity.

```python
assert state.total_enthalpy is None
print(state.thermochemistry_status)
```

## Declare equilibrium

For methane/ethane, declare both phases and the feed in mol. Here temperature
and pressure are fixed; phase amounts and compositions are solved natively.
See the [equilibrium notebook](examples/methane-ethane-phase-equilibrium.ipynb)
and [equilibrium algorithm](theory/algorithms.md#alg-equilibrium) for the
boundary and flash declarations and their numerical checks.

```python
from epcsaft import Mixture, equilibrium as eq

parameters = Parameters.from_packets("data/packets", [
    ("gross-2001-methane-ethane/3", ["methane", "ethane"])])
mixture = Mixture(parameters)
problem = eq.Problem(
    phases=[eq.Phase("L", kind="liquid"), eq.Phase("V", kind="vapor")],
    T=250.0, P=3_750_000.0,
    feed=eq.Amounts({"methane": 0.3, "ethane": 0.7}))
result = eq.solve_equilibrium(mixture, problem)
assert result.temperature is not None and result.pressure is not None
assert result.amounts and result.coordinates
assert max(abs(v) for v in result.classified_residuals) < 1e-8
print(result.temperature, result.pressure)  # K, Pa
```

Inspect the convergence record and component balances as the notebook does.
A finite result and an optimizer termination message have different meanings;
a declared solve does not perform an implicit stability search.

## Read a rejected declaration

A reaction accepts at most one explicit equilibrium target. Omitting it selects
species thermochemistry; supplying two targets is an error.
This small example runs as a documentation doctest and checks the declaration
exception before a thermodynamic solve:

```{doctest}
>>> from epcsaft import ParameterError
>>> from epcsaft.equilibrium import Reaction
>>> try:
...     Reaction({"A": -1.0, "B": 1.0}, equilibrium_constant=1.0, affinity_j_per_mol=0.0)
... except ParameterError as error:
...     print(error)
reaction admits at most one equilibrium_constant, affinity_j_per_mol or correlation
```

For failures after admission, retain the returned cause and residuals; see
[diagnostics](guide/diagnostics.md) and
[termination rules](theory/algorithms.md#alg-termination).

For signatures and result fields use the [model reference](api/eos.rst),
[equilibrium declarations](api/equilibrium.rst) and
[native equilibrium records](api/native/equilibrium.rst).

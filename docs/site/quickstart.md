# First calculations

[Install ePC-SAFT](installation.md) and start Python in an empty directory.
The saturation and flash examples use the packets installed by conda.
For a source installation, run from the extracted source root and use the
packets already in `data/packets`. Skip the NaCl download in that case.

## Pure propane saturation

At 300 K, solve for the liquid–vapor saturation pressure in Pa. Check
`success` before using the returned pressure.

```{doctest}
>>> from pathlib import Path
>>> import epcsaft
>>> from epcsaft import Mixture, Parameters, Pure, equilibrium as eq
>>> packets = Path(epcsaft.__file__).resolve().parent / "data" / "packets"
>>> if not packets.is_dir():
...     packets = Path("data/packets")  # source installation
>>> packet = "feos-gross-sadowski-2001/1"
>>> propane = Pure(Parameters.from_packets(packets, [(packet, ["propane"])]))
>>> saturation = eq.solve_equilibrium(propane, propane.saturation(300.0))
>>> assert saturation.success
>>> print(f"{saturation.pressure:.2f} Pa")
998660.90 Pa
```

At 300 K and 100 kPa, the same propane model has no ideal-gas thermal
record. Total enthalpy is unavailable. See
[caloric data availability](guide/parameters-and-models.md#caloric-data-availability).

```{doctest}
>>> state = propane.state(300.0, P=100_000.0, phase="vapor")
>>> print(state.total_enthalpy, state.thermochemistry_status.name)
None MissingRecord
```

## Methane and ethane flash

Use `eq.tp_flash` at fixed temperature and pressure. Here the feed is
0.3 mol methane and 0.7 mol ethane at 250 K. The unlike interaction
coefficient is assumed to be $k_{ij}=0$. Phase amounts are in mol, in
liquid (`L`) then vapor (`V`) order.

```{doctest}
>>> components = ["methane", "ethane"]
>>> mixture = Mixture(Parameters.from_packets(
...     packets, [(packet, components)], k_ij={("methane", "ethane"): 0.0}))
>>> feed = eq.Amounts({"methane": 0.3, "ethane": 0.7})
>>> flash = eq.tp_flash(mixture, T=250.0, P=3_750_000.0, feed=feed)
>>> assert flash.success
>>> print(f"L: {flash.amounts[0]:.6f} mol, V: {flash.amounts[1]:.6f} mol")
L: 0.835341 mol, V: 0.164659 mol
```

At 100 kPa the two declared phases coalesce. The final residuals meet the
requested tolerance but the declared split is refused. Check `success`
before reading amounts. A declared flash tests the specified phase pair.
Use [phase discovery](guide/stability-and-discovery.md) to search other phase sets.

```{doctest}
>>> refused = eq.tp_flash(mixture, T=250.0, P=100_000.0, feed=feed)
>>> print(refused.success, refused.topology_event.name, refused.amounts,
...       refused.requested_tolerance_met)
False PhasesCoalesced [] True
```

### Export a flash as JSON

Select the fields to retain. Record the input conditions so they are
available even when the split is refused. The phase names below follow
`tp_flash`'s liquid/vapor declaration. The native message explains the refusal.

```{doctest}
>>> import json
>>> def flash_json(result, *, component_ids, temperature_k, pressure_pa):
...     return json.dumps({
...         "components": list(component_ids),
...         "phases": ["L", "V"],
...         "temperature_k": temperature_k,
...         "pressure_pa": pressure_pa,
...         "amounts_mol": result.amounts if result.success else [],
...         "success": result.success,
...         "refusal_cause": None if result.success else result.message,
...     }, allow_nan=False)
>>> encoded = flash_json(flash, component_ids=mixture.component_ids,
...                      temperature_k=250.0, pressure_pa=3_750_000.0)
>>> accepted = json.loads(encoded)
>>> print(accepted["components"], accepted["phases"], accepted["success"])
['methane', 'ethane'] ['L', 'V'] True
>>> assert accepted["amounts_mol"] == flash.amounts
>>> assert accepted["refusal_cause"] is None
>>> encoded = flash_json(refused, component_ids=mixture.component_ids,
...                      temperature_k=250.0, pressure_pa=100_000.0)
>>> rejected = json.loads(encoded)
>>> print(rejected["success"], rejected["amounts_mol"])
False []
>>> assert rejected["refusal_cause"] == refused.message
```

## Aqueous NaCl activity

```{include} ../../README.md
:start-after: <!-- nacl-example-start -->
:end-before: <!-- nacl-example-end -->
```

## Where to go next

- [Units and independent quantities](guide/states-and-properties.md#choose-the-independent-quantities)
  for temperature, pressure, density and composition.
- [Equilibrium stopping criteria](theory/algorithms.md#alg-termination)
  for interpreting residuals and `Problem(row_tolerance=...)`.
- [Warm starts](guide/declared-equilibrium.md#warm-starts-for-ordered-states)
  for `solve_equilibrium(..., start=...)`.
- [Tracing](guide/diagnostics.md#observe-work-explicitly)
  for `epcsaft.trace(diagnostics=True)`.

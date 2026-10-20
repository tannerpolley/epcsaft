# Examples

These notebooks read canonical inputs from a repository checkout and execute
with its prepared Engine environment. The wheel alone does not contain the
packet inputs. Run `nice -n 10 bash tools/build-docs.sh` from the repository root
to execute every example freshly and build this site. See the
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

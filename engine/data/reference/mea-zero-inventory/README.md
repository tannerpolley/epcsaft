# Exact zero-carbon inventory reference

These frozen declarations reproduce the fresh-solvent numerical failure in
[Engine #231](https://github.com/tannerpolley/ePC-SAFT/issues/231). They qualify
the Engine's conservation face and finite-addition enthalpy chain; they do not
adopt these application parameters or validate calorimetry.

The files were copied without modification from the Amine-Thermodynamics #147
temperature-heat diagnosis inputs. The original five reaction laws and source
references remain declared at both endpoints. Component order is CO₂, MEA,
water, MEAH⁺, MEACOO⁻, HCO₃⁻, CO₃²⁻, H₃O⁺, OH⁻. The initial liquid contains
1 mol MEA and 7.911062165006594 mol water at 313.15 K and 1,734,000 Pa, with
exactly zero added carbon. The final liquid adds 0.11 or 0.45 mol CO₂; pure CO₂
enters on the vapor density branch at the same temperature and pressure.

| Declaration | SHA-256 |
|---|---|
| `parameters.json` | `7e9e3f02723700f8b6577681173249c0e7bd2c8c2a098f38a4e0a7e4611852d5` |
| `reactions.json` | `c1cb5c74b2d71ce56e9c0a5d86116f972adbb2d635a0e9ed944152f25bd5fc43` |
| `ideal-calorics.json` | `c94d504e68c6dc6e85d3d950cc51422499846729b1a06994b362040be991924a` |

The declaration reader
[`zero_inventory.py`](../../../../analyses/2026-addition-enthalpy-fitting/scripts/zero_inventory.py)
sets the diagnosis's explicit Born shell and dielectric scales to one. It
exports the native parameter document, original reactions and ideal caloric
coefficients with `--native-input`. The native reader
[`zero_inventory.cpp`](../../../../analyses/2026-addition-enthalpy-fitting/scripts/zero_inventory.cpp)
(compiled and run by `zero_inventory.py --native-run`) owns coefficient conversion, equilibrium re-solves,
finite-addition values and all derivative-verification arithmetic. Python
declares the model and retains native results; it never imports downstream
source. The source bundle's historical qualification labels remain part of
the frozen input and are not a new Engine qualification claim.

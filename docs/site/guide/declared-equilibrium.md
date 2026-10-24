# Declared equilibrium

Declare the phases and species that may participate, the feed basis, and the
fixed or free intensive quantities. The native Engine solves the resulting
balance and chemical-equilibrium equations through one shared driver.
The [methane/ethane notebook](../examples/methane-ethane-phase-equilibrium.ipynb)
shows a bubble point and a flash. The
[equilibrium formulation](../theory/equilibrium.md) and
[equilibrium algorithm](../theory/algorithms.md#alg-equilibrium) define the rows.

## Phases, supports and guesses

A `Phase` has a name, kind, amount declaration and optional species support.
Support means the species with coordinates in that phase; omitted species
have structural zero amount. It differs from an initially small fraction.
Use explicit supports for restrictions that belong to the physical problem.
The compiler determines the admitted default support when none is supplied.

A liquid can contain charged species only with phase electroneutrality.
A neutral vapor support does not acquire ions from their presence in another
phase. Composition guesses initialize a solve; they neither prescribe the
final composition nor prove that the chosen phases are stable.

Use component-name mappings when preparing guesses to avoid confusing input
order with native support order. Returned phase compositions can be read with
`phase_mole_fractions(phase)` in the model's component order. See the
[declaration reference](../api/equilibrium.rst) and
[native result records](../api/native/equilibrium.rst).

## Feed basis and amount constraints

| Declaration | Supplied basis |
|---|---|
| `Amounts` | species amounts, mol |
| `MoleFractions` | dimensionless mole fractions |
| `MassFractions` | dimensionless mass fractions; molar masses are required |
| `Molalities` | solute species mol/kg, declared solvent set and solvent mass in kg |

A salt formula feed and an explicit-ion feed have different total species
amounts. Keep formula-unit molality distinct from the sum of ionic molalities;
see [formula and explicit-ion feeds](../examples/index.md#formula-and-explicit-ion-feeds).

`Free()` leaves a phase amount or intensive quantity to be solved; for temperature or pressure, `Free(guess)` additionally supplies an initial value in K or Pa, respectively. `Pinned(value)` fixes a
phase amount condition. With an absolute-amount feed, phase amounts are mol;
with a normalized fraction feed the amount constraints follow that normalized
basis. Do not treat a phase fraction as an absolute amount.

```{doctest}
>>> from epcsaft.equilibrium import Free, Pinned
>>> temperature = Free(250.0)
>>> temperature.guess
250.0
>>> Pinned(0.0).value
0.0
```

## Boundary and flash specifications

A fixed-T/P flash declares both intensives numerically and lets the phase
amounts change. A bubble declaration pins the incipient vapor amount at zero;
a dew declaration pins the incipient liquid. `bubble_point` and `dew_point`
construct those declarations at fixed temperature with pressure free.
`Pure.saturation(T)` constructs a pure-fluid liquid/vapor boundary problem.
These are specifications for the shared driver, not separate Python solvers.

Check the returned temperature (K), pressure (Pa), amounts (mol on an amount
feed), compositions, balances and classified residuals. Nonempty state arrays and finite temperature/pressure identify an available
candidate; the achieved convergence record determines whether
it meets the requested tolerance. A declared solve performs no unrequested
tangent-plane trials. Explicit tests belong to
[stability and discovery](stability-and-discovery.md).

## Continuation and caloric specifications

`continue_equilibrium` traces from a starting declaration toward a target,
using native prediction and correction. It records accepted states, retries,
reached path and the stopping cause. An earlier accepted state does not mean
the requested endpoint was reached. Ordinary path fractions differ from
logarithmic envelope coordinates; inspect the continuation records.

Admitted neutral PH and PS problems keep pressure fixed and temperature free,
with declared free phase amounts. Supply one extensive target:
`enthalpy_j` in J or `entropy_j_per_k` in J/K, and thermal records for every
consuming species. Charged and source-referenced PH/PS declarations are
unavailable; UV is deferred. These limits follow the
[caloric formulation](../theory/equilibrium.md) and are collected in
[limits](../limits.md); the TP notebook does not qualify caloric specifications.

## Pure stoichiometric solids

A `Solid` declares its formula, source equilibrium constant and reference,
positive molar volume in m³/mol, and amount condition. The native solve
handles absent/present complementarity together with fluid balances.
See [halite complementarity](../examples/index.md#halite-complementarity)
for the one-salt numerical example. See
[solid equations](../theory/algorithms.md#pure-stoichiometric-solids) and the
[reaction guide](reacting-equilibrium.md) before changing reference bases.

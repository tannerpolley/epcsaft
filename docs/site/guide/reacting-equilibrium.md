# Reacting equilibrium

Supply the reaction set and its thermodynamic standards explicitly. The Engine
has no application-specific reaction list or equilibrium-constant default.
The [ammonia notebook](../examples/ammonia-synthesis-equilibrium.ipynb) shows
one homogeneous vapor reaction with published pure parameters; the
[halite notebook](../examples/halite-solubility.ipynb) shows a dissolution
reference and solid complementarity. Their equations belong to the
[equilibrium formulation](../theory/equilibrium.md) and
[reference algorithm](../theory/algorithms.md#alg-reference).

## Declare stoichiometry and the conserved inventory

A `Reaction` maps species to signed stoichiometric coefficients: negative for
reactants and positive for products. Every species must belong to the model.
The reaction set must conserve charge. Supply a physically conserved feed;
check element or moiety balances independently when interpreting a solution.
The native compiler derives the conserved coordinates from the declared
stoichiometry. It does not infer missing chemistry from species names.

Declare at most one explicit target for each reaction: an equilibrium constant,
an affinity in J/mol, or a supported temperature correlation. With no target,
the compiler derives the reaction standard from the model's species
thermochemistry at the fixed temperature and pressure. Every participating
species needs an interval at that temperature; a missing interval remains
unavailable. Caloric conditions and source-reference transfer for this choice
are deferred. A positive
constant K is dimensionless on its stated standard; it is not ln K.
The [reaction declarations](../api/equilibrium.rst) and
[native reference records](../api/native/equilibrium.rst) define the fields.

```{doctest}
>>> from epcsaft.equilibrium import Reaction
>>> reaction = Reaction({"N2": -1.0, "H2": -3.0, "NH3": 2.0}, equilibrium_constant=1.0)
>>> reaction.stoichiometry["NH3"]
2.0
```

This is only a declaration example. The numerical constant 1 above is not an
ammonia correlation; use the sourced constant and component names in the
ammonia notebook for its calculation.

## Keep the source basis with K(T)

`ReactionLogPolynomial` evaluates

$$
\ln K(T)=a+b/T+c\ln(T/T_{ref})+dT.
$$

Temperature and the declared reference temperature are K. The coefficients
have units that make ln K dimensionless: b is K, c and a are dimensionless,
and d is K⁻¹. Supply the correlation's admitted temperature interval and its
standard-state identity. A base-ten reciprocal-temperature source must be converted explicitly to natural-log coefficients before constructing `ReactionLogPolynomial`. Do not pass log10 K as ln K.

A constant or correlation on a source molality or infinite-dilution basis
also needs a `ReactionReference`. Declare the standard concentration or
molality, solvent composition, pressure and ionic gauge where consumed.
Pure-liquid, pure-gas and infinite-dilution species standards are different
physical references. The Engine computes the standard-state transfer through
native thermodynamics; it does not silently reinterpret a source constant.

## Separate fixed-temperature chemistry from caloric consistency

See the [ammonia fixed-temperature standard](../examples/index.md#ammonia-fixed-temperature-standard)
for an explicit scalar conversion.
A fixed scalar does not reproduce the temperature derivatives of independent
thermal source records. Supplying inconsistent reaction and caloric standards
is refused rather than choosing one silently.

For temperature-dependent or caloric work, supply compatible K(T), ideal
thermal records and reference anchors. Residual chemical potentials alone
do not establish reaction enthalpy or total entropy. Read the
[caloric algorithm](../theory/algorithms.md#alg-calorics) and
[limits](../limits.md) before requesting reacting thermal outputs.

## Read the result as a bounded prediction

Inspect the balances, achieved residuals and reference basis using the
[declared-equilibrium guide](declared-equilibrium.md),
[diagnostics](diagnostics.md) and
[ammonia and halite interpretation](../examples/index.md#ammonia-and-halite-interpretation).

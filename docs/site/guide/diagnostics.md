# Diagnostics

Read the cause, available values and achieved residuals together. Declaration
errors, unavailable candidates and finite candidates outside a requested
tolerance describe different outcomes. The
[termination algorithm](../theory/algorithms.md#alg-termination) defines these
meanings. The [water notebook](../examples/parameters-and-water-state.ipynb),
[caloric notebook](../examples/methane-caloric-properties.ipynb) and
[discovery notebook](../examples/water-isobutanol-discovery.ipynb) show pressure
closure, missing thermal outputs and finite-search evidence respectively.

## Start with the declaration

A `ParameterError` identifies an invalid or incomplete declaration, such as
an unknown component, an unsupported unit or specifying both P and density.
An `EosError` can identify a thermodynamic evaluation failure. Use the
[Python declarations](../api/epcsaft.rst) and
[native state records](../api/native/eos-state-reference.rst) to distinguish
exceptions from returned availability causes.

Keep the smallest physical input that reproduces the cause: source selection,
component order, T in K, P in Pa or density in mol/m³, composition, support,
reference and requested observable. A raw optimizer message without those
inputs does not identify a thermodynamic failure.

## Inspect candidates and convergence

For equilibrium, check that temperature and pressure are present and the
state arrays are populated before consuming the candidate. Then inspect
`convergence`, `message`, raw and classified residuals and warnings.
A finite candidate can miss the requested tolerance; retain its achieved
residual instead of reporting a successful physical result from the optimizer
termination alone. No finite candidate means there is no state to substitute
into a downstream calculation.

Balances use the feed's declared basis. A pressure residual has units of Pa;
a reported scaled/classified row uses its own normalization. Compare a row
with the tolerance defined for that row rather than mixing raw residuals with
unrelated dimensionless thresholds. The
[native equilibrium reference](../api/native/equilibrium.rst) describes the
result fields.

A numerical floor estimates attainable accuracy from the reported inner-solve
and input errors. It is not automatically a rigorous error bound or a reason
to silently relax a requested tolerance. Conditioning, residual refinement
and an independent identity provide different evidence.

## Diagnose the density branch

A state at pressure retains its density convergence record and branch
information. `density_diagnostics` requests a bounded report over declared
density nodes, with evaluation order, tolerances, time and evaluation budgets.
It can return samples, pressure crossings, spinodals, declared brackets and
unresolved points. Samples use mol/m³, pressure Pa and pressure-density slope
Pa m³/mol. Read unavailable sample causes before using their fields.

The [density algorithm](../theory/algorithms.md#alg-density) defines the
finite search and connected-cell meanings. A finite scan does not enumerate
all roots by proof. The water notebook checks its selected branch's pressure
closure; it does not demonstrate complete root enumeration.

## Preserve partial availability

A missing ideal thermal record leaves Total calorics unavailable while valid
residual work remains readable. A singular action can leave a central value
available. Inspect `thermochemistry_status`, component-property causes and
per-action causes separately; never replace an unavailable value with zero.
The [state guide](states-and-properties.md) and
[action guide](solved-state-derivatives.md) explain those distinctions.

## Observe work explicitly

`epcsaft.trace(diagnostics=True)` is a context manager for explicit native
observations. Run the calculation inside it and read the captured `spans` and
`events` afterwards. A span can retain its convergence record; an event names
its iteration, values and reason. Model counters record derivative work;
continuation and discovery have their own accepted-state/trial records.
Tracing is an explicit observation, not permission to rerun a hidden search
or to change tolerances. The
[trace definitions](../api/native/eos-state-reference.rst) and
[equilibrium algorithm](../theory/algorithms.md#alg-equilibrium) describe the
available records.

When reporting a problem, include requested tolerance, achieved residuals,
physical inputs and the exact unavailable cause. State whether the evidence
is a source-parameter prediction, numerical verification or an independent
physical comparison. A convergence record alone establishes none of the
last category.

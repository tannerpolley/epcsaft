# Stability and phase discovery

A converged declared split satisfies the equations for the phases you allowed.
It does not establish that another phase cannot lower the Gibbs energy.
Request a stability test or fixed-T/P phase discovery explicitly. The
[water/isobutanol notebook](../examples/water-isobutanol-discovery.ipynb)
compares a declared split with discovery and inspects finite trials. See the
[stability algorithm](../theory/algorithms.md#alg-stability) and
[equilibrium formulation](../theory/equilibrium.md) for trial definitions.

## Test a declared result

`stability_test` evaluates a solved snapshot at its solved temperature and
pressure with the actual standards frozen. Trial supports, starting points
and work limits determine which possibilities are examined. The driver
performs no implicit trial searches during ordinary declared solves,
continuation, parameter fitting or ASL evaluation.

Tangent-plane distance measures reduced molar energy above the candidate's
conserved supporting plane. Inspect trial feasibility, density availability,
charge balance and energy descent together. An optimizer message alone does
not establish a feasible lower-energy phase. The native records distinguish
unresolved trials and numerical uncertainty from an admitted descending split.
The [equilibrium reference](../api/equilibrium.rst) and
[native stability records](../api/native/equilibrium.rst) expose these calls
and evidence fields.

## Discover at fixed temperature and pressure

A `Problem` with `phases=Discover()` asks the native search to consider topology
changes at fixed T/P. Supply a feed and explicit search budgets. Optional trial
declarations restrict the starts and supports examined. The search uses the
shared equilibrium driver for admitted subproblems and retains its work and
stopping causes.

`Discover(mode="fast")` runs a smaller search: it stops at the first resolved
split and solves one two-phase declaration, without a final stability search.
Its record lists the start families and topologies it left unsearched. Use the
default full mode when a fast result must be checked for a further split.

See [water/isobutanol discovery](../examples/index.md#waterisobutanol-discovery)
for the bounded numerical comparison.

Once a topology is found, it may seed later declared calculations within a
qualified domain. Reuse preserves the phase supports and reference choices;
a previous discovery does not remove the need to check a new state's
residuals, branch identity and requested stability evidence.

## Charged trial admission

Charged trials must be electroneutral and satisfy the permittivity model's
domain. Density-dependent forms with the vacuum limit and density-independent
condensed-phase forms have different trial admission rules. Under the latter,
a charged trial counts only at a stable root on the liquid side of its own
composition's mechanical loop. An off-root, vapor-side or loop-less trial can
remain unavailable with an explicit cause.

See [electrostatic domains](../theory/eos.md) and the
[stability algorithm](../theory/algorithms.md#alg-stability). The neutral
notebook does not qualify every charged or reacting search. A refused trial
is not evidence of the physical absence of that phase.

## What a finite search establishes

A resolved feasible Gibbs descent is evidence against the candidate's
stability, within the reported numerical uncertainty. Finding no such descent
in a finite set of local trials is weaker: retain the number of trials,
unresolved attempts, tested supports and exhausted budgets.

Local stability searches are not global proofs. Positive mechanical slopes
and local curvature also do not exclude a distant lower-energy split. A
critical-equation candidate satisfies necessary equations; it is not by itself
a stability certificate. Use the wording supported by the actual search and
read the [limits](../limits.md) before promoting a phase prediction.

## Measure coverage on a declared composition grid

Call `stability_grid(compiled, parent, phases, denominator)` for an accepted
fixed-temperature, fixed-pressure parent. Each `PhaseSpec` declares a support
and a liquid or vapor density side. The Engine enumerates mixtures of the pure
neutral species and neutral opposite-sign ion pairs with integer vertex weights
summing to the denominator. Exact boundary nodes use their positive species
support. The default ceiling of 2,048 raw attempts applies across the request
before deduplication. Set `ceiling` explicitly to allow a larger schedule.

For each support and side, `h` is the largest species mole-fraction L1 distance
from any planned grid node to its nearest admitted raw node. It is zero when
every node is admitted and unavailable when the side has no admitted nodes.
Wrong-side or unresolved pressure roots, refused evaluations, and unvisited
nodes can increase the distance. Shared faces belong to each relevant support.
The rounding bound `min(2, k / denominator)`, for `k` vertices, covers rounding
from the represented continuous composition domain to the planned grid. Add
that bound, `h`, and the separate mapping allowance to bound distance to an
evaluated composition. A sum of at least 2 gives no improvement over the
simplex diameter. The record identifies omitted species and retains the causes
of individual gaps. Each raw trial retains its density and EOS work, including
refused evaluations. The raw reduced molar TPD estimate treats the parent
plane as fixed. Donor-slope and finite-withdrawal errors are separate quantities.

Coverage measures composition distance only. It supplies neither a bound on
unsampled tangent-plane distances nor a count of all density branches, so a
finite search without a resolved Gibbs decrease cannot certify stability.
Optional local refinement and finite withdrawals can supply witnesses without
changing raw grid coverage. Ordinary stability tests, phase discovery, fits
and flashes leave this diagnostic unused.

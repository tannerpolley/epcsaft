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

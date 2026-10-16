---
status: accepted
date: 2026-09-06
---

# General platform with one equilibrium driver

Kind: computational-architecture

The owner selected a general fluid thermodynamics platform whose first required
acceptance uses are the three official applications. The owner explicitly
requires one equilibrium driver and selected Git history as the sole retention
location for obsolete setup and execution documents.

The [project map](../CONTEXT.md) defines the current responsibilities and scope.
All admitted equilibrium specifications use one driver and shared mathematical
rows; separate application/problem-specific equilibrium solvers are excluded.
This replaces the proposal to allow independently selected equilibrium drivers.
It preserves lower-level density/association kernels and the distinct search
and continuation operations that reuse equilibrium calculations. The one bounded neutral recovery
after a failed ordinary L+V solve ([ALG-EQUILIBRIUM](../../../engine/docs/science/algorithms.md#neutral-recovery), #565)
reuses this driver and the straight continuation controller.

The existing native ownership principles remain: one Python-free native shared library for Python and native programs,
one extension and an ASL adapter with its private native copy from the same PIC
objects (#485 owner decision, 2026-10-14; supersedes archive-only linkage), no public C++ ABI or iterative Python callbacks, explicit
physical inputs, common parameter identities and owner-supplied derivative
chains. General platform scope does not introduce compatibility support,
external catalogs, speculative extension frameworks or a second installation.

On 2026-09-07 the owner specified one Ipopt driver with exact curvature for
all admitted equilibrium specifications, including bubble/dew. Root or
limited-memory equilibrium replacements are not authorized. This requirement
does not qualify the existing Hessian: pinned/free-intensive formulation,
derivative correctness and demanded-case cost remain work in the implementation
plan. Lower-level density/association root kernels remain distinct. Material limitations require an owner
decision rather than a silent second driver. Solids remain a separately
justified extension. Runtime budgets retain their current numeric owner.

Obsolete stage designs, briefs, worker records and superseded architecture
ADRs are preserved in Git at `2985cd4`; current scientific evidence remains
with its owner. No historical worker roles, launch commands or per-stage
budgets govern future work. The installed CSE workflow and repository guidance
own execution.

**Amendment 2026-10-13 (#476, owner).** For eligible calls entering `solve_compiled`, the driver first
applies its exact-Jacobian Newton iteration to the shared rows from the same start, without constructing
Ipopt.

**Eligible calls.** Square feasibility cells with one liquid and one vapor phase, where:

- exactly one of T and P is free;
- one phase is pinned at exactly zero amount and the other is free;
- there is no phase discovery, no reaction, no ion, and no caloric, critical-equation or envelope condition.

These are saturation, bubble and dew points.

**Callers.** Public solves, fit observations and the declared solve inside `cooling_onset` enter
`solve_compiled`. The following keep the Ipopt driver alone:

- declared starts (`solve_equilibrium_starts`);
- continuation correctors (`solve_equilibrium_from_coordinates`);
- stability trial solves;
- caloric seed sub-solves.

**The iteration.** It is the driver's existing polish: the CppAD row Jacobian, a full-pivot LU, damping
1, 1/2 and 1/4 with raw-row contraction, and the existing step, correction and iteration limits.

**Lead amendment 2026-10-13 (#476).** Once the existing residual criterion is met, apply exactly
one more full exact Newton correction, with a fresh Jacobian and solve at that point, and return
the corrected point. It must lie inside the same explicit candidate box and still meet the
residual criterion; otherwise the original start falls back to Ipopt. Retain the correction's
evaluation and Jacobian work. No tolerance, comparator or runtime policy changes.

**Lead amendment 2, 2026-10-13 (#476).** At every Newton
candidate box check, a log-composition coordinate above its upper bound of zero by at most
four ulps of the physical mole-fraction bound of one is set exactly to zero. The allowance
is computed as `log(xmax)`, where `xmax` is four successive `nextafter` steps above one
(8.881784197001248e-16 in the log coordinate). Re-prepare the snapped point and recheck
the unchanged residual criterion; a failure falls back to Ipopt. Keep the additional exact
correction and its retained work. Every other bound, tolerance and comparator is unchanged.

**When its point is returned.** Only when all of these hold:

- the iteration ends before its iteration cap;
- every coordinate lies inside the bounds the driver gives Ipopt;
- the classified rows meet the unchanged caller criterion, without a floor allowance;
- the unchanged topology, local phase stability and density-branch admission checks accept it.

**Fallback.** In every other outcome, including an exception, the same start runs the Ipopt driver
unchanged. Ipopt and its checks remain the only source of failures, floor-limited states and refusals.

**What stays the same.** Newton on a square system uses no Lagrangian Hessian, so the exact-curvature
requirement is unchanged wherever Ipopt runs. No finite-difference, secant, quasi-Newton or derivative-free
step is admitted. Solved-state actions read the same coordinates and rows.

**Scope.** #476's equivalence corpus is the qualification. Widening the eligible calls needs a new owner
decision; examples are ions, vapor-fraction pins, and other callers.


**Amendment 2026-10-15 (#545, owner).** Ordinary `solve_compiled` calls may
also attempt the same exact-Jacobian Newton iteration on neutral, nonreacting,
fixed positive T/P cells with exactly two free-amount liquid/liquid or
liquid/vapor phases. Both supports must equal the retained positive-feed
support and contain at least two components; every model charge number must
be zero. Source-reference conversion, caloric, critical, envelope, discovery,
diagnostic-options and optimizer-transfer calls are excluded. The #476
predicate, iteration and admission behavior remain unchanged. Declared starts,
continuation correctors, stability trials and caloric seed sub-solves retain
Ipopt alone.

The existing seed owner constructs the start once. A copy receives Newton;
the original remains the Ipopt fallback. Only this free-amount attempt checks
finite coordinates, the physical-coordinate primal box and distinct,
non-disappeared seed topology before preparing iterations, and full numerical
Jacobian rank and finite corrections before trials. It reuses the actual
Ipopt bounds with the amount chart converted to physical mol and the existing
four-ulp composition snap. Only this attempt permits a trial already meeting
the caller's raw-row criterion past the previous-contraction-ratio test.
The exact CppAD Jacobian, full-pivot LU, damping and iteration/correction
limits, compulsory full correction and fresh classified criterion remain.

Return requires unchanged positive-amount/topology admission, exact reduced
Lagrangian-Hessian positive curvature on the constraint tangent space and
density-branch admission. Equal compositions can occur at a multicomponent
azeotrope; rank and strict reduced curvature must still resolve the amount
split. No global stability or uniqueness certificate is implied. Every rejected
attempt, admission failure or exception falls back once to unchanged Ipopt
from the original seed, retaining all attempted work and the caller's shared
budget. Ipopt remains the sole failure, refusal and floor-limited-state owner.
Solved-state actions use the same accepted coordinates, rows and derivatives.
#545's numerical equivalence and timing evidence qualify this bounded widening;
further widening requires an owner decision.

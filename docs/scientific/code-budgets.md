# Engine code budgets

The owner requested a hypothetical floor as the design goal, with soft and hard
ceilings to constrain duplication and unnecessary machinery. Here **floor
target** means an engineering estimate of a lean complete implementation. It is
not a proven minimum, an allocation to spend or a lower acceptance bound.

The single numeric table is `ARCHITECTURE.yaml:code_budgets`. Run
`python3 tools/diffstat.py --budgets` to display current counts against it.
Both normal local check commands enforce it. `check-affected.sh` provides the
documentation-only route without a wheel build; `check-all.sh` always builds
and installs after readiness checks. Numerical evidence is a separate
acceptance condition.

## Counting and scope

Count physical lines, including comments, blank lines and type stubs, in the
complete working `engine/native/` and `engine/src/epcsaft/` trees. Count tracked
and untracked source alike, regardless of diff attributes. Report native layers
and Python separately, plus their total. Reject unknown native layers, runtime
file types and symlinks rather than silently dropping them from the count.
Tests, documentation, source data, build outputs and third-party dependencies
are outside this runtime budget. They remain subject to scope/ownership review;
none may hide runtime logic to evade a ceiling.

Physical lines are a reproducible size proxy, not a measure of scientific
complexity. Keep conventional readable formatting. Removing useful comments,
packing statements onto lines, hiding code in generated strings, relocating
runtime code, or splitting/renaming files to evade a budget is unacceptable.
Spare lines never justify a second thermodynamic or solver owner.

## Initial estimate and uncertainty

The September 6 baseline is 8,807 lines in 37 source files. The target allocation
starts with the existing Core, Contributions and admission owners. Properties allows
roughly 675 additional lines for reference thermochemistry, property chains
and compact derivative actions. Equilibrium allows roughly 1,820: historical
planning estimated about 1,400 for stability, with the balance assuming reuse
of the compiler, rows and factorizations for continuation and derivatives.
This historical equilibrium allowance includes search work now deferred by the
implementation plan; it is not permission to build that work or spend the spare
lines. All numeric ceilings stay unchanged. Before each slice, inspect its
current/old owners and justify only the bounded repair, adaptation or missing
capability. Verify-only slices begin with zero new runtime. Reuse decisions and
measured equivalent-work costs, not line count alone, justify replacement.

Regression's target follows the retained small fitting-wrapper and four-node
observation design. Binding/Python growth assumes replacement of interface
pieces as all applications migrate; there is no compatibility layer allowance.
The allocation covers the currently documented fluid transfer, including its
required actions, not solids or unrestricted dense derivative tensors. These
estimates are not evidence that the missing science already fits or works.

## Floor, soft ceiling and hard ceiling

- Work toward the floor target. Fewer lines never prove completeness; do not
  add code just to reach it.
- Above a soft ceiling, the command warns. Review each added responsibility,
  reused owner and removable duplication; retain a short reason in the change
  description before acceptance. The soft ceiling grants room to finish and
  simplify, not automatic budget growth.
- Above any layer or total hard ceiling, readiness fails. Refactor/reduce or
  obtain an explicit owner-approved budget amendment with the scientific
  requirement, alternatives and revised numeric table. No override flag or
  automatic ratchet follows the current line count.
- A temporary over-budget working state may be simplified and repaired. Do not
  freeze uncompilable code because an intermediate write exceeded a number.
  The completed candidate must satisfy the hard ceiling before acceptance.

Never satisfy a budget by dropping required physics, domain checks, derivative
terms, failure diagnostics or numerical falsifiers. If the smallest correct
method does not fit, report that evidence and revise the design/budget openly.
Budget changes must not silently change scientific scope.

## Duplication constraints

One common Helmholtz definition supplies properties and their derivatives.
Association elimination, density closure, reference transformations, equilibrium
row assembly and physical-input resolution each retain their accountable owner.
Python maps inputs/outputs; applications supply chemistry and specifications as
data. All equilibrium specifications use one driver; association, density,
fitting and discovery retain their distinct numerical operations using the
shared owners. Existing architecture checks cover known
include/constant/call patterns; size and regex checks cannot establish absence
of semantic duplication. Scientific/code review must inspect the actual calls.

This budget is a local readiness rule, not CI/CD and not a numerical accuracy or
runtime-speed threshold.

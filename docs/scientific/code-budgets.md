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

Count physical lines, including comments and blank lines, in the complete
working `engine/native/` and `engine/src/epcsaft/` trees, except the type stub
`engine/src/epcsaft/_core.pyi`. The stub declares the compiled extension's types
and holds no Python logic; `tools/check-engine-architecture.py` checks that its
top-level names and class members match the installed extension. Count tracked
and untracked source alike, regardless of diff attributes. Report native layers
and Python separately, plus their total. Reject unknown native layers, runtime
file types and symlinks rather than silently dropping them from the count.
Tests, documentation, source data, build outputs and third-party dependencies
are outside this runtime budget. They remain subject to scope/ownership review;
none may hide runtime logic to evade a ceiling.

Each total entry may be tighter than the sum of its layer entries, never looser:
a layer ceiling bounds that layer, and the total independently bounds the whole.

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
owning GitHub issue; it is not permission to build that work or spend the spare
lines. All numeric ceilings stay unchanged. Before each slice, inspect its
current/old owners and justify only the bounded repair, adaptation or missing
capability. Verify-only slices begin with zero new runtime. Reuse decisions and
measured equivalent-work costs, not line count alone, justify replacement.

The parameter regression target follows the retained small fitting-wrapper and
four-node observation design. Binding/Python growth assumes replacement of interface
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

## Amendments

One line per approved ceiling change, newest first (hard ceilings unless marked
soft). The issue or PR holds the requirement, alternatives and measured
accounting; Git keeps the former per-slice tables. Entries marked
**inherited** come from [the `a0adb1cd` architecture comments](https://github.com/tannerpolley/ePC-SAFT/blob/a0adb1cd5ff5f4ee32b3608c78a02c756b06ebd9/ARCHITECTURE.yaml#L31);
their approval assertions were not reverified in the #221 reconciliation.

- 2026-10-07, #341 owner decision: the Python layer and the total no longer count the 1,225-line `_core.pyi` type stub; the architecture check compares its names and class members with the extension instead. The owner then lowered both hard ceilings by those 1,225 lines (`_core.pyi` at base `5cc1001d`): Python **2,892 → 1,667** and total **25,000 → 23,775**. Floor and soft targets are unchanged.
- 2026-10-06, #224 slice B (owner): Equilibrium 8,800 → 8,843, Python 2,850 → 2,892, the measured counts.
- 2026-10-06, #245 (owner): Regression 2,478 → 2,550.
- 2026-10-06, #283 (owner): Core 917 → 1,100, Contributions 1,547 → 1,900, Properties 4,800 → 5,300, total 24,000 → 25,000; soft Core 850 → 1,000, Contributions 1,300 → 1,700, Properties 2,900 → 5,000, Equilibrium 6,000 → 8,000, total 14,850 → 19,500 (the layer soft sum).
- 2026-10-02, #220 envelope raise (owner): Properties 4,600 → 4,800, Equilibrium 7,575 → 8,800, Bindings 1,940 → 2,100, Python 2,700 → 2,850, total 22,000 → 24,000.
- 2026-10-01, #221 finding 5 closeout: Equilibrium 7,650 → 7,575.
- 2026-10-01, #223 (owner): Bindings 1,880 → 1,940, Python 2,620 → 2,700.
- 2026-10-01, #219 slice B closeout: Regression 2,480 → 2,478.
- 2026-10-01, #219 slice B (owner): Regression 2,200 → 2,480.
- 2026-10-01, #219 slice A closeout: Core 1,000 → 917, Contributions 1,612 → 1,547.
- 2026-10-01, #220 roadmap envelope (owner): Properties 4,258 → 4,600, Equilibrium 7,000 → 7,650, Regression 1,619 → 2,200, Bindings 1,726 → 1,880, Python 2,374 → 2,620, total 19,100 → 22,000.
- 2026-09-30, #210 (owner): Properties 4,223 → 4,258, Regression 1,484 → 1,619, Bindings 1,725 → 1,726, Python 2,357 → 2,374, total 18,830 → 19,100.
- 2026-09-30, #211 (owner): Regression 1,463 → 1,484, Python 2,355 → 2,357, total 18,805 → 18,830.
- 2026-09-30, #205 (owner): Regression 1,440 → 1,463, Python 2,350 → 2,355, total 18,790 → 18,805.
- 2026-09-30, #187 (owner): Regression 1,300 → 1,440, Bindings 1,708 → 1,725, Python 2,313 → 2,350, total 18,556 → 18,790.
- 2026-09-29, #129/#193 (owner): Properties 4,211 → 4,223, Regression 1,200 → 1,300, total 18,365 → 18,556.
- 2026-09-26, #172 (owner): total 18,353 → 18,365.
- 2026-09-24, #152 (owner): Bindings 1,700 → 1,708, Python 2,300 → 2,313, total 18,145 → 18,353.
- 2026-09-23, #142 (**inherited**): Properties 4,185 → 4,211, total 18,109 → 18,145.
- 2026-09-23, #140 (**inherited**): Properties 4,183 → 4,185, total 18,102 → 18,109.
- 2026-09-23, #138 (**inherited**): Properties 4,150 → 4,183, total 17,903 → 18,102.
- 2026-09-23, #28 (owner): Contributions 1,500 → 1,612, Properties 4,145 → 4,150, total 17,760 → 17,903.
- 2026-09-23, #84 (**inherited**): Properties 4,100 → 4,145.
- 2026-09-23, #105 (**inherited**): total 17,790 → 17,760.
- 2026-09-23, #87/E4 (**inherited**): total 17,760 → 17,790.
- 2026-09-23, #107–#109 (**inherited**): total 17,601 → 17,760.
- 2026-09-22, #84 (**inherited**): Properties 3,470 → 3,490.
- 2026-09-22, #66 (owner): total 17,600 → 17,601.
- 2026-09-22, #66 (owner): total 17,450 → 17,600.
- 2026-09-21, #89 (owner): Properties 3,460 → 3,470.

Pre-envelope changes without an entry are listed in the #221 reconstruction
comment (5942662477).

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

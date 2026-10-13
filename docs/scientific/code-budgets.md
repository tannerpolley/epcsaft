# Engine code budgets

A hard line ceiling per runtime layer, and one for the total, constrains
duplication and unnecessary machinery. The single numeric table is
`ARCHITECTURE.yaml:code_budgets`. Run `python3 tools/diffstat.py --budgets` to
display each count, ceiling and band. Both normal local check commands enforce
it. `check-affected.sh` provides the documentation-only route without a wheel
build; `check-all.sh` always builds and installs, before readiness checks.
Numerical evidence is a separate acceptance condition.

## Counting and scope

Count the lines of each runtime file as one pinned formatter lays it out,
without changing the file: C++ with `clang-format-18` and the style
`tools/line-count.clang-format`, Python with the `.venv` ruff (`uv.lock`
version) at 120 columns, where a trailing comma does not choose the layout. The
counter's `formatted_lines` holds the exact commands. Comments and blank lines
count as the formatter keeps them; comments are not reflowed. Count the
complete working `engine/native/` and `engine/src/epcsaft/` trees, except the
type stub `engine/src/epcsaft/_core.pyi`. The stub declares the compiled
extension's types and holds no Python logic; `tools/check-engine-architecture.py`
checks that its top-level names and class members match the installed
extension. Count tracked and untracked source alike, regardless of diff
attributes. Report native layers and Python separately, plus their total.
Reject unknown native layers, runtime file types and symlinks rather than
silently dropping them from the count. Tests, documentation, source data, build
outputs and third-party dependencies are outside this runtime budget. They
remain subject to scope/ownership review; none may hide runtime logic to evade
a ceiling.

The total ceiling may be tighter than the sum of the layer ceilings, never
looser: a layer ceiling bounds that layer, and the total independently bounds
the whole.

Formatted lines are a reproducible size proxy, not a measure of scientific
complexity. The formatter splits statements packed onto one line, and the
counter refuses a runtime file holding a formatter-disable directive
(`clang-format off`/`on`, `# fmt: off`/`on`/`skip`, `# yapf: disable`/`enable`),
naming its file and line. Packing still gains lines where the formatter keeps
the packed form: comma-joined C++ declarations such as `int a = 1, b = 2;`, and
a C++ braced list, which shares lines without a trailing comma and takes one
line per element with it. Removing useful comments, hiding code in generated
strings, relocating runtime code, or splitting/renaming files to
evade a budget is unacceptable. Spare lines never justify a second
thermodynamic or solver owner.

## Ceilings and the band

The band of a layer, or of the total, is 1% of its current count, rounded up.

1. **Raise to count plus band.** Above a hard ceiling, readiness fails.
   Refactor or reduce, or raise the ceiling with owner approval to exactly the
   count plus the band, recording the raise below in the same change. No
   override flag or automatic upward ratchet exists.
2. **Lower after consolidation.** When a ceiling sits more than one band above
   the count, readiness fails with the ceiling to lower to (count plus band).
   Lower it in the same change; lowering needs no approval. A ceiling up to one
   band above the count needs no table change.
3. **Approve once per issue.** The owner approves a line allowance per layer
   when authorizing an issue, and the issue body records it. A merge whose
   raises stay within that allowance needs no second approval. Growth beyond
   it needs a new owner decision with the scientific requirement and
   alternatives.

A temporary over-budget working state may be simplified and repaired. Do not
freeze uncompilable code because an intermediate write exceeded a number. The
completed candidate must satisfy both band limits before acceptance.

Never satisfy a budget by dropping required physics, domain checks, derivative
terms, failure diagnostics or numerical falsifiers. If the smallest correct
method does not fit, report that evidence and revise the design/budget openly.
Budget changes must not silently change scientific scope.

## Amendments

One line per ceiling change, newest first (hard ceilings). The issue or PR
holds the requirement, alternatives and measured accounting; Git keeps the
former per-slice tables and the retired floor and soft ceilings. Entries before
#436 count physical lines. Entries marked **inherited** come from
[the `a0adb1cd` architecture comments](https://github.com/tannerpolley/ePC-SAFT/blob/a0adb1cd5ff5f4ee32b3608c78a02c756b06ebd9/ARCHITECTURE.yaml#L31);
their approval assertions were not reverified in the #221 reconciliation.

- 2026-10-13, #435 (rule 2, consolidation; the field-list macro accepted 2026-10-13 as the recommended option while the owner was away, confirmed at merge): Bindings **2,506 → 1,926** and Python **2,050 → 1,989**, the counts 1,906 and 1,969 plus their bands, after `EPCSAFT_READONLY`/`EPCSAFT_READWRITE`/`EPCSAFT_VALUES` (`bindings/fields.hpp`) replaced one line per same-named field or enum value, one aggregate constructor and a typed continuation argument replaced binding glue, and five declaration classes became slotted dataclasses (−576 Bindings, −73 Python, −6 parameter Regression, which stays under its ceiling). Total **24,895 → 24,254**, the sum of the layer ceilings, below the count 24,064 plus its band.
- 2026-10-13, #439 (owner allowance Bindings +1, Python +15; the Python excess accepted 2026-10-13 as the recommended option while the owner was away, confirmed by the owner's merge): Python **2,025 → 2,050**, the count 2,029 plus its band after name-keyed compositions, `molar_masses`, `phase_mole_fractions`, five result reprs and clearer messages (+1 Bindings, +23 Python, +6 Properties, +1 Equilibrium); `mass_basis` required only for `phase_density` (+2), the refusal of unread top-level document fields (+8) and the review corrections (a read-only `molar_masses`, named guesses on the declared support in native component order) bring Python to +35, all within the same ceiling. Its +1 Bindings line fits under the #73 slice 2 ceiling 2,506. Total **24,870 → 24,895**, the sum of the layer ceilings.
- 2026-10-13, #73 slice 2 (owner allowance, rule 1): Bindings **2,283 → 2,506**, the count 2,481 (main 2,283, plus +198 for the `epcsaft_pressure_property` ASL callback and `property_asl_record`) plus its band 25 (+198 Bindings, +1 Python). Total **24,486 → 24,870**, the sum of the layer ceilings, below the count 24,675 plus its band. The #73 cumulative is +313 of about +317 (slice 1 +114, slice 2 +199).
- 2026-10-13, #73 slice 1 (owner): Properties **4,560 → 4,672**, the count 4,625 after merging #441 (main 4,514, plus +111 for the component ln(phi_i) pressure actions) plus its band 47 (+111 Properties, +3 Bindings). The owner raised the slice allowance from +108 to +114 for the dipolar-term and absent-ion status guards; the #73 total is then +316 of about +317. Total unchanged (24,476 / 24,486).
- 2026-10-12, #441 (owner allowance, rule 1): Properties **4,511 → 4,560**, the count plus its band after the dielectric-coefficient families, their admission guard and the reference parameter scale (+23 gross, +15 net formatted lines) on main `42f86c5a`, whose #423 slice 1 used 13 lines of the former headroom. The owner's #441 allowance is Properties 30 plus the scale correction.
- 2026-10-12, #433 (rule 2, consolidation): Equilibrium **9,542 → 9,538**, parameter Regression **3,211 → 3,199**, total **24,499 → 24,486** (the sum of the layer ceilings, below count plus band), the counts plus their bands after one source-reference equality replaced three comparisons.
- 2026-10-12, #436 (owner): re-baseline to formatter-normalized lines, each ceiling the count plus its band; floor and soft ceilings deleted. Core 1,077, Contributions 1,276, Properties 4,511, Equilibrium 9,542, parameter Regression 3,211, Surface 577, Bindings 2,283, Python 2,025, total 24,499.
- 2026-10-10, #391 (owner): Equilibrium **8,843 → 8,917**, the measured count with the explicit cooling-onset dew operation (+95 Equilibrium, +1 Bindings, +2 Python). Round 2 (owner, small bump): **8,917 → 8,921** for the vapor-side reference admission (+4 Equilibrium). Total, floor and soft targets are unchanged.
- 2026-10-09, #380 (owner): parameter Regression **2,550 → 2,815**, the measured count of the exec-started parallel fit workers after the round-3 parent-death guard. It covers the round-1 approval (up to 2,664) and the round-2/3 approval (up to 2,825). Floor and soft targets are unchanged.
- 2026-10-07, #341 owner decision: the Python layer and the total no longer count the 1,225-line `_core.pyi` type stub; the architecture check compares its names and class members with the extension instead. The owner then lowered both hard ceilings by those 1,225 lines (`_core.pyi` at base `5cc1001d`): Python **2,892 → 1,667** and total **25,000 → 23,775**. Floor and soft targets are unchanged.
- 2026-10-06, #224 slice B (owner): Equilibrium 8,800 → 8,843, Python 2,850 → 2,892, the measured counts.
- 2026-10-06, #245 (owner): parameter Regression 2,478 → 2,550.
- 2026-10-06, #283 (owner): Core 917 → 1,100, Contributions 1,547 → 1,900, Properties 4,800 → 5,300, total 24,000 → 25,000; soft Core 850 → 1,000, Contributions 1,300 → 1,700, Properties 2,900 → 5,000, Equilibrium 6,000 → 8,000, total 14,850 → 19,500 (the layer soft sum).
- 2026-10-02, #220 envelope raise (owner): Properties 4,600 → 4,800, Equilibrium 7,575 → 8,800, Bindings 1,940 → 2,100, Python 2,700 → 2,850, total 22,000 → 24,000.
- 2026-10-01, #221 finding 5 closeout: Equilibrium 7,650 → 7,575.
- 2026-10-01, #223 (owner): Bindings 1,880 → 1,940, Python 2,620 → 2,700.
- 2026-10-01, #219 slice B closeout: parameter Regression 2,480 → 2,478.
- 2026-10-01, #219 slice B (owner): parameter Regression 2,200 → 2,480.
- 2026-10-01, #219 slice A closeout: Core 1,000 → 917, Contributions 1,612 → 1,547.
- 2026-10-01, #220 roadmap envelope (owner): Properties 4,258 → 4,600, Equilibrium 7,000 → 7,650, parameter Regression 1,619 → 2,200, Bindings 1,726 → 1,880, Python 2,374 → 2,620, total 19,100 → 22,000.
- 2026-09-30, #210 (owner): Properties 4,223 → 4,258, parameter Regression 1,484 → 1,619, Bindings 1,725 → 1,726, Python 2,357 → 2,374, total 18,830 → 19,100.
- 2026-09-30, #211 (owner): parameter Regression 1,463 → 1,484, Python 2,355 → 2,357, total 18,805 → 18,830.
- 2026-09-30, #205 (owner): parameter Regression 1,440 → 1,463, Python 2,350 → 2,355, total 18,790 → 18,805.
- 2026-09-30, #187 (owner): parameter Regression 1,300 → 1,440, Bindings 1,708 → 1,725, Python 2,313 → 2,350, total 18,556 → 18,790.
- 2026-09-29, #129/#193 (owner): Properties 4,211 → 4,223, parameter Regression 1,200 → 1,300, total 18,365 → 18,556.
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

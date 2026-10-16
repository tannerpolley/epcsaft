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
   the count, readiness fails with the ceiling to lower to (count plus band;
   for the total, no more than the sum of the layer ceilings).
   Lower it in the same change; lowering needs no approval. A ceiling up to one
   band above the count needs no table change.
3. **Approve once per issue.** The owner approves a line allowance per layer
   when authorizing an issue, and the issue body records it. A merge whose
   raises stay within that allowance needs no second approval. Growth beyond
   it needs a new owner decision with the scientific requirement and
   alternatives.

   Delegation (owner, 2026-10-16, #577): for the issues selected from the
   2026-10-15 audits (#563–#579), the lead approves each issue's allowance once
   its independent readiness review returns Ready, using the reviewer's
   recommended allowance, and records it in the issue body. Their combined net
   growth may not exceed +2,000 formatted runtime lines over the 26,886 counted
   at `517a3450`. Beyond that total, and for any build of the scope-change
   issues #573–#576, the owner decides.

A temporary over-budget working state may be simplified and repaired. Do not
freeze uncompilable code because an intermediate write exceeded a number. The
completed candidate must satisfy both band limits before acceptance.

Never satisfy a budget by dropping required physics, domain checks, derivative
terms, failure diagnostics or numerical falsifiers. If the smallest correct
method does not fit, report that evidence and revise the design/budget openly.
Budget changes must not silently change scientific scope.

## Amendments

- 2026-10-16, #566 (lead allowance under the #577 rule: Properties +243/318, Bindings +16/31, Python +14/24 net/gross, other layers 0): request-only derived State properties, partial molar properties and Henry constants use net +137/+11/+12 and gross 142/23/16 against `16fccc32`, re-measured after merging main `5dca2229`. Properties **6,004 → 6,186** and Bindings **2,178 → 2,205**, the counts 6,124 and 2,183 plus their 1% bands; Python 2,192 stays within its ceiling. Total **27,364 → 27,681**, the sum of the layer ceilings, below count 27,418 plus band 275. The response reuses the existing T/V Hessian columns; Henry constants are a thin wrapper over the unchanged reaction-reference owner.
- 2026-10-16, #571 (lead allowance under the #577 rule: Equilibrium +45/50, Bindings +3/3, Python +4/4 net/gross, other layers 0; the lead approved a 2-line in-place edit of the parameter Regression fit worker serializer, net 0): the direct azeotrope cell and predicted-anchor envelope steps use net +36/0/+2/0 and gross 49/2/2/2 (Equilibrium/Bindings/Python/parameter Regression) against `fd14723c`. Equilibrium **9,688 → 9,796**, count 9,699 plus its band 97; the total ceiling stays 27,364 above count 27,258.
- 2026-10-16, #567 (lead allowance under the #577 rule: Core +50/65, Contributions +100/150, Properties +30/45, parameter Regression +15/20, Surface +100/135, Bindings +5/10, Python +65/100 net/gross): temperature laws use net +23/+100/+30/0/+95/+2/+57 and gross 23/114/32/9/105/8/73 against `16fccc32`. Core **1,327 → 1,359**, Contributions **1,276 → 1,389**, Surface **577 → 673**, Python **2,143 → 2,200**, each count plus its 1% band; total **27,066 → 27,364**, the layer-ceiling sum, below count 27,220 plus band 273. Properties, parameter Regression and Bindings stay within their bands.
- 2026-10-15, #556 (lead allowance Properties +181 net / 240 gross, Bindings
  +35 / 45, Python +65 / 90 against `5d69b1d9`): measured growth is
  +181 / +13 / +32 formatted runtime lines. Properties **5,810 → 6,004**,
  Bindings **2,151 → 2,178**, Python **2,107 → 2,143**, each the measured
  count plus its 1% band. Total **26,809 → 27,066**, the layer-ceiling sum after #557.
  One all-species solvent contraction replaces both one-salt contractions;
  the one-salt declaration supplies its particle ratio without a rounding
  round trip through mole fractions. The lead raised the Properties net
  allowance from +180 to +181 for the two delivered-review failure-reporting
  corrections; gross and all other allowances, including total +280 / 375,
  are unchanged. The independent −2 net Parameter Regression change in #557
  is excluded from #556 growth.

- 2026-10-15, #557 (rule 2): Parameter Regression **3,874 → 3,873**, the
  formatted count 3,834 plus its band 39 after replacing the blanket overlap
  refusal with ordered coefficient composition (22 gross additions, −2 net).
  Total **26,810 → 26,809**, the sum of the layer ceilings. Other runtime
  layers are unchanged.
- 2026-10-15, #545 (owner): approved Equilibrium allowance +100 net / 160
  gross formatted lines relative to `0bec2bdf42e31bf37e8d4679169c8174b580408f`,
  with zero growth in other runtime layers. The measured change is +67 net
  formatted lines (9,596 → 9,663); Equilibrium's 9,688 and total's 26,810
  ceilings remain within their existing bands and need no raise.
- 2026-10-15, #543, tape optimization: one added runtime line in Properties
  against `0bec2bdf`; existing ceilings remain in force. Zeroth-order reuse
  retains its existing owner and implementation.

- 2026-10-15, #541, zeroth-order reuse: the selected prototype adds 16 and
  removes 6 runtime lines in Properties against `7029597d`. Tape-local reuse
  replaces the three unconditional order-zero sweeps. Existing ceilings remain
  in force; the separate tape-optimization lever remains unqualified.

- 2026-10-15, #513 (lead decision, R1 corrected after merging #497): **Properties allowance raised by the lead from +200 to +220 within the unchanged +238 total**; gross additions remain at most 300 Properties and 356 total. Against merged main `54a6d798`, the authored change adds Properties +220 and total +238. Main's Properties count is 5,532 after #497; the corrected combined count is 5,752. Properties **5,535 → 5,810**, count plus band 58; total **26,473 → 26,810**, the layer-ceiling sum, below combined count 26,584 plus band 266. The fixed design base remains `5ad73fb9`; separately accepted #505 and #497 additions are excluded from the authored allowance. Both main amendments below are retained. Shared sampled-contradiction detection, local mesh subdivision, partition agreement and cause preservation reuse the existing scan, fold bisection and Newton owners. R1 preserves a refined boundary crossing unless the ordinary partition owns its exact-zero node; it introduces no tolerance. Other layer allowances remain those in the accepted design.

- 2026-10-15, #497 (lead allowance: parameter Regression +60 net / 100 gross,
  Properties +5 net / 10 gross, other runtime layers zero): parameter Regression
  **3,812 → 3,874**, count 3,835 plus its band 39, for opt-in process CPU timing
  of native fit evaluations. Properties gains two trace observations within its
  existing ceiling; the total ceiling stayed 26,083 on the fixed base
  `5ad73fb9`. Integration with #505 retains its later total ceiling.

- 2026-10-15, #505: the accepted design's net allowances are Core +25,
  Equilibrium +380, Properties +10 (entry scopes), other layers zero;
  total +415, gross additions at most 510. Relative to `5ad73fb9`, the
  candidate adds 22 Core, 263 Equilibrium and 2 Properties formatted lines,
  287 total (including the rejected-result tuple reset). Core rises from 1,304 to 1,327, Equilibrium from 9,422 to
  9,688, and total from 26,083 to 26,473. The layer ceilings are the count plus the 1% band; the
  total is limited to the sum of layer ceilings.
  Final delivered review and adoption remain with the lead.

- 2026-10-15, #495 (lead-set allowances under the owner's overnight authority, 2026-10-15): Bindings **2,050 → 2,151**, the count 2,129 plus its band 22 after the nanobind port and delivered-review input corrections, including #486 solver options. The final count is 2,141, within one band of the ceiling. The issue allowance is +115 net / 1,100 gross Bindings lines, raised by the lead from +100 to restore the pybind11 input conversions the delivered review found missing; other runtime layers have zero additions. The total ceiling remains 26,083.

- 2026-10-15, #512: lead-set allowances under the owner's overnight authority, 2026-10-15; Properties +50 net / +140 gross formatted lines, every other runtime layer zero. The bracket change adds 1 net formatted line: Properties 5,517 (band 56), total 25,905 (band 260). No ceiling raise is needed: the existing Properties ceiling 5,535 and total ceiling 26,083 remain within their bands.

- 2026-10-15, #486: lead-set allowances under the owner's overnight authority, 2026-10-15; Core 1,229 → 1,304, Equilibrium 9,316 → 9,422, Bindings 1,923 → 2,050, Python 2,081 → 2,107, each measured count plus its 1% band; total 25,749 → 26,083 (sum of layer ceilings).

- 2026-10-15, #510 (owner net allowances): Properties **5,395 → 5,535**, parameter Regression **3,279 → 3,812**, Python **2,019 → 2,081**, each changed ceiling the formatted count plus its 1% band. Total **25,014 → 25,749**, the sum of the layer ceilings. Current counts are 5,480, 3,776 and 2,060: growth +139/+507/+61 against the fixed `ae6a8f11` base, within the owner allowances +260/+850/+70. Bindings grows +15 within its existing ceiling and +100 allowance. Equilibrium grows +10 from the separately merged #504 main changes; #510 changes no Equilibrium code. Analysis scripts have the separate owner allowance of 700 net lines.

One line per ceiling change, newest first (hard ceilings). The issue or PR
holds the requirement, alternatives and measured accounting; Git keeps the
former per-slice tables and the retired floor and soft ceilings. Entries before
#436 count physical lines. Entries marked **inherited** come from
[the `a0adb1cd` architecture comments](https://github.com/tannerpolley/ePC-SAFT/blob/a0adb1cd5ff5f4ee32b3608c78a02c756b06ebd9/ARCHITECTURE.yaml#L31);
their approval assertions were not reverified in the #221 reconciliation.

- 2026-10-14, #478 (owner allowances, post-#485 base `f64bf08f`): Core **1,077 → 1,229**, Properties **4,672 → 5,395**, counts 1,216 and 5,341 plus their 1% bands. Corrected Taylor density closure and conditional structured search use net +141/+12/+694 and gross +157/+124/+767 in Core/Contributions/Properties, within the authorized 145/20/725 net and 165/135/900 gross. Total **24,139 → 25,014**, the layer-ceiling sum, below count plus band. Other runtime layers are unchanged.

- 2026-10-13, #484 (owner confirmed 2026-10-13; merged in PR #490): parameter Regression **3,199 → 3,279**, Bindings **1,921 → 1,923**, Python **1,989 → 2,019**, the measured counts 3,246, 1,903 and 1,999 plus their 1% bands after the fixed-pivot reaction enthalpy direction, compact fitted-law publication and native point adapter (+74/+5/+33 gross; +64/+2/+30 net). Total **24,027 → 24,139**, the layer-ceiling sum, below count 23,922 plus band 240. The lead permitted development under +140/+12/+35 gross.
- 2026-10-13, #435 (rule 2, consolidation; the field-list macro accepted 2026-10-13 as the recommended option while the owner was away, confirmed at merge): Bindings **2,505 → 1,921** and Python **2,050 → 1,989**, the counts 1,901 and 1,969 on main `543dc58c` plus their bands, after `EPCSAFT_READONLY`/`EPCSAFT_READWRITE`/`EPCSAFT_VALUES` (`bindings/fields.hpp`) replaced one line per same-named field or enum value, one aggregate constructor and a typed continuation argument replaced binding glue, and five declaration classes became slotted dataclasses (−579 Bindings, −73 Python, −6 parameter Regression, which stays under its ceiling). Total **24,672 → 24,027**, the sum of the layer ceilings, below the count 23,826 plus its band.
- 2026-10-13, #434 (rule 2, consolidation): Equilibrium **9,538 → 9,316** and Bindings **2,506 → 2,505**, the counts 9,223 and 2,480 plus their bands after the numerically identical audit set A2–A8, A10 and A15 (−253 Equilibrium, −2 Bindings on main `9b39db7f`); the review then restored 8 lines of pre-refusal diagnostics (Equilibrium 9,231, within the band). Total **24,895 → 24,672**, the sum of the layer ceilings, below the count 24,464 plus its band.
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

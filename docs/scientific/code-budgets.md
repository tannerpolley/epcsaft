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
counter's `formatted_lines` holds the exact commands. Format original source
before excluding documentation or comments; comments are not reflowed.
For Python, mask only first-statement module, class, sync-function and
async-function string docstrings using their exact AST UTF-8 byte spans, then
count the union of lines occupied by retained tokens. Comments, whitespace,
layout and blank lines do not count. Other string expressions, assigned strings
(including blank-looking multiline data), annotations and declarations count;
code sharing a docstring line remains counted. Attribute prose belongs in its
class's NumPy Attributes section.
For C++, scan literals and preprocessing tokens after backslash-newline
continuation: exclude comments outside literals and whitespace-only lines.
Preserve ordinary/character/raw strings, raw-string data, header names and digit
separators. Declarations, braces, directives, macros and inactive branches count;
no preprocessing or compilation changes the count. Unterminated or unrecognised
lexical input fails visibly. This is formatted source size, not executed
instructions. The original-source formatter-disable refusal still runs first.
Count the
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

- 2026-10-29, #886 integration with `94b4de381`: retain both amendment histories. Bindings **2,420 → 2,457** combines the route’s existing ceiling with main’s 37-line #906 increase; the measured combined count is 2,433 (band 25). Main’s other layer ceilings and total 31,163 suffice for the integrated count 31,078. This combines the existing approved allowances; no new runtime capability or numerical rule is added. The #886 allowance remains 750 gross / +610 net.

- 2026-10-28, #886 after merging `fb27d881f`: total **30,372 → 30,394**, the measured formatted runtime count. Main contributes 30,177 lines and the saturation work contributes 217 net Engine runtime lines. Main’s Contributions, Properties, Surface and Python ceilings and the saturation work’s Equilibrium and Bindings ceilings remain sufficient. This increase uses the existing owner allowance of at most 750 gross / +610 net, including campaign support. The final count is 666 gross / 75 removed / +591 net; 34 added lines prepare the stricter campaign by reusing the existing runner. No numerical rule changes.

- 2026-10-28, dedicated pure-fluid route after integration with `f0e0de987`: total **30,253 → 30,276**, exactly the measured formatted runtime count. Main alone counts 30,059; the accepted route adds 217 Engine runtime lines. All layer ceilings pass; the old total alone was short by 23 after retaining main’s regression changes. The complete task count, including campaign support, is 575 gross / +500 net against current main, within the owner’s 750 gross / +610 net allowance. Shared cold enumeration remains main’s single admission owner and receives no new credit here.

- 2026-10-27, #886 (owner allowance 750 gross / +610 net against `7929f13be`):
  Equilibrium **11,950 → 12,083** is count 11,963 plus band 120. The increase
  retains observations from direct pure-state evaluation, one rejected native
  result and original-start continuation while both entries share the reduction.
  Bindings **2,325 → 2,420** is count 2,396 plus band 24, for one bulk snapshot
  with optional records, actual phase pressures and the owning native result.
  Total **30,025 → 30,253** is the layer-ceiling sum; total count is 30,035.
  Equations, precision, physical admission and general continuation remain
  their existing owners. Campaign code is included in the issue accounting.
- 2026-10-28, #906 (owner decision: 1,000 gross / +900 net formatted runtime lines; Equilibrium 500, Parameter Regression 360, Bindings 60, Python 80): the measured use against `f4566b3d0` is 756 added / 72 removed / +684 net. The affine loading coordinate, seeded grid, composed corrections, native comparison and capped-enthalpy composition reuse the existing driver and implicit actions. Equilibrium **11,950 → 12,427** (count 12,303, band 124); Regression **3,955 → 4,207** (count 4,165, band 42); Bindings **2,325 → 2,362** (count 2,338, band 24); Python **2,287 → 2,312** (count 2,289, band 23); Total **30,372 → 31,163** (count 30,861, band 309). Layer ceilings use count plus band; the total is bounded by their sum. Numerical evidence and change accounting are retained in `build/issue-906/build-r2/`.

- 2026-10-28, #896 (lead allowance: 200 gross formatted runtime lines): the declared solvent mole-weighting relation, membership/domain admission and fixed-parameter transport use 124 gross additions against `3cf744eed`. Contributions 1,273 → 1,289; Properties 6,346 → 6,428; Surface 785 → 815; Python 2,246 → 2,287; total 30,025 → 30,372. Area ceilings use their measured count plus the 1% band; the total is capped by the layer-ceiling sum. Other area ceilings suffice. Numerical qualification remains separate.

- 2026-10-27, #863 (lead allowance after Ready): parameter Regression **3,777 → 3,955**, the measured 3,915 plus band 40, for existing rank/whitening/hold decisions and selected first causes through predictions and worker replies, plus minimax's shared Ipopt observations. The approved task allowance was amended to 270 gross / +170 net. Actual runtime growth is 255 gross / +168 net against `7929f13be`; no numerical method or calculation budget changes. Other layer and total ceilings remain sufficient.

- 2026-10-28, #865 (owner decision: merge as a numerically neutral cleanup, with the speed result unresolved): the interim prepared-call package adds +40 net runtime lines. Total **30,025 → 30,059** equals the measured count after integration with main `3cf744eed`. Every layer stays within its existing ceiling.
- 2026-10-27, #862 owner scope change (rules 1 and 2): the allowance remains 1,100 gross / +860 net against `bd81c3ea4`; measured use is 957 added / 299 removed / +658 net formatted runtime lines. Diagnostic allocation-failure containment, incomplete-record fallbacks and controls are removed. Ordinary exception transfers copy the first record, preserving the anchor partition and association consumers' sources. Core **1,382 → 1,323** is count 1,309 plus band 14; Bindings **2,327 → 2,325** is count 2,301 plus band 24. Total **30,086 → 30,025** is the lowered layer-ceiling sum. Other counts/ceilings are Contributions 1,260/1,273, Properties 6,339/6,346, Equilibrium 11,841/11,950, parameter Regression 3,767/3,777, Surface 777/785 and Python 2,223/2,246. Total count is 29,817. Every ceiling remains within one band of its count. Shared first causes, named ranks, chronological recovery and bounded Ipopt observations retain their existing owners. Round-7 evidence is retained in `build/issue-862/build-r7/`.

- 2026-10-26, #855 (lead allowance: 160 gross formatted Equilibrium lines): lazy exact trial Hessians, identical first-point seed-root reuse and removal of the repeated initial Liquid compilation use 17 added / 11 removed formatted runtime lines against `761f14f43`, +6 net, all Equilibrium. Equilibrium counts 10,944 and total 28,028; both remain within the existing rule-1 ceilings and bands. Every other runtime layer and every ceiling is unchanged. Numerical and cost evidence is retained under `build/issue-855/build-r1/`.

- 2026-10-26, #869 (lead allowance, rule 1): Properties **6,129 → 6,346** and Equilibrium **11,399 → 11,530**, the final formatted counts 6,283 and 11,415 plus their 1% bands. Total **29,049 → 29,445**, the smaller of the layer-ceiling sum and total count 29,153 plus its 292-line band. Gross additions stay within the issue's 640-line allowance and its per-area limits; retained accounting is in `build/issue-869/build-r1/runtime-lines.json`.

- 2026-10-26, #853 round 4 (lead-approved Regression allowance: 44 gross formatted additions; Equilibrium 10, Bindings 0): pre-exec descriptor placement and closure bring parameter Regression to 3,739 lines. Rule 1 raises its ceiling **3,728 → 3,777**, count plus band 38. Other layer ceilings and the total remain unchanged. Qualification is retained in `build/issue-853/build-r4/`.

- 2026-10-26, #861 round 4, r5 scope and unchanged r4 allowance (985 gross: Core 230, Properties 190, Equilibrium 455, Bindings 70, Python 40), base `761f14f43`: final gross additions are 200/125/407/64/22 (818 total), 65 removals and +753 net. Source-reference error propagation moves to #869; dependent estimates report explicit unavailability. Frozen endpoint operands and reduced-radius projection reuse the native estimate/comparison owners. Rule-1 final ceilings are Core **1,264**, Properties **6,129**, Equilibrium **11,399**, Bindings **2,244** and Python **2,246**, each final count plus its 1% band. Relative to round 3, Core rises from 1,263, Equilibrium from 11,387, Bindings from 2,216 and Python from 2,244; Properties is unchanged. Total **29,006 → 29,049**, the layer-ceiling sum, below count 28,775 plus band 288. Unchanged layer ceilings remain within their bands.

- 2026-10-26, #861 round 3 preparation, r4 allowance (985 gross: Core 230, Properties 190, Equilibrium 455, Bindings 70, Python 40): after merging `fefb2f1ed` and main `761f14f43`, the runtime tree is unchanged from round 2. Recounting against the new base gives 776 gross additions and 65 removals, **+711 net**, correcting round 2's reported +722 net. The per-layer additions remain Core 199, Properties 125, Equilibrium 395, Bindings 37 and Python 20. Counts and ceilings remain unchanged at total 28,733 / 29,006. Source-reference error production remains incomplete; a forecast above the remaining Properties allowance stops further additions. This is a partial source snapshot, not final numerical qualification.

- 2026-10-26, #861 round 2, r3 allowance (985 gross: Core 280, Properties 180, Equilibrium 415, Bindings 70, Python 40): the partial candidate adds Core 199, Properties 125, Equilibrium 395, Bindings 37 and Python 20 gross formatted runtime lines against `5108b55a1`, 776 total and +722 net. Rule 1 sets Core **1,058 → 1,263**, Properties **6,003 → 6,129**, Equilibrium **11,034 → 11,387**, Bindings **2,179 → 2,216**, Python **2,201 → 2,244**, and total **28,242 → 29,006**, limited by the layer-ceiling sum below count 28,733 plus band 288. Source-reference numerical-error transport remains incomplete; these measured ceilings do not establish scientific acceptance.

- 2026-10-26, #854 round 2 (lead allowance: 950 gross formatted runtime lines, Core 45, Properties 200, Equilibrium 655, Bindings 30, Python 20): measured candidate use against unchanged main 37c58923e is Core 33, Properties 200, Equilibrium 601, Bindings 8 and Python 0 gross additions, 842 total and +536 net. Rule 1 sets Properties **5,893 → 6,003** (count 5,943 plus band 60), Equilibrium **10,861 → 11,034** (count 10,924 plus band 110), and total **27,959 → 28,242**, capped by the layer-ceiling sum below count 28,004 plus band 281. Round 3 retains these ceilings: current counts are Core 1,051, Properties 5,943, Equilibrium 10,927 and total 28,011, all within their bands. Gross additions against the same base are Core 37, Properties 200, Equilibrium 604, Bindings 8 and Python 0, 849 total and +543 net. The named two-density error policy and mathematical unit names add seven formatted lines over round 2; Properties remains at its 200-line gross allowance. Round 4 applies the accepted double-polish stopping rule: Equilibrium counts 10,936 and total 28,020, within the retained ceilings and bands. Gross additions are Core 37, Properties 200, Equilibrium 613, Bindings 8 and Python 0, 858 total and +552 net. Properties is unchanged. The existing 109-row ancillary check passes without a tolerance change; full numerical and cost qualification are retained in build/issue-854/build-r4/. Round 5 preserves the shared work-limit cause when bounds preparation stops: Equilibrium counts 10,938 and total 28,022, within the retained ceilings and bands. Gross additions are Core 37, Properties 200, Equilibrium 615, Bindings 8 and Python 0, 860 total and +554 net. The two added formatted lines change failure reporting only; the regular states remain bitwise identical to round 4. Round-5 evidence is retained in build/issue-854/build-r5/.


- 2026-10-24, #824 (owner allowance Properties +22 / Bindings +6; rule 1): the State density-branch fields `requested_density_branch`, `returned_density_branch` and `density_branch_mismatch` use Properties +12 and Bindings +4 formatted runtime lines, with zero in every other layer. Bindings **2,155 → 2,179**: count 2,157 plus band 22. Properties 5,801 stays within its ceiling of 5,826. The total count of 27,472 stays under the unchanged ceiling of **27,614**.
- 2026-10-24, #823 (owner allowance gross +515/net +275 across Properties and Contributions): the analytic solvent insertion reuses the scalar Born, association and pressure-Taylor owners. Contributions **1,235 → 1,254**, count 1,241 plus band 13. Rule-2 reductions: Core **1,107 → 1,058**, Properties **5,981 → 5,826**, parameter Regression **3,745 → 3,728**, after removing the finite-reference calculation and its seed/diagnostic consumers. Core and parameter Regression equal count plus band; Properties stays within one band of its count. Total **27,816 → 27,614**, capped by the layer-ceiling sum. Other ceilings stay fixed. Final task accounting is retained under `build/issue-823/` against `4fd40fe0f`.


- 2026-10-23, #814 (rule 2): removing the one-anion count requirement adds
  zero runtime lines and removes two from Properties against `e4d87ac8`.
  Lower Properties **5,983 → 5,981** (count 5,921 plus band 60) and total
  **27,818 → 27,816** (the layer-ceiling sum after merging #811). The equation is unchanged.

- 2026-10-23, #811 (owner approval, 2026-10-23): rule 1 records Core
  **1,100 → 1,107** (measured 1,096 plus band 11), Equilibrium
  **10,596 → 10,607** (10,501 plus band 106), and total **27,800 → 27,818**,
  the layer-ceiling sum below count 27,591 plus band 276. Against `e4d87ac8`,
  the approved allowances are Core +5 net / 6 gross and Equilibrium +160 / 280;
  every other runtime layer has zero allowance. The trial polish log-step limit
  uses Core +5 net / 6 gross and Equilibrium +10 net / 11 gross formatted
  lines. The numerical policy
  remains in Core; the existing trial polish remains its sole consumer.

- 2026-10-22, #791 (owner approval, 2026-10-22): Equilibrium **10,450 → 10,596**,
  the measured count 10,491 plus band 105, and total **27,654 → 27,800**, the
  layer-ceiling sum, below count 27,576 plus band 276. The approved task allowance
  is Equilibrium 100/140, Bindings 20/30 and Python 16/24 net/gross, zero elsewhere,
  against integrated main `955ff6d0`. Inherited #793 additions are excluded.
  Fixed-volume molarity rows and the existing ideal activity reader share one
  concentration expression; balance and solver owners are unchanged.

- 2026-10-22, #793 (owner approval, 2026-10-22): Properties **5,842 → 5,983**
  and total **27,384 → 27,654**, the measured Properties count plus its 1% band
  (5,923 + 60) and the layer-ceiling sum, below total count plus band
  (27,495 + 275). The explicit reaction activity reader shares
  the finite-state activity and source-reference owners. The owner allowance is
  +240 net / 342 gross: Properties 120/180, Equilibrium 90/120, Bindings 24/32,
  Python 6/10; all other runtime layers zero. Final measured accounting is in
  the #793 pull request. Other ceilings remain within their bands.

- 2026-10-22, #776 (owner approval, 2026-10-22): the ideal mole-fraction amount
  chart, shared feasibility lifecycle and native effective-water map use +259 net /
  293 gross formatted Equilibrium, +24/30 Bindings, +30/47 Python and -2/1
  Properties against `83d6755d`; total +311/371. After merging main `fc247023`
  with #766 and the species-thermochemistry path, the combined counts are
  Equilibrium 10,213, Bindings 2,138, Python 2,179 and Properties 5,784.
  Rule 1 raises Equilibrium 10,054 to 10,316 (band 103), Bindings 2,130 to
  2,160 (band 22), Python 2,168 to 2,201 (band 22) and total 26,837 to
  27,384: min(27,115 plus band 272, layer-ceiling sum 27,384).
  Rule 2 lowers Properties 5,844 to 5,842 (5,784 plus band 58).
  The #766 amendment below remains as merged. This approval covers the allowance
  and measured ceiling amendments; it does not grant merge authority.

- 2026-10-22, #766 (owner approval, 2026-10-22): the allowance permits
  Core 8/10, Contributions 6/8, Properties 85/125, Equilibrium 145/215,
  parameter fitting 4/8, Surface 30/45, Bindings 14/21 and Python 12/20
  net formatted runtime lines / gross executable additions against `83d6755d`.
  Measured counts require Properties 5,768 → 5,844 (5,786 + band 58),
  Equilibrium 9,936 → 10,054 (9,954 + band 100), and Surface 759 → 785
  (777 + band 8), exactly rule 1. The total ceiling remains 26,837.
  The delivered change uses 223 net / 310 gross.

- 2026-10-22, #781 (lead allowance: Properties at most 2 gross additions, net growth at most 0; every other runtime layer 0): removing the neutral-component count uses 2 added / 3 removed / −1 net formatted runtime line against `2c52c61a`. Properties counts 5,710, so its ceiling lowers **5,769 → 5,768** under rule 2. Total counts 26,580 and its ceiling stays 26,837.
- 2026-10-22, #771 (owner approval, 2026-10-22): Equilibrium allowance
  +135 net / 196 gross and Python +6 / 14 against `83d6755d`; every other
  runtime layer has zero issue-owned growth. This includes the scoped
  normalized declared-root check in the existing charged-Vapor trial polish.
  After the admission-order correction, unused result-member removal and
  current main merges, the issue contributes +134 net / 193 gross Equilibrium
  lines and +6 / 14 Python lines. Equilibrium counts 10,347 against
  the 10,450 ceiling, within its 104-line band. Total counts 27,274, within main's
  unchanged ceiling 27,384. These ceilings record source size and owner
  authority; numerical qualification is established separately by the
  manufactured recovery and the remaining falsifiers.

- 2026-10-20, #700 / #695 (owner decision, 2026-10-20): exclude Python AST-span docstrings, comments and empty layout lines, and C++ comments outside literals and empty lines, while retaining the pinned formatters and all bypass refusals. This accounting change adds or removes no Engine calculation. All eight layer ceilings lower to the measured count plus the rounded-up 1% band under rule 2; the total lowers from 29,191 to 26,837 for 26,571 counted lines, capped by the layer-ceiling sum. Per-layer before/after measurements are retained in the issue's PR.

- 2026-10-20, #692 (owner approval): shared physical-mol recovery has an allowance of +120 net / 166 gross formatted Equilibrium lines against `31460744`, with zero growth in every other runtime layer. Equilibrium **10,582 → 10,794**, the count 10,687 plus its rounded-up 1% band of 107 under rule 1. Total count **29,172** remains under the unchanged **29,191** ceiling. One recovery beside `run_driver` owns eligibility, invocation, admission and work/cause retention across ordinary solves, discovery, start sets and continuation.

- 2026-10-19, #664 after merging main `6385f895` (#663): the full amount-weighted stability correction and stored-update/representation accounting use +40 net / 56 gross formatted Equilibrium lines against both the original `4bf2378a` base and the post-#663 baseline, within the lead allowance of 42 net / 60 gross. The merged count is 29,029; Equilibrium is 10,544 / 10,582. The total remains below main's 29,191 ceiling, so rule 1 needs no further raise. All layer ceilings remain unchanged; #663 growth is recorded separately below.

- 2026-10-19, #663 after merging main `5aa2d3cb` (#658/#661): the approved Equilibrium allowance is +60 net / 100 gross formatted runtime lines; the change uses +22 net / 40 gross. The merged count is 28,989, seven above the previous total ceiling of 28,982. Under rule 1 the total ceiling rises to 29,191: min(28,989 + rounded-up 1% band 290, layer-ceiling sum 29,191). Equilibrium is 10,504 / 10,582; all layer ceilings remain unchanged.

- 2026-10-19, #658: lead-approved 220 gross / +180 net allowance. Contributions rises from 1,424 to 1,451, Surface from 791 to 814 and Python from 2,421 to 2,473: each exceeded ceiling becomes its measured count plus the rounded-up 1% band. The measured counts are 1,436, 805 and 2,448; total is 28,962. The total ceiling stays 28,982. Other ceilings are unchanged.

- 2026-10-19, #652 (lead allowance under the owner's 2026-10-17 delegation: Equilibrium 90 gross / 80 net, plus parameter Regression 1 gross / 0 net for `PhaseSpec` worker serialization): after delivered-review R1 and merging main `45713fae0` (#648), the consumed declared trace-guess cap uses 54 gross / +47 net Equilibrium and 1 gross / 0 net parameter Regression. Equilibrium is 10,477; the lead approved **10,475 → 10,582**, exactly count plus its 105-line band, through the delegated worker's native budget question on 2026-10-19. Total 28,839 remains below 28,982; parameter Regression is 3,962 / 4,000. The one seeding owner excludes declaration guesses replaced by valid previous phases and preserves supplied and calculated provenance through native fit workers.
- 2026-10-19, #648 after merging main `96728cb3` (#647), under the same owner delegation: the combined runtime count is **28,792**, within the delegated 28,919 count cap. Keep #648’s Properties **6,175 → 6,241** and Core **1,375 → 1,390** ceilings and #647’s Equilibrium **10,367 → 10,475** ceiling. Total **28,782 → 28,982**, rule 1 capped by the sum of the layer ceilings: min(count 28,792 + band 288, ceiling sum 28,982). The #648 change remains +129 net against merged main.
- 2026-10-19, #648 (lead allowance under the owner's 2026-10-17 delegation, after independent Ready: at most +137 net runtime lines): the fixed-temperature boundary pressure start adds Core +7. Core **1,375 → 1,395**, count 1,381 plus band 14. The superancillary adds Properties +63 and Equilibrium +11; combined with the loop, Equilibrium is +59 and the total is +129 against main `50df83ee`. Properties **6,175 → 6,241**, count 6,179 plus band 62. After merging main `50df83ee` (#646), Core is 1,376 and its ceiling lowers to **1,390** (count plus band 14); total 28,692 stays below its 28,782 ceiling.

- 2026-10-19, #647 (lead allowance under the owner's 2026-10-17 delegation, after independent Ready: Equilibrium +95 net / 105 gross, Core +5 / 5, with rule 1 raises): the native Michelsen tangent-plane trial minimizer for all-neutral, non-reacting systems uses net +95/+5 and gross 104/5 (Equilibrium/Core) against main `50df83ee`. Equilibrium **10,367 → 10,475**, count 10,371 plus its band 104. Core 1,374 stays within its ceiling 1,375. The total ceiling stays main's **28,782**, above the count 28,663 by less than its band 287.
- 2026-10-18, #641 (lead allowance under the owner's 2026-10-17 delegation, after independent Ready: +28 net / 39 gross runtime lines; design request Core +2/3, Contributions +8/10, Properties +3/4, Equilibrium +1/2, parameter Regression 0/1, Surface +14/18, Python 0/1): the opt-in charged-molecule ion declaration uses net +1/+9/+3/+1/+1/+11/0 and gross 2/12/4/2/2/14/1 (Core/Contributions/Properties/Equilibrium/parameter Regression/Surface/Python), +26 net and 37 gross, against main `cc63a135`. Surface **780 → 791**, count 783 plus its band 8. After merging main `b72de903` (#605) the total ceiling stays main's **28,782**, above the count 28,557 by less than its band 286.
- 2026-10-18, #605 (owner allowance 2026-10-18: +30 net runtime lines, with rule 1 raises): one Ipopt linear-solver owner (MA57 with its own scaling when the user's HSL library loads, else MUMPS) and its read-back use net +27/+1/+1 and gross 31/3/7 (Equilibrium/parameter Regression/Bindings) against main `cc63a135`. Equilibrium **10,249 → 10,367**, count 10,264 plus its band 103. Total **28,514 → 28,782**, the sum of the layer ceilings, below count 28,531 plus band 286.
- 2026-10-17, #594 (lead allowance after independent Ready: Properties about −20 net ±5, Core 0 to +1, Bindings 0 to −1; rule 1 lowering): one terminal-pair limit owner for the ion reference, its parameter and T/P columns and the reaction-reference transfer uses net −18/0/0 and gross 92/5/1 (Properties/Core/Bindings) against main `4d9b7bdb`. Properties **6,186 → 6,175**, count 6,113 plus its band 62. The total ceiling stays **28,514**, above the count 28,432 by less than its band 285.
- 2026-10-17, #608 (lead allowance under the owner's 2026-10-17 delegation, after independent Ready: Python +110 net / 120 gross, other layers 0): `parameter_document` and `Parameters.from_packets` build one parameter document from several packets with set k_ij/l_ij, using net +107 and gross 108 Python against main `ee6c5cd8`. Python **2,303 → 2,421**, count 2,397 plus its band 24. The total ceiling stays **28,514**, above the count 28,426 by less than its band 285.
- 2026-10-16, #575 (owner allowance: parameter Regression +50/70 with its ceiling raised to the final count plus band, Bindings +3/5, Python +12/20 net/gross, total +65): the mixed-composition osmotic observation uses net +42/0/+10 and gross 44/4/10 (parameter Regression/Bindings/Python), re-measured against main `87c48fa1`. Parameter Regression **3,957 → 4,000**, count 3,960 plus its band 40. After merging main `a5e9ecbc` (#595) no layer exceeds its ceiling, and the total ceiling stays main's **28,514**, above the count 28,319 by less than its band 284.
- 2026-10-16, #595 (owner allowance: Equilibrium +30 net / 45 gross, Properties +8 / 12, other layers 0): reacting equilibrium Cp actions request only the reference columns they read, and one ladder per solved state serves them, using net +19/+5 and gross 45/12 (Equilibrium/Properties) against `56216880`, unchanged after merging main `87c48fa1`. Equilibrium 10,166 and Properties 6,131 stay within their ceilings 10,249 and 6,186 and their bands; total 28,267 of 28,514. `reaction_references` and the model's reference store remain the one ladder owner.
- 2026-10-16, #573 (owner allowance: Equilibrium +360, Properties +20, Bindings +40, Python +40, cap +460 net): pure stoichiometric solids, with delivered review r1's complementarity classification and absent-solid curvature row, use net +269/0/+7/+30 (Equilibrium/Properties/Bindings/Python) and parameter Regression +1 (the fit-observation refusal, within its ceiling), +307 total, against main `56216880` (after #565), re-measured after merging main `30b84185` (#574). Equilibrium **9,977 → 10,249** and Python **2,271 → 2,303**, the counts 10,147 and 2,280 plus their 1% bands; Bindings 2,226 stays within its ceiling. Total **27,966 → 28,514**, the sum of the layer ceilings, below count 28,239 plus band 283.
- 2026-10-16, #574 (owner-approved allowance: Contributions +42, Surface +12, Core +4, Properties +3, Python +2, Bindings 0; total +63 net): the Gross (2005) quadrupole term uses net +34/+12/+2/+2/+2/0 and gross 64/21/3/4/2/0 (Contributions/Surface/Core/Properties/Python/Bindings) against `5dca2229`, plus a net-0 one-line in-place edit of the parameter Regression fit worker serializer for the two new parameter vectors. Measured after merging main `56216880`: Contributions **1,389 → 1,424**, count 1,409 plus its band 15; Surface **768 → 780**, count 772 plus its band 8; Core **1,359 → 1,375**, count 1,361 plus its band 14 (#565 had filled Core to its ceiling). The total ceiling stays main's **27,966**, above the count 27,932 by less than its band 280.
- 2026-10-16, #565 (lead allowance under the #577 rule: Equilibrium +270 net / 270 gross, Core +10 / 10, other
  layers 0, net over the merge base): measured growth against main `8fb04d01` (after #568, #569 and #570) is
  Equilibrium +154 net and Core +9 formatted runtime lines. Equilibrium **9,796 → 9,977** (count 9,878 plus band
  99). Core 1,359 sits exactly at its ceiling 1,359 and is not raised. The total ceiling stays main's **27,966**,
  above the count 27,880 by less than its band.
- 2026-10-16, #569 (lead allowance under the #577 rule after Ready: Surface +40, Bindings +10, Python +50, other runtime layers 0): measured growth against main `fe47081d` is Surface +39, Bindings +2, Python +50 formatted runtime lines (total 27,626 → 27,717). Surface **729 → 768** and Python **2,200 → 2,271**, each the count plus its 1% band; Bindings 2,219 stays within its ceiling 2,240. The total ceiling is **27,681 → 27,966**, the layer-ceiling sum. One native cross-association admission owner and one native NASA9 gas-constant rescale replace no existing runtime code.
- 2026-10-16, #570 (lead allowance under the #577 rule: Surface +55, Bindings +35, Python +8 net, other layers 0; Surface includes +1 for delivered review r2's law-term field guard): fitted-parameter write-back and form-tagged ASL thermochemistry records use net +55/+34/+6 against `16fccc32`, re-measured after merging main `35cbe62a`. Surface **673 → 729** and Bindings **2,205 → 2,240**, the counts 721 and 2,217 plus their 1% bands; Python 2,198 stays within its ceiling. The total ceiling stays main's **27,681**, above the count 27,626 by less than its band. One native inverse beside `build_parameters` writes fitted slots, temperature-law coefficients included, back into the start document, and one form-tagged interval writer and reader replace the NASA9-only record.
- 2026-10-16, #572 (lead allowance under the #577 delegation: parameter
  Regression +180 net / 200 gross formatted runtime lines, every other runtime
  layer zero, against `517a3450`): measured growth is 103 added, 20 removed,
  +83 net in `fit.cpp`/`fit.hpp`. Parameter Regression **3,873 → 3,957**, the
  count 3,917 plus its band 40. After merging `4302fc8f` the total
  ceiling stays main's **27,681**, above the count 27,531 by less than its band.
  The +83 counts toward the delegation's +2,000 cap. The joint measurement
  covariance replaces the blanket correlated-declaration refusal within the
  existing measurement chart and evaluator.
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

- 2026-10-26, #854 (lead allowance, rule 1): Properties **5,826 → 5,893**, Equilibrium **10,607 → 10,861**, the partial two-density candidate counts 5,834 and 10,753 plus their 1% bands. Gross formatted additions are 33 Properties, 430 Equilibrium and 2 Bindings (465 total, net +252), within the authorized 150/750/30 and 950-total allowances. Total **27,614 → 27,959**, the layer-ceiling sum, below total count plus band. Numerical qualification remains incomplete: the retained 375.13 K and 375.139 K comparisons exceed the 1e-10 equivalence limit.

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

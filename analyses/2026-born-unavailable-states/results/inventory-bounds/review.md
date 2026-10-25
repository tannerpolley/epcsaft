# Independent delivered review: #237

**Supported for the stated scientific capability.** GPT-6.1 Sol, xhigh, through
Pi independently reviewed the final functions, callers, retained numbers and
installed-wheel probes. Final integrated recheck: `bcbda83f` plus the two prose
updates against `9bc6e902`. No required #237 correction remains.

The first inspection found two defects: HiGHS dropped a tiny nonzero transformed
inventory coefficient, and ordinary boxes requested off-support objectives.
The final inspection verified the equivalent fraction-column scaling and
original-chart witness validation at $10^{-10}$ and $10^{-12}$ mol imbalances,
plus support-union box requests. Empty finite sets, zero inventory, genuine
suprema, structural zero, amount/fraction units, objective restoration and LP
failure diagnostics retain their distinct meanings. One shared HiGHS owner
replaces the embedded box implementation; no second bound solver remains.
The final caller inspection confirmed that main's chart retirement is preserved:
`global_scan.py` declares the copied problem and reads native bounds only. The
reviewer's helper execution reproduced all four retained maxima/logarithms
exactly; the three remaining upstream `linprog` entries are removed.

Born errors are $1.11\times10^{-16}$ for fractions, $6.22\times10^{-15}$ for
logarithms and $3.33\times10^{-16}$ for feasibility. Analytic values, explicit
floors, sevenfold feed scaling and the existing trace-inventory cases support
the stated inputs. Keep the extended analytic, native support and trace checks:
they protect distinct numerical risks. The reviewer confirmed wheel identity
`2ca3b44d…` across bound, analytic and final branch-timing records and verified
copied-input hashes.

The final reacting MEA best-of-three timing increase of **1.09% is flagged**.
Inputs, states, residuals and Ipopt work are identical; overlapping timing
ranges and changing signs across prior runs do not demonstrate a change-caused
slowdown. No speedup or tight sub-percent guarantee is supported.

Measured shipping runtime growth is **199 added / 21 removed / +178 net**,
21,523 total. Full executable accounting is 339 added / 67 removed, including
140 analysis additions and 46 removals: native-consumer conversion, bounded inventory replay and a
cost/provenance reader. Those readers reuse the scientific input owners and
native calculations; they hide no solver. The reviewer justified this growth
and found no remaining duplication. The sprint's 260-line figure is a target,
not a hard gate; no ceiling was amended.

The corrected-tree checks passed 115 Python tests, six fast native programs,
all three mapped slow programs (boundary, equilibrium and regression_fit), Ruff,
mypy and architecture. The queued HEAD check that returned no changed paths is
not numerical evidence; committed changes were checked against `4d1210c1^`.
The earlier IDAES duplicate was removed by main's #254 cleanup, present in
`9bc6e902`. The remaining local duplicate was identical main timing provenance:
`cost-pre-review-main.json` now points to `cost-main.json` with its SHA-256,
while the timing samples remain separate. Pointer hash verification and the
82-analysis layout check pass. Final `bash tools/check-all.sh` and
`bash tools/check-affected.sh HEAD` both passed through the capped slot wrapper
after this correction. The full gate passed 115 Python tests, six native programs,
Ruff, mypy, architecture and the clean 82-analysis layout. The earlier committed
native-change check also passed all three mapped slow programs. No delivery
blocker remains for this slice.

Qualification is numerical verification of conserved-inventory ranges, not
physical validation, thermodynamic equilibrium/stability, derivatives of maxima,
or exhaustive conditioning coverage. #234 batch 15 owns notebook consolidation.

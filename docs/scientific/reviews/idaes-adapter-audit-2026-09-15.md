# IDAES adapter audit - 2026-09-15

## Scope and inspected state

This audit reconciled the neutral IDAES/Pyomo ASL adapter sprint with the current
committed Properties owner. The repair branch was rebased onto
`codex/properties-qualification` at `6795d1e`; `origin/codex/properties-qualification`
was fetched and resolves to the same commit. The saved Properties checkout has
uncommitted P4 source, test, roadmap, and review changes. Those changes were not
captured, altered, or treated as a rebase target.

Local branch ownership is:

- `idaes/adapter/historical-prototype` preserves the original prototype at `c995e39`.
- `idaes/adapter/neutral-current-engine` preserves the pre-rebase sprint at `0194d4e`.
- `idaes/adapter/qualification-repair` contains the reviewed repair and retained evidence.
- `idaes/adapter/properties-integration` identifies the accepted Properties-based tip.

## Findings and disposition

1. Existing retained evidence could be deleted before Ipopt validation. The runner
   now requires a fresh output directory and refuses a nonempty one without changing it.
2. Failure and interruption summaries lost wheel, ASL, packet, source, and Ipopt
   identity. Identity is now written before numerical work and merged into every
   atomic checkpoint.
3. Independent reference rows were written only after the full reference stage.
   Each row is now retained before acceptance; nonfinite values and a false aggregate
   result are rejected.
4. Conditioning covered only 45 grid evaluations. `domain.csv` now owns those 45
   evaluations plus the nine derivative states, all with pressure, slope, RT-scaled
   slope, source, and pass status.
5. Whole-NLP validation followed four solves. The process NL files are now written
   first, `process-0.nl` is validated by `AslNLP`, and only then are the exact same
   retained files passed to Ipopt.
6. Native pressure-Hessian unit metadata labeled both gradient components as Pa.
   The temperature and volume entries now use `Pa/K` and `Pa m3/mol`.
7. The analysis reimplemented packet fingerprinting. It now validates and consumes
   the packet materialization record owned by the repository.
8. Solver row/column maps were duplicated. One verified canonical map is retained;
   each distinct NL remains retained and hash-linked from its process row.
9. The README mixed superseded sprint history with the current claim. It now states
   one current bounded qualification and its non-goals.

Independent scientific/native review found no remaining material issue. Independent
evidence review found the nonfinite-reference gate defect above; after correction,
its focused injection rejected NaN and infinities and the clean-source campaign
passed. `bash tools/check-all.sh` then passed packet validation, nine Python tests,
four native suites, Ruff, mypy, analysis layout, and architecture.

## Accepted result and limits

The retained campaign was produced from clean source commit `be200c6` using wheel
SHA-256 `5b47511b6d6863aef304aaef2d384b4cdde37bdbdd5f33720048f771592869c8`.
It passed 54 conditioning evaluations, seven independent references, 72 scalar
derivative ladders, 216 whole-NLP directional ladders, 27 Hessian matrices, four
optimal KKT cases, and forced-timeout cleanup. SIGTERM and SIGINT probes preserved
completed numerical prefixes, campaign identity, and terminated child groups.

This accepts only the issue #10 neutral propane adapter analysis. The issue-owned
use of `epcsaft._core` is a comparison exception, not a public B1/B2 derivative API.
P4, P5, association, mixtures, electrolytes, phase stability, PROMMIS, and downstream
application qualification remain open.

## Exact next execution order

1. Finish and review the uncommitted Properties P4 work in its saved checkout; do
   not fold it into the adapter branch as an incidental commit.
2. Commit the accepted P4 slice on `codex/properties-qualification`, then rebase
   `idaes/adapter/properties-integration` onto that new committed tip.
3. Rebuild and install one non-editable wheel, run `bash tools/check-all.sh`, and
   replay the bounded adapter campaign into a fresh candidate directory.
4. Compare the new wheel, ASL library, packet, NL, domain, reference, derivative,
   KKT, and failure identities before replacing retained evidence.
5. Complete the P4 observable-coordinate-order matrix and its numerical actions;
   keep P4 open until that broader owner is independently accepted.
6. Complete P5 current-wheel integration qualification across P1-P4.
7. Design and qualify the thin public B1/B2 action interface from accepted P4/P5
   contracts; do not expose analysis-only `_core` objects.
8. Only then adapt and qualify the real IDAES/PROMMIS consumers in their independent
   downstream repositories using an identified immutable Engine wheel.

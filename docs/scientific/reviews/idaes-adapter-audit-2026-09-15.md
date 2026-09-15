# IDAES adapter audit - 2026-09-15

## Scope and inspected state

This audit reconciled the neutral IDAES/Pyomo ASL adapter sprint with the completed
Properties owner. The repair branch was rebased onto
`codex/properties-qualification` at `67a7026`; local and origin resolve to the same
commit. The saved Properties checkout's unrelated user changes were preserved.

Local branch ownership is:

- `idaes/adapter/historical-prototype` preserves the original prototype at `c995e39`.
- `idaes/adapter/neutral-current-engine` preserves the pre-rebase sprint at `0194d4e`.
- `idaes/adapter/qualification-repair` contains the rebased repair and retained evidence.
- `idaes/adapter/properties-integration` is advanced to the final rebased tip during
  closeout.

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

The retained campaign was reproduced from clean source commit `48ee5cf` using wheel
SHA-256 `f8ba51165bb6773ae060c499d723ea3a66313bd6a1d10dd122423c89c06d4691`.
It passed 54 conditioning evaluations, seven independent references, 72 scalar
derivative ladders, 216 whole-NLP directional ladders, 27 Hessian matrices, four
optimal KKT cases, and forced-timeout cleanup. SIGTERM and SIGINT probes preserved
completed numerical prefixes, campaign identity, and terminated child groups.

The stable derivative, domain, whole-NLP and independent-reference tables are
identical to the pre-rebase accepted tables after line-ending normalization. The
adapter source SHA-256 is unchanged. The rebuilt wheel and ASL library identities
changed because the branch now inherits qualified P4/P5.

Commit `3460a53` retains the promoted clean-source campaign, corrected equation
status and interruption-check semantics. The final interruption replay ran from
that clean commit and records its source identity in both retained summaries.

## Claude Fable review and disposition

Claude Fable 5.1 reviewed the rebased implementation and evidence with no workspace
tools. It confirmed the adapter/Properties seam, pressure chain rule, CppAD reverse
indexing, mutex handling and unchanged Contributions. Its `REVISE` verdict identified
three record-integrity corrections:

1. Promoted evidence was not yet committed. Commit `3460a53` now retains it, and
   this follow-up records the clean interruption replay.
2. The `pressure_state_hessian` change note still called the implemented derivative
   slice missing. The canonical equation and generated index now state the bounded
   implemented and tested scope.
3. SIGTERM and SIGINT both reported return code `-2` without explaining that the
   SIGTERM handler deliberately raises `KeyboardInterrupt`. The checker now asserts
   and records that mapping; both paths passed again from clean commit `3460a53`.

This accepts only the issue #10 neutral propane adapter analysis. The issue-owned
use of `epcsaft._core` is a comparison exception, not a public B1/B2 derivative API.
P4 and P5 are inherited as qualified Properties work. Association support in this
adapter, mixtures, electrolytes, phase stability, PROMMIS, and downstream application
qualification remain open.

## Exact next execution order

1. Publish the rebased IDAES branch separately on top of `codex/properties-qualification`;
   do not merge adapter history into Contributions or Properties.
2. Keep the ASL bridge in native bindings and the exact pressure-curvature
   calculation in Properties; Contributions remains unchanged.
3. Design and qualify the thin public B1/B2 action interface from accepted P4/P5
   contracts; do not expose analysis-only `_core` objects.
4. Only then adapt and qualify the real IDAES/PROMMIS consumers in their independent
   downstream repositories using an identified immutable Engine wheel.

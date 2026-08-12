# CAP execution map

GitHub Issues in `tannerpolley/ePC-SAFT` is the sole tracker. The canonical
program is [root issue #10](https://github.com/tannerpolley/ePC-SAFT/issues/10);
native subissues and blocked-by relationships define the live frontier. Read the
selected issue and its comments with `gh` before implementation.

## CAP-09 promotion

Governance D-030 completes the D-029 protected transition. The user-approved
`1e-4` lineage-only comparison passes all eight Khudaida Figure 2 feeds while
the stricter physical and derivative gates remain unchanged. The accepted
receipt is
`validation/campaigns/2026-khudaida/results/receipts/cap-09-phase-equilibrium-promotion.json`.
`phase_equilibrium` is the sole public nonreactive route; omitted fixed-T,P
topology runs the retained internal HELD/HELD2 discovery algorithm and declared
topology never falls back. Issue #18 closes after the final checks, independent
reviews, clean push, and main synchronization complete.

## CAP-10 admission decision

Governance D-031 defers omitted-topology discovery, declared fixed-T,P phase
splits, and fixed-T coexistence boundaries from the Regression Observation
Family catalog. Regression remains EOS-only and no CAP-10A issue is created.
CAP-11 issue #20 is therefore the direct native frontier after issue #19.

## CAP-11 completion

Governance D-032 accepts the exact CAP-11 implementation tree and final
retirement receipt. Superseded live designs, the executable Cython experiment,
split-header packaging remnants, literal-component HELD coordinates, and dead
legacy switches are removed. The internal HELD/HELD2 numerical core remains
the canonical omitted-topology discovery implementation behind
`phase_equilibrium`. Issue #20 and root issue #10 close after the receipt commit
is synchronized. Issue #6 remains a separate post-cutover Validation campaign;
issues #30 and #31 are detached post-cutover reactive-observation follow-ons
and retain the dependency chain #20 -> #30 -> #31 without extending root #10.

## Private development feedback

Governance D-033 makes wheel iteration build-first and latency-bounded. Setup
builds and installs before tests. `tools/check-all.sh` is the parallel fast
installed-wheel gate with a warm target below ten seconds; it does not run slow
Engine stress cases, Validation campaigns or figure reproduction, full native
clang-tidy, or sanitizers. Those remain explicit owner commands selected by the
scientific claim or native risk, not automatic completion ritual.

## Task strategy

Run each unblocked CAP leaf in a fresh user-visible Codex task attached to this
project. Start from the accepted `main` commit, create one `codex/cap-*` branch,
and keep one accountable writer for the full Data–Engine–Validation slice. The
same task completes implementation, verification, review, commit, push, issue
closure, and synchronization of the canonical checkout. Use a worktree only for
explicitly independent parallel work. If a CAP task becomes oversized, take
over that leaf into one successor rather than carrying the whole program in one
conversation.

Accepted plans and ADRs carry design context; GitHub carries lifecycle state.
Task transcripts and handoff summaries are navigation aids only. This project
has no CI/CD.

## Model and review map

Every CAP primary task uses `gpt-5.6-sol` with `xhigh` reasoning. Pass both
settings explicitly during task takeover or creation; if the exact combination
is unavailable, stop instead of substituting it.

| CAP leaf | Primary task | Read-only review |
| --- | --- | --- |
| CAP-01 | `gpt-5.6-sol` / `xhigh` | one focused reviewer |
| CAP-02 | `gpt-5.6-sol` / `xhigh` | one focused reviewer |
| CAP-03 | `gpt-5.6-sol` / `xhigh` | two reviewers: derivative science and code surface |
| CAP-04 | `gpt-5.6-sol` / `xhigh` | two reviewers: Born science and parameter semantics |
| CAP-05 | `gpt-5.6-sol` / `xhigh` | one focused reviewer |
| CAP-06 | `gpt-5.6-sol` / `xhigh` | two reviewers: ownership boundary and regression numerics |
| CAP-07 | `gpt-5.6-sol` / `xhigh` | two reviewers: transport completeness and topology |
| CAP-08 | `gpt-5.6-sol` / `xhigh` | two reviewers: implicit derivatives and conditioning |
| CAP-09 | `gpt-5.6-sol` / `xhigh` | two reviewers: phase science and solver architecture |
| CAP-10 | `gpt-5.6-sol` / `xhigh` | two reviewers: admission evidence and scope control |
| CAP-11 | `gpt-5.6-sol` / `xhigh` | two reviewers: retired surface and final evidence |

Reviewers inspect only stable material and do not own Git or write files. Each
listed review is required; two-reviewer rows use independent passes over their
named concerns.

## Resolution gate

Apply the CAP issue resolution standard in the canonical hard-cutover plan.
Each changed scientific claim gets only its narrowest durable numerical
falsifier; zero new tests is correct when the issue adds no distinct claim.
Readiness checks, schemas, imports, paths, prose, and API shape never become
scientific tests.

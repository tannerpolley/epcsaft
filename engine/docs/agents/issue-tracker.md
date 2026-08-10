# CAP execution map

GitHub Issues in `tannerpolley/ePC-SAFT` is the sole tracker. The canonical
program is [root issue #10](https://github.com/tannerpolley/ePC-SAFT/issues/10);
native subissues and blocked-by relationships define the live frontier. Read the
selected issue and its comments with `gh` before implementation.

## CAP-09 correction

The closed state of issue #18 records the rejected scratch replacement and does
not authorize retirement of `tp_flash`. Governance D-029 and Equilibrium ADR
0002 protect the exact HELD/HELD2 incumbent until a HELD-derived
`phase_equilibrium` candidate completes actual constrained one-to-three-phase
solving and passes the reviewed as-good-or-better promotion gate. Do not infer
current sole-owner capability from #18 or advance a dependent phase capability
on that basis.

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

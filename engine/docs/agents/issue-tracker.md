# Engine issue routing

GitHub Issues in `tannerpolley/ePC-SAFT` is the sole tracker. Read the selected
issue and its comments before implementation; plans and ADRs define design but
do not replace current issue state.

## Historical programs

CAP programs, D-numbered decisions, issue completion and old interface examples
are provenance. They do not establish availability in the rebuilt Engine.
Use the sole accepted ADR sequence under `docs/scientific/adr/` for decisions,
the scientific context for permanent mathematics, and `docs/scientific/implementation-plan.md`
for current implementation/evidence gaps. Read live issue state when issue
work is requested; do not copy it into another status ledger.

## Private development feedback

Build and install the immutable wheel before tests. `tools/check-all.sh` is the
parallel installed-wheel gate with a warm target below ten seconds. Slow Engine
stress cases, analysis reproductions and figures, full native
clang-tidy, and sanitizers run only when the scientific claim or native risk
requires them.

## Task strategy

Use one accountable writer for each Data–Engine–analysis slice and independent
read-only reviewers where the risk warrants them. Persistent tests must
falsify numerical claims; readiness checks, schemas, imports, paths, prose, API
shape, and coverage targets do not become scientific tests. This project has
no CI/CD.

# Engine issue routing

GitHub Issues in `tannerpolley/ePC-SAFT` is the sole tracker. Read the selected
issue and its comments before implementation; plans and ADRs define design but
do not replace current issue state.

## Completed CAP program

Root issue #10 and CAP-01 through CAP-11 are complete under the D-030
through D-033. `solve` is the sole public equilibrium route;
omitted fixed-T,P topology uses the retained internal HELD/HELD2 discovery
algorithm. D-034 admits the issue-#30 fixed-T,P homogeneous reactive
Observation Family, and D-035 admits the issue-#31 fixed-topology reactive
bubble/VLE Observation Family. Nonreactive phase-equilibrium observations
remain unadmitted.

Historical CAP plans, receipts, results, reviews, and source identities remain
provenance. They are not a live execution queue. Issue #6 was superseded when
D-038 established MEA-Thermodynamics as the official downstream application.
Issues #30 and #31 are completed detached post-cutover capabilities with
hash-bound admission receipts.

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

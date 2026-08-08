# Governance issue routing

The sole tracker is GitHub Issues in `tannerpolley/ePC-SAFT`. Infer the current
remote from the attached checkout and use authenticated `gh` operations.

Governance issues own doctrine, manifests, topology, readiness, or an explicit
external operation. Engine production work belongs under its capability root;
[root issue #10](https://github.com/tannerpolley/ePC-SAFT/issues/10) owns the CAP
program. Data, Validation, and Lab are directory owners inside the same branch,
not separate queues or resident workers.

Native subissue and blocked-by relationships define the executable frontier.
Run each unblocked leaf in one fresh task with one accountable writer. The same
task carries implementation, verification, review, commit, push, and issue
closure. Read-only reviewers may inspect stable diffs; parallel writers may not
touch the same CAP slice.

GitHub issue state, Git commits, reviews, integration, and current task-attached
checkouts are lifecycle truth. Plans and ADRs define durable design. This project
has no CI/CD, and labels, milestones, comments, or task summaries grant no
authority by themselves.

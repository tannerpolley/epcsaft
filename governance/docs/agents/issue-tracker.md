# Governance issue routing

The sole tracker is GitHub Issues in `tannerpolley/ePC-SAFT`. Infer the current
remote from the attached checkout and use authenticated `gh` operations.

Governance issues own doctrine, manifests, topology, readiness, or an explicit
external operation. Engine production work belongs under its capability root.
Root issue #10 and CAP-01 through CAP-11 are complete; their issues, plans,
receipts, and reviews are retained provenance rather than a live queue. Data
and Validation are directory owners inside the same branch, not separate queues
or resident workers.

Current issue and blocked-by relationships define the executable frontier. Use
one accountable writer per slice. Read-only reviewers may inspect stable diffs;
parallel writers may not touch the same files.

GitHub issue state, Git commits, reviews, integration, and current task-attached
checkouts are lifecycle truth. Plans and ADRs define durable design. This project
has no CI/CD, and labels, milestones, comments, or task summaries grant no
authority by themselves.

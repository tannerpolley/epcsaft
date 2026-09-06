# Project Context

`tannerpolley/ePC-SAFT` is the private, single-repository home of the active
ePC-SAFT scientific project. Agents
operate only in the checkout or worktree attached to their Codex task; no
absolute filesystem path grants project or Git authority.

The repository contains one runtime distribution, `epcsaft`, under `engine/`.
EOS, equilibrium, and regression are modules of that distribution, not
separate packages. Data packets, production code, installed-wheel analysis
sets, and independent reference cases retain distinct directory ownership while sharing one
branch, one root environment, and one issue tracker.

Current priorities are rapid scientific revision, end-to-end vertical slices,
and low-friction Codex worktrees. Durable release baselines, CI/CD, public ABI
compatibility, and legacy package routes are out of scope.

Developer feedback is build-first and latency-bounded: setup installs the
current wheel before tests, the default installed-wheel gate targets less than
ten seconds on a warm checkout, and publication-grade campaigns or exhaustive
native analysis run only through explicit owner commands.

The current implementation gap plan is `docs/greenfield/TRANSFER.md`.
Permanent mathematics and methods are owned by the documents mapped in
`docs/scientific/CONTEXT.md`. Prior stage documents, issue programs and task
summaries retain historical context, not current scientific authority.

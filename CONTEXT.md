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

The current Engine program is the greenfield rebuild on branch `greenfield`,
mapped in `docs/greenfield/TRANSFER.md`; the earlier Hard-Cutover Plan and root
issue #10 are provenance. Issues #185–#196 hold the scope each stage absorbs;
the stage documents and accepted ADRs define the design. Conversation transcripts and task
summaries are handoff context, not durable authority.

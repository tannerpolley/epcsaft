# Project Context

`tannerpolley/ePC-SAFT` is the private, single-repository home of the active
ePC-SAFT scientific project. Governance doctrine revision 15 applies. Agents
operate only in the checkout or worktree attached to their Codex task; no
absolute filesystem path grants project or Git authority.

The repository contains one runtime distribution, `epcsaft`, under `engine/`.
EOS, equilibrium, and regression are modules of that distribution, not
separate packages. Data packets, production code, installed-artifact campaigns,
and independent oracles retain distinct directory ownership while sharing one
branch, one root environment, and one issue tracker.

Current priorities are rapid scientific revision, end-to-end vertical slices,
and low-friction Codex worktrees. Durable release baselines, CI/CD, public ABI
compatibility, and legacy package routes are out of scope.

The current Engine program is the
[General Scientific Capability Hard-Cutover Plan](engine/docs/plans/2026-08-07-general-scientific-capability-hard-cutover.md).
[GitHub root issue #10](https://github.com/tannerpolley/ePC-SAFT/issues/10), its
native CAP subissues, and their blocked-by relationships are lifecycle truth;
accepted plans and ADRs define the design. Conversation transcripts and task
summaries are handoff context, not durable authority.

# Project Context

`tannerpolley/ePC-SAFT` is the private, single-repository home of the active
ePC-SAFT scientific project. The local canonical checkout is
`~/Workspaces/Engineering/ePC-SAFT-project`.

The repository contains one runtime distribution, `epcsaft`, under `engine/`.
EOS, equilibrium, and regression are modules of that distribution, not
separate packages. Data packets, production code, installed-artifact campaigns,
and independent oracles retain distinct directory ownership while sharing one
branch, one root environment, and one issue tracker.

Current priorities are rapid scientific revision, end-to-end vertical slices,
and low-friction Codex worktrees. Durable release baselines, CI/CD, public ABI
compatibility, and legacy package routes are out of scope.

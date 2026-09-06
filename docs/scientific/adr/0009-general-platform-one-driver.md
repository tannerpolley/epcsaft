---
status: accepted
date: 2026-09-06
---

# General platform with one equilibrium driver

Kind: computational-architecture

The owner selected a general fluid thermodynamics platform whose first required
acceptance uses are the three official applications. The owner explicitly
requires one equilibrium driver and selected Git history as the sole retention
location for obsolete setup and execution documents.

The [project map](../CONTEXT.md) defines the current responsibilities and scope.
All admitted equilibrium specifications use one driver and shared mathematical
rows; separate application/problem-specific equilibrium solvers are excluded.
This replaces the proposal to allow independently selected equilibrium drivers.
It preserves lower-level density/association kernels and the distinct search
and continuation operations that reuse equilibrium calculations.

The existing native ownership principles remain: one Python-free native
archive, one extension, no public C++ ABI or iterative Python callbacks, explicit
physical inputs, common parameter identities and owner-supplied derivative
chains. General platform scope does not introduce compatibility support,
external catalogs, speculative extension frameworks or a second installation.

The choice does not qualify the existing Ipopt Hessian or select an untested
replacement method. Numerical method selection and demanded-case acceptance
remain work in the implementation plan. Material limitations require an owner
decision rather than a silent second driver. Solids remain a separately
justified extension. Runtime budgets retain their current numeric owner.

Obsolete stage designs, briefs, worker records and superseded architecture
ADRs are preserved in Git at `2985cd4`; current scientific evidence remains
with its owner. No historical worker roles, launch commands or per-stage
budgets govern future work. The installed CSE workflow and repository guidance
own execution.

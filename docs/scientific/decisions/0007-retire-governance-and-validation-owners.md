---
status: accepted
date: 2026-09-01
---

# Scientific directory ownership

Kind: computational-architecture

Scientific inputs and calculations are organized by their actual responsibility,
with one owner for each value and requirement. The
[project map](../CONTEXT.md) defines current directories. Execution instructions
live in `AGENTS.md`; native topology and source budgets in
`ARCHITECTURE.yaml`; reusable packets and locks in `data/`; scientific
studies in `analyses/`; equations and algorithms in `engine/docs/`.

Retain one-question analysis sets with their actual inputs, scripts, results,
figures and interpretation. Preserve the established
`analyses/reproduction/<short-id>/` layout. Do not create a miscellaneous
administrative directory that duplicates those owners, GitHub state or the
GitHub issue roadmap.

The retired directory names, decision crosswalk and migration evidence remain
in Git history. Current setup has no aliases, redirects or competing decision
sequence. This ownership choice changes no numerical value, equation, parameter,
unit, standard state or scientific acceptance criterion.

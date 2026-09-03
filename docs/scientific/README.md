# Scientific records

This directory is the CSE record owner for `tannerpolley/ePC-SAFT`: the map of the
scientific question, the owners of every kind of authority, the repository vocabulary, and
the accepted decisions. It copies no equation, tolerance, or result; each stays with the
owner `CONTEXT.md` names.

| file | what it owns |
|---|---|
| `CONTEXT.md` | the question and intended use, the owner map, the glossary with sources, the avoid table, the default deliverable |
| `adr/` | the sole accepted decision-record sequence, `0001` onward; superseded records stay and name their replacement |
| `writing/STYLE.md`, `writing/RULES.md` | house style and the deterministic prose rules for notebooks and manuscripts |
| `numerical-policy-audit.md` | the 2026-09-01 audit of every tolerance in the old tree; the input `docs/greenfield/CORE.md` absorbed |

The rebuild in progress is recorded under `docs/greenfield/`, starting at
`docs/greenfield/README.md`. Its stage documents are the current design authority; the
ADRs here that it supersedes say so in their front matter.

Layout follows CSE issue 55: two mandatory Markdown files, optional records as Markdown,
LaTeX only for manuscript sources under `docs/scientific/latex/` when one exists. The
earlier `formulation.tex`, `methods.tex` and `evidence.tex` maps were removed at this
commit; every path they pointed at had been deleted by the rebuild, and their surviving
content is the owner map in `CONTEXT.md`.

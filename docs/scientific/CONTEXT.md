# Scientific context

default deliverable: capability

## Question and intended use

The engine evaluates ePC-SAFT states, solves phase and reaction equilibrium from a
degrees-of-freedom specification, and fits parameters to declared observations, through one
installed `epcsaft` wheel. The intended use is generic scientific computation for the
formulations the equation registry admits; application chemistry, parameter adoption and
downstream claims stay with the applications that consume the wheel.

The tree is being rebuilt on branch `greenfield` in eight gated stages
(`docs/greenfield/TRANSFER.md`). Four are closed. Until the final check merges,
`main` is the old tree and `greenfield` is the design authority.

## Owners

| authority | owner |
|---|---|
| governing equations, symbols, units, validity domains | `engine/docs/science/full-eos.md` |
| stable equation identifiers, each tied to its `// EqID:` annotation in source | `engine/docs/equations_registry.yaml`; rendered view `engine/docs/equations.md` |
| the rebuild: stages, gates, roles, carry-over, decisions taken while writing | `docs/greenfield/` |
| what every layer may assume of the core: declared choices, convergence records, the trace | `docs/greenfield/CORE.md` |
| the equilibrium formulation and its count rule | `docs/greenfield/EQUILIBRIUM.md`; derivation in `engine/docs/equilibrium/formulation.md` |
| structural rules the checker enforces | `tools/check-engine-architecture.py`; stated in `docs/greenfield/CANON.md` |
| acceptance data: measured rows only | `engine/tests/oracle/oracle.json` |
| reusable source data with provenance and hashes | `data/packets/`; contract `data/docs/packet-contract.md` |
| evidence produced by asking a question of the engine | `analyses/`; reproductions under `analyses/reproduction/<year>-<author>` |
| bibliography snapshot, generated from Zotero, for manuscripts only | `engine/docs/latex/references.bib` via `engine/scripts/docs/sync_bibliography.py`; the equation registry does not read it |
| accepted decisions | `docs/scientific/adr/` |
| runtime semantics, concrete lifetimes | `engine/native/` and `engine/src/epcsaft/` |

Tests and analyses own numerical claims. Documents and architecture views are derived
evidence: they expose structural drift and prove nothing numerical.

## Glossary

Every term names its source. A term with no source is retained vocabulary from the old
tree, kept only until the rebuild retires the text that uses it.

| Term | Meaning | Source | Scope |
|---|---|---|---|
| stage | one unit of the rebuild: one brief, one worktree, one gate, one acceptance | `docs/greenfield/TRANSFER.md` §Execution order | `docs/greenfield/**` |
| brief | the pair `stage-NN.brief.{md,json}`: sources with markers, budgets with ceilings, gates with mutations | `docs/greenfield/ORCHESTRATION.md` | `docs/greenfield/briefs/**` |
| gate | an executable claim a stage must make green | `docs/greenfield/CANON.md` | `docs/greenfield/**` |
| mutation | a named edit that must turn one gate red; acceptance applies every one | `docs/greenfield/ORCHESTRATION.md` | `docs/greenfield/**`, `tools/preflight-brief.py` |
| oracle | measured rows only; a fact about the world no implementation may fail | `docs/greenfield/TRANSFER.md` §Oracle | `engine/tests/oracle/**` |
| canon, core, surface, engine, analyses | the five categories a file belongs to, and the commit categories | `docs/greenfield/CATEGORIES.md` | repository |
| declared choice | a selector with no unspecified member; absence is carried by `declared()` | `docs/greenfield/CORE.md` §4 | `engine/native/**` |
| convergence record | achieved residual, criterion, floor, amplification and a named stop, returned by every solver | `docs/greenfield/CORE.md` §5 | `engine/native/**` |
| floor | the accuracy a nested solve's inputs allow, amplified through their conditioning | `docs/greenfield/CORE.md` §3 | `engine/native/**` |
| derived tolerance | a tolerance constructed from an upstream one by a factor with a stated reason; never tighter than its input | `docs/greenfield/CORE.md` §3 | `engine/native/core/numerical_policy.hpp` |
| scope, trace, capture | one span per instrumented call; one root per bound function; the lexical capture around a call | `docs/greenfield/CORE.md` §6 | `engine/native/**`, `engine/src/epcsaft/**` |
| count rule | rows − unknowns = k − f; square iff f = k | `docs/greenfield/EQUILIBRIUM.md` §2 | `engine/native/equilibrium/**` |
| like-ion dispersion convention | u_ij = 0 for cation–cation and anion–anion pairs, declared at the surface | `docs/greenfield/EOS.md`; Held 2008, Held 2014, Figiel 2025 | `engine/native/**`, `data/packets/**` |
| packet | a versioned, hash-bound, source-traceable data unit | `data/docs/packet-contract.md` | `data/**` |
| reproduction | an analysis leaf that re-derives a published result against its source | `analyses/README.md` | `analyses/reproduction/**` |
| certificate | — | none | `engine/CONTEXT.md` |
| fingerprint | — | none | `engine/CONTEXT.md`, `data/AGENTS.md` |
| admission | — | none | `engine/CONTEXT.md` |
| typed failure | — | none | `engine/CONTEXT.md` |
| claim boundary | — | none | `engine/CONTEXT.md`, `AGENTS.md` |

The five sourceless terms live in `engine/CONTEXT.md`'s glossary, written for the old tree
and not yet retired. Retiring or sourcing them is a Design item for the stage that rewrites
that file, not a Setup edit.

## Avoid

| Avoid | Prefer | Meaning | Scope | Exceptions | Evidence |
|---|---|---|---|---|---|
| saturation solver, flash routine, bubble-point route | equilibrium specification, the cell | one driver serves every degrees-of-freedom cell; a named route is the special case the rebuild removed | `engine/**`; `docs/greenfield/**` | historical text under provenance | `docs/greenfield/TRANSFER.md` §Execution order, owner decision 2026-09-02 |
| HELD, HELD2, GREPE | phase stability, equilibrium | no proper nouns in identifiers; sources are citations | `engine/native/**`; `engine/src/**` | citations, `analyses/**` results, ADR text | `engine/ARCHITECTURE.yaml` naming; owner rule |
| golden file | oracle row, closed form, property | a number this code produced is not evidence this code is right | `engine/tests/**` | — | `docs/greenfield/CANON.md` |
| iteration limit reached | named stop from `classify` | a stop is one of the named outcomes, never hand-set except left_domain, singular, infeasible | `engine/native/**` | — | `docs/greenfield/CORE.md` §5.2 |

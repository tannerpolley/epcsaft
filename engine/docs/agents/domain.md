# Domain Docs

This is a single-context repository. Engineering skills should consume its
domain documentation before exploring the codebase.

## Before exploring, read these

- `CONTEXT.md` at the repository root.
- Relevant decisions under `docs/adr/`.

If either is absent, proceed silently. Domain-modeling workflows create
documentation only when terminology or decisions are actually resolved.

## File structure

```text
/
├── CONTEXT.md
├── docs/adr/
└── src/
```

## Use the glossary's vocabulary

Use terms as defined in `CONTEXT.md` in issue titles, refactor proposals,
hypotheses, and test names. Avoid synonyms the context explicitly rejects.

If a needed concept is absent, reconsider whether the term belongs to this
domain or note the gap for domain modeling.

## Flag ADR conflicts

Surface any conflict with an existing ADR explicitly rather than silently
overriding it.

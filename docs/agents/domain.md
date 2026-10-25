# Domain Docs

How the engineering skills should consume this repo's domain documentation when exploring the codebase. Single context: one root glossary.

## Before exploring, read these

- **`GLOSSARY.md`** at the repo root: project vocabulary in Matt's format, shared with the CSE plugin, which requires a `Source:` line per term.
- **`docs/scientific/decisions/`**: read the accepted decisions that touch the area you're about to work in.
- **`docs/scientific/CONTEXT.md`**: the CSE scientific map (authority, owners, review rules). CSE owns it; read the section the task needs.

## Use the glossary's vocabulary

When your output names a domain concept (in an issue title, a refactor proposal, a hypothesis, a test name), use the term as defined in `GLOSSARY.md`. Write the glossary's bold term rather than a synonym it lists under `_Avoid_`.

If the concept you need isn't in the glossary yet, either you're inventing language the project doesn't use (reconsider) or there's a real gap: add the term through `/domain-modeling` with its `Source:` line.

## Flag decision conflicts

If your output contradicts an accepted decision, surface it explicitly rather than silently overriding:

> _Contradicts decision 0009 (one equilibrium driver), but worth reopening because…_

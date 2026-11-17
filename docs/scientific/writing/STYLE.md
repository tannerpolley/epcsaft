# Scientific writing style
- Lead with the scientific quantity, mechanism, comparison, or result.
- Use active voice and name the actor, method, or calculation.
- State equations, conditions, units, and uncertainty close to the result they govern.
- Prefer precise verbs and concrete nouns.
- Vary sentence length while keeping logical relations explicit.
- Distinguish observation, calculation, interpretation, and speculation.
- State limitations where they constrain a claim.
- Apply the manuscript or analysis-notebook profile while preserving project
  scientific meaning. Follow a publication class only when the user supplies
  that requirement.

## Required

- Cite every source-dependent equation, parameter, dataset, and substantive literature claim.
- Give exact source locators for inputs that drive a calculation.
- Report units, basis, tested domain, metric, and uncertainty or tolerance where relevant.
- Keep code verification, solution verification, model validation, and uncertainty distinct.

## Forbidden patterns

- Throat-clearing and canned transitions.
- Vague declarations that a result is important without naming its consequence.
- False binary contrasts used for drama.
- Business jargon, filler adverbs, stacked sentence fragments, and rhetorical setups.
- Em dashes.
- Claims beyond the tested domain or source support.

## Manuscript deliverables

| Scope | Standing owner decision | Authority |
| --- | --- | --- |
| Project-wide | Deliver manuscripts only as final products. | Owner, 2026-10-15, “manuscripts are only final products ever.” |

## ePC-SAFT engine manuscript

Sources in the private validation checkout are `validation/manuscript/main.tex`,
`validation/manuscript/supplement.tex`, and `validation/manuscript/figures/`.

| Scope | Standing owner decision | Authority |
| --- | --- | --- |
| Main article and supporting information | Deliver the main article and supporting information as separate PDFs. | Owner, 2026-10-28, recorded in issue #37. |
| Supporting information | Put the long figure-comparison table in its own file. The supplement points readers to `supplementary-comparisons.csv`. | Owner decision recorded in the 2026-11-16 review in issue #37. |

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
| Audience and framing | Write the paper that introduces a released open-source thermodynamics package, titled as a package. The same text is the first chapter of the owner's dissertation, so it must read plainly and set the example for later papers. | Owner, 2026-11-19, dissertation chapter 1 and plain-language package paper; owner, 2026-11-27, "easy enough to read and not get lost in the details." |
| Model papers | Measure the structure against the FeOs package paper and other papers that present a full release of a thermodynamics package. | Owner, 2026-11-27, "papers that are similar to this and the feos one." |
| Proportion | Show each capability at the depth a reader needs. Put the detail of each analysis in the supporting information. | Owner, 2026-11-27, "the most bang for buck without diving into the details and nitty gritty of every bit of analysis." |
| Results display | Each section that runs inputs to outputs shows them in a table, plot or equation. Prose keeps only the numbers its claim rests on. Phase equilibria of three-component systems appear on a ternary diagram. | Owner, 2026-11-27, "I know that AI agents like dense text but us humans do not." |
| Tables | Keep table cells minimal and readable at a glance. Put details in the caption or the supporting information. | Owner, 2026-11-27, "easy to read at a glance." |
| Figures | Draw calculated model results as dashed lines and measurements as markers such as dots, crosses or circles. | Owner, 2026-11-27, "model fit data points should almost always be with a dashed line while measured should be the dots, x's or circles." |
| Code listing | Show short routine front-end calculations, not a complicated special case. | Owner, 2026-11-27, the listing was "too long and too hard to read" and "uses a complicated unique example." |
| Verification | Present verification results, not the investigation that produced them. | Owner, 2026-11-19, report the result, not the path; owner, 2026-11-27, "super long text for very little gain." |
| Referee-requested detail | Put detail requested in review in a table row, a figure or the supporting information. Main text grows only when a main claim needs the detail. | Owner, 2026-11-27, approved the retrospective's fixes: "yes do all 1-7." |
| Captions | Lead with the result. Main-text captions run about 60 words or fewer; methods belong in the text or the supporting information. | Owner, 2026-11-27, approved the retrospective's fixes: "yes do all 1-7." |
| Precision | The abstract and main text print values at the precision that changes a reader's conclusion. Roundoff-scale residuals and certificates go to the supporting information. | Owner, 2026-11-27, approved the retrospective's fixes: "yes do all 1-7." |

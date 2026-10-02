# Core readiness and existing implementation review

Date: 2026-09-09. Scope: read-only assessment of C1–C4 by the requested
Fable 5.1 High reviewer, followed by lead inspection of the cited owners and
callers. At that snapshot, all four rows were Pending. This is not implementation
or numerical qualification of Core.

| Row | Finding | Required evidence before qualification |
|---|---|---|
| C1 — Constants and numerical policy | `gibbs_minimizer.cpp` applies `STABILITY_COMPOSITION_FLOOR` to an absolute phase amount although its policy description names a composition bound. Scaling a free-phase feed below that absolute threshold can mark its phases vanished. `FLOOR_MARGIN` and `STALL_WINDOW` explicitly have guessed sources, contrary to the header's claim that unjustified values are absent. | Establish the dimensional meaning and scaling of the phase-vanishing criterion, qualify the chosen policy against justified cases, and resolve unsupported source claims. |
| C2 — Parameter coordinates and scaling | Six active parameter families already exist in `helmholtz_tape.hpp`. Validation checks index range and duplicate indices, but a single pair coordinate can omit its symmetric matrix entry. Existing callers can assemble symmetric pairs manually; the common coordinate/scaling contract is not qualified. | Define the shared identity and scaling contract and falsify incorrect pair mapping and derivative transformations. Do not invent a new public API solely to match the roadmap row. |
| C3 — Diagnostics and availability | `ALG-TERMINATION` requires nonfinite and singular outcomes. `Stop` has no nonfinite state, and `classify` does not reject NaN/nonfinite input before comparing it. General unsuccessful Ipopt statuses enter the same classification route. Density refinement can throw on its iteration limit even after classification reports convergence. | Exercise nonfinite/domain/singular/infeasible outcomes and boundary iterations through their actual producers; retain valid outputs and distinct causes. Existing finite-input classifier checks do not establish these cases. |
| C4 — Core qualification | The installed diagnostics suite passes, but code budgets and existing finite-input cases are not acceptance of the complete Core contract. | Complete C1–C3 evidence, cost review and independent acceptance before qualifying C4. |

Historical source anchors preserved in baseline `8b6835c`: [numerical policy](https://github.com/tannerpolley/ePC-SAFT/blob/8b6835c58a85773059dadb7d1e089e0d7ec07380/engine/native/core/numerical_policy.hpp),
[phase amount use](https://github.com/tannerpolley/ePC-SAFT/blob/8b6835c58a85773059dadb7d1e089e0d7ec07380/engine/native/equilibrium/gibbs_minimizer.cpp),
[active coordinates](https://github.com/tannerpolley/ePC-SAFT/blob/8b6835c58a85773059dadb7d1e089e0d7ec07380/engine/native/model/helmholtz_tape.hpp),
[coordinate validation](https://github.com/tannerpolley/ePC-SAFT/blob/8b6835c58a85773059dadb7d1e089e0d7ec07380/engine/native/model/helmholtz_model.cpp),
[termination doctrine](https://github.com/tannerpolley/ePC-SAFT/blob/8b6835c58a85773059dadb7d1e089e0d7ec07380/engine/docs/science/algorithms.md#alg-termination),
[classifier](https://github.com/tannerpolley/ePC-SAFT/blob/8b6835c58a85773059dadb7d1e089e0d7ec07380/engine/native/core/convergence.hpp) and
[density refinement](https://github.com/tannerpolley/ePC-SAFT/blob/8b6835c58a85773059dadb7d1e089e0d7ec07380/engine/native/model/density.cpp).

The review also checked a suspected gap in `packing_terms`: the omitted ionic
sigma/dispersion derivatives are consistent with the separate ionic packing
diameter. It is not a defect established by this review. The existing D3
readiness case already checks rejection of a tightening derived tolerance;
removing its duplicate Planned declaration does not qualify C1.

Current plans and prerequisites belong to the [GitHub issue roadmap](https://github.com/tannerpolley/ePC-SAFT/milestones).
No constants, tolerances, algorithms or Core
runtime files were changed by this review. The existing F5 disagreement remains
Failed and is not resolved by infrastructure acceptance.

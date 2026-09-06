# Current transfer gap plan

Status: incomplete. Baseline inspected: `greenfield@ff614901be22727095947e71fb4626ef2dae7a7c`.
Documentation reconciled 2026-09-06; no numerical implementation changed by that
reconciliation. The [dated audit](reviews/transfer-audit-2026-09-06.md) retains
exact source anchors, commands, measured results and independent review. Its
F1–F5, E1–E12 and S0–S9 identifiers below refer to that snapshot. Requalify a row
when its code changes; a worker report or equation registration is not acceptance.

## Settled scope and authority

All three official applications migrate to the new general Engine interface.
Retain accepted scientific/runtime thresholds with their case/domain; missing
thresholds require owner decisions. Solids remain a separately justified
extension. Application repositories consume a non-editable wheel identified
by SHA-256; no sibling-source imports or compatibility implementation.

Permanent authority is the [scientific owner map](../scientific/CONTEXT.md):
[equations](../../engine/docs/equations.md),
[EOS conventions](../../engine/docs/science/full-eos.md),
[equilibrium formulation](../../engine/docs/equilibrium/formulation.md), and
[algorithms](../../engine/docs/science/algorithms.md). Prior stage documents,
briefs, launch instructions, performance targets and dated reviews are history,
not competing live doctrine. Supported audit corrections are requirements;
unresolved method selections remain proposals.

## Code and evidence map

| Requirement | Existing owner / implementation at baseline | Gap and numerical falsifier |
|---|---|---|
| Common residual energy and exact derivative chains; ALG-ASSOCIATION | `engine/native/physics/helmholtz.hpp`, `association.hpp`, `electrolyte.hpp`, `born.hpp`; `engine/native/model/helmholtz_tape.hpp` | Common tape exists. F2: weighted association elimination fails exact-zero associating support; E2 must separate regular mass action from the genuinely nonsmooth DH insertion Hessian. |
| EqID: `residual_internal_energy`, residual/total caloric identities | `engine/native/model/helmholtz_model.cpp` | F1: missing T factor in residual u; E1 identity and independent temperature derivative required. Total h/u must not be substituted by residual fields. |
| ALG-ADMISSION | `engine/native/surface/parameter_builder.cpp`, `engine/src/epcsaft/` | F3: neutral subset, charge defaults, source domains and ordinary electroneutrality incomplete; E2 plus explicit selected-support/domain cases. |
| ALG-DENSITY | `engine/native/model/density.cpp` | Gradient scan and corrected reduced slope exist. F4: public branch/count not preserved; E3/E5 root multiplicity, disappearance and branch tests. |
| Equilibrium rows and ALG-EQUILIBRIUM | `engine/native/equilibrium/request_compile.cpp`, `gibbs_minimizer.cpp`, `phase_boundary.cpp` | Compiler/shared minimizer exist. F5 exact pinned Hessian and callback objective-factor weighting are incorrect. E12 directional weighted-gradient ladder at nonzero multipliers, obj_factor 0, 0.5 and 1, including free T/P away from seeds; zero curvature is not exact. Boundary seeds are not continuation. |
| EqID: `implicit_sensitivity`; ALG-ACTIONS | `engine/native/equilibrium/implicit_sensitivity.cpp` | A1 path exists with silent unavailable cases; E3 needs reasons/rank/linear residual and complete mixed chains. A2/A3 required by downstream uses are absent or unqualified; E10 independent action ladder. |
| EqID: `standard_pressure_density_conversion`; ALG-REFERENCE | `engine/native/model/infinite_dilution_reference.cpp`; no complete common thermochemistry owner | E4/E7/E9: consistent fixed-P0 conversion, typed K(T), total properties and reference derivatives. Reaction gauge alone cannot supply absolute calorics. |
| ALG-FIT | `engine/native/regression/` placeholder | General observation algebra/fitting absent; E9 property fits, fixed-coordinate zero derivative, domain failures, chain rules and conditional covariance. Extend active coordinates beyond the six present EOS families only with complete defining chains. |
| Stability/discovery | No qualified current discovery owner | Discover mapping incomplete/rejected. E8/E11: instability witnesses, declared-vs-discovered comparisons, charge/support and genuine bound meaning; local inner minima and exhausted search are not global proof. |
| Installed-wheel/application evidence | `tools/verify-installed-wheel.py`, `tools/check-all.sh`; three downstream repositories | Verifier blocked by missing regression. Ten public tests passed in audit, architecture scan passed but does not establish math or actual artifact arguments. Native configured graph was stale. All-application transfer acceptance remains open. |

The table names intended numerical falsifiers, not newly executed passing tests.
Existing compact checks live under `engine/tests/`; exact existing/proposed
case dispositions and E1–E12 are in audit §7. Do not infer implementation from
an EqID annotation or a source file's existence.

## Reference and application evidence

The 98 retained reference rows have different meanings (audit §7): Hamer–Wu
24 evaluated coefficients; NIST ethane 10 reference-EOS values; Cai 6 rows at
two pressures; speciation 5 Jakobsen and 2 Boettinger rows with differing
reported quantities; reactive VLE 6 Jou, 2 Hilliard and 2 Xu rows with distinct
total/partial-pressure conventions; 6 Amundsen density rows with pressure/basis
assumptions; 8 Khudaida measured tie endpoints with constructed midpoint feeds;
27 Figiel digitized published-model curve points. They are not all measurements.
A reported zero without a detection bound must not be converted into a known
censoring bound or exact mathematical zero.

MEA-Thermodynamics needs reactive equilibrium, caloric observations and fits.
MEA-Absorption-Column needs reference-consistent total properties and derivatives;
the conservative residual's first exact Jacobian already requires A2 actions,
even with a limited-memory outer Hessian. Lithium_Extraction needs coupled
reactive LLE, feed/extract observables and declared/discovered phase evidence.
Retained acceptance domains, commands and missing criteria are in audit §6.
No application was qualified by the documentation cleanup.

Missing May2015 packet inputs, MEA revision-4 R1–R3 source correlations,
Khudaida association data, source calorics and unresolved diluent identity
remain scientific input gaps. No zero/default coefficient is authorized.

## Completion sequence

These are bounded **proposals**, not instructions to launch workers. Keep one accountable writer for each slice across its exact scientific owners; application changes remain separately reviewable in their repositories. No new governance registry, framework, compatibility layer or CI system is needed. Each slice gets independent scientific/numerical review against its falsifier before implementation under the normal CSE route.

| Priority / dependency | Scientific question and smallest change/deletion | Inputs and falsifier | Command / downstream acceptance / exit |
|---|---|---|---|
| **S0 — honest executable baseline** | Can the installed artifact be identified and the actual current numerical checks execute? Correct existing verifier's stage-ahead module/old entrypoint names, preserve all provenance requirements; make build graph/stamp and checker artifact arguments refer to the actual current build. Exclude nonexistent regression honestly until implemented; do not add stubs. | Current wheel, source/build metadata; fail on wrong origin/editable/missing actual binding; include F1 as a failing numerical case before correction. | `bash tools/ensure-engine-wheel.sh`, N:physics/model/diagnostics/specification as the proposed fast subset; N:equilibrium/boundary only as explicit owner-run evidence after F5 qualification; existing root gate after repair. Exit: current build and actual executed checks reproducible; no “fully complete Engine” claim. |
| **S1 — correct EOS/caloric and support admission** depends S0 | Do public states obey registered identities and evaluate valid support limits? Fix F1 at model producer; repair association zero-support elimination; select packet declarations consistently and reject missing physical input. Carry domains/physical-vs-extrapolation semantics without fake defaults. | E1–E3, complete MEA/Gross inputs and boundary assumptions. | N:model/physics, P:eos; one MEA liquid and one ionic/disjoint-support limit. Exit: independent identities/step ladders pass; typed genuine singularity remains distinguishable. |
| **S1b — qualify equilibrium Hessian/method** depends S0 | Does the driver provide the Hessian it advertises? Correct objective/constraint weighting; choose an honest pinned equality method using shared rows, removing the zero-as-exact claim. If choosing limited-memory, revise the exact-only architecture rule and unconditional uses_exact_hessian capability claim in the same slice, preserving exact-Jacobian and scientific checks. | F5 weighted-gradient directional ladder at nonzero multipliers and factors 0/.5/1; existing boundary and finite-phase reference cases. | Native focused callback check, then explicit N:equilibrium/boundary. Exit: correct advertised derivative or explicitly qualified alternative method, same residual/branch criteria and measured cost; no new named bubble/dew owner. |
| **S2 — retain actual cell evidence** depends S0; independent of most S1 implementation | Does every executed numerical claim have inputs, state/failure and matching deviations? Finish the corrected existing harness with one runtime recording boundary; remove unnecessary duplicate schema/value copies. Keep numerical assertions unchanged except separate S1 fixes. | Dirty brief preserved/reviewed; 26 registrations, each actual solve including failures, finite signed deviations and bars. | Existing native checks + proposed summary/claims validation after those tools exist. Future isolated meaningful mutations only: omit one of two solves; emit finite wrong deviation with true; omit required record. Exit: complete cell accounting and honest result interpretation, not source-name presence. |
| **S3 — source and equation completeness** depends input review, can proceed alongside S1/S2 | Which retained data can support the required physical claims? Resolve missing May packet, MEA R1–R3, Khudaida association, dielectric coefficients/domains and HBTA/TOPO gaps. Reconcile registry entries/actual owners; relabel oracle classes. | Retained source CSVs/locators; new acquisition only through separately authorized user-operated source flow. Never fill gaps from an application name or guessed constant. | `data/tools/validate_packets.py` plus actual packet model construction; E6–E8 input audit. Exit per packet: complete typed input or explicit scientific block. A blocked chemistry dataset need not block unrelated numerical fixes. |
| **S4 — reference thermochemistry and general specifications** depends S1/S1b and relevant S3 | Are g°, K(T), h/cp and balances consistent in the same declared reference? Transfer later old thermochemical behavior into model/property/shared-row owners, not old Python orchestration; add explicit total/residual outputs and thermal correlations, energy specification and sensitivity availability. | Old PR199–203/3375ead evidence, E4/E5/E7, declared reference data and conservation basis. | N:model/specification/equilibrium/boundary, existing public tests; MEA finite-dose heat and incipient K(T) slope calculation. Exit: analytic reference conversions and caloric/balance identities; thermal sensitivity conditioning/failures visible. |
| **S5 — demanded higher actions** depends S1/S4; may precede final data fits | Can exact derivatives consumed by the absorber be computed on the same solved branch? Add compact native A2 equilibrium and second caloric actions first; preserve old ordinary A3 domain with a focused replay. Reuse linear factorization/directional work; no general dense tensor API by default. | E3/E10, old higher-state-actions contract, complete directions/units and reference cp slopes. | Native small action falsifiers + installed-wheel absorber node/callback checks. Exit: public action requests/results carry order, complete canonical coordinate identities/units, direction vectors, output identities/units, and typed per-action/caloric failures preserving valid central values; mixed step ladders, conservation/site recurrences and factorial conventions pass; singular/unsupported thermal actions do not erase valid values. Full reactive/caloric A3 scope remains explicitly qualified. |
| **S6 — light fitting** depends S1 and the property types each observation consumes; S4/S5 for caloric derivatives | Can declared weighted observations fit identifiable EOS/reaction/reference coefficients without public evaluator compilation? Implement actual Ceres wrapper, minimal expression algebra, typed active/fixed/shared fields, bounds/transforms and local prediction reuse. Delete planned fixed-coordinate rejection and unconditional covariance-scale rule. | E9; real E6/E7/E11 data only after S3, full coupled chain rule and failures. | P:regression after real module exists; update verifier to current real fit binding; explicit MEA/lithium fits. Exit: identifiable manufactured recovery, J and weighting correct, redundant parameters/bounds/failures honest; real fit adequacy reported separately. |
| **S7 — stable/discovered fluid phases** depends S1/S1b/S4 and complete relevant S3 data | Can declared solutions be tested for lower-energy phases and can missing phases be found within a stated search domain? Implement one controller using shared minimizer/charge/support and scalar density owner; retain distinct LP/global/local methods only where evidence justifies them. Correct Discover mapping; add disappearance/continuation evidence. | E5/E8/E11; exact source equations, physical phase/support domains, clear local/global/search-budget semantics. | Explicit `tools/run-slice.sh khudaida`, native small stability representative after timing; lithium declared-vs-discovered comparison. Exit: found witnesses honored, supports/charge conserved, budget exhaustion never labeled global proof. |
| **S8 — full application qualification** depends required S3–S7 paths | Does each application reproduce its accepted calculation using the new interface and wheel? Migrate all three; retain application transport/chemistry/fit studies locally. Remove superseded adapters only when no consumer remains. | Accepted criteria and independent holdouts; exact wheel identity, canonical input/source data; E7/E10/E11. | Each repository's own integration check and explicit analysis commands after mapping; no sibling-source imports. Exit: all three accepted uses meet retained criteria, or named scientific input gaps remain a truthful block to final transfer completion. |
| **S9 — cost, prevention and closeout** after matching correctness exists; narrow reuse can occur earlier | Did old root causes stay removed and did consumed calculations get faster? Remove duplicate site solve in `state()`, reuse callback primal/derivatives, reuse references within point/fit, then benchmark equivalent work. Apply only the scoped prevention corrections in audit §10. | audit §8 matched cases and all preserved numerical outputs; no timing bar for failed/changed work. | Existing affected/all checks and explicit benchmark; application runtimes against accepted budgets when known. Exit: measured comparison, no uncovered required capability, docs describe implemented/evidenced status; cleanup audit complete. |

Sequence is a dependency graph, not a demand to complete every source acquisition before any code correction. Harness accounting and input recovery can proceed independently; caloric fittings wait for properties; absorber exact Jacobian waits for A2; stable-fluid claims wait for stability evidence. An initial direct-EOS regression subset may land before all caloric work if its limited scope is explicit. The old stage-06 completion alone cannot complete the transfer.


## Code size and consolidation

[Code budgets](../scientific/code-budgets.md) constrain the completed runtime
with a floor target, soft ceiling and hard ceiling. These requested engineering
estimates replace the historical stage ceilings as current size authority.
Neither spare lines nor a smaller implementation establish scientific readiness.
Historical stage work is reconciled in the [sanity review](reviews/reconciliation-sanity-2026-09-06.md).

## Open decisions and completion discipline

Pinned/free-intensive numerical method, discovery controller/global bound,
higher-action storage/algorithm, missing physical inputs and missing application
accuracy/runtime limits remain explicit decisions. Qualify exact pinned
curvature, a Jacobian root method or a limited-memory method against the same
mathematical requirements before choosing. No current exact-only claim may
survive selection of limited memory. Preserve ordinary accepted A3 domains;
broader reactive/caloric A3 requires independent evidence and measured cost.

For every slice, update or explicitly revalidate equation/algorithm IDs in the
same change as code and evidence, then revise this map with the inspected
commit/wheel and actual outcome. Inherited thresholds keep their original
observable/domain; a new threshold needs owner justification. Runtime comparisons
must perform equivalent accepted work. Documentation checks establish only
readiness. Full transfer completion requires all three applications' accepted
uses and the required installed-wheel evidence, with no hidden missing outputs.

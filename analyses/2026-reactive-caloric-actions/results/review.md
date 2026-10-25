# Independent delivered review

Target: GPT-6.1 Sol, reasoning effort `xhigh`, service tier `default`.
T3 task: `issue225-slice1-delivered-review-b02d855e`, request `5c7e9e40-ce2c-4bda-8749-f2201ef53432`.
This is the single independent read-only delivered review requested for slice 1.

**Supported for issue #225 slice 1.** The delivered Properties code and final retained F8/F9 numbers support the bounded reference derivatives on the stated manufactured cases. I found no required defect or unsupported claim within this slice.

Reviewed commit `b02d855ea750aa3a4a8d0ce5d46fd26133bbd06f` against `5a595e9276942ee9532e6a0489ae3b13e6f552e0`, including the complete live [#225 body](https://github.com/tannerpolley/ePC-SAFT/issues/225), all four comments, required instructions, scientific owners, complete changed functions and direct calling paths. The issue body’s normalized SHA-256 matches the retained provenance: `d7a283323f8eeb3e15173c9f3e73b14baad875dc2de42e496b9ce22af9418cf7`.

The implementation and numerical evidence agree:

- [reference_taylor](/home/tnnrpolley21/.t3/worktrees/ePC-SAFT/t3code-225-slice-1-build-properties-reference-6cadaf6a/engine/native/properties/helmholtz_model.cpp:1185) separates the fixed-volume amount insertion from the pressure-closed base path. The degree-seven/nine half-difference excludes repeated insertion terms from the extracted coefficients. The recurrence retains \(V_1\) through \(V_4\), and the existing implicit site calculation advances through degree nine.
- [reference_partial and reaction_references](/home/tnnrpolley21/.t3/worktrees/ePC-SAFT/t3code-225-slice-1-build-properties-reference-6cadaf6a/engine/native/properties/infinite_dilution_reference.cpp:568) restore the factorial and physical-coordinate scales correctly. New columns retain independent admission and derivative/reaction/coordinate causes. Fixed-source system-pressure derivatives are analytic zeros; source-pressure derivatives remain calculated and tested.
- [Native F8/F9 producers](/home/tnnrpolley21/.t3/worktrees/ePC-SAFT/t3code-225-slice-1-build-properties-reference-6cadaf6a/engine/tests/native/properties_checks.cpp:4738) provide the stated independence: F8 uses existing pure Properties caloric and pressure-volume derivatives; F9 third-order comparisons use the independent E2/E3 chain, and fourth-order comparisons difference the separately qualified third-order method at fixed literal composition.

The [retained CSV](/home/tnnrpolley21/.t3/worktrees/ePC-SAFT/t3code-225-slice-1-build-properties-reference-6cadaf6a/analyses/2026-reactive-caloric-actions/results/comparisons.csv) confirms:

| Check | Retained result |
|---|---:|
| F8 pure identities | 48; maximum error/target **0.22121550084594097** |
| Fixed-source system-P zeros | **13 exact zeros** |
| F9 resolved comparisons | 204; maximum error/target **0.4029424761033838** |
| Consecutive scaled derivative changes | 68; maximum **0.0028632302974229192 J/(mol K)** |
| Doubled-terminal changes | 34; maximum **1.3184593974450377e-5 J/(mol K)** |
| Completed Cp/Cv changes | 8; maximum **1.0125233984581428e-10 J/(mol K)** |

All final comparisons pass their declared criteria. The 36 unresolved coarse observations remain identified as refinement rows. The arithmetic diagnoses support pressure-stiffness balancing and balanced independent P4 increments; no tolerance was weakened.

The [notebook](/home/tnnrpolley21/.t3/worktrees/ePC-SAFT/t3code-225-slice-1-build-properties-reference-6cadaf6a/analyses/2026-reactive-caloric-actions/notebook.qmd) accurately reports the retained numbers. Its HTML agrees with the source, its rendered SVG is identical, and the supplied PNG preview clearly presents the derivative errors. The anchored completion uses the existing equilibrium driver only for the solved center composition and correctly limits its claim to ladder-change propagation.

This adds the bounded higher-order reference calculation while retaining the lower-slope owner and E3/E4/P4 calculations still consumed by current Equilibrium. Equilibrium files are untouched, and Properties has no upward runtime dependency. Independent line counting confirms **4,509 Properties lines**, **237 added / 2 removed**, below the **4,600** hard ceiling; aggregate use is **20,368**.

Source and wheel hashes match provenance. Existing logs confirm `check-all` passed with 113 Python tests and six native programs, and `check-affected HEAD` additionally passed equilibrium, boundary and parameter-fit programs. **I performed no build, installation or native replay, and changed no files, Git state, issue content or environment.** The branch remains clean at the reviewed commit.

Qualification is **numerical verification**, limited to the retained states and EOS coordinates. It supplies no universal infinite-dilution error bound, independent ionic-ratio-path evidence, global reacting stability, physical validation, or runtime completed caloric/action delivery. Slices 2–5 remain outside this verdict. No optional suggestion or scope expansion is a completion prerequisite; this review is complete.


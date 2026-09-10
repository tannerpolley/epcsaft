# Quarto notebook hard-cutover record

Date: 2026-09-09. Scope: close out the read-only Quarto notebook boundary.
No Engine equations, algorithms, parameters, numerical criteria, or scientific
milestone status changed. The prior custom presentation/execution surface was
7,022 physical lines. The replacement notebook-specific surface is 623 lines;
including the 59-line ordinary pytest reporter, the touched replacement surface
is 682 lines.

Independent review found routing, report ownership, evidence parsing, freshness,
CSV encoding, native case attribution, comparison validation, and source-identity gaps. The corrections
were independently rechecked and accepted with no remaining material finding.

## Executed evidence

| Check | Evidence result |
|---|---|
| Installed gate | Final `bash tools/check-all.sh` passed in 9.06 seconds after preparing and installing the non-editable wheel: 9 pytest cases and four native suites passed, with ruff, mypy, analysis-layout, packet, and architecture checks also passing. |
| Wheel identity | `build/environment-wheel/preparation.json` records build fingerprint `edfa4eed528b4aad52ab7f0b0a678cb781ed138fcb9eba5d842ba0400c25b596` and wheel SHA-256 `decbb3c191dd4a3406469aacc7bed91cea42657553f130836c0e522c2297bf1a`; the wheel's computed digest matches. |
| Direct numerical falsifiers | The direct F1 executable passed 11 complete comparison rows. The direct F5 executable retained its expected failure (`regular_ok && pinned_ok`); its weighted-curvature disagreement remains unresolved. |
| Quarto watch boundary | The active preview re-rendered changes to the implementation plan, `ARCHITECTURE.yaml`, and stable result resources (`pytest.xml`, `native.xml`, `f1.log`, `f5.log`) after temporary edits were restored. These are Quarto resource triggers, not test or calculation triggers. |
| Render process boundary | An `execve` trace of `quarto render --no-execute` showed only Quarto's Sass, Pandoc, Git metadata, and standard helper processes; it launched no pytest, Meson, native check, Engine build, or scientific calculation. |
| In-app inspection | The rendered Map, Plan, and Results pages were inspected in the Codex in-app browser. Navigation, current-owner links, retained F1/F5 details, stale/failure wording, and downloadable tables were present. |

## Correction closeout

Native runs now own unique Meson logs and atomically publish either that
invocation's JUnit report or an explicit error report. The reader exposes current
pytest/native XML, preparation provenance, representative comparisons, and valid
16-column per-report CSVs; structured CSV fields round-trip to the raw XML.
Newer ordinary results are `Unverified`, older results are `Outdated by
modification time`, and neither label implies qualification. Missing, invalid,
or nonfinite comparison operands emit `Evidence unavailable`. Source identity is
re-enumerated before publication, and affected checks route Engine data and
subproject dependency changes through numerical tests.

Historical F1/F5 records remain separate and outdated. No scientific milestone
was qualified; bounded C1 design is the next queue item.

The final reader review verified all 1,343 native comparison rows across 30 case IDs,
including `density-monotonic/0` under `native.contributions.p2`, and Python
testcase fallback. The installed gate above was rerun in the attached checkout;
the later case-attribution correction affects presentation only.

## Limits

This record verifies the presentation boundary and preservation of existing
outcomes, not thermodynamic correctness or milestone qualification. F1 retains
its stated input limits; F5's weighted-Hessian failure remains a scientific
defect. The watch check does not prove every filesystem event or reboot, and the
`execve` trace covers one inspected render. The notebook remains read-only: it
does not collect, execute, retry, rebuild, or re-certify tests.

[Quarto-native view hard-cutover completion cells](../implementation-plan.md#test-system-preparation)
link this dated evidence to tooling completion. No scientific milestone is
qualified by this record.

# Belov and Aristova test 1: source and scope

DOI 10.5541/ijot.1317496. Zotero parent `5HRFS379`, article PDF `X8VYH7WX`, existing transcription `FDMTPLAW` (SHA-256 `501d9bdb5dfd89cf5584e050cfcd9cc6580fac5395f809cb91d16d81d1bd1f38`). Section 7 links the authors' Julia examples and distinguishes homogeneous gas and heterogeneous calculations.

`test1.dat` was read from `gvbelov/Heterogeneous-Equilibrium`, commit `c74b87545a3418262e8e38a7c7a2e31e1b12966a`, path `test1.dat`, Git blob `5490caab3f7a3f8cdea0d1cb883cce4323902657`. Its fetched, CRLF byte SHA-256 is `2f24904e5db64e96aa2ee96c53bd80863c9518b6c2ddd914e266181fe7459f84`, matching the original source declaration in `source.json`. Git normalizes this text to LF on storage; the stored LF SHA-256 is `8d48f0581f8e181e988b5991e48e772f8a72b0885ebe500a0fff4233da41ecaf`. The two identities describe the same source lines with different line endings. Lines 2–11 give ten species and their standard Gibbs energies divided by RT and carbon/oxygen contents; lines 13–14 give totals 2 and 1 mol; the final lines specify 0.1 MPa and 2000 K.

The published input includes graphite and diamond. `source.json` deliberately omits both and declares an eight-species homogeneous gas restriction. Its expected amounts come from a previously retained 80-decimal element-potential solve, not from a table in the article. Its original amount comparison is therefore numerical verification against that independent restricted calculation, not a reproduction of the heterogeneous published equilibrium. The historical result has no current execution identity and is not current acceptance evidence.

At inspected Engine base `4c0c0c0ea`, `engine/native/surface/parameter_builder.cpp` admits only `pc-saft` as the base model. Exact ideal-gas admission without invented SAFT parameters waits on Engine #766, PR #786 (open on inspection). Do not approximate ideal gas with small diameters or a dilute limit. Once available, the restricted comparison should preserve all eight positive amounts, use the declared absolute logarithmic-amount tolerance 2e-6, material-balance tolerance 1e-10 mol and relative ideal pressure tolerance 1e-10. Proposed checks require independent readiness before adding a producer.

For the full heterogeneous source answer, recover the authors' computed reference and the solid standard-pressure/volume convention; `test1.dat` alone does not supply molar solid volumes required by the current Engine `Solid` declaration. Do not invent them or call the restricted gas result the full published test. No new numerical result is claimed here.

Upstream update: PR #786 merged on 2026-10-22, and remote main now admits `ideal-gas` without SAFT records and derives reaction targets from species thermochemistry. This removes the upstream #766 capability prerequisite, but does not execute this gas-only benchmark or settle the full heterogeneous source inputs. The task checkout now includes that merged source and has its own prepared editable Engine. Independent readiness is still required before any gas-restriction producer.

A reduced eight-gas numerical-verification design is retained in MAIN at
`build/issue-789/2023-belov/design-r1.md`. It uses the existing independent
restricted reference and explicitly proposes a one-state source-Gibbs encoding,
without asserting source calorics or a heterogeneous result. The article
Section 6 identifies the relevant IVTANTHERMO standard pressure as 1 atm;
the separate NIST 1 bar convention is not substituted. No producer has run.

# Generic Parameter Fit Cutover Execution Ledger

Status: active

Canonical design:
[`2026-08-06-generic-parameter-fit-hard-cutover.md`](2026-08-06-generic-parameter-fit-hard-cutover.md)

## Goal

Implement every phase 0--10 of the canonical hard-cutover plan in the current
monorepo worktree. A Luna 5.6 Max worker performs each bounded phase packet. The
main task integrates and validates the result. A fresh Sol 5.6 High reviewer
must then approve that phase before the next phase begins.

## Attached checkout

- Repository root:
  `~/Workspaces/Engineering/ePC-SAFT-project`
- Branch: `codex/reduce-tests-incremental-wheel`
- Baseline HEAD: `1f06ef488c41a732840b3469f5444205e26057fa`
- Remote: `https://github.com/tannerpolley/ePC-SAFT.git`
- Worktrees: one attached canonical worktree
- Existing state: broad pre-existing monorepo, test-reduction, Equilibrium,
  wheel-routing, Validation, and documentation changes; preserve them unless a
  phase explicitly classifies and supersedes an overlapping Regression change.

## Phase gate

For each phase:

1. Main task freezes the phase packet and exact write scope.
2. One Luna 5.6 Max worker performs the packet and reports changed files,
   validation, and unresolved risks.
3. Main task inspects and integrates the phase result, resolving overlap with
   the attached worktree without reverting unrelated changes.
4. Main task runs the phase's narrowest authoritative checks.
5. One fresh Sol 5.6 High reviewer audits the stable phase result against the
   canonical plan, current files, and validation evidence.
6. Reviewer findings are remediated and re-reviewed as necessary.
7. Main task records approval evidence here and only then starts the next
   phase.

No phase is complete because code exists or a narrow test passes. Its exit gate
must be proved against the current worktree. Final completion additionally
requires the full Phase 10 acceptance audit.

## Integration policy

- The main task is the only Git/index writer.
- Workers receive disjoint or explicitly bounded write scopes.
- Reviewers are read-only and fresh for every phase.
- No compatibility aliases, dual public routes, legacy readers, source-specific
  runtime paths, or unrelated cleanup may be introduced.
- Scientific claims require numerical, physical, and source evidence at the
  owner named by Governance NVDD.
- Literature breadth belongs in installed-wheel Validation; compact exact
  derivative and solver invariants belong in Engine.
- Historical plans and receipts remain provenance and are not normalized.

## Tooling warning

The repository profile resolver identifies the root as `scientific-computing`.
The current Chemical Engineer `health_check.py` still expects the retired
`ePC-SAFT-project/ePC-SAFT-engine` checkout and rejects the valid monorepo root.
Until that external skill is updated, use live repository ownership, equation
registry checks, root affected/all checks, installed-wheel verification, and
focused source tracing as the authoritative route.

## Phase state

| Phase | Scope | Worker | Reviewer | Status | Evidence |
|---:|---|---|---|---|---|
| 0 | Freeze contract and machine-readable inventory | Sagan `019fd982-8c28-76a1-8c27-e72e84f8a3dc` | Kuhn `019fd9e8-fef9-74d3-9d84-3538cb097880` | approved | 615 unique entries; JSON shape and summaries, all-entry Git state, 98-file locked-packet closure, Gross 2001/2002 campaign closure, and live API/native/build searches passed. Fresh Sol approval followed five bounded rejection/remediation cycles. |
| 1 | Canonical Parameters and topology | Halley `019fd9ee-1090-7600-9d0f-1d52be81939c` plus primary integration | Confucius `019fda67-9799-7df3-b2f3-13a0dcbdda3f` | approved | Strict source-neutral mapping/JSON document, fixed/selectable separation, model-family metadata, `site_role`, presets, explicit direct/combining topology, noncontinuous derived outputs, exact derivative-source records, JSON-only loading, and old-loader removal integrated. All 11 canonical documents construct; 27 focused tests pass with 5 deselected; packet/preflight/artifact, Ruff, mypy, compile, wheel, diff, and three adversarial combining-rule smokes pass. Fresh Sol verdict: `ship`, no findings or residual core risk. |
| 2 | EOS object and pure-saturation ownership | primary integration | Singer `019fda82-204b-7f02-99b0-08556df29af3` | approved | Strict `Pure`/`Mixture` share `State`; `Pure.saturation` and its exact NLP derivatives call `NativeModel` directly from EOS-owned sources; the Equilibrium source, bindings, exports, fingerprint allowlist, and provider bridge are removed. The installed wheel builds; 31 focused numerical tests, Ruff, mypy, clang-format, YAML, diff, and runtime ownership checks pass. Full Engine status is 70 passed, 9 deselected, with 10 pending-Regression failures confined to deleted `Parameters.from_dictionary`. Fresh Sol verdict: `ship`, no findings or Phase 2 residual risk. |
| 3 | ActiveParameterSet and exact derivatives | Noether `019fda89-7975-78a2-890d-614b8d3d9ed9` plus primary integration | Planck `019fdabe-e8e0-77d0-8423-1deb11d8a0fb` | approved | EOS-owned immutable `ActiveParameterSet` and capability/result records; source-neutral native identity bridge; exact base, pair, association, polar-with-fixed-moment, Born, permittivity, solvation, and suppression derivatives; exact implicit density, association, fugacity, and pure-saturation sensitivities; generalized association coverage for pure, neutral-mixture, electrolyte, and temperature-dependent-diameter contexts. The installed wheel, native CMake build, equation registry, Ruff, mypy, 14-case derivative matrix, and all 45 EOS numerical tests pass. Full Engine status before the planned Regression deletion is 83 passed, 9 deselected, with 10 displaced old-Regression failures confined to deleted `Parameters.from_dictionary`. Fresh Sol verdict: `ship`; no findings. Residual: the bypassed legacy SDK binding remains solely for Phase 8 deletion. |
| 4 | ParameterFit schema, sources, and resolution | Pascal `019fdad6-db4d-7ed2-bef2-7b211a2d46ac` plus primary integration | Aristotle `019fdb6f-c456-72c1-9a0c-fb4f804f2605` | approved | Immutable strict `ParameterFit`, one registry-generated JSON Schema, exact selector/tie/coordinate resolution, source fingerprints, canonical unit conversion, typed table and correlation sources, EOS-domain containment, and converged deterministic quadrature over both objective and exact Jacobian. Exact pressure and molar/mass-density derivatives, all admitted segment-diameter correlation coefficients, Ruff, strict mypy, native/wheel builds, 45 EOS numerical tests, and adversarial schema/unit/domain/overflow smokes pass. Ten full-Engine failures remain only in the Phase 8 deletion set. Fresh Sol verdict: `ship`; no findings or residual risk. |
| 5 | Generic residual engine and Ceres execution | Pasteur `019fdb7b-4ff4-7410-a95d-40edf14991af` (no patch) plus primary integration | Euler `019fdb9c-c901-71c0-9c1d-1a1733b666a7` | approved | One callback-free native pure-saturation residual owner—the only currently admitted Observation Family—binds directly to EOS candidate values and exact sensitivities. It uses scalar coordinate blocks, per-state vector residual blocks, log-ratio and scaled-difference semantics, explicit multipliers, deterministic Halton multistart, single-threaded Ceres, safely serialized Ipopt/MUMPS evaluations, non-evaluable backtracking, complete start retention, and exact SVD/rank/conditioning/bound/basin diagnostics. The installed wheel, Ruff, strict source mypy, 47 focused numerical tests, and a full Engine run of 86 passed/9 deselected with only the 10 planned Phase 8 legacy-test failures pass. Fresh Sol verdict: `ship`; no findings. Residual: delete the old source-specific extension in Phase 8. |
| 6 | FitResult and transfer | Euclid `019fdba1-5319-7562-a5d5-80009e78c271` (no patch) plus primary integration | Galileo `019fdbb2-2d85-7551-aab8-71327316a41f` | approved | New immutable `FitResult`/`FitStart` retain the complete resolved fit, full input Parameters document, row identities, every start and failed evaluation, exact physical/affine Jacobians, block contributions, SVD/rank/conditioning/bounds/basins, deterministic accepted-start selection, and artifact/package fingerprints. Canonical mapping/hash and explicit JSON export normalize nonfinite failed-start evidence to JSON null. `Parameters.to_mapping()` and unit-aware `with_values()` create a new fitted document with fit provenance; input values and files remain untouched. The installed wheel, Ruff, mypy over 24 source files, and three numerical fit/result tests pass. Fresh Sol first rejected incomplete nonfinite-vector normalization; remediation and targeted re-review returned `ship`, no findings or residual risks. |
| 7 | Migrate required scientific slices | Beauvoir `019fdbb8-dcd5-7d20-896b-724d4f96a8b7` (classification, no patch), Mencius `019fdbc3-435c-73d0-b1a9-c65d610ba9f3` (read-only audit), plus primary integration | Hypatia `019fdbdd-b5b2-72a3-a2cd-70daec6d20da` | approved | The same native `EosCost` now executes Pure Saturation State and EOS-owned fixed-state pressure observations. Constant neutral binary `k_ij`/`l_ij` use the exact CppAD volume--parameter Hessian derivative `-RT H[V,p]`, strict positive normalized composition, one canonical pair coordinate, and fail-closed charge/cardinality admission. A source-neutral methane/ethane slice recovers `k_ij=0.03` from four fixed states and matches independently rebuilt central differences; the compact EOS matrix covers base, pair, association, polar-fixed-moment, and electrolyte fixed-topology derivatives. Correlated pair coefficients, direct permittivity calibration, and specialized Born/water/aqueous routes are explicitly unsupported and proceed to deletion rather than wrappers. The installed wheel and origin check, Ruff, strict mypy over 24 source files, 49 affected numerical tests, and invalid-state smokes pass. Full Engine status is 88 passed, 9 deselected, with exactly 10 displaced old-Regression failures scheduled for Phase 8. Luna found and remediation closed zero-fraction, nonpositive-pressure, and charged-pair gaps. Fresh Sol verdict: `ship`; Phase 8 may proceed. |
| 8 | Delete displaced implementation | Zeno `019fdbdf-35f0-7003-866c-9e05b59da6b8` plus Gibbs `019fdbfb-26d1-7573-9582-84d66cb03ee2` remediation and primary integration | Helmholtz `019fdc32-c742-7d10-bdbd-da56e9bcfb31` | approved | Deleted the complete source-specific Regression lifecycle: public object graphs and exports, packaged datasets, native fit owners and extension, Regression-facing SDK/header/callback route, obsolete tests, CMake sources, and live specialized docs. Replaced the public SDK with one private in-wheel `native_transport` used only by current Engine modules and linked generic `fit.cpp` directly into `_core`; removed the final unused pure-parameter exports after review. The rebuilt installed wheel and origin verifier, source/wheel/binary/dead-export audits, equation registry, Ruff, strict mypy, shell syntax, `git diff --check`, all 88 retained Engine tests (5 deselected), and a fresh isolated Clang ASan/UBSan run of 11 tests in 270.76 s pass. Cleanup audit is clean. The bounded Regression cutover removes 37,895 lines while preserving historical plans, receipts, results, and source evidence as provenance. Fresh Sol verdict: `SHIP`; no Phase 8 blocker remains. |
| 9 | Scientific evidence and installed-wheel Validation | Leibniz `019fdc3a-062b-7232-aa76-22d9275da600` (interrupted without patch) plus primary integration | Dewey `019fdc5c-a6f0-7521-a286-9fd0e0971dbc` | approved | Replaced the legacy campaign object graph with one installed-wheel `fit` campaign using three independent Glos 2004 propane saturation states and simultaneous `epsilon/k`, `m`, and `sigma` coordinates. All four low-discrepancy starts are usable in one basin; the accepted solution is rank 3, interior, and has maximum absolute log residual about 0.0021 against the declared 0.005 engineering gate. A deterministic GCC Release crash under repeated Ceres saturation evaluations was traced to per-call Ipopt lifecycle corruption and removed by consolidating pure saturation on one bounded exact-Jacobian Newton solver for empty and nonempty association topology. Public saturation anchors, exact sensitivities, all 88 Engine tests (5 deselected), the installed-wheel campaign, retained Figiel source-value evidence, Ruff, strict mypy, JSON/TOML/hash checks, wheel build/origin verification, and the four-start Release repro pass. Unsupported May VLE and specialized Figiel Regression executables and the obsolete campaign adapter/copies were deleted while source evidence and historical receipts remain provenance. Fresh Sol verdict: `SHIP`; residual association, mixture-VLE, permittivity/electrolyte, predictive, and globality evidence is explicitly deferred outside the approved EOS-only boundary. |
| 10 | Atomic public cutover and completion audit | McClintock `019fdc62-344a-70d3-81f6-7000b0f0ee54` (no patch) plus primary integration | Newton `019fdc70-940e-72a0-a03b-e1f0b12db123` | approved | The public Regression surface is exactly `ParameterFit`, `FitResult`, and `fit`; the quickstart, Python API, Engine README, and architecture describe one generic EOS-only route, immutable parameter transfer, installed-wheel Validation ownership, and explicit electrolyte/phase-equilibrium deferrals. The obsolete Provider Slice 1 executable and live Perdomo split-wheel Regression instructions were removed while receipts and source evidence remain historical provenance. The retained FitResult test now asserts only numerical parameter transfer without input mutation. The rebuilt installed wheel and origin audit, packet validation, all 88 normal and 5 slow Engine tests, focused generic-fit and installed-wheel propane campaign evidence, Ruff, strict mypy, Clang/clang-tidy, fresh ASan/UBSan, YAML/JSON/TOML/hash and dead-surface audits, preflight, artifact-store, diff, and cleanup checks pass. Root affected/full verification stops only at the two explicitly deferred Khudaida electrolyte-flash cases (`native transport is incomplete`). Fresh Sol initially returned `RETHINK` for the three stale surfaces; direct remediation and targeted re-review returned `SHIP`, with no blockers. |

## Execution scope decision — 2026-08-06

The user directed the implementation to prioritize the total Regression rewrite
and hard cutoff. Pre-monorepo Khudaida data-lock compatibility and the legacy
general EOS literature-handoff runner no longer gate Phases 1--10. Completed
results and receipts remain provenance. Phase 9 will create or rewrite only the
minimal installed-wheel campaigns required to prove the new generic Regression
workflow; broader phase-equilibrium campaign maintenance is a separate
follow-up.

## Phase 0 packet

Objective: classify every current Regression-facing runtime, native, data,
test, documentation, Validation, build, and public-export surface before any
new implementation begins.

Required artifact:
`2026-08-06-generic-parameter-fit-inventory.json`

Every entry must record:

- path;
- category;
- current role;
- public/native/build/test/data/documentation surface;
- source-specific or source-neutral identity;
- current Git state;
- disposition: retain, migrate, replace, delete, preserve-as-provenance, or
  unrelated;
- target owner/path when applicable;
- plan phase;
- overlap risk with pre-existing work; and
- evidence supporting the classification.

Phase 0 exit gate: every current surface is classified, the inventory validates
as JSON, no public or native owner is omitted, and a fresh Sol reviewer approves
the inventory against the live tree and canonical plan.

# Independent input and output measurement errors

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. Both campaign
modes refuse changed retained truth on the current Engine; the original normal draws,
truth, tolerances and fit endpoints remain frozen.

This synthetic experiment asks whether the augmented native Ceres estimator reduces signed parameter bias in both the frozen MEA and water/Na/Cl cases. Issue [#510](https://github.com/tannerpolley/ePC-SAFT/issues/510) and [the frozen rules](inputs/campaign-freeze.md) own the method and decision. The supplied MEA parameters are synthetic truth here; this is numerical verification, with no claim of physical validation.

The retained results are the requalification on the current Engine (issue [#624](https://github.com/tannerpolley/ePC-SAFT/issues/624)). After #594 made the ionic reference the √m limit of the 1e-12/1e-14 mol/kg rungs, `prepare` regenerated the truth, sigmas, noisy outputs, `verification.csv` and `rho.csv` with the retained normal draws and bootstrap streams, and both estimators were refitted under the unchanged rules. The comparison rule C1–C7 was frozen in `inputs/campaign-freeze.md` before any fit; `results/requalification.csv` records old (2026-10-14, `3fbccd954`), new and Δ for every compared quantity. C1–C6 hold, and the decisions are unchanged: MEA accepted, ionic inconclusive. C7 `pressure`: NOT MET. The ionic density-root pressure component rose by a ratio of 2.334 for fixed inputs (3.802e-9 → 8.874e-9) and 2.383 for input and output errors (1.390e-6 → 3.313e-6), against the frozen band [0.5, 2]. Cause: the γ± reference now charges the 1e-14 mol/kg root with weight 10/9; at 298.15 K and 1e5 Pa that root stops at a pressure residual of −3.39e-5 Pa (2.5e-13 reduced, inside the 1e-12 density-root tolerance), against −1.39e-6 Pa for the 1e-12 mol/kg root, so the γ± reference pressure term rose about 27×. The driver builds the 1e-14 mol/kg composition exactly as the Engine does. The term is about 3e-5 of the ionic numerical allowance, so it cannot affect the verdict. The band was a prediction for a term the design marked unmeasured; it is not changed after the fact, and the miss is recorded, not waived. The requalification therefore does not pass every frozen check. The after-the-fact note is a separate field, `outcome_note_after_fit`, in `inputs/campaign-freeze.json`; the frozen rule text is unchanged.

The first requalification `run` fitted nothing. Since #486, `FitControls` begins with the six-entry `solver_options` array, so the driver's positional braces filled tolerances and left the iteration and elapsed-time limits at 0; every fit was refused before its first iteration. The driver now assigns the frozen controls by name, and the inputs were regenerated and refrozen before any fit. The refused records are in the main checkout's `build/issue-624/build-r1/failed-run-results.tar.zst`; none was read as a result. The completed run took 514 s wall and 3 h 48 min CPU; the cgroup `memory.peak` was 8.99 GB (systemd's unit summary reported 5.9 G for the same unit); the longest fit took 122.5 s (MEA).

`scripts/campaign.cpp` owns native truth generation, normal draws, both estimators, endpoint error propagation, the four tighter replays and paired bootstrap statistics. Python mechanically declares the existing MEA inputs and reads retained results. All EOS, equilibrium, derivative and fit calculations run in C++.

[`scripts/diagnosis/`](scripts/diagnosis/) retains the native reference-ladder probes and small numerical tables; their original output directory was `build/ionic-reference-allowance/`, while `scripts/diagnosis/present.py` only reads retained tables for the notebook. The probes reproduce the retained tables at main `4ff4e4f5b` (#510), before the √m-limit reference (#594, merged `44e6b416`). The patch still applies through `2666b5296`, but it fails against the oracle rewritten by #623 (`35fd2c9f`) and is not ported.

The historical declaration source is `../2026-fitting-state-reuse/inputs/mea-initial.json`; its SHA-256, component order and first sixteen state identities are retained in `inputs/provenance.json`. `scripts/prepare.py` generates the corresponding native declaration header under `build/issue-510/`. The campaign’s preparation and fit commands are retired from this record.

The runner uses niceness 10 and single-thread numerical libraries alongside process-level replicate parallelism. Existing `draws.csv` files are read; they are never regenerated, clipped or redrawn. The bootstrap index streams are likewise retained. `statistics` recomputes only native statistics from retained fit endpoints. `prepare` verifies exact columns against rebuilt-state differences at the frozen steps; the ionic cancellation probe uses the existing binary128 contribution-kernel verification facility.

Each case has `inputs/<case>/truth.csv`, `draws.csv`, `verification.csv`, `rho.csv` and `bootstrap-indices.csv`. Git retains each case's `summary.csv`, `decision.csv`, `numerical-components.csv`, `paired-errors.csv` and `bootstrap-means.csv`, together with one `endpoints.csv`. The latter has 264 rows: 256 main fits and eight tighter replays. Its first three columns identify replicate, method (`fixed` or `eiv`) and run (`main` or `replay`); the remaining columns preserve the original endpoint fields exactly: regular-fit flag, signed optimizer-scaled error and nine propagated error components in the order of `numerical-components.csv`. No fit or replay is omitted. `results/requalification.csv` holds the C1–C7 comparison with the 2026-10-14 round.

These are frozen finite-ladder results. After #823, analytic references have no measured terminal change; the current allowance reader refuses that missing estimate. Reproducing the historical numerical allowance requires its recorded Engine and campaign source, rather than assigning zero error to the new limit. No frozen experiment was rerun for #823.

The row-level prediction error table separates pressure roots, reference values, fit-own/cold discrepancies and rounding. The action table separates checked linear solves, reference actions, fit-own/cold discrepancies, rounding and the frozen verification estimate. The same-point cold evaluation is performed at every true state and endpoint. MEA consumes laddered neutral reaction references. Its hidden reference-trace root residuals are represented by terminal changes and observed discrepancies, without a separate invented term. The unestimated cross term involving the derivative of the output/reference sensitivity remains a limit of the first-order allowance; the declared factor of two covers it only nominally.

Signed bootstrap means and their four simultaneous percentile intervals precede the absolute-value transformation. Numerical qualification requires every main fit to be regular, interior and full rank, and every numerical term and replay to be available. The strict decision is `U_Delta < -epsilon_num` in both cases. Finite-sample bootstrap coverage is approximate. These conditions and all failures remain part of the result even if a point bias decreases.

`allowance` reconstructs each retained fitted point and its existing SVD, then evaluates the cold observation and error components natively; it never refits. The historical `allowance` and `statistics` modes read the original endpoint directories preserved in the raw archive below. The reporting corrections use the complete augmented column norm and apply the ordinary fit's relative weights to its cold columns. The latter defect was diagnosed after the first bias summary was read. Earlier values remain in the archived `before-complete-norm-*` and `before-cold-scaling-*` files and local commits; fitted values and signed-mean bootstrap results are unchanged. `scripts/calculation_ids.py` adds or updates a source-record path on tables with identical content or an existing identifier, preserving their numerical fields. Line endings are LF.

`scripts/present.py` reads the native summary, decision and numerical-component CSVs to prepare HTML tables and the signed-mean interval figure; it also tabulates `results/requalification.csv`. `notebook.qmd` includes only those retained presentations and renders with execution disabled. The matched A environment is a clean main archive: its Engine build predates the final main refresh, which changed only the accepted corpus inventory tool. The source identities and this equality are recorded under `build/issue-510/`.

The separate matched-main fixed-input verification driver and both results are retained at local commit `f3e8a53b` and under `build/issue-510/`. Its verification-only `fixed` mode is absent from the final campaign entry; no test helper is placed in an unrelated runtime or test layer to meet the analysis allowance.

## Raw record archive

The requalification's detailed records from `7af85124` are stored outside Git at `/home/tnnrpolley21/Workspaces/Engineering/ePC-SAFT/build/archive/issue-624-campaign-raw.tar.zst` in the main checkout. A [raw-file index](results/raw-manifest.tsv) records each archived file's path relative to this leaf, byte count and SHA-256, plus the archive's size, SHA-256 and source commit. The archive preserves the `results/<case>/<replicate>/<fixed|eiv>/<main|replay>/` directories: predictions, residuals, exact columns, singular vectors, covariance, iteration costs, work, timing, typed refusals, reconciliation and numerical-error diagnostics. It also contains the endpoint files and the true-state error records under `results/<case>/truth/`.

The 2026-10-14 round's archive, `build/archive/issue-510-campaign-raw.tar.zst` (records from `e9a9850e`, with earlier allowance versions and the entire partial `results/attempt-1/` round), and its raw-file index are retained at `3fbccd954`.

To verify and restore those paths, run from the repository root:

```bash
leaf=analyses/2026-errors-in-variables
archive="$(git rev-parse --show-toplevel)/build/archive/issue-624-campaign-raw.tar.zst"
printf '%s  %s\n' "$(awk -F '\t' '$1 == "# archive_sha256" {print $2}' "$leaf/results/raw-manifest.tsv")" "$archive" | sha256sum --check
nice -n 10 zstd -t "$archive"
nice -n 10 tar --zstd -xf "$archive" -C "$leaf"
```

The repository's `.gitignore` keeps restored raw directories outside Git. Rendering `notebook.qmd` needs only the retained summaries, figure and six HTML table includes; it does not need this archive or execute the model.

## Historical campaign modes (#720)

`scripts/campaign.cpp` prepare and run are **historical, not rerunnable**.
At #696’s measured commit `335707d59c67cf1f92df12b5d25530bf31d916a6`, both
`native-2026-errors-in-variables-campaign.cpp-prepare` and
`native-2026-errors-in-variables-campaign.cpp-run` stopped with
`retained true output changed` (report lines 700–701). Prepare was compiled
with the required `EPCSAFT_EIV` definition and generated MEA declarations.
The retained-truth guard refused the current image in both modes before the
campaign proceeded. Frozen truth, normal draws and tolerances were preserved;
these failures are not completed fitting evidence. The retained #624 values
above remain evidence for their original image only.

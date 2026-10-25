Independent delivered review by Claude Opus 5.5, high effort, task
`c5d6705e-81e9-4acd-9503-ce4802071aad`. The review below is retained verbatim
as reviewer evidence; its suggested gates and optional work do not amend the
owner's accepted scope.

> # Delivered review: ePC-SAFT #221 finding 5, lever-derivative correction
> 
> **Verdict: Supported.**
> **Reviewed:** commit `09bcb72f` (branch `issue-221-finding5`) against merged main `bdde4003`. This was a read-only review. I did not build or run anything; I inspected the final source, the retained analysis files and the logs in `/tmp/finding5-correction/`.
> **Next stage:** stop (closeout). It is gated on the owner's cost acknowledgment in limit 1.
> 
> ## Result
> 
> `CompiledSystem::lever_jacobian` now gives the correct first derivatives of the equilibrium rows along the temperature, pressure, feed-amount and active-parameter directions, at a fixed solved state. The old dense assembly never differentiated an ion's packing diameter (`PackingDiameter`). Because volume at fixed packing fraction is the packing coefficient divided by the packing fraction, and that coefficient contains the ion's diameter cubed, chemical and pressure parameter columns were missing the volume response.
> 
> At the accepted aqueous electrolyte equilibrium (300 K, 100 kPa, feed 0.98/0.01/0.01 mol, central classified residual 1.76e-14), the cation packing-diameter column changes as follows. Units are per ångström.
> 
> | Row | Old | Corrected | Centered, step 1e-4 Å | Centered, step 1e-5 Å |
> |---|---:|---:|---:|---:|
> | 0 | 0.305912 | 0.0724851948 | 0.0724851949 | 0.0724851948 |
> | 1 | 7.145186 | 7.0080692539 | 7.0080692552 | 7.0080692538 |
> | 2 | 0.191408 | 0.0328172189 | 0.0328172190 | 0.0328172177 |
> | 3 (pressure) | 0.303791 | 0.0720755439 | 0.0720755440 | 0.0720755440 |
> 
> The old values were wrong by up to 83 % relative. The corrected values match both step sizes to within 1.35e-9 per Å. This is numerical verification of a derivative. It is not physical validation and does not adopt any parameter.
> 
> ## Correction mathematics and code inspection
> 
> - **Assembly.** `lever_jacobian` prepares one center. For each of the 2 + C + A unit columns (C components, A active parameters), it takes `row_series(...).rows[1]`. That function and its recorded geometry are the same ones production uses: `sensitivity()` in `implicit_sensitivity.cpp:173` already calls `row_series` the same way.
>   - The recorded geometry builds the diameter from the parameter-dependent copy of the parameters, so the ion diameter is differentiated. This matches the equation tag `packing_coefficient_parameter_derivative`.
>   - The volume chain then reaches the chemical rows through −ln v and the residual-gradient Hessian action, and the pressure rows through 1 − Pv/RT and the −v·A_v term.
> - **Free temperature or pressure.** For a free T or P, the recorded geometry takes T or P from the state vector, so the control entry is ignored. The new code instead sets a unit direction on that state coordinate.
>   - The indices match `prepare` and the recorded geometry: T at `primal_size − free_intensive_count`, P at `primal_size − 1`.
>   - This keeps the old column definition, a partial derivative at fixed other coordinates. The free-T/P, pinned, scaled-feed case agrees to ≤4.7e-15 relative in columns 0 and 1.
> - **Padding.** The point is padded to `unknown_count()` after `prepare`, so `prepare` still rejects points shorter than `primal_size` as before. The padded multiplier and electric-potential entries get a zero direction, so first-order coefficients do not depend on them. Padding is needed only to keep the recorded geometry's electric-potential and multiplier indices in range for primal-only callers such as S9's `fixed_start`.
> - **Feed scaling, pins and reactions.**
>   - Balance (moiety) rows divide by a fixed `feed_total`. A pinned phase amount is `fraction × Σ n_feed` and is differentiated. This gives −B/N_feed + fraction·B·x/N_feed, the same as the removed `pinned_amount_feed_term`.
>   - The control has only T, P, feed and parameter entries, so the reaction-target branch of `row_series` is inactive, as in the old code.
>   - Charge-balance and normalization rows stay zero.
>   - The non-reference standard state fills its pressure and parameter slopes with zeros, so chemical-row pressure and parameter columns gain no new standard-state term.
> - **Refusal, shape and order.** The reference-reaction refusal still comes first, with the same message. Width is 2 + C + A and the layout is `row × width + column`, as before. The `equilibrium_lever_jacobian` profiling scope is kept. Nested `equilibrium_row_series` scopes now appear once per column; they are profiling entries only.
> - **Removal.** Only `PackingTerms::parameter_derivative` and the block in `packing_terms` that produced it are removed.
>   - The removed neutral `SigmaConstant` branch (d·d/σ) and `DispersionEnergyOverK` branch (0.36σ/T·e^{−3ε/T}) are algebraically identical to differentiating `segment_diameter`.
>   - Ions take d = `packing_diameter_angstrom`, so excluding ions from those two families was correct. The only defect was the missing `PackingDiameter` branch.
> - **Other writer's file.** `git diff bdde4003 09bcb72f -- engine/native/equilibrium/gibbs_minimizer.cpp` is empty.
> 
> ## Falsifier and tolerance
> 
> - `boundary_checks packing-derivative` compares the four chemical and pressure rows at steps of 1e-4 Å and 1e-5 Å (8 comparisons) with an absolute tolerance of 2e-9 per Å. The largest measured error is 1.348e-9 per Å.
> - On the base build it fails at its first comparison. `base/falsifier.log` records model 0.30591249238 against reference 0.07248519485 (error 0.2334), verdict false. On the candidate, all eight pass.
> - The stencil holds the state coordinates fixed (log x, log η, amounts), so it tests exactly the fixed-packing-fraction derivative the lever defines.
> - The defect is about 1e8 times larger than the tolerance, so the test separates "chain missing" from rounding cleanly.
> 
> ## Evidence and counts, checked against the retained files
> 
> - **Matrix comparison.** `columns.csv` has 230 entries (56 + 126 + 48). Recomputing it:
>   - Exactly 4 entries differ by more than 1e-12 relative: the ionic case, column 5, rows 0–3.
>   - Of the other differences, 12 are signed zeros and 18 are non-zero, with absolute difference ≤2.22e-16 and relative difference ≤4.7e-15.
> - **Native records.** `comparison.json` matches `summary.json`: 6267 records (128 + 4865 + 675 + 520 + 79), 5679 judged, no false judgments, 7 numeric differences, maximum 1.96e-16 relative.
> - **Khudaida slice.** I re-ran `cmp` on base against candidate for `stage-results.csv`, `engine-residuals.csv` and `qualification-comparison.csv`. All three are byte-identical and their SHA-256 hashes match `summary.json`.
> - **Python tests.** `candidate.log` and `affected.log` show 113 passed, 1 deselected, and "All checks passed". The parametrized `test_referenced_continuation_follows_the_cold_solves` (two pressure cases) is not marked slow or expensive, so it is in the default collection. The logs record only the count, not test names.
> - **Reproduction.** `compare.cpp` calls the compiled `lever_jacobian`; it contains no duplicate assembly and no stencil code.
> - **Line budget.** Executable code is 91 lines added and 92 removed. That is 74 analysis-script lines plus 17 runtime lines added, against 92 runtime lines removed, giving −75 net in Equilibrium (6129 → 6054). The Equilibrium hard ceiling is lowered from 7650 to 7575, which keeps the 1521-line shared roadmap allowance. Lowering a ceiling needs no approval.
> - **Lithium.** This was a source inspection only. Neither inspected reference builder sets active parameters, so no packing-diameter column is exposed there. This is stated as limited to those sources.
> 
> ## Limits and open items (not required corrections)
> 
> 1. **Cost is above the band the owner accepted.**
>    - The owner accepted a 1.5–1.9× slowdown, from the interleaved study (1.53–1.91×).
>    - The final measurement used separately installed builds on a shared, changing host load. Its ratios of median times are 1.92×, 2.35× and 2.25×, and single repetitions reach 3.27×.
>    - The notebook discloses this honestly. However, the `code-budgets.md` amendment cites only the accepted 1.53–1.91× range.
>    - The absolute cost is about 0.1–0.4 ms per call, and Lithium's live-derivative path calls `lever_jacobian`.
>    - The owner should confirm that the higher measured ratio is accepted.
> 2. **Not all parameter families were compared before and after.** The 230-entry matrix covers Kij, segment count and ion packing diameter, plus the T, P and feed columns. Neutral `SigmaConstant` and `DispersionEnergyOverK` lever columns were not compared. Their removed formulas are algebraically identical to the shared diameter chain, which I checked by inspection, not by numbers.
> 3. **The tolerance is empirical.** 2e-9 is 1.48× the largest measured error. Rounding error at the 1e-5 Å step reaches 1.12e-9 on row 2, so a compiler or platform change could narrow that margin. The separation from the defect is decisive regardless.
> 
> ## Optional suggestions (not work orders)
> 
> - A step-scaled or relative stencil tolerance would be more robust than a fixed 2e-9.
> - Recording the geometry once and sweeping all directions forward would recover much of the cost. The per-column recording is the cost source the notebook names.
> - The notebook calls the 18 non-zero differences "last-bit", but one is about 24 units in the last place (4.7e-15 relative). The absolute ≤2.22e-16 statement is correct; "rounding-level" would be more accurate wording.
> - Branch commit `74f1370e` renames "Regression" to "Parameter regression" in past budget amendments. It is a separate, wording-only commit, outside the scientific slice.

# Damped-trial diagnosis and negative design result

No tested predicate separates the rejected ionic/Khudaida trials from all required productive trials while meeting the ionic work limit. Disposition proposed to the accountable lead: **optimization not adopted**. This is a completed diagnosis/design handoff, not a new solver implementation or a claim that all possible predicates are impossible.

## Conditions and recording verification

Unchanged polish source: `ec0bafa71b5022821f1774b2b9fffb06622be87c`. Assigned branch after authorized main merge: `cbef363598704af19ece9146ec5aae97cf11a16e`. Pristine main wheel SHA-256: `f5b38aaadb19f2339c8215e7a321f373f03e2940dc1d1ce12cabc5e843ed48a2`. Recording-only wheel SHA-256: `44a42e3a648ffa944cdbb3d89c95413e7f62a565f2d814624d3eb7063e644adc`.

The recording compares 15 returned states exactly: all coordinates, amounts, fractions, densities, intensive values, residuals, multipliers and work match. All four existing comparison checks pass with exactly identical retained rows. The replay retains 220 trials. It runs main's original Newton solve, damping and admission; no shortcut is enabled. Engine and fitting calculations stay native C++.

Cold ionic LLE and the sequential eight-point Khudaida replay use the declared water/ethanol/isobutanol/Na+/Cl− packet at 293.15 K and 100000 Pa. The four protective check families reuse the exact current declarations in the trace-isomerization, native action transport and parameter-fit tests. Masking and small-phase declarations reuse the prior retained inputs. Input-owner hashes are in `provenance.json`.

## Numerical question

For the rounded trial displacement h = x_trial − x_current, predict raw rows as r_pred = r + J h. The native recorder accumulates this in long double. Predicted merit decrease is max|r| − max|r_pred|. The assembly rounding estimate for each row is u Σ|t|, u = 2^-53, with each t a consumed assembly operand. Norm and active-row comparisons are separate. The native association-gradient/evaluation floor is retained separately. Neither estimate is a rigorous total thermodynamic error bound.

Coordinates are compared by distance between ordered binary64 representations. Both the number of changed coordinates and the smallest/largest changed distance are retained, as are identical-current and identical-previous-trial comparisons. All predicate values depend only on the current prepared state, Jacobian, computed direction and representable trial, before preparing/evaluating that trial.

## Cause-level findings

All five rejected cold-ionic trials and all thirteen rejected sequential-Khudaida trials are inside the physical composition cutoff, in nominal polish steps before the extension guard. Each fails strict raw contraction; all dominant rejected rows are chemical. Seventeen dominate at row 2 (aqueous isobutanol), and one Khudaida trial dominates at row 9 (organic sodium). Every one has a positive linear predicted merit decrease. Thus the rejection follows recomputation of the thermodynamic rows, rather than a duplicate representable coordinate or predicted absence of descent.

For cold ionic step 1, full trial, the controlling chemical row is predicted near 4.20e-15 but recomputes to 1.98e-13. Its own assembly estimate is 1.17e-15. A final assembly estimate therefore underestimates this discrepancy. The maximum native association estimate is about 2.15e-10, which is too broad to distinguish productive polish steps. These comparisons do not isolate a particular association or electrolyte operation as the cause.

A bounded precision alternative was checked: extended-precision summation of the already evaluated operands still rejects all eighteen ionic/Khudaida rejected trials. It also reverses the accepted ionic quarter-step merit comparison. This does not recompute the equation of state in higher precision and does not justify changing the accepted merit or summation method.

## Predicate verdict

Exact-current identity, identity with the preceding trial, and all coordinates within one ulp catch zero rejected ionic and Khudaida trials. They cannot meet the work target. The rejected ionic trials change up to 81–325 ulps; rejected Khudaida trials change up to 36–702 ulps.

Predicted decrease below the norm assembly estimate catches three rejected ionic trials but also its productive quarter trial. It catches six rejected Khudaida trials but also three productive ones. It catches eighteen productive fit trials (thirteen damped), two productive isomerization trials and a productive masking correction. Comparing only the controlling merit row catches two rejected ionic trials and zero Khudaida rejects, but twelve productive fit trials. Comparing each row's predicted change with its own estimate catches zero ionic/Khudaida rejects, but nine productive fit trials and a productive isomerization half step. Increasing the estimate to include the native association floor catches every rejected ionic/Khudaida trial while discarding productive ionic, Khudaida and two-phase action trials.

The protected action half step changes up to 44 ulps and improves the raw maximum from 1.3766765505351941e-14 to 6.439293542825908e-15. Its predicted decrease is 6.546231967915157e-15, versus a 2.378664671080309e-15 assembly estimate, but well below the native association estimate. A productive fit quarter step changes up to seven ulps and has predicted decrease 4.440892098500626e-16 versus assembly estimate 2.043609448574777e-15; raw merit actually falls from 2.664535259100376e-15 to 2.220446049250313e-15. These are direct counterexamples to treating sub-resolution predicted gain as nonprogress.

A predicate restricted to trials after a rejected full step can remove at most three rejected ionic evaluations while retaining the accepted quarter step, leaving 50 evaluations, above the required 49. The historical failed-full norm-band saving also omitted the productive quarter; it is not evidence for a stop restricted to nonproductive trials.

## Work, preservation and line estimates

No rule is proposed. Expected unchanged work: ionic 16/53/17; Khudaida 66/119/74 with all eight cases meeting the requested row criterion; MEA 43/57/44; free T and free P each 7/13/8; masking 6/11/7; unloaded MEA 14/20/15, including its fourth correction; small positive phase 4/9/5.

The three unchanged budgets are observable-relative |ΔT|/T, |ΔP|/P and |Δρ|/ρ ≤ 1e-11; trace-relative |Δn_i|/n_i ≤ 1e-11 for species below 1e-6 of feed; phase-relative |ΔN_i,p|/N_p ≤ 1e-11 for other amounts. The recording has exactly zero changes on all compared states. Those budgets do not loosen row, extensive stationarity, derivative-action, masking or parameter-fit criteria.

New runtime line estimates: Core 0, Equilibrium 0, all other runtime layers 0. Proposed runtime ceiling: 0, because optimization is not adopted. Existing unqualified amendment-2 runtime remains untouched. The presentation-only addition is 39 lines, with a proposed 60-line tooling ceiling; no persistent numerical test is added. Disposable native recording is 105 header lines plus six inserted source lines, retained under `build/issue-272/damped-stop/`.

## Evidence and limits

`trials.md` is the complete per-trial diagnosis table; `trials.csv` owns its numerical values. `coordinates.csv` retains every coordinate; `rows.csv` retains before/after/predicted rows and per-row estimates. `predicates.csv` tests every candidate on every case. `work.csv`, `states.csv`, `checks.csv` and `provenance.json` retain the unchanged values, comparison criteria, source/input hashes and wheel identities. The existing research notebook reads these results through its presentation-only script.

This is numerical verification and local merit prediction on a fixed set of declared states, with no experimental validation or universal impossibility proof. No new readiness or delivered review was run: the assigned endpoint is this design/negative-result handoff. The accountable lead owns issue disposition and removal of the rejected shortcut in a separately authorized change.

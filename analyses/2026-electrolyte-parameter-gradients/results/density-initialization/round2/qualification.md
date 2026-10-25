| Quantity | Count | Maximum relative difference | Limit | Outcome |
|---|---:|---:|---:|---|
| density column | 1 | 0 | 1e-10 | Pass |
| density | 5 | 1.70919409e-15 | 1e-10 | Pass |
| pressure | 5 | 7.44919816e-11 | 1e-10 | Pass |
| slope | 5 | 1.26453027e-14 | 1e-10 | Pass |
| sites | 8 | 6.33893517e-15 | 1e-10 | Pass |
| columns and observables | 28 | 5.18748376e-13 | 1e-10 | Pass |
| physical chart | 39 | 0 | 1e-10 | Pass |
| observables | 204 | 1.07616692e-10 | 1e-10 | Fail |
| physical columns | 2496 | 1.37000329e-12 | 1e-10 | Pass |

Near zero the error is divided by max(abs(reference), 0.01), corresponding to the 1e-12 absolute floor. The native density column and held-trial values/columns include exact optimizer columns. Full-fit physical columns cover the positive-weight observations.

| Complete call | Cold median (s) | Seeded median (s) | Median paired saving | Paired range | Paired median absolute deviation | Pairs | Outcome |
|---|---:|---:|---:|---:|---:|---:|---|
| F2 full | 13.653682 | 8.823412 | 40.190% | 24.507% to 48.054% | 7.864 percentage points | 7 | Pass for timing |
| B4 full | 9.616262 | 9.150722 | 9.759% | -7.602% to 18.183% | 5.633 percentage points | 7 | Pass for timing |

Positive paired saving means a shorter seeded call. F2 requires at least 5%; B4 requires (owner amendment, 2026-10-06, replacing at least 5% for B4: no slower beyond noise) a nonnegative median or a negative median within the paired median absolute deviation. Each pair returns matching dimensions, solver outcomes and trial-failure counts.

| Key cold solve | Cold iterations / evaluations / Hessians | Seeded iterations / evaluations / Hessians | Maximum relative state difference |
|---|---:|---:|---:|
| case-6/figure-2-point-7/cold | 16 / 53 / 17 | 16 / 53 / 17 | 0.0 |
| case-7/flash-5918.5286-Pa/cold | 43 / 57 / 44 | 43 / 57 / 44 | 0.0 |

The manufactured-recovery variance factor is 3.47653699266e-27, versus 3.47493504088e-27; its absolute difference is 1.60195178e-30, against 6.66666667e-21. The floor is twice the final-cost criterion divided by the three residual degrees of freedom.

The 5 independent pressure stencils pass 1e-6 max(1, abs(reference)); the largest error is 0.000702397446 of its allowance. Density, sites, slope and EOS pressure are retained for the ethane closure, the 45 K natural case, first ionic-reference rung, pure solvent and finite activity closure.

| Invalidation or fallback | Result |
|---|---|
| Missing branch witness after a parameter change | Cold selected state preserved |
| Liquid density offered for vapor, on the wrong side of the witness | Cold vapor selected |
| Irregular retained closure | Offer rejected; cold result preserved |
| Density outside the anchor window | Recorded left-domain cause, offered density and convergence operands; cold recovery matches |
| Empty cold domain | Anchor and final cold refusal retain the same cause |
| Same liquid request, new upper branch | Cold upper density selected |
| Failed evaluator trial, then prior successful-point replay | Predictions, chart and columns preserved exactly |
| Changed-composition anchor and unresolved fold | Existing P2 numerical checks pass |

The failed F2 osmotic coefficient at the manufactured optimum is 0.983825866820915, versus cold 0.983825866926791, for KCl at 0.002 mol/kg, 298.15 K and 100000 Pa. The absolute difference is 1.05876086e-10; the relative difference is 1.07616692e-10. Every physical-chart entry is bit-identical. Replaying with and without hydration holds gives bit-identical cold predictions, so re-solving holds does not explain this difference. The scalar pressure-closed-property cause remains unresolved; neither a new method nor a changed tolerance has been selected.

All other retained equivalence rows pass. Both full B4 starts pass their same-point cold replay. The fit timing evidence does not close the failed observable comparison. These calculations provide numerical verification and measured runtime, with no physical validation, global-optimum claim, parameter adoption, delivered review or manuscript promotion.

**Superseded on 2026-10-06 by round 3.** The owner's equivalence standard (CODING_STANDARDS.md, 2026-10-06) is 1e-10 relative with absolute floors for near-zero quantities. The practical osmotic coefficient is phi = -ln a_w / s with s = M_w sum(nu_i m_i); at low molality ln a_w is a near-zero difference of order-one terms, so the 1e-12 floor belongs to ln a_w. Carried to phi, the per-row allowance is max(1e-10 abs(phi), 1e-12/s), which grows as 1/m. For KCl at 0.002 mol/kg, s = 0.0180153 x 2 x 0.002 = 7.206e-5 and the allowance is 1.388e-8; the observed 1.0588e-10 is 0.0076 of it, equivalent to a ln a_w difference of 7.6e-15. The row passes, and the "Fail" outcome and "unresolved" cause above are superseded. The delivered review of 5fec07f8..ab797752 found no wrong-branch or stale-seed path and attributed the residual differences to liquid stiffness at the same converged state. The B4 timing gate is a nonnegative median paired saving, or a negative median within the paired median absolute deviation; B4 passed it at 9.759%. The round-3 record holds the regenerated comparisons.

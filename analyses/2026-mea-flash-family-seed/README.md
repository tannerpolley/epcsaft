# Reacting MEA flash family before and after PR #370 (#350 slice 2, #356)

**Question.** Does the true-species split start with the amount step rule on free phase amounts (PR #370)
recover points of the reacting MEA flash family that the K-factor start lost, without losing any, changing their
states, or slowing the key flash?

**Declaration.** The case-7 reacting flash of `analyses/2026-tpd-phase-discovery/scripts/baseline.py` (MEA, water,
CO₂ and their true species; one free liquid, one free vapor; fixed T and P; row tolerance 1e-10) with
3 temperatures (313.15, 318.15, 333.15 K) × 8 loadings (1e-5 to 0.5 mol CO₂/mol MEA) × 3 pressures, plus the key
flash (318.15 K, 5918.5286 Pa, loading 0.25558): 73 points (`scripts/family-points.txt`). Feed 9.7292821 mol MEA,
76.9751925 mol water, CO₂ = loading × MEA.

**Result (numerical verification only).**

| | before (main 7fa511b6) | after (main 6bbaaf1e, rerun 2026-10-09) |
|---|---|---|
| converged points | 58 of 73 | 72 of 73 |
| lost / new | — | 0 / 14 |
| key flash work (iterations/evaluations/Hessian calls) | 43/57/44 | 19/42/20 |
| work summed over the 58 common points | 2778/4388/2835 | 1866/3479/1922 |

The 58 common states agree within 1.5e-11 relative (worst: loading 1e-5, 318.15 K, a mole fraction). The rerun on
current main is identical to the retained PR #370 branch run at all 73 points (success, work, amounts, mole
fractions, densities). Current main cannot select the K-factor start, so "before" is the retained run from #356.
The notebook computes every number from `results/comparison.csv`.

**Limits.**
- Numerical verification on one declaration; no physical validation and no global-stability claim.
- Row 66 (loading 0.5, 313.15 K, 5061.6 Pa) still fails: "phase L vanished" before, "phase V vanished" after (#358, handed to Amine-Thermodynamics #172 as a model decision).
- 12 of the 58 common points need more of some work component; worst 49/69/50 → 145/446/144 (loading 0.5,
  318.15 K). The amount step rule has no upper bound on added iterations.
- Wall time is not analyzed here; PR #370's timing gives relative changes only, on one shared host.

**Files.**
- `scripts/family.py`: runs the 73 points through public `solve_equilibrium`, one native solve per process.
- `scripts/compare.py`: joins the three runs into `results/comparison.csv` with the frozen #350/#356 state metric.
- `scripts/render.py`: `figures/family.png` from the CSV.
- `results/family-before-7fa511b6.jsonl`: retained main-before run (#356 build, `family-K-7fa511b6.jsonl`).
- `results/family-pr370-branch.jsonl`: retained fix-round run of the #370 branch (#356 fix, `family-fix.jsonl`;
  bitwise identical to a60d9d9f, `family-identity.txt`).
- `results/family-main-6bbaaf1e.jsonl`: rerun on current main.

Raw evidence: `build/archive/issue-356-slice2.tar.zst` in the main checkout. Issues #350, #356, #358; PR #370.

## Named-release rerun (#680)

The [research notebook](notebook.qmd) reads the new results under
[results/main-8298c602/](results/main-8298c602/). The earlier retained files remain
under their original identities. The new provenance records the downloaded
wheel, its verified hash, the actual installed Engine identity and the declared
producer. Unavailable portions retain their stopping states and causes.
Independent delivered review and owner promotion remain open; the manuscript
is unchanged.

## Public release v1.0.0 rerun (#845)

The [notebook](notebook.qmd) reports the current-method family on the public
release and compares it with the earlier recovery and the unpromoted
`main-8298c602` run. [Retained results](results/v1.0.0/summary.json) own the counts,
state comparisons and numerical checks; [provenance](results/v1.0.0/provenance.json)
and [source identity](results/v1.0.0/source-identity.json) identify the detached
release image, DOI and public tree. No Engine or criterion changed.
Independent delivered review and owner promotion remain open. The manuscript
is unchanged.

# Conserved-inventory bounds: numerical qualification

The native amount bound and aggregate fraction bound share the HiGHS program
used by optimizer boxes. They impose conserved rows, nonnegative supported
amounts and the requested fraction floor, not equilibrium, energy or additional
charge constraints. No derivative of the nonsmooth maximum is supplied.

## Bound comparisons

`comparison.csv` retains the historical SciPy maxima/logarithms alongside native
values, errors and residuals. The two copied Born requests are unchanged; input
hashes, wheel SHA-256 and both types of witnesses are in `witnesses.json`.
The floor is $10^{-20}$. These aggregate fractions equal the liquid fractions
only for the original sole finite liquid with zero-amount vapor.

| Case | Bicarbonate maximum | Hydronium maximum |
|---|---:|---:|
| 1 | 0.01716948693975162 | 0.45247530902624283 |
| 2 | 0.07204451382562489 | 0.47991282246917943 |

Maximum absolute fraction error is $1.11\times10^{-16}$ (criterion $10^{-10}$),
logarithm error $6.22\times10^{-15}$ (criterion $10^{-8}$), and normalized
feasibility residual $3.33\times10^{-16}$ (criterion $10^{-10}$).
`global_scan.inventory_limits` now reads these native bounds; its `linprog`
implementation and three corresponding allowlist locations are removed.
Main's #234 retirement removes the diagnostic chart code while preserving its
historical results. The standalone `global_scan.py --case 1/2` now declares the
copied problem and reads only native inventory bounds.

`analytic.csv` retains the manufactured $n=(2-z,1-z,z)$ checks: amount maxima
$(2,1,1)$ mol, fraction maxima $(2/3,1/3,1/2)$, and floored maxima
$(0.6633333333333333,0.32666666666666666,0.485)$ at $f=0.01$. Scaling the feed
by seven multiplies amount maxima by seven and leaves fraction maxima unchanged.
All these value errors are zero in the retained binary64 calculation; maximum
normalized residual is $1.11\times10^{-16}$, against $10^{-10}$ criteria.

Signed inventory $n_0-n_1=0$ has unbounded amounts but fraction maxima $1/2$.
For $n_0-n_1=1$, $f=1/2$ gives an empty finite set even though the transformed
point $(y_0,y_1,\tau)=(1/2,1/2,0)$ exists. At $f=0$, the finite witness
$(1,0,3)$ establishes a nonempty set and the species-1 supremum is $1/2$.
`analytic-diagnostics.json` retains these witnesses and LP causes. A returned
zero-$\tau$ maximizer does not itself exclude other finite maximizers, as the
zero-inventory case demonstrates. Native specification checks separately cover
infeasible support and structural zero.

Delivered review found that HiGHS dropped the transformed coefficient for
$n_0-n_1=10^{-10}$ or $10^{-12}$ mol: an empty half-floored inventory was returned
as feasible. `small-balance-before.csv` retains that reproduction. The shared
owner now solves with the exactly equivalent scaled variable $q=s\tau$, then
validates the reconstructed witness against the original equations. Both small
nonzero cases now return an empty finite set, zero residual and no maximum;
the zero-inventory case remains feasible with maximum $1/2$. This corrects a
coefficient-resolution defect, not a change to the inventory constraints or
tolerance. The reviewer also found unnecessary objectives for off-support
species; the ordinary box now requests only the compiled support union.

The existing #172 trace-inventory checks pass at $10^{-9}$ and $10^{-12}$ mol.
Recovered $x_{s2}$ are $7.499999992499985\times10^{-10}$ and
$7.500000035072494\times10^{-13}$: relative errors about $2.1\times10^{-15}$
and $4.68\times10^{-9}$, below $10^{-6}$. Their temperature-action relative
errors are below $4.68\times10^{-9}$, against $5\times10^{-6}$.

## Ordinary-solve cost

The explicit fraction-feasibility program never runs in ordinary amount boxes.
The shared amount call replaces, rather than supplements, main's species LPs.
The frozen #223 declarations are reused without altered starts or acceptance:
Khudaida Figure 2 point 7 ionic LLE, and MEA flash at 318.15 K, 5918.5286 Pa.
`cost-main-source.txt` identifies main `40e49c5e`; its runtime matches `edf7672e` and is unchanged
in the later merged main `f593e803`. Measurements used one installed non-editable
wheel at a time, `nice -n 10`, OMP and OpenBLAS threads both one. Every cost
sample has the same Ipopt iteration/evaluation/Hessian counts and identical
phase amounts, compositions and residuals between main and branch.

| Frozen declaration | Main best of 3 (s) | Branch best of 3 (s) | Ipopt iterations, main/branch | Change |
|---|---:|---:|---:|---:|
| Ionic LLE | 0.116668683 | 0.093580790 | 16 / 16 | -19.79% |
| Reacting MEA flash | 0.167978461 | 0.169816676 | 43 / 43 | +1.09% (flagged) |

`cost-{main,branch}.csv/json` retain the final interleaved measurements, literal
inputs, reactions, packet and wheel hashes, and numerical outcomes. An initial
noninterleaved run showed MEA +6.65% and LLE -14.34%; its samples are retained as
`cost-initial-*`. The first interleaved run, retained as `cost-pre-review-*`,
had MEA -1.39%; after the required coefficient/support corrections the final
interleaved best-of-three has MEA +1.09%, explicitly flagged above. Main MEA
samples span 0.168–0.327 s and branch samples 0.170–0.218 s. Such variation
prevents attributing that 0.00184 s difference to this change or establishing
a speedup or precise sub-percent cost bound. Ipopt work is identical. The
review-required support correction removes the extra off-support objectives;
no additional LP remains in the ordinary path. No reproducible change-caused
slowdown is demonstrated, but the observed final MEA increase is not hidden.

These checks establish numerical verification of inventory bounds and preserved
ordinary solves at the stated declarations. They establish neither physical
validation nor thermodynamic equilibrium/stability from an inventory bound.
Notebook/prose consolidation remains with #234 batch 15.

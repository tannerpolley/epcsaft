# Declared phase equilibrium and explicit tangent-plane qualification

The corrected wheel passes the supplied neutral, reacting, permittivity-form,
refusal and 100 mol extensive-margin checks and reproduces all **4,374 frozen
numerical rows exactly**. Its SHA-256 is
`8470b49a9283528b83dbbdab46dc62dee0fc70e5e2149f8f218b823dd85812cb`.
The [notebook](notebook.qmd) presents the retained numerical evidence. The initial delivered review returned
**Not supported** for an intensive/extensive addition-margin defect; the native
amount-weighted correction and exact 100 mol falsifier now pass.
**The same reviewer returned Supported, with no required defects, at
`0f9bbe59b13bedc1b0884908ce7ba55be5194224`.** This verdict covers the qualified
cases and limits below. The notebook remains unpromoted for manuscript use;
manuscript work was not requested and is not a tool-qualification gate.
The [review record](../../docs/scientific/reviews/issue-223-delivered-review.md)
states the final verdict and separate declared-conditioning scope. The 100 mol
public-reuse failure comes from the lead and retained ordinal-5 diagnostic;
the reviewer did not independently inspect that reuse evidence. These results
verify calculations and compare implementations; they establish neither
physical validation nor global stability. No admitted charged trial has been
qualified as resolved descent followed by ionic phase addition; the manufactured
reacting fixture stays one phase.

The [owner-adopted permittivity admission](https://github.com/tannerpolley/ePC-SAFT/issues/223#issuecomment-5947638210)
and [Ready review](../../docs/scientific/reviews/issue-223-permittivity-admission-readiness.md)
allow charged trials only for declared density–temperature or dipolar-perturbation
permittivity, whose defining equations approach vacuum permittivity as density
tends to zero. This is a necessary search-domain condition, not proof of physical
applicability over the entire domain. Bülow–Ascani–Held Part I §2.1 and the
[Ascani–Held 2021](https://d-nb.info/1232858226/34) framework (pp. 3–4, Eqs. 10–13)
supply no density or branch exclusion inequality for charged trials. Density-independent choices raise
`ParameterError: permittivity model has no low-density limit; declare phases`
before numerical trial or discovery-driver work. Neutral-only caller families
remain available on ionic models. Declared calculations retain their physics.

## Current numerical evidence

| Calculation | Retained rows | Timing samples | Named indeterminate trials |
|---|---:|---:|---:|
| [Neutral cases 1–4](results/qualification-8470b49a/neutral.json) | 31,677 | 57 | 26 |
| [Engine comparison with refined FeOs](results/qualification-8470b49a/feos.json) | 3,359 | 9 | 0 |
| [Manufactured reacting reference](results/qualification-8470b49a/reacting.json) | 2,726 | 21 | 4 |
| [Permittivity forms and limits](results/qualification-8470b49a/permittivity.json) | 1,650 | 39 | 3 |
| [Ionic declarations, refusals and neutral-only searches](results/qualification-8470b49a/ionic.json) | 5,434 | 48 | 0 |
| [100 mol default discovery and margin](results/qualification-8470b49a/gibbs-scaling.json) | 2,231 | 3 | 1 |

Each numerical CSV preserves round-trip decimal values and exact binary64 bits,
including signed zero. Compact JSON retains requests, diagnostics, provenance
and named limitations. Current refusal metadata preserves both physical-source
and request annotations. Byte-identical current reacting/permittivity
rows link the two 0b297e6f CSVs kept from that run; cost references were rebased with all
observed operands checked unchanged. Original JSON/cost hashes identify the
offline retention changes. Executed-script hashes remain those recorded at run time.

The 177 timing samples cover 59 operations. [Median-three times and work](results/presentation-8470b49a/cost-summary.md)
are derived from the six `*-costs.csv` tables by `scripts/present.py`, which
checks identical work counts across all three samples. Sample-zero native work
links the same report operands. Declared density roots are existing
`solve_density` spans; density evaluations are `evaluate` spans with that
ancestor. Explicit native work includes initialization, references, branch
admission, declarations and trials. Trace counting adds no EOS calls. Host load
contributes to elapsed times; these are observations, not an isolated speed
benchmark or a universal wall-time threshold.

The refusal rows count observed driver callbacks, not all EOS evaluations.
Discovery owns a private trace and compiles before admission, including a
generic-point Jacobian and EOS work that the outer capture cannot see. Source
inspection shows admission before the first declaration driver, density root
and trial; zero callback columns do not mean zero compile/EOS work. Refusal
wall times include that compilation.

On 8470b49a, case-6 neutral-only discovery returned G/RT=-8.662995081104183 mol,
above the declared ionic split's -8.663320942099162 mol. The restricted family
returns an ion-free second liquid; its finite no-witness outcome is limited to
those caller supports and does not reproduce the ionic split's lower energy.

## Frozen declared calculations

The original wheel is from `443d828908f4ad5127ab14384f97635dca30bd30`, SHA-256
`298a360e532f6f4c5e7ef70c5be45e8b9852310d1c388d6225c504714c3dc410`.
[Original identities](results/baseline-identity.json) and [final replay metadata](results/after-8470b49a/)
identify wheel, dependencies, physical inputs, declarations, script hashes and
source modification times. No equality tolerance is applied. Identical replay
CSV rows link the originals rather than being retained twice.

- [Declared rows](results/before/declared.csv): 1,841 rows for cases 1, 3, 6 and 7,
  including the reactive MEA bubble, flash and liquid states; all available
  first/second temperature actions, residuals, extensive objectives and counts.
- [Continuation rows](results/before/continuation.csv): 838 rows for the accepted
  methane/ethane quarter-step bubble path from 180 to 181 K and the half-step
  branch-loss boundary, including the last attempted point.
- [Fitting rows](results/before/fitting.csv): 217 rows for #219 held-NaCl observations
  and fit at 298.15 K, 100000 Pa and molalities 0.1, 1 and 4.5 mol/kg; predictions,
  exact columns, residuals, weighted least-squares cost, fitted coordinates,
  covariance and work counts.
- [ASL rows](results/before/asl.csv): 1,478 rows for the live Khudaida anchor and
  `start-0.9`/`start-1.1` feed points, ten phase-component amounts and first/second
  four-feed derivatives through public actions and ASL.

All four replay captures contain **zero tangent-plane trials**. Values, actions,
objectives, residuals, text and exposed iteration/evaluation/Hessian counts are
exact. Amounts use the declared mol feed basis, compositions mol/mol, densities
mol/m³, T in K, P in Pa and extensive G/RT in mol. ASL feeds and amounts use mol/s;
first derivatives are dimensionless and second derivatives have units s/mol.
The fitting cost is the native weighted least-squares objective.

The parent-owned [final ASL hidden-trace comparison](results/asl-native-after-8470b49a/summary.json)
finds all 396 work rows and all 1,478 numerical rows identical to the
[original native capture](results/asl-native-before/summary.json): 30 ASL calls,
9 driver finalizations, 66 iterations, 130 evaluations, 75 Hessian calls and
zero trials. The current library SHA-256 is
`12f4a434a00d824da6065dba579c0b07b88f300a646067167b92edb95fda2ba1`.
The static ASL archive has a separate hidden trace pointer. ROOT inspected the
current Linux/glibc ELF pair at `0x6ebe50`, offset `0x28`, and Scope instructions
at `0xa4500`, using unchanged Trace/Convergence headers and helper. The shim
restores the pointer after every call. This diagnostic is specific to the
identified library and C++ layout; its offset is not portable to another build.

All four full public traces are retained once as `.csv.gz` under `results/after-0b297e6f/`. Their complete decompressed bytes match the 8470b49a captures exactly.
[Trace identities](results/after-8470b49a/trace-retention.json) retain compressed
and complete decompressed CSV hashes: 190,868,907 CSV bytes in 32,415,807 compressed
bytes. No span, event or operand was filtered. No original-wheel public traces
exist, so their `trace_comparison` entries remain unavailable. The separate ASL
work comparison is available and exact.

## Reacting reference and failure diagnosis

At 298.15 K and 100000 Pa, the existing density–temperature selection supplies
water/Na/Cl in that order. A single liquid state at x*=(0.98,0.01,0.01) defines
frozen g_i=-RT[ln(rho*x_i*)+mu_i,res/RT] and K=exp(-nu·g/RT). The
[original frozen operands](results/diagnostics/reacting-reaction-mapping/reacting.csv)
remain unchanged across all later attempts and wheels. For feed (0,1,1) mol and
nu=(2,-1,-1), the final declared amount error is 4.884981308350689e-15 mol against
1e-8 mol. Discovery finds one phase, classified residual 1.0441577530190186e-14
against 1e-10, zero moiety errors and delta G/RT=-1.2434497875801753e-14 mol.
The returned declaration recompile/test/re-solve has scaled Gibbs difference
1.7763568394002505e-15 against 1e-9. The pinned/free-pressure variant uses the same
g, nu=(1,-1,-1), amount fraction 50/99 of the 2 mol feed, expected amounts
(98/99,1/99,1/99) mol and amount error 3.469446951953614e-18 mol. No free-pressure
discovery was requested. Independent affinity is checked only where every
reaction participant has positive fraction and the phase has positive amount;
structural-zero phases have a named unavailable affinity, while shared rows,
charge and moieties remain checked. The largest measured affinity/RT magnitude
is 1.2073675392798577e-13, below the unchanged 1e-10 criterion.

The [6e4b4e29 initial discovery failure](results/diagnostics/reacting-initialization-6e4b4e29/reacting.json)
and [diagnosis operands](results/diagnostics/reacting-initialization-6e4b4e29/diagnosis-operands.csv)
remain retained. The ideal-speciation seed had literal chloride zero and log
minus infinity. Supplying the existing positive charge-neutral base before the
same ideal-speciation step recovered the defined reference. ROOT changed only
the explicit reacting discovery initializer, indexed by the admitted support;
ordinary declared initialization, default trial starts, budgets and tolerances
were preserved. The [diagnosis provenance](results/diagnostics/reacting-initialization-6e4b4e29/diagnosis-provenance.json)
pins source, mapped inputs, frozen g/K and the bounded probe. Duplicate prefix
rows in the two failed analysis attempts are linked through `row-retention.json`;
reconstruction recovers their original CSV bytes exactly. Distinct attempted
elapsed times and original JSON provenance remain unchanged.

## Earlier evidence and current scope

The [historical MEA descent diagnosis](results/diagnostics/phase-domain/diagnosis.md)
and [operands](results/diagnostics/phase-domain/operands.csv) retain the selected
solvent-only mathematical model's charged descent at density
2.3214665376575057 mol/m³, epsilon=71.48980191703282 and Born molar Helmholtz/RT
=-322.0579048293252. Direct delta G/RT=-1.9951025594702504e-5 mol, with zero errors
in all nine component inventories; donor-pressure reclosure preserves descent.
This establishes neither a physical ionic phase nor three-phase equilibrium.
Fresh charged requests on that form now refuse under the adopted rule; the
historical decrease is preserved. Earlier case-6 charged-search passes on
384d608/76b13fd3 apply to those wheels and the former admission domain only.

The original failed bubble initializer and the exact composition-derived fixture
correction are retained in [density-bubble operands](results/diagnostics/density-bubble/operands.csv)
and [provenance](results/diagnostics/density-bubble/provenance.json). The accepted
fixture bubble has P=3094.0685042628947 Pa and liquid density
55795.196278370924 mol/m³; this source probe does not alter frozen baseline starts
or establish general cold-start behavior. Earlier association-request,
balance-basis, returned-declaration and FeOs coordinate failures remain in
`results/diagnostics/` with their genuine logs.

The [final native provenance](results/native-8470b49a/provenance.json) links the
92 comparisons/1,306 rows owned once by the 6e4b4e29 CSV. All 90 asserted numerical
inequalities pass their original tolerances; two entries describe optimizer
topology. Current charged-first stationary-dual defect is
6.414777798225091e-15 against 1e-10, not the historical 384d608 value.
[Verification identities](results/verification-8470b49a.json) associate the final
wheel with the lead's required numerical and source checks.

The [refusal counting diagnostic](results/diagnostics/refusal-scope-count-0b297e6f/scope-count.json)
retains an analysis-only stopping point. Native `solve_equilibrium` starts its
trace Scope at line 1042 before discovery dispatch at line 1043. The refusal
therefore retains one enclosing span, but zero finalized drivers, declaration
scopes, roots, trials, iterations, evaluations and Hessian calls. The checker
now distinguishes the enclosing span from a driver invocation; the exact typed
reason and zero-driver/trial requirements are unchanged. The private compile/EOS
work described above remains outside this outer capture.

## Extensive margin and separate declared-conditioning limitation

The native addition comparison now uses the extensive estimate
sum_p N_p epsilon_p, in mol for G/RT, with the same per-mole phase estimates.
[The 100 mol falsifier](results/qualification-8470b49a/gibbs-scaling.json)
uses case-1 fractions (0.6,0.4), T=293.15 K, P=100000 Pa, original parameters,
row tolerance and default native options. It makes only three default Discover
calls and reads the existing 1 mol neutral reference. Maximum composition
error is 1.6653345369377348e-15 against 1e-9, phase-amount scaling error
2.4158453015843406e-13 mol and component-inventory error
7.105427357601002e-15 mol against 1e-12 times the 100 mol feed. Scaled Gibbs
error is 6.698679190032706e-15 against 1e-10; relative density difference
2.1094237467877974e-15 is reported without another tolerance.

The retained first candidate amount is exactly 100 mol. Its positive per-mole
estimate is 5.312046591948305e-11; native and independently assembled extensive
estimates both equal 5.312046591948305e-9 mol. Relative difference is zero
against 1e-13. The old unweighted counterfactual differs by 0.99 and fails the
same criterion. Final-state scaling alone would not falsify the missing amount
factor. Median elapsed time is 0.15303935300471494 s, with 17 driver invocations,
14 trial invocations, 162 iterations, 302 evaluations and 179 Hessian calls;
25 density roots/1152 root evaluations and 1592 energy evaluations are counted.
The final near-pure-1 trial remains indeterminate, with slope
-4.387731590963851e-15 against estimated error 8.331709354552604e-11.
It is the 34th named indeterminate trial across current qualification.

The [original cold attempt](results/diagnostics/100-mol-cold/gibbs-scaling.json),
and three failed cost observations are preserved unchanged; the executed
diagnosis sources under `scripts/diagnostics/` remain at record commit `21d7a112`. The cold declaration reached a coincident stationary point, with
classified row residual 1.2816303435237685e-15 and G/RT=-54.058413564546925 mol,
then was rejected and its public state cleared. Each attempt used
55 iterations, 98 evaluations and 56 Hessian calls. The separately authorized
[ordinal-5 probe](results/diagnostics/gibbs-scaling-ordinal5/provenance.json)
also found that public reuse of the discovered declaration and supplied state
coalesced after 11 iterations, 18 evaluations and 12 Hessian calls. Neither is
labelled bad input, physical nonexistence, successful reuse or a fixed path.
Zero reduced curvature in those rejected records is an uncomputed default,
not a stability calculation.

[Owner issue #243](https://github.com/tannerpolley/ePC-SAFT/issues/243) resolved
that declared-conditioning question: Ipopt now measures each free
phase amount in units of its start amount (ALG-EQUILIBRIUM), and the 100 mol cold
declaration and the returned-declaration reuse keep the split
(`engine/tests/equilibrium/test_feed_scale.py`); the records above are kept as
the pre-#243 evidence. The original #223 reuse
cases and all frozen no-overhead rows pass; those checks do not promise universal
reuse or cold-start convergence at every feed scale. This separate limitation
does not block the accepted #223 extensive-margin falsifier.

[Ordinal-6 full operands](results/diagnostics/gibbs-scaling-ordinal6/operands.csv)
and [provenance](results/diagnostics/gibbs-scaling-ordinal6/provenance.json)
retain 1,127 numerical operands from one independent default 100 mol discovery,
its exact source and log, and original JSON/selected-CSV hashes. The final
three-sample qualification is separate from that diagnostic. Ordinal 5's
recorder stopped at cleared arrays; no lost margin values were invented.
#341 removed the superseded 384d608, 411c8, 6e4b4e29, 76b13fd3 and 0b297e6f
snapshots, their presentations and the HTML renders; they remain at record
commit `21d7a112`.

## Reproduction and presentation

Read the full current #223 design and source identities first; check branch,
dirty state and wheel/source modification times before every heavy call. Use one
analysis process, the existing installed wheel and a new output directory.
`qualify.py --stage neutral` precedes the FeOs Engine comparison, then `--stage
reacting --reacting-reference .../reacting-reaction-mapping/reacting.csv`,
`--stage permittivity` and `--stage ionic`. FeOs comparison uses `--reference
results/qualification-76b13fd3/feos.json`: both independently refined references
and achieved residuals stay owned once; no FeOs import or rerun is required.
After the four frozen baseline sections, `--stage gibbs-scaling --neutral-reference
.../qualification-8470b49a/neutral.csv` performs only default 100 mol discovery.

```sh
timeout 300 nice -n 10 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python analyses/2026-tpd-phase-discovery/scripts/qualify.py --stage neutral --output analyses/2026-tpd-phase-discovery/results/new-run
timeout 180 nice -n 10 env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python analyses/2026-tpd-phase-discovery/scripts/baseline.py --section declared --output analyses/2026-tpd-phase-discovery/results/new-replay
```

The baseline sections are `declared`, `continuation`, `fitting`, `asl`.
`scripts/replay_compare.py BEFORE AFTER OUT.json` compares two replays of all four
sections under #243's quantity tolerances; `scripts/physical_comparison.py` runs
the #243 42-row feed-scale matrix and `scripts/neutral.py` the 70-case neutral
methane/ethane sweep. `tools/check-equilibrium-corpora` runs all three on two wheels.
`scripts/present.py --qualification results/qualification-8470b49a --output
results/presentation-8470b49a` formats retained values only. Completed duplicate
or empty logs link their retained originals in
[log links](results/completed-log-links-8470b49a.json); genuine failed logs remain.
Render the notebook without execution from its retained tables into `build/`.

## Named-release rerun (#680)

The [research notebook](notebook.qmd) reads the new results under
[results/main-8298c602/](results/main-8298c602/). The earlier retained files remain
under their original identities. The new provenance records the downloaded
wheel, its verified hash, the actual installed Engine identity and the declared
producer. Unavailable portions retain their stopping states and causes.
Independent delivered review and owner promotion remain open; the manuscript
is unchanged.

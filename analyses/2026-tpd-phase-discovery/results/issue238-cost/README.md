# Declared equilibrium starts and default-solve cost

The installed wheel SHA-256 is `4e00f1536cddea911baa5cf157892a8fc6def636aceeb6aa9b5d74ff1fa1540d`, built from runtime source at `daa98f1f`, integrating main `875526d0` density diagnostics. Manufactured and default-solve evidence below was regenerated on this wheel. Born and Uyan diagnosis remains scoped to its recorded `624881f1` wheel; it is not a new-wheel application acceptance claim. The final full wheel gate passed 114 Python tests, six native targets, ruff, mypy, analysis layout and architecture; the affected check reported no changed paths. [Independent review](../../../../docs/scientific/reviews/issue-238-delivered-review.md) is Supported with no required corrections.

The manufactured exchange case at 200 K and 5 MPa selects the second caller start after rejecting negative packing; the third is unattempted. Both phase amounts are 1.3 mol within 4.44e-16 mol, composition errors are at most 1.11e-16, and exchange ratios are 1.9999999999999996 and 0.49999999999999883. Original rows reach 5.773159728050814e-15 against 1e-10. A perturbed start with one aggregate optimizer iteration retains residual 1.224413540303225 and no selected state. See [manufactured.csv](manufactured.csv) and [work diagnostics](manufactured-work.json). Coordinates, previous results, partial packing, pin conflicts and free-pressure previous-result reuse are checked numerically by `engine/tests/equilibrium/test_declared_starts.py`.

## Default ordinary solves, not multistart calls

The unchanged harness declares two frozen inputs and runs each three times, single-threaded at nice 10. All 1,176 retained state rows, including binary64 values, and all objective/residual/work counts agree exactly. Evaluations include Ipopt objective callbacks and exact polish evaluations; they are not isolated objective-callback counts.

| Case | Main iterations/evaluations/Hessians | Branch | Main best (s) | Branch best (s) | Change |
|---|---|---|---:|---:|---:|
| Ionic LLE, figure 2 point 7 | 16/49/17 | 16/49/17 | 0.104208528 | 0.065512062 | -37.13% |
| Reacting MEA flash, 5918.5286 Pa | 43/57/44 | 43/57/44 | 0.159753735 | 0.112479235 | -29.59% |

Both wall differences exceed the approximately 10% reporting threshold. They are flagged as shared-host load variation, **not an algorithmic speedup**. No changed-work or induced-slowdown evidence occurs in these two cases. This does not establish timing for other fluids or configurations. Every run is in [main-cost.csv](main-cost.csv) and [branch-cost.csv](branch-cost.csv); input/wheel provenance is alongside them. Byte-identical main/branch states are retained once in [main-states.csv](main-states.csv); the comparison records the common SHA-256. [comparison.json](comparison.json) records the exact agreement and timing flags. The main wheel SHA-256 is `da56e55eab2dac8b75c2497da07e4895d838d82d2d78b305e4a30b81bfef3950`.

The elapsed deadline controls state availability, separately from numerical convergence. No post-deadline state is published; measured rows retain the common nonfinite/converged/floor/budget classification, native optimizer termination, controller cause and visible overrun. The review request to force a budget numerical stop was withdrawn because it contradicted the accepted #238 design and ALG-TERMINATION.

Uyan's attempted reference conversion did not preserve its historical values and is excluded from this delivery. That starts-only handoff was superseded by the source-ledger native conversion; the allowlist is empty. Its [current notebook](../../../reproduction/2015-uyan/notebook.qmd) retains changed predictions without attributing them solely to the convention. See the [Uyan diagnosis](../../../reproduction/2015-uyan/README.md). Born's twelve declared attempts are retained without claiming nonexistence; see the [Born analysis](../../../2026-born-unavailable-states/README.md).

Runtime growth against active main `875526d0` is 374 added, 50 removed, +324 net lines; integrated total 22,013 (Core 917, Contributions 1547, Properties 4763, Equilibrium 7170, regression 2398, Surface 532, bindings 1945, Python 2741). The documented 275-line task target was exceeded to complete partial/pinned admission, between-kernel deadlines and measured failure diagnostics. No architecture ceiling was amended.

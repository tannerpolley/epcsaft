# Exact residual-polish termination

**On merged main the amendment-2 shortcut passes every affected check.**
The [physical comparisons](results/merged-main/physical-state-comparison.csv),
[work](results/merged-main/work.csv), [timing](results/merged-main/timing/summary.json),
[criteria](results/merged-main/criteria.json) and
[qualification provenance](results/merged-main/qualification.json) retain the
requalification under `CODING_STANDARDS.md`, PR #311. The
[four checks](results/equivalence/four-checks.csv) record the earlier base.
The two native checks that failed at the rounding floor now use derived criteria:
the action forward bound includes the computed Jacobian's backward error, and the
wide-ionic oracle checks a volume accuracy target of 2^-100.

The earlier candidate failed the former criteria, and the separate damped-trial
predicates did not meet the work and productive-step constraints.
`results/resumed/` and `results/damped-stop/` retain those historical findings.
Their original conclusions apply to their original criteria and inspected snapshots.

Can a fresh, sufficiently small Newton correction avoid damped trial work while
preserving the physical states of unchanged polishing? Issue
[#272](https://github.com/tannerpolley/ePC-SAFT/issues/272) owns the accepted
amendment-2 design and qualification budgets. [ALG-EQUILIBRIUM](../../engine/docs/science/algorithms.md#alg-equilibrium)
owns scaling, rank, forward-error admission and best-point selection.

`results/equivalence/` retains the fresh installed candidate and main comparisons.
`results/resumed/` retains the earlier comparisons:
returned state rows, physical differences in the three approved bases, work
counts, detailed polish traces, declarations, wheel identities and paired timing.
`notebook.qmd` reads retained values without calling the Engine. Before rendering it,
run `python3 analyses/2026-equilibrium-polish-masking/scripts/present_polish.py`, which
formats two of its included tables into `build/notebook-data/`.
The numerical-policy values are read from the installed candidate wheel.

The masking declaration comes from
`engine/tests/equilibrium/test_trace_reliability.py`. The free-temperature and
free-pressure declarations reuse `test_cold_vapor.py` and `test_continuation.py`.
The small positive vapor amount lies above the phase-disappearance cutoff;
its trace distinguishes scaling by its own amount from scaling by the feed.
The native diagonal-system falsifier has an analytically known forward correction
that exceeds the local threshold despite a tiny backward residual.

The earlier norm-wide residual-band stop is rejected. `results/rows.csv`,
`operands.csv`, `steps.csv`, `trials.csv`, `work.csv` and
`physical-state-comparison.csv` retain its trace-inventory counterexample.
`results/correction-size/` retains the disposable fresh-versus-previous-direction
comparison. `results/amendment-2/` preserves the paused, unqualified build.
These historical records are evidence for the design, not current executors.

After installing each corresponding non-editable wheel, retain states with:

```bash
.venv/bin/python analyses/2026-equilibrium-polish-masking/scripts/qualify_polish.py baseline build/issue-272/resumed/main-wheel/epcsaft-0.2.0.dev0-cp313-cp313-linux_x86_64.whl
.venv/bin/python analyses/2026-equilibrium-polish-masking/scripts/qualify_polish.py candidate build/issue-272/resumed/candidate-wheel/epcsaft-0.2.0.dev0-cp313-cp313-linux_x86_64.whl
```

With the candidate wheel installed and both retained state sets present,
`scripts/delivery-evidence.sh` compares the states, runs `tools/time-key-solves`
and executes the required affected and full checks. The timing tool installs
both identified wheels in its own temporary environments and preserves the
checkout environment.

Evidence is numerical verification on the declared local states, without
experimental validation, a global stability certificate or a general state-error
bound. Independent delivered review and investigator promotion remain separate.

The amendment-2 baseline build log and the equivalence `check-affected`/`check-all` logs are at record commit `21d7a112`; `results/amendment-2/change-accounting.json` lists the files as they stood when that amendment was recorded.

## Historical declaration replay (#720)

On `ca719012`, the unchanged producer solves all 15 declared states and passes
its declaration-identity assertion. The generated
[inputs](results/resumed/inputs.json) have SHA-256
`577bd90565f348a5652cab8d4492a2cc646c4412e4121559a56ceed07a4dd06a`,
exactly the stored reference declaration hash, also at audit base `335707d5`.
The archived #696 log reports a mismatch, but its differing declaration file
was not retained in that archive. The exact changed field and cause in that
execution remain unresolved; no declaration was edited or reference replaced.

After verifying that exact identity, the original comparison against the
stored baseline fails its unchanged 1e-11 preservation limit. The largest
relative trace-amount difference is 3.821850521611541e-7 for unloaded-MEA
carbonate; the largest phase-relative difference is 6.259862029612324e-9
for the small vapor amount. [Current states](results/issue-720/issue720-states.csv),
[polish traces](results/issue-720/issue720-trace.csv),
[comparison rows](results/issue-720/physical-state-comparison.csv), and
[maxima](results/issue-720/budget-maxima.csv) retain the numerical evidence.
That replay did not establish the reason for those historical-baseline differences.
The subsequent diagnosis and authorized reference refresh are recorded below.

The producer was called through its existing `run` function with `OUT` set to
`build/issue-720/polish-base/` and an exact copy of the stored `inputs.json`,
so the historical comparison files were never overwritten. The command and logs
are `build/issue-720/polish-base.py`, `polish-base.log`, `compare-polish.py`
and `polish-compare.log`. Copied wheel metadata preserves those original output
paths and hashes. Its installed wheel is the same image as the other #720 runs.

## Current-equation reference refresh (#754)

The [accepted diagnosis](https://github.com/tannerpolley/ePC-SAFT/issues/754#issuecomment-6328215094)
attributes the MEA shift to the intended infinite-dilution reference correction
in [#353](https://github.com/tannerpolley/ePC-SAFT/issues/353), and brackets the
masking trace-inventory correction at [#633](https://github.com/tannerpolley/ePC-SAFT/issues/633).
The #353 attribution is an inference supported by the source change and numerical
probes: it cannot be isolated with the published wheels because the earliest
available wheel already contains the shift. The read-only diagnosis is
`/home/tnnrpolley21/Workspaces/Engineering/ePC-SAFT/build/issue-754/diagnosis.md`,
SHA-256 `3e518992ff71c9c5aa0a7a49a7a4fe60abba282900c42ae460c555787ded9140`.

The original reference remains frozen in
[amendment-2/baseline-states.csv](results/amendment-2/baseline-states.csv),
from reference source `d7739584458aa43a265dc7cc8b6f0e1a097c1480` and wheel
SHA-256 `ae3d960f3b310d88c83468e361c390392a8b58b6e5dafbead662f42e0e7fd8d3`.
Its [wheel record](results/resumed/baseline-wheel.json), declarations, comparisons
and all other historical rows retain their original bytes and meanings.

The unchanged producer at main `8ffd627f0f7e4e68fdf9db634773ff00e058d85e`
generated a fresh reference and a separate repeat using the same installed wheel.
[Qualification provenance](results/issue-754/qualification.json) and the
[reference](results/issue-754/baseline-wheel.json) and
[repeat](results/issue-754/candidate-wheel.json) records identify the wheel,
declarations, states and traces by SHA-256. Each fresh state and trace table was
byte-identical to the #720 numerical owner, so the records reuse those files
instead of keeping duplicate rows. The comparison against the historical
reference also reproduced the #720 rows byte for byte.

Both fresh runs solve all 15 states and meet every requested row tolerance.
The [residual maxima](results/issue-754/residual-maxima.csv),
[physical comparisons](results/issue-754/physical-state-comparison.csv) and
[comparison maxima](results/issue-754/budget-maxima.csv) retain the results;
the notebook computes the before/after table from these files. The original
comparison code passes every assertion with the unchanged **1e-11 comparison
limit**. This is a repeat on the accepted current equations, not a new isolated
experiment with the polish shortcut disabled.

The small-positive-phase discrepancy is a rounding-amplification note, not a
loosened bound: recovering a tiny phase from differences of order-one inventories
amplifies binary64 rounding. The diagnosis supports refreshing the reference;
it supplies no universal forward-error bound. This evidence is local numerical
verification, without physical validation or a global stability claim.

To regenerate with a wheel built from the identified producer commit, use this
checkout's `.venv` and leave `EPCSAFT_HSL_LIBRARY` unset. The output override keeps
the historical directories frozen; it changes no input or solve operation:

```bash
nice -n 10 env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 agent-heavy bash .codex/environments/setup.sh
nice -n 10 env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 agent-heavy bash tools/build-engine.sh --install-wheel
nice -n 10 env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 PYTHONDONTWRITEBYTECODE=1 agent-heavy .venv/bin/python - <<'PY'
from pathlib import Path
import runpy
import sys
sys.path.insert(0, 'analyses/2026-equilibrium-polish-masking/scripts')
import qualify_polish as q
q.OUT = q.ROOT / 'analyses/2026-equilibrium-polish-masking/results/issue-754'
q.OUT.mkdir(exist_ok=True)
wheel = Path('build/environment-wheel/epcsaft-0.2.0.dev0-cp313-cp313-manylinux_2_39_x86_64.whl')
q.run('baseline', wheel)
q.run('candidate', wheel)
runpy.run_path(str(Path(q.__file__).with_name('compare_polish.py')))
PY
```

The producer writes complete tables on replay. Before retaining them, compare
their bytes with the existing numerical owners, remove identical copies, and
update the wheel records' paths to the single retained owner, as in this delivery.
Historical qualification of the shortcut remains tied to its original builds;
independent delivered review and investigator promotion remain separate.

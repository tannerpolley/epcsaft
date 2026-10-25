# Density root benchmark

Shared numerical verification and warm timing for [#478](https://github.com/tannerpolley/ePC-SAFT/issues/478).
The committed [states](inputs/cases.json), [native oracle](results/oracle.csv),
[application values](results/application-reference.json), [comparators](results/competitors.csv)
and [main reference](results/main-0225e03b.csv) are identical across prototype checkouts.
Existing declarations and packets own the inputs; calculations run in the installed Engine.

```bash
python3 analyses/2026-density-root-benchmark/scripts/run.py \
  --python /path/to/checkout/.venv/bin/python --label prototype
python3 analyses/2026-density-root-benchmark/scripts/run.py --compare A.csv B.csv
```

Outputs are `build/density-benchmark/<label>/rows.csv`, `summary.md` and retained
returns. Three fresh processes alternate case order, with nice 10, single-thread
pools, five warm-ups and the declared N per process. Each Tier A/B primary sequence
contains 32 distinct pressures, fixed T and composition, and a native oracle per
state. All states retain the original stable-root count and mechanical branch;
metastable sequences stay between model saturation and the appropriate spinodal.
Six cheap Tier C applications use 32 temperatures from T to T + 0.031 K. Other
applications are labeled **repeated identical input (rewards caching)**. The
same label marks every changing sequence's secondary baseline timing. The full
fit retains the origin start, full controls and observations, one warm-up and
one timed fit per process, as authorized.

The Tier A/B primary metric is the **marginal cost of a pressure request**:
`state(T,P)` minus `state(T,rho)` at the same T, composition and returned pressure-state density.
Both raw call times are reported. This subtraction includes differences in the
complete state requests and can amplify timing noise. FeOs `State` and Clapeyron
`volume` are raw call comparisons with the raw Engine pressure request, with
different returned quantities; they are never compared with the subtraction.
Preparation, checks and trace capture are outside timing. Medians and p10–p90
pool warm calls across three processes; FeOs includes Python dispatch overhead.

Every warm-up and timed return is retained and checked after the timer stops.
Density requests require relative density agreement at most 1e-10, positive pressure slope, the requested branch and frozen stable-root count, and calculated reduced pressure residual at most 1e-12. The current criterion has no additive floor operand. Earlier retained files preserve their historical criteria and numerical-floor values. A sole stable root returns `Unspecified` and count 1, including the absent-vapor request. These checks provide numerical verification of the calculated equations, not physical validation or global phase stability.

Applications compare every
retained field with main at the `time-key-solves` quantity scale, 1e-10; TPD
uses absolute 1e-10 in RT units. Finite stability searches are not certificates.

Work is replayed on fresh models with the same warm-up and changing/repeated
history as the measured mode; every replay return is checked too. Density-solver
spans and their descendant evaluations remain diagnostics. The total counts
**every `evaluate` span in the native trace, independent of its ancestors**;
this trace exposes no derivative order, so the total is unclassified. Stability
also retains its native `energy_evaluations`. Counts are means per primary call;
individual counts remain beside each replay return.
Every Helmholtz or pressure evaluation, by any kernel, must appear as one `evaluate` span in the native trace; prototypes must preserve this requirement.

`density_diagnostics(adaptive=False)` scans 16,384 declared log intervals to the
native packing bound and refines with log tolerance 1e-14. Its crossing detection
and local `refine_crossing` are independent of `solve_density`'s `scan_roots`,
`refine` and `anchor`; it shares `at_density`, the EOS evaluator and the packing
bound. Native refinement bisects density; exact floating pressure zeros may stop
before the bracket reaches the log tolerance. The CSV records stop and bracket
width. The [56 original roots](results/oracle-reproduction.csv) are reproduced;
1,792 roots cover all changing states. Finite sign scans cannot certify completeness
or detect every double root. Cases avoid exact spinodals.

Comparator densities carry their individual oracle differences and a separate
1e-10 accuracy result. Pure-water sigma(T) is evaluated by the canonical native
`segment_sigma` template from the Held packet. FeOs cannot represent C3's sigma
screening choice; neither conversion represents MEA, and Clapeyron cannot return
the sole liquid-side root for the absent-vapor request. Shared-host timing noise,
comparator conventions and any recorded accuracy differences limit selection.
This is numerical verification, without physical validation.

The current density/spinodal producer replays the retained run declaration with
native calculations and requires bitwise reproduction of every canonical state
record. Use the verified named-wheel environment:

```bash
nice -n 10 /path/to/named-wheel/venv/bin/python -B \
  analyses/2026-density-root-benchmark/scripts/pressure_residual.py
```

Add `--check` to check the retained values without calculating. Both commands
leave retained evidence unchanged. The per-state JSON files under
`results/main-8298c602-pressure-residual/` own the raw returns, diagnostic
crossings and endpoint evaluations. The notebook reduces those records directly;
the earlier aggregate return, bracket and nested-grid exports are historical
copies. The original producer hash in `provenance.json` identifies the original
calculation; the PR records the consolidated producer and its reproduction run.
`scripts/check.py` checks the earlier timing/application records and uses the
same density-return mutations. The [research notebook](notebook.qmd) interprets
retained values only.

## Named-release rerun (#680)

The [research notebook](notebook.qmd) reads the new results under
[results/main-8298c602/](results/main-8298c602/). The earlier retained files remain
under their original identities. The new provenance records the downloaded
wheel, its verified hash, the actual installed Engine identity and the declared
producer. Unavailable portions retain their stopping states and causes.
Independent delivered review and owner promotion remain open; the manuscript
is unchanged.

The #683 preparation attempt stopped at the first state because its historical
floor diagnostic was unavailable. Its stopping record and notebook account
remain unchanged. The current producer uses the direct residual criterion above.

## Current producer diagnosis (#750)

[The current diagnosis](results/current-main/diagnosis.md) records the unchanged failure,
the responsible Engine changes, retained comparisons, repaired reruns and remaining
limits. New results are separate from the original retained results.

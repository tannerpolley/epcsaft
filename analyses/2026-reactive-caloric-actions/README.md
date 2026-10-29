# Anchored Gibbs calorics and reference derivatives

The [notebook](notebook.qmd) reads retained numerical evidence for issue #225
slices 1–3 and central-temperature extraction scaling. The manufactured E/I/U
and E/N bubble cases test values, exact actions, reference completion and
limiting admission under unchanged criteria. They supply numerical verification,
not physical calorimetry, reaction chemistry validation, parameter adoption or
global stability. The apparent-Cp build has a separate retained numerical record.

## Apparent-Cp build

`results/slice4/` retains the native comparisons and source/wheel identity for
Equilibrium's apparent isobaric heat capacity, its first actions and exact fitting
columns. Independent normalized H values qualify the first actions; rebuilt
native Cp predictions qualify the fitting columns. The two quadrature identities,
nonreacting reduction, bounded synthetic recovery and enthalpy preservation use
the issue's unchanged criteria. The notebook computes every reported number from
those retained rows. Key-solve work and counterbalanced wall times remain under
the raw build path identified by the provenance file.

From the repository root, after preparing the installed wheel:

```bash
bash analyses/2026-reactive-caloric-actions/scripts/apparent-cp-qualify.sh
.venv/bin/python analyses/2026-reactive-caloric-actions/scripts/apparent_cp_transport.py
```

The first runs native calculations; the second reads
native result rows and renders the comparison figure.
Delivered review returned Supported with no findings (PR #310); the investigator
promoted the result for #228 on 2026-10-09. The manufactured E/I/U
branches and fixed-source recovery do not physically validate reaction chemistry
or establish general coefficient identifiability.

## Current-main replay

`results/main-2026-10-06/` retains the replay after integrating main `d7739584`.
Candidate source `2428781b` preserves the density publication correction,
duplication removals, trace payload construction inside `Scope::observe` and
the separate exact-inventory incipient extension. The binary128 verification
callers use the current parameter record without narrowing epsilon motion.
No runtime calculation or numerical criterion is added by integration.

The native tables retain every F6 selector, affected prescribed calorics,
P4/ionic/Hessian/pressure/volume/F8/F9/Born comparisons, the lower-admission
ladder, the original-order previous-failure row and the binary128 pin.
`row-overlap.tsv` identifies the exact 67 shared row movements from the common
recording; all have fresh independent comparisons. The movement-table owner named in the current provenance
retains all moved operands, causes, methods, criteria and verdicts. The current
replay regenerated that earlier table byte-for-byte, so its numeric rows are
retained once; fresh raw logs establish current-source requalification.
`provenance.json` identifies the exact source, wheels, commands, costs, checks,
raw logs and limits. Independent delivered review remains with the lead.
No notebook promotion for manuscript use is implied.

Raw outputs, native reference sources, baseline source archive, separate
non-editable wheels and the initial compile and stale-baseline failures were
written under `build/current-calorics/main-2026-10-06/` and are archived in
`build/archive/issue-225-slice4.tar.zst`. The exact current-main archive
has its own environment; the timing tool reads only its installed wheel.

From the repository root, the current native producers are:

```bash
bash tools/setup.sh
native_build_dir=$(<build/environment-wheel/meson-build-dir)
"$native_build_dir/properties_checks"
"$native_build_dir/equilibrium_checks" reference-caloric-identities
"$native_build_dir/equilibrium_checks" reactive-caloric
"$native_build_dir/regression_checks" reactive-caloric-F6-addition
```

## Frozen dated reference programs (#720)

The dated `prediction-reference.cpp`, `fit-row-references.cpp` and
`surface-state-references.cpp` under `scripts/main-2026-10-06/` are
**historical, not rerunnable**. Record snapshot: commit
`ca719012c31e8b41edcb98ae82a388117b58d343`; their comparison source is
`2428781b` and the cited tables remain in `results/main-2026-10-06/`.

At #696’s measured commit `335707d59c67cf1f92df12b5d25530bf31d916a6`:

- `supplement-native-prediction-reference` failed with `failed: found`
  (report line 703).
- `supplement-native-fit-row-references` failed with `map::at`
  (report line 704).
- `supplement-native-surface-state-references` failed with `map::at`
  (report line 705).

All three require the historical generated table
`build/current-calorics/main-2026-10-06/native-rows/native-rows.tsv`, which is
not retained in current main. The latter two failures are subsequent row
lookups into that absent table. No replacement values were invented. The
former staging and `references.sh` recipe is withdrawn; the cited comparison
and qualification tables, sources and original wheel identities remain intact.
The historical oracle source was `2428781bd:engine/tests/oracle/native-rows.tsv`
(SHA-256 prefix `610f72bf…`); the replay’s selector list remains at record commit
`21d7a112`. The independent reference sources also include
`surface-properties.cpp` and `surface-equilibrium.cpp`, compiled with
`tools/native-driver`. The historical `classify.cpp` reduced `comparisons.tsv`
to `qualification.tsv`. These method identities do not supply the absent table.
The current caloric producers above retain their separate scope.

The notebook and plot read retained results and do not run Engine calculations.

## Earlier evidence

Each earlier table retains its own source/wheel identity and claim limits:

- `results/diagnosis.csv`, `provenance.json`, `review.md`:
  original reference qualification and arithmetic diagnosis.
- `results/slice2-provenance.json`, `slice2-review.md`: the earlier prescribed
  caloric replay and its separate retention review.
- `results/current-source-provenance.json`, `current-source-review.md`: the scoped Supported prescribed-T/P runtime
  review at `14c32d1af98dbf362eb0f48e34129ee08814b3ee`.
- `results/lower-reference-causes.csv`, `lower-reference-provenance.json`:
  the measured admission boundary, available sibling column and owned cause
  transported to the consuming caloric completion.
- `results/temperature-extraction/`: corrected temperature conditioning,
  both reference orders, native movements and equivalent-work timing.
- `results/main-integration/`: the earlier density/scaling integration,
  independent references and two warmed timing rounds at `f242f778`.

The previous scoped reviews support their inspected snapshots and stated uses.
They do not establish delivered acceptance of the current integrated source.
#341 removed the superseded full comparison tables (`results/comparisons.csv`,
`current-source-comparisons.csv`, `main-integration/` and
`temperature-extraction/comparisons.tsv`), raw timing state tables, earlier
figures with their plotting and collection scripts, and the per-run
reproduction copies; they remain at record commit `21d7a112`.

## Named-release rerun (#680)

The [research notebook](notebook.qmd) reads the new results under
[results/main-8298c602/](results/main-8298c602/). The earlier retained files remain
under their original identities. The new provenance records the downloaded
wheel, its verified hash, the actual installed Engine identity and the declared
producer. Unavailable portions retain their stopping states and causes.
Independent delivered review and owner promotion remain open; the manuscript
is unchanged.

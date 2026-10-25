# Figiel 2025 Figure 9 interior comparison

Current scope is Figure 9's 409 interior mixed-solvent observations at 298.15 K
and 100000 Pa with unchanged `figiel-yu-held-2025/4`. The retained table also
accounts for all 57 excluded endpoints. This is source-limited physical
comparison: current errors, RMSE and maxima are descriptive, with unknown
measurement uncertainty and no physical pass threshold. The legacy-derived
threshold file applies only to its original historical results/population.

Run from the repository root using this checkout's own editable Engine:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 .venv/bin/python analyses/reproduction/2025-figiel/scripts/reproduce.py
cd analyses/reproduction/2025-figiel
QUARTO_PYTHON=../../../.venv/bin/python nice -n 10 quarto render notebook.qmd
mkdir -p ../../../build/issue-789/2025-figiel/render
mv notebook.html ../../../build/issue-789/2025-figiel/render/
```

`results/figure-09-comparison.csv` retains every source file/hash/row identity,
coefficient or cause, descriptive error and numerical diagnostic. Native
reference terminal change must be within 5e-5 and relative pressure closure
within 1e-8, with finite positive coefficient, density and pressure slope.
`results/figure-09-run-identity.json` identifies all source and packet bytes and
the loaded native library. The notebook computes its coverage, descriptive
statistics and comparison plots from retained files only. It neither reruns the
Engine nor applies historical thresholds. Independent delivered review remains
open.

Source: Figiel, Yu and Held (2025), DOI 10.1021/acs.iecr.5c00475, Figure 9;
exact source rows remain under `figures/figure-09/source/`. The observation and
parameter ledgers remain in `inputs/source/`. The source checker is rerunnable:

```bash
.venv/bin/python analyses/reproduction/2025-figiel/scripts/check_figiel_2025_source_packet.py
```

Figures 2–8 remain uncompleted under the retained input/comparison limits;
Figure 1/Tables 1–5 remain source material. The recovered Figure 5 LiBr
extraction uncertainty is preserved for later separate work, and is not
measurement uncertainty. The 36 missing-primary organic NaBr ledger records
belong to the single-solvent figures, not this Figure 9 population. Older
Provider rows remain historical. No literal paper-curve parity, refit or
independent prediction is claimed.

## Historical figure reproduction (does not run on the current wheel)

`figures/figure-0N/results/data/figure-0N-series.csv` hold literature rows and
Provider rows separately. The Provider rows are a historical record: they were
generated at 298.15 K and 1 bar on wheel SHA-256
`4b6385ab4f4a4ba1bba2c69a807e6c105120e2b9460f84fd2c4006e963540901` with the
parameter bundle `data/packets/figiel-2025-validation-parameters/1/parameters`
(its byte-identical leaf copy was removed by #341)
(`figiel-2025-validation-parameters`, version 1). That campaign ended `failed`
against the legacy-derived gates in `inputs/docs/legacy-derived-gates.json`
(10% margin over the legacy model): Figure 4 (RMSE 2.206 > 1.945 kJ/mol, max abs 3.549 > 3.350)
and Figure 7 (RMSE 0.03095 > 0.02252, max abs 0.0994 > 0.0671) exceeded them; Figures 1–3, 5, 6, 8 and 9
passed. The per-figure metrics recomputed from the retained series match that
campaign result exactly. Its full result file, `results/receipts/figiel-2025-full.json`, was removed with
the retired result machinery in commit `9a7fc51a`; read it at `9a7fc51a^`.

`scripts/figiel_2025_full.py`, the method record (grids, carrier salts, panel
conventions, gate metrics), was removed by #341 (record commit `21d7a112`). Its last revision pinned wheel
`48f3a75c9fc16ba71616aa703b526f41c2dcf89a7e00eebe23f75fcb8fa24594`, reads packet
`figiel-2025-reference-electrolytes/2`, and called the retired
`epcsaft.unit_registry` and `Mixture.reference(solvent_x=…, ion_component_ids=…,
reference_convention=…)` keywords, so it did not run on the current wheel and did
not produce the retained series. `scripts/figiel_2025_miac.py` (removed at the same commit) was historical in the
same way; its retained aqueous result `results/data/figiel-2025-miac.csv` came from
wheel `758327312c1a7f84a4e836e37b6798b1a1452506994da84429095247c3ce1474`.
The other `results/*.json` files record their own wheel identities.

## Tables

Tables 1–5 are source transcriptions, not Provider results. Each has one owner,
`results/tables/table-00N/results/data/table-00N.csv`; no script writes them.
The former specialized aqueous `k_ij`, solvation, and permittivity parameter-fit
replays have been removed; those observation families are outside the generic
parameter-fitting boundary.

The figure plotting helpers live in `scripts/shared/` (moved from `inputs/shared/` by #341).
Figure 3 reuses `figures/figure-02/source/dielc-salts-in-water.csv`; #341 removed its byte-identical copy.

# Perdomo strong-electrolyte phase stability

Current native comparison at 298.15 K: NaCl–water solvent-only-vapor bubbles
for two explicit dielectric declarations, Figure 1a topology stress, Figure 3
LiCl–water–1-butanol feeds, and Table 5 liquid splits and pressure sensitivity.
The [notebook](notebook.qmd) reads retained numbers and distinguishes numerical
closure, comparison between different models, and unresolved declared branches.
Independent delivered acceptance and manuscript promotion remain pending.

Run from this checkout's root with its own editable Engine environment:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 env EPCSAFT_BUILD_JOBS=6 \
  .venv/bin/python analyses/reproduction/2025-perdomo/scripts/reproduce.py
```

`--case <exact ID from inputs/current/cases.json>` runs a single declaration;
`--read-results` refreshes comparison CSVs from retained records without solving.
Each calculated case binds its declaration, selected parameter file, criteria and
applicable diagnostic seed to its producer and Engine identity. Comparison refresh
preserves those identities and refuses execution-input mismatches; a partial run
updates only the executed case, preserving other cases’ Engine identities.
The producer preserves every case, achieved native rows, full search diagnostics,
bounded native diagnosis attempts, errors and Engine identity under
[`results/current/`](results/current/summary.csv).

Render the notebook as a standalone presentation without changing the shared
Quarto site configuration or rerunning thermodynamics:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 env EPCSAFT_BUILD_JOBS=6 \
  bash analyses/reproduction/2025-perdomo/scripts/render.sh
```

The standalone [rendered notebook](../../../build/issue-729/render/notebook.html) is
retained for review. Root-site registration is a lead-owned gate because this
build excludes changes to `analyses/_quarto.yml`.

Perdomo's paper uses SAFT-γ Mie; these calculations use ePC-SAFT. Their source
comparisons are numerical verification between different models, with no measured
data or physical-validation claim. The source tables retain the temperature,
nominal molality, vapor-volume, pressure-header and tie-line-count discrepancies.
Source uncertainty remains unknown. The pressure display, refinement, balance,
charge, chemical-potential and objective criteria are frozen in
[`inputs/current/criteria.json`](inputs/current/criteria.json).

Bubble vapor contains water only. The separate explicit charged-vapor inputs
remain pending [Engine #771](https://github.com/tannerpolley/ePC-SAFT/issues/771).
No trace vapor-ion prediction or salt chemical-potential equality is claimed.

**Conditional LiCl diagnostic: explicit zero assumptions for the three
Li-containing diameter corrections; historical inputs retained; comparison
between different models.** The complete model and provenance are in
[`inputs/current/licl-conditional.json`](inputs/current/licl-conditional.json)
and the [selected input table](inputs/current/parameter-selection.csv).
Finite discovery and bounded declared-start attempts do not prove global stability
or physical nonexistence of an unresolved second liquid. No refit or reusable
parameter adoption is made.

The original bundle, source ledger, Figure 1a sweep, Figure 3 fits and rendered
comparisons remain historical evidence with their original bytes and software
identities. They grant no current parameter validation or executable replay.
The Figure 3 historical transfer-scale and predictive limits are described
[beside its sources](figures/figure-03/README.md).

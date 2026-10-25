# Historical MEA speciation parameter regression

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The pre-D042
parameter-fitting runner and its GREPE interfaces were removed; retained fit inputs and
the declarations imported by later studies remain available.

This directory preserves the completed pre-D-042 Regression campaign, its
inputs, generated values, and review record. It is immutable provenance, not a
current replay target. The former IPOPT/`regression.fit` runner was removed by
the native owner-neutral cutover; no command here should be used to recreate
that lifecycle.

This fixed campaign fits three shared ePC-SAFT parameters to 131 admitted
speciation measurements at 44 homogeneous MEA states. It uses the installed
Engine wheel, public GREPE and Regression APIs, and one simultaneous sparse
full-space IPOPT solve. The parameter and observation inputs and complete
result are retained under `inputs/` and `results/`.

The completed installed-wheel campaign evaluated all 131 targets without a
drop or rejection using the simultaneous sparse full-space backend. It fitted
`k_ij = 0.10159970424983134`, `sigma_MEAH+ = 4.365905218305267 angstrom`, and
`sigma_MEACOO- = 5.1809945813168845 angstrom` with objective
`0.0343599102370016` in 26 IPOPT iterations and 95.89 seconds.

## Historical CO2 partial-pressure comparison

`scripts/run_pressure_validation.py`, a historical method record removed by #341
(record commit `21d7a112`), used
retired `ActiveParameterSet`, `general_reactive_equilibrium_problem_from_mapping`
and `equilibrium.solve` interfaces. Its retained result identifies wheel SHA-256
`405e79b7ea100e9bafaced88acfa6b1a7b7819d91d3a4bb176c651417f2331bd`.
This rejected speciation-fitted candidate supports no current capability claim;
current MEA qualification is #61 and
[`2026-mea-recovered-state-actions`](../2026-mea-recovered-state-actions/README.md).
The plot command below reads the historical table without Engine evaluation.

The independent pressure check uses 58 source-backed 30 wt% MEA observations
from Hilliard (2008) and Jou et al. (1995) at 40, 60, and 80 °C:

```bash
.venv/bin/python analyses/2026-mea-grepe-regression/scripts/plot_pressure_validation.py
```

Each in-domain row is a fixed-pressure, single-finite-phase public GREPE
request. The predicted liquid CO2 fugacity is compared with the measured CO2
partial pressure; the specialized reactive-bubble calling route is not used. The
four-process replay took 29.77 seconds and evaluated 55 rows with zero solver
or certificate failures. Three Jou rows at 500, 600, and 1000 kPa remain in
the table as `outside_parameter_domain` because the fitted parameter document
declares a 300 kPa upper pressure bound.

The pressure result is a failed predictive check, not another fit: median
absolute `log10(model / observed)` is 1.024 and the maximum is 1.211. The
revised speciation-fitted set therefore overpredicts pressure by roughly one
order of magnitude and must not be described as pressure-validated.

Hilliard and Jou are observation sources, not separate models. The plot uses
points for their observations and marker-free lines for the one Engine model
evaluated at each source state. The complete interpretation, retained-input
links, and source-by-source recalculation are in
[`notebook.qmd`](notebook.qmd). The fit input
supplies reaction constants at 20, 40, 60, and
80 °C, so this check does not extrapolate them to the available 100 or 120 °C
observations.

The local observation snapshot is materialized from
`Canonical_VLE_Observations.csv` (SHA-256
`9e7d9ba5fead8bfa83a311dad341e3e2e8df1806d5249642a23562e99a72cb73`)
and `pco2_metrology_manifest.csv` (SHA-256
`0d14803873a60534ec5d7df382cfbd0ae03e4aaeba68bb5d54be7e4def8397cc`).

## Historical fitted parameters

`mea-epcsaft-parameters.zip` and the retained parameter document preserve the
nine-species candidate used in this comparison. They are not current MEA
adoption instructions. The historical problem supplied feed, reactions and
source-standard-state correlations through `ChemicalEquilibriumProblem` and
`GeneralReactiveEquilibriumProblem`; those interfaces have been removed.
The producer projected rounded source molar masses onto its declared
reaction-conservation subspace before constructing that problem. These
historical transformations remain in the method record.

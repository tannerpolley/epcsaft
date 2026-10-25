# NaCl activity in ethanol–water

Current native Figiel predictions are compared with the complete Esteso Table I
activity grid. The [notebook](notebook.qmd) reads the retained results, shows every
solvent series and reports each printed-precision outcome. Disagreements remain
visible with undetermined physical causes; numerical resolution is assessed
separately. Tables II–V, reference coefficients and measured voltages are source
records only. Esteso's fitted Debye–Hückel/Pitzer, neutral-pair and forward cell
voltage calculations are outside this reproduction's scope.

Run from this checkout after its normal environment setup:

```bash
nice -n 10 env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 .venv/bin/python analyses/reproduction/1989-esteso/scripts/reproduce.py
```

`--compare-only` refreshes the comparison from retained native records after
checking their input hashes; it preserves their original execution identity.
The producer uses only `figiel-2025-reference-electrolytes/4`, formula-unit
molality per kg mixed solvent, salt-free solvent composition, 298.15 K and
100000 Pa. No parameter is refitted. Native C++ owns reference construction,
density closure and activity conversion. The actual activity-call traces supply
convergence evidence; the binding does not expose its consumed State.

The absolute printed-precision criterion is 0.0005 and the separate numerical
allocation is 0.000005. Neither describes experimental uncertainty. See the
[source ledger](inputs/source/source-ledger.json),
[verified Table I](inputs/source/esteso-1989-table-I.csv),
[native results](results/native-results.json.gz) and
[row comparisons](results/comparisons.csv). The notebook distinguishes physical
comparison from numerical verification and preserves the inconsistent source
reference statements without recalibration.

Source: Esteso et al. (1989), *Activity coefficients for NaCl in ethanol-water
mixtures at 25°C*, [DOI](https://doi.org/10.1007/BF00652989).
Model: Figiel, Yu and Held (2025), [DOI](https://doi.org/10.1021/acs.iecr.5c00475),
[accepted packet](../../../data/packets/figiel-2025-reference-electrolytes/4/).
The earlier comparison's producer was not retained. Its values and figures
remain at record commit `ca719012c31e8b41edcb98ae82a388117b58d343`; earlier
method text remains at `21d7a112`. The legacy source input is preserved, and its
values are checked against the full transcription. Historical results supply
no live prediction or qualification claim.

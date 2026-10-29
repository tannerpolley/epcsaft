# Measured NaCl dilution enthalpy

The [notebook](notebook.qmd) compares the declared Held strategy-2 `/2` caloric predictions with the Messikomer–Wood measurements reproduced by Partanen–Partanen.
It reports every selected replicate, signed nominal errors, nonnegative interval gaps, finite-concentration apparent heat-capacity differences and exact thermal actions.
Disagreement means exceeding the study-chosen 0.003RT benchmark, without an uncertainty-complete experimental rejection claim.
The [source notes](inputs/sources.md) preserve the partly audited caloric fitting history, unverified historical permittivity identity, secondary transcription and higher-temperature extrapolation.

Calculation and presentation are separate:

```bash
EPCSAFT_BUILD_JOBS=6 nice -n 10 .venv/bin/python analyses/2026-nacl-dilution-enthalpy/scripts/run.py
nice -n 10 .venv/bin/python analyses/2026-nacl-dilution-enthalpy/scripts/present.py
quarto render analyses/2026-nacl-dilution-enthalpy/notebook.qmd --no-execute
```

Use this checkout's prepared editable Engine environment.
The [run identity](results/run-identity.json) records exact native library, input and producer hashes.
The [comparison CSV](results/comparison.csv), [chord intervals](results/chords.json), [native readouts](results/native-states.json) and [failures](results/failures.json) own the retained numerical values.
The notebook reads generated tables and a figure; rendering never runs the model.
Independent delivered review and investigator promotion remain outstanding.

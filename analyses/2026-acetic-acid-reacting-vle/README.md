# Aqueous acetic acid: coupled bubble pressure

This leaf compares native reacting bubble pressures with the water-majority
Bernatová–Aim–Wichterle measurements at 323.15 K, using fixed published parameters.
The [notebook](notebook.qmd) reports pressure errors, chemistry and permittivity
controls, numerical verification and conditional input-interval comparisons.
Historical calibration overlap is unresolved. This case does not test speciation.

The exact accepted-main wheel and its installed environment are retained in the
main checkout under `build/issue-856/build-r1/`; [executable.json](results/executable.json)
identifies that environment, wheel hash, clean Engine source and loaded libraries.
The native reference diagnostic is retained beside that environment; its exact
compiler command is `build/issue-856/build-r1/reference-build-command.json`
in the main checkout and links the installed wheel's native library.
From this checkout, regenerate with that qualified environment:

```bash
nice -n 10 python3 analyses/2026-acetic-acid-reacting-vle/scripts/extract_measurements.py
nice -n 10 /home/tnnrpolley21/Workspaces/Engineering/ePC-SAFT/build/issue-856/build-r1/environment/bin/python analyses/2026-acetic-acid-reacting-vle/scripts/run.py
nice -n 10 /home/tnnrpolley21/Workspaces/Engineering/ePC-SAFT/build/issue-856/build-r1/environment/bin/python analyses/2026-acetic-acid-reacting-vle/scripts/reference_refinement.py
nice -n 10 python3 analyses/2026-acetic-acid-reacting-vle/scripts/summarize.py
nice -n 10 .venv/bin/python analyses/2026-acetic-acid-reacting-vle/scripts/present.py
```

Rendering reads retained displays and never runs thermodynamic calculations:

```bash
cd analyses
nice -n 10 quarto render 2026-acetic-acid-reacting-vle/notebook.qmd --no-execute
```

[Source notes](inputs/source-notes.md) record E1–E5, original locators,
activity bases and source-access limits. Harned–Ehlers was added to Zotero
as metadata only through a nonbrowser operation. The original PDF page could
not be retrieved; its indexed primary text supplies the retained table check.
No local fit or Engine runtime change is part of this analysis.

The notebook also retains the pressure diagnosis: pure-endpoint errors, an exact
endpoint-line decomposition, bounded interaction and association sensitivities,
and ranked explanations with their source and numerical limits.
The interaction sensitivities do not identify a physical cause or a correction;
pure-parameter influence through mixture nonideality remains unresolved.
Sources and extracts are in `inputs/source/`, numerical outputs and generated
tables in `results/pressure-diagnosis/`.
The earlier executable identity is retained alongside the fresh-run identity.

Prepare this checkout's own `.venv` using the root README Checks, then run from
this study directory:

```bash
EPCSAFT_BUILD_JOBS=6 nice -n 10 ../../.venv/bin/python scripts/probe.py
EPCSAFT_BUILD_JOBS=6 nice -n 10 ../../.venv/bin/python scripts/pure_neutral_controls.py
nice -n 10 ../../.venv/bin/python scripts/decompose.py
```

`probe.py` performs 39 native equilibrium solves and 12 direct state evaluations;
`pure_neutral_controls.py` separately checks the two pure endpoints and three
neutral bubble states. Both require this checkout's loaded Engine.
`decompose.py` reads retained values, checks pressure-decomposition closure and
generates the diagnosis JSON, CSV and notebook tables without model execution.

# Optional Gross (2005) quadrupole–quadrupole term: qualification (#574)

**Question.** Does the quadrupole term of ADR 0014 equal the Gross (2005)
equations, and with the paper's own parameters does it reproduce the paper's
vapor-pressure improvement for CO2 and benzene?

**Answer.** Yes, within the limits in the [notebook](notebook.qmd).
- **Code identity.** The term, saturation pressures and the CO2/N2 bubble point agree
  with Clapeyron 0.6.29 `QPCPSAFT` far inside the 1e-10 (term) and 1e-9 (solved)
  criteria. FeOs 0.10.1 agrees within its 1e-6 criterion where every m ≤ 2.
- **CO2.** With the quadrupole, the vapor-pressure AAD against Span–Wagner over
  220–300 K is about 0.42 %, against 3.0 % for the nonpolar Gross–Sadowski 2001 set.
- **Benzene.** The quadrupole gives about 0.28 % against DIPPR-101 over 290–540 K.
  That needs the uncapped chain length: the dipole's m ≤ 2 cap would give 3.0 %.

The exact values are in [summary.json](results/summary.json) and
[engine-rows.csv](results/engine-rows.csv); the figure is
[saturation-deviation.png](figures/saturation-deviation.png).

## Reproduce

```bash
CACHE=~/.cache/epcsaft-comparators/222
JULIA_DEPOT_PATH=$CACHE/depot $CACHE/julia-1.13.1/bin/julia --startup-file=no --project=analyses/2026-competitor-cross-check/inputs analyses/2026-polar-quadrupole/scripts/clapeyron_reference.jl
$CACHE/feos/bin/python -I analyses/2026-polar-quadrupole/scripts/feos_reference.py
.venv/bin/python analyses/2026-polar-quadrupole/scripts/compare.py
.venv/bin/python analyses/2026-polar-quadrupole/scripts/plot.py
```

The Clapeyron and FeOs rows feed `engine/tests/eos/test_quadrupolar.py`.

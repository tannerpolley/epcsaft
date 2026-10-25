# HDA two-flash ePC-SAFT inputs and Engine-side checks (#73)

[Issue #73](https://github.com/tannerpolley/ePC-SAFT/issues/73) asks whether the IDAES two-flash
hydrodealkylation (HDA) simulation and optimization can take four-component (H2, CH4, benzene,
toluene) phase properties with exact first and second derivatives from one Engine wheel, while
IDAES keeps its own `SmoothVLE`, bubble/dew and fugacity-equality rows. This leaf holds the
Engine-side evidence: the input checks, the benzene/toluene binary path and the pair-by-state
domain matrix. The IDAES model and its campaign live downstream (issue #73, slice 3 design, D1);
`results/campaign/` holds a hashed copy of the campaign's files and [the notebook](notebook.qmd)
reports it. Every full-package and recycle-gas state carries hydrogen and methane, so at least the
hydrogen/methane pair, fitted to cryogenic data, is extrapolated there; the F102 states carry only
benzene and toluene and lie within that pair's source temperatures, with its pressure range unstated
(`results/campaign/auxiliary_domains.csv`, `results/domain_matrix.csv`). Nothing here is a predictive
process claim.

## Inputs

- `data/packets/hda-pc-saft/1`: one nonassociating Gross-Sadowski PC-SAFT family. H2 is Ghosh,
  Chapman and French (2003) Table 1; CH4, benzene and toluene are the full-precision Gross-Sadowski
  records. Its `parameters/parameters.json` holds all six k_ij with their sources.
- `data/packets/nasa-glenn-2002-hda-species/1`: NASA/TP-2002-211556 gas NASA9 records at 1 bar.
  `hda_inputs.thermochemistry` rescales every coefficient by 8.314510 / R_engine, so the source
  cp, h and s are unchanged. Formation-based stream enthalpies carry the reaction enthalpy, so the
  tutorial's separate reaction heat is not used.
- `data/packets/hydrogen-binary-vle/2`: the transcribed H2/benzene and H2/CH4 rows, each with its
  source class.
- `inputs/hydrogen_partitions.csv`: this analysis's training, holdout and stress partition and
  the fixed engineering discrepancy scales `a_*` for each packet row.
- `inputs/yu_rounded_pure.csv`: Yu et al. (2026) Table 1 printed benzene and toluene pure values.

### Partitions and discrepancy scales

The partitions follow the issue's immutable row rules:

- Park 323.2 and 423.2 K rows are training, and its 373.2 K rows are holdout.
- Benham-Katz and Sagara ordinary rows with 0.01 <= x_H2 < 1 and 0 < y_H2 < 1 are training.
- Hong-Kobayashi ordinary rows are holdout.
- Pure endpoints, critical markers and the six ordinary rows with x_H2 < 0.01 are stress.

The scales are not Gaussian standard uncertainties. Each comes from the literal accuracy wording
quoted in the packet's `provenance/source-origin.toml`; the 0.0005 absolute composition floor is
the issue's engineering floor, not a source number.

- Park: a_T = 0.1 K ("0.1 K in temperature"), a_P = 60000 Pa ("approximately 0.06 MPa"), and
  a_x = 0.001 ("less than 0.001 in composition").
- Benham-Katz: a_T = 5/9 K (the ±1.0°F low-temperature bound), a_P = 0.005 P ("±0.5% of the true
  pressure"), and a_z = max(0.005 z, 0.0005) for z = x_H2, y_H2 ("±½% of the component composition").
- Sagara: a_T = 0.1 K (±0.1°C), a_P = 10132.5 Pa (±0.1 atm), and a_z = max(0.01 z, 0.0005)
  ("within 1 per cent").
- Hong-Kobayashi: a_T = 0.01 K ("better than 0.01 °C"). a_P is 0.001 times the smallest of the
  500, 2000 and 6000 psia gauge ranges containing P ("0.1% of the full scale reading"). a_z is
  max(0.005 z, 0.0005), or max(0.083 z, 0.0005) when z <= 0.0051 (the 8.3% error at the lowest
  calibrated hydrogen concentration).

## Checks

The three Engine scripts read the packets through `scripts/hda_inputs.py`, declare their inputs in
`declarations()`, and leave every density, state and equilibrium calculation to the Engine.
`domain_matrix.py` imports no Engine and is carried unchanged from `acf63fd2`.

```bash
.venv/bin/python analyses/2026-hda-two-flash/scripts/nasa9_reference_check.py
.venv/bin/python analyses/2026-hda-two-flash/scripts/rounding_check.py
.venv/bin/python analyses/2026-hda-two-flash/scripts/binary_path.py
.venv/bin/python analyses/2026-hda-two-flash/scripts/domain_matrix.py
```

The three Engine scripts exit nonzero when their falsifier fails, and every script writes its retained
result:

- `nasa9_reference_check.py` writes `results/nasa9_reference_check.{csv,json}`. It fails if any of
  the issue's 16 rows (four species at 298.15, 600, 999.999 and 1000.001 K) misses the table by
  more than 1e-6 J/(mol K) in cp or s or 1e-3 J/mol in h, on the independent evaluator or the
  Engine ideal-gas path. It also fails if Delta_rH(298.15 K) of toluene + H2 -> benzene + CH4
  misses -41890 J/mol by more than 1e-2 J/mol.
- `rounding_check.py` writes `results/rounding_check.{csv,json}`. At equimolar benzene/toluene,
  105 kPa and 350, 370 and 400 K, it fails if density or total enthalpy differ by 1e-4 relative or
  more, or bubble or dew by 0.01 K or more. It also fails if the largest ln(phi_i) difference leaves
  the accepted band [4.75e-4, 4.85e-4).
- `binary_path.py` writes `results/binary_path.{csv,json}`. On the 5 K grid from 300 to 450 K plus
  the two crossings, it fails if any branch lacks two stable roots, if a crossing's recomputed
  fugacity residual exceeds 1e-8, or if a crossing moves more than 1e-6 K from the retained
  bubble and dew.
- `domain_matrix.py` writes `results/domain_matrix.csv`: the below/within/above classification of
  each pair against the feed, heater/reactor, F101 and F102 states. It overwrites the table without
  checking it; the separate check is a byte comparison with the `acf63fd2` table
  (`git show acf63fd2:analyses/2026-hda-two-flash/results/domain_matrix.csv | cmp - results/domain_matrix.csv`), which passed for this
  delivery.

The binary path's central difference of h only verifies the exact cp. The rounding check fails
the issue's ln(phi_i) limit on the liquid branch. The owner adopted the full-precision Gross set
with k_ij = -0.003 on 2026-09-26 and accepted that deviation; no refit followed. Each JSON
records the Engine it ran on (`epcsaft.provenance()`).

## IDAES campaign (slice 3b)

The downstream branch `idaes/hda-epcsaft-73` of the IDAES examples fork builds the two-flash
flowsheet with three property packages joined by IDAES `Translator` blocks (owner decision,
2026-10-13): the full four-component package for the feeds, M101, H101, R101 and F101; a
vapor-only package for the recycle gas (S101, C101, P103); and a benzene/toluene package for
F102, P101 and P102. Every density, enthalpy and ln(phi_i) comes from the installed wheel's
`epcsaft_pressure_property` callback. The `to_aromatics` inlet's H2 and CH4 overall fractions
are fixed at zero and their two defining rows dropped (decision 1, 2026-10-13). The campaign runs
the binary path, the fixed simulation and the three fixed starts at smoothing widths 1, 0.1 and
0.01, and decides each criterion per attempt. Each start is square-solved at its fixed decisions;
when its F102 temperature lies outside F102's bubble-dew band, that coordinate moves to the band's
midpoint and is square-solved again before release (decision 2, 2026-10-13). From the downstream checkout root, with the pinned
non-editable main wheel installed:

```bash
PYTHONPATH=. python -m idaes_examples.mod.hda.hda_epcsaft_campaign
```

`scripts/collect_campaign.py` copies `hda_epcsaft_results/campaign/` byte for byte into
`results/campaign/` and adds table logic only: `auxiliary_domains.csv` (A10 classes of every
audited row), `near_active_rows.csv` (each near-active KKT row and its identity) and
`measure_reference.csv` (the 3b-measure attempts the reformulated optimum is compared with). It
also copies the ideal-package control's summary (`ideal_control.csv`; diagnostic, not acceptance;
produced by `build/73-s3b/ideal-control/control.py` in the main checkout).
`provenance.json` records each copied file's SHA-256, the downstream provenance record and the
packet locks. It recomputes no acceptance:

```bash
python3 analyses/2026-hda-two-flash/scripts/collect_campaign.py <downstream>/hda_epcsaft_results/campaign \
  <3b-measure s4a dir> <ideal-control summary.csv>
```

The diagnostic probes of the failed attempts (trajectories, width and bound isolation, Ipopt's
multipliers) are not retained here; the notebook states their outcomes. The notebook reads only
these retained files; it runs no Engine or IDAES calculation. Render it
into `build/` with Quarto (`nbclient`, `ipykernel` and `pyyaml` in the checkout's `.venv`):

```bash
mkdir -p build/hda-two-flash-notebook
QUARTO_PYTHON="$PWD/.venv/bin/python" quarto render analyses/2026-hda-two-flash/notebook.qmd --execute \
  --to html --output - > build/hda-two-flash-notebook/notebook.html
```

## Hydrogen-pair fit (historical)

`results/hydrogen_pair_fit.json`, `results/hydrogen_pair_rows.csv`,
`results/hydrogen_pair_profile.csv` and `figures/hydrogen_pair_fit.png` are carried byte-identical
from branch `issue-73-hda-two-flash` (`acf63fd2`) and are not rerun (slice 3 design, D2). They came
from the native extension `_core*.so` with SHA-256 (`core_sha256` in the fit JSON)
`54f51f98206584f0eed0123c426c3cf4119fd94efd2477470360a6c24ffab908` and
`hydrogen-binary-vle/1` rows SHA-256
`fab15358e4fb7e1bd5f82cd4871e6aee341d853bf6123c25644289ec84443f01`, which carried the partitions and
scales inline. Version 2 joined with `inputs/hydrogen_partitions.csv` reproduces those rows
string-identically. The only exception is the source class of the six x_H2 < 0.01 rows, which
was `dilute_below_0.01` and is now `ordinary`, still in the stress partition. The fitted
H2/benzene and H2/CH4 k_ij in `hda-pc-saft/1` are that fit's values. The H2/CH4 data are
cryogenic and the H2/benzene data stop at 423.2 K, so hotter HDA states are extrapolations.

## Limits

- Density roots come from a deterministic finite scan, not an exhaustive certificate. Branch
  identity uses the lowest and highest mechanically stable roots.
- The NASA9 check tests record transcription and the NASA9 forms, not calorimetric data.
- No global optimum or global phase-stability claim.

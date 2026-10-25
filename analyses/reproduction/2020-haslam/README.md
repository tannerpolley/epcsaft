# ePC-SAFT comparison with Hamer–Wu; Haslam Table 8 values are not recomputed

Current native ePC-SAFT predictions are compared with Hamer–Wu's evaluated
osmotic and molality-scale mean ionic coefficients. Haslam's printed Table 8
deviations appear alongside these physical comparisons. The [research
notebook](notebook.qmd) reports deviations, rounding diagnostics, reference
convergence and the source-average arithmetic audit. Independent delivered
scientific review and investigator promotion remain open.

The physical comparison uses existing Figiel parameters without refitting,
water with fully dissociated monovalent ions, formula-unit molality in mol/kg water,
298.15 K and 101000 Pa. This is a declared 1% pressure extrapolation from the
packets' 1-bar domain; its pressure error has not been bounded. Mean ionic activity
coefficients use the molality scale. Osmotic coefficients use pure liquid water
at the same temperature and pressure as solvent reference.

The fixed printed-rounding diagnostic is 0.0005 absolute for each dimensionless
coefficient. The native reference-limit diagnostic must be at most 5×10⁻⁵.
Neither tolerance is experimental uncertainty. Reference diagnostics and effective activity settings are retained in the run JSON.
Even after admission, the bound activity result does not expose its per-point
density residual or stopping reason; no second density solve supplies them.

Run from this checkout, after its own environment setup:

```bash
nice -n 10 env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 \
  .venv/bin/python analyses/reproduction/2020-haslam/scripts/haslam_2020_table8.py
```

The command retains predictions and any unavailable values with their causes.
It never supplies a missing model choice. Plot the retained values without
rerunning the Engine:

```bash
nice -n 10 .venv/bin/python analyses/reproduction/2020-haslam/scripts/plot_haslam_2020_table8.py
```

The source-only check is:

```bash
.venv/bin/python analyses/reproduction/2020-haslam/scripts/check_haslam_2020_source_packet.py
```

[Point records](results/data/haslam-2020-table8.csv),
[salt/property accounting and Haslam's printed deviations](results/data/haslam-2020-table8-aad-comparison.csv)
and [run identity, diagnostics and tolerances](results/data/haslam-2020-table8.json)
retain the comparison. The
[corrected source ledger](inputs/source/haslam-2020-table8-source-ledger.yaml)
links the accepted transcriptions. Haslam's parameters are source records only.
Hamer–Wu supplies evaluated, thermodynamically consistent reference values;
they are not independent raw measurements of both properties at every row.
No held-out-data or physical-validation claim is made.

Haslam's SAFT-γ Mie/mean spherical approximation predictions are **not recomputed:
model absent from the Engine**. Exact Table 8 deviations additionally require
unresolved source-row selections. Matching counts establish no exact membership.
Haslam's printed deviations remain descriptive, with no agreement tolerance
against a different model. The printed osmotic average has unresolved aggregation;
it must not be silently replaced by a mean inferred from the saltwise entries.

The accepted design selected superseded packet versions /2. This run uses
`figiel-2025-reference-electrolytes/4` and `figiel-2025-aqueous-iodides/3` by the
owner's round-2 instruction. Both supply `ionic_dispersion`; reference /4 also
declares the Wolbach–Sandler cross-association size factor. No packet was edited
or refitted. Neither current packet contains a source-backed Li⁺/I⁻ pair, so
LiI remains unevaluated.

Historical values and figures are preserved under [results/historical](results/historical/).
They came from EOS commit `53b0b74bf5de86d459c49eb2cc7e75b7cfd4ce3e`, wheel SHA-256
`961156a641435e33746a65d07d97e07a40243e0cfdf870932bd79961783afa96`.
The identity is retained in Git at
`9a7fc51a^:analyses/reproduction/2020-haslam/results/receipts/haslam-2020-table8.json`.
Their original source-membership labels are superseded by the corrected ledger;
they support no current Engine capability or exact Table 8 reproduction claim.

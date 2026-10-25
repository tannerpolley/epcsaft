# MEA parameter estimation

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The scorecard
records a blocked reactive fit and historical physical comparisons; the CO2/water
producer now refuses the required anchored vapor root, as detailed below.

The sole human-facing scientific scorecard and decision record is
[results/analysis/README.md](results/analysis/README.md). Bounded commands are
listed there; retained failures are immutable and are not replay defaults.

The 123-state, 210-target fit observations live once in
[`inputs/observations.json`](inputs/observations.json); candidate parameter and
active-set differences are recorded in [`inputs/fit-candidates.json`](inputs/fit-candidates.json).

## Historical physical CO2/water producer (#720)

`scripts/physical_co2_water_run.py` is **historical, not rerunnable**.
At #696’s measured commit `335707d59c67cf1f92df12b5d25530bf31d916a6`
(report line 712), row `kiepe-313.20-13` raised `RuntimeError` after a returned
`SUCCESS`: phase vapor had no locally unique anchored root. The finite search
refused the vapor branch because its sole stable root was
`17505.97517342151 mol/m3`, above the negative-slope sample
`9740.5196072389117 mol/m3`. This is a local admission failure, not proof that
a physical vapor phase cannot exist. All original physical-comparison values
remain unchanged.

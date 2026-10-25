# Khudaida Figure 2 comparison

The current reproduction covers Figure 2 only: eight source feed declarations
at 293.15 K and 100000 Pa with nominal 5 wt% NaCl, unchanged packet
`khudaida-2026-figure-2-electrolyte-lle/7`, and the exact reversed Table 3 pairing.
All auxiliary feed mappings and actual declared two-liquid states satisfy the
accepted numerical criteria. The computed [notebook](notebook.qmd) compares
sixteen phase rows and all sixty-four formula coordinates. The grand mean
absolute deviation meets the approved 0.0161 paper-correlation-scale limit;
separate standard-uncertainty comparisons do not establish agreement within
measurement uncertainty. Independent delivered review is pending.

[Native results](results/figure-2-diagnostics.json) retain the original source
feed, pairing, one-mole explicit-species feed, auxiliary and actual states,
raw/classified equation rows, projections and their causes. Every projection
uses validated `project_prediction` on its corresponding admitted result.
Exact state preservation and unchanged solve-span counts verify this reuse.
Formula ratios use sodium as one NaCl formula unit with a separate chloride
equality check; salt-free ratios use the three-solvent denominator. The
[comparison table](results/figure-2-comparison.csv) is computed from these
retained results and source measurements when rendering the notebook.
[Run identity](results/figure-2-run-identity.json) binds packet, source, native
library and analysis reader hashes. The native reader is confined to this leaf.

Run from the repository root with this checkout's editable Engine and notebook
dependencies in its own `.venv`:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 bash analyses/reproduction/2026-khudaida/scripts/reproduce.sh
cd analyses/reproduction/2026-khudaida
QUARTO_PYTHON=../../../.venv/bin/python nice -n 10 quarto render notebook.qmd
mkdir -p ../../../build/issue-789/2026-khudaida/render
mv notebook.html ../../../build/issue-789/2026-khudaida/render/
```

The auxiliary one-liquid solve is a native normalization/readout vehicle,
without a physical single-phase stability claim. Actual solves reuse the
source-independent composition guesses of the declared settler branch.
Measured endpoints never seed a solve. This is a fixed-model physical
comparison on a local branch; global stability and the causes of individual
coordinate discrepancies remain unestablished. Feed traces are manufactured
overall declarations, not measured inventories. Paper traces retain unresolved
pixel uncertainty and receive no new threshold.

Figures 3–7 remain unexecuted pending source-backed temperature/salt-domain
and comparison decisions. All thirty-nine original feed identities remain
retained. The older split-wheel HELD2 evidence and plotted outputs keep their
historical identities; they are not this run or a current discovery claim.
Their record snapshot is `ca719012c31e8b41edcb98ae82a388117b58d343`, and the
original packet materializations remain under `data/packets/`. This delivery
runs neither that historical campaign nor an IDAES process model.

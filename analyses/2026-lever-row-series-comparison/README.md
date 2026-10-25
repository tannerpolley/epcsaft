# Lever row-series comparison

[The notebook](notebook.qmd) qualifies the correction of the ionic packing-diameter
volume chain in `CompiledSystem::lever_jacobian`. The existing row-series owner
replaces its dense derivative assembly. Four parameter-column entries change;
the other 226 retained entries and the S9/S10/E5 comparisons remain within the
owner's 1e-12 relative requirement. The Khudaida numerical CSVs reproduce byte
for byte against `bdde4003`.

`results/columns.csv` and `centers.csv` retain the complete matrix comparison and
exact centers. `stencils.csv` retains the two ionic refinements, and
`native-comparisons.csv` the S9/S10/E5 values. `cost.csv` distinguishes the
initial interleaved proposal measurements from the final installed-call timings
measured in separate builds. `summary.json` identifies wheels, comparison
counts, file hashes, scope and downstream source inspection.

The one persistent falsifier is `boundary_checks packing-derivative`: eight
chemical/pressure rows at 300 K, 100 kPa, cation packing diameter active, with
steps 1e-4 and 1e-5 angstrom and an absolute 2e-9 per-angstrom criterion.

Reproduce the matrix and cost measurement by linking the driver against the native
build:

```bash
mkdir -p build/lever-row-series-comparison
bash -c '"$(tools/native-driver boundary_checks analyses/2026-lever-row-series-comparison/scripts/compare.cpp)" > build/lever-row-series-comparison/comparison.log'
```

The command writes a new log under `build/`. It calls the current shared Engine
calculation; the prospective duplicate assembly has been removed. The baseline
values are retained from the `bdde4003` build. Timing changes with shared-host
load and is not expected to reproduce byte for byte. These are numerical
verification and bounded cost measurements, with no physical validation or
parameter adoption claim.

# Fixed-depth native association comparison

This analysis retains one canonical installed-wheel campaign comparing the
fixed-depth log-Newton association candidate with the implicit mass-action
reference at depths 8, 16, 24, and 32. It covers the 14 named controls plus
64 development and 32 held-out generated states across eight strata, including
weak/strong association, disconnected/nonassociating components, seeded
topologies, active association energy and dimensionless kappa coordinates, and
trace/exact-zero support.

The seeded mixture set uses `std::mt19937`; each association energy/volume
matrix is generated with symmetric tied entries and the actual seed is retained
in every corresponding comparison input. Generated input
manifests retain the actual arrays, kappa/Delta distinction, base-state
diameters, packing fraction, topology provenance, and ordinary/stress/support/
boundary regime labels.

Prepare the canonical environment and rerun the retained campaign with:

```bash
bash tools/setup.sh
.venv/bin/python analyses/2026-association-native-comparison/scripts/run.py
```

The script records the native runner's homogeneous comparison table once at
[`results/association-comparison.csv.gz`](results/association-comparison.csv.gz)
(gzip; the recorded SHA-256 is of the decompressed CSV),
one keyed input table at [`results/inputs.csv`](results/inputs.csv), a
selector/status manifest at
[`results/campaign-manifest.json`](results/campaign-manifest.json), and compact
provenance at [`results/summary.json`](results/summary.json). The native
campaign is capped at 600 seconds and 4 GiB address space per process. Use
`--selector pilot` for the dense 24-site resource pilot or
`--selector case:<id>` to recover one affected state; missing state/depth rows
remain explicitly `unattempted` in the manifest.

The centered finite-difference set is frozen to pure 2B, pure 3B, strong,
heterogeneous-cross, active-energy, active-kappa, trace-support,
exact-zero-support, and one ordinary-dense generated state from each partition.
It checks both extensive energies at two steps, including mixed T/V and
state/active entries. State/composition steps use relative sizes `1e-4` and
`3e-5`; active-parameter steps use `1e-3` and `3e-4` to avoid cancellation in
the K/kappa coordinates. One-sided stencils are used when a nonnegative
amount or parameter boundary admits no negative perturbation; unresolved steps
remain recorded as unavailable/failed evidence.

This is a bounded comparison result, not qualification of every mixture,
pressure-closed density derivative, equilibrium specification, or production
solver replacement.

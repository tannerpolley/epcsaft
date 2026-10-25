# Historical fixed-pressure anchor analysis

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The original
producer was removed by #341 because its fixed-pressure anchor and density-diagnostic
interfaces were retired.

The retained study ran on Engine commit
`f55b217ceef9be78253885f90caf91a8f5f6c491`, wheel SHA-256
`687481a96890443023ef52388fd5b45fa5ff001d45bd5e4b9bb65cf16337658e`.
`results/summary.json` and `results/calls.csv` retain its inputs, comparisons,
numerical diagnostics and timings. They are historical numerical verification,
not qualification of current anchored pressure requests.

`scripts/run.py`, the method record, was removed by #341 (record commit `21d7a112`). It used the retired
`epcsaft.unit_registry`, quantity-valued state outputs, `density_anchor` keyword
and `State.density_diagnostics` view and did not run on the current wheel.
The checker and notebook read retained values without running the model.
Current reactive MEA anchored/continuation evidence belongs to
[`2026-mea-recovered-state-actions`](../2026-mea-recovered-state-actions/README.md)
(#131/#61), which does not requalify this older density-anchor study.

# Native diagnostic trace

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The original
capture-preservation and cost probes require a special fault-injection build and
generated headers excluded by #696.

Issue #284 qualifies optional capture of equilibrium polishing, density scans,
association iteration and derivative failures, anchors and solved Jacobians.
The revised design is supplied by the accountable lead. Numerical qualification
uses the installed non-editable Engine wheels from the same scientific base.

The pinned cold-scan input and 50-iteration reference originate in #275 at
`19a0c825`; they describe a cryogenic interruption, not physical validation.
The masking and MEA endpoint declarations (`scripts/trace_declarations.py`,
`scripts/endpoint_declarations.py`; under `inputs/` until #341) are frozen copies of
their sources at `9fdfe7ce` (`engine/tests/equilibrium/test_trace_reliability.py`
and `2026-addition-enthalpy-fitting`); #341 replaced the original source slicing
with these copies, field-for-field identical, so `results/provenance.json`
hashes name the originals. #272 owns the masking sequence interpretation. The notebook reports capture preservation and measured cost.

Raw timing rows (`results/timing*/pairs.csv`, and `results/timing/{starts,states,work}.csv`) are at record commit `21d7a112`; each timing folder keeps its `summary.json`, and `results/timing-r2-final/pairs.csv` stays.

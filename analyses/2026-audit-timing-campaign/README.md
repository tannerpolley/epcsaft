# Audit timing campaign (#578)

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The producers
require the frozen `results/image.json`, its guarded Engine image and an external
persistent Julia service; they are not portable current-main commands. The
external audit probe scripts they load (`build/audit-2026-10-15/`) were lost on
2026-10-28, so the producers can no longer run; the retained results stand.

What do the four 2026-10-15 audits' speed questions (#551–#554, merged as T1–T14)
cost at current main, against FeOs 0.10.1 and Clapeyron 0.6.29 at identical
parameters? The [notebook](notebook.qmd) answers each kept question from the
retained files in `results/` and states, per audit question, whether the finding
still holds, was fixed by a named issue, or is new. The accepted design is the
main checkout's `build/issue-578/design.md`; by the owner's 2026-10-16 decision
the timed commit is current main instead of `517a3450`.

Every Engine interpreter loads the frozen image in `results/image.json` through
`scripts/guard.py`, which forbids rebuilds and process launches; native fit
workers start through the same guard (`scripts/fit_launch.py`). Clapeyron runs
in one persistent Julia service for the whole campaign. Results are shared-host
software measurements and numerical verification, not physical validation.

The original campaign used the guarded preparation and measurement methods
named above; its executable rerun recipe is withdrawn.

Each section (`T7`, `T3`, `T1-verify`, `T1`, `T2`, `T9-serial`, `T9-fits`, `T6`,
`T4`, `T8-verify`, `T8`, `T10`, `T10-cohorts`, `T11`) writes its own
`results/raw/<section>.jsonl.gz` and refuses to overwrite one; the Clapeyron
sections need the `epcsaft-578-clapeyron` service running
`scripts/clapeyron_server.jl` on `build/issue-578/clapeyron.fifo`.

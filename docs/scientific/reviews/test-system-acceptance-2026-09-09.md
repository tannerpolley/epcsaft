# Test-system implementation acceptance

Date: 2026-09-09. Current full Fable 5.1 High disposition: Not ready;
interruption recovery, worker liveness and nested recording corrections are in progress.

Scope: the five preparation steps in the
[implementation plan](../implementation-plan.md#test-system-preparation), against
the [test-system design](../test-system-design.md). This qualifies recording,
applicability and display behavior only within the checks actually exercised.
The earlier restart check did not establish reboot persistence, and the integrity
harness was not individually registered in the notebook. The implementation plan
tracks their correction and a fresh full review. All 35 scientific sublayers remain Pending.

## Executed evidence

| Exercise | Observed result |
|---|---|
| Shared tooling harness, `python3 tools/check-test-system.py` | 17 passed, none failed or skipped. Includes real pytest collection/outcomes, criteria and incomplete evidence, record corruption, input/artifact invalidation, shared input groups and locking. Synthetic native protocol and worker transition cases are explicitly labelled simulated. |
| Installed gate, `bash tools/check-affected.sh HEAD` | The expanded gate passed in 66.44 seconds, including 41.90 seconds for recorded automatic evidence. The normal direct integrity command took 29.90 seconds and recorded all 17 individual passes. The warm-wheel target below 10 seconds was not met. |
| All seven installed native suites through `tools/run-numerical-evidence.py --native-suite` | 49 Passed, one Failed, no Error. Contributions 8, Properties 9, diagnostics 6, specification 9, boundary 7 and Properties floor 1 passed; equilibrium had 9 passes and the established F5 failure. |
| Numerical preservation | All 24 retained task-start F1/F5 operand pairs preserved. F1 Passed; F5 Failed with 990 complete, explicitly judged Hessian rows. Its Jacobian comparisons are also retained. All 36 task-start Engine runtime file digests were unchanged. |
| Real two-tab browser check | Both tabs Synced with complete inventory; all 17 integrity outcomes were visible, including four explicitly simulated checks. Contributions and Properties categories appeared, all details were initially collapsed, and there were no browser errors. |
| Real source edit and restoration | Adding a temporary computational input made retained F5 Outdated. One automatic selection completed with exit 0; a second worker and the idle watch did not duplicate completed cases. Removing the input restored the generation and completed the restored automatic selection. |
| Real managed-service restart | Completed run IDs remained unchanged after restarting the viewer and worker. |
| Shell/service artifact stability | After canonical preparation, managed-worker preparation and the installed gate, retained F5 was still Current. Locale-independent fingerprint sorting corrected the previously observed unnecessary wheel replacement. |

The generated inventory contains 136 collected declarations (132 Implemented,
4 Planned), with no collection errors; the public snapshot groups parameter cases
into 92 declarations. Parameter cases retain their identities;
the notebook groups their comparisons. Actual shared inputs are retained once
per group, with per-row overrides and explicit null values preserved. Case-definition
metadata is separate from scientific inputs. On the inspected live snapshot,
JSON was approximately 3.14 MB. After individual integrity registration, two concurrent requests completed in
1.287 and 2.492 seconds, both Synced with complete inventory. The final native
replay left F5 Current/Failed. These are observations, not performance guarantees.

Browser transition checks also exercised Updating, expiry, Offline, recovery,
search focus/cursor preservation, owner failure badges and blank stale result
cells. Those transition cases used controlled snapshots; the table above identifies
the separate checks against the actual running service.

## Review and limits

The current correction has a persistent enabled systemd user service with
lingering and an explicit startup PATH; Docker is excluded. A computer reboot
was not performed. The final restart check preserved all 86 completed run IDs, the successful
terminal worker record and F5 Current/Failed. A regression check also preserves
failed terminal details without rerunning unchanged cases. The stale Queued
status found during acceptance was corrected. The requested Fable 5.1 High
full review then found additional interruption-recovery, dead-worker visibility
and nested-record reconciliation defects; corrections and closure are pending. Same-process registry
add, rename, remove and same-size/restored-mtime edits were detected correctly;
missing harness source was rejected.

The requested independent Fable Medium review initially returned Not ready.
Corrections addressed migrated criterion mismatches and omitted comparisons,
shared input retention, public failure reasons, live-response cost and artifact
preparation stability. Closure review returned Ready. The additional dependency
version/search-path cache and corrupt-record applicability fixes passed the
17-check harness, installed gate and two-tab check; their focused independent
review also returned Ready. All seven native suites were replayed again against
the final computational inputs with the same 49 Passed / F5 Failed outcomes.

Nonblocking review notes: schema-style coordinate labels in some rows coexist
with actual coordinate values in inputs; the cache regression exercises invalidation
but does not assert a hit for an unchanged context. Runtime-probe timeouts can exceed
the browser request timeout, in which case the browser fails closed; final concurrent
live requests completed within that timeout. None of these observations qualifies science.

No tolerance, reference value or thermodynamic algorithm was changed to make a
test pass. Existing F5 disagreement remains a scientific limitation. Planned
FEOS/Clapeyron cases were neither installed nor executed. Slow analysis replays,
full clang-tidy and sanitizers were outside this infrastructure validation.

Run records and numerical tables remain rebuildable under
`build/numerical-evidence/`; this document retains the bounded acceptance conclusions.
The implementation plan remains the sole delivery-status owner.

The separate [Core readiness review](core-readiness-2026-09-09.md) retains C1–C4
gaps without changing scientific qualification. The performance target remains
unmet; reducing recorder work must preserve detection of external edits during a run.

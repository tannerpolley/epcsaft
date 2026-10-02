# Historical test-system acceptance observations

Date: 2026-09-09. Scope: bounded recording, applicability and display checks
of the retired integrity harness and worker system. These are historical
observations, not current execution instructions or scientific qualification.
The [test-system design](../test-system-design.md) owns current viewer behavior.

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

## Review and limits

A computer reboot was not performed. The final restart check preserved all
86 completed run IDs, the terminal worker record and F5 Current/Failed.
All seven native suites were replayed against the final computational inputs
with the same 49 Passed / F5 Failed outcomes. These checks do not qualify science.

No tolerance, reference value or thermodynamic algorithm was changed to make a
test pass. Existing F5 disagreement remains a scientific limitation. Planned
FEOS/Clapeyron cases were neither installed nor executed. Slow analysis replays,
full clang-tidy and sanitizers were outside this infrastructure validation.

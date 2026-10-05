# Declared equilibrium starts: delivered review

**Supported for #238 starts-only; no remaining required corrections.** Independent read-only review inspected runtime `daa98f1f` and final wheel SHA-256 `4e00f1536cddea911baa5cf157892a8fc6def636aceeb6aa9b5d74ff1fa1540d`. Evidence/prose/deduplication is retained at `4de9fbe5`; runtime is unchanged. The reviewer executed no numerical probes.

The manufactured exchange at 200 K and 5 MPa closes original rows to $5.773159728050814\times10^{-15}$ against $10^{-10}$. Phase amount error is at most $4.440892098500626\times10^{-16}$ mol; composition error is at most $1.1102230246251565\times10^{-16}$. A one-optimizer-iteration start retains original residual 1.224413540303225, no selected state and no later attempted start. The shared polish, supplied partial/pinned seed admission and measured-work corrections are resolved.

Elapsed exhaustion governs state availability separately from numerical classification. The request to force a numerical budget label after deadline expiry is **withdrawn**: it conflicts with the accepted design and ALG-TERMINATION. Measured rows retain the common nonfinite/converged/floor/budget priority. Success is false, no state is published, and controller cause/overrun records elapsed rejection. Unmeasured NaNs remain nonfinite.

The main-integrated wheel logs confirm 114 Python tests, six native targets, ruff, mypy and architecture passed. The full gate initially failed only on duplicate main/branch state tables; the identical rows are now retained once. Final gate `d7e3279e` completed with all 114 Python tests, six native targets, ruff, mypy, analysis layout and architecture passing. The affected check reported no changed paths. The reviewer confirmed Supported again after reading that log, with no required corrections.

Two frozen default solves preserve iterations/evaluations/Hessians (16/49/17 ionic LLE; 43/57/44 reacting MEA), objectives and residuals, with 1,176 binary64-identical state rows against the retained main baseline. Wall differences above 10% are flagged shared-host variation, not algorithmic speedup. Evaluations include exact polish. [Retained operands and provenance](../../../analyses/2026-tpd-phase-discovery/results/issue238-cost/README.md) own the exact comparisons.

Growth of 374 added / 50 removed / +324 net runtime lines is justified by admission, between-kernel deadline checks and measured failure accounting under #248's paused ceilings. The final integrated total is 22,013; no architecture ceiling is amended.

Born/Uyan evidence remains scoped to its recorded `624881f1` wheel. Uyan historical agreement fails and reference reconstruction remains with #234; its old generator/allowlist is restored. Finite Born start failures establish neither physical nonexistence nor global stability. This review qualifies the declared-start capability on the stated cases, not Uyan chemistry, application acceptance or global root completeness.

Review provenance: Pi `openai-codex/gpt-6.1-sol`, xhigh, existing independent reviewer task `issue238-pi-review-1`, final terminal result 2026-10-02 (reviewer turn 6).

# Independent delivered review: #240

## Verdict and inspected snapshot

**Supported — scientific parameter-continuation capability.** Reviewed candidate
`853574a5` against integrated main `875526d0`. The global-scalar selection and
literal Born-column tolerance findings are closed; no further required
scientific correction was identified. Performance qualification is separate.

The lead checked the findings against source, the executed qualification logs
and retained values before accepting this bounded verdict.

## Scientific ownership and integration

The existing continuation controller, exact row/sensitivity actions, shared
scalar/vector parameter access, native fitting chart and trial-model lifetimes
remain the numerical owners. No second solver, Python numerical continuation or
changed thermodynamic formulation was introduced.

Upstream density changes extract the existing slope search and reject nonfinite
derivatives. Diagnostic negative-slope roots do not enter ordinary stable-branch
selection. Post-merge continuation evidence exercises density admission. Stop
naming retains its native meanings. Budget and allowlist conflict resolutions
preserve the upstream delivery and remove the converted continuation entry.

## Independent checks actually executed

On the current installed wheel, the reviewer executed seven global-chart
comparisons and six component-target rejection probes; 24 strict native Born
comparisons; 83 manufactured-exchange comparisons; and 4,500 historical-point
comparisons in memory. All numerical criteria passed. The reviewer checked the
six paths' termination causes, native work and failed residuals, accepted raw and
classified residuals, unchanged first-action ladders, package payload hashes,
per-layer accounting, timing cohorts, source diffs and whitespace. Analysis
layout passed with 82 leaves.

The supplied post-merge logs establish successful full/affected gates: 115
Python checks, six native programs, and affected Equilibrium, Boundary and
Parameter-regression-fit programs. Native/test sources are unchanged between
the gated merge commit and the reviewed candidate.

No reviewer writes, builds, installations, timing experiments or expensive
six-path/full-gate reruns were performed. Numerical values, criteria, physical
inputs and wheel identity are retained in the neighboring CSV/JSON records and
[qualification](summary.md).

## Growth and duplication

Independently verified 103 gross added, 60 removed and +43 net runtime lines:
21,689 to 21,732. All per-layer counts match the code-budget closeout. No numeric
ceiling changed. Growth is below the 270-line task target and 24,000-line
envelope; the corrections preserve shared storage, parameter actions and one
native continuation controller rather than duplicating numerical methods.

## Performance qualification

Current-main/current-branch inputs, returned states and iteration/evaluation/
Hessian work are identical. Best wall times from three matched pairs are:

- Ionic LLE: main 0.056040324 s, branch 0.056819409 s; 16/16 iterations.
- Reacting MEA: main 0.099776606 s, branch 0.100137661 s; 43/43 iterations.

These are increases of 1.390% and 0.362%. **Strict no-wall-time-increase is not
supported.** Three shared-host samples cannot establish causation, a speedup or
universal performance behavior. All earlier cohorts remain retained. This
caveat does not negate the bounded scientific-capability verdict; the issue's
numerical acceptance specification contains no strict timing threshold.

## Claim limits

All six Born paths terminate `infeasible`; neither c=0 nor epsilon=2 was reached.
Supported evidence covers finite accepted prefixes, not folds, endpoint
nonexistence, global stability or physical validation. Propagated residual-floor
estimates are not rigorous derivative-error bounds.

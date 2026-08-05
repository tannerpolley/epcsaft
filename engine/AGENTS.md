# Unified Engine Repository Agent Contract

Repository Profile: scientific-computing

Read doctrine revision 10 in `../ePC-SAFT-governance/GOVERNANCE.md` and the
machine-readable topology in
`../ePC-SAFT-governance/manifests/management.json` before changing this
repository. Local policy may narrow that doctrine but must not contradict it.

This repository owns the single `epcsaft` runtime distribution and its three
public module surfaces:

```text
epcsaft
epcsaft.equilibrium
epcsaft.regression
```

Its end-to-end scientific lifecycle is:

```text
materialized Data packet -> Parameters -> RegressionResult.to_parameters()
                         -> Mixture/State -> Equilibrium result
```

Keep one root Python build, one CMake build, one version, one wheel, and one
installed module origin. Do not add compatibility imports, alternate package
identities, installed-provider discovery, sibling source discovery, mutable
wheel paths, hidden physical defaults, or duplicate equation owners.

Data remains the owner of reusable scientific packets. Engine tests may retain
exact byte-for-byte packet fixtures only when a materialization receipt binds
the Data repository, commit, packet path, version, fingerprint, and every file
hash. Runtime code and tests must never discover a sibling Data checkout.

Validation remains the black-box evidence repository. Installed-artifact
campaigns must consume one immutable Engine wheel and hash; they must not import
this source checkout or private implementation seams.

Preserve historical plans and receipts as provenance. Their former repository,
package, symbol, and artifact names are historical facts, not live setup
instructions.

Use Numerical Verification-Driven Development from Governance section 7.
Retain tests only for numerical values, identities, derivatives, limits,
solver or equilibrium certificates, manufactured solutions, metamorphic
relations, and source-backed scientific benchmarks. Do not add or retain API
compatibility, legacy, mock, exact-prose, snapshot, solver-label, iteration
count, or coverage-only tests. A private seam earns a test only when it uniquely
exposes a numerical invariant such as a production derivative tensor.

Classify Engine verification into three lanes. The default `Test` environment
action owns compact persistent numerical evidence and excludes the `campaign`
marker. Expensive replays use the explicit campaign route and normally belong
in Validation. The `Check`, native-check, sanitizer, build, typing, lint,
environment, and legacy-path commands are readiness checks rather than
scientific tests.

For each retained parameterized matrix, use one module-level evidence header
or adjacent manifest to satisfy Governance section 7's claim, oracle,
quantity, units or scaling, tolerance-basis, and failure-meaning record. An
expected value produced by this implementation is change detection only.
Copied production equations, generated values, or shared kernels do not count
as independent verification and must be declared when retained.

This private repository has no CI/CD. Keep reproducible checks in
`.codex/environments/environment.toml` and repository scripts.

Current issue tracking and runtime authority belong to
`ePC-SAFT/ePC-SAFT-engine`. Resolve the currently admitted immutable wheel and
full SHA-256 from Governance's artifact-store records; do not encode a mutable
`dist/` path or a release version in this instruction file.

After changes, run the applicable named actions in
`.codex/environments/environment.toml`, including the installed-wheel Engine
suite and readiness checks affected by the change, then run the repository
cleanup audit.

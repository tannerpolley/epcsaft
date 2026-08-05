# Unified Engine Repository Agent Contract

Repository Profile: scientific-computing

Read doctrine revision 9 in `../ePC-SAFT-governance/GOVERNANCE.md` and the
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

Current issue tracking and runtime authority belong to
`ePC-SAFT/ePC-SAFT-engine`. The accepted `v0.2.0.dev0` Release wheel and its
Governance artifact-store object are the only supported installed runtime.

After changes, run focused scientific checks, the retained Engine suite from an
installed wheel, static legacy-path searches, and the repository cleanup audit.

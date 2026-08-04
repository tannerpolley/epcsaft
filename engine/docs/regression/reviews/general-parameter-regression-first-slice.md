# General Parameter Regression First-Slice Evidence

Status: `FIT_READY_AUTHORITY_NEUTRAL`

Date: 2026-07-27

## Implemented contract

The first general pair families independently fit one shared, dimensionless
`k_ij` or `l_ij` for a
caller-supplied neutral, nonassociating binary model and source-bound
fixed-composition VLE rows. Each training row owns two lifted log-volume
coordinates and four scaled residuals:

```text
(P_liquid - P_observed) / pressure_scale
(P_vapor - P_observed) / pressure_scale
(mu_1,liquid/RT - mu_1,vapor/RT) / chemical_potential_scale_1
(mu_2,liquid/RT - mu_2,vapor/RT) / chemical_potential_scale_2
```

The exact `4R x (1 + 2R)` Jacobian consumes the installed Provider's
row-major value/gradient/Hessian contract in
`(n1,n2,V,active_pair_parameter)`. No third
derivatives, density roots, copied EOS equations, numerical production
derivatives, or Equilibrium dependency are present.

## Exact subjects

- Provider commit:
  `1e571ab0a84603a51ed6994b14286f683fb12b88`
- Provider tree: `910e01c1571d4d6128db60a30d6e1948f5c8ac4c`
- Provider wheel:
  `~/Workspaces/Engineering/ePC-SAFT-project/artifacts/provider-general-parameter-sdk-v1/1e571ab/epcsaft-0.1.0.dev0-cp313-cp313-linux_x86_64.whl`
- Provider wheel SHA-256:
  `6536edc63adaa13c5c6c67c185d82c9ae232048e99dc3dc3be502708eea4410f`
- Installed Provider header SHA-256:
  `b667379c2d7106d012c6b57f96b6f32dd23ef305fe6f15a87c22ab20029008f8`
- Regression commit:
  `da7a44ce093201022aec2f3514d4e4fd9d8d2929`
- Regression package tree: `2f4f82fbb5751317b8314309ff4246f0ebfde7ee`
- Regression wheel:
  `~/Workspaces/Engineering/ePC-SAFT-project/artifacts/regression-general-parameter-v1/da7a44c/epcsaft_regression-0.1.0.dev0-cp313-cp313-linux_x86_64.whl`
- Regression wheel SHA-256:
  `3a59d2233fec51f949a7784937b54f7f66beae2476fa8e33976672a480b67137`

## Verification

The exact installed wheels passed:

```text
90 full installed-artifact tests passed; 3 campaigns deselected in 228.92 s
2 explicit all-17-row pair-family campaigns passed in 0.36 s
115 Provider native-SDK/association/EOS tests passed in 309.47 s
56 Validation tests passed in 0.88 s
```

The explicit campaign consumes
`evidence/may-2015-methane-ethane-vle.csv`, SHA-256
`5cd1e74925a3c6504f5106dcf911f2cae2d6e99a5133fccc20454d8991bdbc7f`.
It reports:

```text
rows                         17
residuals                    68
variables                    35
Ceres termination            CONVERGENCE
full lifted rank             35
projected parameter rank     1
fitted k_ij                  -0.00843032298906253
fitted l_ij                  -0.002774426668544412
active parameter bound       none
declared confirmation starts 2
```

Independent scientific and code reviews caught and closed false convergence
reporting, missing Provider-domain preflight, unbounded fixed-field reads,
post-Ceres diagnostic loss, incomplete result provenance, Provider capability
superset handling, per-iteration Regression allocations, legacy native-model
loading, and accepted pure-component numerical parity.

## Interpretation

This evidence establishes reusable, source-bound, exact-derivative
fixed-composition neutral-binary `k_ij` and `l_ij` runtimes. It does not
establish predictive validity, global uniqueness, uncertainty, catalog
authority, or support for association, `k_hb_ij`, polar, dielectric, reactive,
or other parameter families. The May rows are all training data. The
historical pressure-closure miss remains immutable negative evidence and is
not relaxed or reinterpreted here.

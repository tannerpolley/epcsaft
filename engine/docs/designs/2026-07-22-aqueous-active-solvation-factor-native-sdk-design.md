# Aqueous Active Solvation-Factor Native SDK Design

Status: implementation candidate for Project Truss Provider issue #22. This
document owns only the appended Provider callback; Regression's staged target,
residual, bounds, starts, ranks, and acceptance remain in
`ePC-SAFT-regression/docs/science/figiel-aqueous-kij-regression.md`.

## Capability

For the exact ordered Figiel `(water, sodium-cation, bromide-anion)` model at
`298.15 K`, `100000 Pa`, and formula-unit molality `0.001-6 mol/kg`, the
callback accepts one finite positive trial water solvation factor and returns:

- `ln(gamma_pm^m)` on the existing molality convention;
- exact total fixed-pressure
  `d ln(gamma_pm^m) / d f_water`;
- terminal reference molality and value/derivative convergence diagnostics;
- the bound parameter fingerprint; and
- structured status/error output.

The ABI appends a 280-byte size-gated result and two 64-bit table fields at
offsets 208 and 216. The complete table becomes 224 bytes; every prior field
and offset is unchanged.

## Derivative owner

`f_water` is one CppAD independent coordinate in the existing scalar-generic
Born SSM+DS contribution. The existing association implicit solve, pressure
root, fugacity construction, and four-molality electrolyte reference sequence
remain the only owners. At each finite or reference composition, Provider
forms the exact fixed-pressure derivative by

```text
rho_f = -P_f / P_rho
d ln(phi_i)/df |_P = partial_f ln(phi_i)
                       + partial_rho ln(phi_i) * rho_f
```

and returns

```text
d ln(gamma_pm^m)/df_water
  = 0.5 * [(d ln(phi_c^finite)/df + d ln(phi_a^finite)/df)
           - (d ln(phi_c^reference)/df + d ln(phi_a^reference)/df)].
```

The molality conversion is independent of `f_water`. No third derivative or
numerical production derivative is required.

## Trial-bundle admission

Installed-catalog NaBr remains fingerprint-bound. A user-provided staged trial
is admitted only when resolution proves corrected Figiel SSM+DS topology,
ordered water/+1/-1 components from Li+/Na+/K+ and Cl-/Br-, exact
`figiel-yu-held-2025` record/domain provenance (plus the exact NIST water molar
mass identity), and the frozen `298.15 K`, `100000 Pa`, ion-fraction domain.
Within that boundary, the existing active-Born callback is exposed only for its
five declared tracer selections, the active-`k_ij` callback for the six aqueous
salts, and this callback only for NaBr. Reordered, incomplete, unknown-domain,
wrong-purpose, wrong-source, or arbitrary models fail closed. Generic phase
transport for user bundles remains independently available and unchanged.

## Evidence gate

At `f_water=1.35` and `1.65`, all 21 audited NaBr molalities are checked against
centered callback values at steps `1e-4` and `5e-5` using

```text
abs(J_exact - J_h/2)
  <= max(1e-8, 20*abs(J_h - J_h/2), 2e-8*abs(J_exact)).
```

The catalog coordinate `1.5` must reproduce the existing public reference
value and diagnostics. Header layout, prefix offsets, malformed result size,
fingerprint, state, nonfinite input, trial-bundle negative space, complete
Provider suite, installed wheel, linkage, and cleanup must pass.

Excluded: regression targets/residuals/Ceres, bounds/scales/fit decisions,
catalog mutation, generic overlays, arbitrary solvents or ions, copied EOS or
reference equations, density alternatives, finite-difference runtime, second
or third parameter derivatives, prediction, promotion, and authority change.

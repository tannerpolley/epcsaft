# Explicit modular baseline on the current Engine

## Question and input

Can the declared Figiel model and explicit-charge neutral packets reproduce
their retained fixed-density calculations after admission is made fail-closed?
The canonical equations and physical meanings remain in
[full-eos.md](../../../engine/docs/science/full-eos.md). This review uses the
single installed noneditable wheel identified in the
[manifest](modular-baseline-admission-2026-09-22.json), never Engine source imports.

The Figiel input is the immutable figiel-yu-held-2025/3 parameter document
(SHA-256 3d37235cd8217e056f55005bfdf2b65dcf7c4a8621e239bdd40a70c273a9be8c),
selected in order as water, sodium-cation, chloride-anion. The state is
298.15 K, 56000 mol m^-3, and mole fractions (0.98, 0.01, 0.01).
The 23 retained donor scalars include all five Helmholtz contributions,
pressure, residual calorics, the pressure density derivative, two site
fractions, fugacity coefficients, and residual chemical potentials.

The neutral control selects each of 82 components separately from the
version-3 Gross/DEEA packets at 300 K and 10 mol m^-3. Five scalars per
component cover pressure, its density derivative, residual Helmholtz energy,
residual internal energy, and residual isochoric heat capacity. These
version-3 packets add explicit zero charge without changing physical
coefficients; historical version-2 packets and results remain intact.

## Numerical result

[The CSV](modular-baseline-admission-2026-09-22.csv) retains one row per
quantity, with the donor value, current-wheel value, absolute difference,
and error normalized by max(1e-12, 1e-12 max(|current|, |donor|)).
The donor control is commit
45fbef5ca3f4f9ff7a72c31f17b49400034f0879; the current packet and wheel
hashes are in the manifest.

| Control | Rows | Largest normalized error | Result |
|---|---:|---:|---|
| Figiel fixed density | 23 | 0.024670 | all within criterion |
| Neutral charge | 410 | 0.0008884 | all within criterion |

An independent symmetric water 2B reduction gives bond density
536.45118457612296 and solves B X² + X - 1 = 0. The current installed
general-site model reproduces both site fractions and association
Helmholtz energy to relative 1e-11 and absolute 1e-12 in the Engine
numerical check. Changing document/family display IDs leaves the sampled
pressure, association, residual energy, pressure derivative, site
fractions and chemical potentials invariant to relative/absolute 1e-12.

Duplicate water records cannot select between physically different
coefficients silently. In a single-record sensitivity calculation, reducing
water dispersion energy by 10% changes the sampled pressure from
17.1188294 MPa to 94.3738520 MPa. The duplicate-record input is rejected
before subset selection; association none with supplied sites, missing
neutral charge, and malformed component/endpoint identities also fail by
named diagnostic.

## Interpretation and limits

These are equation and admission checks, not source-domain validation.
The fixed-density Figiel pressure is far above the packet's reported
1-bar conditions. Its source evidence is at 298.15 K and reaches ion mole
fraction about 0.38; no wider applicability or fit transfer follows.
The neutral rows test numerical preservation at their listed dilute
states, not source uncertainty. The Figiel packet retains sources,
units and transformations, but a complete fixed-versus-fitted partition
and calibration-observation attachment remain open in
[issue #25](https://github.com/tannerpolley/ePC-SAFT/issues/25).
Previous IDAES and pure-saturation results stay bound to their recorded
version-2 inputs and older wheel. Their current runners have version-3
inputs; those campaigns were not replayed here.

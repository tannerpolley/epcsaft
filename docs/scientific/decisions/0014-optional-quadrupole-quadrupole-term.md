---
status: accepted
date: 2026-10-16
---

# Admit an optional Gross quadrupole–quadrupole residual term

Kind: scientific model admission (#574). This partially supersedes ADR 0006 for
the quadrupole–quadrupole term only and extends ADR 0011.

## Context

ADR 0006 retired the residual point-multipole family; ADR 0011 readmitted the
Gross–Vrabec dipole–dipole term. CO2, N2 and aromatics lack the standard PCP-SAFT
quadrupole treatment. Gross (2005, AIChE J. 51:2556) gives it with tabulated,
independently measured moments and no extra adjustable parameter.

## Decision

- **Opt-in choice.** The `polar` kind gains the choice `quadrupole-quadrupole`.
  `none` stays the default.
- **Inputs.** Every selected component supplies `quadrupole_moment` in
  `debye * angstrom` and `quadrupole_count` (n_Q, dimensionless), finite and
  non-negative; zero is allowed. They are rejected as unconsumed under any other
  choice, and the dipole records are rejected under this one.
- **Equations.** Gross (2005) Eqs. 8–10 and 12–19 with Table 2, as in
  [equations](../../../engine/docs/equations.md#quadrupolar_helmholtz). The
  two-body third-order integral is zero (Eq. 14). ε_ij = √(ε_iε_j) without k_ij;
  m_ij = √(m_im_j) and m_ijk = (m_im_jm_k)^{1/3} with **no chain-length cap**. The
  source prints none; the m ≤ 2 cap is the Gross–Vrabec (2006) dipole convention
  and stays dipole-only. The moment bridge is the dipole's debye constant, which
  uses the CODATA 2018 vacuum permittivity.
- **Where the equation lives.** Contributions evaluates both terms in one
  multipole function that keeps the dipole's operation order, sharing η, σ_i, ε_i,
  the Padé form and the J-polynomial helper on the same tape. The helper no longer
  caps; the dipole caps its own chain arguments. The result member `dipolar` carries
  whichever polar term the model declares.
- **Unchanged models.** An undeclared family evaluates no multipole code, and the
  dipole values are bit-identical to those before this decision.
- **Comparisons.** Polar models are compared only with their own parameters
  (ADR 0011).
- **Still retired.** The Vrabec–Gross dipole–quadrupole cross term. Its sources
  disagree on the third-order term, and `polar` takes one choice, so no model
  declares both dipolar and quadrupolar components.

## Validity and verification

The constants were fitted to two-centre Lennard-Jones plus point-quadrupole
simulations over elongations 0 ≤ L* ≤ 0.8 (m up to 1.8447) and
0.5 ≤ T/T_c ≤ 0.99, with (Q*²)² ≤ 4 in the 2CLJ reduction, for linear symmetric
molecules. Larger m (benzene, m = 2.2463) is the source's own extrapolation. The
source finds its mixing rule unsuitable for real mixtures with more than one
quadrupolar component (benzene + CO2, Figure 20 and Conclusions); the CO2/N2
bubble-point check is a code-identity test of that rule, not validation.

No claim is made for ions or electrolyte combinations, or for composition
derivatives above first order along an absent quadrupolar species (the dipole's
A2-only limit).

Numerical evidence lives in
[analyses/2026-polar-quadrupole](../../../analyses/2026-polar-quadrupole/README.md)
and `engine/tests/eos/test_quadrupolar.py`.

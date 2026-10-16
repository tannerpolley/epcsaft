---
status: accepted
date: 2026-09-23
---

# Admit an optional Gross–Vrabec dipole–dipole residual term

Kind: scientific model admission (#28). This partially supersedes ADR 0006 for
the dipole–dipole term only.

## Context

ADR 0006 retired the residual point-multipole family. Issue #28 asks whether one
targeted new physical description improves a predeclared deficiency. The owner
authorized one term on 2026-09-23: the Gross and Vrabec (2006) dipole–dipole
contribution, known as PCP-SAFT (the polar part of ePPC-SAFT).

Nonpolar PC-SAFT parameters and their k_ij were fitted to absorb polar effects.
Adding a polar term on top of them therefore double-counts those effects, and a
polar model is only meaningful together with its own parameters.

## Decision

- **Opt-in family.** `model_families` gains an optional `polar` kind. Its choices
  are `none`, which is also the default when the kind is absent, and
  `dipole-dipole`; [ADR 0014](0014-optional-quadrupole-quadrupole-term.md) adds
  `quadrupole-quadrupole`.
- **Inputs.** `dipole-dipole` requires two explicit records for every selected
  component: `dipole_moment` in `debye` and `dipole_count` (n_μ, dimensionless).
  Both must be finite and non-negative; zero is allowed for nonpolar species.
  With `none`, these records are rejected as unconsumed.
- **Inputs are kept separate from permittivity.** The residual dipole inputs are
  distinct from the #27 relative-permittivity records `dipole_scaling`,
  `polarizability_scaling` and `correlation_integral_parameter`. Neither set
  feeds the other equation, and they are never translated into each other.
- **Where the equation lives.** The Contributions layer owns it, in
  `multipole_helmholtz`, shared with the quadrupole term since ADR 0014, using the equations cited in [equations](../../../engine/docs/equations.md#dipolar_helmholtz).
  It shares η, σ_i and ε_i with dispersion on the same tape.
- **Combining rules.** ε_ij = √(ε_iε_j) without k_ij (source p. 1196). k_ij
  remains a dispersion correction only.
- **Nonpolar models are unchanged.** A model that does not declare the family
  evaluates no polar code and is bit-identical to the pre-#28 Engine.
- **Comparisons.** A polar model is compared only with its own published or
  documented-refit parameters, against a nonpolar model with its own parameters
  on the same data. It is never added to nonpolar-fitted parameters.
- **Still retired.** The Vrabec–Gross dipole–quadrupole cross term remains retired
  under ADR 0006 until a separate owner approval; the quadrupole–quadrupole term is
  admitted by ADR 0014.

## Validity and verification

The J-integral constants were fitted to two-centre Lennard-Jones plus
point-dipole simulations. These cover elongations 0 ≤ L* ≤ 1 (hence the caps
m_ij, m_ijk ≤ 2) and 0.55 ≤ T/T_c ≤ 0.99. The dipole moment is aligned axially,
and n_μ = 1 for small single-dipole molecules.

No claim is made for ions, for combining the term with an electrolyte family, or
for composition derivatives above first order along the direction of an absent
dipolar species.

Numerical evidence and the model/parameter comparison live in
[analyses/2026-polar-dipole](../../../analyses/2026-polar-dipole/README.md) and
`engine/tests/eos/test_dipolar.py`.

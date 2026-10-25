---
status: accepted
date: 2026-10-16
---

# Pure stoichiometric solids in the one equilibrium driver

Kind: computational architecture.

The owner [adopted option (a), this record and the build](https://github.com/tannerpolley/ePC-SAFT/issues/573)
on 2026-10-16 after an independent readiness review of the #573 design.

## Setting

- ADR 0009 and `ARCHITECTURE.yaml` excluded solids, so the Engine could not
  express salt solubility or precipitation.
- Lithium extraction from concentrated brines needs saturation indices and
  precipitated amounts. The salting-out LLE and solid saturation share one
  aqueous phase.

## Alternatives

- (a) A native pure-solid phase in the existing driver.
- (b) Couple PHREEQC (Pitzer) or Reaktoro, with the Engine supplying activities.

## Choice: (a)

- A pure stoichiometric solid has one coordinate, its amount n_s ≥ 0.
- Its standard Gibbs energy is g_s = Σ ν_s,i g_i + RT ln K_s^ρ(T, P°) + V_s (P − P°)
  (EqID `solid_standard_gibbs`). The existing reaction-reference transfer converts
  the sourced dissolution constant to the Engine density basis at the declared
  reference pressure P°, which is required; unset is refused.
- The objective gains n_s g_s/RT, the material balance gains Bν_s n_s, and the solid
  has no chemical row in the Ipopt problem. The amount bound then gives the
  complementarity n_s ≥ 0, SI_s ≤ 0, n_s SI_s = 0.
- One classification at the Ipopt endpoint makes the declaration square for polish
  and actions: a present solid gets the row g_s/RT = λ·Bν_s; an absent solid gets
  n_s = 0. Violations after the polish flip once; a remaining violation is a typed
  failure.
- Lead decision, 2026-10-16 (delivered review r1): a solid is present when
  n_s/N_feed exceeds its bound multiplier z_s = −A_s/RT, the complementarity
  indicator. The first rule, the scaled affinity against the row criterion,
  started a stable solid absent at small excess feeds, and with a second free solid
  the single flip round then failed.
- A `Pinned(0)` solid reports its saturation index without precipitation.

[ALG-EQUILIBRIUM](../../../engine/docs/science/algorithms.md#alg-equilibrium) owns the
method and [formulation](../../../engine/docs/equilibrium/formulation.md) the physical
meaning.

## Evidence

- One activity model per calculation and one driver.
- Exact derivatives: the solid term is linear in n_s; ∂g_s/∂P = V_s at fixed T.
- PHREEQC takes no external activities. Reaktoro would be a second solver
  without exact curvature.
- Halite in water at 298.15 K with the Held 2014 strategy-2 NaCl set
  (`engine/tests/equilibrium/test_pure_solid.py`): solubility 6.2989 mol/kg against
  Clarke–Glew 1985's 6.1476 (+2.46 %, tolerance 3.0 %); saturation index within
  0.032 of PHREEQC `pitzer.dat` at 5–6 mol/kg.

## Validity boundary

- Pure solids at fixed T and P.
- One pure neutral solvent with a charge-neutral dilution path.
- The Newton attempts are excluded.
- Discovery, stability trials, continuation, fit observations, ASL records, PH/PS,
  critical, envelope and azeotrope conditions are refused with solids.
- Solid solutions and kinetics stay outside the Engine.
- Each mineral and parameter set needs its own qualification, which belongs to
  the issue that adopts it.

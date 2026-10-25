# Current Held-2014 replay (#750)

Decision **(a)**: the changed printed extrema follow the intended infinite-dilution correction in commit `488451a320735a4882375c367647365960d61166` ([PR #619](https://github.com/tannerpolley/ePC-SAFT/pull/619)). The same MgCl2 state at 2.25 mol/kg has gamma 1.5623036395978573 on published wheel `517a3450` and 1.5622972483698117 on `cc63a135` and current main. [The bounded comparison](wheel-comparison.json) retains the values and wheel hashes.

The maximum printed gamma errors change from 0.085606 to 0.085605 decade for MgCl2 and from 0.089353 to 0.089351 for CaCl2, both at the original 2.25 mol/kg extrema. Only those diagnosed stored strings changed. The osmotic extrema, concentration checks, energy-coordinate invariance, Debye checks, positive osmotic coefficients, and every-row 0.02-decade criterion retain their original limits. LiCl and NaCl pass, and both divalent salts still fail the physical-data criterion. The [rerun log](held-rerun.log), [current extrema](sigma-replay/maximum-errors.csv), and [identity](sigma-replay/run-identity.json) retain the numerical result. The original sigma replay remains unchanged.

This is numerical verification and a model comparison against the retained evaluated aqueous data at 298.15 K and 100000 Pa. It supplies no new experimental validation, refit, or parameter adoption.

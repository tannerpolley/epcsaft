# Perdomo Figure 3 and Table 5

This directory preserves historical source rows, fitted values, replay
outputs, and results for an old split-wheel charged-parameter workflow. Its
specialized Regression executable was removed by the monorepo hard cutover.
The current Engine supports electrolyte parameter regression and phase observations. This
historical executable is not replayed: its inputs lack a confirmed Li transfer
concentration-scale convention and its two one-row fits have no held-out
evidence. The current conditional diagnostic retains the coefficients without
refitting; see the [live notebook](../../notebook.qmd).

The regression targets combine Figiel's retained water x-process single-ion
solvation values with published water-to-1-butanol transfer Gibbs energies.
The historical transform used a declared c-to-x conversion. The Li+ source
does not state its concentration-scale convention explicitly, so its c-process
classification remains a working assumption rather than a source fact.

The fitted values are:

- `k_ij(1-butanol, Li+) = -1.1524901988171197`
- `k_ij(1-butanol, Cl-) = 0.8430300746192005`

The two one-row fits had no held-out observations and confer no practical
identifiability, predictive, or authority claim. The results preserve their
exact historical split-wheel identities and differing EOS payloads; those
identities are provenance, not current installation instructions.

This material does not admit the candidate parameter bundle. The authoritative
Figure 1(a) result remains a separate NaCl/water namespace.

## Historical results

All three Table 5 endpoint-midpoint feeds return accepted two-phase states.
Their mean endpoint absolute errors on the explicit-species mole-fraction
basis are 0.03747, 0.02982, and 0.01114 for 4.58, 4.95, and 5.74 mol/kg,
respectively.

Of the 33 Figure 3 author-model tie-line midpoints, nine return accepted
two-phase states. Eighteen return a one-phase result with solver, numerical,
and physical checks passed; the final six are explicit flash-error
non-admissions. Every result retains `globality_certificate=not_guaranteed`.
Because the paper curve uses
SAFT-gamma-Mie GC and this campaign uses ePC-SAFT, this is a descriptive
cross-EOS challenge, not a same-EOS or predictive reproduction.

The synthetic midpoint feeds are only basin requests constructed from source
endpoints; they are not experimental measurements. The source curve itself is
also the authors' model output rather than experiment.

## Retained layout

- `source/`: source tie lines, metadata, and historical fit targets.
- `../../results/data/figure-03/`: historical row outputs.
- `../../results/rendered/figure-03/`: historical plots.
- `../../results/`: immutable historical execution records.

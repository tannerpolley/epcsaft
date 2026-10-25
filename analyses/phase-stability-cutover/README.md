# Native tangent-plane cutover

**Frozen record — not rerunnable as a current-main calculation.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The historical
trial declaration places a charged component in vapor support and is refused by the
current zero-charge vapor-support rule; later tangent-plane leaves own current
qualification.

The current native results, conditions, changed claims and reproduction command are in
[the research notebook](notebook.qmd). The prose below preserves the historical
`c74a91f6` result and its interpretation; “current Engine” in that record refers
only to that historical wheel. Its charged stability claims are not re-established
by the native replay. The original result files and wheel identity remain unchanged.

# Phase count and stability of the retained phase-stability cases

The superseded project's issue 154 compared a HELD2 phase-discovery baseline
(`87c93528`) with a proposed discovery cutover (`07d99449`) on four cases. It
recorded a supported negative because the cutover returned one phase instead of
three for a neutral ternary and exhausted its budget on Khudaida row 7. Neither
formulation exists in the current Engine. Phase discovery is future scope
(`engine/docs/science/algorithms.md`), and an undeclared topology is rejected
with "phase count unknown; declare phases". This record evaluates the same four
cases with the current Engine as declared topologies. A tangent-plane distance
(TPD) search tests their stability.

## Result

[`scripts/current_engine.py`](scripts/current_engine.py) produced the
[rows](results/current-engine.csv) and [summary](results/summary.json) with
wheel `c74a91f6…`, which includes the test-point fix below. TPD minima come from a composition grid plus a local
Nelder–Mead refinement, on the liquid and vapor density roots. A minimum that
is zero to roundoff is the reference phase itself.

| case | T, P, feed | stable result | min TPD, liquid trials | min TPD, vapor trials |
|---|---|---|---:|---:|
| methane/ethane | 300 K, 1 bar, 50/50 | 1 vapor | `-1.3e-15` | `-1.3e-15` |
| water/Na+/Cl- | 298.15 K, 1 bar, 0.98/0.01/0.01 | 1 liquid | `-7.5e-10` | `+2.94` |
| neutral ternary fixture | 300 K, 1 bar, 1/1/1 | 1 vapor | `+1.87` | `-3.4e-15` |
| neutral ternary fixture | 300 K, 1 bar, 3 declared liquids | converged, 9 iterations, Δμ/RT `2.2e-10`; metastable against vapor | `-2.1e-10` | `-2.04` |
| neutral ternary fixture | 300 K, 10 MPa, 3 declared liquids | **3 liquids**, converged, 9 iterations, Δμ/RT `2.4e-13` | `-1.6e-12` | one root only (the liquid) |
| Khudaida row 7, packet v6 | 293.15 K, 1 bar | 2 liquids, Δμ/RT `8.9e-10` | `-3.7e-11` | `+1.84` |

The neutral ternary fixture is a vapor at its retained conditions. Each pure
component has liquid ln φ of `+0.87` to `+1.19` at 1 bar, against about `-0.02` as vapor, so each is above its
boiling point. The liquid states have large negative vapor TPDs: one liquid
`-2.07` and two liquids `-2.06`. Per mole, the single vapor is `2.05 RT` below
one liquid (`G/RT` = `-1.112` vs `+0.935`). The baseline's three declared liquids
were therefore metastable, not the equilibrium. The cutover's one-phase answer
had the correct phase count, although the retained record does not say whether
that phase was a vapor. The first reason for the old negative was wrong. The
second reason, the cutover's Khudaida budget failure, belongs to a formulation
that no longer exists. The current Engine solves row 7 in 15 iterations to a
stable two-liquid state.

## Fixed: declared three-phase problems were refused

Before this change, the Engine refused the declared three-liquid ternary before
iterating: "rank-deficient declared equations; nonlinear solve not run". Its
structural rank check evaluates the Jacobian at `CompiledSystem::generic_point`
(`engine/native/equilibrium/request_compile.cpp`). That point set the phase-p
compositions proportional to `1 + 0.37 i + 0.11 p`, which is affine in the phase
index, so any three or more phase-composition vectors were coplanar and the
material-balance block over phase amounts dropped to rank 2. On this fixture the
determinant of the phase-composition matrix was `-6.9e-19` and the Jacobian rank
was 17 of 18. The rejection came from the test point, not from the equations. It
applied to every declared problem with three or more phases and free amounts,
vapor–liquid–liquid equilibrium included.

Phases 3 and higher now add `0.05 i^p`, a Vandermonde column that is independent
for up to as many phases as components. Phases 1 and 2 are unchanged, so every
one- and two-phase start is bit-identical. The structural lattice in
`engine/tests/native/specification_checks.cpp` now requires full generic rank for
every accepted square cell that is not all-pinned; that check failed before the
change and passes after it. The three-liquid rows above are the result. At 1 bar
the fixture's liquids are metastable against its vapor. At 10 MPa only one
density root exists (a vapor request returns the liquid root), and the
three-liquid state passes the tangent-plane search on it: the minimum TPD is
zero to roundoff at one of its own phases.

## Inputs and reproduction

The two fixtures are the retained issue-154 inputs, with fields that the current
schema requires made explicit. The ternary now declares `charge_number = 0`.
The electrolyte now declares `ionic_dispersion = exclude-same-sign-ion-pairs`,
the former implicit behavior. In the lineage record at commit `4b990757`,
implicit packet v5 and explicit packet v6 agreed to `2.4e-10` in mole fraction. Methane/ethane uses packet `gross-2001-methane-ethane/3`, which states
the implicit neutral charges of version 2. Khudaida row 7 uses packet v6, the
current successor of v5.

Current reproduction uses the native call and retained output path in the notebook.
The historical grid/refinement implementation has been removed.

Claim limits: declared topologies with a bounded TPD search on liquid and vapor
roots; not a global-stability certificate; no solid phases.

The owner-approved rounded-point comparison uses `scripts/evaluate_historical_trials.py` and retains its rows in `results/historical-supplied-round3/`. Native supplied trials use branch-selected density roots without packing guesses. The notebook reads all comparison values from those retained rows. Historical candidate snapshots and charged trial comparisons remain unverifiable for the reasons retained in `unverifiable.csv`; no global certificate follows. Batch 14 delivered review is Supported in [PR #307](https://github.com/tannerpolley/ePC-SAFT/pull/307); manuscript promotion is not established.

## Historical supplied-trial producer (#720)

`scripts/evaluate_historical_trials.py` is **historical, not rerunnable**.
At #696’s measured commit `335707d59c67cf1f92df12b5d25530bf31d916a6`
(report line 716), it raised `epcsaft._core.DeclarationError`: vapor phase
`phase-0` support contains charged component index `1` (`z=1`); vapor support
must contain only zero-charge components. The current rule refuses this
historical declaration before any solve is accepted. The earlier rounded-point
tables and their unavailable comparisons remain frozen.

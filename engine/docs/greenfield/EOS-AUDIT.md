# EOS layers audit — findings (2026-09-02, greenfield@c860395)

Evidence tags: **measured** (run on this tree today), **derived** (follows from reading the
code or the equations), **guessed** (stated as such). Scripts under this directory reproduce
every measured number: `ab_dispersion.py`, `slope_probe.cpp`, `libr_curve.py`.

## A. The three defects — verdicts

1. **Like-ion dispersion not excluded — confirmed, measured.** At fixed density
   (KCl, 298.15 K, ρ = 55975.4012 mol/m³) the greenfield dispersion term equals the superseded
   engine run with `ionic_dispersion = all-ion-pairs` to twelve digits (−8.852592255706 at
   1 mol/kg; −8.754778549632 at 2 mol/kg) and differs from its `exclude-same-sign-ion-pairs`
   branch (−8.849957737809; −8.744923772855). Against Figiel's digitized LiBr model curve
   (27 points, 2.16–5.85 mol/kg): greenfield 49.19 % AARD, max 61.25 %; the superseded engine
   under the exclusion 0.27 % AARD, max 0.56 %, every point within 1.10× its digitization
   uncertainty; the superseded engine under all-ion-pairs reproduces greenfield's numbers to
   four decimals. Not overstated. It is a carried kernel minus its resolution: `helmholtz.hpp`
   diffs against `9eccbfd:_native/helmholtz.hpp` only in includes, namespace, constants and
   the deleted helpers; the pair loop is byte-identical.
2. **`PressureSlope` is a frozen-site derivative — confirmed, measured** with a native probe
   linking the archive: Φ_VV ratio PressureSlope/Hessian = 1.0000 ethane, 14.51 associating
   binary (313.15 K, 3.3e4), 0.679 aqueous electrolyte, 0.807 nine species; the centred
   difference of the gradient agrees with the Hessian path to 1e-8 relative in all four.
   NaCl (Figiel/2) at 280 K, 5e5 Pa raises `density root refinement did not converge`
   through the public API. Not carried: the old solver called the full `(1+S)` Hessian per
   iteration (`CARRY-OVER.md` §2); the 4-sweep path and its missing correction are Stage 3's.
3. **Anchor bypasses branch selection — confirmed, measured.** Ethane at 280 K/1e6 Pa and
   250 K/5e5 Pa: `phase="liquid"` with a vapour anchor returns the vapour root (475.51,
   257.27 mol/m³) and `phase="vapor"` with a liquid anchor returns the liquid root. Not
   carried: `9eccbfd:_native/density.hpp` has no anchor path at all (it has a 360-interval
   scan per phase); the anchor and its `stable_root_count = 1` are newly written.

## B. Sweep: registered case splits and their selectors

85 registry entries; 55 distinct `// EqID:` annotations under `native/physics` + `native/model`
(58 lines); 83 in the old tree. Every entry whose `latex` has `\begin{cases}` or whose
description/change note uses selector language, with the state of its selector:

| eqid | split | selector in greenfield | evidence |
|---|---|---|---|
| `same_sign_ionic_dispersion` | `cases` on `z_i z_j > 0` **and** a declared choice | **missing** — annotated line is unconditional | measured (A.1) |
| `debye_finite_size_chi` | `cases` on `y_ref < 1e-3` | implemented (`electrolyte.hpp:331`, `SCREENING_SERIES_SWITCH`), mirrored in the tape cache key — **under a second policy name** (C.1) | derived |
| `ionic_packing_diameter` | charged → declared packing diameter, else eq 3 | implemented (`helmholtz.hpp:92-101`); consistent with `u_ii = 0` only through packet data (0.88 σ) | derived |
| `born_radial_charging_work` | limits selected by `c_shell`, `c_dielectric` data | implemented; builder requires both under Born | derived |
| `component_mixed_bulk_permittivity` | with/without a neutral component | implemented in-kernel (`electrolyte.hpp:164-180`) | derived |
| `log_fugacity_coefficient` (old-only) | defined only for Z > 0 | implemented unannotated (`helmholtz_model.cpp:524`) | derived |
| `residual_heat_capacity_cp` (old-only) | conditioned on dP/dρ > 0 | implemented unannotated (`:540`) — its selector input is the defect-2 quantity in the density solve, correct here (Hessian path) | derived |
| `practical_osmotic_coefficient` | single active solvent | implemented (solvent-index detection, duplicated — C.6) | derived |
| `reference_molality_sequence` (old-only) | terminal-pair selection | implemented as `reference_limit_terminal_pair` (unregistered — C.3) | derived |
| `declared_solvent_source_reference_transfer`, `neutral_reference_stoichiometric_basis` (old-only, `cases` on `z_i = 0`) | not yet written — owed by Stage 6's standard-state transform | derived |

**One registered case split has no selector, and it is the one that failed.** The permittivity
three-way choice (ion-fraction / ion-specific / solvent-only) is the healthy counter-example:
Stage 2 merged three equations into one and Stage 4's builder resolves the choice into
per-component data (`parameter_builder.cpp:364-375`) — the same "resolution wraps kernel"
shape, done right, because α is expressible as data with no diagonal problem.

Carried-kernel diffs against `9eccbfd` (measured with `git diff`): `helmholtz.hpp` 17+/94−
(all mechanical: includes, namespace, constants, deleted helpers; the equations are identical);
`electrolyte.hpp` 57+/110− (the permittivity merge, scalar dispatch, policy names — the
Debye–Hückel and permittivity equations are identical); `born.hpp` 10+/10− (names only);
`association.hpp` vs `a248de6` 62+/95− (Eigen `solve_linear`, policy names, iteration
out-parameter). No other carried kernel had a consumer-side selector in `_resolve.py`: `l_ij`
optional-default-0 and `k_ij` required-for-unlike-pairs are reproduced by the builder;
the domain limits (`T`/`P` range, `ion_fraction_max`, `k_ij(T) ≤ 1` over the range) and the
`missing-opposite-charges` / `enabled-without-ions` checks were dropped without a record —
validation decisions, not physics, listed under C.

## C. Borrowed and newly written crud, by file and line

1. **Duplicate policy entry, two names one meaning** — `numerical_policy.hpp:136`
   `ELECTROLYTE_SERIES_THRESHOLD` (1e-3, added in the Stage 3 commit) and `:242`
   `SCREENING_SERIES_SWITCH` (1e-3, added in `155fa63`). `helmholtz_model.cpp:295` keys the
   tape cache on the first; `electrolyte.hpp:334` branches on the second. Diverging either
   reuses a tape across the series switch with no error. Newly written; exactly the #186
   defect class. Checker rule 5 cannot see it (matches `constexpr`, not `inline const
   Tolerance`). Measured (grep), consequence derived.
2. **Dead double-path evaluation** — `helmholtz_model.cpp:324-338, 461`: `electro`, `born`,
   `direct`, `direct_value` computed on every `evaluate` and discarded (`(void)direct_value`).
   Three kernel evaluations per call, ~41 per cold density scan. Newly written. Derived.
3. **Annotations without registry entries** — `ion_suppressed_bulk_permittivity`
   (`electrolyte.hpp:183`) and `reference_limit_terminal_pair`
   (`infinite_dilution_reference.cpp:374`) are not in `equations_registry.yaml`; the registry
   still carries their predecessors (`ion_fraction_…`, `ion_specific_…`,
   `reference_molality_sequence`). The registry sync (`engine/scripts/docs/sync_equation_registry.py`)
   does not check the source tree, so nothing flags it. Measured.
4. **Unreachable permittivity forms and five dead parameter arrays** —
   `PermittivityForm::DensityTemperature` and `::DipolarPerturbation` (`parameters.hpp`,
   `electrolyte.hpp:140-162`) can never be selected: `parameter_builder.cpp:163-171` maps
   packet choices to `None`, `IonSuppression`, `ComponentMixing` only, and no entry field
   fills `dielectric_volume_m3_per_mol`, `dielectric_temperature_k`, `dipole_scaling_debye2`,
   `polarizability_scaling_angstrom3`, `correlation_integral_parameter`. Carried
   (`CARRY-OVER.md` §1 kept "four published forms"), reachable only from a hand-built struct.
   Derived. Either the builder learns the two forms or the forms and arrays go.
5. **Hidden physical default in Python** — `eos.py:175` `stoichiometry: tuple[int, int] = (1, 1)`
   on `Mixture.reference`. A 2:1 salt with the default silently computes a 1:1 reference.
   Newly written. Derived. (`parameters.py:120`'s `charge_number = 0` when absent is a schema
   convention exercised by 82 of 104 neutral components across the packets — measured — and
   matches the old `_dictionary_io` optional `fixed` block; acceptable, noted.)
6. **Duplicated solvent-index detection** — `infinite_dilution_reference.cpp:408-416` and
   `:477-485`, same nine lines twice. Newly written. Measured (diff).
7. **Same-sign k_ij records required by the builder** — `parameter_builder.cpp:346-355` requires
   a `k_ij` record for every unlike pair including same-sign ion pairs, so selecting two
   cations from `figiel-yu-held-2025/2` (21 absent) is rejected by name today; the old loader's
   `not_applicable` set is the missing decision. Carried-minus-resolution, same root as A.1.
   Derived.
8. **Validation decisions dropped without a record** — the old `_resolve.py` rejected a
   neutral component carrying `packing_diameter`/`born_diameter` (`ionic-component-required`),
   an electrolyte choice with no ions, and enforced the source domain limits; the builder
   accepts the first two silently (a neutral's packing diameter is stored and never read —
   `segment_diameter` takes the σ branch) and has no domain concept. Not physics; listed so
   the decision is on record. Derived.
9. **Policy entries sourced to the old tree only** — `REFERENCE_CONVERGENCE_TOLERANCE 5e-5`,
   `REFERENCE_MOLALITY_*`, `REFERENCE_CONVERGENCE_WITNESS`, `COMPOSITION_SUM_TOLERANCE`
   cite "electrolyte_reference.hpp" as their source — a file, not a justification. Category 3
   in form, category 4 in substance. Carried. Guessed as to whether a paper value exists.
10. **Stale numbering** — `helmholtz_tape.hpp:16` "TRANSFER.md, Stage 2" (now Stage 3);
    `density.cpp:300` "Source: density.hpp on main". Cosmetic.
11. **Phase concept in `model/`** — `DensityBranch` is the one sanctioned place;
    `infinite_dilution_reference.cpp:364,394,421,448` hard-code `Liquid`, which is the MIAC
    reference's definition (a liquid solution) rather than a route flag. No action.
12. **Branches on a problem kind** — none found in `model/` beyond `charge_number != 0`
    (physics) and the declared enums. `parameters.py:203` rejects `family == "k_hb_ij"` by
    name: a dropped decision (`cross_association_energy_correction`, registered, no
    successor); no packet carries the family as a coefficient (the `ascani-2022` hits are
    provenance locator text). Measured.

Retracted during the audit: an early read that the ascani packet is unloadable because of
`k_hb_ij` — wrong; see 12.

## D. Repair design — see `stage-04.repair.brief.md` §2

Placement, in one line each: the selector is an enum in `physics/parameters.hpp` and a branch
on the annotated line of `physics/helmholtz.hpp` (the registry says the split is on ε_ij, so
that is where it lives, and every consumer — fixtures included — inherits it); the
declaration is required by the model constructor (no default) and mapped by the builder from
a packet key on the electrolyte family record (where the old schema and the registry both put
it); the slope fix reuses the S×S factorisation the Hessian order already forms; the anchor
is honoured only with a declared branch and a checked far side; eight packets gain a
revision declaring the exclusion.

**The `k_ij = 1` mechanism is challenged, with evidence.** It zeroes the term exactly, but it
must also zero the temperature slope (`helmholtz.hpp:237-241`) and it leaves a same-sign
`k_ij` activatable by `ParameterFamily::Kij` with a non-zero derivative, so a fit could
"discover" like-ion dispersion. The branch on ε_ij has neither hole: an excluded pair's `k_ij`
and slope are unread and its derivative is identically zero. Derived.

**Default: none.** Every source in the tree is the Held lineage and would declare the
exclusion, so a default would be right nine times today — and wrong the first time a Ji 2012
set arrives, silently. The cost of no default is eight packet revisions, which Stage P will
migrate anyway.

## E. Stage numbering

One brief, `briefs/stage-04.repair.brief.{md,json}`, branch `stage-04-repair`, re-accepting
Stages 2, 3 and 4 **together**: their gates are one command set run on one tree, the physics
property checks cannot see the change (they hold for any u_ij), and only Stage 4's oracle gate
sees the convention. `TRANSFER.md`'s execution table gets a "repaired in 04.repair" note on
rows 2 and 3.

## F. Prevention — recommended: a registry-to-selector rule in the checker

Registry entries whose `latex` contains `\begin{cases}` gain a `selector:` field naming the
symbol that chooses the branch (`ionic_dispersion`, `SCREENING_SERIES_SWITCH`, `charge_number`).
Checker rule 8: for every such entry, the file carrying its `// EqID:` annotation must
reference the selector symbol within the annotated function (±40 lines is enough here).
Companion: the set of EqIDs annotated at `source_commit` minus the set annotated in the
target tree must be listed in `CARRY-OVER.md` under "no successor" — today 30 are missing
without a record. Cost: ~35 lines in `check-engine-architecture.py`, one field on 2–4
registry entries, one line in `sync_equation_registry.py`'s `FIELDS`, ~15 lines for the
successor rule. It would have gone red at Stage 2 on `same_sign_ionic_dispersion`.

What it would not catch: (i) a split written in prose rather than `\begin{cases}`
(`ionic_packing_diameter`'s dependence on `u_ii = 0` is a change note); (ii) a branch on the
right symbol choosing the wrong pairs — presence, not correctness; (iii) defects 2 and 3,
which are not case splits: the derivative-order defect is caught only by a finite-difference
check on **every** `DerivativeOrder` on a fixture **with sites** (12 lines, in the brief), and
the anchor defect only by a branch check with an adversarial anchor (in the brief). Both belong
in `TESTING.md` as standing rules: every derivative order against a centred difference on an
associating fixture; every warm start against its cold answer.

Costed alternatives: a required diff-against-source reading step is what `CARRY-OVER.md`'s
audit already was, and it missed this — not mechanical, not recommended as the mechanism. A
fixture that flips every enum and asserts the output moves (`same_sign_selector_live`, ~40
lines, in the brief) catches a dead selector but only once the choice is reified as an enum —
it is the companion of rule 8, not a substitute.

## Missing facts, named

- Figiel's own reported per-salt deviations (the published AARD) — still missing; the
  digitized curve stands in for the model, with its own 0.45–0.8 % uncertainty.
- Whether a paper value exists for `REFERENCE_CONVERGENCE_TOLERANCE` (C.9).
- Whether Stage P's provenance form should carry `ionic_dispersion` for the eight revisions
  instead of `repository-migration-v1` — a contract decision for the coordinator.

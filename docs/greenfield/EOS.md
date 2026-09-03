# EOS — the equation kernels, the tape, the state and the density closure

Stage 4's design document. Four documents merged here on 2026-09-02, each keeping its own
section: the audit that swept the layer, the two measured defects, and the performance floor.
Evidence is **measured**, **derived** or **guessed** throughout.

---


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
   under the exclusion 0.27 % AARD, max 0.56 %, every point within 1.097× its digitization
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
   iteration (`CARRY-OVER.md` §2); the 4-sweep path and its missing correction are Stage 4's.
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
| `declared_solvent_source_reference_transfer`, `neutral_reference_stoichiometric_basis` (old-only, `cases` on `z_i = 0`) | not yet written — owed by Stage 5's standard-state transform | derived |

**One registered case split has no selector, and it is the one that failed.** The permittivity
three-way choice (ion-fraction / ion-specific / solvent-only) is the healthy counter-example:
Stage 4 merged three equations into one and Stage 3's builder resolves the choice into
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
   `ELECTROLYTE_SERIES_THRESHOLD` (1e-3, added in the Stage 4 commit) and `:242`
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
10. **Stale numbering** — `helmholtz_tape.hpp:16` "TRANSFER.md, Stage 4" (now Stage 4);
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
Stages 3 and 4 **together**: their gates are one command set run on one tree, the physics
property checks cannot see the change (they hold for any u_ij), and only Stage 3's oracle gate
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
successor rule. It would have gone red at Stage 4 on `same_sign_ionic_dispersion`.

What it would not catch: (i) a split written in prose rather than `\begin{cases}`
(`ionic_packing_diameter`'s dependence on `u_ii = 0` is a change note); (ii) a branch on the
right symbol choosing the wrong pairs — presence, not correctness; (iii) defects 2 and 3,
which are not case splits: the derivative-order defect is caught only by a finite-difference
check on **every** `DerivativeOrder` on a fixture **with sites** (12 lines, in the brief), and
the anchor defect only by a branch check with an adversarial anchor (in the brief). Both belong
in `CANON.md` as standing rules: every derivative order against a centred difference on an
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

---


Status: findings, 2026-09-02. Written because Stage 3's activity gate passed while the
engine's mean ionic activity coefficients drifted low at high molality, and the first
account of that drift — recorded in `STAGE-LOG.md` and corrected below — was wrong in two
places.

Evidence is **measured** (run on this tree today), **derived** (follows from the equations),
or **guessed**.

## 1. The equation, and where the number comes from

Every ePC-SAFT electrolyte term reads a bulk relative permittivity that falls as salt is
added. The engine spells that with one equation (`EqID: ion_suppressed_bulk_permittivity`,
`physics/electrolyte.hpp`):

```text
ε = ε_solvent-mixture / (1 + Σ_i α_i x_i)
```

`α_i` is `ion_suppression_coefficient[i]`, and the sum runs over ions in ion mole fraction.
This is Zuber et al. 2014 Eq. 9 (Fluid Phase Equilibria 376, 116–123). **One fitted `α` per
ion per solvent** — Table 2 gives eighteen for water, Table 3 eight for methanol, Table 4
four for ethanol, and Table 9 a general cation and anion value per solvent for ions the
paper does not fit. Now in the tree as
`data/packets/zuber-2014-electrolyte-permittivity/1`.

Two published values are negative (Ca²⁺ −0.53, SO₄²⁻ −0.66 in water) and are retained as
fitted: the paper explains that the dielectric constant of water *rises* with MgSO₄ to a
point before falling, which Eq. 9 can only represent with a negative parameter.

**The engine has always supported per-ion values.** `parameter_builder.cpp` accepts the
`ion-specific-suppression` choice and requires a coefficient for every charged component
(`:374`); the `ion-fraction-suppression` choice is the one-number form, read from
`model_coefficients` and broadcast to every ion (`:366-370`). What was missing was the data,
not the capability.

## 2. Why Figiel uses 7.01 for every ion, and why that is not a free parameter

`figiel-2025-reference-electrolytes/2` declares `ion-fraction-suppression` with a single
`7.01`. That is not an arbitrary placeholder to be swapped. **measured**, from
`analyses/reproduction/2025-figiel/scripts/check_figiel_2025_source_packet.py:275-293` and
the parameter-provenance ledger it checks: the Born diameters of Li⁺, Na⁺, K⁺, Cl⁻ and Br⁻
are classified `fitted-in-Figiel-2025` and are staged **after** the dielectric suppression
stage, whose single reported value is `7.01`.

So the ion diameters in that packet were regressed **with `α = 7.01` held fixed**. The two
are one calibration. Changing `α` alone leaves diameters that were fitted to absorb its
error.

> **The rule.** `born_diameter` and the dielectric suppression coefficient are **paired**.
> A packet's Born diameters are only valid with the suppression coefficient they were
> regressed against, and this packet's are only valid with `7.01`. Never evaluate
> `figiel-2025-reference-electrolytes` with any other suppression value, and never
> substitute a suppression coefficient into a packet without re-fitting the diameters in
> the same step. Any consumer that changes one must change both, and a fit that varies one
> while holding the other must say so.

## 3. The experiment

Per-ion Zuber values substituted into the Figiel packet, nothing else changed, mean ionic
activity coefficient against the digitized literature of Figiel Figure 5, `m ≤ 2 mol/kg`,
21 points per salt (**measured**, `analyses/reproduction/2025-figiel/figures/figure-05/`):

| salt | AARD, α = 7.01 | AARD, Zuber per-ion |
|---|---|---|
| LiCl | 6.23 % | 7.02 % |
| LiBr | 6.83 % | 6.97 % |
| KCl | 5.29 % | 9.11 % |
| KBr | 7.54 % | 10.62 % |
| NaCl | 3.34 % | 6.89 % |
| NaBr | 2.18 % | 4.58 % |

**Physically correct per-ion coefficients make the fit worse on every salt.** That is the
expected outcome of §2 and is evidence for it, not against per-ion support: the diameters
are absorbing the difference. Using Zuber's values properly means re-regressing the Born
diameters against them, which is Stage 6's kind of work, on Stage 6's machinery.

## 4. The mechanism, and the direction of the error

A mean ionic activity coefficient falls steeply from 1 at low molality — Debye–Hückel — then
passes through a minimum and rises. The rise is what the Born term does when the permittivity
itself depends on composition. The Born contribution to a species chemical potential goes as
`−(z_i² e²/8πε₀ d_i)(1 − 1/ε)`; as salt is added, `ε` falls, `(1 − 1/ε)` falls, the Born
energy becomes less negative, and `γ±` is pushed back up. **Stronger suppression means a
stronger upturn.** `derived`, and confirmed by the experiment: lowering K⁺'s coefficient from
7.01 to Zuber's 2.43 weakens the upturn and moves KCl *further* below the data, not toward it.

Against the digitized literature, every deviation is one-signed and grows with molality
(**measured**, α = 7.01):

| salt | 0.001 m | 0.1 m | 0.5 m | 1.0 m | 2.0 m |
|---|---|---|---|---|---|
| KCl | −0.0 % | −0.5 % | −3.0 % | −7.8 % | −19.3 % |
| KBr | −0.0 % | −1.0 % | −5.2 % | −11.4 % | −24.9 % |
| LiBr | −0.0 % | −1.5 % | −4.6 % | −9.7 % | −23.6 % |
| NaCl | −0.0 % | −0.9 % | −1.9 % | −4.5 % | −12.4 % |

Exact at infinite dilution, and the upturn is too weak. Whether the cause is the Born
diameters, the suppression, or something else is a question for the planner; §5 says what has
already been ruled out.

## 5. Ruled out: the mapping

Every field this model reads is mapped and assigned by Stage 3's builder, checked line by
line: `ion_suppression_coefficient` per component (`:227`) and from the model-level record
(`:366-370`), `ionic_region_relative_permittivity` (`:292`), `born_c_shell` and
`born_c_dielectric` (`:297-298`), each rejected by name when the declared choice requires it
and it is absent (`:359-363`). The physics layer reads all four (`physics/born.hpp:42-76`).
The drift is not a transcription loss between packet and array.

## 6. Correcting the first account

`STAGE-LOG.md` first recorded this as the engine being 1.8–2.9× worse than the superseded one
and missing a qualitative feature. Both halves were wrong, for one reason: that comparison
used the six oracle rows against
`analyses/reproduction/2025-figiel/results/data/figiel-2025-miac.csv`, and **that file's rows
carry a different `parameter_fingerprint` per salt** — it is a per-salt fitted comparison, not
one model. Against the one consistent superseded run in
`figures/figure-05/results/data/figure-05-series.csv`, over 21 points per salt:

| salt | superseded | this engine |
|---|---|---|
| LiCl | 6.09 % | 6.23 % |
| LiBr | 7.10 % | 6.83 % |
| KCl | 6.10 % | 5.29 % |
| KBr | 5.26 % | 7.54 % |
| NaCl | 6.58 % | **3.34 %** |
| NaBr | 8.35 % | **2.18 %** |

Comparable or better on four of six. And at low molality the superseded engine returns
**exactly 1.0000** for NaCl and NaBr below 0.05 mol/kg — no electrostatics at all, up to 14.7 %
wrong — where this engine matches the literature to 0.02–0.34 %. The "missing upturn" claim
was also wrong: LiCl and LiBr are the only salts with an upturn inside 2 mol/kg and this
engine reproduces both, too weakly. KCl and KBr fall monotonically in the *data* as well;
their minimum is near 2.5–3 mol/kg, outside the oracle's range.

## 7. What the gate measures, and what it should

The gate asks whether the engine reproduces **the measured activity coefficients** within the
parameter source's own adequacy criterion. That is a fact about the world and a legitimate
oracle claim. It is not the question the reproduction is trying to answer, which is whether
the engine reproduces **the model Figiel published**.

Nothing in the tree currently answers the second. The digitized files under
`figures/figure-05/source/` are literature points; `figure-05-series.csv`'s model rows are the
superseded engine's own output, and `inputs/docs/legacy-derived-gates.json` keeps only scalar
per-figure bounds derived from it (figure 5: RMSE ≤ 0.151, max ≤ 1.276) with the model series
deliberately not retained. **Missing fact, named:** Figiel's own reported deviations on these
salts, per salt, from the paper. With them, the reproduction claim becomes "our AARD is within
a stated margin of the published AARD" — a comparison against the paper rather than against a
superseded implementation.

Two things follow for the briefs, neither of which is a change to Stage 3:

1. An average over a molality range cannot see a systematic drift that is exact at one end.
   A criterion of the form "AARD ≤ x" should be paired with a **shape** condition — the sign
   of the deviation must not be constant across the range, or the highest-molality row must
   be within its own bound.
2. Per-ion suppression is data plus a re-regression, not a substitution. The packet is in the
   tree; the fit that would make it usable belongs to Stage 6.

---


Status: findings, 2026-09-02, from a measured sweep of 630 cases (5 fixtures × 7 temperatures
× 6 pressures × 3 branches) plus three real packets. Independently verified on the two claims
that carry the rest. Evidence is **measured**, **derived** or **guessed**.

## 1. The central defect: `PressureSlope` returns a frozen-site derivative

`DerivativeOrder::PressureSlope` sets

```cpp
result.volume_second_derivative = second[2 * 1 + 1];        // helmholtz_model.cpp:385
```

which is the raw tape second derivative at **frozen site fractions**, while the `Hessian`
order sets

```cpp
result.volume_second_derivative = result.hessian[1 * total_coordinates + 1];   // :460
```

which is the Schur-corrected total derivative `base(1,1) − cross_Vᵀ H_X⁻¹ cross_V`, including
`dX/dV`. **Verified by reading both branches.** Central differences on `P` confirm which is
which — associating binary at 313.15 K, ρ = 33000: finite difference 70463.61706, Hessian
path 70463.61705, `PressureSlope` path 40494.85915.

`solve_density` consumes `PressureSlope` for both its Newton derivative and its
mechanical-stability filter, so with association active it is iterating on the wrong slope.

**Consequences, measured:**

- The Newton degenerates to a linear iteration at rate `|1 − exact/frozen|`. Predicted 0.7358
  for the 313.15 K liquid root; the trace measures 0.735 per step over 60 steps. Steps to
  reach 1e-12: 90 associating, 99 aqueous, 118 nine-species — against a
  `DENSITY_MAX_ITERATIONS` budget of 60.
- **102 of 630 sweep cases throw** `density root refinement did not converge`. **Verified on a
  real packet through the public API:** `figiel-2025-reference-electrolytes/2`, NaCl at
  280 K and 5 × 10⁵ Pa, raises `EosError: density root refinement did not converge`.
- The stability filter tests a quantity up to 368× wrong. Its *sign* never disagreed across
  471 converged cases, so no returned density is currently wrong from this — luck, not a
  guarantee.
- With the exact slope the same brackets converge in 4–8 evaluations.

**Why no gate caught it.** `model_checks.cpp:227` finite-difference-checks
`State::pressure_density_derivative`, which comes from the *Hessian* path and is correct. The
only fixed-pressure native test (`:276`) uses `model_support::ethane()` — **non-associating**,
where the frozen and exact slopes are identical by construction. `PressureSlope`'s
`volume_second_derivative` has no finite-difference check anywhere. The suite is not weak; it
tests the path that works and never exercises the one that does not.

## 2. It can miss the liquid root

At **ethane, 280 K, 10⁴ Pa** the reference (1200-interval bisection) finds three roots:
4.299 (dP/dρ = +2324), 9567.8 (−461), 11163.5 (+645). The two dense roots are a factor 1.167
apart and one scan interval is 1.985× wide, so both fall inside a single interval, the two
sign changes cancel, and neither is seen. `phase="liquid"` returns the **vapour** root with
`stable_root_count = 1` and no diagnostic. 200 intervals are needed to find it.

`DENSITY_BRACKET_INTERVALS = 40` is inherited rather than tuned; the miss count against the
reference is non-monotone in granularity — 10 → 26 misses, 20 → 2, 30 → 4, 40 → 1, 200 → 0 —
and polish cost is flat above 20 intervals, so a finer scan buys root coverage and nothing
else.

## 3. The anchor path bypasses branch selection

Requesting `Liquid` with a vapour-side anchor returned the **vapour** root at **57 of 57**
two-root conditions, labelled `branch = Liquid` with `stable_root_count = 1` — a count it
never computed. Example: associating binary, 450 K, 10⁵ Pa, liquid root 29013, anchored
answer 27.05.

This is why threading more anchors through the callers (§5) must wait for the fix: more warm
starts would mean more silently wrong branches.

## 4. Cost

Per cold solve, `OMP_NUM_THREADS=1`:

| fixture | evaluations | scan | polish | wall |
|---|---|---|---|---|
| ethane | 46 | 41 | 5 | 273 µs |
| neutral binary | 46 | 41 | 5 | 392 µs |
| associating binary | 99 | 41 | 58 | 2007 µs |
| aqueous electrolyte | 90 | 41 | 50 | 1915 µs |
| nine species | 86 | 41 | 45 | 6437 µs |

**This corrects `EOS.md` Mechanism 4**, which attributes the cost to the scan.
On neutral fluids polishing is 11 % of evaluations; on associating fluids it is **60–66 %**.
The scan is not the dominant term on the fluids that matter — the broken Newton is.

Through Python: nine-species cold `Mixture.state` is 8024 µs at 108 density evaluations;
with an exact anchor, 460 µs at 1 evaluation. `Reference.activity(1.0)` is 4.04 ms — one
fully cold solve per activity call.

## 5. Ranked, with what each changes

| # | change | saving | changes an answer? |
|---|---|---|---|
| 1 | Give `PressureSlope` the true `dP/dρ` — the same Schur correction the `Hessian` order already forms | polish 27–79 → 4–8 evaluations; nine-species solve 7.6 → 3.2 ms | **Yes.** Turns 102 throws into answers and changes what the stability filter accepts |
| 2 | Scan at `DerivativeOrder::Gradient` — the scan reads only `residual`, which needs 2 sweeps, and asks for 4 | 1.36 ms per solve; a further 43 % after #1. Answers bit-identical on 20 of 20 rows | No |
| 3 | Make the anchor path respect `branch` and stop asserting a stable-root count it never computed | correctness only | **Yes**, explicitly |
| 4 | Thread anchors through `Pure.state`, `electrolyte_activity`, and the two cold `state_at_pressure` sites in `infinite_dilution_reference.cpp` | 2.9–4.7× end to end | No — **but only after #3** |
| 5 | Seeded two-branch Newton (`P/RT` for vapour, 0.6 × the hard-sphere bound for liquid) with the scan as fallback | 41 + 27–79 → about 9–15 evaluations; 75 of 75 at two-root conditions | **Yes** — the scan enumerates roots and can report ambiguity; two seeds cannot prove there are only two. Argue separately |

Changes 1–4 need no new policy entry and no literal. Change 5 needs one seed fraction,
expressible as a multiple of `CLOSE_PACKING_FRACTION` through the existing
`molar_density_upper_bound`; it is a seed, not a tolerance, so it does not sit under
`DENSITY_ROOT_RESIDUAL` in the nesting order.

## 6. What this says about `CORE.md` §3.3

The `Convergence` record sketched there wants `amplification = 1/|dr/dρ|` — **which is exactly
the quantity defect 1 corrupts.** It cannot be built correctly until that is fixed.

And the floor is **not** currently the binding constraint: with a correct Newton the residual
reaches 1e-16 to 5e-13 everywhere tested, and over the 471 converged cases the achieved
reduced residual is p50 1.17 × 10⁻¹³, max 9.97 × 10⁻¹³. `converged_to_floor` would never
fire here. The 102 failures are a real non-convergence, not a mislabelled floor.

`DENSITY_MAX_ITERATIONS`'s own source string — "5-8 steps from an anchor, 20-40 for a cold
bracket" — is falsified by measurement today and becomes true after change 1. The value needs
no change.

## 7. Delivered accuracy, stated

The tolerance is on the **reduced** residual, so the delivered relative pressure error is
`|r|/z`, and on a dense liquid root `z ≈ 7e-5…7e-4`. A reduced 1e-12 therefore buys only about
1e-8 relative. Measured worst cases: 1.23 × 10⁻⁸ (MEA nine-species), 1.70 × 10⁻⁸ (Figiel
NaCl), 6.6 × 10⁻¹³ (methane–ethane). Neutral fluids sit at p50 3 × 10⁻¹⁶; associating and
electrolyte fluids at p50 3.3–3.7 × 10⁻¹². This is worth stating in the policy header's
rationale, because "1e-12" reads tighter than what a caller actually receives.

## Not measured

- The association-noise contribution to the residual floor (`CORE.md` §3.3's
  `V · floor_grad[V]`): `solve_site_fractions` does not expose its converged residual, so the
  floor was bounded from above but not measured.
- The isolated wall cost of the S×S correction in change 1 — bounded below one sweep of noise
  by difference, not timed directly.
- Composition. The sweep varies temperature and pressure at one composition per fixture, so
  the two-in-98 real-packet failure rate is a **lower bound**.

---


> **The scripts that produced this document are deleted, 2026-09-02.**
> `profile_entry_points.py` and `capture2.py` both import `epcsaft.equilibrium`,
> `epcsaft.regression`, `unit_registry` and `_core._instrumentation`, none of which exist in
> this tree, and `profile_entry_points.py` hard-coded an absolute path into the superseded
> checkout. They ran against `perf/instrumentation`, a branch of the tree being replaced.
> **The measurements below are facts and stay**; the scripts are recoverable from history,
> and `measure_floor.py` — rewritten at Stage 3 against the current surface — is the one that
> still runs. A canon file that cannot run is a rule that has quietly stopped checking.


Measured 2026-09-02 on branch `perf/instrumentation` (`f0abace`, `main @9eccbfd` plus counters only), wheel `287ad813144acc856381d64ca43177e061531f28cec605b7a036a3943b34f57b`, `OMP_NUM_THREADS=1`, gcc 13.3.0, AMD Ryzen 5 5500, Python 3.13, `ldd` in `floor-measurements.json`. Script: `measure_floor.py`. Counters: tape recordings and recording time, tape operations (`ADFun::size_op`), Jacobian/Hessian/Forward calls, sweep time, association solves and iterations, exposed as `_core._instrumentation()`.

## 0. The premise, tested first

The estimate this replaces assumed re-taping dominated the 33 ms reactive evaluation. It does not.

| case | wall | recording | sweeps (Jacobian + Hessians + Forward) | Hessian calls | assoc. iterations |
|---|---|---|---|---|---|
| neutral binary, `Mixture.state(rho)` | 660 µs | 17.8 µs (2.7%) | 49 µs (7.5%) | 1 | 0 |
| associating binary, `Mixture.state(rho)` | 1480 µs | 37 µs (2.5%) | 791 µs (53%) | 5 | 126 |
| aqueous electrolyte (3), `Mixture.state(rho)` | 1006 µs | 31 µs (3.1%) | 397 µs (39%) | 3 | 14 |
| nine-species electrolyte, `Mixture.state(rho)` | 8494 µs | 100 µs (1.2%) | 7315 µs (86%) | 7 | 44 |
| nine-species, GREPE per-evaluation path (`evaluate_parameterized_phase`) | 4815 µs | 124 µs (2.6%) | 4572 µs (95%) | 8 (two tapes) | 44 |
| reactive bubble point, whole solve (4 cases) | 10.9–11.4 s | 160–180 ms (1.5%) | 10.6–11.0 s (97.5%) | 11,041–11,547 | 69,649–74,020 over ~1,600 solves |
| neutral binary stability decision | 671 ms | 68 ms (10%) | ~185 ms (est., small-tape sites not timed) | 4,632 | 0 |

Two consequences. First, one cached tape per topology is worth 1–10% of an evaluation, not the factor the estimate assumed; the premise was wrong. Second, the reactive solve reports 183–357 "EOS evaluations" but performs 1,700–1,830 tape recordings and 11,000–11,500 Hessian sweeps: about 10 tapes and 63 dense Hessians per counted evaluation. The cost is the Hessian count, and the count is set by two multipliers stacked: `(S + 1)` Hessians per tape to assemble the implicit association second derivative, and roughly 10 tapes per counted evaluation from re-evaluating the same state across the objective, gradient, Hessian, vapour-closure, and derivative-audit callbacks.

The `Mixture.state` path adds about 590 µs of Python and pybind marshalling on top of 67 µs of native work for the neutral binary; that is the 50× the earlier note attributed correctly to marshalling.

## 1. Arithmetic floor for one residual Helmholtz value

The operation count is measured, not estimated: CppAD records one tape operation per elementary scalar operation, and `size_op` reports it.

| case | tape ops (value + mass-action residuals) | transcendental ops (counted from the kernels) | association fixed-point iterations | double-precision floor |
|---|---|---|---|---|
| n = 2 neutral | 301 | 2 `exp` (diameters), 4 `sqrt` (pair energies), 3 `log` (hard sphere, 2 chain contacts) | 0 | 0.4–0.8 µs |
| n = 2 associating, S = 4 | 601 | as above + 16 `exp` (site-pair strengths) + 4 `log` | 126 (damped substitution, 0.5) | 1.5–4 µs incl. the solve |
| n = 3 electrolyte, S = 2 | 572 | + 1 `sqrt` (κ), 2 `log` (χ), permittivity correlation | 14 | 1–2 µs |
| n = 9 electrolyte, S = 6 | 2164 | 9 `exp`, 81 `sqrt`, 11 `log`, 36 `exp` (site pairs), 6 `log` (χ) | 44 (1,584 site-pair MACs per iteration) | 5–12 µs |

Scaling read off the kernels (`helmholtz.hpp:215-357`): the pair loop is `n²` with a `sqrt` per pair and about 20 flops; moments and chain terms are `n` with a `log` per component; the seven-term integrals and `C1` are constant (~120 flops); association strengths are `S²` with an `exp` each and the residual is `S²` per iteration; Debye–Hückel is a `sqrt` plus a `log` per ion; Born is `O(n_ion)`.

Throughput assumptions (line items):

| assumption | value used | sensitivity |
|---|---|---|
| scalar double throughput, compiled, no vectorisation (the loops are short and data-dependent) | 1 flop / 0.5–1 ns | linear in the floor |
| transcendental cost (`exp`, `log`, `sqrt` via libm) | 15–40 ns each | the 9-species value has ~140 transcendentals: 2–6 µs of its 5–12 µs; dominant at n = 9 |
| cache | everything (parameters, ~50 KB tape) resident in L1/L2 | if not, ×2–3 on the upper end |
| association solve | Newton with `H_X` reaches 1e-13 in 3–5 iterations from the previous point (today 14–126 damped-substitution iterations) | today's 126 iterations at S = 4 are ~15 µs on their own; the floor assumes ≤5 |

CppAD's interpreted forward sweep runs the same 2164 ops at 14 ns per op measured (below), so a value through the tape costs ~30 µs against 5–12 µs compiled.

## 2. Derivative multiplier

Measured per-op cost on the sweep path: nine species, 7 Hessians in 7,315 µs over a 2164-op tape of width `w = 2 + n + S = 17` (`T, ρ, x, X_sites`). A dense CppAD Hessian is one forward-0 plus `w` (forward-1 + reverse-2) pairs, `1 + 2w = 35` sweeps; 7 × 35 = 245 sweeps in 7,315 µs → 30 µs per sweep → **14 ns per tape op per sweep**. Recording: 100 µs / 2164 ops → 46 ns per op recorded. Neutral binary (301 ops, w = 4): 9 sweeps ≈ 38 µs, consistent with the 49 µs measured for Jacobian + Hessian + Forward.

Relative to one value sweep (30 µs at n = 9):

| quantity | sweeps | multiple of the value | today at n = 9, S = 6 |
|---|---|---|---|
| value | 1 | 1 | 30 µs |
| gradient (reverse, one dependent) | 1 + 1 | 2 | ~60 µs; today the Jacobian has `1 + S` dependents |
| dense Hessian, `w` independents | `1 + 2w` | 23 at `w = 11` (sites eliminated), 35 at `w = 17` | 1.0 ms each, 7 of them: 7.3 ms |
| parameter Jacobian, P active | included in the same Hessian when parameters are tape inputs | +2P sweeps | — |

The implicit association second derivative today: `(S + 1)` dense Hessians of the full tape, then `X_qr` assembled by linear solves. With one Hessian of the Lagrangian `a + λ·H` (λ from one linear solve of `H_X`) the same tensor costs one Hessian. Multiplier: `(S + 1)` → 1 (7× at S = 6, 5× at S = 4), and if the sites leave the tape as independents (they are eliminated by the Lagrangian), `w` drops from 17 to 11: 35 → 23 sweeps, another 1.5×. Combined: **245 sweeps → 23 sweeps per derivative-rich evaluation, about 10×**, from 7.3 ms to ~0.7 ms at the same 14 ns/op.

Floor for value + gradient + Hessian + P-parameter Jacobian, nine species:

- **0.7 ms** with the tape interpreted by CppAD as today (23 + 2P sweeps × 30 µs);
- **60–100 µs** if the tape is compiled to machine code (CppADCodeGen or equivalent), where a sweep runs at ~1.5 ns per op;
- textbook "3–5× the value" applies to the gradient only; a dense Hessian of `w` variables is `O(w)` sweeps and cannot be brought below that without sparsity, and this tape is dense.

## 3. Tape cost, and whether caching is the whole story

Measured: recording is 1.2–3.1% of a single evaluation, 1.5% of the reactive solve, 10% of the stability decision (small tapes, cheap sweeps). Caching one tape per topology removes at most that. The premise that caching is the main lever is falsified by these numbers. The main levers, in order:

1. Hessian count per evaluation, `(S + 1)` → 1 (7×), and site elimination (1.5×).
2. Tapes per counted evaluation, ~10 → 1, by evaluating each state once per Ipopt iteration and serving objective, gradient, constraint Jacobian, and Hessian callbacks from it (today each callback and the incipient-vapour closure re-evaluate). Worth ~5–10× on the reactive solve on its own, with no change to the derivative code.
3. Evaluation count in discovery (section 4).
4. Recording (≤10%).
5. Python marshalling on the `Mixture.state` path only (590 µs per call), irrelevant to native solves.

## 4. Evaluation-count floor for a binary stability decision

A tangent-plane stability decision for a binary at fixed T, P in the `(x, v)` space: minimise `TPD(x, v)` from a start set that covers every basin the EOS can have. For a binary there are at most two composition basins per density branch; a defensible start set is the two pure-component corners on each of the liquid and vapour branches plus the feed on the other branch, 4–6 starts. Newton with the exact 2×2 Hessian from the tape converges in 6–12 iterations from those starts. The dual bound loop adds one LP (no EOS work) per major iteration and 1–3 major iterations for a one-phase conclusion. Floor: **40–150 evaluations** (lower end: 4 starts × 8 Newton steps plus 8 for bound verification; upper end: 6 starts × 12 steps, 3 major iterations, and a 20-point coarse sweep to seed the starts).

Measured: 4,631 evaluations, 42 optimizer solves. The excess is spent on (a) the DIRECT-L global search with a work budget of 6,500 (`held2_step2.cpp`, `step2_search_work_budget`), which samples the domain rather than descending from a start set, and (b) the pressure-envelope scan (64 log-volume intervals × 3 points, subdivided to depth 8, up to 120 bisections per bracket) used to find density roots at each trial composition, where a bracketed Newton from the previous root needs 3–6. Both are algorithmic, not implementation, choices. Each of those 4,631 evaluations also records a fresh 301-op tape and takes a dense Hessian it never uses for a root scan.

## 5. Whole-problem floor

Composition: evaluations × per-evaluation cost + KKT linear algebra. A 30-variable KKT system factorises in ~9,000 flops (≈10 µs); over 30 Ipopt iterations that is 0.3 ms and never sets the floor.

**(a) Nine-species reactive bubble point** (finite electrolyte liquid, incipient neutral vapour, two reactions). Unknowns: 9 amounts + 1 packing coordinate for the liquid, 3 for the incipient vapour, P; with multipliers ~25–30. Ipopt iterations: 15–30 with an exact Hessian (today the finite solve reaches its answer in 11–61 iterations per the GREPE plan; the bubble closure adds an outer Newton of 5–18 steps today, which folds into the same system when P is a row). Evaluations: 100–200 liquid evaluations with Hessian plus 50–100 vapour evaluations with gradient.

| endpoint | evaluations | per evaluation | total |
|---|---|---|---|
| lower: compiled tape, one Hessian, 100 + 50 evaluations, 15 Ipopt iterations | 150 | 60–100 µs liquid, 10 µs vapour | **10–15 ms** |
| upper: interpreted CppAD tape cached, one Hessian, 200 + 100 evaluations, 30 iterations | 300 | 0.7 ms liquid, 0.1 ms vapour | **150–220 ms** |
| today | 183–357 counted, 1,700–1,830 tapes, 11,000–11,500 Hessians | 33 ms per counted evaluation | 7.7–11.4 s |

The range is 50–1,000× below today. The upper end needs only the two derivative changes (one Hessian, one evaluation per state per iteration); the lower end needs a compiled tape.

**(b) Neutral binary stability decision.**

| endpoint | evaluations | per evaluation | total |
|---|---|---|---|
| lower: 40 evaluations, compiled tape, Hessian only at Newton steps (w = 4, 9 sweeps at 1.5 ns/op) | 40 | ~5 µs | **0.2 ms** |
| upper: 150 evaluations, interpreted tape (38 µs Hessian + 5 µs value), plus density roots by bracketed Newton at ~4 value-only sweeps | 150 | ~60 µs | **9 ms** |
| today | 4,631 | 145 µs incl. recording and an unused Hessian | 671 ms |

## 6. Assumptions and what falsifies them

| assumption | used | falsified by |
|---|---|---|
| 14 ns per tape op per sweep (interpreted) | measured at n = 9 and n = 2 | a sweep counter that shows more sweeps per Hessian than `1 + 2w` (CppAD may use more for `Hessian(x, l)` with many dependents); measure with `size_forward_bool`/timing per call |
| 1.5 ns per op compiled | literature for CppADCodeGen on this size of tape; not measured here | compile one 2164-op tape and time a sweep; if > 4 ns/op the lower endpoints double |
| `(S + 1)` → 1 Hessian yields the identical tensor | algebra of the implicit function theorem | the #189 equivalence table (elementwise 1e-10) |
| 10 tapes per counted evaluation are redundant re-evaluations of one state | inferred from 1,759 recordings / 183 evaluations | instrument the callback sequence: if the extra tapes are distinct states (line search, vapour closure) the 5–10× from state caching shrinks to 2–3× |
| Newton association in ≤ 5 iterations | standard for a well-scaled mass-action system from a warm start | the iteration counter after the change; today 14–126 |
| 40–150 evaluations suffice for a binary stability decision | algorithm analysis, section 4 | the eight LLE feeds and the neutral test discovering the wrong phase count from the reduced start set |
| Ipopt 15–30 iterations for the bubble point | GREPE plan's measured 7–61 with median 14, plus the pressure row | the iteration counter on the eight cases after the boundary owner lands |
| KKT algebra negligible | 30³/3 flops | a KKT larger than ~200 variables, which no admitted problem has |

## 7. Acceptance criteria, as written into #189 and #188

**#189 (Helmholtz tape and Hessian), quantity: derivative-rich evaluation cost.** Case: the nine-species MEA model (`analyses/2026-mea-parameter-estimation/inputs/neutral-selected-parameters.json`) at 313.15 K, ρ = 3.0e4 mol/m³, the composition in `measure_floor.py`, through the native per-evaluation entry with one active parameter. Threshold: value + gradient + Hessian + parameter Jacobian in ≤ 1.0 ms (interpreted tape; the upper endpoint with 30% margin) and ≤ 23 + 2P sweeps and ≤ 1 Hessian call per evaluation on the counters; the associating binary at 376.88 K in ≤ 5 association iterations per solve. Measurement: `OMP_NUM_THREADS=1 EPCSAFT_WORKTREE=<worktree> python measure_floor.py <out>` (the script in `oracle/`; the counters must remain available as `_core._instrumentation()`), median of 20 calls after 3 warm-ups. Environment record: wheel SHA-256, commit, `ldd` of `_core*.so`, compiler, CPU model, `OMP_NUM_THREADS`, as `measure_floor.py` writes them. A compiled-tape implementation is judged against 100 µs on the same case.

**#188 (phase stability), quantity: evaluations and wall time per stability decision.** Case: the neutral binary fixture (`analyses/2026-fixed-pressure-anchor/inputs/fixtures/neutral_binary.json`) at 300 K, 1e5 Pa, feed (0.5, 0.5), which is one phase. Threshold: ≤ 150 EOS evaluations and ≤ 150 tape recordings per decision (today 4,631), wall ≤ 10 ms on the interpreted tape; the eight water/ethanol/isobutanol/NaCl feeds each in ≤ 600 evaluations. Measurement: `discovery_case()` in `measure_floor.py` and the eight-feed comparison, same environment record. A decision that meets the count but changes the phase count or compositions beyond the Table 3 uncertainties fails.

**Reactive bubble point (composite, #187 + #189):** the four `vle_obs_*` cases in `measure_floor.py` in ≤ 250 ms each (upper endpoint) with ≤ 300 recordings and ≤ 300 Hessian calls per solve.

## 8. Does the floor justify the greenfield rebuild?

On performance alone, no. The measurements say that three bounded changes inside the existing tree reach most of the gain: one Lagrangian Hessian instead of `(S + 1)` (7×), one state evaluation per Ipopt iteration shared across callbacks instead of ~10 (5–10× on the reactive solve, pending the callback instrumentation named above), and Newton-from-starts instead of a 6,500-work global sample plus envelope scan in discovery (30–100× on evaluation count). Tape caching, the change the plan had at the front, is worth 1–10%. Together these take the bubble point from 11 s to the 150–220 ms upper endpoint and the stability decision from 671 ms to under 10 ms without a new tree; the remaining 10× to the lower endpoints needs a compiled tape, which is orthogonal to the rebuild. The rebuild is still justified by what the audit found elsewhere: six problem classes, fifteen definitions of R, the Python numerics, and 46k removable lines. #189 should therefore be reordered so the Hessian change and the per-iteration state cache come first and tape caching last, and its acceptance is the table above, not "one cached tape".

---

# Revision 2 (2026-09-02): every public entry point

Sweep of 39 entry points on `perf/instrumentation` `0db4d46` (adds tape width, range, and density-iteration counters), wheel `324eac4f…`, `OMP_NUM_THREADS=1`, gcc 13.3.0, Ryzen 5 5500. Full table with every counter: `entry-point-profile.md` / `.json`; script `profile_entry_points.py`. Counters are reset before each call; state calls are medians of 5–20, solves single calls after one warm-up.

Counter gap, stated: sweep time is measured only at the two full Helmholtz tapes (`thermodynamic_derivatives`, `fixed_topology_parameterized_mixture_phase_block`); the five small tapes (pure phase, packing, inverse packing, mixture phase) report recordings and Hessian calls but not sweep time, so "sweeps 0%" for saturation, discovery, and the declared routes means unmeasured, not zero. Ipopt iteration counts for saturation and discovery are not exposed. Python marshalling is obtained only by subtraction.

## The four mechanisms, per call

Shares are of wall time; `rec` = tape recordings; `H/rec` = Hessian calls per recording; `assoc/solve` = successive-substitution iterations per association solve; `dens` = density-root iterations.

| entry point | wall ms | rec | rec % | sweep % | H/rec | assoc/solve | dens | dominant mechanism |
|---|---|---|---|---|---|---|---|---|
| `Mixture.state` ρ, neutral | 0.46 | 1 | 3 | 8 | 1 | 0 | 0 | Python marshalling (~85% of wall) |
| `Mixture.state` ρ, associating | 0.99 | 1 | 3 | 53 | 5 | 126 | 0 | (S+1) Hessians; association 126 iterations |
| `Mixture.state` ρ, electrolyte 3 | 0.69 | 1 | 3 | 39 | 3 | 14 | 0 | Hessians |
| `Mixture.state` ρ, nine species | 5.25 | 1 | 1 | 87 | 7 | 44 | 0 | (S+1) Hessians, width 17 |
| `Mixture.state` P liquid, neutral | 7.5 | 158 | 19 | 67 | 1 | 0 | 141 | Hessian per density iteration; 141 iterations; re-taping 19% |
| `Mixture.state` P vapor, neutral | 13.4 | 321 | 21 | 70 | 1 | 0 | 304 | 360-interval scan; re-taping 21% |
| `Mixture.state` P liquid, associating | 42.5 | 80 | 3 | 89 | 5 | 381 | 64 | 5 Hessians × 64 iterations; 30,515 association iterations |
| `Mixture.state` P vapor, associating | 150 | 321 | 4 | 95 | 5 | 27 | 304 | scan × 5 Hessians |
| `Mixture.state` P liquid, nine species | 382 | 77 | 2 | 98 | 7 | 44 | 61 | 7 Hessians of width 17 per density iteration |
| `Mixture.state` P vapor, nine species | 1,691 | 322 | 2 | 98 | 7 | 28 | 306 | scan × 7 Hessians |
| `Mixture.reference`, electrolyte 3 | 108 | 387 | 6 | 92 | 3 | 13 | 307 | four fixed-P states, each a scan |
| `ElectrolyteReference.state` | 47.6 | 158 | 6 | 88 | 3 | 13 | 124 | one fixed-P state |
| `observable_block` fixed-state, neutral | 0.72 | 2 | 4 | 11 | 1 | 0 | 0 | marshalling |
| `observable_block` fixed-pressure-fugacity, neutral | 7.1 | 159 | 18 | 64 | 1 | 0 | 141 | as fixed-P state |
| `Pure.saturation`, propane | 454 | 5,422 | 11 | 15 (partial) | 1 | 0 | 0 | 21 Ipopt starts × ~258 callbacks; ~75% outside tape work |
| `Pure.saturation` with start | 1.3 | 6 | 8 | 8 | 1 | 0 | 0 | already near floor |
| `Pure.saturation`, associating pure | 531 | 1,700 | 4 | 24 (partial) | 3 | 22 | 0 | starts × Hessians |
| declared fixed-T,P one phase | 4.6 | 126 | 24 | n/m | 1 | 0 | 0 | 48-interval pressure scan; re-taping 24% |
| declared fixed-T,P three phase | 10.1 | 45 | 7 | n/m | 1 | 0 | 0 | near floor for its route |
| fixed-T coexistence, associating | 98 | 240 | 5 | 90 | 5 | 212 | 0 | 5 Hessians; 50,840 association iterations |
| fixed-P coexistence, associating | 1,559 (7,166 on the first run; variance unexplained) | 4,082 | 5 | 90 | 5 | 211 | 0 | outer Python bisection × inner solve; 863,025 association iterations |
| fixed-T coexistence, charged 4-comp | 63 | 260 | 8 | 87 | 3 | 23 | 0 | Hessians |
| discovery, neutral binary | 478 | 4,632 | 10 | n/m | 1 | 0 | 0 | evaluation count (DIRECT-L 6,500 budget + envelope) |
| discovery, charged electrolyte 3 | 870 | 3,626 | 8 | 87 | 3 | 16 | 0 | evaluation count × 3 Hessians |
| GREPE finite LLE, ternary | 99 | 1,430 | 17 | 67 | 1 | 0 | 1,028 | 1,028 density iterations inside the start construction |
| GREPE finite + incipient, solved P | 438 | 2,213 | 7 | 78 | 2.3 | 20 | 527 | Hessians; density iterations in vapour closure |
| GREPE homogeneous reaction | 13 | 0 | 0 | 0 | — | — | 0 | Python only (no EOS evaluation) |
| GREPE solved-T bubble, binary | 1,678 | 4,082 | 5 | 90 | 5 | 211 | 0 | same path as fixed-P coexistence |
| GREPE reactive bubble, nine species | 7,416 | 1,759 | 2 | 97 | 6.5 | 46 | 572 | 11,485 Hessians of width 16 |
| `regression.least_squares`, saturation 3 params | 11,705 | 140,562 | 11 | 15 (partial) | 1 | 0 | 0 | 26 cold saturation solves × 5,422 |
| `regression.curve_fit`, same | 880 | 10,848 | 10 | 15 (partial) | 1 | 0 | 0 | 2 cold solves |
| `regression.least_squares`, equilibrium owner | 227 | 795 | 5 | 52 | 1.9 | 15 | 0 | 8 solved-P evaluations |

**Mechanism 1, re-taping per iteration.** Recordings per unit of useful work: 158–322 per fixed-pressure state, 387 per reference, 5,422 per cold saturation, 4,632 per binary stability decision, 4,082 per fixed-P boundary, 140,562 per three-parameter saturation fit. Share of wall: 1–4% on wide tapes (nine species, associating), 6–8% on the electrolyte-3 tape, 18–25% on the small neutral tapes (width 2–4). So tape caching is worth up to a quarter on the neutral fixed-P, declared, and discovery paths and 1–3% on the wide-tape paths. It moves up for the EOS entry points but is nowhere the largest term.

**Mechanism 2, Hessians where a first derivative suffices.** The shared entry is `NativeModel::evaluate` → `thermodynamic_derivatives` (`derivatives.hpp:62`), which always builds the Jacobian and `(1+S)` dense Hessians. `solve_density` (`density.hpp:109`) and `polished_pressure_density` (`eos.cpp:105`) call it once per iteration to read `P` and `dP/dρ`, so every fixed-pressure state pays 61–306 full Hessian sets: 67–98% of its wall. A value-plus-`dP/dρ` path costs one forward-0 sweep, one reverse-1 sweep for the gradient, and one forward-1/reverse-2 pair for `∂²a/∂ρ²`: 4 sweeps against `1 + 2w` per Hessian times `(1+S)`; at width 17 with S = 6 that is 4 against 245, 60×; at width 4 neutral, 4 against 9, 2×. The same entry serves `Pure.saturation` through `pure_phase_block` with an extra third-order sweep per callback.

**Mechanism 3, association convergence.** Iterations per solve today: 126 (associating state at fixed ρ), 381 (fixed-P liquid associating), 212 (coexistence), 44–46 (nine species), 13–27 (electrolyte 3), 22 (associating pure saturation). Distribution: the fixed-P associating liquid ran 30,515 iterations over 80 solves; the fixed-P boundary 863,025 over 4,082. Newton on the site fractions with `H_X` converges in 3–6 from the previous state (the residual is 2e-13). Cost share: at S = 4 an iteration is ~150 flops, so 863k iterations are ~0.1–0.2 s of the boundary's 1.6–7 s (5–10%); on `Mixture.state` at fixed ρ the 126 iterations are ~20 µs of 0.99 ms. This mechanism is real, isolated, and second-order in wall time; it is being fixed separately and is assumed landed below.

**Mechanism 4, evaluation count per solve.** Saturation: 5,422 recordings for a three-variable problem, from 21 Ipopt starts (7 vapour × 3 liquid seeds, `pure_saturation.cpp:317-330`) each running ~258 callbacks (max_iter 50, no callback caching); with a start it takes 6. Floor without a start: 2–3 seeds from a Clausius–Clapeyron estimate × 6–8 Newton steps = 15–25. Discovery: 4,632 against 40–150 (section 4 above). Fixed-P boundary: 4,082 recordings from a Python temperature bisection (64 iterations cap) around an inner fixed-T solve; with T as an unknown row, one Newton of 10–15 steps over two phases: ~30 evaluations. Fixed-pressure states: 141–306 density iterations from the 360-interval scan; with an anchor or a bracket from the previous root, 5–8 (the `density_anchor` path already does this: `Pure.saturation` with a start is 6 recordings). GREPE finite LLE: 1,028 density iterations in start construction against ~40. Regression on saturation: 140,562 recordings for 26 residual evaluations because each is a cold 21-start solve; with continuation, 26 × 6.

## Per-entry-point floors

Lower endpoint: compiled tape (1.5 ns/op/sweep), one Hessian, Newton counts; upper endpoint: interpreted CppAD tape cached (14 ns/op/sweep), one Hessian, Newton counts; both after association Newton. Marshalling floor for a `py::class_` binding is taken as 20–50 µs per call.

| entry point | today | floor (compiled … interpreted) | what sets the endpoints |
|---|---|---|---|
| `Mixture.state` ρ, neutral | 0.46 ms | 30–110 µs | marshalling 20–50 µs + 9 sweeps × 301 ops |
| `Mixture.state` ρ, nine species | 5.25 ms | 0.1–0.7 ms | one Hessian at width 11 (23 sweeps × 2164 ops) |
| `Mixture.state` P, neutral, with anchor | 7.5 ms | 0.1–0.2 ms | 6 iterations × 4 sweeps |
| `Mixture.state` P, neutral, no anchor | 7.5–13 ms | 0.3–0.5 ms | 20–40 value-only scan points + Newton |
| `Mixture.state` P, nine species | 0.38–1.7 s | 1–5 ms | 6–40 iterations × (value + 3 sweeps) |
| `Mixture.reference` | 108 ms | 2–5 ms | four fixed-P states with anchors |
| `Pure.saturation`, cold | 454 ms | 0.3–1 ms | 15–25 evaluations, 3-variable Newton, no Ipopt |
| `Pure.saturation`, associating pure | 531 ms | 1–3 ms | as above with one Hessian |
| declared fixed-T,P three phase | 10 ms | 2–5 ms | 20–40 evaluations at width 4 |
| fixed-T boundary, associating | 98 ms | 2–5 ms | 10–15 Newton steps × 2 phases, one Hessian |
| fixed-P boundary | 1.6–7.2 s | 3–8 ms | T as an unknown row; no outer bisection |
| discovery, neutral binary | 478 ms | 0.2–9 ms | 40–150 evaluations |
| discovery, charged 3 | 870 ms | 15–30 ms | 150 evaluations × one Hessian at width 6 |
| GREPE finite LLE, ternary | 99 ms | 3–6 ms | 30–60 evaluations × 0.1 ms; anchored density |
| GREPE solved P, 4-comp charged | 438 ms | 20–50 ms | 150 evaluations × 0.3 ms |
| GREPE reactive bubble, nine species | 7.4 s | 10–15 ms … 150–220 ms | section 5 |
| `regression.least_squares`, saturation | 11.7 s | 30–50 ms | 26 evaluations × warm saturation |
| `regression.least_squares`, equilibrium owner | 227 ms | 10–20 ms | 8 evaluations × 2 ms |

## Revised ordering for #189 (per family, not one list)

- **EOS fixed-pressure states, reference, and every route that solves density** (`Mixture.state(P)`, `observable_block` fixed-pressure, `Mixture.reference`, GREPE start construction): (1) a derivative entry with an order argument so density iterations get value and `dP/dρ` from 4 sweeps, not `(1+S)` Hessians (2–60×); (2) Newton from an anchor or the previous root instead of the 360-interval scan (20–50× on iterations); (3) tape caching (10–25% on width ≤ 4, 1–3% otherwise); (4) one Hessian where a Hessian is wanted.
- **Saturation and saturation fits**: (1) warm continuation seeds and a 3-variable Newton in place of 21 Ipopt starts (5,422 → 15–25 evaluations; the fit from 140,562 to ~160); (2) drop the third-order sweep from callbacks that do not need it; (3) tape caching (11%).
- **Associating and electrolyte states, coexistence, discovery (charged)**: (1) one Lagrangian Hessian instead of `(1+S)` (3–7×); (2) evaluation count (discovery 30–100×; boundaries 100×); (3) tape caching (4–10%).
- **Reactive bubble point**: unchanged from revision 1: Lagrangian Hessian, one state evaluation per Ipopt iteration shared across callbacks, tape caching last.

## Acceptance criteria (revised; posted on #189)

Held to on `profile_entry_points.py` with the environment record it writes, medians of the stated repetitions, after association Newton lands:

| case | quantity | threshold |
|---|---|---|
| `Mixture.state` P liquid, nine species, with a density anchor | recordings ≤ 8, Hessian calls ≤ 1, wall ≤ 5 ms | today 77 / 539 / 382 ms |
| `Mixture.state` P vapor, neutral, no anchor | recordings ≤ 40, wall ≤ 0.5 ms | today 321 / 13 ms |
| `Mixture.reference`, electrolyte 3 | recordings ≤ 40, wall ≤ 5 ms | today 387 / 108 ms |
| `Pure.saturation`, propane 300 K, cold | recordings ≤ 30, wall ≤ 1 ms | today 5,422 / 454 ms |
| `regression.least_squares`, saturation, 3 parameters | recordings ≤ 300, wall ≤ 50 ms | today 140,562 / 11.7 s |
| `Mixture.state` ρ, nine species | Hessian calls = 1, sweeps ≤ 23 + 2P, wall ≤ 1 ms | today 7 / 245 / 5.3 ms |
| GREPE reactive bubble, nine species | ≤ 300 recordings, ≤ 300 Hessian calls, ≤ 250 ms | today 1,759 / 11,485 / 7.4 s |
| discovery, neutral binary | ≤ 150 recordings, ≤ 10 ms (#188) | today 4,632 / 478 ms |

## Verdict, from the complete measurement

Unchanged in conclusion, sharpened in detail. Every entry point is dominated by one of four bounded defects inside the current tree: a derivative entry that always builds `(1+S)` Hessians (fixed-pressure states, reference, GREPE start construction, saturation callbacks), a 21-start Ipopt saturation with no continuation (saturation and every fit on it), evaluation-count policies in discovery and the fixed-P boundary, and `(1+S)` Hessians per evaluation in the reactive solve. Tape caching is the fifth lever and worth 1–25% depending on tape width; it is not the plan's first item on any path. None of the four requires a new tree; the rebuild stands on ownership and size, and #189's acceptance is now the per-entry table above rather than one number.

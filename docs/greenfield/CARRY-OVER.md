# Carry-over manifest

What is copied from `main @9eccbfd`, what is written fresh, and what has no successor.

The default is **written fresh**. A file carries only if it is (a) correct, (b) covered by
a test or an oracle row, and (c) free of the structures the rebuild exists to remove. Three
files meet all three; a handful more carry as algorithms rather than as text.

A carried file arrives with its `// EqID:` annotations intact — 57 of them across the
physics headers — because the equation registry (`equations_registry.yaml`,
`sync_equation_registry.py`) is a working system that the rebuild has no quarrel with.

---

## 1. `physics/` — mostly carried

| target | source | lines | how |
|---|---|---|---|
| `physical_constants.hpp` | `_native/physical_constants.hpp` | 14 → 30 | carried, extended to hold every category-1 constant. Today 15 definitions of `R` live across the tree; the checker forbids the 16th. |
| `helmholtz.hpp` | `_native/helmholtz.hpp` | 357 → 200 | carried. Hard chain, dispersion, packing-volume closed form, 26 `EqID`s. Trimmed where it reaches for solver state it should not have. |
| `association.hpp` | `_native/association.hpp` **on `perf/association-newton` @a248de6**, not `main` | 265 → 150 | carried from the branch. The hybrid Newton/successive-substitution solver with the analytic Jacobian `∂F_a/∂X_b = X_a ρ_N x_{c(b)} n_b Δ_ab + δ_ab (1 + ρ_N B_a)`. Measured: 264→7, 486→8, 31,056→645 iterations, agreement 1.5e-14. `main` still has the damped-substitution version at 1,000 iterations — do not copy that one. |
| `electrolyte.hpp` | `_native/electrolyte.hpp` | 399 → 260 | carried with one merge: five permittivity formulations become four, because Figiel is Zuber with `α_i = C` for every ion, and solvent-only is ion suppression with `α = 0`. The four that remain are distinct published equations — ion suppression `ε_sf/(1 + Σ α_i x_i)`, linear mixing (Ascani), density–temperature (Schreckenberg), dipolar perturbation (Rueben) — selected by the coefficient set in the parameter record rather than by a formulation enum. |
| `born.hpp` | `_native/born.hpp` | 87 → 60 | carried nearly verbatim. |
| `numerical_policy.hpp` | — | 0 → 250 | **written**. Every tolerance as a `Tolerance` entry with its category. See TRANSFER.md Stage 4 for the HELD2 drift resolved here. |

## 2. `model/` — algorithms carry, code does not

| target | relationship to today | lines |
|---|---|---|
| `helmholtz_tape.hpp` | **written**. `derivatives.hpp` (238) is the `(1+S)`-Hessian site the rebuild exists to remove; it is not a starting point. The tape is designed around the Lagrangian Hessian and site elimination from the first line. | 350 |
| `helmholtz_model.{hpp,cpp}` | **written**. Replaces `eos.cpp` (2,033) + `eos.hpp` (514) + `model.hpp` (59) + `eos_access.cpp` (598). Parameter validation happens once, in the constructor; the three passes in `_dictionary_io.py`, `_resolve.py`, and the native side become one. | 750 |
| `density.cpp` | **algorithm carries** from `density.hpp` (255): the bracketed scan plus Newton on the pressure row is right. What does not carry: 360 intervals (→ anchor-first, 5–8 iterations), and calling `thermodynamic_derivatives` per iteration for a full `(1+S)` Hessian set to read `P` and `dP/dρ` — a value-plus-`dP/dρ` path is 4 sweeps against 245 at nine species. | 250 |
| `infinite_dilution_reference.cpp` | **written**, from `eos_source_reference.cpp` (904) + `electrolyte_reference.hpp` (221, 5 `EqID`s). The dilution sequence and MIAC/osmotic definitions carry as equations; the four-fixed-P-states-each-a-scan structure (108 ms per call) does not. | 450 |

`phase_block.hpp` (1,499, 20 `EqID`s) has no successor as a file. Its equations are
absorbed into the tape and the minimizer; its `py::dict` marshalling is not written.

## 3. `equilibrium/` — the formulation carries, no file does

Nothing in `engine/cpp/equilibrium/` (34,095 lines, 30 files) is copied.

What carries is the objective, already written down and registered:

    G = Σ_p [A_p(T, V_p, n_p) + P V_p] + Σ g° n      over (ln n_p, ln η_p)

This is `reactive_finite_phase.cpp`'s objective (3,723 lines) and Pereira 2012's, and it
is documented with 15 registered equations in
`engine/docs/equilibrium/formulation.md` (PR #194, branch `docs/equilibrium-formulation`)
— **merge that branch before Stage 5 begins**; it is the specification the stage is
written against.

| what carries as an algorithm | from | note |
|---|---|---|
| the LP master problem (problem OP) | `held2_step4.cpp`, 60 lines of model assembly | follows the paper exactly, feed cut included. Its 90-line certificate does not carry: HiGHS's `kOptimal` with feasible primal and dual status is the same statement. |
| the electroneutral elimination | `held2_step1.cpp` and `chemical_equilibrium.cpp` — **two implementations of one idea** | one chart is written. The 850 lines become 60. |
| balance/reaction rank by SVD | `phase_direction_compiler.cpp` (430), and again in `_grepe.py:616`, and again in `_phase_reactions.py:134` | one native implementation; the two Python versions are not written. |
| the Gibbs minimization itself | six `Ipopt::TNLP` subclasses across `grepe.cpp` (4,487), `reactive_finite_phase.cpp` (3,723), `phase_equilibrium.cpp` (1,051), `coexistence.cpp` (1,815), `held2_step5.cpp` (1,984), `held2_step8_nlp.cpp` (2,137) | one `GibbsMinimizer`. All six solve the same problem with different row sets; the row set is now data. |
| the dual bound loop and TPD search | HELD2 stages I–II, ~6,480 lines | written to Pereira 2012 per `PHASE-STABILITY.md`. HELD v1 (`held_stage_{i,ii,iii}.cpp`, 3,051) has no successor at all — it is the superseded implementation of the same algorithm. |

## 4. Python — three modules survive in substance

| target | source | lines |
|---|---|---|
| `units.py` | `units.py` | 86 → 80, carried |
| `__init__.py` | `__init__.py` | 85 → 60, rewritten to 11 exports |
| `equilibrium/loader.py` | `equilibrium/_feed.py` (364) + parts of `_grepe_io.py` (956) | → 120. The mapping from a document to a `Problem` is worth keeping; the mass-fraction and molality *arithmetic* moves native into `request_compile.cpp` (Python does no arithmetic), and the request-family branching is not written. |

No successor: `_grepe.py` (2,770), `_api.py` (2,423), `eos.py` (2,412), `_dictionary_io.py`
(1,839), `records.py` (1,171), `_resolve.py` (1,165), `_grepe_io.py` (the rest),
`_phase_reactions.py` (328), `regression/_api.py` (292), `_canonical.py` (58),
`_trials.py` (53). The new `eos.py` (600), `parameters.py` (400), `reference.py` (100),
`problem.py` (450), and `regression.py` (200) are written against the bindings, not
derived from these.

The Python `_fingerprint` (25 sites) and the C++ canonical-JSON hash become one native
`parameter_hash`. `sha256.hpp` carries.

## 5. Tests — the oracle carries, the suite does not

`engine/tests/` is 4,404 lines across a `numerical/` / `behavior/` split that the rebuild
drops: tests are organised by capability, because a test is written against a user-facing
claim.

Carried, as the acceptance contract rather than as test code:

| oracle case | rows | what it pins | stage |
|---|---|---|---|
| `aqueous_alkali_halide_activity` | 24 | LiCl/LiBr/KCl/KBr mean ionic activity coefficients; the infinite-dilution log fugacities and convergence errors were computed columns the oracle rebuild removed | 4 |
| `light_alkane_saturation` | 10 | ethane saturation pressures and saturated-liquid mass densities against NIST SRD 69 — the propane and n-butane rows, and the vapour molar volumes, were computed columns the oracle rebuild removed | 7, 10 |
| `reactive_vle_baseline.json` | 4 | nine-species MEA–CO₂–H₂O reactive bubble point: liquid molar densities 46,425–47,148 mol/m³, and today's 7.4–8.2 s at 225–252 evaluations. Reported as a change indicator, never pinned — they are this engine's own output | 12 |
| `floor-measurements.json`, `entry-point-profile.json` | 39 entry points | every counter and wall time on `main`, with the environment record | every stage |

`reactive_equilibrium_reachability` (2 rows) is not a target — it is a record that both
downstream repositories cannot construct a model against `main` at all. It is the
evidence for Stage 8, not a gate.

One synthetic fixture carries: the associating binary, which is Baygi & Pahlavanzadeh
2015's MEA–water model and already lives in `tests/native/model_support.hpp`. The
immiscible ternary does **not** carry (revised 2026-09-02): its two-liquid cell was never
verified to exist, and `EQUILIBRIUM.md` §4.2 replaces it with a published
fixture, the hypothetical `A + B ⇌ C` of Ascani et al. 2023 (Table 3, `K_a = 2.25`),
whose miscibility gap and reactive tie lines are documented in that paper. Everything
else is written against the oracle rows, on the user's rule: real mixtures and real data,
no route-specific or downstream-specific gates or flags.

## 6. Outside the engine — untouched

`data/packets/`, `analyses/`, `docs/scientific/adr/`, `docs/latex/`,
`equations_registry.yaml`, `scripts/docs/sync_bibliography.py`, `scripts/docs/
sync_equation_registry.py`, and the Zotero export path are not part of the transfer.
`analyses/` is migrated at Stage 8 when the field names change.

## 7. Two branches to land before their stage

- **`perf/association-newton` @a248de6** — the source for `physics/association.hpp`.
  The solver was uncommitted working-tree state in a scratch worktree until
  2026-09-02; `f0abace`, the branch tip before that, carries only the
  instrumentation counters.
  Either merge it to `main` first or copy from the branch; do not carry `main`'s
  damped-substitution solver.
- **`docs/equilibrium-formulation` (PR #194)** — `formulation.md`, 15 registered
  equations, the written specification for Stage 5. Merge before Stage 5.

The saturation-continuation change (`eos_regression_evaluator.cpp`, fit 3,509 ms → 330 ms
with identical fitted parameters `[2.001999999999, 3.618400000002]`) is on the same perf
work; it does not carry as code — Stage 5 writes continuation into `phase_boundary.cpp` —
but it carries as the evidence that the 21-start grid is unnecessary.

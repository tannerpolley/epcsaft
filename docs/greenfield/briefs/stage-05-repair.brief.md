# Stage 05 repair — eight audit findings on the merged equilibrium

Branch `stage-05-repair`, based on the greenfield commit the coordinator pins. Machine-readable
claims are in `briefs/stage-05-repair.brief.json`; that file governs where the two disagree.
`EQUILIBRIUM.md` stays the specification; this brief names seven defects the independent
audit (b74ae3b → f6efa0b) found in the accepted tree, each verified by reading, and the
check that turns red on each. Nothing else moves. Stage 6 depends on the temperature and
feed sensitivities, so this repair merges before Stage 6 is accepted.

## The seven, and what is true

**1. The feed sensitivity omits the pinned amount's dependence on the feed total.** A pinned
phase `q` has `N_q = β_q · N_feed` with `N_feed = Σ_j n_feed,j`. Balance row `r` is
`(Σ_{p free} N_p (B x_p)_r + β_q N_feed (B x_q)_r − (B n_feed)_r) / N_feed`. Its derivative in
`n_feed,i` at a converged point is

    ∂row_r/∂n_feed,i = −B[r,i]/N_feed + Σ_{q pinned} β_q (B x_q)_r / N_feed

The second term is the same for every feed column, vanishes when `β_q = 0` (which is why
attempt 4's S9 passed on boundary cells), and is what the audit measured missing at
`β = 0.2` (analytic −1,361,623 against −1,748,683 Pa/mol by central difference, 22 % off).
`lever_jacobian` gains it on the balance rows; S9 gains the fixed-fraction flash at
`β = 10⁻²` with the feed column.

**2. The temperature column of the chemical rows omits the standard term's derivative, and
the standard term is frozen at the guess temperature.** The chemical row carries
`s_i = g°_i / (RT)`, and compile evaluates `standard_gibbs` once at `temperature_guess`
(`request_compile.cpp:161`). Two consequences on a free-`T` cell: the term's own
`∂s_i/∂T` (for a constant `g°_i` in J/mol it is `−g°_i/(RT²)`; for a constant `K` written as
`s = −ln K` it is zero; for a `K(T)` correlation it is what the Properties stage's
`standard_gibbs(T, order)` will return) is absent from the `T` column of `lever_jacobian`,
of the free-`T` Jacobian and of the objective gradient; and the value itself is stale away
from the guess. The repair: the standard term is a function `standard_term(T)` of the
declared source returning value and `∂/∂T` exactly, evaluated at the iterate's `T`
everywhere the row is assembled, and its derivative enters every `T` column. S9 gains a
reactive free-`T` cell with a nonzero `standard_gibbs_j_per_mol` offset (the audit's
±1,000 J/mol at 300 K measured a ±0.0013 mismatch).

**3. Michelsen's deficit divides by `RT` twice.** The rows are dimensionless (§2.2), and
`gibbs_minimizer.cpp:117` forms `exp(−r/RT)`. The deficit is `1 − Σ_i x_i exp(−r_i)`. S15
asserts the reported deficit equals that expression recomputed from the returned residual
rows, so the diagnostic and the number the failure carries are one.

**4. The equilibrium record never reaches the trace.** The driver builds
`result.convergence` and never calls `scope.finish(convergence)`; `density.cpp` does. Rule
9's promise, a solver's span carrying its record, is broken for `solve_equilibrium`. S1
asserts, under `core::Capture`, that the `solve_equilibrium` span's record equals
`result.convergence`.

**5. `compile_problem` returns a `CompiledSystem` holding a raw model pointer with no
`keep_alive`**, so Python can collect the model under a live compiled object. The binding
gains `py::keep_alive<0, 1>()`; a Python test asserts the model's reference count rises by
one while the compiled object lives and falls when it is deleted.

**6. The surface drops unknown species keys and unknown amount statuses silently.**
`_vector` reads `values.get(component_id, 0.0)`, so a misspelt species key is ignored
(a missing key is legitimately zero: a pure feed in a binary); `_phase_from_mapping` turns
any status other than `"pinned"` or `"free"` into free. Both reject by name:
`unknown species <name>`, `unknown amount status <status>`.

**7. Two physical defaults on the surface, unseen by the checker.**
`Reaction(equilibrium_constant=1.0)` and `solvent_mass_kg=1.0` (two sites) are rule-12
defaults; the checker scans only `src/epcsaft/*.py` (`surface.glob`), not the `equilibrium/`
subdirectory. The checker moves to `rglob`, and the defaults become required arguments.
The checker change lands with the fix: the architecture gate is red on this branch between
the two writes and green after both, and a mutant that restores one default is red.

Canon, one sentence for the coordinator, not a worker write: S9's row says `−1/RT` on the
pressure row; the assembly has `−v_p/(RT)`, written as `−c_p/(η_p RT)` with `c_p` the packing
coefficient, which is the derivative of `Z_p − P v_p/(RT)` in `P`.

**8. The assembly has no active-parameter columns in `F_θ`.** `EQUILIBRIUM.md` §12 promised
`dy/dθ` for the model's active parameters from `∂²A/∂y∂θ`, and attempt 4's `lever_jacobian`
carries `T`, `P` and the feed only, so a regression over the six EOS families has no analytic
sensitivity. The tape already carries the block: `HelmholtzDerivatives.gradient` and `hessian`
are laid out over `(T, v, x_1..x_C, θ_1..θ_A)` with `A = active_.size()` (`helmholtz_model.cpp`,
`total_coordinates = 2 + components + active_.size()`). `lever_jacobian` becomes
`rows × (2 + C + A)` with the active columns after the feed columns, in the model's active order;
`sensitivity(compiled, point, lever)` takes a weight vector of that length. Per row, column
`2 + C + k`: a chemical row of species `i` in phase `p` is `hessian[(2 + i) · width + (2 + C + k)]`,
the θ-derivative of the residual chemical potential over `RT` at fixed `(T, v, x)`, and nothing
else, because the log-coordinate transforms and the standard term carry no θ (the ions'
standard-state transfer's θ-dependence belongs to the Properties stage with
`standard_gibbs(T, order)`); the pressure row `Z_p − P v_p/(RT)` is
`−v_p · hessian[1 · width + (2 + C + k)]`, the θ-derivative of `∂A/∂v` scaled as the row is;
balance, normalisation and charge rows are zero. The `k = 0` path is the same block: the
stationarity rows are the chemical rows scaled by `N_p`, so the driver's KKT `F_θ` reuses it.
Token `active_parameter_column`. S9 gains: the neutral binary `{V,L}` flash at 180 K, fixed `T`
and `P`, with the methane–ethane `k_ij` and methane's segment number active; `dy/dθ_k` from
`sensitivity` with the unit weight on column `2 + C + k` against a centred difference of two full
solves on models rebuilt at `θ_k ± δ` (`HelmholtzModel(parameters_perturbed, active)`),
`δ = 10⁻⁴` relative, 1e-6 relative; and the exact zero on a balance row's θ-column. Mutation:
the active column multiplied by zero.

**Amendment after the first stop (coordinator, 2026-09-03).** Rule 12 is `SURFACE.md`'s: every
function-signature default under `src/epcsaft/` is the constant `None`. Widening the checker to the
subpackage therefore governs every default in `equilibrium/problem.py`, not only the two physical
ones: the physical ones become required arguments; structural empties become `None` with an
explicit replacement at the top of the function. That is within this repair's writes and ceilings
now. A replacement may be an empty structure, never a number: `support`, `composition_guess`,
`reactions`, `standard_gibbs_j_per_mol` and the molality solvent list become empty tuples or lists;
`Free(guess=None)` and `packing_guess=None` cross as *absent*, so the start rule seeds them and the
trace's `start_rule` event says so (a zero sentinel that the native side turns into 300 K or 1e5 Pa
is a physical default in disguise). If `FreeValue` needs an explicit absent state, that is a
two-line native change under the same finding (`CORE.md` §4.3), not scope creep. And `sensitivity()` lives in `implicit_sensitivity.cpp`, which is now a listed write.

**S9-T-standard, replacing the free-`T` solve (planner's ruling after the fourth stop).** Compile the
reactive fixture twice with `standard_gibbs_j_per_mol = {1000, −1000, 0}`: (i) `L₂` pinned at
`β = 10⁻²`, `T` free, `P` fixed at 1 bar; at the start the start rule produces (`boundary_start`),
the `T` column of `compiled.jacobian` against a centred difference of `evaluate` in the `T`
coordinate, step 1e-3 K, 1e-6 relative; (ii) the same problem with `T` fixed at 298.15 K; at its
start, the `T` lever column of `lever_jacobian` against a centred difference of `evaluate` over the
system recompiled at `T ± 1e-3` K, 1e-6 relative. No solve in either. Mutation unchanged: the
standard term's derivative multiplied by zero reddens (i) and (ii). The comparison at any interior
point falsifies exactly finding 2's defect; whether the free-`T` two-liquid solve converges is a
separate finding, recorded in `STAGE-LOG.md`, not this repair's.

## What you inherit and may not change

Everything else under `equilibrium/`: the rows, the driver, the start rule, the coincidence
rule, the tolerances, `classify`. No route name, no second sensitivity, no third derivative.
Read the tree through Serena before writing: `get_symbols_overview` on
`equilibrium/request_compile.cpp` and `equilibrium/gibbs_minimizer.cpp`; `find_symbol` for
`lever_jacobian`, `standard_state`, `michelsen_deficit`, `solve_equilibrium`;
`find_referencing_symbols` for `standard_gibbs` and `lever_jacobian`. Record what it found.

## Budgets

Ceilings only (`CANON.md` "Budget"). `request_compile.cpp` sits at 848; its ceiling moves
to 960 for finding 2's function. A breach is `budget_ceiling`, measured as first written.

## Checks, named

| id | asserts | mutation (must go red) |
|---|---|---|
| S9-feed-pinned | on the neutral-ternary fixed-fraction flash at `β = 10⁻²`, `dy/dn_feed,i` from `sensitivity` against a centred difference of two solves, every `i`, 1e-6 relative | the pinned-amount term multiplied by zero |
| S9-T-standard | on the reactive fixture with `standard_gibbs_j_per_mol` offsets of ±1,000 J/mol and `T` free at fixed `P`, `dy/dT` against a centred difference, 1e-6 relative | the standard term's `∂/∂T` multiplied by zero |
| S15-deficit | the failure's deficit equals `1 − Σ_i x_i exp(−r_i)` recomputed from the returned rows, 1e-12 | `exp(−r/2478)` in place of `exp(−r)` |
| S1-record-in-trace | the `solve_equilibrium` span's record equals `result.convergence` | `scope.finish` removed |
| keep-alive | `sys.getrefcount(model)` rises by one across `compile_problem` and falls after `del` | `keep_alive` removed |
| reject-by-name | a misspelt species key and an amount status `"fixed"` each raise with the offending name | the message loses its name |
| no-surface-defaults | the architecture checker is green with `rglob` and the defaults removed | one default restored |

## Stopping

`budget_ceiling`, `gate_impossible`, `source_contradicted`, `environment_broken`,
`mutation_invalid`. Do not accept your own work; the coordinator rebuilds clean, runs
every gate and every mutation, then reads.

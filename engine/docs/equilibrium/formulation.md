# Equilibrium formulation

Status: scientific equation authority for phase, chemical, and coupled reactive-phase equilibrium. Equation identifiers below are registered in `engine/docs/equations_registry.yaml`; runtime source is the executable authority and this document creates no second implementation.

## 1. One problem

Phase equilibrium, chemical equilibrium, and coupled reactive-phase equilibrium are one constrained minimization. Equality of chemical potential across phases and vanishing reaction affinity are both first-order stationarity conditions of that minimization under different constraint rows. Nonreactive equilibrium is the case with no reaction rows; single-phase chemical equilibrium is the case with one phase; the coupled problem has both constraint sets active. Nothing below is derived separately for the three cases.

Sources: \cite{Pereira2012} (HELD: Helmholtz free Energy Lagrangian Dual) for the Helmholtz-coordinate objective, the tangent-plane criterion in volume–composition space, and the Lagrangian dual bound loop; \cite{Michelsen1982} for the tangent-plane stability criterion; \cite{MitsosBarton2007} for the dual extremum principle; \cite{Perdomo2025} for electroneutral reduced coordinates; \cite{Gross2001} for pressure and chemical potential from the residual Helmholtz energy.

## 2. Free energy and coordinates

The EOS supplies the residual Helmholtz energy `a^r(T, ρ, x)` and its exact derivatives. For a phase `p` with amounts `n_p` (mol) and volume `V_p` (m³) at temperature `T`:

```text
A_p(T, V_p, n_p) / (R T) = Σ_i n_{p,i} [ ln(n_{p,i} / V_p / ρ_ref) − 1 ] + N_p a^r(T, N_p / V_p, n_p / N_p),   N_p = Σ_i n_{p,i}
```

`EqID: phase_helmholtz_amount_volume`. Source: \cite{Gross2001} Appendix A with the ideal-gas mixing term; `ρ_ref = 1 mol m⁻³` contributes a species-linear constant only.

The system objective at fixed `T` and `P` is the Gibbs energy expressed in Helmholtz coordinates:

```text
G(T, P, {n_p, V_p}) / (R T) = Σ_p [ A_p(T, V_p, n_p) / (R T) + P V_p / (R T) ] + Σ_p Σ_i n_{p,i} g°_i / (R T)
```

`EqID: system_gibbs_helmholtz_coordinates`. Source: \cite{Pereira2012} Eq. 2–4; exact derivation from `G = min_V [A + P V]`. `g°_i` is the standard-state chemical potential offset supplied by a reaction record (zero when no reaction record references species `i`).

**Helmholtz versus Gibbs.** `G(T, P, n) = min_V [A(T, V, n) + P V]` for a mechanically stable phase. Minimizing `A + PV` jointly over `(n_p, V_p)` and minimizing `G` over `n_p` with `V_p` eliminated by the density solve are the same problem; the joint form makes the pressure-equality condition `∂A_p/∂V_p = −P` a stationarity row instead of an inner root solve and needs no density branch selection inside the optimizer. This is the formulation of \cite{Pereira2012} and of the current `reactive_finite_phase.cpp` objective (`evaluate_phase`: `Φ + P V / (R T) + Σ g° n`). Every part of the target uses it: the equilibrium minimizer, the phase-stability inner problem, the phase boundary, and pure saturation. "Gibbs minimizer" in the architecture text means this object. The dual in \cite{Pereira2012} is taken over the material-balance multipliers of the same objective (section 6); its inner problem is the single-phase form of the same objective with a linear term, so the discovery kernel reuses the minimizer with fixed multipliers. The reuse is sound because the inner problem differs from the one-phase minimization only by the term `−λ · n` and by dropping the balance rows.

Coordinates used by the solver are logarithmic amounts `ln n_{p,i}` on the phase support and the logarithmic packing fraction `ln η_p` (`EqID: phase_packing_volume`: `V_p = (π/6) Σ_i n_{p,i} m_i d_i³ / η_p`, exact derivation from the packing-fraction definition in `full-eos.md` §4). Bounds `0 < η_p < η_max` keep every trial inside the EOS domain; the logarithm keeps amounts positive without inequality constraints. For a charged phase the amount coordinates are the electroneutral chart of \cite{Pereira2012}-style modified fractions (`EqID: electroneutral_amount_chart`; \cite{Perdomo2025} §3, C−2 independent coordinates per phase); a neutral phase uses the ordinary chart.

## 3. Variables

| symbol | meaning | present when |
|---|---|---|
| `n_{p,i}` | amount of species `i` in phase `p`, `i ∈ S_p` (support) | always |
| `V_p` (as `ln η_p`) | phase volume | every EOS phase; absent for an ideal-solution phase |
| `T`, `P` | temperature, pressure | one may be unknown |
| `x_inc` | composition of an incipient phase (zero amount) | boundary problems |
| `Π` | phase count | unknown when topology is omitted |

Reaction extents are not separate variables: reactions enter through the balance matrix `B` (conserved moieties), whose nullspace is the permitted reaction span, and through `g°`. Extents are recoverable as `ξ = ν⁺ (n − n_feed)`.

## 4. Constraints

| id | constraint | rows |
|---|---|---|
| `EqID: material_balance` | `Σ_p B n_p = B n_feed` (`B = I` when no reactions) | rank `B` |
| `EqID: phase_electroneutrality` | `z · n_p = 0` for every phase with charged support | one per such phase |
| `EqID: phase_support` | `n_{p,i} = 0` for `i ∉ S_p` (structural, removes the variable) | — |
| `EqID: incipient_normalization` | `Σ_i x_{inc,i} = 1` | one per incipient phase |
| `EqID: pressure_equality` | `−∂A_p/∂V_p = P` | one per EOS phase (a stationarity row, see §5) |

Well-posedness of the declared rows is proved before the solve by `compile_phase_direction_system` (`engine/cpp/equilibrium/phase_direction_compiler.cpp`): singular-value rank of the declared direction rows against `rank_tolerance · max(1, ‖D‖)`, conservation of every balance row and per-phase charge along each direction, and the orthonormal invariant complement. Redundant reactions, non-conserving directions, or a non-neutral reference phase fail with the named deficiency. No second rank routine exists in the target; the two Python `_matrix_rank` functions are deleted.

## 5. Stationarity and KKT

With multipliers `λ` on `material_balance` and `φ_p` on `phase_electroneutrality`, the Lagrangian stationarity conditions are

```text
∂G/∂n_{p,i} = μ_{p,i} = λ · B_i + z_i φ_p        for i ∈ S_p          (EqID: chemical_potential_stationarity)
∂G/∂V_p     = −P_p(T, V_p, n_p) + P = 0                                (EqID: pressure_equality)
```

Sources: \cite{Pereira2012} Eq. 5–8; \cite{Gross2001} Eq. A.33 for `μ_i^r`. Consequences that are not imposed a second time:

- **Phase equilibrium.** For a species in two supports, `μ_{p,i} − μ_{q,i} = z_i (φ_p − φ_q)`: equal chemical potential for neutral species, equal electrochemical potential `μ̃_i = μ_i + z_i F Δφ` for ions (`EqID: interphase_potential_equality`).
- **Chemical equilibrium.** For any reaction `ν` in the nullspace of `B`, `Σ_i ν_i μ_{p,i} = λ · (B ν) + φ_p (z · ν) = 0` because `B ν = 0` and `z · ν = 0` (`EqID: reaction_affinity`). Affinity vanishes as a consequence of the balance multipliers; a problem with `B = I` has no reaction and no affinity row.
- **Coupled.** Both hold with the same `λ`, `φ`.

Second-order local condition: the reduced Hessian `Z^T ∇²L Z` over the constraint nullspace `Z` is positive definite (`EqID: reduced_hessian_stability`; today `reduced_hessian_minimum` in `reactive_finite_phase.cpp`). This certifies a local minimum, not a global one.

Residuals reported with every result (dimensionless): balance rows scaled by feed contribution, per-phase charge, relative pressure, chemical-potential equality over `R T`, affinity over `R T`, KKT stationarity, phase distinctness. Their thresholds are category-3 entries of the numerical-policy header; a breach is a named stop, not a physical conclusion.

## 6. Stability

For a candidate state with phase `0` at `(x⁰, v⁰)` and reduced objective `g(x, v) = a(T, v, x) + P v / (R T)` per mole, the tangent-plane distance in the extended space is

```text
TPD(x, v) = g(x, v) − g(x⁰, v⁰) − ∇_x g(x⁰, v⁰) · (x − x⁰)      (EqID: tangent_plane_distance)
```

The state is stable if and only if `TPD ≥ 0` for every trial `(x, v)` in the EOS domain (\cite{Michelsen1982} in `x`; \cite{Pereira2012} in `(x, v)`; \cite{Perdomo2025} in electroneutral coordinates). A feasible trial with `TPD < 0` is a one-sided instability witness: the candidate is not the global minimum. Completing a finite search without a witness is not the converse proof; results record `global = None` unless the bound loop closes.

Lagrangian dual (\cite{Pereira2012}; \cite{MitsosBarton2007}): `q(λ) = min_{x,v} [ g(x, v) − λ · x ] ≤ G*/N ≤ UBD`, where `UBD` is the value of the linear master problem over the sampled points and `q(λ)` is evaluated by the inner minimization. The inner minimization is the one-phase equilibrium minimizer with the linear term `−λ · n` and no balance rows (`EqID: dual_inner_problem`). The loop terminates when `UBD − q(λ)` is below the policy tolerance or a budget ends it (named stop). The kernel is specified in `phase-stability.md` once designed.

## 7. Specification patterns

The caller declares which variables are specified. The active equation set follows; nothing else selects it.

| specified | unknown | equations added | count check |
|---|---|---|---|
| `T`, `P`, feed, `Π` phases with supports | `n_p`, `V_p` | balance, charge, pressure equality, stationarity | closed system at two intensive variables: determined (Duhem) |
| `T`, reference-phase composition, one incipient phase | `P`, `x_inc`, `V_ref`, `V_inc` | pressure equality (2), `μ_inc,i = μ_ref,i` on the incipient support, normalization | unknowns `1 + C_inc + 2`, equations `2 + C_inc + 1`: square |
| `P`, reference composition, incipient phase | `T`, `x_inc`, `V`'s | same with `T` free | square |
| `T` and `P` both fixed with an incipient phase | — | — | over-specified by one: fails "free T or P" |
| `T` and `P` both free with an incipient phase | — | — | under-specified by one: fails "fix T or P" |
| `T`, `P`, feed, reactions, one phase | `n` | balance (rank `C − R`), charge, pressure equality, stationarity | determined; affinity follows |
| `T`, `P`, feed, reactions, `Π` phases | `n_p`, `V_p` | all rows | determined |
| `T`, `P`, feed, topology omitted | `Π`, then as row 1 | phase stability supplies `Π` | discovery precedes the minimizer |
| feed on an affine coordinate `q`, one target equation (e.g. `P = P_target`) | `q` plus row 1 unknowns | target row | square when one coordinate per target |

Charged species add the `phase_electroneutrality` row and the electroneutral chart; nothing else changes. Pure-component saturation is the second row with one component and two phases (`P` unknown at fixed `T`, or `T` at fixed `P`); the density solve is the pressure-equality row alone for one phase at fixed composition. A request whose unknown count and independent equation count differ fails naming the count and the missing or extra specification; it never defaults.

Cases this specification cannot express: solids (different state variables), critical points (require second-order conditions as equations), and kinetics. They need their own formulation, not a hidden branch.

## 8. Implicit sensitivity

Let `F(y, θ) = 0` collect the active stationarity and constraint rows in the unknowns `y` (amounts, volumes, the free intensive variable, multipliers) with parameters `θ` (the active EOS coordinates). Then

```text
dy/dθ = − F_y⁻¹ F_θ                                                   (EqID: implicit_sensitivity)
```

`F_y` is the KKT matrix assembled from the Helmholtz tape's second derivatives (`∂²A/∂n∂n`, `∂²A/∂n∂V`, `∂²A/∂V²`) and the constraint Jacobians; `F_θ` from the mixed derivatives `∂²A/∂n∂θ`, `∂²A/∂V∂θ`, and `∂g°/∂θ`. Rank is judged by singular values against `ε · max(1, dim) · σ_max` and the condition number is reported; a deficient or ill-conditioned system makes the derivative unavailable and leaves the values available (today `populate_sensitivities`, `solved_state_sensitivity.cpp`). Output rows (pressure, amounts, mole fractions, densities, chemical potentials) are linear or rational in `y` and carry their own chain rule.

## 9. Registry entries

The identifiers above are registered with `implementation_refs` pointing at the current owners: `reactive_finite_phase.cpp` (`evaluate_phase`, `evaluate_multiphase`, `lagrangian_gradient`, `reduced_hessian_minimum`, `populate_sensitivities`), `chemical_equilibrium.cpp` (`make_amount_chart`), `phase_direction_compiler.cpp` (`compile_phase_direction_system`), `grepe.cpp` (`incipient_closure`), `held2_step2.cpp` (`evaluate_held2_tpd`), `held2_step1.cpp` (`make_held2_coordinates`), `solved_state_sensitivity.cpp`. The rewrite moves the refs to `gibbs_minimizer.cpp` and `phase_stability.cpp`; the identifiers do not change.

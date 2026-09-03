# Equilibrium formulation

Status: scientific equation authority for phase, chemical, and coupled reactive-phase
equilibrium in the greenfield engine. Equation identifiers below are registered in
`engine/docs/equations_registry.yaml`; runtime source is the executable authority and
this document creates no second implementation. Re-landed 2026-09-02 from PR #194 with
the revisions `engine/docs/greenfield/EQUILIBRIUM.md` §6 lists; that document
owns the activation rule, the rejection table and the checks, and is cited here rather
than repeated.

## 1. One problem

Phase equilibrium, chemical equilibrium, and coupled reactive-phase equilibrium are one
constrained minimization. Equality of chemical potential across phases and vanishing
reaction affinity are both first-order stationarity conditions of that minimization
under different constraint rows. Nonreactive equilibrium is the case with no reaction
span; single-phase chemical equilibrium is the case with one phase; the coupled problem
has both. Nothing below is derived separately for the three cases, and nothing below
names a problem kind.

Sources: \cite{Pereira2012} (HELD: Helmholtz free Energy Lagrangian Dual) for the
Helmholtz-coordinate objective, the tangent-plane criterion in volume–composition space,
and the Lagrangian dual bound loop; Michelsen 1982 for the tangent-plane stability
criterion; Mitsos and Barton 2007 for the dual extremum principle; Perdomo et al. 2025 for
electroneutral reduced coordinates in the stability search; \cite{Gross2001} for pressure
and chemical potential from the residual Helmholtz energy; \cite{Ascani2022} and Tsanas et
al. 2019, 2022 for the electroneutrality rows of multiphase electrolytes.

## 2. Free energy and coordinates

The EOS supplies the residual Helmholtz energy `a^r(T, ρ, x)` and its exact derivatives.
For a phase `p` with amounts `n_p` (mol) and volume `V_p` (m³) at temperature `T`:

```text
A_p(T, V_p, n_p) / (R T) = Σ_i n_{p,i} [ ln(n_{p,i} / V_p / ρ_ref) − 1 ] + N_p a^r(T, N_p / V_p, n_p / N_p),   N_p = Σ_i n_{p,i}
```

`EqID: phase_helmholtz_amount_volume`. Source: \cite{Gross2001} Appendix A with the
ideal-gas mixing term; `ρ_ref = 1 mol m⁻³` contributes a species-linear constant only.

The system objective at fixed `T` and `P` is the Gibbs energy expressed in Helmholtz
coordinates:

```text
G(T, P, {n_p, V_p}) / (R T) = Σ_p [ A_p(T, V_p, n_p) / (R T) + P V_p / (R T) ] + Σ_p Σ_i n_{p,i} g°_i / (R T)
```

`EqID: system_gibbs_helmholtz_coordinates`. Source: \cite{Pereira2012} Eqs. 2–4; exact
derivation from `G = min_V [A + P V]`. `g°_i` is the standard-state chemical potential
offset supplied by a reaction record (zero when no reaction references species `i`).

**Helmholtz versus Gibbs.** `G(T, P, n) = min_V [A(T, V, n) + P V]` for a mechanically
stable phase. Minimizing `A + PV` jointly over `(n_p, V_p)` and minimizing `G` over `n_p`
with `V_p` eliminated by a density solve are the same problem; the joint form makes the
pressure-equality condition `∂A_p/∂V_p = −P` a stationarity row instead of an inner root
solve and needs no density branch selection inside the solver. Every part of the engine
uses it: the minimizer, the phase-stability inner problem, every boundary, and pure
saturation. The dual in \cite{Pereira2012} is taken over the material-balance multipliers
of the same objective (§6); its inner problem is the single-phase form of the same
objective with a linear term, so the discovery kernel reuses the minimizer with fixed
multipliers.

**Coordinates.** Each phase carries its **amount** `N_p`, its **composition** `x_p` on
the support `S_p` (as `ln x_{p,i}`), and its **volume** as the logarithmic packing
fraction `ln η_p` (`EqID: phase_packing_volume`: `V_p = (π/6) Σ_i n_{p,i} m_i d_i³ / η_p`,
exact from the packing-fraction definition). The amount is separate from the composition
so that pinning it — to zero for an incipient phase, to a fraction for a fixed-fraction
flash — is a specification that removes one variable rather than a limit that breaks one;
`n_{p,i} = N_p x_{p,i}` is a linear change of variables in log space. Bounds
`0 < η_p < η_max` keep every trial inside the EOS domain; the logarithm keeps compositions
positive without inequality constraints. The EOS is evaluated at unit amount, `(T, v_p,
x_p)`, and extensive results are scaled by `N_p` afterwards, so that a pinned phase is a
plain zero in every sum and no `ln 0` is ever formed. Charged phases use the same
coordinates; the electroneutral chart of Perdomo et al. 2025 belongs to the stability
search (`phase-stability.md`), not to the minimizer.

## 3. Variables and their statuses

| symbol | meaning | status the specification may give it |
|---|---|---|
| `x_p` | composition of phase `p` on `S_p` | always free |
| `N_p` | amount of phase `p` | free, or pinned to `β_p N_feed` (`β_p = 0`: incipient) |
| `V_p` (as `ln η_p`) | phase volume | always free |
| `T`, `P` | temperature, pressure | each fixed or free |
| `Π` | phase count | declared, or omitted (discovery, §6) |

Reaction extents are not variables: reactions enter through the balance matrix `B`
(conserved moieties), whose nullspace is the reaction span, and through `g°`. Extents are
recoverable as `ξ = ν⁺ (n − n_feed)`. There is no separate variable kind for an incipient
phase.

## 4. Rows

| id | row | count |
|---|---|---|
| `EqID: chemical_potential_stationarity` | `μ_{p,i} − λ·B_i − z_i φ_p = 0`, `i ∈ S_p`, written intensively | `Σ_p |S_p|` |
| `EqID: pressure_equality` | `−∂A_p/∂V_p − P = 0` | one per phase |
| `EqID: material_balance` | `Σ_p N_p B x_p = B n_feed` (`B = I` when no reactions) | `rank B` |
| `EqID: phase_electroneutrality` | `z · x_p = 0` | one per ion-containing phase except one phase of nonzero amount; all of them if every ion-containing phase is pinned |
| `EqID: phase_normalization` | `Σ_i x_{p,i} = 1` | one per phase |
| `EqID: phase_support` | `n_{p,i} = 0` for `i ∉ S_p` (structural: removes the variable) | — |

The normalisation row is active for every phase, so the coordinates are the same for a
free and a pinned phase and the two differ by exactly one removed variable. The
chemical-potential row is the same equation for both because `μ` is intensive.

**Why one charge row fewer than the ion-containing phases.** The feed is neutral and every
admitted reaction conserves charge, so `z` lies in the row space of `B` and overall
neutrality `Σ_p N_p z·x_p = 0` is a consequence of the balance rows. With `π_ch`
ion-containing phases of nonzero amount, `π_ch − 1` rows are independent, their multipliers
`φ_p` being the electrostatic potential differences between phases (\cite{Ascani2022}
eq. 5; Tsanas et al. 2019 §2.2.1 and 2022 §2.4). A phase pinned at zero amount drops out
of the balance sum, so the row-less phase must have nonzero amount.

Well-posedness of the declared rows is proved before any solve by the compile step
(`EQUILIBRIUM.md` §2.3 and §3.2): the count rule `rows − unknowns = k − f`, the
phase rule in rank form `Π − k ≤ rank B`, the rank of the assembled Jacobian at a generic
point, and the named rejections. No second rank routine exists.

## 5. Stationarity and KKT

With multipliers `λ` on the balance rows and `φ_p` on the charge rows, the rows of §4 are
the KKT conditions of §2's objective. Sources: \cite{Pereira2012} Eqs. 5–8;
\cite{Gross2001} Eq. A.33 for `μ_i^r`. Consequences that are not imposed a second time:

- **Phase equilibrium.** For a species in two supports, `μ_{p,i} − μ_{q,i} = z_i (φ_p −
  φ_q)`: equal chemical potential for neutral species, equal electrochemical potential for
  ions (`EqID: interphase_potential_equality`).
- **Chemical equilibrium.** For any reaction `ν` in the nullspace of `B`, `Σ_i ν_i μ_{p,i}
  = λ · (B ν) + φ_p (z · ν) = 0` because `B ν = 0` and `z · ν = 0` (`EqID:
  reaction_affinity`). Affinity vanishes as a consequence of the balance multipliers.
- **Euler.** `g_p = λ·B x_p + φ_p z·x_p` for every phase: for a pinned phase, the statement
  that it is at the threshold of existing (`TPD = 0`); for a free phase, the equilibrium
  condition.
- **Coupled.** Both hold with the same `λ`, `φ`; and `G_min = λ · B n_feed`.

Second-order local condition: the reduced Hessian `Zᵀ ∇²L Z` over the constraint
nullspace `Z` is positive definite (`EqID: reduced_hessian_stability`). This certifies a
local minimum, not a global one.

Residuals reported with every result (dimensionless): balance rows scaled by feed
contribution, per-phase charge for every ion-containing phase, relative pressure,
chemical-potential rows over `R T`, affinity over `R T`, the reduced Hessian's smallest
eigenvalue, phase distinctness. Their thresholds are category-3 entries of the
numerical-policy header; a breach is a named stop, and a declared phase set that resolves
to fewer distinct phases than declared is a failure with no equilibrium state
(`EQUILIBRIUM.md` §3.4).

## 6. Stability

For a candidate state with phase `0` at `(x⁰, v⁰)` and reduced objective `g(x, v) = a(T,
v, x) + P v / (R T)` per mole, the tangent-plane distance in the extended space is

```text
TPD(x, v) = g(x, v) − g(x⁰, v⁰) − ∇_x g(x⁰, v⁰) · (x − x⁰)      (EqID: tangent_plane_distance)
```

The state is stable if and only if `TPD ≥ 0` for every trial `(x, v)` in the EOS domain
(Michelsen 1982 in `x`; \cite{Pereira2012} in `(x, v)`; Perdomo et al. 2025 in
electroneutral coordinates). A feasible trial with `TPD < 0` is a one-sided instability
witness. Completing a finite search without a witness is not the converse proof; results
record `global = None` unless the bound loop closes.

Lagrangian dual (\cite{Pereira2012}; Mitsos and Barton 2007): `q(λ) = min_{x,v} [ g(x, v)
− λ · x ] ≤ G*/N ≤ UBD`, where `UBD` is the value of the linear master problem over the
sampled points and `q(λ)` is evaluated by the inner minimization. The inner minimization
is the one-phase form of the minimizer with the linear term `−λ · n` and no balance rows
(`EqID: dual_inner_problem`). The loop terminates when `UBD − q(λ)` is below the policy
tolerance or a budget ends it (named stop). The kernel is specified in
`phase-stability.md`. Discovery is the same row system with the normalisation row as a
complementarity rather than an equality (Michelsen 1994): a phase absent at equilibrium
satisfies `Σ_i x_{p,i} exp(−r_{p,i}/RT) < 1`, with `r` its chemical-potential residuals,
where a present phase satisfies `= 1`.

## 7. Specification

The caller declares the status of every variable in §3; the active rows follow from §4 by
one rule and nothing else selects them. With `k` pinned phases and `f` free intensives,
`rows − unknowns = k − f`: the flash is `k = f = 0`; a boundary — bubble, dew, cloud
point, pure saturation at `C = 1`, the fixed-fraction flash — is `k = f = 1`,
distinguished only by which phase is pinned and which intensive is free; a three-phase
point at fixed reference composition is `k = f = 2`. `f ≠ k`, `k = Π`, `Π − k > rank B`, a
dependent reaction set, a charge-violating reaction, a charged feed, and a support that
cannot be neutral are rejected naming the deficiency; the engine never defaults. The full
table, the compile-tier checks and the mutations are `EQUILIBRIUM.md` §2.3, §3.2
and §4.1.

Cases this specification cannot express: solids (different state variables), critical
points (second-order conditions as equations), and kinetics. They need their own
formulation, not a hidden branch.

## 8. Implicit sensitivity

Let `F(y, θ) = 0` collect the active rows in the unknowns `y` (compositions, volumes,
free amounts, the free intensive, multipliers) with parameters `θ` (the active EOS
coordinates). Then

```text
dy/dθ = − F_y⁻¹ F_θ                                                   (EqID: implicit_sensitivity)
```

`F_y` is the matrix the driver factorised at the solution, assembled from the Helmholtz
tape's second derivatives (`∂²A/∂n∂n`, `∂²A/∂n∂V`, `∂²A/∂V²`) and the constraint
Jacobians; `F_θ` from the mixed derivatives `∂²A/∂n∂θ`, `∂²A/∂V∂θ`, and `∂g°/∂θ`. No
third derivative appears. Rank is judged by singular values against `ε · max(1, dim) ·
σ_max` and the condition number is reported; a deficient or ill-conditioned system makes
the derivative unavailable and leaves the values available. Output rows (pressure,
amounts, mole fractions, densities, chemical potentials, activities) are functions of `y`
and carry their own chain rule; a regression consumes `dy/dθ` through them.

## 9. Registry entries

Registered by this document: `phase_helmholtz_amount_volume`,
`system_gibbs_helmholtz_coordinates`, `phase_packing_volume`,
`chemical_potential_stationarity`, `pressure_equality`, `material_balance`,
`phase_electroneutrality`, `phase_normalization` (replaces the retired
`incipient_normalization`, which named a row that is now universal), `phase_support`,
`interphase_potential_equality`, `reaction_affinity`, `reduced_hessian_stability`,
`tangent_plane_distance`, `dual_inner_problem`, `implicit_sensitivity`. Registered by
`phase-stability.md`: `electroneutral_amount_chart`, `modified_fraction_chart`,
`multiplier_potential_link`, `outer_approximation_bound`, `candidate_phase_criteria`,
`helmholtz_lagrangian_dual`, `tunneling_transform`. Owners are the `// EqID:` annotations
in `native/equilibrium/`, added by the stage that writes each row; the registry carries
no owner list.

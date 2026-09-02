# Stage 3 — `engine/native/model/`

You are implementing Stage 3 of a greenfield rebuild in this worktree (branch `stage-2`).
You have no one to ask. Read, in order, before writing anything:

1. `engine/docs/greenfield/TESTING.md` — a hard restriction on what may be tested
2. `engine/docs/greenfield/TRANSFER.md` §"Stage 3" and §"Naming"
3. `engine/docs/greenfield/STAGE-LOG.md` — Stage 1's decisions; you build on them
4. `engine/native/physics/*.hpp` — all eight files, completely
5. `engine/tests/native/physics_checks.cpp` — the check pattern you follow
6. `tools/check-engine-architecture.py` — the five rules your result must satisfy
7. `engine/docs/greenfield/briefs/stage-2.brief.json` — sources, budgets, gates

Verify every source marker and gate command in the JSON against the tree **before**
writing. If a claim is false, stop and say which.

## 0. Non-negotiable

**Layers.** `physics <- model <- equilibrium <- regression <- bindings`. A file under
`native/model/` may include `"physics/..."` and `"model/..."` only, rooted at
`engine/native`, never relative.

**No pybind11, no Python** outside `native/bindings/`. `_core` does not exist until Stage 3.

**Constants and tolerances.** Category-1 constants come from `physics/physical_constants.hpp`.
Every threshold, budget, step factor or fraction you introduce is a `Tolerance` entry in
`physics/numerical_policy.hpp`, with category and source, in the existing style. No numeric
literal acting as a tolerance may appear in `model/`. The checker rejects any `constexpr`
name defined twice — give each a unique name.

**Naming.** No proper nouns in identifiers. Published equation names Stage 1 already uses
(`debye_huckel`, `born`) are fine; algorithms, authors, products and stage numbers are not.
Cite sources in comments and `// EqID:` annotations only.

**Testing.** Before writing a check ask: would it survive a legitimate reimplementation of
the module from its equations? Only closed forms, exact symmetries, limits, derivative
agreement and real data qualify. No golden files, no assertions on internals, no iteration
count as a correctness claim. Counts and timings are reported, never pass/fail.

**`NDEBUG`.** `buildtype=release` with `b_ndebug=if-release` compiles `assert` to nothing.
Use the `CHECK` macro pattern from `physics_checks.cpp` (fprintf + `std::abort`).

**Environment is not yours to repair.** You may edit `engine/meson.build` only to add your
sources to `model_sources` and to add the two test executables. You may **not** touch its
dependency block, `default_options`, the virtual environment, `engine/build/subprojects`,
or how you were invoked, in order to make a build succeed. That is the orchestrator's
problem: stop with `environment_broken`.

**Physics is closed.** Two edits permitted, no others: (a) *adding* entries to
`numerical_policy.hpp`; (b) an optional trailing `std::size_t* iterations = nullptr` on
`solve_site_fractions` in `association.hpp`. Do not change any kernel's arithmetic, any
existing tolerance value, or any existing signature.

## 1. Stop conditions — stopping is a first-class outcome

- `source_contradicted` — a source lacks its stated marker, or says something other than
  this brief claims. Name source and marker.
- `gate_impossible` — a gate cannot run as written, or cannot be made to fail by its named
  mutation. Name the gate and why.
- `budget_ceiling` — a file would exceed its ceiling. Stop *before* exceeding it.
- `environment_broken` — configure, compile or a precondition fails for a reason outside
  the tree. Quote the command and its first error line. Do not work around it.

An implementer that quietly absorbs a wrong plan is the failure. A correct stop is not.

## 2. Design

### 2.1 The tape — `helmholtz_tape.hpp`

One `CppAD::ADFun<double>` per (model, active-parameter set, recording-time predicate
signature), reused across points.

**Coordinates.** Independents `y = (T, V, n_1..n_C, theta_1..theta_P)`: temperature K,
volume m^3, amounts mol, active parameter values. Width `w = 2 + C + P`. Dependent 0 is the
extensive residual Helmholtz `Phi = A_res/(RT)`, dimensionless. Site fractions are **not**
independents: they enter as CppAD *dynamic parameters* (`CppAD::Independent(y, dynamic)`,
then `f.new_dynamic(X)` per point), so the tape records once and is re-pointed.

**Site elimination.** With `rho = n_tot/V`, `x = n/n_tot`, `m_a = x_{c(a)} nu_a`,
`rho_N = rho * NUMBER_DENSITY_PER_ANGSTROM3`, `Delta_ab` the physics association strength,
`B_a = sum_b m_b X_b Delta_ab`, and `H_a = X_a(1 + rho_N B_a) - 1` the physics
`association_residuals`, define the association potential

    Q(y; X) = sum_a m_a (ln X_a - X_a + 1) - (rho_N/2) sum_a m_a X_a B_a

At the physics solution `X*(y)` (where `H = 0`): `Q = a_assoc` and `dQ/dX_a = 0`. Both
follow from `a_assoc - Q = 0.5 sum_a m_a H_a`, which lets you build `Q` from physics calls
alone: `residual_helmholtz(...).association - 0.5 sum_a m_a H_a`. Then

    Phi(y; X)   = n_tot * [a_hc + a_disp + a_dh + a_born](T, rho, x) + n_tot * Q(y; X)
    F_a(y; X)   = dPhi/dX_a = -(n_tot m_a / X_a) * H_a

Dependents are `(Phi, F_1..F_S)`. With `a(y) := Phi(y; X*(y))`:

    grad a  = Phi_y                                  (envelope: Phi_X = 0)
    hess a  = Phi_yy - Phi_yX * Phi_XX^-1 * Phi_Xy   (Schur complement)
    Phi_XX  = -n_tot [ delta_ab m_a / X_a^2 + rho_N m_a m_b Delta_ab ]   (analytic, S x S)

This is TRANSFER's "Lagrangian a + lambda.H" with the multiplier that makes it a potential.
The tape never carries `lambda` or `X_y`.

**Sweeps — never call `ADFun::Hessian`.** It discards the `Forward(1)` outputs you need,
forcing `S` extra reverse sweeps. Drive the sweeps yourself; count each `Forward`/`Reverse`
as one:

| order | sweeps | how |
|---|---|---|
| value | 1 | `Forward(0, y)` |
| gradient | 2 | + `Reverse(1, e_0)` -> `Phi_y` |
| pressure slope | 4 | + `Forward(1, e_V)` (gives `dPhi/dV` **and** `Phi_XV` for all S dependents) + `Reverse(2, e_0)`; `a_VV = Phi_VV - Phi_VX Phi_XX^-1 Phi_XV` |
| hessian | `1 + 2w` | `Forward(0)`; per `j`: `Forward(1, e_j)` (column j of `Phi_Xy`), `Reverse(2, e_0)` (row j of `Phi_yy`); then one `S x S` solve and the Schur product |

At nine species, `S = 6`, `P = 1`: `w = 12`, hessian = 25 sweeps and one 6x6 solve, against
`7 * 35 = 245` today. Report `sweeps`, `hessian_calls`, `recordings` via a counters struct
on the model — the model layer owns evaluation accounting.

**Active parameters.** `ActiveParameter { ParameterFamily family; std::vector<std::size_t>
indices; }`, indices being flat positions in that family's array sharing one value (a
symmetric `k_ij` pair has two). Families at Stage 3: `segment_count`, `sigma_constant`,
`dispersion_energy_over_k`, `k_ij`, `association_energy_over_k`, `association_volume`.
Implement with a `TapeParameters` mirror of `Parameters` whose arrays for those families are
`std::vector<AD>`; the physics kernels are templated on `ParameterSet` and accept it.
`electrolyte_contribution` and `born_contribution` take `const Parameters&` and get the
double record. Other families are not activatable at Stage 3 because physics reads them as
`double` — say so in a comment, do not extend physics.

**Cache key — read this twice.** Two predicates in `physics/electrolyte.hpp` are evaluated
from `double` inputs *while recording* and baked into the tape: the `finite_size_chi` series
switch per ion (`reference_kappa * 1e-10 * d_i < 1e-3`, line 332 — pass the evaluation point
as the reference) and the `reference_charge_square_sum == 0.0` branch (line 242). Both flip
inside Stage 3's own scope: the dilute end of the molality ladder takes the series, a
pure-solvent state inside an electrolyte model takes the no-ion branch. A tape keyed on
topology alone is silently wrong across them, with no error raised. Key the cache on
(active set, predicate signature) and re-record when it changes. Check 3.14 exists for this.

**Order per point.** Always run the double path first — `segment_sigma`,
`segment_diameter`, `solve_site_fractions` (where domain errors throw) — then
`new_dynamic(X*)` and the sweeps. `evaluate` is not thread-safe (an `ADFun` holds sweep
state); guard with a `std::mutex` member and say so.

### 2.2 The model — `helmholtz_model.{hpp,cpp}`

```cpp
enum class DerivativeOrder { Value, Gradient, PressureSlope, Hessian };
struct HelmholtzDerivatives {           // in (T, V, n, theta); w = 2 + C + P
  double value;                         // Phi = A_res/RT
  std::vector<double> gradient;         // w, empty unless requested
  std::vector<double> hessian;          // w*w row-major, empty unless Hessian
  double volume_second_derivative;      // a_VV, set for PressureSlope and Hessian
  std::vector<double> site_fractions;   // S
};
struct SweepCounters { std::size_t forward, reverse, hessian_calls, recordings; };

class HelmholtzModel final {
 public:
  explicit HelmholtzModel(Parameters, std::vector<ActiveParameter> = {});  // validates once
  HelmholtzDerivatives evaluate(double T, double V, const std::vector<double>& n,
                                DerivativeOrder) const;
  State state(double T, double molar_density, const std::vector<double>& x) const;
  State state_at_pressure(double T, double P, const std::vector<double>& x,
                          DensityBranch, std::optional<double> anchor_density) const;
  double packing_fraction_coefficient(double T, const std::vector<double>& x) const;
  double phase_packing_fraction(double T, const std::vector<double>& n, double V) const;
  double molar_density_upper_bound(double T, const std::vector<double>& x) const;
  SweepCounters counters() const;  void reset_counters() const;
};
```

`State` carries the five term values, `bulk_relative_permittivity`, `residual_helmholtz`,
`compressibility_factor`, `pressure`, `residual_chemical_potential_over_rt[]`,
`log_fugacity_coefficient[]` (absent when `Z <= 0`), `residual_entropy`,
`residual_internal_energy`, `residual_enthalpy`, `residual_isochoric_heat_capacity`,
`residual_isobaric_heat_capacity` (optional), `molar_density`,
`pressure_density_derivative`, `site_fractions[]`. `state(T, rho, x)` is
`evaluate(T, 1/rho, x, Hessian)` with `n_tot = 1`, mapped by (carry the EqIDs from eos.cpp):

    Z = 1 - V Phi_V              P = RT (1/V - Phi_V)          mu_res,i/RT = Phi_{n_i}
    ln phi_i = Phi_{n_i} - ln Z  S_res/R = -Phi - T Phi_T      U_res/RT = -T Phi_T
    H_res/RT = -T Phi_T + Z - 1  Cv_res/R = -2T Phi_T - T^2 Phi_TT
    dP/drho|_{T,x} = RT (1 + V^2 Phi_VV)        (dP/dT)_V = P/T - RT Phi_TV
    (dP/dV)_T = RT (-1/V^2 - Phi_VV)            Cp_res from Cv_res and the two P derivatives

Every second derivative is the *reduced* one. Two identities you will check: Euler
homogeneity `V Phi_V + sum n_i Phi_{n_i} = Phi`, and `sum x_i mu_res,i/RT = a_res + Z - 1`.

**Validation, once, in the constructor** (`std::invalid_argument` naming the field): array
lengths against `C`, `C^2`, `S`, `S^2`; `segment_count > 0`, `sigma_constant > 0`,
`dispersion_energy_over_k >= 0`, `molar_mass > 0`, all finite; `k_ij`, `l_ij` and the three
association matrices symmetric; `site_component` in range; `association_active` in {0,1}
with energy and volume `>= 0` where active; every `sigma_term_component` and
`permittivity_term_component` in range; charged components have positive `packing_diameter`
and `debye_huckel_diameter`, and under `DebyeHuckelBorn` also positive `born_diameter` and
`ionic_region_relative_permittivity`; a `permittivity_form` other than `None` requires at
least one dielectric solvent; `ion_suppression_coefficient` sized `C`. Point validation in
`state`/`evaluate`: `T > 0`, `V > 0`, `n_i >= 0`, `sum x = 1` within a new policy entry, all
finite. Do **not** require electroneutrality — later stages evaluate charged trial states.

Helpers, from `main`'s `helmholtz.hpp` lines 155-198, using greenfield constants
(`NUMBER_DENSITY_PER_ANGSTROM3`, not `PHYSICAL_AVOGADRO`): `packing_fraction_coefficient`
= `packing_fraction_from_geometry(1.0, x, p, d(T))` after checking sigma and d positive and
finite; `phase_packing_fraction(T, n, V)` = eta at `rho = n_tot/V`; `molar_density_upper_bound`
= `CLOSE_PACKING_FRACTION / coefficient`.

### 2.3 Density — `density.cpp`

```cpp
enum class DensityBranch { Unspecified, Vapor, Liquid };
struct DensityRoot { double molar_density, pressure_residual, pressure_density_derivative;
                     std::size_t iterations, stable_root_count; DensityBranch branch; };
DensityRoot solve_density(const HelmholtzModel&, double T, double P,
                          const std::vector<double>& x, DensityBranch,
                          std::optional<double> anchor_density);
```

Residual in reduced units, as `DENSITY_ROOT_RESIDUAL` already states:
`r(rho) = Z(rho) - P/(rho RT)`, converged at `|r| <= DENSITY_ROOT_RESIDUAL`. Each evaluation
is `evaluate(T, 1/rho, x, PressureSlope)` — four sweeps, never a Hessian.

1. **Anchor.** If given and in domain, safeguarded Newton on `r(rho)` with `dr/drho` from
   `dP/drho`. Accept when converged and mechanically stable (`dP/drho > 0`, strict — a
   constraint, not a tolerance). On leaving the domain, stalling, or `dP/drho <= 0`, fall to 2.
2. **Cold.** Scan `DENSITY_BRACKET_INTERVALS` (40, exists) points in `log rho` over
   `[rho_low, rho_high]`, `rho_low` a fraction of the ideal-gas density `P/(RT)`, `rho_high`
   a fraction of `molar_density_upper_bound` — two **new** policy entries. Each sign change
   of `r` brackets a root; refine with bisection-safeguarded Newton to
   `DENSITY_ROOT_RESIDUAL` within `DENSITY_MAX_ITERATIONS` (exists). Keep roots with
   `dP/drho > 0`.
3. **Branch.** `Liquid` -> highest-density stable root; `Vapor` -> lowest; `Unspecified` ->
   the unique stable root, else throw naming the count. When exactly one stable root exists,
   both `Liquid` and `Vapor` return it — a deliberate change from `main`, which threw; a
   sweep crossing the critical region must not fail. Record it as a deviation. No root: throw.

Not written: the 360-interval scan, the 64 x depth-8 envelope, `polished_pressure_density`'s
`nextafter` walk, the `std::string` branch names, the certified-root-index continuation API.
Association non-convergence at a scan point is skipped as `main` did; other exceptions propagate.

### 2.4 Infinite-dilution reference — `infinite_dilution_reference.cpp`

The ladder `1e-6, 1e-8, 1e-10, 1e-12`, the witness `1e-14`, and the `5e-5` threshold move
into `numerical_policy.hpp` (category 3, source `electrolyte_reference.hpp` on `main`). At
each molality: `formula_composition` (EqID `formula_unit_composition`), then a `Liquid`
`state_at_pressure` **anchored on the previous root** — composition moves by 1e-6 between
points, so the anchor is nearly exact and each point costs a few four-sweep evaluations.
`convergence_error` = max over the two ions of `|delta ln phi|` between the last two ladder
points; above threshold, evaluate the witness once and fail if still above. Pure-solvent
`ln phi` when there is one solvent. Activity per `reference.py` 426-490:
`ln gamma_i = ln phi_i(m) - ln phi_i^inf`; `ln gamma_pm^x = sum nu_i ln gamma_i / sum nu_i`;
`gamma_pm^m = gamma_pm^x / (1 + M_s m sum nu)`; osmotic
`phi = -[ln x_s + ln phi_s(m) - ln phi_s^pure] / (sum nu * m * M_s)`.

Not written: the four-fixed-P-states-each-a-scan structure, the neutral-basis rank
machinery, the request/result/status types.

## 3. Checks — `model_checks.cpp`, meson test `model`

Scaled coordinates: scale `T` by `T`, `V` by `V`, `n_i` by 1, `theta_k` by `|theta_k|`, so
every Hessian entry is `O(Phi)`; then `close(a,b,tol) := |a-b| <= tol*max(1,|b|)`. FD uses
relative step 1e-5 (absolute 1e-5 for amounts), centred, on strictly positive compositions.
Run on all four fixtures unless stated.

1. **Tape reproduces the double path.** `evaluate(...).value` equals `n_tot` times the sum
   of the five physics terms on the double path, to 1e-14. Pins `Q = a_assoc` at `X*`.
2. **Gradient vs FD**, every coordinate including theta, 1e-6 scaled.
3. **Hessian equivalence — the check that licenses the 7x.** On the associating binary and
   the nine-species case with one active parameter (water `dispersion_energy_over_k`): the
   model's reduced Hessian equals the `(S+1)`-Hessian reference tensor computed in
   `model_support.hpp` elementwise to **1e-10 scaled**. Write the reference independently of
   `model/`: independents `z = (T, V, n, theta, X)` (width `w+S`); dependents
   `G_0 = n_tot*[a_hc + a_disp + a_dh + a_born + a_assoc(X)]` using the **standard**
   association form straight from `residual_helmholtz` with X on the tape, and `G_a = H_a`
   from `association_residuals`; `J = Jacobian(z)`; `X_y = -H_X^-1 H_y`; `S+1` calls to
   `ADFun::Hessian(z, k)`; `hess a[l][r] = d_l^T hess G_0 d_r + sum_a A_X[a] * X_lr[a]` with
   `d_l = (e_l ; X_y[:,l])` and `X_lr = -H_X^-1 (d_l^T hess G_a d_r)_a`. Also compare its
   gradient `A_y + A_X X_y` to the model's, 1e-10. This is evidence rather than a
   restatement because the two sides share only the Stage-1-checked kernels and CppAD: they
   differ in formulation (standard `a_assoc` with sites as independents and the implicit
   function theorem, versus potential plus Schur with sites eliminated) and in every line of
   assembly. Both are exact, so 1e-10 is a real threshold, not FD noise.
4. **Hessian vs FD of the gradient**, 1e-6 scaled, all fixtures. Guards the failure both
   exact formulations share: a quantity that varies on the double path recorded as a tape
   constant.
5. **Hessian symmetric**, 1e-12 scaled.
6. **Euler homogeneity** `V Phi_V + sum n_i Phi_{n_i} = Phi`, 1e-12; and from `state`,
   `sum x_i mu_res,i/RT = a_res + Z - 1`, 1e-12.
7. **Pressure slope.** `pressure_density_derivative` vs centred FD of `pressure` in rho,
   1e-6, all fixtures — exercises the 4-sweep path with its Schur correction.
8. **Permutation.** Relabel the nine-species components, sites, pairs and edges: `a_res`,
   `P`, `Z` unchanged to 1e-12; `mu_res` permuted accordingly.
9. **Zero-density limit.** `a_res -> 0`, `Z -> 1` on the neutral binary.
10. **Parameter derivative.** `dPhi/dtheta` vs FD over a rebuilt model with `theta +- h` (all
    indices of the family moved together), 1e-6, for a `k_ij` pair on the neutral binary and
    an association energy edge on the associating binary.
11. **Density roots against NIST.** Ethane at 200 K and 240 K, `P` = the oracle
    `saturation_pressure_pa` for that row: `Liquid` and `Vapor` give two distinct roots, each
    `|P(rho*) - P|/P <= 1e-9` with `dP/drho > 0`, liquid denser, and `1/rho*` within 5% of
    the row's molar volumes. `Unspecified` throws there. Anchored solves from 1.05x and 0.95x
    each cold root return it to 1e-10 relative. Ethane at 320 K, 1e6 Pa (supercritical): all
    three branches return the same single root. Nine-species `Liquid` at 313.15 K, 1e5 Pa
    converges with `eta < CLOSE_PACKING_FRACTION`.
12. **Reference: the Debye-Huckel tail.** Aqueous NaCl at 298.15 K, 1e5 Pa, solvent fractions
    (1,0,0), ions (1,2), stoichiometry (1,1). The recorded `convergence_error` equals
    `A_gamma (sqrt(m3) - sqrt(m4))` to **1e-3 relative**, where m3, m4 are the last two ladder
    molalities and `A_gamma = e^2/(8 pi eps0 eps_r k T) * sqrt(2 e^2 N_A rho_w/(eps0 eps_r k T))`
    with `rho_w` the model's own pure-water liquid density at (T,P) times `M_w`, and `eps_r`
    the physics `bulk_relative_permittivity` of pure water there. This is the exact `kappa -> 0`
    limit; every other term is analytic in ion fraction and contributes ~1e-10. On `main`'s
    oracle the recorded error is 1.0628346e-05 and this formula gives 1.06287e-05 at
    `rho_w = 997`. Then `electrolyte_activity(m = 1e-8).mean_ionic_log_molality_scale =
    -A_gamma(sqrt(1e-8) - sqrt(m4))` to 1e-2 relative. Swapping cation and anion order changes
    `ln gamma_pm` by < 1e-12.
13. **Helpers.** `packing_fraction_coefficient(T,x)*rho == packing_fraction_from_geometry` and
    `phase_packing_fraction(T,n,V) == eta(n_tot/V, n/n_tot)`, 1e-14; the constructor throws
    `std::invalid_argument` on a mismatched array length, an asymmetric `k_ij`, and a
    `site_component` out of range.
14. **Cache invariance.** A model that has evaluated point A (NaCl at m = 1e-2), then the pure
    solvent, then m = 1e-12, then A again, agrees with a fresh model at A to 1e-14 on value,
    gradient and Hessian. This is the check for the recording-time predicates in §2.1.

## 4. Performance — `model_floor.cpp`, meson test `model_floor`, suite `floor`

Reported, never pass/fail. `env: {'OMP_NUM_THREADS': '1'}`, timeout 600. Median of 20 after
3 warm-ups, `std::chrono::steady_clock`. One table row per case with `wall_us`, `sweeps`,
`hessian_calls`, `recordings`, and where relevant `density_iterations`,
`association_iterations`. Print `__VERSION__`, the `/proc/cpuinfo` model name, `OMP_NUM_THREADS`.

| case | what | today | floor |
|---|---|---|---|
| A | nine-species, 313.15 K, rho = 3.0e4, one active parameter (water `dispersion_energy_over_k`), `evaluate(Hessian)` | 5.25 ms, 245 sweeps, 7 Hessians | <= 1.0 ms, <= 23 + 2P sweeps, 1 Hessian |
| B | same point, the (S+1) reference tensor from `model_support.hpp` | — | the measured multiplier |
| C | associating binary, 376.88 K, rho = 1e4, x = (0.1385, 0.8615), `state` | 126 assoc. iterations | <= 5 |
| D | nine-species `Liquid` at 313.15 K, 1e5 Pa, cold and anchored | 382 ms | <= 5 ms |
| E | neutral binary `Vapor` at 300 K, 1e5 Pa, cold | 13.4 ms | <= 0.5 ms |
| F | `infinite_dilution_reference`, aqueous NaCl, 298.15 K, 1e5 Pa | 108 ms | 2-5 ms |

A miss is a conversation, not a failure. Report numbers exactly as measured.

## 5. Fixtures — `model_support.hpp`

Functions returning `epcsaft::Parameters`, fully populated (see `physics_checks.cpp::neutral_binary`
for required fields). Transcribed from the sources in the JSON; confirm, do not re-derive.

**Ethane**: m = 1.6069, sigma = 3.5206, eps/k = 191.42, M = 0.03007. **Methane**: 1.0,
3.7039, 150.03, 0.016043. Neutral binary = [methane, ethane], k_ij = 0, FD composition (0.4, 0.6).

**Associating binary** [solute, solvent]: solute m = 4.5354, sigma = 2.6019, eps/k = 204.0438,
M = 0.06108; solvent 2.1945, 2.229, 141.66, 0.01801528; k_ij = -0.052. Sites
[solute a (nu=2), solute b (1), solvent a (2), solvent b (2)]; active symmetric edges (energy K,
volume): solute a-solute b (2383.4744, 0.118488); solute a-solvent b (2093.8222,
0.1540467874409226); solute b-solvent a the same; solvent a-solvent b (1804.17, 0.2039). No
electrolyte, no permittivity. Composition (0.1385, 0.8615).

**Aqueous NaCl** [water, sodium-cation, chloride-anion]: water m = 1.2047, eps/k = 353.95,
M = 0.0180153, sigma(T) = 2.7927 + 10.11*exp(-0.01775 T) - 1.417*exp(-0.01146 T)
(`sigma_constant` 2.7927, two sigma_term entries on component 0), `relative_permittivity`
78.09 constant, solvation 1.5, charge 0, sites [a (1), b (1)] with edge (2425.7, 0.04509).
Na+: z = +1, m = 1, sigma = 2.8232, eps/k = 230.0, packing = DH diameter 2.484416, born 3.445,
solvation 1.0, M = 0.02298977. Cl-: z = -1, m = 1, sigma = 2.756, eps/k = 170.0, packing = DH
2.42528, born 4.1, solvation 1.0, M = 0.0354527. k_ij: water-Na+ -0.3, water-Cl- -0.3,
Na+-Cl- 0.8. `permittivity_form = IonSuppression` with `ion_suppression_coefficient` =
(0, 7.01, 7.01); `electrolyte_form = DebyeHuckelBorn`, `ionic_region_relative_permittivity` = 8.0,
`born_c_shell = born_c_dielectric = 0`. FD composition (0.98, 0.01, 0.01).

**Nine-species MEA**, order [water, monoethanolamine, carbon-dioxide, protonated-monoethanolamine,
carbamate-anion, bicarbonate-anion, carbonate-anion, hydronium-cation, hydroxide-anion],
charges (0,0,0,+1,-1,-1,-2,+1,-1):

| i | component | m | sigma | eps/k | M | packing = DH d | born d | alpha_i | solvation |
|---|---|---|---|---|---|---|---|---|---|
| 0 | water | 1.2047 | sigma(T) as NaCl | 353.95 | 0.01801528 | — | — | 0 | 1.5 |
| 1 | monoethanolamine | 3.0353 | 3.0435 | 277.174 | 0.06108 | — | — | 0 | 1.0 |
| 2 | carbon-dioxide | 2.079 | 2.7852 | 169.21 | 0.04401 | — | — | 0 | 1.0 |
| 3 | protonated-MEA | 1.0 | 3.48508556586 | 232.687201645 | 0.06209 | 3.0668752979568 | 3.53322927146 | 2.6 | 1.0 |
| 4 | carbamate | 1.0 | 3.53543525721 | 453.265244384 | 0.10408 | 3.1111830263448 | 3.54107030822 | 7.01 | 1.0 |
| 5 | bicarbonate | 1.0 | 2.9296 | 70.0 | 0.0610168 | 2.5780480000000003 | 3.0 | 7.89 | 1.0 |
| 6 | carbonate | 1.0 | 2.4422 | 249.26 | 0.06001 | 2.149136 | 3.0 | 7.89 | 1.0 |
| 7 | hydronium | 1.0 | 3.4654 | 500.0 | 0.01902 | 3.049552 | 1.218 | 9.55 | 1.0 |
| 8 | hydroxide | 1.0 | 2.0177 | 650.0 | 0.01701 | 1.775576 | 3.08107689400404 | 13.96 | 1.0 |

Permittivity: water `78.41043478578715 - 0.36133766233691506*(T - 298.15) +
0.00076555618295*(T - 298.15)^2` (`relative_permittivity[0]` = 78.41043478578715,
`permittivity_reference_temperature_k[0]` = 298.15, two terms of power 1 and 2 on component 0);
MEA 32.0; CO2 1.4122; ions 0. `permittivity_form = IonSuppression`,
`electrolyte_form = DebyeHuckelBorn`, `ionic_region_relative_permittivity` = 8.0,
`born_c_shell = born_c_dielectric = 0`. k_ij (symmetric, else 0): MEA-water -0.07352749874985018;
carbamate-protonated-MEA -0.00201813457644; carbonate-water -0.25; hydronium-water +0.25;
hydroxide-water -0.25. l_ij = 0. Sites, all nu = 1: [MEA a, MEA b, water a, water b, CO2 a,
CO2 b]. Active symmetric edges: MEA a-MEA b (2586.3, 0.03747); water a-water b (2425.7, 0.04509);
CO2 a-water b (1212.85, 0.04509); CO2 b-water a (1212.85, 0.04509); MEA a-water b and
MEA b-water a by arithmetic energy and geometric volume: (2506.0, 0.0411037990944876). Floor
composition (0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0). FD/permutation composition:
replace the zeros by 1e-3 each, add 3e-3 to hydronium, renormalise.

## 6. Build

Add the four `.cpp` to `model_sources`. Add two executables in the style of `physics_checks`,
each `link_with: epcsaft_native`, `dependencies: [eigen_dep, cppad_headers_dep]`,
`include_directories: native_inc`; do **not** list `cppad_temp_file` again. `test('model',
model_checks)`; `test('model_floor', model_floor, suite: 'floor', env: {'OMP_NUM_THREADS': '1'},
timeout: 600)`. Then `meson compile`, `meson test physics model`, `meson test --suite floor`,
`python3 tools/check-engine-architecture.py` — all clean, worktree free of stray files.

## 7. What to report

A JSON object: `files[]` {path, lines from `wc -l`}; `gates[]` {id, command, passed} for every
gate in the JSON, run by you; `floor` — the §4 table verbatim plus the environment lines;
`deviations[]` — every departure from this brief and why; `policy_entries_added[]`;
`stopped_because` in {none, budget_ceiling, gate_impossible, source_contradicted,
environment_broken} with `detail`. Do not run the mutations; acceptance runs them. Do not
claim a number you did not measure.

Deferred deliberately, do not write: `parameter_hash` (first consumer is Stage 3), a
tape-side ideal-gas Helmholtz (Stage 4), component and site names (Stage 3 maps names; the
model works on indices).

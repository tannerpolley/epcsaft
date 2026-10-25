# FeOs 0.10.1 controls and returned quantities

The analysis uses the CPython 3.13 wheel, numpy 2.5.3 and si-units 0.11.5,
pinned in `requirements-feos.lock`, installed with:

```bash
uv venv --python 3.13 ~/.cache/epcsaft-comparators/222/feos
uv pip install --python ~/.cache/epcsaft-comparators/222/feos/bin/python \
  --only-binary :all: -r analyses/2026-competitor-cross-check/inputs/requirements-feos.lock
```

Numerical calls run serially without the retired heavy-process lock, with a 180-second limit,
one allowed CPU and priority 10. `set_num_threads(1)` caps FeOs;
`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` and
`JULIA_NUM_THREADS` are 1 before imports. The lead owns affinity, hardware
records and Julia's 6 GB memory cap. Caches remain through delivered review,
per the [owner decision](https://github.com/tannerpolley/ePC-SAFT/issues/222#issuecomment-5956426820).

## Native solve controls

All source locators refer to FeOs tag `v0.10.1`. Default values are preserved
where the selected callable exposes no control.

| Calculation | Entered controls or source defaults | Source locator |
|---|---|---|
| `PhaseEquilibrium.pure` | Entered `tol=1e-12`, `max_iter=1000`; automatic initialization | [pure signature](https://github.com/feos-org/feos/blob/v0.10.1/py-feos/src/phase_equilibria.rs#L59) |
| `PhaseEquilibrium.bubble_point` | Entered `tol_inner=tol_outer=1e-12`, `max_iter_inner=max_iter_outer=1000`; automatic initialization | [bubble signature](https://github.com/feos-org/feos/blob/v0.10.1/py-feos/src/phase_equilibria.rs#L317) |
| `PhaseEquilibrium.tp_flash` | Entered `tol=1e-12`, `max_iter=1000`; automatic native stability initialization; one mole feed | [flash signature](https://github.com/feos-org/feos/blob/v0.10.1/py-feos/src/phase_equilibria.rs#L133) |
| Native density iteration | Verbatim `let (abstol, reltol) = (1e-12, 1e-14);`, `let maxiter = 50;`; pressure-error stopping test in reduced units | [density iteration, lines 110–125 and 223](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos-core/src/density_iteration.rs#L110) |
| Neutral association | Defaults `max_eta: 0.5`, `max_iter_cross_assoc: 50`, `tol_cross_assoc: 1e-10`, `dq_variant: DQVariants::DQ35` | [PC-SAFT options, lines 34–42](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos/src/pcsaft/eos/mod.rs#L34) |
| Revised electrolyte association | Defaults `max_eta: 0.5`, `max_iter_cross_assoc: 50`, `tol_cross_assoc: 1e-10` | [ePC-SAFT options, lines 37–44](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos/src/epcsaft/eos/mod.rs#L37) |
| Exact vapor-pressure property | Defaults `MAX_ITER_PURE: usize = 50`, `TOL_PURE: f64 = 1e-12`; property uses `Default::default()` | [pure constants](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos-core/src/phase_equilibria/vle_pure.rs#L13), [gradient](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos-core/src/ad/properties/vapor_pressure.rs#L35) |
| Exact bubble-pressure property | Defaults `MAX_ITER_INNER: usize = 5`, `TOL_INNER: f64 = 1e-9`, `MAX_ITER_OUTER: usize = 400`, `TOL_OUTER: f64 = 1e-10`; property uses `Default::default()` | [bubble constants](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos-core/src/phase_equilibria/bubble_dew.rs#L16), [gradient](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos-core/src/ad/properties/bubble_point_pressure.rs#L52) |

Association includes analytic paths and a cross-association Newton path;
the latter tests the mass-action vector norm against its tolerance
([association, lines 99–222 and 299](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos/src/association/mod.rs#L99)).
These stopping rules do not imply the issue's agreement targets are achieved.
Python state/equilibrium objects expose no achieved residuals, iteration
counts, site fractions or association-strength matrix. The runner identifies
these absent diagnostics. Pressure, exact pressure-density derivative and
component fugacities allow the lead's observable closure comparisons.

## Parameters, units and exact coordinates

`cases.json` owns ordered components, complete unlike interactions, packet
hashes and constants. `PureRecord` receives molar mass in g/mol (the shared
table is kg/mol), sigma in Å and energies divided by Boltzmann's constant
in K. Returned `model_serializations` retains `Parameters.to_json_str()`
contents for each constructed selection. At version 0.10.1 each associating
species uses `association_sites=[{kappa_ab,epsilon_k_ab,na:1,nb:1}]`;
obsolete top-level association fields are not used
([Python records](https://github.com/feos-org/feos/blob/v0.10.1/py-feos/src/parameter/mod.rs)).
Cross-volume mixing is geometric and cross-energy mixing arithmetic
([combining rule](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos/src/pcsaft/parameters.rs#L163));
the native strength includes `(sigma_i sigma_j)^(3/2)`
([strength](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos/src/pcsaft/eos/mod.rs#L182)).

C1 exact `Property.vapor_pressure_derivatives` uses
`EquationOfStateAD.PcSaftNonAssoc`, canonical array `[m,sigma,epsilon_k,mu]`,
and seeds those four names. First three columns map to shared
`(m,sigma,epsilon)` in `(Pa,Pa/Å,Pa/K)`; `mu` stays in `gradient_extra`.
C2 `PcSaftFull` uses component 1's
`[m,sigma,epsilon_k,mu,kappa_ab,epsilon_k_ab,na,nb]`, then component 2's
same eight values, then `k_ij`. Selected derivative names are shared-order
`[k_ij,epsilon_k_ab1,kappa_ab1]` in `(Pa,Pa/K,Pa)`. One self-association
parameter owns both donor–acceptor directions. Order follows
[ParametersAD](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos-core/src/ad/mod.rs#L17),
[pure build](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos/src/pcsaft/eos/pcsaft_pure.rs#L180)
and [binary build](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos/src/pcsaft/eos/pcsaft_binary.rs#L43).
Pressure and gradient return together; no finite differences run. Native
property failure produces a boolean, NaN pressure and zero-filled gradient
storage; that storage is not an available derivative.

States return SI molar density and pressure, dimensionless `ar=A_res/(nRT)`,
`Z`, `lnphi` and `mu_res=mu_res/(RT)`. Native molar energies and chemical
potentials are divided by installed `RGAS*T`; the runner records quantity
constants `NA`, `kB`, `R`. Hard sphere and chain shares are summed once to
match Engine hard-chain. C4 liquids sort by water fraction; beta maps the
native mole-based `vapor_phase_fraction` to the water-rich phase
([fraction getter](https://github.com/feos-org/feos/blob/v0.10.1/py-feos/src/phase_equilibria.rs#L479)).

## Electrolyte mapping and timing

C3 uses revised ePC-SAFT, Born off, the Held base water sigma activating
native temperature-dependent water, one donor/acceptor and every unlike
interaction. Native ion packing and screening diameters are both 0.88 sigma
([diameter](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos/src/epcsaft/parameters.rs#L365),
[Ionic call](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos/src/epcsaft/eos/mod.rs#L139)).
Sigma-screening requests are explicitly skipped. Matched rows apply the
issue's complete electrostatic factor conversion using shared constants;
literal rows keep 78.09. Equal `ExperimentalData` values at 298/299 K
enter constant permittivity. Ionic uses `EPSILON_0=8.85416e-12`,
`QE=1.602176634e-19`, `BOLTZMANN=1.380649e-23`
([constants](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos/src/epcsaft/eos/ionic.rs#L10)).
Literal `kcl=0` rows reproduce the omitted water–chloride interaction as a
distinct diagnostic selection; no packet changes.

C5 retains literal Debye inputs by default. Diagnostic requests with
`dipole_convention=matched` multiply only the entered FeOs dipole moment by
`sqrt([(1e-21/c)^2/(4*pi*epsilon0_Engine*angstrom^3)]/1e-19
*kB_FeOs/kB_Engine)`, with exact `c=299792458 m/s`, `angstrom=1e-10 m`
and the shared constants. This preserves the reduced dipole coupling;
the factor is `1.0000000002722094`, so 2.88 D becomes
approximately 2.880000000783963 D. `dipole_mu_factor` and the parsed model
records retain the actual computational inputs. Engine's conversion is in
[physical_constants.hpp](../../../engine/native/core/physical_constants.hpp)
and [dipolar_helmholtz](../../../engine/native/contributions/helmholtz.hpp);
FeOs uses `mu^2*1e-19/kB` in
[parameters.rs, lines 236–243](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos/src/pcsaft/parameters.rs#L236).
This small conversion correction does not change the polynomial coefficients.
Engine retains the printed Gross–Vrabec coefficients, while FeOs carries
additional digits in
[polar.rs, lines 11–33](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos/src/pcsaft/eos/polar.rs#L11).
Literal and mapped native values remain distinct; the lead's source-digit
diagnostic owns any coefficient correction and the raw misses remain misses.

The default zero-ion reference and pure solvent have identical T/P and `x=[1,0,0]`;
one native state supplies both. `C3-diagnosis-finite` changes only the ionic
reference to 1e-12 mol/kg; its pure-water reference remains zero ion. Reference requests evaluate each entered
finite-ladder point independently. Activity rows return Julia helper fields
`lnphi`, `reference_lnphi`, `pure_lnphi`, `m`, `MW` (kg/mol), `xwater`, plus
native reference inputs/quantities. Python performs no activity/osmotic
contraction. The [NaCl dilute-cancellation diagnosis](https://github.com/tannerpolley/ePC-SAFT/issues/222#issuecomment-5946042190)
remains a returned-value comparison target; no method or tolerance changes.

An integer third command argument enables only requests with a nonempty
`timing` field. Otherwise the runner collects values. Model construction
and imports precede the first call. Warm runs contain five unrecorded
warm-ups and 25 samples. T1 includes native State construction and
`ar`, `Z`, `lnphi` extraction. Neutral T2 includes complete native solve and
phase pressure/composition/density extraction, plus C4 mole fraction.
C3 T2 includes fresh solution/reference calls and extraction, JSON write/read
and resident compiled Julia contraction in every sample. Helper startup is
outside warm work; its first request is inside cold first-call work.
Each row retains `transport_ns` within its complete elapsed time. Optional
retained-reference timing is omitted.

C1 T3 includes native exact pressure/gradient call and extraction. C2 T3
first runs native automatic bubble with entered controls for its own pressure
guess, then calls the exact property routine. Both solves and array preparation
are inside every sample. No Engine state seeds a comparator. Full shares,
chemical potentials, slopes and nonfinite diagnostics are extracted outside
timing. Nonfinite fields retain state and an absence reason, never substituted
zero. Negative-Z fixed states retain energy, chemical potential and slopes;
`lnphi` is outside its real-valued domain because the native expression uses
`log(Z)` ([native producer](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos-core/src/state/residual_properties.rs#L348)).
The lead owns qualification, CSV retention and three serial fresh-process
repetitions; these are complete case-output costs, not equal-algorithm costs.

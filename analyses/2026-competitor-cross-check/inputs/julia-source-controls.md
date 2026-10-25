# Julia comparator conventions and calls

The comparison uses native Julia 1.13.1, Clapeyron 0.6.29, JSON 1.10.0 and
ForwardDiff 1.4.6. `Project.toml` declares exact direct versions;
`Manifest.toml` retains dependency versions and source tree hashes. The cache
is `~/.cache/epcsaft-comparators/222/`, with Julia in `julia-1.13.1/` and packages
in `depot/`. The archive is
[julia-1.13.1-linux-x86_64.tar.gz](https://julialang-s3.julialang.org/bin/linux/x64/1.13/julia-1.13.1-linux-x86_64.tar.gz),
SHA-256 `0f2e18c8dea60a2c8711d089cd9612f7a8df394b89e5208096fe5867647e0908`,
checked against Julia's [version index](https://julialang-s3.julialang.org/bin/versions.json).

Run from the checkout root. Replace the two paths at the end for a selected
input file and output file. Omit a repetition argument for qualification;
adding `1`, `2` or `3` enables the selected warm timing blocks.

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 \
  JULIA_DEPOT_PATH="$HOME/.cache/epcsaft-comparators/222/depot" \
  timeout -k 5 1800 nice -n 10 taskset -c 0 \
  systemd-run --user --scope -p MemoryMax=6G \
  "$HOME/.cache/epcsaft-comparators/222/julia-1.13.1/bin/julia" \
  --startup-file=no --project=analyses/2026-competitor-cross-check/inputs \
  analyses/2026-competitor-cross-check/scripts/clapeyron.jl \
  analyses/2026-competitor-cross-check/inputs/cases.json build/competitor-cross-check/Clapeyron.json
```

Julia first compilation uses the lead-amended 1800 s cap; the original 300 s
cap was exhausted during C2 inference/compilation. Normal optimized qualification
and three fresh-process timing repetitions subsequently completed. The earlier
`--compile=min` probe is not benchmark evidence. Post-compilation case cap remains
180 s and total warm campaign 1200 s; keep cold and warm costs separate.
Keep one Julia task, one BLAS thread and one precompile task, with CPU 0 affinity
and at most two computational processes. The retired global flock is not used.
Temporary native outputs live under `build/competitor-cross-check/`; retained
measurements and host identity are in the analysis results.

## Bounded C4 method diagnosis

`--association-probe INPUT NATIVE_CLAPEYRON OUTPUT` changes only documented
association options outside qualification/timing: baseline, `implicit_ad=true`,
and tighter atol/rtol=1e-30, each at Float64 and 256-bit precision. It retains
native state quantities and actual association-matrix/compressed dimensions.
The [compressed initializer](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.29/src/models/SAFT/association.jl#L490)
uses a hardcoded 1e-12 threshold; initializer success returns before the general
solver's option-controlled route at lines 584–588. These probes do not replace
primary results, change physical inputs, or claim that higher precision fixes it.
A separate public `volume(...;vol0=...)` call starts from each library-returned
C4 liquid volume to test pressure closure; it is a diagnostic reseed workload,
not an automatic-initialization benchmark or a runtime workaround.

## State quantities and explicit parameters

`cases.json` owns inputs, sources, units, component order and requested points.
All library state calls use one mole, molar volume `1/rho` in m³/mol, T in K,
rho in mol/m³ and pressure in Pa. `a_res` returns dimensionless Aʳ/(nRT).
`pressure` supplies Z=P/(rho RT), and `VT_chemical_potential_res` supplies
ln(phi)=muʳ/(RT)-ln(Z). Values at nonpositive Z retain muʳ/(RT), Aʳ/(nRT),
pressure and Z, with the real-valued ln(phi) explicitly absent.
The producer is [VT.jl](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.29/src/methods/VT.jl#L516).
`p∂p∂V` supplies the density slope by the exact chain dP/drho=−(dP/dV)/rho².

C1, C2 and C4 use `PCSAFT`; C5 uses `PCPSAFT`, the
[Gross–Vrabec dipole model](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.29/src/models/SAFT/PCSAFT/variants/PCPSAFT.jl).
The input nmu=1 is the selected case; its dipole in debye is passed explicitly.
Molar masses are converted from kg/mol to the library's g/mol parameter.
No default physical parameter is selected by an application name.

The association topology uses one `H` donor and one `e` acceptor for each
associating component, with explicit self H–e/e–H and unlike H–e/e–H entries.
The cross energy is arithmetic and cross kappa includes the
Wolbach–Sandler size factor specified in #222. Both entries are supplied with
`AssocOptions(combining=:nocombining,atol=1e-13,rtol=1e-13,max_iters=1000)`;
the remaining defaults are dampingfactor=0.5 and implicit_ad=false.
[PCSAFT.jl, Delta](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.29/src/models/SAFT/PCSAFT/PCSAFT.jl#L224)
multiplies dimensionless `bondvol` by sigma_ij³. Its docstring's m³ label is
not the entered basis. The returned association parameter check compares
sigma_ij³ kappa_ij with sqrt(kappa_ii kappa_jj)(sigma_i sigma_j)^(3/2);
it checks construction only. Native `a_assoc` at the two C4 fixed states is
the separate equation falsifier. No native Delta or site-residual value is
invented from a library's unavailable accessor.

C3 explicitly constructs `pharmaPCSAFT` then applies the like-charge dispersion
suppression in the pinned
[ePCSAFT constructor](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.29/src/models/SAFT/PCSAFT/variants/ePCSAFT.jl#L79).
Mixed unlike energies are constructed before ion diagonal epsilon is zeroed;
the latter also produces the native 0.88 sigma ion packing diameter.
The native [water correlation](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.29/src/models/SAFT/PCSAFT/variants/PharmaPCSAFT.jl#L109)
uses the supplied base diameter. The combined `ePCSAFT` object uses `DH` or
`hsdDH` and `ConstRSP`; Born is absent. This avoids the convenience constructor's
CSV-only parameter-path argument while retaining its numerical model.
`iondata`, the returned screening diameters, evaluated water sigma, Bjerrum
length and kappa*d expose the actually consumed convention.
Matched relative permittivity compensates both charge and vacuum permittivity
differences exactly as declared in #222; literal 78.09 remains separate.
The constants come from [constants.jl](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.29/src/base/constants.jl).

## Solvers and complete timed blocks

Saturation uses `ChemPotVSaturation(atol=1e-12,rtol=1e-12,max_iters=1000)`;
bubble uses `ChemPotBubblePressure` with the same controls. LLE uses
`MichelsenTPFlash(equilibrium=:lle,K_tol=1e-12,ss_iters=1000,nacc=5,second_order=true)`
and `tp_flash2`. Automatic initialization runs inside every solve. Liquids
are ordered by water mole fraction; beta is the water-rich mole fraction
of total phase amount. No returned split establishes global stability.

The native [volume solver](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.29/src/methods/property_solvers/volume.jl#L15)
has max_iters=100, log-volume convergence `(atol=0,rtol=1e-12)` at line 131,
and the early pressure test `abs(p-p0)<3eps(p0)` at line 81. Its mechanical
and spinodal checks remain unchanged. Library iterations, association residual
and solver reduced residual are absent from the returned values used here;
observable pressure and chemical equalities are computed from native states
by the lead. A reported solver success does not substitute for those residuals.

Automatic LLE initialization uses `K0_lle_init`, then
`tpd2(lle=true,strategy=:pure,break_first=true)`. The tangent-plane solver's
source defaults are tol_trivial=1e-5, tol_equil=1e-10 and it_ss=30;
its minimizer uses f_abstol=1e-12, f_reltol=1e-10 and maxiter=100.
See [tpd.jl](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.29/src/methods/property_solvers/stability/tpd.jl#L165).
The Gibbs minimizer after Michelsen successive substitution uses
f_abstol=0, f_reltol=0, g_reltol=K_tol, g_abstol=1e-8 and maxiter=100,
with the Newton2 bounded line search when second_order=true;
see [Michelsentp_flash.jl](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.29/src/methods/property_solvers/multicomponent/tp_flash/Michelsentp_flash.jl#L447).

T1 calls `a_res`, `pressure`, and `VT_chemical_potential_res`, then extracts
Aʳ/(nRT), Z and each ln(phi); it includes residual chemical-potential automatic
differentiation. It excludes model construction, density solves, contributions,
density slopes and source diagnostics. T2 saturation/bubble includes the
complete native solver and pressure, density and composition extraction.
T2 LLE also includes automatic tangent-plane initialization, phase sorting and
beta extraction. T2 C3 calls BOTH `mean_ionic_activity_coefficient` and
`osmotic_coefficient`, extracts the two coefficients, and includes their repeated
solution/reference solves. These helpers use the 1e-20 mol/kg reference
in [mT.jl](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.29/src/methods/property_solvers/electrolytes/mT.jl#L1).
Reference and pure-solvent states are collected again outside timing for
qualification. Adopted D6 excludes Clapeyron gradient timings and makes no
claim about gradient availability.

Without a repetition argument only cold import, model construction and first
case calls are recorded. With a repetition argument, the script imports and
constructs models first, qualifies initial outputs outside timing, then runs
five warmups and 25 samples for requests marked T1 or T2. Each sample measures
the complete block above with `time_ns`; the lead serializes fresh-process
repetitions and statistics. Case-output costs differ in internal work across
libraries and are not equal-algorithm benchmarks.

## Resident FeOs contraction and source diagnostics

The `--contract` command in `clapeyron.jl` imports JSON and LinearAlgebra, sets BLAS
threads to one, reads one JSON object per stdin line, and emits/flushed one JSON
object per stdout line. It does not import Clapeyron or ForwardDiff. Input keys
are `lnphi`, `reference_lnphi`, `pure_lnphi`, `m`, `MW`, `xwater`; `pure_lnphi`
may be a scalar water value or a vector with water first. Returned keys are
`gamma_x`, `gamma`, `phi`; the ionic finite/zero reference and pure water are
separate inputs. The equations are full-eos §10, applied to native FeOs values,
not a FeOs activity-helper capability. Time native states, input/output
serialization and helper transport completely; transport-only elapsed time
is separately retained by the caller.

`--diagnose INPUT OUTPUT` evaluates C3 source errors outside timing.
The [DH rational expression](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/v0.6.29/src/models/Electrolytes/Ion/DH.jl#L80)
is evaluated with its actual Float64 coefficients lifted into BigFloat.
The analytic expression uses the alternating series below y=0.5, stopping
when the next term is below 1e-60, and log1p above it. ForwardDiff returns
first/second y derivatives and propagates source-minus-analytic errors through
the dimensionless Helmholtz energy, pressure and chemical potentials. Density
refinement moves toward the analytic pressure root with the sign declared in
#222, and retains two successive updates separately from the raw result.
Each actual screening y and both sides of 1.5 are retained. The 256/512-bit
checks precede Float64 transport of property rows; string function/density
values retain high-precision digits. Source-error diagnostics quantify a
specific comparator expression, not Engine properties or physical accuracy.

The diagnostic first resolves the unmodified library pressure root with two
high-precision updates, separately recording its shift from the supplied
Float64 density. Source-error displacement then starts at that root. Fixed
requests keep their supplied density and have only the fixed-state correction.
Mole fractions lifted into BigFloat are normalized to one mole. Activity rows
retain corrected gamma/Phi, signed delta ln(gamma)/delta Phi, and the native
helper minus the refined unmodified-library contraction separately.

`--screening` accepts JSON-lines with `id`, `T`, `rho`, `x`, `diameters` (m),
`epsilon_r` and `constants` (`e`, `epsilon0`, `kB`, `NA`). It emits one row
at each of 256/512 bits with Bjerrum length, kappa*d and analytic/rational
function values and derivatives at the supplied native state. The extra
Float64 source-expression values use the exact subtraction order in
[FeOs ionic.rs](https://github.com/feos-org/feos/blob/v0.10.1/crates/feos/src/epcsaft/eos/ionic.rs#L67);
they diagnose expression cancellation only, not returned FeOs fugacities.
The zero-ion expression is absent because the library takes its explicit zero
ionic-energy branch. The producer is an independent source calculation;
these are not labelled as a native FeOs kappa accessor.

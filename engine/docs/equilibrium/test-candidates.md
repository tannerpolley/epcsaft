# Potential equilibrium numerical tests

Status: candidates, not implemented or accepted tests. Researched 2026-09-07.
The owner requested both matching-model numerical references and comparisons
against other equilibrium models with justified allowances for EOS differences.
No upstream calculation was rerun during this inventory. Numbers below are
upstream assertions or displayed examples; rounded examples need reproduction
before use. This document owns candidate selection only; the current completion
sequence remains in the GitHub issue roadmap. PH/PS/UV and
phase-discovery candidates are retained for future implementation and do not
add current checks or acceptance gates. Reuse applicable declared-phase cases
for the bounded gap before inventing another test campaign.

## Focused falsifiers retained before producer corrections

These two native checks protect the identified producer corrections. They test independent differential identities. F1 passes after the 2026-09-08
caloric correction; F5 passes and runs in `Automatic_fast`. A failure is diagnostic evidence, not a qualified
Engine result. The native checks print their ladder values and maximum errors to
ordinary process output. Run `properties_checks F1` and `equilibrium_checks F5`
from the configured build directory recorded in
`build/environment-wheel/meson-build-dir`. The retained pre-cutover tables
preserve source, wheel, executable, and physical-input provenance; they are
historical evidence, not a current-code pass.

### F1 residual internal energy

`engine/tests/native/properties_p1_m3_checks.cpp` evaluates the Gross 2001 revision-2
ethane packet at `T=300 K`, `rho=10000 mol/m^3`, `x=(1)`. At fixed density and
composition it compares the state field with the independent relation

$$u^{\mathrm{res}}=-RT^2\left(\frac{\partial a^{\mathrm{res}}}{\partial T}\right)_{\rho,\mathbf{x}}=h^{\mathrm{res}}-RT(Z-1).$$

The centered temperature ladder uses `Delta T = (0.1, 0.03, 0.01, 0.003) K`
and evaluates `a^r` through separate state calls. The `2e-7` relative ladder
allowance is provisional and local to this retained state: it is just above
the largest relative discrepancy in the retained historical refinement rows,
with a roundoff margin. The `1e-10` scaled identity allowance is a provisional
same-state detection bar for the direct relation. Neither number is a
cross-model or universal caloric accuracy limit. The historical wheel probe
and its rows are stale evidence; the focused current native run is recorded
separately with its wheel and executable identity.

### F5 weighted Hessian and gradient

`engine/tests/native/equilibrium_s6_e2_checks.cpp` differentiates the full weighted
gradient independently of the callback Hessian for both a named `regular`
two-phase neutral ternary problem and a named `pinned`/free-intensive neutral
binary boundary problem. For every full primal coordinate, including embedded
chemical multipliers and free intensive coordinates where present, it checks

$$\frac{\partial}{\partial\mathbf y}\left[\alpha\nabla f(\mathbf y)+\sum_j\lambda_j\nabla g_j(\mathbf y)\right]$$

against the lower triangle returned by the Ipopt callback. The cases use
nonzero callback multipliers and `alpha=(0, 0.5, 1)`, so objective weighting
and constraint curvature cannot cancel. The callback row ordering follows
the same active-row selection as the driver; optimizer callback multipliers
remain distinct from chemical balance multipliers.

Each objective factor is checked on a step-size ladder of
`(10, 3, 1)*sqrt(epsilon)*max(1,|y|)`, and the check prints the case,
factor, step factor and largest error/scale. The regular case keeps the
existing `1e-6` scaled comparison; the pinned case uses a provisional `1e-5`
allowance because its full coordinate vector includes embedded multipliers.
These are local finite-difference detection allowances, carry no model-
adequacy claim and are not a blanket tolerance. The current source is
expected to fail from the pinned zero-Hessian path and the objective-factor
treatment; retain this failure before any production correction. The retained
CSV records the zero-based lower-triangle row and column for each largest-error
entry, while every checked entry remains subject to the local comparison.

## What a comparison establishes

- Matching equations, parameters, mixing/association rules and references can
  support tight independent-implementation comparisons.
- Different EOS/activity models can support physical similarity comparisons:
  phase compositions/fractions, boundary locations, trends and energy responses.
  Their discrepancy includes model differences, not just numerical error.
- Each Engine solve must independently meet its requested conservation,
  equilibrium and energy residual criteria. Cross-model allowances never loosen
  those criteria. Matching another package alone does not prove correctness.
- Select cross-model tolerances per observable/domain before acceptance. Use
  available measurements, independent model spread and application requirements;
  do not choose a tolerance merely large enough to pass the current result.
  Use absolute composition/fraction tolerances near zero. Compare phase sets
  without depending on phase ordering. Different phase counts near boundaries
  require diagnosis or boundary-curve comparison, not automatic failure or
  an arbitrarily widened tolerance.
- The default replay declares phases. Stability/discovery is a separate explicit
  mode. Borrow inputs and scientific comparisons, not upstream solver routes.

## External execution status

Not run here. On 2026-09-08 the owner selected running these cases with their
scientific slices. Use pinned package installations first; source clones are
needed only when an upstream test or source modification requires them.
See [Clapeyron installation](https://clapeyronthermo.github.io/Clapeyron.jl/dev/installation/)
and [FEOS installation](https://feos-org.github.io/feos/installation.html).
Retain actual outputs and package/version/input identities with the slice;
upstream assertions below must never appear as locally computed results.

## Sources and parameter availability

Pinned inspected source trees:

- [Clapeyron](https://github.com/ClapeyronThermo/Clapeyron.jl/tree/b2a2027def8a03d2107567d4a1e4e38e64675ee5), abbreviated C below.
- [FeOs](https://github.com/feos-org/feos/tree/f1fd55d8a3cd7a254cd8721619e01820a453234e), abbreviated F below.

[C PC-SAFT tables](https://github.com/ClapeyronThermo/Clapeyron.jl/tree/b2a2027def8a03d2107567d4a1e4e38e64675ee5/database/SAFT/PCSAFT)
contain pure parameters for methane, ethane, propane, butane, nitrogen,
cyclohexane, hexane, octane, ethanol, methanol, water and 1-butanol.
Association records are supplied for the associating species. Methanol/cyclohexane
has an explicit interaction value 0.051. Table source locators include
Gross–Sadowski 2001 (10.1021/ie0003887) and 2002 (10.1021/ie010954d).
Hydrogen and oxygen pure records also exist, with distinct cited sources.
The local `data/packets/gross-sadowski-2001/` already contains many hydrocarbons,
nitrogen and carbon monoxide; this is not an exact parameter-match claim.

[F test parameters](https://github.com/feos-org/feos/blob/f1fd55d8a3cd7a254cd8721619e01820a453234e/crates/feos/tests/pcsaft/test_parameters.json)
include propane, butane, hexane and associating nonpolar-model water (`water_np`).
Some numerical values differ from C even for the same species.

Pure-parameter availability is confirmed, but a complete mixture packet still
needs exact binary interactions (including explicit intended zeros), association
rules/sites and enabled contributions. Do not reuse PR/SRK interaction constants
as PC-SAFT interaction constants. Source-range applicability stays with the
caller/study; the Engine does not acquire or police parameter literature.

## Previously identified matching-model candidates

These can become matching-model references only after compatibility checks.
All pressures below are Pa and temperatures K; compositions are mole fractions.

[C equilibrium tests](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/b2a2027def8a03d2107567d4a1e4e38e64675ee5/test/methods/eos_pcsaft.jl):

| Candidate | Inputs | Stored expected quantity |
|---|---|---|
| Ethanol saturation | T=298.15 | P=7972.550405922014 |
| Ethanol boiling | P=100000 | T=351.32529505096164 |
| Methanol/cyclohexane bubble | T=313.15, x=(.5,.5) | P=54532.249600937736 |
| Same, bubble temperature | P=2000000, x=(.5,.5) | T=435.80890506865 |
| Same, dew pressure | T=443.15, y=(.5,.5) | P=1655548.6543884084 |
| Same, dew temperature | P=2000000, y=(.5,.5) | T=453.0056727580934 |
| Same, LLE pressure | T=313.15, one liquid x=(.27,.73) | P=737971.7522006684 |
| Same, LLE temperature | P=100000, one liquid x=(.27,.73) | T=312.9523684945214 |
| Same, azeotrope | T=443.15 / P=100000 | P=2443546.2800998255 / T=328.2431049077264 |
| Same, VLLE | T=313.15 | P=54504.079665621306 |

These LLE compositions specify one coexisting liquid, not an overall feed.
Azeotrope and VLLE inputs require their own well-posed shared-row specifications.
Critical/spinodal assertions in the source are adjacent derivative candidates,
not additions to the current equilibrium scope by this document.

[C multiphase benchmark](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/b2a2027def8a03d2107567d4a1e4e38e64675ee5/test/tpflash/tp_flash_benchmark_cases.jl):

| Candidate | Inputs | Available evidence |
|---|---|---|
| Water/cyclohexane/propane VLLE | T=298.15, P=100000, z=(.333,.333,.334) | Three complete phase compositions and fractions stored upstream |
| Same species, propane absent, LLE | T=298.15, P=100000, z=(.5,.5,0) | Two complete compositions/fractions, exactly zero propane |

Retain the full upstream precision once in a numerical table when adopted;
do not duplicate it here. Compare densities and residuals as well as splits.

[F TP tests](https://github.com/feos-org/feos/blob/f1fd55d8a3cd7a254cd8721619e01820a453234e/crates/feos/tests/pcsaft/tp_flash.rs):
propane/butane at 250 K with pressure/feed constructed from pure saturation;
propane/butane/hexane at 300 K, 120000 Pa, z=(0,.5,.5), compared with the
reduced binary. These assert equilibrium identities/support reduction; a pinned
run is needed for complete reference rows.

[F pure tests](https://github.com/feos-org/feos/blob/f1fd55d8a3cd7a254cd8721619e01820a453234e/crates/feos/tests/pcsaft/vle_pure.rs):
propane at 170, 200, 250, 300, 350 K and separately 0.1, 1, 10, 30, 44 bar.

[F stability tests](https://github.com/feos-org/feos/blob/f1fd55d8a3cd7a254cd8721619e01820a453234e/crates/feos/tests/pcsaft/stability_analysis.rs):
water_np/hexane equimolar liquid at 300 K, 1 bar requires an instability
witness; propane/butane bubble phases at 300 K are searched without a witness.
No-witness outcomes are not global certificates.

[F thermal flashes](https://github.com/feos-org/feos/blob/f1fd55d8a3cd7a254cd8721619e01820a453234e/crates/feos/tests/pcsaft/px_flashes.rs):
manufactured binary PH/PS cases at 50 bar, composition .3, explicit residual and
Joback parameters, including a kij temperature sensitivity. Check every enabled
contribution before adoption, particularly polar terms excluded by our scope.
Parameter availability alone does not establish compatibility.

[C electrolyte LLE](https://github.com/ClapeyronThermo/Clapeyron.jl/blob/b2a2027def8a03d2107567d4a1e4e38e64675ee5/test/methods/electrolytes.jl):
water/ethanol/NaCl at 298.15 K, 1 bar, solvent ratio (.5,.5), salt molality 6.
Stored compositions use SAFTVREMie, so this is a cross-model candidate. Complete
ePC-SAFT ionic/solvent interactions and dielectric choices remain unresolved.

## Additional cross-model candidates

### Thermo: PR, SRK and IAPWS models

[Official examples](https://thermo.readthedocs.io/thermo.flash.html), inspected
as documentation version 0.6.1; pin the executable version/data before replay.

| Candidate | Source inputs and displayed result | PC-SAFT preparation |
|---|---|---|
| Methane/ethane/nitrogen PR flash | 110 K, 1 bar, z=(.965,.018,.017); vapor fraction about .0890 | Pure records available; select PC-SAFT interactions separately |
| Same mixture dew/bubble | Dew T about 133.8 K at 1 bar; bubble P about 515029.6 Pa at 133 K | Compare boundary curves and compositions |
| Same mixture PH/PS inversion | Reconstruct 110 K from its own TP enthalpy/entropy | Match thermal references for direct cross-model comparison; also retain each model's internal inversion |
| Butanol/water/ethanol SRK VLLE | 361 K, 1 bar, z=(.25,.7,.05); source reports three phases | Pure/association records available; interactions unresolved; full splits need replay |
| Pure water IAPWS saturation | 600 K, vapor fraction .5; P=12344824.3572 Pa | Water parameter choice available; expect meaningful model discrepancy |
| Pure water UV | V=.00061 m3/mol, U=55850 J/mol; T about 800.5922 K, P about 10144789 Pa | Requires consistent reference calorics before comparison |

### Phasepy: association-sensitive bubble point with a different model

[Official example](https://phasepy.readthedocs.io/en/latest/): ethanol/water,
liquid composition (.5,.5), 320 K, Peng–Robinson with modified Huron–Vidal and
Dortmund UNIFAC. Displayed vapor composition (.70761727,.29238273), pressure
0.23248584919691206 bar. PC-SAFT pure and association records exist. Resolve
mixture parameters separately. Useful for pressure/composition comparisons and
a later composition sweep; do not use its PR inputs as PC-SAFT inputs.

### BioSTEAM/Thermosteam: phase fractions, LLE and VLLE

[Official examples](https://biosteam.readthedocs.io/en/latest/tutorial/Phase_equilibrium.html):

| Candidate | Inputs / displayed output | Use and limitation |
|---|---|---|
| Water/ethanol VLE | Equal feed amounts; 355 K, 101325 Pa | Stored phase flows; normalize explicitly |
| Same, specified vapor fraction | VF=.5, 101325 Pa gives 353.94 K; VF=.5, 360 K gives 127822 Pa | General free-T/free-P constraints |
| Same, PH pressure change | From the 355 K state, double P at fixed H: displayed T=373.69 K | Reference alignment required |
| Water/1-butanol/octane LLE | Feed amounts (304,30,100), 360 K; displayed P=101325 Pa | Liquid partitioning; pure/association records available |
| Water/octane/ethanol VLLE | Feed (.5,.5,.1), 355 K, 101325 Pa | Complete displayed phase flows; candidate only |
| Same VLLE at VF=.5 | P=101325 Pa; displayed T=358.76 K | Difficult phase-fraction specification |

Pin the resolved property package; do not assume the tutorial's defaults are
immutable. The documentation explicitly calls VLLE experimental and potentially
inaccurate/inconsistent: these rows require an independent corroborating run,
not direct promotion to authoritative expected values. Some other tutorial
cells display conditions inconsistent with their calls; do not adopt them.

### Clapeyron: additional PR/UNIFAC mixture-space cases

The pinned C multiphase benchmark above also provides PR four-hydrocarbon
near-boundary and trace-component cases. Candidates include equimolar
isobutane/n-butane/n-pentane/n-hexane at 1 bar and 282.2, 282.3 and 284.4 K,
and n-butane/n-pentane/n-hexane/n-heptane with strongly unequal trace feeds.
PC-SAFT pure records exist. Recompute phase boundaries: a tiny PR vapor fraction
need not remain tiny or nonzero under PC-SAFT. The UNIFAC water/hexane LLE
case at 303.15 K, 101325 Pa, equimolar feed is another cross-model candidate.

### Cantera: restricted reactive gas equilibrium

[Official methane/air example](https://cantera.org/3.2/userguide/flame-temperature.html)
includes a restricted five-species ideal-gas calculation (CH4, O2, N2, CO2,
H2O) and a larger mechanism. Its HP sweep starts at 300 K and 1 atm,
with equivalence ratio 0.5–2.0 and oxidizer O2:N2=1:3.76.
Pure PC-SAFT records exist for the five-species
subset; use identical species/element balances and compatible standard-state
thermochemistry. This is a candidate for TP reaction equilibrium followed by
HP equilibrium, not a kinetics simulation. Exact reference outputs need replay.
Do not compare a restricted Engine species set to the full mechanism's result.
Radicals in the full mechanism have not been checked for PC-SAFT parameters.
High-temperature model applicability must be assessed in the study, even though
runtime range policing is not the Engine's job. Defer adoption until reference
thermochemistry and reaction inputs are complete.

## Adoption order and retained evidence

Owner selection, 2026-09-08: prepare a small diverse first set comprising
ethanol saturation, methanol/cyclohexane bubble/dew, the zero-component flash,
declared VLLE, Thermo's PR natural-gas and Phasepy's ethanol/water cases.
Use case-specific evidence for cross-EOS limits, independently of numerical
residual tolerances. This selects preparation work, not accepted test values.

1. Establish matching-model ethanol and hydrocarbon cases, then bubble/dew and
   declared LLE/VLLE. These isolate numerical disagreement most directly.
2. Add Thermo's natural-gas and Phasepy's ethanol/water comparisons, with
   independent PC-SAFT inputs and agreed model-discrepancy limits.
3. Add liquid partitioning and phase-fraction cases after parameter review;
   thermal and reactive cases wait for explicit reference inputs.
4. For each adopted case, retain package version/revision, resolved parameters,
   species order, units, phases, initial conditions, physical specifications,
   returned states/failures, residuals and uncertainty/tolerance rationale.
   Review source licensing/attribution when copying code or parameter tables.
5. Put long upstream replays and comparisons under one-question analyses using
   the installed Engine wheel. Promote only compact distinct numerical
   falsifiers to persistent Engine tests; do not require external packages in
   the ordinary Engine test process or reproduce their solver implementations.

No tolerance numbers are approved by this inventory. The owner accepts
cross-EOS similarity as a goal; numerical and model-error allowances remain
separate, explicit decisions supported by the adopted cases.

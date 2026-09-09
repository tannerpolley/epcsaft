# Preparation doctrine record: equilibrium driver, thermal forms and reaction inputs

Date: 2026-09-08

Status: preparation record only. This file records source-backed equation and
convention proposals. It does not qualify a scientific sublayer, repair
production numerics or qualify the owner-selected DIPPR107 form.

## Scope and authority

The governing project records are the [scientific context](../CONTEXT.md), the
[implementation plan](../implementation-plan.md), ADR
[0009](../adr/0009-general-platform-one-driver.md). The canonical equation owner is
[engine/docs/equations.md](../../../engine/docs/equations.md), and the algorithm
owner is [engine/docs/science/algorithms.md](../../../engine/docs/science/algorithms.md).
The generated equation index is
[engine/docs/equations_registry.yaml](../../../engine/docs/equations_registry.yaml).

The changes in this preparation slice are limited to those two Engine
documents, the generated index, and this record. No implementation-plan row,
production source, parameter bundle or application source was changed.

## Equilibrium driver doctrine

Every admitted equilibrium specification uses one Ipopt driver and the shared
thermodynamic energy. The driver must supply the exact weighted Lagrangian
curvature required by EqID [`weighted_lagrangian_hessian`](../../../engine/docs/equations.md#weighted_lagrangian_hessian)
at arbitrary callback points, objective factors and optimizer multipliers:

\[
\nabla_y^2\mathcal L_{\mathrm{opt}}
=\mathrm{obj\_factor}\,\nabla_y^2 f(y)
+\sum_j\Lambda_j^{\mathrm{opt}}\nabla_y^2g_j(y).
\]

The objective factor scales objective curvature only. Pinned nonlinear pressure
and potential rows retain their underlying Hessians when the objective is
constant. A zero multiplier makes the corresponding weighted term zero; it does
not remove the requirement to implement the constraint Hessian. This applies to pinned/free-intensive and
incipient bubble/dew specifications. Root and limited-memory equilibrium
replacements are outside the accepted driver doctrine. A specification that
cannot be represented by the shared exact-curvature formulation is unavailable
with a named diagnostic; it does not acquire a second driver.

The algorithm record now states this requirement and retains the derivative
obligation. Constraint curvature generally needs third derivatives of the
energy, and lifting variables does not remove that requirement automatically.
The wording is a capability requirement, not evidence that the current
implementation already provides the callback or its third derivatives. A
weighted-Hessian falsifier and implementation qualification remain later work
under the plan.

## NASA7 and NASA9 thermal equations

The selected initial thermal forms are polynomial heat capacity (including the
constant case), Shomate, NASA7, NASA9 and DIPPR107. This slice registers NASA7 and NASA9
as bounded empirical relations. Their canonical equations are in
[`cp_ideal_nasa7`](../../../engine/docs/equations.md#cp_ideal_nasa7) and
[`cp_ideal_nasa9`](../../../engine/docs/equations.md#cp_ideal_nasa9).

The primary source locator for the old seven-coefficient and new nine-constant
forms is NASA RP-1311 Part I (1994), section 4.2, Eqs. (4.6)-(4.11), available
through the [official NTRS record](https://ntrs.nasa.gov/citations/19950013764).
The NASA Glenn data report is McBride, Zehe and Gordon, NASA/TP-2002-211556
(2002), under “Empirical Equations for Fitting Thermodynamic Functions,” Eqs.
(1)-(3), [official NTRS record](https://ntrs.nasa.gov/citations/20020085330).
The NTRS record describes the data as least-squares coefficients for a
seven-term heat-capacity form with integrated enthalpy and entropy constants.
The [Cantera species-thermo documentation](https://cantera.org/stable/reference/thermo/species-thermo.html)
provides an accessible technical transcription: NASA7 is the original
NASA/Chemkin form, and NASA9 is the NASA Glenn form with the two inverse-power
terms and an arbitrary number of temperature regions.

The Engine transcription uses the following owner conventions.

| Convention | Prepared rule |
| --- | --- |
| Temperature | `T` is in K and is positive. The dimensionless `tau = T/(1 K)` keeps the source's numerical Kelvin convention explicit. |
| Energy units | `R` is in J mol^-1 K^-1; `h` is in J mol^-1; `Cp` and `s` are in J mol^-1 K^-1. |
| Standard pressure | Each coefficient interval carries an explicit source pressure `P_src`. Entropy on caller pressure `P` includes `-R ln(P/P_src)` for the ideal-gas pressure change. No 1 atm or 1 bar default is inserted. |
| Intervals | Each record declares an interval `I_r`. Evaluate only inside that interval. Preserve source boundary behavior; do not smooth a boundary or promise a derivative at a boundary without source continuity data. |
| Integration constants | NASA7 `a5,a6` and NASA9 `a7,a8` are integration constants in the source form. A record must state whether they provide absolute standard-state `h` and `s` or whether explicit `h_ref` and `s_ref` anchors are also required. No 298.15 K anchor is inferred from the equation kind. |
| Ideal reference assembly | [`ideal_gas_helmholtz`](../../../engine/docs/equations.md#ideal_gas_helmholtz) consumes caller-declared standard-state thermochemistry. Missing formation or entropy references make total calorics unavailable while leaving phase-only requests independent of those inputs. |

NASA RP-1311 describes simultaneous fitting of heat capacity, enthalpy and
entropy and constraints at 298.15 K or common interval points. That source
practice does not authorize an Engine-wide anchor or smoothing policy. The
recorded interval, source pressure, standard state and anchor status remain
input data. The equations satisfy the interior identities
\(\mathrm dh^\circ/\mathrm dT=C_p^\circ\) and
\(T\,\mathrm ds^\circ/\mathrm dT=C_p^\circ\); with source and caller standard pressures held fixed; these identities are algebraic
checks for a later falsifier, not numerical qualification.

## DIPPR107 owner decision

On 2026-09-08 the owner selected DIPPR107 as the fifth initial thermal form,
alongside polynomial (including constant), Shomate, NASA7 and NASA9.
`cp_ideal_dippr_107` remains in the canonical equation owner. Its implementation,
source/domain trace, derivatives and numerical qualification are required; this
scope decision alone supplies none of that evidence.

## Reaction definition trace

The current native input is an evaluated, dimensionless equilibrium constant or
an affinity target. [`ReactionSpec`](../../../engine/native/equilibrium/problem.hpp#L83)
contains a name, products-positive stoichiometry, `equilibrium_constant` and an
optional affinity. It has no native temperature-correlation object. When no
standard Gibbs vector is supplied, the native stoichiometry path maps the
positive constant to \(-RT\ln K\) in
[`stoichiometry.cpp`](../../../engine/native/equilibrium/stoichiometry.cpp#L176),
and the affinity path uses the same mapping at lines 439-452. This means a
correlation evaluator and its source standard-state metadata remain outside
the current Engine input structure. No code change is proposed here.

The read-only absorber draft contains nine species in this order:
carbon-dioxide, monoethanolamine, water, protonated-monoethanolamine,
carbamate-anion, bicarbonate-anion, carbonate-anion, hydronium-cation and
hydroxide-anion. Charges and products-positive stoichiometry are declared in
[`reaction-system.json`](../../../downstream/MEA-Absorption-Column/src/mea_absorption_column/data/epcsaft_datasets/MEA_reactive_epcsaft_bundle/reaction-system.json#L1).
The source-contract owner is
[`chemical_reaction_source_contract.json`](../../../downstream/MEA-Thermodynamics/data/reference/MEA/manifests/chemical_reaction_source_contract.json#L1).

| Reactions | Current draft candidate | Source-contract record | Read-only limit |
| --- | --- | --- | --- |
| R1-R3 | `a + b/T + c ln(T) + dT + standard_state_offset`, with offsets 8.0330699846, 4.0165349923 and 4.0165349923 | Austgen 1991, mole-fraction unsymmetric aqueous infinite-dilution source; source domain 273.15-498.15 K | The draft offsets are rounded representations of a source-to-molality conversion. Their source convention and interval still need qualification. |
| R4 | `a + b/T`, `a=3.3515778177997895`, `b=-1895.3 K`, draft domain 293.15-393.15 K | Tong 2012 via Aroua 1999 candidate, `a=2.151`, `b=-1545.3 K`, molality infinite-dilution water, source domain 293.15-323.15 K | The draft is a fitted application candidate and differs from the source-contract candidate. The 323.15-393.15 K portion of the draft domain has no source support in this record. |
| R5 | `-ln(10) (a_k/T+b+cT)`, `a_k=2597.91 K`, draft domain 293.15-393.15 K | Bates and Pinching 1951 Eq. 7, `a_k=2677.91 K`, `b=0.3869`, `c=0.0004277 K^-1`, source domain 273.15-323.15 K | The coefficient and domain differ between the draft and source contract. The draft range above 323.15 K is extrapolation and is not source-qualified. |

The source contract also records the common aqueous-molality scale, water molar
mass 0.01801528 kg mol^-1, and the source-to-common offsets in
[its standard-state block](../../../downstream/MEA-Thermodynamics/data/reference/MEA/manifests/chemical_reaction_source_contract.json#L215).
Those records are evidence of declared inputs and transformations, not a
qualification of any reaction model.

## Native reaction reparameterization

The selected native representation is

\[
\ln K(T)=a_n+\frac{b_n}{T}+c_n\ln\left(\frac{T}{T_{\mathrm{ref}}}\right)+d_nT,
\qquad T_{\mathrm{ref}}>0,
\]

with an explicit `T_ref` supplied by the caller or record. For a source form

\[
L_s(T)=a+\frac bT+c\ln(T/T_{\mathrm u})+dT+o,\qquad T_{\mathrm u}=1\,\mathrm K,
\]

where `o` is a temperature-independent source-to-common standard-state offset
and coefficients use the numerical Kelvin convention, the exact
algebraic map is

\[
a_n=a+c\ln(T_{\mathrm{ref}}/T_{\mathrm u})+o,
\qquad b_n=b,\qquad c_n=c,\qquad d_n=d.
\]

For the R1-R3 infinite-dilution conversion, products are positive and solvent
water is excluded from the solute stoichiometric sum. With standard molality
\(m^\circ=1\ \mathrm{mol\,kg^{-1}}\),

\[
\frac{m_i}{m^\circ}=\frac{x_i}{x_{\mathrm{H_2O}}M_{\mathrm{H_2O}}m^\circ}\;\longrightarrow\;\frac{x_i}{M_{\mathrm{H_2O}}m^\circ}\quad(x_{\mathrm{H_2O}}\to1),
\qquad
\ln K_m=\ln K_x-\Delta\nu_{\mathrm{solute}}
\ln(M_{\mathrm{H_2O}}m^\circ).
\]

Using the source-contract water mass, \(\ln(M_{\mathrm{H_2O}}m^\circ)
=-4.016534992299479\). Therefore the exact offsets for R1, R2 and R3 are
8.033069984598958, 4.016534992299479 and 4.016534992299479. The draft JSON
stores them rounded to 10 decimal places, so its offset representation differs
from the exact derived values by about \(10^{-12}\). The conversion is a
declared standard-state transformation, not a fitted correction.

For the arithmetic identity check only, choose \(T_{\mathrm{ref}}=313.15\) K,
which is the retained provider-transform temperature. Evaluating R1-R3 with
the exact offsets gives common \(\ln K\) values of
`-31.17540213354659`, `-14.505037067112742` and `-23.53653838436904`; expanding
the native form returns the same values to floating-point roundoff. This
check does not establish the physical standard state, source applicability,
coefficient selection or provider qualification.

In source-form labels throughout this record, `ln(T)` means
`ln(T/(1 K))`; source coefficients use the numerical Kelvin convention.
The two other transformations are direct algebraic cases. For a candidate
`a+b/T`, set \(a_n=a\), \(b_n=b\), and
\(c_n=d_n=0\). For a candidate
\(-\ln(10)(a_k/T+b+cT)\), set
\(a_n=-\ln(10)b\), \(b_n=-\ln(10)a_k\), \(c_n=0\), and
\(d_n=-\ln(10)c\), with any independently declared temperature-independent standard-state offset
folded into `a_n`. The R4 and R5 coefficient disagreements above remain
unresolved by this reparameterization.

Differentiating the native form gives

\[
\frac{\mathrm d\ln K}{\mathrm dT}
=-\frac{b_n}{T^2}+\frac{c_n}{T}+d_n,
\qquad
\Delta_r h^\circ=RT^2\frac{\mathrm d\ln K}{\mathrm dT},
\]

on the interior of the declared source interval with temperature-independent
standard pressures and a consistent standard-state definition held fixed.
Temperature-dependent standard-state conversions add their own derivative.
This derivative identity does
not supply a standard-state definition. The existing
[`reaction_constant_ln_polynomial`](../../../engine/docs/equations.md#reaction_constant_ln_polynomial)
remains the source-form owner with `c ln(T)`; the native `ln(T/T_ref)` form is
an explicit representation change with a recorded `T_ref`.

The separate [`standard_pressure_density_conversion`](../../../engine/docs/equations.md#standard_pressure_density_conversion)
handles pressure-to-density reference changes and contributes its own
temperature derivative. It must not be conflated with the R1-R3
mole-fraction-to-molality conversion. Activity, mole-fraction, molarity and
molality support still needs explicit standard quantities, solvent basis,
pressure or temperature conversion data and later numerical checks. Missing
source information remains unavailable; no standard state is guessed from the
equation shape.

## Verification performed and remaining work

The changed equation definitions were parsed and the generated registry was
regenerated. The focused registry check and source-only repository check are
the appropriate validation for this documentation preparation slice. No Engine
wheel was built, no numerical test was run and no scientific claim was
qualified.

Remaining preparation work includes the disposition of differing R4/R5 draft
candidates and the source and applicability checks
for each reaction basis. The next scientific work must compare source records,
standard-state conversions, coefficient intervals and independent numerical
identities before accepting any thermal or reaction capability.

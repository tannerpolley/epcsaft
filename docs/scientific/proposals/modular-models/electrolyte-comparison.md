# Electrolyte alternatives: equations, parameters and evidence

[Overview](README.md) · [Program issue #23](https://github.com/tannerpolley/ePC-SAFT/issues/23) · [Annotated sources](references.md)

This comparison preserves the scientific detail behind the modular proposal.
Published results below are reported results, not a new common-data reproduction.
Yang 2024 predates Figiel 2025 and compares primarily against its own ePPC-SAFT
predecessors. The papers do not establish a controlled Yang-versus-Figiel ranking.

## Complete formulations first

On the common reduced molar basis $a=A/(nRT)$, the source decompositions are

$$
a^{res}_{Yang}=a^{hc}+a^{disp}+a^{assoc}+a^{polar}
+a^{NAHS}+a^{MSA}+a^{Born},
$$

$$
a^{res}_{Figiel}=a^{hc}+a^{disp}+a^{assoc}+a^{DH}
+a^{Born,SSM+DS}.
$$

Yang uses an unrestricted-primitive MSA description, explicit multipolar
interactions and a nonadditive-hard-sphere correction, applied here to
water–methanol. Figiel uses DH and its modified Born expression. Their solvent
models and parameters also differ. Replacing only an association topology on
Figiel's unchanged parameters creates a new hybrid, not either published EOS.
[YangPair2024, Eq. (1), §2](references.md#yangpair2024),
[Figiel2025, theory](references.md#figiel2025)

| Feature | Bülow ePC-SAFT advanced | Figiel 2025 | Yang 2024 IP-Wertheim |
|---|---|---|---|
| Long-range ionic description | DH | DH | MSA |
| Solvation | Born with composition-dependent permittivity | Fitted Born diameter, shell and dielectric saturation | Born with explicit short-range ion–solvent association |
| Pairing | Part II introduces dissociation into DH | No explicit ion-pair sites | One dedicated pair site per ion |
| Main scientific advance | Dielectric/solvation treatment and extension to organic media | Solvation/transfer energies with useful finite-concentration properties | Separate hydration and pairing capacity; improved NaCl temperature behavior |
| Direct comparator | Earlier ePC-SAFT | Earlier ePC-SAFT advanced | Earlier ePPC-SAFT |

## Association: topology changes within Wertheim mathematics

For site type $A$ on species $i$, multiplicity $M_{iA}$, and unbonded
fraction $X_{iA}$,

$$
a^{assoc}=\sum_i x_i\sum_A M_{iA}
\left(\ln X_{iA}-\frac{X_{iA}}2+\frac12\right),
$$

$$
X_{iA}=\left[1+\rho_N\sum_j x_j\sum_B
M_{jB}X_{jB}\Delta_{iA,jB}\right]^{-1},
$$

$$
\Delta_{iA,jB}=D_{ij}^{3}g^{HS}_{ij}\kappa_{iA,jB}
\left[\exp\left(\frac{\epsilon_{iA,jB}}{k_BT}\right)-1\right].
$$

The contact diameter $D_{ij}$, contact function and cross-volume conventions
must follow the selected source; the generic expression does not make those
definitions interchangeable. Use compatible number-density and volume units.
The current definitions remain in [Engine equations](../../../../engine/docs/equations.md).

In Yang's predecessor, an ion's sites could bind either solvent or counterions,
allowing larger ionic clusters. Yang separates solvation sites from one
additional counterion-pairing site. Like-charge association is excluded. In the
reported topology, chloride has six solvation sites plus one pair site; sodium
has seven plus one. This preserves solvent coordination while restricting the
ion–ion bonding network to pairs, without increasing the number of adjustable
parameters relative to the preceding association model.
[YangPair2024, §2.1, Fig. 1, Table 1](references.md#yangpair2024)

For a single 1:1 salt with equal cation and anion mole fractions $x_s$, one
dedicated site per ion, no competing partner and symmetric strength, the pair
block satisfies $cX^2+X-1=0$, where $c=\rho_Nx_s\Delta_{pair}$. Hence

$$
X_{pair}=\frac{2}{1+\sqrt{1+4\rho_Nx_s\Delta_{pair}}},
\qquad
a_{pair}=2x_s\left(\ln X_{pair}-\frac{X_{pair}}2+\frac12\right).
$$

This is an exact algebraic derivation of that restricted block, not a general
mixture solution. Coupled solvation still requires its internal-state solution;
multiple salts introduce competing pair partners. If $\Delta$ itself depends
on a coupled dielectric state, that dependency must also be solved and
differentiated. The existence of the quadratic does not decouple all physics.

## Wertheim pairing, Bjerrum pairing and screening

Yang's IP-Wertheim uses fitted counterion association energies. Its IP-Bjerrum
alternative substitutes a Coulomb-based association integral. With the paper's
pair-distance convention,

$$
\ell_{B,ca}=\frac{b_{ca}|z_cz_a|e^2}
{8\pi\epsilon_0\epsilon_r k_BT},
\qquad
\Delta_{Bj}^{molar}=4\pi N_A
\int_{\sigma_{HS,ca}}^{\ell_{B,ca}}
r^2\exp(2\ell_{B,ca}/r)\,dr.
$$

Here $b_{ca}$ is fitted. This cutoff convention has a factor of two relative
to the common Bjerrum length; preserve the source definition. Yang's analytical
approximation is bounded by $2.1\le 2\ell_B/\sigma_{HS,ca}\le16$.
The integral requires a physically meaningful interval and its stated domain;
do not extrapolate it into a negative association strength. Convert the molar
strength before inserting it into a number-density kernel.
[YangPair2024, §2.1](references.md#yangpair2024)

Bülow Part II instead introduces a dissociation degree into DH, including

$$
\kappa_D^2=\frac{\rho_Ne^2}{\epsilon_0\epsilon_r k_BT}
\sum_i\alpha_i x_i z_i^2,
$$

and modifies closest approach. That changes the screening population, whereas
Yang's pairing sites alter short-range bonding. Bülow reports modest effects on
the investigated salt solubilities and neglects composition derivatives of
$\alpha$. Reusing that approximation would require explicitly declaring it;
it cannot silently satisfy the Engine's complete derivative requirement.
[BulowII2021, §2.3 and conclusions, with corrigendum](references.md#bulowii2021)

Yang likewise identifies missing feedback from pairing into electrostatics.
Adding one pair model on top of another without reconciling screening and
species definitions is not supported by either paper.

## Corrected Figiel Born model and dielectric response

The following expression is explanatory; the current canonical definition and
correction reconciliation remain in [full EOS §7](../../../../engine/docs/science/full-eos.md#7-canonical-born-source-reconciliation).
With $D_i=d_i^{Born}+\Delta d_i$,

$$
a^{Born,SSM+DS}=-\frac{e^2}{4\pi\epsilon_0 k_BT}
\sum_i x_i z_i^2\left[
\frac{1-\epsilon_{bulk}^{-1}}{D_i}
+(1-\epsilon_{ion}^{-1})\left(\frac{1}{d_i^{Born}}-\frac{1}{D_i}\right)
\right].
$$

The outer region uses bulk permittivity; the inner region uses ionic-region
permittivity. The 2026 correction restores the missing bracket and changes the
second dielectric label to $\epsilon_{ion}$. The original misprint must not
be reintroduced as a distinct selectable model.
[FigielCorrection2026](references.md#figielcorrection2026)

The shell thickness is

$$
\Delta d_i=\frac{f_{mix}-1}{|z_i|}d_i^{Born},
\qquad f_{mix}=\sum_kx_k f_k.
$$

In the investigated systems, ions have $f_i=1$, and fitted solvent values
are 1.5 for water, 1.4 for methanol and 1.6 for ethanol. The ionic-region
permittivity used in the current source-bearing model is 8. These are explicit
source parameters, not universal defaults. Distances in the SI Born prefactor
must be in metres. When the shell vanishes the second term disappears; when
the two permittivities are equal the shell partition collapses to direct Born.

Figiel's bulk relation is

$$
\epsilon_{bulk}=
\frac{\sum_s w_s^{salt-free}\epsilon_s}{1+7.01x_{ion}}.
$$

The numerator uses salt-free solvent **mass** fractions, whereas $x_{ion}$
is the total ion mole fraction. The relation was correlated using water and
methanol data and checked against ethanol; the reported study centers on
298.15 K. It does not establish arbitrary solvent or temperature transfer.
[Figiel2025, Eqs. (8)–(12)](references.md#figiel2025)

Yang's principal dielectric relation uses solvent partial molar density:

$$
\epsilon_r=1+\rho_{solv}d_{V,mix}
\left(\frac{d_{T,mix}}{T}-1\right),
$$

with $d_V$ and $d_T$ averaged using salt-free solvent mole fractions.
Density and volume units must match. Its RSP variants preserve solvent
parameters, fit effective alcohol parameters, introduce salt corrections, or
combine those changes. This explicitly density-dependent relation creates
additional pressure and composition derivative paths.
[YangPair2024, dielectric model](references.md#yangpair2024),
[YangDielectric2024, Eqs. (1a), (2a), (2b)](references.md#yangdielectric2024)

## Calibration and parameter meaning

| Parameter or choice | Figiel | Yang IP-Wertheim |
|---|---|---|
| Solvent association | Water, methanol and ethanol use 2B | Water 4C, methanol 2B, ethanol 3B |
| Ion size | Inherited segment parameters; separately fitted Born size | Separate hard-sphere, MSA and Born diameters |
| Short-range interactions | Effective dispersion and pair parameters | Pseudo-unary association energies, salt-pair and alcohol–ion corrections |
| Born calibration | Primarily infinite-dilution solvation | Inherited solvation-based procedure from preceding work |
| Concentration observations | Aqueous MIAC interaction fits; NaBr used to fit solvent shell factors | Joint MIAC, vapor-pressure and density parameterization |
| Mixed-solvent work | Predictions from the stated single-solvent/ion parameterization | Fits alcohol–ion and dielectric-related parameters to mixed-solvent properties |

Yang's association-energy mixing includes

$$
\epsilon_{w,i}^{AB}=(\epsilon_w^{AB}+\epsilon_i^{AB})/2,
$$

$$
\epsilon_{c,a}^{AB}=(\epsilon_c^{AB}+\epsilon_a^{AB})(1-w_{ca})/2,
\qquad
\epsilon_{alc,i}^{AB}=(\epsilon_{alc}^{AB}+\epsilon_i^{AB})(1-w_{alc,i})/2.
$$

A pseudo-unary ion energy constructs cross interactions; it does not enable
like-ion self-association. Fifteen parameters are fitted jointly for the four
reference chloride systems before extension. “Ion-specific” does not mean
absence of salt-pair parameters. Supporting-information Tables S2–S4 retain
the different solvent parameters and solvent-pair corrections.
[YangPair2024 and SI](references.md#yangpair2024)

Two source examples illustrate effective parameterization. Na–Br and Na–I
pairing energies are only 0.006199 and 0.03311 K, compared with Na–Cl's 949.8 K.
Substantial pairing attraction is not necessary for every fitted salt. Lithium's
hard-sphere diameter is 0.0526 Å, accompanied by a large MSA diameter ratio and
additional dispersion. Transferring that isolated fitted size as a literal ion
dimension would discard the compensating full model.
[YangPair2024, parameter tables](references.md#yangpair2024)

## Reported results and their limits

Yang covers 19 aqueous alkali-halide salts, excluding LiF, and reports overall
deviations of approximately 1.1% for MIAC, 1.7% for vapor pressure and 0.36% for
density in Table 7. These are the paper's aggregate comparisons, not certified
heldout errors from a common comparison with Figiel.

| Aqueous salt | IP-Wertheim MIAC AAD | IP-Bjerrum MIAC AAD |
|---|---:|---:|
| NaCl | 0.56% | 1.2% |
| KCl | 0.54% | 2.7% |
| RbCl | 0.33% | 2.0% |
| CsCl | 0.69% | 0.55% |

IP-Wertheim improves several chloride comparisons, while IP-Bjerrum is slightly
better for CsCl. Yang describes the overall improvement over the preceding
clustering-capable ePPC-SAFT as modest. The particularly useful result is aqueous
NaCl MIAC curvature across 273.15–383.15 K: Fig. 4 captures an observed maximum
better than the predecessor's nearly linear decline, without its previous
MIAC-focused density penalty. Contribution analysis attributes this to the
balance among association, electrostatics and the remaining terms; it does
not isolate MSA as the cause or directly validate caloric quantities.
[YangPair2024, §4, Fig. 4, Table 7](references.md#yangpair2024)

| Mixed-solvent dielectric treatment | Overall reported MIAC deviation |
|---|---:|
| RSP-1: original solvent parameters | 11% |
| RSP-2: effective alcohol parameters | 3.8% |
| RSP-3: explicit salt correction | 6.1% |
| RSP-4: both changes | 3.0% |

Despite RSP-4's smaller aggregate MIAC error, Yang recommends RSP-2 because
RSP-3/4 show questionable dielectric or vapor-pressure extrapolation. The
authors also report pairing decreasing with increasing alcohol content, likely
connected to a Wertheim strength without explicit dielectric dependence, and
identify unresolved electrostatic feedback. Those limitations matter directly
to organic-rich extraction media.
[YangPair2024, Table 8, §4.3 and conclusions](references.md#yangpair2024)

Figiel instead constrains organic-solvent ion interactions with solvation and
transfer information and then predicts mixed-solvent MIACs without additional
mixture-specific fitting. The conclusions highlight useful behavior up to
roughly 2 mol kg⁻¹ in the studied mixtures, with underprediction for concentrated
LiBr and some water–ethanol–NaBr conditions. The study centers on 298.15 K.
Its prediction exercise cannot be ranked by directly comparing its deviations
with Yang's differently fitted aggregates.
[Figiel2025, Figs. 6–9 and conclusions](references.md#figiel2025)

## Why MIAC alone cannot identify the missing physics

For component amounts $n_i$, total moles $n$, and compressibility $Z$,

$$
\ln\phi_i=\left[\frac{\partial(na^{res})}{\partial n_i}\right]_{T,V,n_{j\ne i}}
-\ln Z,
$$

$$
\ln\gamma_{\pm}^{*,x}=
\frac{\sum_{i=+,-}\nu_i[\ln\phi_i-\ln\phi_i^\infty]}{\nu_++\nu_-}.
$$

The infinite-dilution reference must use the specified solvent composition,
temperature and pressure; apply the appropriate molality conversion when
comparing practical molal activities. These expressions explain the
source comparison; the canonical [electrolyte reference definitions](../../../../engine/docs/science/full-eos.md#10-electrolyte-reference-states)
own the project's exact convention.

A large absolute solvation contribution may largely cancel in the reference
subtraction. Concentration-dependent dielectric response, association and
dispersion can compensate. A small MIAC error therefore does not uniquely
establish pair fractions or single-ion solvation. Single-ion quantities also
retain their extrathermodynamic convention; salt observables alone cannot
resolve it. Joint solvation/transfer, density, dielectric and MIAC constraints
are more informative than MIAC alone.

## Additional research branches

**Rueben 2024** supplies perturbation-based dielectric response with solvent
partial density and ionic decrement parameters. Its ionic calibration is at
298.15 K; reported transfer to low-permittivity solvents is limited. The
framework motivates a candidate for comparison, not broad electrolyte
temperature validation. The native form's existence is not public admission.
[Rueben2024, Eqs. (4)–(6), §§2–4](references.md#rueben2024)

**Biazi 2025** couples a Maribo-Mogensen dielectric description to Wertheim
association and separates solvent hydrogen bonding from ion–solvent bonds.
It fits density, MIAC and relative static permittivity, and discusses an
independent-additive approximation for the two association energies. The
NaCl/NaBr MIAC fitting data extend to 573 K and 1000 bar; dielectric data span
much narrower conditions. Density increments are underpredicted near water's
critical region and at high salt concentration. A broad MIAC data range is not
uniform qualification of every property throughout that range.
[Biazi2025, theory, Table 5 and results](references.md#biazi2025)

**Silva 2023** distinguishes ways of introducing concentration-dependent
permittivity in the derivation of DH/Born models. A formulation ID such as
“DH with variable permittivity” needs the actual derivation and reference
convention. Differentiating a post-substitution equation exactly does not
reproduce a different charging-process derivation.
[Silva2023](references.md#silva2023)

**Maghsoodloo 2025** discusses mixed-solvent electrolyte modeling, parameter
degeneracy and dissociation assumptions using eNRTL. It supports retaining the
chemical description and calibration with parameters. It is adjacent modeling
evidence, not a PC-SAFT contribution or proof of SAFT compatibility.
[Maghsoodloo2025](references.md#maghsoodloo2025)

## Consequences for this Engine

| Candidate action | Assessment |
|---|---|
| Labeled ion-solvation/pair sites | Worth a bounded later study; existing explicit site multiplicities and interaction matrix provide the mathematical starting point. |
| Second general Wertheim solver or energy owner | Unnecessary for Yang's topology alone; reuse the existing owner. |
| Yang sites on unchanged Figiel parameters | Research hybrid requiring reassessment and likely refitting, not an accepted source model. |
| DH-to-MSA replacement | Not justified by these results alone; no isolated causal comparison. |
| Polar and NAHS additions | The optional Gross–Vrabec dipole–dipole term is admitted by [ADR 0011](../../decisions/0011-optional-dipole-dipole-term.md). Its combination with electrolyte families and any NAHS addition still require separate justification. |
| Bjerrum pairing | Lower-priority candidate requiring domain, screening and derivative reconciliation; not uniformly better in Yang's comparison. |
| Temperature/density/solvent validation | Required regardless of the eventual physical selection. |

The repository's `Yang2024` citation currently denotes the **dielectric** paper
(`10.1016/j.fluid.2024.114103`), not the ion-pairing paper
(`10.1021/acs.iecr.4c01579`). This proposal distinguishes `YangDielectric2024`
and `YangPair2024` to prevent accidental source substitution.

The present literature gives a credible reason to investigate ionic association.
It does not demonstrate a missing term that should be switched on immediately
in the current Figiel model. Alkali-halide evidence does not establish
multivalent Mg/Ca competition, reactive amines or extraction complexes.

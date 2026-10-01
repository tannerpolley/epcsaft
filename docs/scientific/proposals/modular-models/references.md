# Annotated primary sources and evidence locators

[Proposal overview](README.md) · [Family map](physics-and-families.md) · [Paper comparison](electrolyte-comparison.md)

This reading catalog records the sources used in the discussion and the
proposal. It is a claim-oriented reference list, not a refreshed Zotero/BibTeX
export. No library records or attachments were changed. DOI links identify
publications; Zotero parent/attachment keys identify the locally consulted
sources without embedding private machine paths.

The [expanded family atlas guide](family-tree-guide.md#evidence-behind-the-expanded-branches)
records additional sources and their claim locators, including the 2025
Binding-MSA branch. Its PDF has a unified numbered index for all 56 records.

**Evidence classes.** “Source” means a paper's equations, assumptions or reported
results; “repository” means the inspected canonical owner or retained analysis;
“synthesis” means the proposed compatibility/architecture interpretation.
Reported fits are not relabeled as independent validation. A historical map or
abstract does not supply implementation-ready equations: the selected branch's
complete derivation, parameters, corrections and numerical references must be
checked before future implementation.

## Foundation and neutral/polar branches

### Wertheim background

Wertheim's thermodynamic perturbation theory motivates association and chain
formation. This proposal uses the concrete SAFT/PC-SAFT equations in Chapman
and Gross–Sadowski below; it does not claim a fresh audit of every paper in
Wertheim's original series. The generic free-energy and implicit-derivative
identities in the proposal are explanatory derivations, not transcriptions
of a new selectable EOS.

### Chapman1989

Chapman, Gubbins, Jackson and Radosz (1989), *SAFT: Equation-of-state solution
model for associating fluids*.
[DOI](https://doi.org/10.1016/0378-3812(89)80308-5).
Zotero parent `2P8DY5AX`. Use: the SAFT branch and association/chain background.
Locator: model development. This establishes ancestry, not cross-family
parameter compatibility.

### Chapman1990

Chapman et al. (1990), *New reference equation of state for associating liquids*.
[DOI](https://doi.org/10.1021/ie00104a021).
Repository source ID `Chapman1990`. Use: Wertheim association Helmholtz and
mass-action structure. Locator: association development, as mapped by the
[canonical EOS sources](../../../../engine/docs/science/full-eos.md#1-sources-and-evidence-classes).

### Gross2001

Gross and Sadowski (2001), *Perturbed-Chain SAFT: An Equation of State Based on
a Perturbation Theory for Chain Molecules*.
[DOI](https://doi.org/10.1021/ie0003887).
Zotero parent `KXZN4KLQ`; repository `Gross2001`. Use: hard-chain reference,
dispersion and thermodynamic derivatives. Locator: theory/reference and
dispersion sections. The reference package and its fitted parameters must
remain consistent.

### Gross2002

Gross and Sadowski (2002), *Application of the Perturbed-Chain SAFT Equation of
State to Associating Systems*.
[DOI](https://doi.org/10.1021/ie010954d).
Zotero parent `2AJYEZS6`; repository `Gross2002`. Use: PC-SAFT association and
cross-interaction treatment. Locator: association model and mixing rules.

### SimplifiedPCSAFT2003

von Solms, Michelsen and Kontogeorgis (2003), *Computational and Physical
Performance of a Modified PC-SAFT Equation of State for Highly Asymmetric and
Associating Mixtures*.
[DOI](https://doi.org/10.1021/ie020753p).
Use: historical branch with modified reference/mixing approximations, not a new
additive energy. Publication identity checked through DOI metadata; no claim of
a full implementation-ready transcription or current Engine support is made.

### GrossQuadrupole2005

Library citation key: `grossEquationofstateContributionPolar2005`.

Gross (2005), *An equation-of-state contribution for polar components:
Quadrupolar molecules*.
[DOI](https://doi.org/10.1002/aic.10502).
Zotero parent `EIPM974F`. Use: quadrupolar branch identification. Locator:
quadrupolar perturbation formulation. The proposal does not select a particular
quadrupolar implementation or transfer its parameters into current scope.

### GrossVrabec2006

Library citation key: `Gross2006`.

Gross and Vrabec (2006), *An equation-of-state contribution for polar components:
Dipolar molecules*.
[DOI](https://doi.org/10.1002/aic.10683).
Zotero parent `WGN5SY8E`, PDF `DWISZT5F`. Use: explicit dipolar perturbation,
simulation-based development and application with molecular dipole moments.
Locators: abstract, theory and pure-component/mixture application sections.
Other fitted molecular parameters belong to that polar formulation; a measured
dipole input does not make a nonpolar fit automatically transferable. Online
publication was in 2005; the journal-volume citation is 2006.

### VrabecGrossDQ

Library citation key: `vrabecVaporliquidEquilibriaSimulation2009`.

Vrabec and Gross, *Vapor-liquid equilibria simulation and an equation of state
contribution for dipole-quadrupole interactions*.
[Journal article, 2008](https://doi.org/10.1021/jp072619u).
[Author preprint](https://arxiv.org/abs/0904.4637).
Zotero parent `7H3X4JVA` records the 2009 arXiv version. Use: direct evidence
that a distinct dipole–quadrupole cross contribution exists. Locator: abstract
and interaction/formulation development. The expanded atlas verifies the 2008
journal identity; the preprint date is not the journal publication date.

### Kleiner2007

Kleiner and Sadowski (2007), *Modeling of polar systems using PCP-SAFT: An
approach to account for induced-association interactions*.
[DOI](https://doi.org/10.1021/jp072640v).
Zotero parent `28XJ6NNB`. Use: cross-association for molecules that do not
self-associate. Locator: induced-association model and mixtures with water or
alcohols. “No self-association” cannot be used as a blanket mixture-site rule.

### Esper2023

Library citation key: `esperPCPSAFTParametersPure2023`.

Esper, Bauer, Rehner and Gross (2023), *PCP-SAFT Parameters of Pure Substances
Using Large Experimental Databases*.
[DOI](https://doi.org/10.1021/acs.iecr.3c02255).
Zotero parent `ZYLLGSVM`. Use: formulation-specific parameter estimation for
nonassociating, associating and dipolar substances using vapor pressures and
densities. Locator: model choices and parameter-estimation method. A parameter
library carries its fitting formulation; the proposal does not claim this
library's coverage as current Engine availability.

### Lafitte2013

Lafitte et al. (2013), *Accurate statistical associating fluid theory for chain
molecules formed from Mie segments*.
[DOI](https://doi.org/10.1063/1.4819786).
Use: adjacent SAFT-VR Mie branch with variable repulsive/attractive exponents
and its own perturbation/reference treatment. Locator: abstract and theory.
The historical distinction is supported; a full Mie implementation is outside
this proposal's selected initial scope.

## Electrolyte development

### Cameretti2005

Cameretti, Sadowski and Mollerup (2005), *Modeling of Aqueous Electrolyte
Solutions with Perturbed-Chain Statistical Associated Fluid Theory*.
[DOI](https://doi.org/10.1021/ie0488142).
Zotero parent `CNAR9DGU`. The same-year correction is retained as parent
`VZWJHRBS`, [DOI 10.1021/ie051055i](https://doi.org/10.1021/ie051055i),
PDF `RPMVXCZ2`.
Use: the early DH-based electrolyte PC-SAFT branch. Locator: electrolyte model
development. The expanded atlas inspection confirms that p. 8944 is headed
“Additions and Corrections,” replacing Eq. (13) and providing corrected
parameters in Table 2. These records do not define separate model families.

### Held2008

Held, Cameretti and Sadowski (2008), *Modeling aqueous electrolyte solutions.
Part 1. Fully dissociated electrolytes*.
[DOI](https://doi.org/10.1016/j.fluid.2008.06.010).
Zotero parent `D6FBI6FG`; repository `Held2008`. Use: ion parameterization and
full-dissociation assumptions. Locator: model and Table 2; current repository
records retain a bounded 298.15 K, 1 bar reference domain.

### Held2014

Held et al. (2014), *ePC-SAFT revised*.
[DOI](https://doi.org/10.1016/j.cherd.2014.05.017).
Use: explicit cation–anion dispersion in a revised short-range ionic model,
particularly for concentrated solutions. Locator: abstract and revised-model
discussion. This is the basis for reassessing effective interactions before
adding explicit ionic association.

### Rozmus2013

Rozmus et al. (2013), *Modeling of strong electrolytes with ePPC-SAFT up to high
temperatures*.
[DOI](https://doi.org/10.1021/ie303527j).
Zotero parent `YMIPN9IR`, PDF `X8IHURRC`. Use: polar/electrolyte combined
formulation with MSA and Born. Locator: abstract and full EOS theory. Its
parameterization is not assumed equivalent to DH-based Figiel ePC-SAFT.

### Bulow2019

Bülow et al. (2019), *Incorporating a concentration-dependent dielectric constant
into ePC-SAFT. An application to binary mixtures containing ionic liquids*.
[DOI](https://doi.org/10.1016/j.fluid.2019.03.010).
Zotero parent `WNXMNLKL`; repository `Bulow2019`. Use: concentration-dependent
dielectric treatment in the electrolyte development. Locator: dielectric and
DH model sections. Ionic-liquid results are not automatically alkali-halide
or mixed-solvent validation.

### BulowI2021

Bülow, Ascani and Held (2021), *ePC-SAFT advanced – Part I: Physical meaning of
including a concentration-dependent dielectric constant in the Born term and
in the Debye–Hückel theory*.
[DOI](https://doi.org/10.1016/j.fluid.2021.112967),
[corrigendum](https://doi.org/10.1016/j.fluid.2021.113184).
Parent `L79J5R52`, PDF `S4WEK4G3`, Markdown `8AH6M9CN`;
correction PDF `YHZIPCDI`, separate correction parent `9XTQJ9P7`.
Use: physical interpretation of shared dielectric dependence in DH/Born.
Locator: theoretical development and corrigendum. Local metadata includes a
2020 date, while the journal DOI/volume identifies the 2021 article; the DOI
is the unambiguous publication locator.

### BulowII2021

Bülow, Ascani and Held (2021), *ePC-SAFT advanced – Part II: Application to Salt
Solubility in Ionic and Organic Solvents and the Impact of Ion Pairing*.
[DOI](https://doi.org/10.1016/j.fluid.2021.112989),
[corrigendum](https://doi.org/10.1016/j.fluid.2021.113183).
Parent `QYRA4L6M`, PDF `Y3P62JXY`, Markdown `SITUIV2V`, SI `AGLEWE3W`;
correction PDF `D4M5NHIV`, separate correction parent `5E8Y9396`.
Use: dissociation degree in screening/closest approach, limited salt-solubility
impact and neglected composition derivatives. Locators: §2.3, conclusions and
corrigendum. The derivative omission remains an explicit source approximation.

### AscaniHeld2021

Ascani and Held (2021), *Prediction of salting-out in liquid-liquid two-phase
systems with ePC-SAFT: Effect of the Born term and of a concentration-dependent
dielectric constant*.
[DOI](https://doi.org/10.1002/zaac.202100032).
Zotero parent `KLBK4XDV`; repository `AscaniHeld2021`. Use: direct Born and
combined solvent/ion dielectric rules already represented in source records.
Locator: model definition; see current EOS §§6–7 for its distinction from
Figiel's radial ionic-region permittivity.

### YangParameter2023

Library citation key: `Yang2023`.

Yang et al. (2023), *Systematic evaluation of parameterization approaches for
the ePPC-SAFT model for aqueous alkali halide solutions. II. Alkali bromides,
iodides, fluorides, and lithium halides*.
[DOI](https://doi.org/10.1016/j.fluid.2023.113853).
Zotero parent `AW8TRCCS`. Use: preceding ion-specific parameterization branch
and context for the 2024 comparison. Locator: parameterization discussion.
The roadmap claims no separate numerical reproduction of this paper.

### YangDielectric2024

Library citation key: `Yang2024`.

Yang, Kontogeorgis and de Hemptinne (2024), *Composition-dependence of relative
static permittivity in ePPC-SAFT for mixed-solvent alkali halides*.
[DOI](https://doi.org/10.1016/j.fluid.2024.114103).
Zotero parent `LQMTYLCE`; existing repository ID **`Yang2024`**. Use:
density–temperature dielectric formulation and solvent mixing. Locators:
Eqs. (1a), (2a), (2b), as mapped by current EOS §6.2. It is a different paper
from the ion-pairing source immediately below.

### YangPair2024

Library citation key: `Yang2024a`.

Yang et al. (2024), *Ion Pairing in ePPC-SAFT for Aqueous and Mixed-Solvent
Alkali Halide Solutions*.
[DOI](https://doi.org/10.1021/acs.iecr.4c01579),
[DTU supporting information](https://backend.orbit.dtu.dk/ws/portalfiles/portal/373059111/Supporting_Information.pdf).
Zotero parent `5KAIFLSL`, PDF `S33CTT65`.

| Locator | Retained use |
|---|---|
| Eq. (1), §2 | Full polar/NAHS/MSA/Born ePPC-SAFT boundary. |
| §2.1, Fig. 1, Table 1 | Separate solvation and one-per-ion pairing sites; IP-Wertheim/IP-Bjerrum definitions and domain. |
| Parameterization section and parameter tables | Fifteen-parameter reference fit; separate sizes; pseudo-unary/pair parameters; Na and Li examples. |
| Fig. 4 | NaCl MIAC temperature curvature over 273.15–383.15 K. |
| Table 7 | Aqueous MIAC, vapor-pressure, density aggregates and pairing-model comparison. |
| Table 8 | Mixed-solvent RSP-1 through RSP-4 MIAC deviations. |
| §4.3 and conclusions | Preferred RSP-2 despite aggregate RSP-4 advantage; questionable alcohol/pairing trend; missing electrostatic feedback. |
| SI Tables S2–S4 | Solvent parameters, topology and additional pair corrections. |

These reported comparisons have their own fitting definitions. They do not
establish superiority to the later Figiel paper or validate multivalent and
reactive extraction chemistry.

### Figiel2025

Paul Figiel, Gangqiang Yu and Christoph Held (2025), *Predicting Thermodynamic
Properties of Ions in Single Solvents and in Mixed Solvents Using a Modified
Born Term within the ePC-SAFT Framework*, *Industrial & Engineering Chemistry
Research* 64(18), 9406–9418.
[DOI](https://doi.org/10.1021/acs.iecr.5c00475),
[supporting information](https://doi.org/10.1021/acs.iecr.5c00475.s001).
Zotero parent `YEENX35P`, PDF `WE6I8RN3`, Markdown `2M74BKJG`;
SI PDF `3DKPDJ45`, Markdown `D5NILP7Y`.

| Locator | Retained use |
|---|---|
| Theory, Eq. (7) with 2026 correction | Shell and dielectric-saturation Born formulation. |
| Eqs. (8)–(12) and parameter tables | Shell factors, fitted Born diameters, dielectric averaging/decrement and calibration. |
| Eqs. (17), (19), (20)–(22) | Solvation, transfer and activity/reference-state definitions. |
| Figs. 6–9 and conclusions | Mixed-solvent predictions, concentrated-solution shortcomings and bounded approximately 2 mol kg⁻¹ performance. |
| SI | Direct/SSM/SSM+DS comparisons and supporting parameters/results. |

The study's 298.15 K emphasis does not establish Yang's temperature-range
claim. Retain source-specific convention and calibration differences.

### FigielCorrection2026

Figiel, Yu and Held (2026), correction to the preceding paper.
[DOI](https://doi.org/10.1021/acs.iecr.6c01415).
PDF attachment `QS7QHR5L` under parent `YEENX35P`. Locator: corrected Eq. (7).
Authority: restores the bracket and substitutes ionic-region for the second
bulk-permittivity label. The current canonical equation already adopts it.

## Dielectric coupling and calibration interpretation

### Rueben2024

Rueben et al. (2024), *Permittivity Modeling in Electrolyte PC-SAFT*.
[DOI](https://doi.org/10.1021/acs.jced.4c00347).
Zotero parent `VKRDYYKN`, PDF `R4FIF8WI`, Markdown `UUH34V28`.
Locators: Eqs. (4)–(6), §§2–4 and PT2 ion parameters. Use: partial-density
perturbation dielectric model and ionic decrement, 298.15 K ionic calibration,
and limited low-permittivity-solvent transfer. No broad temperature or
electrolyte coverage is inferred from the framework alone.

### Biazi2025

Biazi et al. (2025), *Modeling aqueous solutions with ePC-SAFT including a
relative static permittivity model connected to Wertheim's Association Theory*.
[DOI](https://doi.org/10.1016/j.molliq.2024.126726).
Zotero parent `FRL3PY4D`, PDF `5AKEMWKV`.
Locators: dielectric/association theory and its additive-independence caveat;
Table 5 fitting-data ranges; density, MIAC and permittivity results/conclusions.
Use: joint hydration/dielectric modeling and its property-specific limits.
The 2025 issue date and 2024 DOI stem refer to the same paper. Wide MIAC
conditions are not uniform dielectric or density validation.

### Silva2023

Silva, Liang and Kontogeorgis (2023), *How to account for the concentration
dependency of relative permittivity in the Debye–Hückel and Born equations*.
[DOI](https://doi.org/10.1016/j.fluid.2022.113671).
Zotero parent `6ZGX63EB`, PDF `9JABVAJE`. Locators: electrostatic/charging
derivations and conclusions. Use: model differences from the stage at which
concentration dependence is introduced. Complete differentiation is necessary
but does not equate different derivations.

### Maghsoodloo2025

Maghsoodloo et al. (2025), *Thermodynamic Modeling of Mixed Solvent Electrolyte
Solutions: Challenges and Practical Guide*.
[DOI](https://doi.org/10.1021/acs.jced.5c00412).
Zotero parent `BS2SL7QC`, PDF `MMRNBNT9`. Locators: parameter estimation,
dissociation assumptions and model-comparison discussion. Use: identifiability
and chemical-model/calibration linkage. This eNRTL work is adjacent evidence,
not a PC-SAFT formulation or a compatibility proof.

## Repository evidence and interpretation boundary

| Retained owner | Use in this proposal |
|---|---|
| [Scientific context](../../CONTEXT.md) | Shared ownership, one wheel/driver, explicit inputs, mathematical versus empirical validity. |
| [Current implementation plan](../../implementation-plan.md) | Current availability and qualification gaps; unchanged by this proposal. |
| [Equations](../../../../engine/docs/equations.md) and [full EOS](../../../../engine/docs/science/full-eos.md) | Canonical equations, units, corrections, source conventions and dielectric forms. |
| [Algorithms](../../../../engine/docs/science/algorithms.md) | Implicit states, derivative and density requirements. |
| [Packet contract](../../../../data/docs/packet-contract.md) | Reusable parameter/observation ownership, provenance and immutable records. |
| [Association study](../../../../analyses/2026-association-native-comparison/notebook.qmd) | Bounded neutral fixed-volume values/first/second-derivative evidence at snapshot `0ae7806`. |
| [ADR 0003](../../adr/0003-explicit-parameter-inputs.md), [ADR 0006](../../adr/0006-retire-residual-polar.md) | Explicit physics and current residual-polar exclusion. |

The complete-model selectable object, compatibility levels, proposed records,
implementation mapping and five-step sequence are **design synthesis** from
these sources and the user's intended workflow. They are not recommendations
quoted from one paper, accepted new equations, or a completed pre-build review.

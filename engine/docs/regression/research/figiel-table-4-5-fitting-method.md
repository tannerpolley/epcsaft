# Figiel (2025) aqueous Table 4/5 fitting method

**Question.** What fitting method for the aqueous `k_ij` values in Figiel,
Yu, and Held (2025) is actually stated by the local primary literature, and
what is the narrowest reproducible reconstruction when details are omitted?

## Bottom line

The 2025 paper states the parameter roles and order, but not a complete
regression recipe.  It says that (i) ion Born diameters were adjusted to
infinite-dilution water solvation Gibbs energies, (ii) one solvent-specific,
ion-independent `f_k` was adjusted to NaBr MIAC data, and (iii) ion--water and
ion--ion `k_ij` values were adjusted to MIAC literature data for different
salts in water; it then performed an iteration, with minimal parameter change
(main paper lines 275--283).  The paper does **not** state
an MIAC objective/residual, row weights, exact row subset/cutoff, optimizer,
starts, bounds, iteration count, or termination tolerance.

The most defensible *reconstruction* is therefore an explicit analogue of
Held et al. (2014), not a claim about the authors' hidden implementation:

```text
r_q(theta) = 1 - gamma_pm,m,calc(q; theta) / gamma_pm,m,exp(q)
min_theta sum_q r_q(theta)^2
```

This mirrors Held et al.'s published unweighted relative-squared osmotic-
coefficient objective, `Σ(1 - phi_calc/phi_exp)^2` (2014 Eq. 20), while using
the 2025 paper's molality-scale MIAC observable.  Equal row weighting is an
inference; no uncertainty weighting is published.  Do not describe this
residual, any bounds, starts, Ceres settings, or a solver choice as Figiel's
reported method.

## Verified source statements

### 2025 main paper

The retained primary-paper transcription (SHA-256
`ce80533925a91bc59d8d0d8056113c40611ca26c2edf04aced76986d50bd4bae`) states:

| Claim | Locator | Evidence status |
| --- | --- | --- |
| All MIACs are reported on the molality scale. | lines 265--271 | verified |
| Born diameters were adjusted to infinite-dilution water solvation Gibbs energies; each solvent's ion-independent `f_k` was adjusted to NaBr MIAC data. | lines 273--277 | verified |
| Ion--water and ion--ion `k_ij` values were adjusted to experimental MIAC literature data for different salts in water. | lines 279--280 | verified |
| Fit order was `d_Born`, `f_k`, then aqueous `k_ij`; an iteration followed and parameter changes were minimal. | lines 281--283 | verified |
| Table 4/5 values are valid at 298.15 K; the Table 5 footnote says ion--solvent cells were fitted and solvent--solvent cells inherited. | lines 315--341 | verified |
| Water `f_k` was fitted only to NaBr in water and then held constant for other calculations. | lines 359--361 | verified |
| For ethanol, the paper explicitly says only MIAC data up to 0.25 mol kg⁻¹ were used for `f_ethanol`; it gives no analogous water cutoff. | lines 384--387 | verified (ethanol only) |

The active 11-cell aqueous tuple used by the current Regression scope is read
directly from Table 4/5 (lines 315--341):

```text
water--Li+  -0.4    water--Na+  -0.3    water--K+  -0.1
water--Cl-  -0.3    water--Br-  -0.3
Li+--Cl-     0.8    Na+--Cl-     0.8    K+--Cl-    0
Li+--Br-     0.5    Na+--Br-    0.65    K+--Br-   -0.35
```

The printed tables contain additional H⁺, I⁻, SO₄²⁻, VO²⁺, and V³⁺ cells;
their fitting data and the treatment of blank cells are outside this 11-cell
scope.

### Supporting Information

The official local SI transcription (SHA-256
`85bd39f727158d5a9d6eea6828c1673f73850e783a655b09660cc9b66d84321a`) only
documents model screening and target tables.  SI S2 says candidate Born-term
approaches were screened against NaBr MIACs in water, retaining approaches
with AARD < 10%, then tested against NaBr MIACs in methanol; the `s-SSM+DS-f`
approach was selected (SI lines 55--80).  SI Table S4/S5
contains vanadium density/osmotic data and water-solvation literature values
(SI lines 181--260).  No objective, weights, starts,
bounds, optimizer, or iteration controls for Tables 4/5 are given.

### Predecessor fitting lineage

* **Held et al. (2014), ePC-SAFT revised** (SHA-256
  `b8b1e46bf870224de5de68b5989f9cb377d17445d87109a5462a94f1efaafbda`):
  data were generally limited to atmospheric pressure, 298.15 K, and salt
  concentrations ≤5 mol kg⁻¹.  The published objectives are a density
  difference-to-pure-water sum of squares (Eq. 19) and
  `OF_OC = Σ(1 − phi_calc/phi_exp)^2` (Eq. 20), summed over all points
  (lines 251--268).  The sequence was sensitivity analysis
  of ion dispersion energies (with `u_ion/k_B ≤ 400 K`), then adjustment of
  water--ion `k_ij`, ion--ion `k_ij`, and ion diameters, followed by successive
  adjustment of other ions.  No optimizer name, starts, or `k_ij`/diameter
  bounds are stated.
* **Held, Cameretti, and Sadowski (2008)** (SHA-256
  `dfcac60a91d6b0ec5de23dfb171c276ee75d290f331a047e7789f971b5043ff6`): ion
  diameter and dispersion parameters were fit to density data at 20--30 °C
  plus 25 °C MIAC data (lines 169--171); 14 salts were
  used for the Na⁺/Li⁺/K⁺ × F⁻/Cl⁻/Br⁻/I⁻/OH⁻ family (LiF omitted), then other
  ions were fit successively (lines 351--356).  No exact
  objective, optimizer, starts, or bounds are stated.
* **Cameretti, Sadowski, and Mollerup (2005)** (SHA-256
  `09828f160ea0089fdd1d923a18d7921e48f77bf4ffc8bcbb14fe140383a68e69`):
  describes a nonlinear least-squares fit; nine alkali-halide salt solutions
  (192 vapour-pressure and 189 density points) were fit simultaneously, with
  binary `k_ij` fixed to zero (lines 223--257).  No
  optimizer, starts, bounds, or explicit objective formula is given.
* **Bülow, Ascani, and Held (2020)** (SHA-256
  `6128216697888c5a6088e15d5dd3c798a63ce5b616bc52fc4506d53d116d2397`):
  confirms that the aqueous ion parameters and `k_ij` values were inherited
  from the 2014 work and that no ion--organic `k_ij` values were added for its
  prediction study (lines 165--180).
* **Bülow, Ascani, and Held (2021), Part II** (SHA-256
  `b8de66d4dae1b58f9e23056962f332afecc7c56c0f59ed2c9f7aef4208bedcf8`):
  states that all ion-related parameter estimation used aqueous salt data
  (lines 348--350); it does not add an objective or solver
  recipe.

## Data, fixed quantities, and reconstruction boundary

* **Aqueous MIAC rows.** The local Validation packet transcribes Hamer and Wu
  (1972), Tables 9, 10, 17, 28, 29 for LiCl, LiBr, NaCl, NaBr, KCl, and KBr:
  164 rows total (LiCl 29, LiBr 29, NaCl 29, NaBr 21, KCl 28, KBr 28),
  0.001 ≤ molality ≤ 6 mol kg⁻¹.  CSV SHA-256 is
  `2f63e13f06a5b0f4e8bca2980b6a8d9d7fb0f839153c43e3a71952daf9796595` and
  metadata YAML SHA-256 is
  `7a1880502848c4d3ab5ad18f69bb9700136517a9d44ae9d5ab57d7b521160345`.
  This is a reproducible local reconstruction of “different salts in water,”
  not a paper-verified row cutoff: the 2025 paper cites Hamer--Wu for Figure 5
  but does not enumerate rows or a concentration limit.
* **Preceding stages.** For a staged replay, hold the five active Born
  diameters (Li⁺, Na⁺, K⁺, Cl⁻, Br⁻) after fitting SI Table S5 reported
  averages, then hold the water `f_k` after fitting NaBr MIACs (the packet has
  21 rows).  The Table 4/5 stage then fits the 11 aqueous cells above against
  the 164-row packet.  These row choices and the five-ion restriction are
  Regression reconstruction decisions, not disclosed Figiel settings.
* **Fixed model state.** Use 298.15 K and 1 bar (the Table 4/5 validity
  condition).  Pure-ion parameters in Table 3 are inherited from the cited
  ePC-SAFT work; solvent--solvent `k_ij` cells are inherited.  No source gives
  uncertainty values or a weighting matrix.

## Unknowns that must remain explicit

The local primary sources do not identify: (1) whether Figiel minimized MIAC
relative, absolute, logarithmic, or another residual; (2) whether density,
osmotic, or other observables entered the Table 4/5 aqueous fit; (3) the exact
MIAC rows and any high/low-concentration exclusion; (4) row/property weights;
(5) optimizer implementation; (6) starts, bounds, scaling, random seeds;
(7) iteration count and stopping tolerance; or (8) whether all Table 4/5
cells were fit jointly or in additional sub-stages.  Any implementation must
record these as Regression choices and must not attribute them to Figiel.

## Source locators and hashes

All paths below are local, read-only primary-paper/SI transcriptions or the
hash-bound local data packet used to make the reconstruction reproducible.

| Source | SHA-256 | Locators used |
| --- | --- | --- |
| `~/Zotero/storage/2M74BKJG/Figiel, Yu, Held - 2025 - Predicting Thermodynamic Properties of Ions in Single Solvents and in Mixe.md` | `ce80533925a91bc59d8d0d8056113c40611ca26c2edf04aced76986d50bd4bae` | main paper lines 265--283, 315--341, 359--361, 384--398 |
| `~/Zotero/storage/3DKPDJ45/Supporting information Predicting thermodynamic properties of ions in single solvents and in mixed solvents using a modified Born term within the ePC-SAFT framework.mmd` | `85bd39f727158d5a9d6eea6828c1673f73850e783a655b09660cc9b66d84321a` | SI lines 55--80, 181--260 |
| `~/Zotero/storage/TZJJENHW/Held et al. - 2014 - ePC-SAFT Revised.md` | `b8b1e46bf870224de5de68b5989f9cb377d17445d87109a5462a94f1efaafbda` | lines 251--268 |
| `~/Zotero/storage/NG7QD2QZ/Held, Cameretti, Sadowski - 2008 - Modeling aqueous electrolyte solutions. Part 1. Fully dissociated.md` | `dfcac60a91d6b0ec5de23dfb171c276ee75d290f331a047e7789f971b5043ff6` | lines 169--171, 351--356 |
| `~/Zotero/storage/NUQK8KKZ/Cameretti, Sadowski, Mollerup - 2005 - Modeling of Aqueous Electrolyte Solutions with Perturbed-Chai.md` | `09828f160ea0089fdd1d923a18d7921e48f77bf4ffc8bcbb14fe140383a68e69` | lines 223--257 |
| `~/Zotero/storage/8AH6M9CN/Bülow, Ascani, Held - 2020 - ePC-SAFT advanced - Part I Physical meaning of including a concentratio.md` | `6128216697888c5a6088e15d5dd3c798a63ce5b616bc52fc4506d53d116d2397` | lines 165--180 |
| `~/Zotero/storage/SITUIV2V/Bülow, Ascani, Held - 2021 - ePC-SAFT advanced – Part II Application to Salt Solubility in Ionic and.md` | `b8de66d4dae1b58f9e23056962f332afecc7c56c0f59ed2c9f7aef4208bedcf8` | lines 348--350 |
| `~/Workspaces/Engineering/ePC-SAFT-project/ePC-SAFT-validation/data/hamer-wu-1972-aqueous-alkali-halides.csv` | `2f63e13f06a5b0f4e8bca2980b6a8d9d7fb0f839153c43e3a71952daf9796595` | metadata: six salts, 164 retained rows, 0.001--6 mol kg⁻¹ |
| `~/Workspaces/Engineering/ePC-SAFT-project/ePC-SAFT-validation/data/hamer-wu-1972-aqueous-alkali-halides.yaml` | `7a1880502848c4d3ab5ad18f69bb9700136517a9d44ae9d5ab57d7b521160345` | source DOI, table mapping, transformation record |

The official PDF byte hashes are retained for audit but are not needed for
line-level locators: main paper `a3c940895c530f72f47f22f8c0f4796b5ad37918a4edd274f85497c6dd81ad1f`,
SI `005b38ed566ec3c09b87e1ca3a9dd6eeafc9ba75e1a30b9322291d770bb93895`,
Held 2014 `dea9aa05e2ee8eb1c675873fa9d8737943312484874d844f01b45613da261acf`,
Held 2008 `20855b3de708fd529b29c7c4d8f6d0a1270802b610408606db1488b32d9c97e8`,
and Cameretti 2005 `2fc7c1ce1674c38a0bf9c9997a2d56d2ef3f319a6ec640901d82d199135897c0`.

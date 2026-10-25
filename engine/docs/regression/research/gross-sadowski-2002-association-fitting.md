# Gross--Sadowski (2002) pure-association parameter fitting

**Question.** Does Gross and Sadowski (2002) fit the two association
parameters independently, or jointly with the three ordinary PC-SAFT pure
parameters, and does the paper retain enough information to reproduce that
fit?

## Bottom line

**Verified:** the paper says that all five pure-component parameters are
adjusted **simultaneously** against vapor-pressure and liquid-density data:

```text
m, sigma, epsilon/k, epsilon^A_iB_i/k, kappa^A_iB_i
```

The five are not presented as a two-parameter association-only regression.
This is stated in the pure-component section (PDF p. 5512, continued on p.
5513; local Markdown lines 100--112) and repeated in the conclusions (PDF p.
5514; local Markdown lines 164--166).  Therefore, fitting only association
energy and association volume while freezing `m`, `sigma`, and `epsilon/k`
would be a **conditional reconstruction**, not the Gross--Sadowski procedure.

**Verified:** the paper does not publish the raw vapor-pressure or liquid-
density rows, row counts, an objective equation for the five-parameter fit,
property weights, residual scales, optimizer, bounds, starts, stopping
criteria, or a data-selection rule.  Table 1 gives only fitted values,
aggregate AAD percentages, temperature ranges, and pointers to three external
property compilations.  Figures show plotted points but are not a numerical
data table.  Consequently the paper alone cannot support a source-faithful
replay, even though it proves the joint-fit structure.

**Engineering implication (inference):** acquire the cited source rows and
document an explicit reconstruction before claiming a Gross--Sadowski fit.
The current Regression pure-2B surface is useful derivative/mechanics
evidence, but it is intentionally a one-parameter-at-a-time, one-density-row
surface and does not implement the paper's five-coordinate objective.

## Primary-source evidence

The source artifacts inspected were:

| Artifact | SHA-256 | Locators used |
|---|---|---|
| `/home/tnnrpolley21/Zotero/storage/UGT3UV2L/Gross and Sadowski - 2002 - Application of the PC-SAFT equation of state to associating systems.pdf` | `e4991a757cad356c2d0f6ec20972adb74458218b2583cd92b19b825ccf1d89bb` | PDF pp. 5510--5515; especially Table 1 on p. 5511, Eqs. 1--3 and §3 on pp. 5511--5512, §4.1 and conclusions on pp. 5512--5514 |
| `/home/tnnrpolley21/Zotero/storage/FZVTDKXU/Gross, Sadowski - 2002 - Application of the PC-SAFT equation of state to associating systems.md` | `dc4695f03a2511f0ac416bfb54923ed2b7b7a9ced8240d10b112b42ad977d732` | lines 27--55 (Table 1), 76--112 (model and pure fit), 124--166 (results/conclusions), 174--207 (references) |

The PDF page numbers above are the journal page numbers, not the local PDF
index.  The equation-heavy p. 5511 table/model page was also visually checked
against the PDF rendering; the Markdown transcription preserves the same
parameter columns and numerical values.

### What was fitted

| Claim | Evidence | Status |
|---|---|---|
| PC-SAFT is written as `Z = 1 + Z^hc + Z^disp + Z^assoc` (Eq. 1). | PDF p. 5511, Eq. 1; Markdown lines 76--86 | verified |
| Two pure association coordinates are association energy `epsilon^A_iB_i/k` and effective association volume `kappa^A_iB_i`. | PDF p. 5511, §2; Markdown lines 88--89 | verified |
| The authors assign two association sites to every associating component (2B). | PDF pp. 5511--5512, §3; Markdown lines 100--110 and Table 1 footnote line 55 | verified |
| The pure fit adjusts `sigma_i`, `m_i`, `epsilon_i/k`, association energy, and association volume—five parameters—for each component. | PDF p. 5512 (text continues at p. 5513); Markdown lines 112 and 273--278 of the current Regression design note | verified |
| The five parameters are fitted simultaneously to vapor pressure and liquid density. | PDF p. 5512; Markdown line 112; conclusion PDF p. 5514 / Markdown line 166 | verified |
| Binary cross-association uses Eqs. 2--3 (arithmetic mean energy and a geometric/diameter correction for volume); no extra cross-association correction is introduced. | PDF p. 5511, Eqs. 2--3; Markdown lines 88--98 | verified |
| In mixture calculations, the only adjusted binary parameter is dispersive `k_ij`; cross-association uses the stated combining rules. | PDF pp. 5511--5514; Markdown lines 98, 124, 156--166 | verified |

The word “simultaneously” is material: the paper does not say “fit the
association pair after fixing the three PC-SAFT parameters.”  It reports one
five-coordinate pure-component identification for each of 18 associating
substances (abstract, PDF p. 5510; conclusion, PDF p. 5514).

### Parameter identities and units

Table 1 (PDF p. 5511; Markdown lines 27--55) identifies the following units:

| Coordinate | Meaning | Unit in Table 1 |
|---|---|---|
| `m_i` | segment number | dimensionless |
| `sigma_i` | segment diameter | Angstrom |
| `epsilon_i/k` | segment energy divided by Boltzmann constant | K |
| `kappa^A_iB_i` | effective association volume | dimensionless |
| `epsilon^A_iB_i/k` | association energy divided by Boltzmann constant | K |

Molar mass `M_i` is tabulated in g/mol but is not described as a fitted
coordinate.  The table footnote fixes `N^site = 2` for all substances.

## Methanol, ethanol, and water

The three relevant Table 1 rows are reproduced below so that the provenance
and the association scheme are unambiguous.  The AAD columns are ordered as
the table heading states: vapor pressure (`P^sat`) then liquid-density/volume
property (`v`).

| Component | `m` | `sigma` (Angstrom) | `epsilon/k` (K) | `kappa^A B` | `epsilon^A B/k` (K) | AAD `P^sat` (%) | AAD `v` (%) | T range (K) | Table 1 refs |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| methanol | 1.5255 | 3.2300 | 188.90 | 0.035176 | 2899.5 | 2.36 | 2.01 | 200--512 | 1, 2 |
| ethanol | 2.3827 | 3.1771 | 198.24 | 0.032384 | 2653.4 | 0.99 | 0.79 | 230--516 | 1, 2 |
| water | 1.0656 | 3.0007 | 366.51 | 0.034868 | 2500.7 | 1.88 | 6.83 | 273--647 | 1 |

**Methanol and ethanol (verified):** both use the same two-site 2B topology,
the same five-coordinate simultaneous fit, and literature data spanning the
ranges shown above.  Their source references are VDI-Wärmeatlas (ref. 1) and
Daubert et al., *Physical and Thermodynamic Properties of Pure Chemicals* (ref.
2); see the Table 1 footnote (PDF p. 5511; Markdown lines 55 and 174--188).

**Water (verified):** water is also assigned the 2B model and is fitted with
the same five-coordinate procedure.  The authors explicitly call this a
simplification: Economou and Tsonopoulos are cited as favoring four sites,
whereas Suresh and Elliott found two sites at least as good; Gross and Sadowski
choose two sites “for simplicity” (PDF pp. 5511--5512; Markdown lines 100--110).
The paper separately warns that cyclic water topologies and dipole--dipole
interactions are omitted; its liquid-density AAD is 6.83% versus 1.88% for
vapor pressure (PDF p. 5512; Markdown lines 124--127).  Water's Table 1 data
pointer is VDI-Wärmeatlas alone (ref. 1).

This is not a separate association scheme for water in the 2002 fit.  The
three rows all use `N^site = 2`; the four-site model appears only as a cited
alternative and comparison context.

## Data actually retained by the paper

Table 1 supplies aggregate outputs, not observations:

* Each row has a temperature range and AAD percentages for vapor pressure and
  liquid density/volume.
* The cited sources are property compilations: VDI-Wärmeatlas (1994), Daubert,
  Danner, Sibul, and Stebbins (1989), and Chao et al. (1990).  The mapping is
  explicit in the Table 1 footnote (PDF p. 5511; Markdown line 55).
* Methanol and ethanol point to refs. 1 and 2; water points to ref. 1.  No row
  count, pressure value, density value, uncertainty, phase-volume convention,
  source table number, or extraction subset is printed for any of the three.
* Figure 1 plots saturated liquid and vapor densities for methanol,
  1-pentanol, and 1-nonanol, but the plotted symbols are not accompanied by a
  machine-readable table (PDF p. 5512; Markdown lines 104--105 and 124--127).

Therefore the exact raw series used by the authors is **absent** from both the
paper and the local Markdown transcription.  The external compilations may be
acquired separately, subject to access/licensing and a documented extraction,
but choosing rows from them would be a reconstruction decision unless the
authors' exact subset can be established.

## Objective, weighting, and numerical method

The following are deliberately marked `unknown`, rather than inferred from
the AAD columns:

| Item | What the 2002 source says | Status |
|---|---|---|
| Objective equation | “Simultaneously fitting vapor pressure data and liquid density data” and reporting “dimensionless absolute average deviations (AAD%)”. No minimization equation is printed in the paper. | unknown |
| Residual definition | AAD implies an aggregate absolute relative deviation may have been reported, but the paper does not define whether the fit minimized absolute, squared, logarithmic, or another residual. | unknown; the AAD wording is not an objective definition |
| Property weighting/scales | No uncertainty, row weighting, pressure/density scale, or relative-vs-absolute weighting is stated. | unknown |
| Algorithm | No optimizer or least-squares algorithm is named for the 2002 five-parameter pure fit. | unknown |
| Bounds and starts | No parameter bounds, initial values, multistart policy, stopping tolerances, or iteration controls are stated. | unknown |
| Data partitions | No training/held-out partition or row-selection rule is stated. | unknown |

The cited 2001 nonassociating PC-SAFT paper does publish a Levenberg--Marquardt
and relative-squared objective for a different task: fitting universal
dispersion-series coefficients (2001 Markdown lines 175--185, Eq. 20).  That
is **not evidence** that the 2002 associating pure-component fit used the same
objective, weights, or algorithm; the 2002 paper does not refer back to Eq. 20
for its five-coordinate regression.

## Historical implementation boundary

The implementation inspected for this research note is preserved at commit
`98a433660af7b889202b489d373c4091773429b2`. It was not a reproduction of the
2002 five-parameter identification and must not be used to claim recovery of
`epsilon^A/k` or `kappa^A B` from a single density anchor.

## Decision boundary

1. **Source-faithful route:** obtain the cited compilation rows (at minimum
   vapor pressure and liquid density over the component-specific temperature
   range), record the extraction and use basis, and establish the paper's
   objective/weights/algorithm from an additional authoritative source or mark
   them as reconstruction choices. Fit all five pure coordinates jointly.
2. **Conditional diagnostic route:** fit only the two association coordinates
   with `m`, `sigma`, and `epsilon/k` fixed to Table 1 values. This can be useful
   for sensitivity or initialization, but must be labeled conditional and not
   attributed to Gross--Sadowski.
3. **Do not** treat a one-row, one-parameter surface as either route;
   it lacks the source series and the five-coordinate objective.

The 2002 paper therefore resolves the original question—association energy and
volume were normally regressed with the other three pure parameters—but does
not, by itself, supply enough numerical detail to execute that fit faithfully.

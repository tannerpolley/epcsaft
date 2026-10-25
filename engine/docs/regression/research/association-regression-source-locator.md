# Association-parameter regression source locator

Date: 2026-07-28

Scope: read-only audit of the local Zotero library and official publisher or
library records for sources that estimate PC-SAFT/PCP-SAFT pure association
energy and association volume. This note distinguishes reproducible regression
inputs from parameter tables and plots.

## Bottom line

- **Verified:** Gross and Sadowski (2002) adjusted all five pure parameters
  (`m`, `sigma`, `epsilon/k`, `epsilon_AB/k`, and `kappa_AB`)
  simultaneously against vapor-pressure and saturated-liquid-density data for
  18 associating substances. It is not an association-only regression.
- **Verified:** The paper supplies fitted parameters, temperature ranges, and
  property AADs, but no row-level fitting data, fitting grid, objective
  equation, starts, bounds, or Supporting Information was located.
- **Verified:** The underlying Gross--Sadowski property sources are the 1994
  *VDI-Wärmeatlas* and the 1989 Daubert--Danner/DIPPR compilation. Both are
  copyrighted reference works; no open machine-readable copy of the required
  rows was located.
- **Verified:** Three later PC-SAFT papers already stored in Zotero repeat the
  same five-parameter joint pattern for alkanolamines and expose the objective
  equation and property correlations more clearly. They still omit the exact
  temperature sampling and optimization starts/bounds.
- **Inference:** Baygi and Pahlavanzadeh (2015) is the cheapest locally
  available source from which to define a new, explicitly source-bound MEA
  correlation campaign, because it prints the target correlations, units,
  range, objective, association schemes, and fitted five-parameter tuples.
  It cannot support an exact replay of the authors' optimization until the
  temperature grid (or number and placement of points), starts, and bounds are
  resolved.

## Local Zotero sources

| Source | What was regressed | Data available in the article | Objective | Starts/bounds | Classification |
|---|---|---|---|---|---|
| Gross & Sadowski (2002), DOI [`10.1021/ie010954d`](https://pubs.acs.org/doi/10.1021/ie010954d) | Five simultaneous pure parameters for 18 associating substances; 2B topology | Table 1 gives fitted parameters, temperature ranges, and separate vapor-pressure/liquid-density AADs. The raw VDI/Daubert rows are absent. No SI was located. | Not printed | Not printed | Canonical method and comparison table; **not a replay-ready dataset** |
| Nasrifar & Tafazzol (2010), DOI [`10.1021/ie901181n`](https://pubs.acs.org/doi/10.1021/ie901181n) | Five simultaneous PC-SAFT parameters for MEA (4C), DEA (4C), and MDEA (3B) | Table 2 gives tuples, ranges, and AADs; plotted/smoothed targets come from Yaws (1999), but row values and sampling grid are absent | Eq. 38: sum of absolute relative errors in saturation pressure plus saturated-liquid density | Not printed | Strong objective/formulation evidence; **not exact replay-ready** |
| Pahlavanzadeh & Fakouri Baygi (2013), DOI [`10.1016/j.jct.2012.12.021`](https://doi.org/10.1016/j.jct.2012.12.021) | Five simultaneous MDEA parameters for 4(2:2,0:0) and 6(4:2,0:0) association schemes | DIPPR-derived vapor-pressure and density correlation equations/constants, 30--170 °C range, fitted tuples and AADs; exact point grid absent | Eq. 8: sum of absolute relative saturation-pressure and density errors | Not printed | Correlation-defined targets and competing topology evidence; **grid/start/bound incomplete** |
| Baygi & Pahlavanzadeh (2015), DOI [`10.1016/j.cherd.2014.07.017`](https://www.sciencedirect.com/science/article/pii/S0263876214003360) | Five simultaneous MEA parameters for 2B, 3B, and 4C schemes | Eqs. 9--10 and Table 1 print DIPPR-derived target correlations/constants; 303.15--443.15 K range; Table 2 prints all three fitted tuples and AADs. Exact point grid absent. | Eq. 8: sum of absolute relative saturation-pressure and density errors | Not printed | **Best local bounded campaign source**, but not an exact author-run replay |
| Pakravesh & Zarei (2025), DOI [`10.1021/acs.jced.4c00390`](https://pubs.acs.org/doi/10.1021/acs.jced.4c00390) | New five-parameter fits are for P-rho-T-SAFT-HR, not PC-SAFT; PC-SAFT tuples are reused from earlier work | Article and free SI expose extensive calculated-property tables and source ranges | P-rho-T objective is printed; no new PC-SAFT objective | Detailed P-rho-T procedure is delegated to earlier work | Useful cross-check and data-source index; **not a new PC-SAFT association fit** |

### Immutable local artifact identities

| Artifact | SHA-256 |
|---|---|
| `/home/tnnrpolley21/Zotero/storage/UGT3UV2L/Gross and Sadowski - 2002 - Application of the PC-SAFT equation of state to associating systems.pdf` | `e4991a757cad356c2d0f6ec20972adb74458218b2583cd92b19b825ccf1d89bb` |
| `/home/tnnrpolley21/Zotero/storage/3G4FGGY4/Nasrifar and Tafazzol - 2010 - Vapor-liquid equilibria of acid gas-aqueous ethanolamine solutions using the PC-SAFT equation of sta.pdf` | `29165d43cf374760cc17a730d92e15ed77a9e969af1f3355a1ca39a275f7110f` |
| `/home/tnnrpolley21/Zotero/storage/YAUNIZGX/Pahlavanzadeh and Fakouri Baygi - 2013 - Modeling CO2 solubility in Aqueous Methyldiethanolamine Solutions by Perturbed Chain-SAFT Equation o.pdf` | `b737d0b83a13c5bb868d36459d1c0474eccb40e7ef700bde50fdaa77de31673a` |
| `/home/tnnrpolley21/Zotero/storage/JWH69DKG/Baygi and Pahlavanzadeh - 2015 - Application of the PC-SAFT equation of state for modeling CO2 solubility in aqueous monoethanolamine.pdf` | `7e8e77577a34bd9867489faee992dd192e8cbbc728c50a26e8264b0e09192365` |
| `/home/tnnrpolley21/Zotero/storage/TTAMHITN/Pakravesh and Zarei - 2025 - Thermodynamic Modeling of Pure, Binary, and Ternary Mixtures of Alkanolamines Using Three Versions o.pdf` | `29b1f1e115eecb11cde903aa12d1f9e6a98d2f85499879a5e05b34f429d6e0ea` |

## Gross--Sadowski property references

1. **VDI-Wärmeatlas: Berechnungsblätter für den Wärmeübergang**,
   7th expanded edition, VDI-Verlag, Düsseldorf, 1994. The
   [KIT catalog record](https://publikationen.bibliothek.kit.edu/77794)
   gives ISBN `3-18-401362-6`; `3-18-401361-8` appears in records for an
   alternate binder/format. There is no DOI for this edition. Later Springer
   editions are not the cited 1994 source.
2. **T. E. Daubert and R. P. Danner, Physical and Thermodynamic Properties
   of Pure Chemicals: Data Compilation**, Hemisphere, 1989, ISBN
   `0-89116-948-2`, OCLC `19513328`. Locate through
   [WorldCat](https://search.worldcat.org/title/physical-and-thermodynamic-properties-of-pure-chemicals-data-compilation/oclc/19513328),
   [Google Books](https://books.google.com/books/about/Physical_and_Thermodynamic_Properties_of.html?id=dt0UngEACAAJ),
   or [Open Library](https://openlibrary.org/books/OL2209327M/Physical_and_thermodynamic_properties_of_pure_chemicals).
   The [NIST catalog description](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication782-1992.pdf#page=60)
   identifies the four-volume loose-leaf compilation and supplements.

These references are the data authority cited by Gross and Sadowski, but
access to a book is not yet a source packet: the selected component rows,
edition/page/table locators, units, transformations, and license/use basis
would still need to be recorded.

## Strong download targets not currently present in Zotero

1. **Albers & Sadowski (2011), “Minimal Experimental Data Set Required for
   Estimating PCP-SAFT Parameters,”**
   DOI [`10.1021/ie2010803`](https://pubs.acs.org/doi/10.1021/ie2010803).
   This is a direct experimental-design source: it reports that five points
   can identify an associating PCP-SAFT pure-component problem and examines
   the required balance and temperature span of vapor-pressure and
   liquid-volume observations. It is method evidence, not a replacement for
   the Gross--Sadowski rows.
2. **Esper, Bauer, Rehner & Gross (2023), “PCP-SAFT Parameters of Pure
   Substances Using Large Experimental Databases,”**
   DOI [`10.1021/acs.iecr.3c02255`](https://pubs.acs.org/doi/10.1021/acs.iecr.3c02255).
   The free official SI includes initial values, bounds, diagnostics, and all
   fitted parameters as PDF, CSV, and JSON. The approximately 551,000 raw
   property rows came from public and commercial databases and are not
   redistributed in that SI. This is the strongest reproducible optimization
   contract located, but it uses PCP-SAFT rather than the exact 2002 PC-SAFT
   model.
3. **Albers, Heilig & Sadowski (2012), “Reducing the Amount of PCP-SAFT
   Fitting Parameters. 2. Associating Components,”**
   DOI [`10.1016/j.fluid.2012.04.014`](https://doi.org/10.1016/j.fluid.2012.04.014).
   It studies reduced association-parameter strategies across homologous
   families. Use it to test whether association parameters can be constrained
   externally, not as evidence that the original five-parameter PC-SAFT fit
   was separable.

## Remaining acquisition order

The Baygi--Pahlavanzadeh paper and the additional association-method papers
listed above are now present in the local Zotero Association collection. They
support the bounded reconstruction design, but they do not recover Baygi's
unpublished grid, optimizer, starts, bounds, or tolerances.

1. Obtain the exact DIPPR/Daubert or VDI records for one Gross--Sadowski
   component and freeze the row-level source packet before claiming an
   experimental-data replay.
2. Obtain Esper et al. (2023) SI for explicit starts/bounds and regression
   diagnostics, while keeping its PCP-SAFT model distinction visible.
3. Use Gross and Sadowski's fitted tuple only as the comparison target until
   the original property rows and unresolved fitting choices are recovered.

No inspected source supplied a justified association-only fit of
`epsilon_AB/k` and `kappa_AB` while holding `m`, `sigma`, and `epsilon/k`
fixed. The verified source-faithful PC-SAFT route remains a joint
five-parameter fit.

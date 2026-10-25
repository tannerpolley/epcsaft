# Observation and source availability inventory

Checked read-only in the Zotero library on 2026-09-23 using Zotero Companion (`doctor --json`, author/title searches, attachment-child inspection) and `data-search`. The articles below contain parameter fits, AADs, data-source citations, or database totals. They do **not** contain the row-level vapor-pressure and saturated-liquid-density observations needed to reproduce the pure-fluid fits. Zotero Companion returned no canonical retained acetone data package. A citation to a database or compilation is not the same as those observations being present in the library.

## Pure-fluid vapor pressure and saturated-liquid density

| Substance | Exact tables present in Zotero | What the attached tables establish | Point-level experimental/evaluated data in Zotero |
|---|---|---|---|
| Acetone | Gross & Vrabec 2006, p.1200, Table 2 (`WGN5SY8E`, PDF `DWISZT5F`); Kleiner & Sadowski 2007, p.15547, Table 3 (`28XJ6NNB`, PDF `TY4A7SB4`) | GV2006 cites data refs 1 and 2 and prints AADs, but leaves the temperature-range cell blank. Kleiner gives a 200–508 K fit range and refs 46 and 49. | No source table of the measurements. GV refs 1/2 resolve to KDB (accessed 2004) and Daubert et al. (1989); Kleiner refs 46/49 resolve to Daubert and KDB. Neither source is an attached Zotero item. |
| Butanone (2-butanone) | GV2006 p.1200 Table 2; Kleiner 2007 p.15547 Table 3 | GV2006 cites refs 1 and 2 and leaves the temperature-range cell blank. Kleiner gives 186–536 K and refs 46 and 49. | No source table of measurements. The cited KDB and Daubert sources are not attached in Zotero. |
| DMSO | GV2006 p.1200 Table 2; Kleiner 2007 p.15547 Table 4 | Both report 291–519 K. GV cites ref 1; Kleiner cites refs 46 and 49. | None in Zotero. Measured rows are now taken from the open NIST ThermoML Archive (see "Open data acquired" below). |
| Dimethyl ether | GV2006 p.1200 Table 2; Kleiner 2007 p.15546, Table 1; Gross & Sadowski 2001, p.1249, Table 2 (continued) (`KXZN4KLQ`, PDF `ZAETGZK9`) | GV and Kleiner give 200–400 K; their data refs are 2/3 and 46/47, respectively. Gross & Sadowski gives a nonpolar PC-SAFT fit over 200–400 K using refs 5 and 12 for both vapor pressure and liquid-volume data. | No source observation table is attached. The named sources are Daubert (1989), VDI-Wärmeatlas (1994), and the Gross & Sadowski reference list; the source books/atlas are not Zotero attachments. |
| Diethyl ether | GV2006 p.1200 Table 2; Kleiner 2007 p.15546, Table 1; Gross & Sadowski 2001 p.1249, Table 2 (continued) | GV gives 220–466 K; Kleiner gives 220–466 K; Gross & Sadowski gives 240–467 K for its nonpolar fit. The data refs are GV 2, Kleiner 46, and Gross & Sadowski 5/12. | No source observation table is attached. The cited Daubert and VDI sources are not Zotero attachments. |
| Acetonitrile | Kleiner 2007 p.15547, Table 4 | Gives 292–529 K and data refs 46 and 51; there is no acetonitrile row in GV2006 Table 2. | No source observation table is attached. Kleiner ref 51 is Antosik, Galka, and Malanowski (2004), *J. Chem. Eng. Data* 49, 11; Zotero search returned no item. |

In GV2006, Table 2 refs 1–4 are KDB (2004), Daubert et al. (1989), VDI-Wärmeatlas (1994), and Hsu & McKetta (1964). In Kleiner 2007, refs 46, 47, 49, and 51 are Daubert (1989), VDI-Wärmeatlas (1994), KDB (accessed 2004), and Antosik et al. (2004). In Gross & Sadowski 2001, the `nonpolar-parameters.tsv` file resolves each Table 2 row's vapor-pressure and liquid-volume reference numbers to the titles printed in that paper's footnote. Those references identify what was fitted; they do not provide the individual measured values here.

Esper et al. 2023 (`ZYLLGSVM`, PDF `AS4LVNHT`, MMD `ZP8IBLXL`) gives database-source summaries in p.15302, Table 2, and data-point counts in Table 3. It states that vapor pressure and liquid density are drawn primarily from DDB, with additional DIPPR 801 and ThermoML data, after filtering. The article says its full PCP-SAFT parameter set is in Supporting Information as PDF/CSV/JSON, but the Zotero item has only the article PDF, HTML snapshot, and MMD child; no SI or substance-level observation rows are attached. Do not treat the article's aggregate tables as pure-component data rows.

## Binary VLE papers and plotted observations

The GV2006 reference numbers and attached-article figures identify these original experimental papers. Author/title searches found none of the five original papers in Zotero.

| System and conditions | GV2006 locator and cited experiment | Original paper in Zotero? | DOI / acquisition identifier |
|---|---|---|---|
| Acetone–ethane, 25 °C (298.15 K) | p.1201, Figure 13; ref. 48: Ohgaki, Sano, and Katayama, *J. Chem. Eng. Data* 1976, 21, 55–58 | No | [10.1021/je60068a016](https://doi.org/10.1021/je60068a016) |
| Acetone–n-butane, 40 °C (313.15 K) | p.1201, Figure 14; ref. 49: Matvienko and YarymAgaev, *Russ. J. Appl. Chem.* 1994, 67, 145–147 | No | No DOI located in the DOI/Crossref search. |
| Acetone–n-pentane, 372.7, 397.7, and 422.6 K | p.1201, Figure 15; ref. 50: Campbell, Wilsak, and Thodos, *J. Chem. Eng. Data* 1986, 31, 424–430 | Yes since 2026-09-29: `R9D3UYNT`, PDF `2LB6R2IN` | [10.1021/je00046a016](https://doi.org/10.1021/je00046a016) |
| Acetone–n-decane, 40 °C and 60 °C | p.1202, Figure 16; ref. 51: Messow et al., *Z. Phys. Chem.* 1977, 258, 90–96 | No | [10.1515/zpch-1977-25815](https://doi.org/10.1515/zpch-1977-25815) |
| DMSO–toluene, 50 °C | p.1202, Figure 18; ref. 52: Nissema and Karvo, *Finn. Chem. Lett.* 1977, 8, 222–226 (author spelling as extracted from GV2006) | No | No DOI located in the DOI/Crossref search. |

Figure 15 is digitized by `scripts/digitize_fig15.py`. The owner added Campbell et al. 1986 to Zotero on 2026-09-29. Table I (p. 426) is transcribed into `inputs/campbell1986-table1.csv`, including the rows right of the azeotrope. The 2026-09-24 search below is kept as the record of open availability.

**Campbell et al. 1986 search (2026-09-24).** No open copy of the tables was found:

| Route | Result |
|---|---|
| Zotero (`zotero-search` for Campbell, Thodos, acetone, pentane) | not in the library |
| ACS PDF and article page | HTTP 403 (paywalled) |
| Unpaywall API | `is_oa: false` |
| NIST ThermoML Archive (v2020-09-30 and trc.nist.gov) | no file; the archive starts around 2003 |
| Author publication page, eng.usf.edu/~campbell | connection refused |
| GitHub code search for the DOI | no hits |
| DDBST free data pages | no acetone–pentane page (404) |

The companion papers at the same isotherms (n-pentane–methanol, n-pentane–ethanol, and the ternaries at 372.7 K) are also paywalled and are not used.

## Other polar-parameter sources searched

- Kleiner & Sadowski 2007 explicitly uses the Gross–Vrabec dipolar model (reference 23) and lists PCP-SAFT fits in Tables 1, 3, and 4. Its values, including acetonitrile, are in `polar-parameters.tsv`; its tables contain no `n_mu` field.
- de Villiers et al. 2011 (`DCI4JHHU`, PDF `C9D9NIZS`), p.176, Table 1, gives `sPC-SAFT-GV` values for acetone, 2-butanone, and diethyl ether. This is a soft-SAFT base and is marked **not transferable** to PC-SAFT. Table 1 prints `n_p`, not `n_mu`.
- Tumakaka & Sadowski 2004, Dominik et al. 2005, and Jog, Sauer, Blaesing & Chapman 2001 were searched by author/title and were not found in Zotero. Their exact DOIs and unverified availability are recorded in `polar-parameters.tsv`.
- Esper et al. 2023's SI parameter file was not attached to its Zotero item. Sauer et al. 2014 (`NSV69I6I`, PDF `I3STB7SM`) has group-contribution tables, not per-substance parameters for the requested compounds.
- Gross 2005 (`EIPM974F`) concerns quadrupoles. Vrabec & Gross 2008/2009 (`7H3X4JVA`) concerns dipole–quadrupole work, not a separate requested pure-component GV-dipole set. The `D6FBI6FG` Held/Cameretti Zotero result is electrolyte ePC-SAFT and did not supply neutral-compound GV-dipole rows.

## Open data acquired (2026-09-24)

- **Reference-EOS uncertainties.**
  - Lemmon and Span 2006: acetone subsection, open text at datapdf.com.
  - Wu, Zhou and Lemmon 2011: NIST-hosted PDF, abstract and conclusions.
  - The diethyl-ether EOS used by CoolProp 8.0.0 (`Thol-IJT-2014`) is Thol, Piazza and Span, *Int. J. Thermophys.* 35 (2014) 783, valid 270–500 K according to CoolProp and the EES documentation. It is paywalled (Springer, Unpaywall `is_oa: false`) and absent from Thol's open dissertation (d-nb.info/1109051441). The owner added it to Zotero on 2026-09-29 (`XRMZRYJI`, PDF `LSLF8W4R`). Its Sect. 3.1 deviation envelopes (pp. 799–801) are the diethyl-ether reference bounds (see README).
- **DMSO saturation rows** come from the NIST ThermoML Archive, `ThermoML.v2020-09-30.tgz` (doi:10.18434/mds2-2422, SHA-256 `231161b5e443dc1ae0e5da8429d86a88474cb722016e5b790817bb31c58d7ec2`). Of the 253 archive files that name DMSO, four supply the rows, extracted by `scripts/dmso_observations.py`:
  - vapor pressure: Fulem, Růžička and Růžička, *Fluid Phase Equilib.* 303 (2011) 205 (doi:10.1016/j.fluid.2011.02.002);
  - vapor pressure: de Pablo et al., *J. Chem. Thermodyn.* 2018 (doi:10.1016/j.jct.2018.04.014);
  - vapor pressure: Zhang et al., *J. Chem. Eng. Data* 2013 (doi:10.1021/je400531a);
  - liquid density: Paknejad, Mohammadkhani and Zarei, *J. Chem. Eng. Data* 64 (2019) 5174 (doi:10.1021/acs.jced.9b00396).
- **What remains unavailable.** Fulem et al.'s recommended vapor-pressure equation (paywalled) and the Jakli–Van Hook 1972 point data would add evidence but are not needed for these rows.

Searches also covered acetone/DMSO parameter terms and named pure-property sources (Daubert, KDB, Antosik et al.); no separate observation package or source table was located. The binary papers above are the only missing original data sources explicitly needed for the listed VLE systems.

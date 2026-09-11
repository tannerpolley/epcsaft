# Reading and reproducing the PC-SAFT family atlas

[Open the five-page vector PDF](../../latex/figures/generated/pcsaft-family-tree.pdf) ·
[Editable TikZ source](../../latex/figures/tikz/pcsaft-family-tree.tex) ·
[Model proposal](README.md)

The atlas expands the original branching diagram into **56 primary publication
records, including corrections**, with clickable paper numbers and DOI links.
It is an A3 landscape document intended for zooming or large-format printing.
Its editable source uses TikZ; the PDF contains vector paths and selectable text.

## Page map

| Page | Contents |
|---|---|
| 1 | Wertheim/Chapman foundations, the PC-SAFT trunk and nearby Huang–Radosz, soft-SAFT, SAFT-VR, Mie, SAFT-γ Mie and eSAFT-VR Mie branches. |
| 2 | PC-SAFT association/reference choices, dipolar/quadrupolar/cross terms, group construction, parameter libraries and critical-region corrections. |
| 3 | DH-based ePC-SAFT, Born/solvent-transfer developments, polar/MSA ePPC-SAFT, and shared dielectric/pairing theory. Includes the later 2025 Binding-MSA work. |
| 4–5 | Numbered author/year/short-title index with exact DOI links. Numbers in diagram nodes jump to this index. |

## What the lines mean

**Solid arrows** show selected model development or explicit reformulation.
The page-2 fan connects PC-SAFT to modeling responsibilities; each column then
contains important publications in that branch. It does not claim that every
box descends from the box immediately above it.

**Dashed arrows** connect shared theory, parameter inputs or related formulation
work. Examples include the two polar contributions informing dipole–quadrupole
work, the pure-component library feeding mixture parameterization, and the
association-based dielectric theory used by Biazi.

**Double-line arrows** denote a published correction, illustrated by Figiel's
2026 correction. The other corrections have explicit paper references in the
relevant node. Corrections are not independent physics families.

The figure is a source-backed **schematic**, not a computed citation network,
an exhaustive catalog of every SAFT paper, a priority ranking, or a compatibility
proof. Dates follow journal publication years where online and volume dates
differ. Vertical distances are not time intervals. Branches can have multiple
antecedents; the diagram does not assert exclusive parentage or unchanged
parameters along a line.

## Scientific reading notes

- Original PC-SAFT changes the reference/dispersion description within the SAFT
  lineage. The nearby SAFT-VR/Mie and soft-SAFT boxes therefore branch from the
  common SAFT background, not from PC-SAFT.
- Simplified PC-SAFT changes reference/mixing approximations. It is not an
  additional molecular attraction. The 2011 dipolar comparison uses alternative
  Jog–Chapman and Gross–Vrabec descriptions with sPC-SAFT.
- DD, QQ and DQ mean dipole–dipole, quadrupole–quadrupole and dipole–quadrupole.
  Distinct cross interactions must be selected explicitly. A residual polar
  contribution is different from a dielectric constitutive model.
- Group-contribution methods and parameter libraries change how parameters are
  obtained and, for heterosegmented descriptions, how molecules are represented.
  They are not another intermolecular force.
- Critical-region corrections address long-wavelength density fluctuations.
  The dashed connection between the two RG examples identifies this shared
  approach without claiming an exclusive paper-to-paper inheritance.
- ePC-SAFT and ePPC-SAFT denote different complete electrolyte formulations in
  this comparison. The latter includes a polar background and MSA/Born choices.
  The original discussion's Figiel and Yang parameterizations cannot be
  interchanged merely because both are descended from PC-SAFT.
- Dielectric papers in the last column are not a sequential chain of new EOSs.
  Zuber is an empirical dielectric correlation; Silva's papers clarify
  derivation/association choices; Rueben and Biazi develop different dielectric
  descriptions within electrolyte modeling.
- The Yang dielectric and ion-pairing papers are separate 2024 publications.
  The original project's `Yang2024` source identifier denotes the dielectric
  paper. The index uses distinct numbers [44] and [45].
- The 2025 BiMSA papers [55–56] are a later coupled electrostatic/pairing branch.
  Their inclusion broadens the source map beyond the earlier Yang comparison.
  It does not qualify their implementation or claim that they resolve every
  mixed-solvent, parameter-identifiability or derivative issue.

## Evidence behind the expanded branches

The [existing annotated references](references.md) retain the detailed source
locators for the original comparison. The table below covers the added or
more precisely identified branches. Bracketed numbers match the PDF index;
the DOI links identify the exact primary records. Publisher abstracts support
the high-level model classification where a complete paper was not inspected.
Such classifications do not constitute implementation-ready mathematical review.

| Papers | Primary evidence and locator | Supported use |
|---|---|---|
| [01–04] Wertheim, 1984/1986 | [Part I](https://doi.org/10.1007/BF01017362), [II](https://doi.org/10.1007/BF01017363), [III](https://doi.org/10.1007/BF01127721), [IV](https://doi.org/10.1007/BF01127722); publisher titles/abstracts and Chapman theoretical background. | Directional attraction, multiple sites and polymerization as the association foundation. |
| [07] Huang–Radosz, 1990 | [Primary paper](https://doi.org/10.1021/ie00107a014); model identity and original SAFT formulation. | A nearby SAFT branch, not a PC-SAFT descendant. |
| [09–10] SAFT-VR and soft-SAFT | [Gil-Villegas 1997](https://doi.org/10.1063/1.473101), [Blas–Vega 1997](https://doi.org/10.1080/002689797170707); original abstracts and theoretical model descriptions. | Different reference/potential choices under the shared SAFT foundation. |
| [12–13] Mie extensions | [Papaioannou 2014](https://doi.org/10.1063/1.4851455), [Selam 2018](https://doi.org/10.1016/j.fluid.2018.02.018); abstracts and Selam's non-electrolyte-model section. | Heteronuclear group construction and a DH/Born electrolyte extension of the Lafitte Mie model. |
| [19] DQ journal identity | [Vrabec–Gross 2008](https://doi.org/10.1021/jp072619u); abstract and journal issue metadata. | The DQ paper is a 2008 journal article; the previously cited arXiv copy was posted in 2009. |
| [20] Dipolar alternatives | [de Villiers 2011](https://doi.org/10.1016/j.fluid.2011.03.025); abstract/introduction comparing JC and GV formulations. Zotero parent `DCI4JHHU`. | Alternative dipolar descriptions used with simplified PC-SAFT. |
| [21] GC-SAFT | [Tamouza 2004](https://doi.org/10.1016/j.fluid.2004.06.038); abstract and parameter-construction description. | Group construction of SAFT parameters across hydrocarbon series. |
| [22] Homo/heterosegmented GC | [Sauer 2014](https://doi.org/10.1021/ie502203w); abstract and the two molecular-model descriptions. Zotero parent `NSV69I6I`. | Distinguishes segment construction; includes a heterosegmented dipolar treatment. |
| [24] Mixture parameterization | [Rehner 2023](https://doi.org/10.1007/s10765-023-03290-3); parameterization methods and conclusions. | Explicitly consumes Esper's library and extends Sauer's GC cross-association and group interactions. |
| [25] Vector GC | [Hemprich 2024](https://doi.org/10.1021/acsomega.4c04867); model development and abstract. | Group/geometry-based dipole construction for PCP-SAFT. |
| [26–27] RG corrections | [Bymaster 2008](https://doi.org/10.1021/ie8001167), [Tang 2010](https://doi.org/10.1021/ie100890d); abstracts and stated crossover/mixture approximations. | Pure-fluid and mixture density-fluctuation corrections, not extra bond energies. |
| [30] Weak electrolytes | [Held–Sadowski 2009](https://doi.org/10.1016/j.fluid.2009.02.015); source model and title. | Chemical dissociation supplements the physical EOS. |
| [40–41] ePPC development | [Ahmed 2018](https://doi.org/10.1016/j.fluid.2017.12.002), [Roa Pinto 2022](https://doi.org/10.1016/j.fluid.2022.113486); abstracts/theoretical introductions; Yang's predecessor discussion. | Mixed solvents and temperature/short-range interaction developments in the Rozmus/ePPC research line. |
| [42] First Yang parameter study | [Yang et al. 2023](https://doi.org/10.1016/j.fluid.2023.113778); exact publication identity and the 2024 paper's predecessor discussion. | Separates the reference-system study from the later Part II salt extension. |
| [46] Association/permittivity | [Maribo-Mogensen 2013](https://doi.org/10.1021/jp403375t); abstract and Biazi's theory section, where this source is reference [1]. | Shared association-linked dielectric theory informing Biazi's later formulation. |
| [47] Ion-specific decrement | [Zuber 2014](https://doi.org/10.1016/j.fluid.2014.05.037); current EOS §6.3 and source Eqs. (8)–(10). | A constitutive correlation rather than a complete selectable EOS on its own. |
| [51] Pairing representations | [Silva 2026](https://doi.org/10.1021/acs.jpcb.5c04209); abstract and SI overview. Zotero parent `5C98GI58`. | Physical/chemical association, explicit pair species and free-energy/activity relationships. No new numerical results are imported into the project. |
| [54] Early ePC correction | [Cameretti 2005 correction](https://doi.org/10.1021/ie051055i); p. 8944, corrected Eq. (13) and Table 2, inspected in Zotero PDF `RPMVXCZ2`, parent `VZWJHRBS`. | Resolves the earlier catalog's uncertainty: this record is an addition/correction, not a second model family. |
| [55–56] Binding-MSA branch | [Raeispour Shirazi et al., first paper](https://doi.org/10.1016/j.fluid.2024.114231), [improved formulation](https://doi.org/10.1016/j.fluid.2025.114396); primary publisher abstracts and model descriptions. | BiMSA within ePPC-SAFT for aqueous/mixed-solvent salts, followed by a more flexible formulation. These papers were not found by the bounded local ion-pairing search; no PDFs were acquired or Zotero records added. |

Several publications have online dates preceding their issue year. The atlas
uses 2006 for Gross–Vrabec, 2008 for the DQ paper, and 2025 for the first BiMSA
article in volume 588, even where a DOI stem or online record contains the
preceding year. Short titles in the index are reading labels; DOI identity is
authoritative. Blas–Vega uses the publisher's `10.1080/002689797170707` record;
some later reference lists give a different digitized identifier for that title.

## Reproduction and artifact ownership

The PDF is a read-only publication artifact under
`docs/scientific/latex/figures/generated/`. The single editable diagram/index
source is under `figures/tikz/`; common visual styles are in `figures/styles.tex`.
No numerical simulation, parameter fitting, library mutation or Engine build
is performed. The current roadmap and accepted physics remain unchanged.

The atlas uses an explicit A3 landscape page with one margin owner (`geometry`),
LuaLaTeX, TeX Gyre Heros, Latin Modern Math and TikZ. This isolated diagram
profile uses the standard `article` class; no manuscript template, Quarto
execution or second scientific environment is introduced. Only exercised TikZ
libraries are loaded; shell escape and externalization are unnecessary.

From the repository root:

```bash
atlas_root="$PWD"
latexmk -cd -lualatex -interaction=nonstopmode -halt-on-error \
  -outdir="$atlas_root/build/pcsaft-family-tree" \
  docs/scientific/latex/figures/tikz/pcsaft-family-tree.tex
cp build/pcsaft-family-tree/pcsaft-family-tree.pdf \
  docs/scientific/latex/figures/generated/pcsaft-family-tree.pdf
```

Compilation products and raster QA pages are disposable under
`build/pcsaft-family-tree/`; the PDF, source and this guide are retained.
The installed CSE Visualize `scripts/latex_visual_qa.py` can be run with
`--profile standalone --dpi 120 --contact-sheet` against the PDF, followed by
page-by-page color and grayscale inspection. The tool supplements visual review;
it does not establish scientific ancestry or mathematical compatibility.

## Checked disposition

Scientific logic: **schematic**. Layout: **keep**, with the expanded content
split into three diagram pages and two reference pages. The selected claims
and arrow meanings were checked against the primary-source classifications
above. This is source/editorial assessment, not independent pre-build review
of a new scientific capability.

LuaHBTeX 1.24.0 / TeX Live 2026 and latexmk 4.88 produced five A3 pages without
final-pass TeX warnings, missing references or overfull boxes. Automated text,
vector-clearance and raster QA reported **zero findings**. All five pages were
visually inspected in both color and grayscale. All 56 publication IDs are
unique, used by the diagram and resolved in its source index. The compiled
PDF's internal destinations and external DOI annotations are also checked.

The inspected source/output fingerprints are:

| Artifact | SHA-256 |
|---|---|
| `pcsaft-family-tree.tex` | `5492fb4284d05608cb2a7d6e40de5cc04364a60fb43ba2e8987753a5abbbaba4` |
| `styles.tex` | `7a9ec0fc8e9b7db0914db88c284aeed13bca0a6474417f0772f5043e00853d73` |
| `pcsaft-family-tree.pdf` | `a779fad88408c97dc182e0628675ecb9fa173354005c1c7c0a4b36ddca4f12b4` |

The retained PDF is byte-identical to the inspected build. Rebuilding can change
PDF timestamps and therefore the PDF hash without changing the drawing.

# Rezaee 2025/2026 Source Packet

This folder contains source-backed Rezaee data and a compact provenance
manifest. The paper PDFs, Markdown conversions, generated model results, and
legacy executables are not retained here; the manifest identifies the former
source packet only. Reproduction must use the current installed artifact and
the canonical source rows in this campaign and in `2025-rezaee`.

- [`combined-source-manifest.csv`](./combined-source-manifest.csv) retains compact provenance for the source packet that was split by paper year.
- [`rezaee-2026-reaction-constants.csv`](./rezaee-2026-reaction-constants.csv): reaction A/B standard Gibbs-energy changes and equilibrium constants from the 2026 supporting information.
- [`rezaee-2026-organic-pcsaft-parameters.csv`](./rezaee-2026-organic-pcsaft-parameters.csv): paper-reported DES, TOPO, RLi, and RNa PC-SAFT parameters.
- [`rezaee-2026-organic-binary-interactions.csv`](./rezaee-2026-organic-binary-interactions.csv): paper-reported organic-phase binary interaction parameters.
- [`rezaee-2026-si-density-tables.csv`](./rezaee-2026-si-density-tables.csv): density tables transcribed from the 2026 supporting information. Table S5 is DES; Table S6 is 10 wt% TOPO in DES.
- [`../../../2025-rezaee/inputs/source/rezaee-2025-headline-extraction-points.csv`](../../../2025-rezaee/inputs/source/rezaee-2025-headline-extraction-points.csv): source-backed headline extraction percentages from the ACS abstract for the 2025 experimental paper.
- [`../../../2025-rezaee/inputs/source/rezaee-2025-doe-extraction-responses.csv`](../../../2025-rezaee/inputs/source/rezaee-2025-doe-extraction-responses.csv): Table 5 designed-experiment Li/Na extraction response rows from the 2025 experimental paper.
- [`../../../2025-rezaee/inputs/source/rezaee-2025-optimum-neighborhood.csv`](../../../2025-rezaee/inputs/source/rezaee-2025-optimum-neighborhood.csv): Table 7 optimum-neighborhood validation rows.
- [`../../../2025-rezaee/inputs/source/rezaee-2025-screening-extraction.csv`](../../../2025-rezaee/inputs/source/rezaee-2025-screening-extraction.csv): Table 3 DES-ratio screening and Table 4 coextractant screening rows.
- [`../../../2025-rezaee/inputs/source/rezaee-2025-real-brine-extraction.csv`](../../../2025-rezaee/inputs/source/rezaee-2025-real-brine-extraction.csv): Table 8 real/synthetic brine composition and cation extraction efficiency rows.
- [`../../../2025-rezaee/inputs/source/rezaee-2025-extraction-equilibrium-mole-fractions.csv`](../../../2025-rezaee/inputs/source/rezaee-2025-extraction-equilibrium-mole-fractions.csv): 26 designed-experiment aqueous and organic phase mole-fraction rows transcribed from the 2025 ACS supporting information Tables S1 and S2.
- [`rezaee-2026-epcsaft-species-dataset.json`](./rezaee-2026-epcsaft-species-dataset.json): analysis-local manifest of the species, source parameter tables, reaction constants, organic binary interactions, and ePC-SAFT option context used by the Rezaee validation scripts.

The 2025 experimental paper is `10.1021/acs.iecr.4c03496`; its supporting-information file is listed by ACS as `ie4c03496_si_001.pdf`. The 2026 PC-SAFT/ePC-SAFT paper is `10.1016/j.fluid.2026.114737`; its supporting-information file is `1-s2.0-S0378381226000786-mmc1`. Rezaee 2026 states that the phase-equilibrium rows come from the 2025 work, so both source years are required for reproduction. The current benchmark set is the 26 source-backed SI equilibrium rows; treat the 2026 main-text statement of 36 equilibrium points as a clerical mismatch unless a source-backed workbook is later supplied. Do not fill source input rows from model output or inferred values. Add only source-backed rows from the 2025/2026 papers, their supporting information, or from a verified digitization workflow.

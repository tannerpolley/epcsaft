# Preparation conversion arithmetic review — 2026-09-08

Status: preparation evidence only. This record performs source-backed
arithmetic checks for the retained Khudaida parameter packets and the retained
absorber R1–R3 reaction correlations. It changes no packet, application,
Engine source, runtime behavior, or numerical tolerance.

## Scope and evidence boundary

The numerical comparisons are retained once in
[`preparation-conversions-2026-09-08.csv`](preparation-conversions-2026-09-08.csv).
The CSV is the canonical table for calculated values; this narrative refers to
row IDs so values are not duplicated in prose. The calculation uses only
already-retained files and source records. No source acquisition or wheel/build
operation was performed for this record.

## Khudaida cross-association arithmetic

Revision 4 materializes cross-association values in
`parameters/association.csv`. Revision 6 carries the same self-association
values and size correlations, while its `topology.edges` entries mark the six
cross edges as `source.kind = combining-rule` with input identities and no
materialized cross values. The three unordered pair calculations are retained
as `CA-*` rows. Each row represents both reciprocal site orientations present
in revision 4.

The retained Gross–Sadowski source evidence gives the rules

\[
\epsilon^{A_iB_j}=\frac{\epsilon^{A_iB_i}+\epsilon^{A_jB_j}}{2},
\qquad
\kappa^{A_iB_j}
=\sqrt{\kappa^{A_iB_i}\kappa^{A_jB_j}}
\left(
\frac{\sqrt{\sigma_i\sigma_j}}
{(\sigma_i+\sigma_j)/2}
\right)^3.
\]

The retained rule also states that no binary correction is added to these
cross-association quantities. The water size is evaluated from its retained
temperature correlation at the packet condition; ethanol and isobutanol use
their retained sizes. The recomputation reproduces the revision-4 decimal
records at their persisted precision. The tiny residuals in the CSV are decimal
rounding diagnostics from retaining more calculation digits, not an adopted
acceptance tolerance. These are conditional combining-rule calculations and do
not recover a source-transcribed cross fit or qualify the model.

Formula authority and provenance are retained in the project review
[`gross-sadowski-2002-association-fitting.md`](../../../engine/docs/regression/research/gross-sadowski-2002-association-fitting.md),
the local source transcription, and the retained source PDF; the project source note links the retained source locators.

## R1–R3 source-to-native representation

The absorber record stores R1–R3 as

\[
L_s(T)=a+b/T+c\ln(T/(1\,\mathrm K))+dT+o,
\]

where the JSON stores a rounded `o`. The retained doctrine defines the native
representation

\[
L_n(T)=a_n+b_n/T+c_n\ln(T/T_{\rm ref})+d_nT,
\]

with

\[
a_n=a+c\ln(T_{\rm ref}/1\,\mathrm K)+o,\qquad
b_n=b,\quad c_n=c,\quad d_n=d.
\]

The `RXN-R1` through `RXN-R3` rows use the source-to-common offsets from
the retained reaction source contract to construct `a_n`. They compare the
exact-offset source expression and native expression over the retained source
interval endpoints and transform temperature. They also compare the JSON's
rounded-offset expression against that same native expression.

The exact-offset comparison is at binary64 roundoff. The rounded-offset
comparison differs by the same temperature-independent offset rounding at all
checked temperatures. This separates representation arithmetic from physical
qualification: the identity check does not establish the source activity
convention, source applicability, coefficient selection, standard-state
validity, or provider neutral-reference behavior. The source contract and
absorber JSON remain exploratory inputs under their existing ownership.

The retained inputs are
[`reaction-system.json`](../../../downstream/MEA-Absorption-Column/src/mea_absorption_column/data/epcsaft_datasets/MEA_reactive_epcsaft_bundle/reaction-system.json),
[`chemical_reaction_source_contract.json`](../../../downstream/MEA-Thermodynamics/data/reference/MEA/manifests/chemical_reaction_source_contract.json),
and the doctrine record
[`preparation-doctrine-2026-09-08.md`](preparation-doctrine-2026-09-08.md).
Their inspected identities are recorded below.

## Reproduction and limits

Recompute the displayed rules from the linked packet self-association values,
size correlation and reaction coefficients. Cross values used high-precision
Decimal arithmetic; reaction values used binary64 math.log. The CSV reports
arithmetic discrepancies and offset rounding separately, without adopting a
numerical acceptance tolerance. It does not test the native provider.

The lead independently recomputed all six cross-energy/volume comparisons from
the packet inputs to the stored precision and checked the reaction-basis
identities separately.

Revision 6 already declares its combining-rule inputs; resolving those rules for
a consumer is implementation work under the existing input policy, not a new
owner decision. Physical convention and provider agreement remain separate
qualification work. No additional owner approval is inferred here.

## Input identities checked by the lead

| Input | SHA-256 |
|---|---|
| `data/packets/khudaida-2026-figure-2-electrolyte-lle/4/parameters/association.csv` | `ea8e7d42383d62ac7f3fa73392112c7d1e3d8fd008e2895a6a78ae928a341161` |
| `data/packets/khudaida-2026-figure-2-electrolyte-lle/4/parameters/single.csv` | `c66a460176920dbab620d69b9c442fbceb5193247e491079759b2ef43414ad9d` |
| `data/packets/khudaida-2026-figure-2-electrolyte-lle/4/parameters/correlations.toml` | `f284b10a36200624fea7f8c36ee0b6db5df891bbbbe0e4eea2ff97928b4cfbf7` |
| `data/packets/khudaida-2026-figure-2-electrolyte-lle/4/parameters/bundle.toml` | `f5d0fb91d5d6fe4837411a0d27d1e53000685b4062ecd7d6398db6e2377b1fa2` |
| `data/packets/khudaida-2026-figure-2-electrolyte-lle/4/packet.toml` | `7bfe5fb47418da30084efe1cde3c4913edac9b42c2a2069ea837fd33d10eb81c` |
| `data/packets/khudaida-2026-figure-2-electrolyte-lle/4/provenance/parameter-origin.toml` | `6f4564227912f49113b9c3908c4343420ddc5d188024eff8058b1e24794aa6f2` |
| `data/packets/khudaida-2026-figure-2-electrolyte-lle/6/parameters/parameters.json` | `05d684f3ec8a58c51a3d86bcc35b808efc9ee0a858086f705cb1e4ea011ab2e7` |
| `data/packets/khudaida-2026-figure-2-electrolyte-lle/6/packet.toml` | `5d1238112b6de4119e3347b4ff9b628a5adcf11cd2eaa37665a4fba3e25bb921` |
| `data/packets/khudaida-2026-figure-2-electrolyte-lle/6/provenance/source-origin.toml` | `aaaab21624b5df8661fdb20009539c75c07b6396b3bf969ce5c8344aecfc33b1` |
| `engine/docs/regression/research/gross-sadowski-2002-association-fitting.md` | `4d2429b1641180d951b2b792fd758cd677301e6d353c15eb3140d9f8836ca0ea` |
| `downstream/MEA-Absorption-Column/src/mea_absorption_column/data/epcsaft_datasets/MEA_reactive_epcsaft_bundle/reaction-system.json` | `810dfec15760cf74451df91743d6e63684cee93ddaf3e1ff4e42bf4a686afe29` |
| `downstream/MEA-Thermodynamics/data/reference/MEA/manifests/chemical_reaction_source_contract.json` | `39db0d7ef972dc7eb41328bdf2ec3f67f62c33fc2bf0fdc7bab471ade9aefb55` |
| `docs/scientific/reviews/preparation-doctrine-2026-09-08.md` | `89730bb8e20aa2c3a7e98de1482ed48acfb81c7edf0e2dc94a26fc3439016f53` |

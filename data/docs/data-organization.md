# Data organization

## Ownership boundary

Data is the canonical store for reusable, source-traceable scientific inputs.
That includes two kinds of packets:

1. **Parameter bundles**: self-contained ePC-SAFT component, site, association,
   binary-interaction, model, and correlation values with their source and
   domain metadata.
2. **Reusable property observations**: source-backed records such as pure-fluid
   saturation pressure and density, with stable component identity, explicit
   units, independent conditions, source locators, and uncertainty fields when
   the source provides them.

Analysis sets own question-specific evidence. Paper- or figure-specific tables, digitized
curves, direct experimental records retained only to reproduce a named paper,
campaign inputs, fitted alternatives, generated results, plots, and receipts
remain under `analyses/`. A numeric table does not become a Data packet merely
because it has been normalized. It must be reusable independently of one
campaign and have a defensible source and use basis.

Engine owns equations, algorithms, and any separately admitted self-contained
runtime copy; its Regression module owns fitting workflows and fitting results.
The installed package may not import from or fetch Data at runtime. Consumers materialize a
packet from an immutable Data commit and record its packet fingerprint.

## Packet layout

Every admitted asset uses the existing packet contract:

```text
packets/<packet-id>/<version>/
  packet.toml
  parameters/       # reusable EOS parameter bundles, when present
  observations/     # reusable source-backed property records, when present
  sources/          # distributable source material, when present
  provenance/       # origin, transformation, and use-basis records
```

Do not create a parallel top-level `species/`, `properties/`, or mutable
catalog directory. Packets are the unit of review, hashing, versioning, and
consumer locking. A packet may contain several species when they share a
source product and license boundary; split packets when source identity,
versioning, or use basis differs.

Property observation files should preserve source values and decimal strings,
use stable component identifiers where the source supports them, keep
condition units explicit in column names or unit fields, retain source row
roles, and carry source locators or source IDs. Missing source uncertainty is
represented as missing; it must not be replaced with zero or converted into a
model-acceptance threshold. Model predictions and fitted values are not
observations.

## Current parameter-bundle inventory

The repository currently contains the project’s reviewed parameter families:

| Packet family | Primary use |
|---|---|
| `ascani-2022-case-study-2-electrolyte-lle` | Ascani electrolyte-LLE source bundle |
| `crespo-2021-deea` | Crespo DEEA source bundle |
| `figiel-2025-aqueous-iodides` | Figiel aqueous iodide parameters |
| `figiel-2025-reference-electrolytes` | Figiel reference electrolyte parameters |
| `figiel-2025-validation-parameters` | Figiel validation-only parameter packet |
| `figiel-yu-held-2025` | Figiel/Yu/Held solvent and ion parameters |
| `gross-2001-methane-ethane` | Gross/Sadowski methane and ethane bundle |
| `gross-2001-propane` | Gross/Sadowski propane bundle |
| `gross-sadowski-2001` | Broad Gross/Sadowski source-complete bundle |
| `held-cameretti-sadowski-2008` | Held/Cameretti/Sadowski electrolyte bundle |
| `khudaida-2026-figure-2-electrolyte-lle` | Khudaida Figure 2 parameter bundle |

Engine runtime materializations are not additional Data ownership locations.

## Current reusable property packet

`nist-srd69-pure-saturation-properties/1` is a reviewed reusable packet with
three pure-component property tables. It currently has no live consumer or
packet lock:

- methane, NIST SRD 69, 100–180 K;
- ethane, NIST SRD 69, 100–280 K; and
- propane, NIST SRD 69 auxiliary reference-equation values, 100–360 K.

The direct Glos 2004 propane experiment is intentionally not in this packet;
it remains reproduction-analysis evidence. The same rule applies to the Figiel,
Khudaida, Crespo, Ascani, and other paper/figure-specific target matrices.
Those records may be promoted later only through a separate source and
ownership review.

Older Khudaida packet versions are immutable historical admissions and must not
be rewritten in place. `data/packet-locks.json` identifies the current packet
version for new analysis materializations.

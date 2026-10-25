# Cameretti, Sadowski and Mollerup (2005): aqueous electrolytes with ePC-SAFT

Reproduction of Ind. Eng. Chem. Res. 44, 3355 and its correction notice (44, 8944) with
the published parameters, solved natively by the Engine (issue #387). Sources, Zotero keys
and PDF hashes are in `inputs/source/source-manifest.csv`. The record with every number is
[notebook.qmd](notebook.qmd).

## Result

- **Vapor pressures (Figure 2) are partly reproduced.** The closest agreement comes from the
  article's original ion table with an ionic hard-sphere diameter equal to σ_j: 248 of 280
  points sampled from the paper's recovered model lines lie within the pre-set
  max(1 %, two pixels) tolerance, with mean deviation 1.30 % and largest 10.00 %. Original/σ
  is the choice for comparing with the recovered original curves; it does not adopt those
  parameters, which the correction declares incorrect, for general calculations.
- **The corrected ion table of the correction notice agrees less closely** under either
  diameter convention (σ_j or 0.88σ_j).
- **Densities (Figures 3 and 4) are not reproduced.** The same model gives liquid densities
  up to about 25 % below the paper's lines. A factor 1 + (ν − 1)x_salt, the error from
  multiplying an ion-counted molar density by a salt-basis molar mass, brings eight of
  nine recovered salts within 2.5 %. This explanation is inferred, not established, and is
  not an Engine density correction.
- **Table 5 (mixed salts) is not reproduced** in either printing. Against its experiments the
  Engine mean |Δp| is 0.74 %, 2.53 % and 2.55 % for original/σ, original/0.88σ and corrected/σ,
  closer than the printed model column (3.11 %); corrected/0.88σ (4.25 %) is not. Read in units
  of 100 kPa, 22 of 50 corrected Table 5 model pressures exceed the vapor pressure of the
  pure-water model, which no measured pressure and no Engine value does.
- **Pure water** meets the paper's claim (average deviation below 0.9 % from IAPWS for
  vapor pressure and liquid density, 278–393 K).
- **The Engine default start solves every bubble-pressure state** of all four variants, as does
  the declared ideal-solution start (`results/default-start-summary.csv`, by molality band; any
  failure would be listed in `results/default-start-failures.csv`). The record before #391 failed
  for 169–327 of 450 states per variant (record commit `cd2b3bc7`).
- **Mean ionic activity coefficients at 25 °C** (not reported by the paper) are poor for the
  Li and Br salts under every variant, as Held, Cameretti and Sadowski (2008) also state.

## Conventions

Debye–Hückel term of Eq. 13 with closest approach a_j = σ_j; dispersion only between
water and water and between water and ions (same-sign ion pairs excluded, cation–anion
k_ij = 1); solvent-only water permittivity from Floriano and Nascimento (2004) at 100 kPa,
interpolated by a cubic through their four tabulated temperatures; no Born term; pure-water
vapor; molality per kg water; x_salt on a salt-plus-water basis. The ionic hard-sphere
diameter is not printed, so d_j = σ_j and d_j = 0.88σ_j (Held et al. 2014, Eqs. 21–22) are
both run, each with the original and the corrected ion table.

## Files and rerun

- `scripts/model.py` declares the parameters and conventions (single owner), citing each input
  to its own source (original or corrected Table 2, the article's Table 1, Floriano–Nascimento).
- `scripts/digitize.py` recovers model lines and experimental crosses from the article
  figures; QA overlays are `figures/figure-0N/qa-overlay.png`.
- `scripts/water_reference.py` evaluates the IAPWS saturation equations (Wagner and Pruss 1993).
- `scripts/run.py` runs every native solve in parallel and writes `results/*.csv` and
  `results/run-identity.json` (Engine commit, solve count, wall clock). Each bubble pressure is
  solved from the Engine default start and from the declared start.
- `inputs/source/table-005-corrected.csv` transcribes the corrected Table 5; the original
  tables stay under `results/tables/`.

```bash
.venv/bin/python analyses/reproduction/2005-cameretti/scripts/digitize.py
cd analyses/reproduction/2005-cameretti/scripts && ../../../../.venv/bin/python run.py
```

The template parameter files of the source packet (`inputs/parameters/`, with uncorrected
σ_j, unrelated solvents and placeholder permittivity and Born fields) were removed;
`scripts/model.py` replaces them. Git history keeps them.

## Limits

Figure values come from raster images and cover a subset of the paper's data; the sulfate
vapor-pressure panel and several dense density lines were not recovered. The permittivity
interpolation is a declared choice. No Engine source was changed. Bubble pressures use a
declared ideal-solution start.
Table 4 comparisons use the original printing. Not yet independently reviewed.

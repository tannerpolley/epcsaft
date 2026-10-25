# Chloride brine ion refit

Na⁺, Mg²⁺ and Ca²⁺ refit jointly on binary and mixture osmotic data in the
Held 2014 and Figiel 2025 forms
([Engine issue 575](https://github.com/tannerpolley/ePC-SAFT/issues/575)).
Neither form meets the held-out rule for the NaCl–CaCl₂ and NaCl–MgCl₂ water
properties at the limits of the
[mixed-chloride study](../2026-mixed-chloride-water-activity/), with
Dinane–Mounir judged by source agreement (owner amendment below). With the
Saad 1980 calcium mixtures in training, the Figiel form passes both
Robinson–Bower buckets and Dinane–Mounir above ionic strength 3; it fails only
Dinane–Mounir I ≤ 3 in the joint printed-rounding realization `r1408`. The
owner qualified that Figiel-form refit on 2026-10-16 with `r1408` as its
stated limit; the Held form stays unqualified. The
[notebook](notebook.qmd) holds the numbers, figure and limits.

- [Frozen fit declaration](inputs/fit-declaration.json), including its dated
  amendments; `results/stage-a/fit-declaration.json` is the version that run
  used.
- Owner amendment, 2026-10-16 (`held_out.amendment` in the declaration):
  Dinane–Mounir is judged by source agreement. Above I = 3 the training
  Rard–Miller source, converted to water activity at the Dinane–Mounir
  compositions, disagrees with Dinane–Mounir by more than the frozen 0.0013
  RMS limit, so no model trained on Rard–Miller can meet it. For each
  Dinane–Mounir bucket and frozen realization the limit is the larger of the
  frozen limit and the RMS and maximum source disagreement on the same rows
  and realization (`results/source-agreement-limits.csv`). `scripts/evaluate.py`
  applies it to the existing optima without refitting; every other criterion
  is unchanged.
- Owner qualification, 2026-10-16 (`held_out.owner_qualification` in the
  declaration): `figiel-2025-chloride-mixture-refit/1` is qualified for
  NaCl–MgCl₂–CaCl₂ chloride brines at 25 °C over the training and held-out
  ionic strengths, with the `r1408` failure as its limit
  (`results/owner-qualification-table.md`). The computed verdicts in
  `selection.json` are unchanged.
- `results/stage-a/`: binaries plus Rard–Miller; `results/stage-a/run-60s/`
  keeps the 60 s run's F files and fit records (its H optima equal the 900 s
  ones, kept once); `results/stage-a-calcium/`: with Saad 1980, the run the
  candidate packets store.
- Packets: [held-2014-strategy2/2](../../data/packets/held-2014-strategy2/2/packet.toml)
  (settled conventions), the owner-qualified
  [figiel-2025-chloride-mixture-refit/1](../../data/packets/figiel-2025-chloride-mixture-refit/1/packet.toml),
  the unqualified
  [held-2014-chloride-mixture-refit/1](../../data/packets/held-2014-chloride-mixture-refit/1/packet.toml),
  and the [Saad 1980 observations](../../data/packets/saad-padova-marcus-1980-nacl-kcl-cacl2-observations/1/packet.toml).
- `inputs/binary-only-transfer.csv` transcribes the readiness-review
  binary-only transfer probe (`build/issue-575/readiness-r1-probes/`, hashes in
  `inputs/binary-only-transfer.sha256`) as a comparator.

Reproduce from the repository root:

```bash
OMP_NUM_THREADS=1 nice -n 10 .venv/bin/python analyses/2026-chloride-brine-refit/scripts/fit.py
OMP_NUM_THREADS=1 nice -n 10 .venv/bin/python analyses/2026-chloride-brine-refit/scripts/evaluate.py
OMP_NUM_THREADS=1 nice -n 10 .venv/bin/python analyses/2026-chloride-brine-refit/scripts/evaluate.py --run stage-a
.venv/bin/python analyses/2026-chloride-brine-refit/scripts/consistency.py
.venv/bin/python analyses/2026-chloride-brine-refit/scripts/present.py
```

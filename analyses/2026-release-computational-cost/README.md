# Release computational cost

Issue [#681](https://github.com/tannerpolley/ePC-SAFT/issues/681) remeasures the
named release `main-8298c6029ebb8b467929badc4a223c580ec591ab` against FeOs and
Clapeyron. This leaf reuses the persistent workers, timer and interleaving loop
from the audit timing campaign, with release declarations and strict quiet-host
admission before each block. Its earlier timing values are not inputs.

The [notebook](notebook.qmd) reports the ratios, absolute times, spread,
agreement and limits. Measurements and numerical agreement checks belong to
`results/`; `statistics.csv`, `ratios.csv`, `checks.csv` and `host.csv` own the
displayed quantities. `results/raw/` preserves all blocks, including the
disturbed saturation block and its complete retry. `deferred-host.jsonl`
preserves load and compiler observations while admission was deferred.

The owner's task response sets admission at one-minute load below four, with
no compiler processes from other tasks before or after each block. All serial
workers use one numerical thread and nice level ten. The package-to-package
controls and paired ranges describe the observed timing spread.

This is a frozen release record. Repeating calculation requires a fresh
detached checkout of the named tag in `build/issue-681/release`, prepared with
its README setup and `EPCSAFT_BUILD_JOBS=6`, a newly frozen image using the
existing campaign's freeze method, and the existing comparator environments.
The guarded workers refuse rebuilds or changed image bytes. The measurement
record is not an instruction to overwrite these retained outputs.

Presentation reads retained results only:

```bash
nice -n 10 .venv/bin/python analyses/2026-release-computational-cost/scripts/present.py
quarto render analyses/2026-release-computational-cost/notebook.qmd --no-execute
```

Independent delivered review and owner promotion remain outstanding.
No manuscript is changed here.

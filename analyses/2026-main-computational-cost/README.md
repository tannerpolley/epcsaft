# Computational cost on current main

This analysis measures warm serial calls after fresh numerical agreement with
FeOs 0.10.1 and Clapeyron.jl 0.6.29. It retains the same physical declarations
as the frozen release comparison and adds a separately prepared propane
saturation call when both comparators match.

The accepted measurement design is issue 929's lead-amended design, with
independent readiness r2. Scripts reuse the existing persistent workers and
native calculations. The Engine image is an ordinary detached clone with its
own environment under `build/issue-929/main`; `results/image.json` identifies
its exact commit and installed bytes. No timing from the release leaf enters
this analysis.

Admission requires physical core 15/31 reserved from user programs, CPU 31
worker affinity, nice 10 and no compiler before or after a block. Load is
recorded as host context. The reservation wrapper restores user slices on exit.
A disturbed operation is retried in full, including both members of an action
comparison; the reducer selects the latest complete operation attempt.
Numerical failures are excluded under the unchanged quantity-specific bounds.

Run `scripts/run.py preflight` through `scripts/reserve-core.sh`, then
`scripts/verify.py`, then `scripts/run.py measure --name measurement-1` through
the same wrapper, using the isolated image's Python. Freeze preparation bytes
with `scripts/provenance.py freeze` before preflight and recheck them afterward.
Reduction and `scripts/present.py` read retained files; presentation never
imports the Engine. Render `notebook.qmd` with `quarto render --no-execute`.

The notebook and retained CSV files report actual selected samples, exclusions,
absolute times and paired ratios. They support numerical verification and a
host-specific hardware measurement, not physical validation or equivalent
global phase-search coverage. Delivered review and investigator promotion
precede manuscript use.

The full fast and affected-area gates pass. The analysis-layout checker permits
only the exact design-required input snapshot paths and hashes. Identical derived
preflight checks and normalized values are stored once, with both attempts mapped
to the retained files and hashes in `results/identical-preflight-copies.json`.
The distinct raw preflights, initial preparation and excluded raw attempt remain.

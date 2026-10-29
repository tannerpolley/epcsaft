# Equilibrium numerical-error estimates

[The notebook](notebook.qmd) records numerical determination and comparison
coverage for the retained pure-fluid and neutral-mixture endpoints, the
source-reference chemical and first-thermal estimates, falsifiers and delivery
checks. Results are numerical verification, with independent pure-fluid and fixed-input mixture density references;
independent delivered review and investigator promotion remain outstanding.

`results/` owns the displayed values and judgments; `inputs/` owns the frozen
physical endpoints and original residual criteria. `scripts/present.py`
regenerates displays without running a model. `scripts/compare_854.py`
replaces the retained pure-fluid comparison reader's local judgment with the
shared native comparator; it accepts the original evidence directory and a
separate output directory. Runtime estimates and all numerical judgments
remain native.

Run `scripts/qualify.py --out <build-directory>` with this checkout's prepared
Python to compile the retained native probes and regenerate endpoint estimates,
manufactured perturbations, comparison boundaries, density-reference coverage
and row-error checks. It uses the exact recorded compiler flags and native
library; the output belongs under `build/`. The original installed/development
corpus is retained under the main checkout's `build/issue-853/build-r1/`.
Its original physical quantities and criteria were compared by the same
production reader with native-lifted coordinates added as diagnostic operands.

The source-reference additions use the accepted Route A design for issue
[#869](https://github.com/tannerpolley/ePC-SAFT/issues/869). Their final native
probe sources, exact compiler commands, logs and timing controls are retained
in the main checkout at `build/issue-869/build-r1/`; their numerical tables and
source hashes are retained here in `results/source-*`. The notebook states the
local estimate's underestimation limits separately from unavailable results.

The solve-observation qualification for [#862](https://github.com/tannerpolley/ePC-SAFT/issues/862)
is retained in `results/solve-*`. Its native/public cause table and numeric
records compare the recovered association, callback, polish, discovery,
continuation, failed Ipopt point and capture-log boundaries against the detached
baseline. Exact probe inputs, source, compile commands and full logs are retained
in the main checkout at `build/issue-862/build-r9/`; `solve-provenance.json`
records their hashes and loaded-library identities. `scripts/present.py` reads
these retained values only. This reporting qualification introduces no new
physical accuracy or validation claim.

Diagnostic retention uses ordinary copies. Out-of-memory failure during record
construction or retention propagates normally; allocation-failure controls and
incomplete-record fallbacks are outside this qualification.

The selected-endpoint mixture refinement for [#864](https://github.com/tannerpolley/ePC-SAFT/issues/864)
is retained in `results/mixture-*`, with the original binary64 declarations in
`inputs/mixture-refinement.txt`. The production association reference isolates
the existing contribution at the declared four-site mathematical input and
compares its amount/volume derivatives with binary128 differentiation. Native
probe sources, instrumentation, compile commands, trace rows and full delivery
logs remain in the main checkout's `build/issue-864/build-r3/`. The notebook
states the unmeasured complete-call cost and the remaining review boundary.

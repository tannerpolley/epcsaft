# Lithium chloride refit on the qualified chloride basis

Stage A of [Engine issue 640](https://github.com/tannerpolley/ePC-SAFT/issues/640)
fits Li⁺ on `figiel-2025-chloride-mixture-refit/1` with every other parameter fixed.
The fit repairs the trace-LiCl bound and Lithium's activity targets, but fails
the Hamer–Wu guard and the held-out Yao mixtures. The unchanged KCl–CaCl₂
baseline also fails its held-out check. Stage A created no parameter packet.
The [notebook](notebook.qmd) holds all values and the
pre-registered diagnosis and stage B trigger.

The [declaration](inputs/fit-declaration.json), its digest and freeze time
precede the fit. Source observation packets retain the checked transcriptions;
the declaration lists their hashes, source PDF hashes and exact row locators.
`results/fit.json` reports every start and both distinct optima. The selected
optimum is the lowest-cost converged start; only one optimum falls in the
declared cost band. `results/parameter-changes.json` records the Li-only changes.

## Frozen first Stage-A fit (#720)

`scripts/fit.py` is **historical, not rerunnable** for the first Stage-A fit.
Record snapshot: commit `ca719012c31e8b41edcb98ae82a388117b58d343`;
`results/fit.json`, its declaration and all cited tables stay frozen.
At #696’s measured commit `335707d59c67cf1f92df12b5d25530bf31d916a6`
(report line 711), it raised
`AssertionError: retain the first frozen run; no implicit reruns` because
`results/fit.json` already exists. That result was preserved rather than
removed to force execution. Evaluation reads the retained fit; the later
cation-mixing qualification for packet `/2` remains separate from this freeze.
The commands below evaluate retained parameters without repeating that fit.

```bash
env -u EPCSAFT_HSL_LIBRARY OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 10 ~/.local/share/mamba/envs/lithium-extraction-reaktoro/bin/python analyses/2026-lithium-chloride-refit/scripts/pitzer.py
env -u EPCSAFT_HSL_LIBRARY OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 10 .venv/bin/python analyses/2026-lithium-chloride-refit/scripts/evaluate.py
```

The report-only trace LiCl supplement in MgCl₂ and CaCl₂ reads the selected
fit and unchanged baseline without fitting or replacing qualification scores:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 ~/.local/share/mamba/envs/lithium-extraction-reaktoro/bin/python analyses/2026-lithium-chloride-refit/scripts/pitzer.py --trace-divalent
env -u EPCSAFT_HSL_LIBRARY nice -n 10 .venv/bin/python analyses/2026-lithium-chloride-refit/scripts/trace_activity.py --divalent
nice -n 10 .venv/bin/python analyses/2026-lithium-chloride-refit/scripts/present.py
quarto render analyses/2026-lithium-chloride-refit/notebook.qmd
```

Run these supplement commands before presentation when reproducing all
results. `results/trace-divalent.csv` retains the eight states and the Song
Li–Mg comparator; the accompanying JSON files retain the post-fit evaluation
identity and Pitzer coefficients. Li–Ca has no comparator. None of these
report-only values changes the frozen declaration, limits or packet decision.

Presentation reads retained values and does not import the Engine. The
readiness design and stage zero are retained in the main checkout's
`build/issue-640/`; stage A execution evidence is in `build/issue-640/stageA/`.
Independent delivered review remains with the lead.

## Joint Li⁺/Na⁺/Mg²⁺ experiment (stage B)

The owner's 2026-10-19 “Re-weight + stage B” decision authorizes the joint
cation experiment in #640. Ca²⁺, K⁺, Cl⁻ and water retain the `/1` coefficients.
The [stage B declaration](inputs/stage-b-declaration.json) owns the expanded
training rows, seven fixed discrepancy scales, ordered starts, numerical
limits and source identities. The stage A declaration and results above stay
intact. The experiment acknowledges that its held-out sources had already
been inspected in stage A and diagnosis.

`scripts/stage_b.py` reuses the #575 training and complete water-property
qualification readers and the stage A Li and K readers. Its `freeze` action
ran before the first fit; the digest and UTC time are retained beside the
declaration. Its `run` action performs the native fit and required evaluations
under the shared fitting and total elapsed budgets, refusing to replace an
existing fit. Reproduction requires the declared source and native-library
identities; a new experimental declaration is a separate authorized action.

```bash
env -u EPCSAFT_HSL_LIBRARY OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 timeout --signal=TERM 2400 nice -n 10 .venv/bin/python -u analyses/2026-lithium-chloride-refit/scripts/stage_b.py run
nice -n 10 .venv/bin/python analyses/2026-lithium-chloride-refit/scripts/stage_b_present.py
```

[Stage B results](results/stage-b/evaluation.json) retain every close optimum's
verdict. Each optimum directory retains training and held-out row predictions,
rounding metrics, activity transfers and all eight divalent trace predictions.
The complete #556 comparison outputs and execution logs are retained under
`build/issue-640/stageB/`. A budget stop is incomplete numerical evidence and
does not establish physical failure or permit a new parameter packet.

The joint fit fails the old chloride rule, both Yao buckets, divalent transfer
preservation and the high-ionic-strength Li–Na reference residual. No new
parameter packet is created. The original executed source is retained in
`inputs/stage-b-executed.py.txt`; the post-fit NaCl descriptive-metadata
correction and numerical verification are linked from the notebook. The
current calculation script fixes the cation-list scope and comparison-module
identity without changing a native fit observation or repeating the fit.

## Additional neutral components (#781)

The #781 comparison is retained as it ran for that PR:
- `scripts/neutral_components_inputs.py` writes `build/issue-781/packet-inputs.hpp` from packet `/2`, plus zero-amount dodecane from `gross-sadowski-2001/3`.
- `scripts/neutral_components.cpp`, compiled with `tools/native-driver` against that header, prints the 113 value and derivative comparisons and writes the occupied-neutral rows.

Its outputs are kept as `results/cation-mixing/neutral-derivatives.csv` and `neutral-occupied.csv`. `scripts/cation_mixing_present.py` writes the notebook paragraph from them.

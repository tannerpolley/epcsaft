# Bülow et al. (2021) sour-gas model comparison

**Historical results below.** Record snapshot:
commit `ca719012c31e8b41edcb98ae82a388117b58d343`. Original wheel identities, source
inputs and retained numerical evidence keep their recorded provenance. The pre-D042
IPOPT/`regression.fit` execution runner was removed; the retained partial reproduction
is not a current replay. The current Figure 1 subset is described below.

This campaign is retained as historical evidence from the pre-D042
Regression lifecycle. Its former IPOPT/`regression.fit` execution runner has
been removed; the recorded inputs, results, and results remain immutable
evidence and are not a current validation command.

This campaign is the source-faithful reproduction target for **Predicting
Vapor–Liquid Equilibria for Sour-Gas Absorption in Aqueous Mixtures of
Chemical and Physical Solvents or Ionic Liquids with ePC-SAFT** (Bülow,
Ince, Hirohama, Sadowski, and Held, *Ind. Eng. Chem. Res.* 2021, 60,
6327–6336, DOI `10.1021/acs.iecr.1c00176`).

The paper declares one reacting liquid and one vapor phase. The vapor support
is restricted to water, MDEA, sulfolane, CO2, and H2S; charged species remain
in the liquid. This campaign therefore uses a declared phase topology only.
It does not run phase-candidate discovery or make a globality claim;
composition remains constrained to the declared supports.

## Implemented sequence

1. Bind the article and 45-page supporting information by SHA-256 and retain
   the exact reaction constants, pure parameters, interaction coefficients,
   and SI ARDs.
2. Replay the neutral sulfolane–CO2 and sulfolane–water binaries using the
   declared liquid/incipient-vapor topology and approximate Figure 1/3
   digitizations.
3. Retain the former sulfolane–CO2 Regression diagnostic and its generated
   result as provenance; the pre-D042 generic Regression runner is no longer
   a live executable route.
4. Assemble the eight-species aqueous MDEA–CO2 reaction system with the
   paper's original Born term, solvent-only permittivity, liquid-only ions,
   and neutral water/MDEA/CO2 vapor support.
5. Evaluate one positive-loading Kuranov state at each of 313, 333, 373, and
   393 K through the generic declared-topology GREPE bubble solve and retain
   the exact plotted table and liquid/vapor states.

This is a **retained historical partial reproduction**, not a claim that the
complete paper has been reproduced. Figures 4–10 do not publish their numerical point
arrays, so the reactive checks use adjacent Kuranov Table 2 data at about
19 wt% MDEA. The article Figure 1/3 arrays are explicit approximate
digitizations. Source-reported ARDs remain separate from generated results.

## Retained inputs

- `source/paper-result-summary.csv`: exact headline conditions and ARDs
  transcribed from the article captions and text.
- `source/reaction-network.csv`: the five liquid reactions and phase support
  stated by the article.
- `parameters/paper-parameters.csv`: Tables 1–3 values transcribed from the
  article; these are not silently converted into a current Engine packet.
- `results/data/bulow-2021-binary-plotted-data.csv`: exact neutral-binary
  paper-parameter predictions used by the generated figure.
- `results/data/bulow-2021-sulfolane-co2-fitted-data.csv`: exact forward
  pressure replay at the generic Regression result used by the generated figure.
- `results/data/bulow-2021-reactive-mdea-co2-plotted-data.csv`: exact coupled
  reactive predictions, phase compositions, and densities used by the figure.
- `results/rendered/bulow-2021-sour-gas-reproduction.{png,svg,pdf}`: generated
  model-vs-data comparison.

## Interpretation

The aqueous reactive route was evaluated with the source-shaped original-Born
model. All four reactive states evaluated, with a generated mean total-pressure
ARD of 43.07% (individual values 34.13–48.75%). This is close to the scale of
the paper's reported CO2 partial-pressure errors but is not the same dataset or
measurement basis, so it is diagnostic model-reproduction evidence rather
than a validation or a regressed MEA parameter set.

The sulfolane–water binary evaluated for all 31 digitized points with 14.74%
mean ARD. The published sulfolane–CO2 coefficient evaluated 18 of 24 points
but gave 395.42% mean ARD. The former fixed-pressure GREPE Regression run
admitted all 24 sulfolane–CO2 rows and recorded an owner result for every
fitted target. Its recorded Ceres cost, bound activity, and forward-pressure
replay remain truthful diagnostic evidence, not a current executable result,
reproduction claim, or promotable parameter. A missing source-formulation
detail is still indicated; loosening numerical tolerances would not repair it.

The existing `analyses/reproduction/2021-bulow/` campaign is a different paper
(Bülow, Ascani, and Held, *ePC-SAFT advanced – Part II*, salt solubility,
2021). It is not reused for this sour-gas campaign.

Current source limits for #789 are recorded in
`inputs/source/reference-limits.md`. Source article/SI identities remain
available, but the approximate figure arrays lack calibrated uncertainty and
Kuranov's total-pressure rows are not the paper's CO2 partial-pressure rows.
The retired interface is not an Engine capability diagnosis. Historical numerical results remain unchanged; the current subset uses the accepted calibrated references.

A calibrated eight-point Figure 1 model-line subset is retained in
`inputs/source/figure-1/`, with a source image, pixel transform, CSV, uncertainty
metadata and QA overlay. This is a reference preparation result. The current producer below uses these accepted references; historical output retains its separate identity.

## Current Figure 1 subset

The live producer covers only the eight calibrated sulfolane/CO2 model-line
points and four source-precision probes per point. It uses the literal printed
parameters and unchanged pressure allowance; every mismatch is retained.
[comparison.csv](results/comparison.csv),
[source-precision-probes.csv](results/source-precision-probes.csv) and
[diagnostics.json](results/diagnostics.json) contain actual native values.
The [computed notebook](notebook.qmd) distinguishes source agreement from
numerical closure. Other figures and all old results above remain historical
or uncompleted; the old aggregate errors are not current evidence.
Independent delivered review is pending.

From the repository root with this checkout's editable Engine:

```bash
env -u EPCSAFT_HSL_LIBRARY nice -n 10 .venv/bin/python analyses/reproduction/2021-bulow-sour-gas/scripts/reproduce.py
cd analyses/reproduction/2021-bulow-sour-gas
QUARTO_PYTHON=../../../.venv/bin/python nice -n 10 quarto render notebook.qmd
mkdir -p ../../../build/issue-789/2021-bulow-sour-gas/render
mv notebook.html ../../../build/issue-789/2021-bulow-sour-gas/render/
```

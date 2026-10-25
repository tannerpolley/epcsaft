# Native density diagnostics: delivered review

Independent delivered review of `b3e9b070` against integrated main `8999d6a2`
on 2026-10-02: **Supported; no remaining required findings**.

## Capability and numerical evidence

Properties' existing density owner reports pressure crossings of either mechanical
sign, fixed-composition spinodals ($P_\rho=0$), connected declared cells,
unavailable samples and finite-evaluation brackets. Pressure evaluation, exact
slope, adaptive slope search, packing admission and ordinary stable closure
remain shared. Negative-slope diagnostic roots do not supply ordinary pressure
inverse actions. Python declares inputs and reads the native results.

The independent reviewer reproduced these comparisons on the installed,
non-editable wheel:

| Quantity | Achieved | Required limit |
|---|---:|---:|
| Cold-vapor maximum density error | $1.4188117348\times10^{-10}$ mol/m³ | $10^{-9}$ mol/m³ |
| Spinodal $\lvert P_\rho\rvert/(RT)$ | $1.3100631691\times10^{-14}$ | $10^{-10}$ |
| Fine pressure-density action relative error | $3.2292210383\times10^{-6}$ | $5\times10^{-6}$ |
| Coarse pressure-density action relative error | $1.2917561274\times10^{-5}$ | $2\times10^{-5}$ |
| Reciprocal action relative error | $0$ | $10^{-12}$ |

Cold-vapor pressure-root counts are 1/3/3 and spinodal counts 2/2/2. Ethane
control root counts are 1/3/3/0. Near-coalescence shifted declared-grid counts
remain 3/1/1/3, distinct from the adaptive counts 3/3/3/3. The finite-domain
bracket maximum relative width is $2.13\times10^{-16}$ against $10^{-13}$.
Maximum dilute pressure and slope errors are $3.77\times10^{-10}$ and
$7.54\times10^{-10}$ against $10^{-6}$, decreasing with density.

- [Cold-vapor rows](../../../analyses/2026-cold-vapor-branch/results/pressure-roots.csv)
  retain the independent reference densities, residuals, signs and brackets.
- [Ethane controls](../../../analyses/2026-density-closure-diagnosis/results/root-controls.csv),
  [shifted grids](../../../analyses/2026-density-closure-diagnosis/results/scan-grids.csv),
  [pressure actions](../../../analyses/2026-density-closure-diagnosis/results/pressure-actions.csv),
  [before/after comparison](../../../analyses/2026-density-closure-diagnosis/results/before-after.csv)
  and [provenance](../../../analyses/2026-density-closure-diagnosis/results/replay-summary.json)
  retain the inputs, historical values and final comparisons.
- [Ladder refinement](../../../analyses/2026-density-closure-diagnosis/results/ladder-refinement.csv)
  diagnoses amplification of liquid root-bracket uncertainty when subtracting
  nearby densities. Tighter verification-only brackets correct that error;
  acceptance limits remain unchanged. The relative metric uses the expected
  action magnitude, not a dimensionful unit floor.

Qualified wheel SHA-256:
`79fca1a03770a9feeb0a2ff0c48f1bb9e45e29935309d14625b1149077adaec5`.
Source provenance is `72b3c80b`; `b3e9b070` adds the regenerated results, not
runtime changes. `check-all.sh` and `check-affected.sh HEAD` passed, including
114 Python tests, six native suites, architecture, mypy and Ruff. The reviewer
inspected the execution evidence and performed no-write numerical probes.

## Required corrections checked

The initial review identified three omissions, corrected in `6bb46d73`:

1. Resolved roots now accumulate rather than being replaced by an incomplete
   adaptive pass. Crossing ownership uses connected bracket enclosure or exact
   numerical zeros, never approximate density distance. Three ethane roots
   survive budgets 650 and 780; an exact 100 mol/m³ node root survives budget 3.
2. Shared adaptive slope-search failure retains its active bracket and cause.
   Budget 5 retains the 5500–7000 mol/m³ bracket at 308.8 K and $10^7$ Pa.
3. Control and grid tables retain the cold-closure cause at 200 K/$10^{10}$ Pa:
   `no stable density root (count 0)`, separate from geometric domain exit.

Sorting applies on partial returns. Thirty additional bounded ownership probes
found no over-counts. Unavailable association samples retain their causes and
break connectivity; unstable anchors supply no ordinary pressure actions.

## Growth, integration and limits

Runtime growth is **402 gross, 40 removed, +362 net**, total **21707**.
Properties is 4763 (294/40/+254), Bindings 1922 (36/0/+36), and Python including
stubs 2707 (72/0/+72); other layers are unchanged. The independent reviewer
justified the two-line gross excess over the original 400-line ceiling under
the owner-directed #248 pause: required result/bracket preservation adds 17 net
native lines. No numeric ceilings changed. The shared pressure/slope owners
remain single implementations; batch 12's Python solver loops and seven
allowlist entries are retired. [Code-budget closeout](../code-budgets.md)
retains per-layer counts and the scientific reason for growth.

Integration commit `72b3c80b` removes four byte-identical batch-7 table copies,
retaining their canonical numerical rows once through the existing retained-file
mapping. The reviewer verified both aliases, canonical LF hashes and original
collected CRLF hashes. This changes no numerical value or process/cancellation
qualification.

Qualification covers the stated finite searches and reference cases. It is not
exhaustive root enumeration, global stability or physical validation. Ordinary
closure retains its documented one-loop and finite-search limitations. Neutral
packing/finite-boundary equality does not qualify an electrolyte boundary.
Batch 15 retains notebook/prose consolidation responsibility; hosted issue
closeout and publication remain with the accountable lead.

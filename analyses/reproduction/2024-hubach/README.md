# Hubach 2024 reduced S1/S6 comparison

Current scope: four nonzero-salt native S1/S6 rows, all numerically admitted.
Retained values and direct checks are in `results/subsystem-comparison.json`;
`results/subsystem-run-identity.json` records input, driver and library hashes.
The notebook computes separate source differences and mean absolute errors.
They are descriptive under the declared benchmark choices; the causes of
remaining source discrepancies are undetermined. Delivered review is pending.

From the repository root, using this checkout's editable Engine:

```bash
nice -n 10 agent-heavy bash analyses/reproduction/2024-hubach/scripts/reproduce.sh
cd analyses/reproduction/2024-hubach
QUARTO_PYTHON="$OLDPWD/.venv/bin/python" nice -n 10 quarto render notebook.qmd --to html
mkdir -p ../../../build/issue-789/2024-hubach/render
mv notebook.html ../../../build/issue-789/2024-hubach/render/
```

The comparison covers the four nonzero-salt rows of SI Tables S1/S6.
It uses Lithium's benchmark assumptions: neutral `[emim][tcb]`, 1 kg aqueous
solution plus 1 kg pure solvent, and two declared liquids at 294.15 K and
101300 Pa. The tables state 100000 Pa; the pressure difference and unknown
source feed ratio preclude exact source parity. All assumptions and their
Lithium file paths, commit and hashes are retained in
[the declaration inputs](inputs/lithium-declarations.json).

Lithium's benchmark assumption retains the Engine dispersion rule
$u_{ij}=\sqrt{u_i u_j}(1-k_{ij})$, with the printed coefficients unchanged.
Hubach Eq. (6) puts $(1-k_{ij})$ inside the square root (PDF CBJZ36VC p. 3,
article p. 2257). This difference also precludes source parity.

S1 compares an aqueous emim-cation-equivalent mass fraction calculated from
the neutral pseudo-component; it does not introduce a charged emim species.
S6 compares the organic TBP mole fraction. Measured and paper-calculated
columns remain separate descriptive references. The two nominal zero-salt
rows are retained but excluded; no trace-salt replacement is proposed.

Full extraction remains **unavailable: charged exchange species and the
source feed ratio are not decided in Lithium** ([Lithium #237](https://github.com/tannerpolley/Lithium_Extraction/issues/237#issuecomment-6355298159)).
S11 is retained source evidence only. This comparison does not fit parameters.

The review draft is `build/issue-789/2024-hubach/design-lithium.md` in the main
checkout. The sole required readiness correction, the dispersion disclosure, is retained
above and in the declaration; execution follows the task authorization.

The earlier numerical files, where present, remain frozen historical evidence
at their recorded Engine identity. The source assessment in
`inputs/source/reference-limits.md` predates the Lithium declaration above;
its source-parity limits remain. Neither retained source rows nor this
assessment establish numerical verification or physical validation.

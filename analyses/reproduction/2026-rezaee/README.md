# Rezaee 2026 density-only comparison

Current scope: seven S5 pure-DES and eleven S6 DES–TOPO density rows.
All rows pass the fixed numerical solution checks. Density residuals and
separate table AARDs remain descriptive; their causes are undetermined.
Retained rows, direct state diagnostics, exact thermal actions and local
uncertainty scales are in `results/density-comparison.json`. The notebook
computes all displayed values and its discrete model-versus-measurement
figure from retained files. Input, driver and loaded-library identities
are in `results/density-run-identity.json`. Delivered review is pending.

From the repository root, using this checkout's editable Engine:

```bash
nice -n 10 agent-heavy bash analyses/reproduction/2026-rezaee/scripts/reproduce.sh
cd analyses/reproduction/2026-rezaee
QUARTO_PYTHON="$OLDPWD/.venv/bin/python" nice -n 10 quarto render notebook.qmd --to html
mkdir -p ../../../build/issue-789/2026-rezaee/render
mv notebook.html ../../../build/issue-789/2026-rezaee/render/
```

The comparison uses seven literal SI S5 pure-DES rows and eleven
SI S6 rows for 90 wt% DES + 10 wt% TOPO, at Lithium's benchmark pressure of
101325 Pa. Only the published Table 4 DES and Table 6 TOPO tuples are consumed.
[The declaration inputs](inputs/lithium-declarations.json) retain Lithium's
benchmark assumptions, exact paths, commit and hashes: the rounded F/3 mass
0.20748 kg/mol, TOPO mass 0.3866349 kg/mol, DES donor/acceptor self association,
one self-paired TOPO site, no cross association and density-stage dispersion
`k_ij=0`. These are Lithium's benchmark assumptions, not recovered author defaults.

The S5 272.04–330.37 K and S6 293.15–343.15 K grids stay separate. Which of
the article's range or S5 is erroneous is not decided in Lithium. Measurement
uncertainties remain 0.01 K and 0.1 kg/m³; they are separate from numerical
solution checks. No curve, temperature or parameter is fitted to the targets.

Extraction remains **unavailable: the published A/B reference conversion is
not decided in Lithium** ([Lithium #237](https://github.com/tannerpolley/Lithium_Extraction/issues/237#issuecomment-6355298159)). The extraction-fitted
DES–TOPO `k_ij=0.0623`, complexes, aqueous strategies and Lithium's fitted
reaction coordinates are excluded from the density calculation.

Readiness is recorded in `build/issue-789/2026-rezaee/readiness-lithium-r1.md`.
The native liquid request picks the highest stable root in its finite search;
where only one stable root exists the Engine reports the actual branch as
`Unspecified`. This does not establish global stability. Local propagated
uncertainties assume independent temperature and density errors as a Lithium
benchmark assumption; input/topology/model uncertainties remain unquantified.
The paper fitted Tait-smoothed data, so this literal-SI comparison is neither
exact optimizer replay nor independent physical validation.

The earlier numerical files, where present, remain frozen historical evidence
at their recorded Engine identity. The source assessment in
`inputs/source/reference-limits.md` predates the Lithium declaration above;
its source-parity limits remain. Neither retained source rows nor this
assessment establish numerical verification or physical validation.

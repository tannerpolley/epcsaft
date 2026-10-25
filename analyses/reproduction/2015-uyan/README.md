# Uyan et al. 2015 explicit-ion MDEA ePC-SAFT

This is the source-located transfer-reconstruction lane for Uyan et al., *Fluid Phase
Equilibria* 393 (2015) 91--100, DOI `10.1016/j.fluid.2015.02.026`.

The paper declares one finite reacting liquid and one vapor containing only
the neutral volatile species water, CO2, and MDEA. The historical campaign passed
that topology directly to its then-installed equilibrium interface. It does not run
phase-candidate discovery or make a global-stability claim.

## Historical evidence (does not run on the current wheel)

The former generator and its retained 30-row results were recorded in
commit `566d57de5ec6b2bd36d1e1f0460b62a837e70e0e`. No wheel SHA-256 was
recorded with those rows; this exact record commit is the available execution
provenance, not a verified wheel identity. That generator used the retired `epcsaft.unit_registry`,
`SourceReferenceDeclaration`, `source_reference_transfer` and
`GeneralReactiveEquilibriumProblem` / `equilibrium.solve` interfaces. Git retains
its method identity. The current generator uses native declarations and writes
a separately identified source-ledger result directory.
Current reacting MEA qualification is owned by #61 and
[`2026-mea-recovered-state-actions`](../../2026-mea-recovered-state-actions/README.md).

The local Uyan source contains the four reactions, pure and binary parameter
tables, reaction-correlation coefficients, source ranges, figure captions,
and reported aggregate ARD values. Exact experimental observations were
materialized separately from Ermatchkov et al. 2006 Tables 1--3. The retained
source CSV contains all 101 rows and their reported uncertainties; the
executable target is the 30-row lane with CO2-free solvent MDEA mass fraction
between 0.18 and 0.20, corresponding to Uyan Figure 6.

The installed public wheel evaluated all 30 target rows with a finite reacting
liquid and incipient neutral vapor. No solver or physics tolerance was
weakened. The replay used eight campaign workers and completed in 315.95 s;
one row used the campaign's next declared pressure seed. Every result had a
denser liquid than vapor. The exact plotted table and the SVG, PNG, and PDF are
retained under `results/`.

The retained reconstruction’s CO2-partial-pressure ARD is 118.98% (180.07% at
313.15 K, 79.14% at 353.15 K, and 97.72% at 393.15 K). It demonstrated that wheel’s
specified-topology numerical path; it did not reproduce Uyan’s model accuracy.
The reconstructed temperature trend is much too weak at 353.15 and 393.15 K. Uyan's 19.22% Figure 5 value and 19.82% abstract
value remain provenance only.

The remaining model-fidelity gap is explicit. Uyan did not publish author
code, solver details, or an unambiguous water--MDEAH+ interaction statement;
the article gives conflicting zero/inherited values. This campaign follows
Table 5 and the VLE-results text and inherits the water--MDEA interaction. It
also reconstructs the paper's Hsieh pure-MDEA permittivity dependency from
the six official ThermoML values and retains Uyan's extrapolation beyond their
measurement range. Water and MDEA are declared independently as dielectric
solvents, their correlations are evaluated at the solve temperature, and CO2
is excluded from the solvent mixing rule. Uyan does not publish the evaluated
Floriano--Nascimento water correlation, so the campaign explicitly labels its
Archer--Wang water polynomial as a transfer approximation. The historical installed
Engine used its revised Debye--Huckel and source-reference
implementation, so this is a transparent
model transfer rather than an exact author-code replay.

## Declared-start conversion diagnosis

The #238 candidate at `624881f1` changed the retired reference path to current
neutral EOS reference declarations. That reconstruction did **not** preserve
this historical lane. The exact 30-row, 0.18–0.20 mass-fraction replay returned
26 numerical roots and four failed rows; 189 of 238 historical comparisons
missed. Maximum relative pressure/density error was 59.589566888229804 against
1e-7, and maximum absolute vapor-composition error was 0.5851705927164594
against 1e-7. Accepted roots closed original rows to at most
4.440892098500626e-14 against 1e-10. Small equation residuals therefore do not
establish historical agreement. See [comparisons](results/data/native/comparisons.csv)
and [diagnosis](results/data/native/comparison-summary.json).

An initial unfiltered 101-row run returned 81 roots and 20 failures, with the
same 189 misses on the historical subset; its [table and summary](results/data/native/diagnosis-101/)
are diagnosis evidence, not acceptance. The bounded
[initializer comparison](results/data/native/reference-path-probe.json) gives
CO2 partial pressure 2.42482383379317 kPa with ordinary initialization and
2.4248238337931682 kPa with one declared start for t1-r01, versus the historical
0.7893425273754909 kPa. This rules out the ordered-start controller as the cause
of that row's discrepancy under the current reference declaration; it does not
identify every historical model/reference difference. The old transfer,
reference-declaration and chemical-standard interfaces are absent from the
installed wheel, SHA-256
`84dbd5d062ccff3e6b539c4fe1e6498f54002c67823dc8cc15baadf253233bbe`.

The earlier starts-only diagnosis remains historical evidence; it does not
define the current source convention. The source ledger states water at its
pure-component reference and all other species at infinite dilution. Uyan's
text accompanying Eqs. 22–23 states the same, while Eq. 18 labels the CO2
coefficient with a pure-component superscript. The historical generator's
default pure-CO2-vapor lane therefore must not be attributed unambiguously to
the paper. The retained table's own lane labels remain unchanged; the historical
wheel was not identified. #274 added per-species pure-gas references (PR #313),
but no value of that default pure-gas lane is replayed exactly.

## Source-ledger calculation

`generate_data.py` declares `NeutralReference(solvent="water", ...)` with a
charge-neutral dilution path and `RawMoleFractionInfiniteDilution` for all four
source reaction correlations. It retains the source mole-fraction constants,
the 100 kPa reference pressure, physical parameters and feed declarations.
The polynomial's reference-temperature shift changes its coordinates only,
preserving $\ln K=c_1+c_2/T+c_3\ln T+c_4T$. No molality conversion, Python
activity-scale inversion, or pure-gas CO2 approximation is used.

`solve_equilibrium_starts` owns the row-specific pressure start followed by the
de-duplicated 100/300/500 kPa nominal list, including selection and termination causes.
Each physical pressure guess is the larger of its nominal list entry and the
observed CO2 partial pressure plus 100 kPa; retained selected seeds report
that declared physical guess, with nominal entries recorded separately.
For the four literal 313.15 K source rows t1-r09 through t1-r12, each pressure
start also receives the owner-approved fixed
[t1-r08 composition and packing guesses](inputs/starts/313-15-k-t1-r08.json).
The pressure order, equations, feeds, work limits and tolerances are unchanged.
The [26-root pressure-only replay](results/data/native/source-ledger/pressure-only/uyan-2015-ermatchkov-2006-grepe-summary.json)
retains its separate identity. The fixed declarations add no Python result
selection or continuation.
Native `PartialPressure` observation evaluation supplies CO2 partial pressure
after the selected result declares the reporting problem's physical guesses.
Python reads outputs and computes comparison statistics. The original-row
criterion remains $10^{-10}$; historical pressure/density and composition
comparison criteria remain relative and absolute $10^{-7}$, respectively.

The calculation command for the 30-row lane is:

```bash
.venv/bin/python analyses/reproduction/2015-uyan/scripts/generate_data.py --minimum-mdea-mass-fraction .18 --maximum-mdea-mass-fraction .20 --workers 1
```

Numerical runs use the checkout's own non-editable Engine
wheel. The reference wheel and earlier retained results keep their own identities.

The historical pure-CO2-vapor declaration is now expressible through
`NeutralReference.species_reference`. Its separate diagnostic uses raw source
constants and the same native starts:

```bash
.venv/bin/python analyses/reproduction/2015-uyan/scripts/qualify_species_reference.py
.venv/bin/python analyses/reproduction/2015-uyan/scripts/present_species_reference.py
```

The [mixed-reference results](results/data/native/species-reference/summary.json)
retain the actual historical-table misses. The
[reference qualification](results/data/native/species-reference/qualification-summary.json)
retains native pure-state identities, the decreasing-pressure limit, refined
derivatives, default-map reduction and unchanged key-solve work. This expresses
the generator at `e12a5ff8`; it makes no historical-equivalence claim. The paper's
documented all-solute infinite-dilution calculation remains the default.

The [notebook](notebook.qmd) reports the historical/current calculation differences and
qualification. The final replay returns all 30 roots with maximum original-row
residual $5.062616992290714\times10^{-14}$ against $10^{-10}$. CO2 partial-pressure
average relative deviation is 698.100321% on the same 30 feeds, versus the
historical 118.976509%; the numerical roots do not reproduce the authors'
accuracy. Current [values](results/data/native/source-ledger/uyan-2015-ermatchkov-2006-grepe-plotted-data.csv),
[attempts](results/data/native/source-ledger/attempts.csv),
[original rows](results/data/native/source-ledger/rows.csv),
[reference diagnostics](results/data/native/source-ledger/references.csv) and
[all before/after quantities](results/data/native/source-ledger/reference-convention-comparison.csv)
retain their separate identity. Historical results and rendered figures remain
under their original paths. Changed values report the difference between historical reconstruction and
current source-ledger calculation. The [notebook](notebook.qmd) states the
historical-identity limit. This is numerical verification, not physical
validation or a reproduction of the authors' accuracy.

Compare and render the retained historical and source-ledger tables without Engine evaluation:

```bash
.venv/bin/python analyses/reproduction/2015-uyan/scripts/compare.py
.venv/bin/python analyses/reproduction/2015-uyan/scripts/render.py
```

No private Engine import is used.

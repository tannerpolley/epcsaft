# Finite-addition enthalpy fitting

The signed closed-system enthalpy increment per mole of added feed is
qualified on the manufactured neutral and ionic cases in
[issue #210](https://github.com/tannerpolley/ePC-SAFT/issues/210).
[The notebook](notebook.qmd) interprets the retained values; this is numerical
verification, with no calorimetric data, physical validation or parameter adoption.

`results/comparisons.csv` owns the F1–F8 comparisons, refinement values,
last-three mixed-slope ladder values, admission changes, failure messages and
Jacobian timings. `results/summary.json` identifies their wheel and inputs.
`results/regularity.csv` and `results/provenance.json` retain the earlier
endpoint regularity calculation, including its original free-pressure guess
of 1.5 MPa. The native check uses that same guess; a draft 3 MPa guess converged
to coincident phases and was removed. The F4 check uses the six fitted
coordinates of F3, including tied epsilon/k, rather than an extra individual
species slot. No tolerance changed.

Inputs and executable numerical checks are owned by
`engine/tests/native/equilibrium_s1_s2_checks.cpp`,
`engine/tests/native/regression_reaction_checks.cpp` and
`engine/tests/equilibrium/test_reference_calorics.py`.
After preparing the current installed wheel and native programs with the
repository environment command, the native `equilibrium_checks S2` selector,
one `regression_checks` call each for `addition-enthalpy`, `ionic-addition`,
`inventory-face` and `neutral-electrostatics`, and the installed-wheel
`test_reference_parameter_enthalpy_actions_match_resolves` check emit the
comparison rows. `scripts/retain.py` imports their logs.
`scripts/plot.py` reads retained comparisons and renders the figure without
running the thermodynamic model.

The #231 extension retains its exact zero-inventory face, neutral electrostatic
restriction and fresh-solvent heat/actions in `results/zero-inventory.csv`,
with separate wheel/input/source provenance and original failure diagnostics in
`results/zero-inventory.json`. The frozen Engine declarations live in
`engine/data/reference/mea-zero-inventory/`; `scripts/zero_inventory.py`
exports those declarations (`--native-input <file>`), and `--native-run <log>`
compiles the native readout `scripts/zero_inventory.cpp` with `tools/native-driver`
and runs it on them. Python does not calculate this extension's thermodynamics or derivative
verification. The notebook records the numerical domain and remaining delivered
review; these results are not calorimetric validation or
application parameter adoption.

The #269 incipient liquid/vapor extension keeps the #231 tables unchanged.
`results/incipient-wheel.csv` owns the installed-wheel boundary and positive
inventory states; `results/incipient-native.csv` owns their numerical checks
and the separate #231 replay. The paired JSON files identify source, wheels,
inputs, endpoint rank/conditioning and refusal causes. `scripts/incipient.py`
only declares cases and reads installed native results; the native
`inventory-boundary` selector and `scripts/zero_inventory.cpp` qualify the solved
endpoints. `scripts/retain.py --zero-inventory --stem incipient-native --base
8f5f5685` keeps these new rows separate from historical results.
The `incipient-cp`, `incipient-provisional` and `incipient-cost` tables retain
#261/#265 preservation and the existing #238 speed controls, with exact
supplied compositions in `inputs/incipient-guesses.json`. Their JSON provenance
records producer identities and the inherited unaccepted flash limits.

The F5 historical trace-Lithium replay retains 1080 native fraction-action rows
in `results/lithium-trace-actions.csv` and 38 native endpoint rows in
`results/lithium-trace-endpoints.csv`. `results/lithium-trace.json` owns their
wheel/source lineage, units, criteria, raw-output hashes and limits. The paired
current/original endpoint readouts agree exactly; action differences remain
observations without a new comparison margin. This is reconstructed numerical
preservation evidence for the literal historical scalar reaction targets.

`inputs/lithium-trace/` retains the byte-exact ASL record. The
`probe.py`/`campaign_probe.py` declaration/readout snapshots, which ran only on
their identified wheel, were removed by #341 (record commit `21d7a112`). Original native JSON/CSV readouts are losslessly compressed under
`build/issue-231/lithium-trace/`, with raw differences, log and restoration
records identified in the compact provenance. `scripts/plot.py --lithium-trace`
reads the retained endpoint table without rerunning thermodynamics.

## Named-release rerun (#680)

The [research notebook](notebook.qmd) reads the new results under
[results/main-8298c602/](results/main-8298c602/). The earlier retained files remain
under their original identities. The new provenance records the downloaded
wheel, its verified hash, the actual installed Engine identity and the declared
producer. Unavailable portions retain their stopping states and causes.
Independent delivered review and owner promotion remain open; the manuscript
is unchanged.

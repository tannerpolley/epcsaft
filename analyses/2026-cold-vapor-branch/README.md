# Cold-start incipient vapor density

The #204 failure comes from a dense-root seed for the incipient vapor, not
from a narrow admission window. At the exact manufactured methane/ethane
pseudo-mixture input, the shared initializer now recovers the retained J1
state at 180 K and 1,596,859.1113027716 Pa. Its maximum scaled equilibrium
residual is 1.8122e-15 against the existing 1e-10 criterion. This is numerical
verification of the declared bubble state, not physical validation of an
isomerizing methane/ethane fluid or a global stability calculation.

[The notebook](notebook.qmd) interprets [densities.csv](results/densities.csv),
[pressure-roots.csv](results/pressure-roots.csv) and
[summary.json](results/summary.json). The exact inputs and original failure
remain in [cold-start.json](../2026-reaction-constant-fitting/results/cold-start.json).
That original record is the single input and J1 reference for the new
previous-failure check in `engine/tests/equilibrium/test_cold_vapor.py`.

At the equimolar 180 K feed, the vapor spinodal maximum is
1,595,392.9545692347 Pa, 1,466.15673354 Pa below the requested pressure.
Only the dense stable pressure root remains: 19,698.5002963 mol/m³, packing
0.367749292482. The density owner correctly returns it without a phase label.
The initializer incorrectly used that same root for both fugacity estimates,
making every K equal to one and seeding identical phases. The rejected
corrector reached 201.101220717 K and vapor-labelled density
15,112.6344208 mol/m³, packing 0.253523129240. The declared vapor root there
is 1,143.15865409 mol/m³, packing 0.0191771435167, with window
[1,131.84025158, 1,154.59024063] mol/m³. The log displacement is 2.58173593042
against 0.00995033085317. The dense root has positive pressure slope, but
it is the wrong declared phase; positive slope and small equilibrium residual
do not make it an admissible vapor. Admission correctly publishes no state.

`start_compositions` now uses the ideal-vapor fugacity limit, ln phi = 0,
when the feed search has one unlabelled root. This replaces its self-ratio;
actual phase seeds still use shared EOS pressure closure at the estimated
compositions. Caller-guessed starts and two-root estimates retain their rules.
The corrected cold start seeds vapor at 1,282.74550558 mol/m³ and converges
to J1 vapor density 1,288.45507053 mol/m³, packing 0.0203906816053. No density
window, numerical tolerance, equilibrium equation or admission rule changed.
Properties and bindings have zero line growth; equilibrium grows by three
physical runtime lines, from 6,018 to 6,021 (hard ceiling 7,000; soft ceiling
6,000). The three lines prevent a sole-root self-ratio without another solver.

Build and install the sole non-editable wheel before replaying:

```bash
bash tools/build-engine.sh --install-wheel
.venv/bin/python analyses/2026-cold-vapor-branch/scripts/run.py
.venv/bin/python -m pytest -q -s engine/tests/equilibrium/test_cold_vapor.py
engine/tools/run-native-tests.sh 'equilibrium/*' boundary
bash tools/run-slice.sh khudaida
bash tools/check-all.sh
```

The native diagnostic invokes the existing internal corrector with the former
K=1 composition seeds, inspecting coordinates before public branch admission.
It is not a second equilibrium method or a published physical result. Its
temperature and compositions reconstruct the original rejected iterate to the
reported digits. The pressure-root probe uses the installed wheel for pressure
and analytic slope through native density diagnostics with declared log-density
cells and native refinement; packing comes from native readouts. The
[delivered comparison](../../docs/scientific/reviews/issue-239-delivered-review.md)
retains the historical root and spinodal comparison criteria. It distinguishes the
negative-slope middle crossing from the two stable crossings. Its sampled
root list is not proof of global completeness. The Engine density admission
retains the finite-search limits in ALG-DENSITY.

The new numerical falsifier failed with the original branch-window message on
the baseline wheel and passes on the repaired wheel. The explicit native
equilibrium and boundary runs pass with all 7,197 and 3,813 unique comparison
ID/model values exactly unchanged; their counts and common hashes are in
[verification.json](results/verification.json). The Khudaida fixed-settler candidates pass qualification against
the retained reference, with no retained values edited. The full fast gate and
independent review are recorded in the final delivery section of the notebook.
These checks bound preservation to their exercised states. No physical data,
application result, critical-point passage or phase discovery was qualified.

## Current replay and refused equimolar seed (#720)

The dense-corrector endpoint above is historical. The #719 diagnosis on
`ca719012` found that the explicit [0.5, 0.5] phase seeds now return
`phases_coincide` with no admitted state; the exact first revision that changed
that endpoint was not bisected. Reading the absent density vectors caused the
probe's segmentation fault. The probe now checks success, optional temperature
and pressure, and vector sizes before printing final values. Its relocated
executable uses the current shared native-library directory; shared linkage
replaced archive-only linkage under #485.

The [current results](results/issue-720/summary.json) retain the refusal,
[seed and admitted densities](results/issue-720/densities.csv), and
[pressure-root comparisons](results/issue-720/pressure-roots.csv). No final
physical values are emitted for the refused seed. The public cold start still
recovers J1 at the original tolerances, and every pressure-root assertion passes.
The earlier dense endpoint remains in the original results; it is not a new
current-corrector result.

```bash
env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 nice -n 10 \
  .venv/bin/python analyses/2026-cold-vapor-branch/scripts/run.py \
  --output-dir build/issue-720/cold-vapor
```

The notebook retains its exact base content. Its historical source-line counts
have no retained CSV/JSON source in this leaf, so the #720 repair notes and
current-result links remain here rather than extending that notebook.

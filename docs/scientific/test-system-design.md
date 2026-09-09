# Test system and live companion design

Design reviewed Ready, 2026-09-08. This specifies the required capability. The
implementation plan owns delivery status and scientific qualification; the
[acceptance record](reviews/test-system-acceptance-2026-09-09.md) records the
implementation checks and their limits. No thermodynamic equation,
accepted tolerance, phase specification or numerical qualification changes here.

## Required behavior

A test is declared once beside its executable definition. Running it records its
own outcome and useful numerical observations. The Tests page renders those
records, with no test-specific numerical parser, copied input, equation catalogue
or editable result in the presentation. A new registered test appears without a
frontend change. Planned checks are visible without pretending they ran.

The guarantee is **no unverified result presented as current**. There is a bounded
observation delay, not instantaneous execution: an active browser checks every two
seconds; request timeout is five seconds, and current badges expire after ten
seconds without a successful coherent refresh. These are presentation defaults,
not scientific tolerances. On return from a suspended browser, invalidate badges
before refreshing. An unavailable service or broken inventory cannot look green.
Expensive studies need not run on every edit; they must visibly become outdated.

Use one local Quarto website: Map, Plan and Tests. The companion is a read-only
viewer. Rendering, refreshing and starting the service never run scientific
calculations or tests; results come from the normal build, check and study
commands. Automatic execution on source edits was retired on 2026-09-09 (see
[the cleanup brief](project-map-design.md)). No Docker requirement, database,
remote service, CI pipeline or new testing framework is needed. Use the existing pytest, native executables, wheel
preparation, standard-library files and local preview.

The installed local service must have persistent startup configuration rather
than a transient session-only unit. Enable startup for the user's service manager;
boot operation before login requires lingering. Validate installed/enabled state
and service restart separately from an actual computer reboot, and state when
the latter has not been exercised.

## Ownership and structure

| Existing owner | Target responsibility |
|---|---|
| `engine/tests/` | Executable numerical checks, colocated declarations and small reporting helpers |
| Native test support and executable registration | Discoverable native test IDs, metadata, selection and individual outcomes |
| pytest collection and reporting hooks | Collected Python cases, parameters, skip/error/outcome records and numerical observations |
| `analyses/` | Explicitly run studies; same result format where a useful table exists |
| `data/` | Reference measurements, source-model inputs and scientific conventions |
| `tools/run-numerical-evidence.py` | Shared execution/recording entry point, replacing its F1/F5 special cases |
| `tools/check-all.sh`, `tools/check-affected.sh`, slice runners | Invoke that recording path; never create a second untracked result path |
| `tools/numerical-inventory.py` | Consume runner discovery, planned declarations and completed records; replace source-name guessing |
| `tools/companion-data.py` | Validate and project useful scientific fields; never assign an inferred pass |
| `docs/companion/` | Generic Quarto view, sorting, search, filtering and formatting |
| `build/numerical-evidence/` | Rebuildable compact run manifests and homogeneous numerical CSV/TSV tables |
| Implementation plan | Capability requirements, delivery order and reviewed qualification; links to test IDs when concrete |

Keep these existing owners and the capability-based directory organization; kinds
are metadata, not separate numerical/package directory trees. A shared record reader/writer helper is justified by
native and Python producers; it is tooling, not an Engine runtime interface.
Internal fingerprints and error logs remain private to the machinery. Hashes,
raw diagnostics and secrets do not appear in the notebook. Source paths are
useful navigation in the Map view and milestone details; scientific result
tables stay focused on quantities and their meanings.

## What a test means

Distinguish an execution test from a comparison. A single executable test can
produce several comparisons and many state points. Its declaration supplies the
test identity; each comparison supplies a stable local key and scientific meaning.
One notebook row represents one comparison or one check without numeric output.
Parameter cases and derivative step ladders expand inside that row. Never split a
single state into repeated cards, or discard additional quantities to fit a label.

Use independent attributes rather than one overloaded status:

| Attribute | Values and applicability |
|---|---|
| Kind | Numerical; Package; Tooling. Numerical is the default scientific work. Package covers an existing installed-wheel/binding/integration check with a concrete purpose. Tooling covers existing discovery, rendering and readiness checks. These labels do not authorize more tests. |
| Category | Core, Contributions, Properties, Equilibrium, Regression, Surface, Bindings/Python, or Applications; a study additionally identifies itself as a study |
| Property group | Ideal, Residual, Total, only for relevant quantities |
| Comparison basis | Relationship, Experimental, Matching-model, Cross-model sanity; absent when irrelevant |
| Method | Identity, derivative, limit, conservation, equilibrium observation or another concise scientific method; only when useful |
| Coordinates | T, rho, composition direction, parameter direction or equilibrium coordinates, when relevant |
| Execution policy | Automatic fast or Explicit; existing slow/expensive markers and studies are Explicit |
| Lifecycle | Planned or Implemented |
| Execution | Not run, Queued, Running, Passed, Failed, Skipped, Error or Interrupted |
| Freshness | Current, Outdated or Unverified; independent of execution outcome |
| Evidence | Numerical rows, Outcome only, or Unavailable after an interrupted/invalid record |
| Judgment | Criterion satisfied, Criterion failed, or Observation only; applies to a comparison, not automatically the whole test |

Compatibility is not a default kind or a required test layer. An actual approved
compatibility obligation can later be a Package check with that scope. Do not
preserve legacy behavior or introduce API/snapshot/coverage tests to fill labels.
Physical-input admission checks already present are classified honestly; this
design neither silently removes them nor declares them thermodynamic validation.

A Package check can emit numerical observations. Kind describes its purpose;
evidence describes its output. An automatic fast test can belong to any category.
A solver converging is an execution fact, not agreement with an experimental value.

## Declarations and numerical reporting

Python uses one registered pytest marker with keyword metadata; native tests use
one small registration entry beside their function. Registration drives native listing and
selection as well as execution; do not maintain a second test list for HTML.
Existing pytest parametrization remains responsible for executable Python cases.

Required implemented-test metadata: stable ID, human scientific title, kind,
category and execution policy. Add group, basis and comparison metadata only where
applicable. Every emitted scientific numerical comparison must declare its basis
and concise method/claim; package/tooling observations need not claim a physical basis. Native and Python declarations normalize to the same record shape.
A missing declaration remains visible as **Unclassified** and makes inventory
validation incomplete; do not silently omit it or guess from its filename.
Duplicate IDs or duplicate parameter-case identities are inventory errors.

At the actual comparison, one reporting call accepts the computed model value,
reference value, actual input object, units, comparison key and criterion (or explicitly observation-only mode). It writes
the observation and evaluates that same criterion for the test assertion. Do not
have one tolerance for the assertion and another for the notebook. Ordinary pytest
outcomes without a reporting call still appear honestly as Outcome only. Numerical
checks intended to support scientific qualification must emit their relevant values.

The test/reference helper owns the relationship and its LaTeX label together;
the frontend receives the label. Inputs are captured from the actual invocation,
including resolved state, composition, branch, parameters and reference convention
as scientifically needed. No independent hand-entered T/rho/x display fields.
If a test internally sweeps or transforms inputs, emit the effective inputs for
each resulting observation rather than repeating the top-level request.

For a scalar comparison, retain at least reference, model, unit and criterion.
Store input/parameter vectors and multidimensional observations with explicit
coordinates, not only their largest error. Homogeneous sweeps belong once in CSV
or TSV; JSON contains compact metadata, artifact identities and outcomes. The reader
may derive display errors from the retained operands. Preserve full precision;
format scientific notation and smaller parenthesized units only in HTML.

The criterion is explicit: absolute/relative tolerance and its exact rule, norm,
conservation bound, derivative-ladder requirement, or an observation-only judgment.
Keep existing criteria during conversion. No universal tolerance, silently adopted
threshold, or inference from printed precision. For zero references, relative error
is undefined unless a scientific scale was explicitly declared; display a dash and
retain absolute error. Nonfinite values, domain exclusion and missing values are
explicit states, never coerced to zero. A criterion-free comparison is **Observed**,
never Passed. Remove the ambiguous **Compared** label.

A planned case has a lightweight declaration in its future test/study owner with
lifecycle Planned and no callable. Native listing and Python metadata collection
include these declarations, but neither executes a fake skipped test. Required
planned metadata is only ID, purpose/category and the concrete unmet prerequisite.
Do not require hypothetical inputs or three comparison types. High-level capability
requirements remain in the plan; concrete declarations reference that requirement
instead of duplicating its status. Implementing a case replaces its planned entry
under the same ID, leaving one visible row. FEOS/Clapeyron remain planned with their
scientific slices; this work neither installs nor runs them.

### Concrete record contract

Use three small versioned records; do not store separate frontend DTOs on disk.
A generated inventory is disposable and never an editable source of declarations.

| Record | Required fields | Conditional fields |
|---|---|---|
| Declaration | schema_version, test_id, title, kind, category, lifecycle, execution_policy, evidence_expectation | property_group, requirement_id, unmet_prerequisite for planned work; expected comparison keys for numerical evidence |
| Run | schema_version, run_id, test_id, parameter_case_id, input_generation, artifact_identity, execution, started_at, finished_at and artifact digests for terminal states | outcome_reason and comparison records; unavailable identities/artifacts use an explicit unavailable state, never an omitted field |
| Comparison | comparison_id, quantity, unit, inputs, judgment, evidence_state; basis and method/claim for scientific numerical evidence | relationship_latex, coordinates, criterion, observations artifact, observation-only reason, unavailable_reason, internal reference identity |

Numerical observation rows contain case/coordinate values, model and reference
operands, plus per-row criterion outcome where defined. Uncertainty, scales and
norms are retained where the scientific criterion needs them. A comparison's
criterion applies uniformly to its table or names explicit per-row criteria;
never aggregate away failing points. Inputs can refer to one shared case record
inside the run, avoiding repeated composition/parameter vectors in every CSV row.
The input record used to invoke the calculation is serialized once by the reporter.

Python authoring surface: a `pytest.mark.evidence(...)` declaration and
an `evidence` fixture. The fixture offers `compare(..., criterion=...)` for
a judged comparison and `observe(...)` for a finite observation with no acceptance
claim. `observe(...)` cannot produce Passed. Custom scientific criteria are
evaluated by the test helper and passed to `compare(...)` with their recorded
operands, rule and verdict; simple scalar criteria can be evaluated by the reporter. Native registration and reporting expose
the equivalent two operations through existing test support. Keep custom derivative
or equilibrium criteria in their test helpers; do not build an expression language
or move scientific algorithms into the reporter. Helpers provide a human equation
label and serialize their actual thresholds. A fixture/hook failure is a recording
error, not permission to publish a passed test without its expected observations.

An execution may pass all its assertions while containing observation-only
comparisons. Conversely, a failed numerical criterion must fail its owning test.
A recorded comparison failure paired with a passed owner is an invalid record and
must be exposed as Error. A test that declares numerical comparisons but emits
none cannot be accepted as complete numerical evidence; preserve its execution
outcome separately and mark its evidence unavailable. Classify a deliberate
outcome-only check explicitly instead of inferring a numerical claim.

### Identity, criteria and completeness invariants

`test_id` is unique across the checkout's supported inventory. Each collected
parameter case has one `parameter_case_id`; use the explicit value `single` for an
unparameterized test. A run records exactly one collected pair. Parameter cases
must use stable declared keys or canonical serialized parameter contents, never
collection order. Multiple comparisons use unique keys within that pair; a sweep
uses distinct state/coordinate keys within its comparison table. Duplicate records
are errors. Removed tests disappear from the current inventory; orphaned retained
results cannot become current entries. Planned declarations do not produce runs.

`input_generation` is the digest of a canonical sorted manifest of pre-build
source/input paths and content hashes, declarations, data, execution configuration,
lockfiles and expected environment/toolchain configuration. It excludes generated
wheels, binaries and generated inventory artifacts. Recompute it before/after
build and before/after execution using one shared input identity routine.

`artifact_identity` identifies the wheel, native binaries, actual compiler and
native/Python dependencies, provider versions and generated collected inventory
used by the run. These are recorded after build/collection and verified before
and after execution; they are not recursively included in `input_generation`.
The build/collection record explicitly binds its artifacts to the pre-build input
generation, which must still match after building. Preserve both manifests
internally to explain invalidation; neither belongs in the frontend.

At build/run time a case is executed only when its artifact/environment identity
verifies before and after execution. At display time a run is shown Current only
when its input generation still matches the checkout and the recorded wheel and
binaries are still the files on disk; the compiler, runtime and dependency probes
are not repeated on refresh. Display freshness therefore means "sources and
artifacts unchanged since this run", not a re-verified environment. Every terminal run requires both identity fields, its case key and artifact
digests. A build failure or unidentifiable dependency records an explicit unavailable
identity/artifact state and cannot qualify as Current evidence. Failed, skipped,
interrupted, error and incomplete records require concise reasons. Unavailable
identity does not hide the observed execution failure; it limits its applicability.

A scalar criterion records its nonnegative finite `atol`, `rtol`, scale and combining
rule explicitly. For example, the existing Python symmetric closeness rule uses
`abs(model-reference) <= max(atol, rtol*max(abs(model),abs(reference)))`.
Other existing rules retain their exact definition: do not translate a native
normalized error, summed tolerance or derivative-ladder rule into that example.
A custom criterion records its helper identity, exact human mathematical rule,
operands, scale/norm, thresholds, tolerance rationale and verdict. The reporter
validates completeness and propagates that verdict; it does not interpret a new
expression language. Changing a criterion invalidates the prior judgment.

`evidence_expectation` is either Outcome only or Numerical. Numerical declarations
name the expected comparison keys; parameterized/swept helpers generate these keys
from their actual case definition during collection/setup, before evaluation.
Adaptive studies declare required named outputs and a coverage/termination rule;
they need not predict every internal solver iterate. A failure cannot retroactively
shrink the expected evidence to the rows that happened to survive.
Reconcile expected and emitted keys and state points on completion. Missing output
has `evidence_state=Unavailable`, `judgment=Unavailable` and an `unavailable_reason`;
it never becomes Observed merely because an artifact exists. Nonfinite operands,
excluded domains, missing ideal/reference inputs, provider or reporter errors are
represented explicitly. Partial valid rows may remain inspectable but do not
satisfy the missing comparison. Only finite, complete criterion-free observations
can be Observed. A prior complete run cannot fill gaps in a failed replacement.

Experimental, Matching-model and Cross-model records retain the source/provider
identity and version, source units, domain and state/convention basis, plus supplied
uncertainty or the explicit fact that none was supplied. Missing uncertainty does
not forbid a useful observation, but cannot be silently interpreted as zero.
Matching-model records also identify the compared equations, parameters and
conventions and record whether matching was established. A mismatch cannot retain
a matching-model acceptance claim. Reference lineage stays internal; notebook text
contains only useful science. Old same-engine values and timings cannot be relabelled
as experimental or independent matching-model validation. Do not introduce historical
snapshots as acceptance references merely to populate a table.

## Scientific comparison rules

Relationships check the stated identity or derivative chain within this model;
they do not establish experimental accuracy. Perturbations use the same admitted
branch and composition chart. Composition derivatives use explicit zero-sum
directions. Retain perturbation size, directional coordinates and actual state
points. Do not force T, rho and x tests where the quantity or available evidence
does not justify them.

Experimental comparisons use the observable actually measured. Compare total
caloric quantities with consistent ideal contributions and reference conventions;
do not invent experimental residual quantities. Unavailable ideal caloric inputs
block that comparison alone. Derived observables such as isothermal compressibility
must transform model and reference consistently with explicit density basis and
units; singular transformations must report their domain limitation.

Matching-model comparisons require matched equations, parameters and conventions.
FEOS or Clapeyron may serve this role where their implemented model matches ours.
Different EOS or activity-model equilibrium outputs are Cross-model sanity checks:
record discrepancies without making exact agreement a gate. Adopt a quantitative
criterion only when justified for that use. No model is selected by application
name and no phase-analysis requirement is added to a case that does not need one.
Reference notes may contain concise scientific conditions/uncertainty; source
lineage remains in retained data, excluded from the notebook as requested.

## Execution and freshness

1. Discover declarations and planned cases; build an inventory generation. Python
   uses actual pytest collection, including parametrized cases and collection
   errors. Native listing uses the registered executable built from current inputs.
   If a build or collection fails, show Inventory unavailable/incomplete and retain
   previous entries as Unverified. Never claim the old inventory is exhaustive.
2. Changed inputs make recorded results Outdated at the next refresh. Additions,
   removals and renames count. A presentation edit only rebuilds the view; it
   does not rerun thermodynamics.
3. A single shared lock covers check-all, affected, focused and slice entry
   points from wheel preparation through publication; nested calls reuse the
   owning context and cannot reacquire it or start independent builds. Existing
   wheel verification is mandatory before Engine execution. The viewer never
   takes this lock to execute anything.
4. The automatic selection (`--automatic`) runs when a normal gate invokes the
   runner, including registered native fast cases as well as Python cases.
   Native registration alone currently does not make a suite part of the fast
   gate; the migration must explicitly connect its permitted selection. The automatic set is the collected implemented cases explicitly declared
   Automatic fast, intersected with the existing permitted fast selection. Slow,
   expensive, campaign and unclassified cases are excluded. Python fast candidates
   come from the existing not-slow/not-expensive selection; native candidates must
   be explicitly classified and measured before admission, never inferred from a
   suite name. The initial native conversion admits F1 only after confirming its
   bounded runtime; F5 remains an explicit diagnostic until classified against the
   actual fast budget. Explicit tests show Outdated or Not run with their scientific
   prerequisite. A failure is recorded once; retry only for a new generation or an
   explicit request. Bound runtime with declared timeouts; stop only owned processes.
5. Capture `input_generation` before build, bind the resulting `artifact_identity`
   to those inputs, and check both before and after execution. Edits during build or execution
   invalidate that run. Do not install from one source generation while labelling
   results with a later generation. A running generation can finish but cannot
   overwrite results of a newer generation.
6. Publish each completed test record and its table atomically. Write tables in a
   private run directory, then replace the manifest pointer only when validation
   succeeds. Under the shared lock, atomically replace one current manifest indexing
   all completed case records; each entry retains its own generation and can be
   labelled outdated independently. Readers never concatenate whatever F1/F5 files
   happen to exist. Readers cannot see half-written CSV/JSON pairs. Retain per-test outcome,
   not a suite outcome copied onto its children. Native registration brackets each
   test; on abort, the active check is Error/Interrupted and unstarted checks are
   Not run. Partial observations are identifiable and never a completed pass.
7. Serve plan, inventory, current applicability and records as one coherent snapshot
   with an opaque revision. Validate the relevant input identity before and after
   snapshot assembly; retry once or return Updating/Unavailable. A browser replaces
   its snapshot atomically. Never mix old numerical rows with new conditions.

Initially fingerprint the whole computational input set in `input_generation`, rather than attempting a
fragile per-function dependency graph. Include native/Python source, all tests and
shared helpers, declarations, build files, lockfiles, environment configuration,
recording code and relevant materialized data/reference contents.
Record the resulting installed wheel, binaries, actual compiler/native dependency
identities and Python package environment in `artifact_identity`. For studies
include scripts, configuration, seeds and materialized reference inputs in the
input generation; record actual external provider identities with the artifacts. An unidentifiable external dependency makes that
result Unverified. Content identity, not only modification time or Git commit,
controls applicability; uncommitted edits and input deletion must invalidate it.

Conservative invalidation can rerun more fast tests than necessary. Narrow it only
when a complete declared input closure and measured benefit justify doing so.
Recompute the input generation on every result request so a stale badge cannot
survive an edit. Never silently access changing remote reference data: materialize it through its
existing data owner first. Environmental identity is operational evidence, not a
claim to bitwise reproducibility or universal numerical correctness.

Normal checks and explicit slice invocations publish through the same reporter.
Uninstrumented direct binaries or foreign reports may be discovered as Unverified;
no machinery can retroactively attest a run that bypassed it. At migration close,
all repository-supported entry points must record outcomes. The local service
starts only the viewer.

## Notebook behavior

Retain collapsed scientific categories and collapsed rows. Keep the four visible
columns Quantity/Check, Inputs, Relationship/Comparison and Status. A small kind
label in the first column and a Kind filter distinguish scientific, package and
tooling work; tooling is collapsed and visually secondary. No extra Properties tab.

Examples of unambiguous status: **Passed**, **Failed**, **Observed**, **Outdated ·
last passed**, **Unverified · last failed**, **Running**, **Planned**, **Skipped**.
Outcome-only rows need no empty numerical table. Show their scope and concise
reason on expansion. Fresh numerical rows show the actual inputs, model/reference,
errors and applicable criterion. Old observations can be inspected on expansion
with Outdated/Unverified labels; their current-result cells remain empty. Never
substitute a previous pass for a failed or incomplete current run.

Replace Live with **Synced** (inventory/view current), **Updating**, or **Offline**.
Synced does not mean tests passed. If any recorded check is queued/running, say so. Planned and explicit work are not
synchronization failures. Keep selected scientific comparisons visible only once.
If one execution test emits several comparisons, show its aggregate execution
outcome only in context; do not duplicate a generic test row beside every comparison.
A partially emitted test also retains its interruption or omitted-evidence state.

## Delivery and deletion

Use one bounded local capability task with these ordered acceptance steps. Split
into subissues only if separate delivery or ownership actually becomes useful;
parallel trackers and a new remote are unnecessary for this design.

1. Define the shared declaration/result schema and inventory validation. Convert
   F1, F5 and one existing Python numerical check as vertical examples;
   preserve values, criteria and the known F5 failure. Current Engine collection
   has nine unparameterized tests; exercise parametrized collection in the tooling
   acceptance harness without inventing additional scientific cases. Demonstrate missing metadata,
   collection failure, multiple observables and partial-run handling.
2. Convert all actual native/Python registrations and existing package/tool checks;
   route fast/focused/affected/slice entry points through the same recorder. Add
   useful planned declarations and connect existing studies without replaying them.
   Reconcile collected IDs against the discovered inventory: each is classified,
   deliberately planned, or explicitly reported as unavailable. No silent omissions.
3. Implement identities and atomic publication. Verify edit-during-build/run,
   dependency/data changes and result corruption. Existing known failures remain
   failures.
4. Replace special-case frontend/backend projections with the shared records and
   independent statuses. Verify notebook refresh, expiry, offline behavior, search,
   math, collapsed rows and the no-provenance presentation.
5. Remove F1/F5 log parsers, native regex discovery, category-name heuristics,
   frontend property/equation catalogue, stale-status string parsing and competing
   result writers. Remove obsolete generated assets. Update current owner docs;
   qualify this infrastructure only after the acceptance evidence below passes.

Do not keep a permanent legacy importer. During conversion old results remain
Unverified; once converted, preserve useful scientific evidence with its owner and
remove the obsolete display route. No numerical acceptance is earned by conversion.
This is test infrastructure preceding scientific implementation, not completion of
Core or any thermodynamic layer. Use the task-configured implementation agents and independent reviewer. Model
settings belong in task configuration, never in test metadata.

## Acceptance evidence

Use a small end-to-end harness for this tooling behavior, separate from scientific
falsifiers. These checks verify record integrity and display truth, not thermodynamic
accuracy, and do not add a general API/coverage test program.
Declare its individual checks once beside their executable definitions, record
their own outcomes through the shared tooling route, and display them in Tests.
The installed fast gate invokes the harness; synthetic scenarios remain labelled
as such. Contributions and Properties are the supported scientific suite names;
internal directory and binary identifiers do not set presentation vocabulary.

| Exercise | Required observed behavior |
|---|---|
| Add/rename/remove a registered native or parametrized Python test | Inventory updates without editing the website; build/collection failures show incomplete inventory |
| Run one test normally and through a focused command | Same stable ID, inputs, numerical rows, criterion and outcome; no duplicate notebook row |
| Run a test with several quantities and a derivative ladder | Every comparison and coordinate is retained once; expandable rows group state points correctly |
| Reuse existing F1/F5 evidence | F1 values preserved; F5 disagreement remains visible; no criterion or success invented |
| Fail, skip, abort or time out part of a suite | Individual outcomes correct; unstarted tests never inherit Passed |
| Edit source, helper, declaration, tolerance, data or environment | Prior results cease to be current before the replacement run completes |
| Edit during build/run, corrupt an artifact or interrupt publication | No falsely current or mixed-generation result; failed replacement cannot resurrect a pass |
| Open multiple tabs, restart the viewer, repeat an unchanged failure | No test process starts; unfinished work reconciles honestly |
| Break collection, disconnect server, suspend/resume browser | No stale current badges; inventory incompleteness and offline state visible |
| Add a planned or explicit external comparison | Visible without a fake run, forced input fields or automatic installation |
| Edit only presentation or plan wording | View refreshes; numerical evidence is not rerun or silently requalified |
| Duplicate IDs, missing declarations, dynamic/parametrized collection or missing expected comparisons | Explicit inventory/evidence failure; no guessed classification or silently accepted empty numerical run |
| Concurrent gate execution, crash at manifest replacement and failed rerun after prior pass | One shared lock and coherent manifest; no mixed artifacts, duplicate work or resurrected pass |
| Criterion/owner disagreement, zero reference, nonfinite operand or missing reference identity | Invalid/unavailable evidence handled explicitly; no false Passed/Observed/Current judgment |
| Inspect compact UI and internal records | Values/units/criteria agree; provenance stays internal; no mandatory empty comparison classes |

Design review, 2026-09-08: independent Luna High scientific-basis and execution
reviews both returned Ready after corrections to reference/criterion completeness,
expected evidence, concurrent recording and the input/artifact identity split.
The lead checked these findings against the current test entry points and this
contract. The design is Ready; the implementation plan owns current delivery status.
No scientific milestone is qualified by this review.

# Reference-limit diagnosis (#750)

Decision **(a)**: the old coordinate vectors were computed before the neutral reaction-reference correction in `bd7b6e388f97f24893f74346ae31e4c35d951d04` ([PR #385](https://github.com/tannerpolley/ePC-SAFT/pull/385)). The correction returns the square-root-molality Richardson limit instead of a finite dilution rung. This is a reviewed removal of the reference bias, not evidence that the physical Born formulation has no solution. [The wheel comparison](wheel-comparison.json) evaluates exactly the same two old coordinate vectors. The retained pre-correction wheel from main `026b3031` gives row maxima 2.4868995751603507e-13 and 2.984279490192421e-13. The post-correction `e9b1aed5` wheel gives about 1.85795e-7, as do the published main wheels and current main. The release catalog begins after this correction, so the earlier pair uses explicitly hashed retained local wheels; no historical Engine was compiled.

The turning replay now uses the old points as guesses for the existing native ordered-start driver at each unchanged physical parameter. All eleven roots are corrected and meet the original 1e-10 row criterion; the worst retained row maximum is 4.82502926502093e-13. [The corrected roots](turning/turning-corrections.json), [phase values](turning/turning-verification.csv), [bounded correction probe](turning-correction-probe.json), and [rerun log](born-turning-rerun.log) retain the calculation. No old archive member was rewritten, and no solve controls, bounds, alternate starts, or tolerance were introduced.

The parameter-action verifier also passes. Its smallest finite-difference scaled errors are 8.963037445800033e-10 and 5.408773657762672e-9, within the existing 1e-5 limit. All 4500 common-point comparisons between the two retained historical/native path tables pass, with maximum error 3.33574570376877e-13. That input-table check is not a new continuation on the corrected Engine; the first-action checks use fresh native solves. The old audit's empty-comparison failure does not reproduce in this checkout, where the retained path inputs exist. The comparison now reads the retained path inputs independently of `--out`; a separate output directory no longer makes that input set empty. A current [custom-output rerun](born-parameter-custom-output.log) also verifies all common points. Provenance collection uses the loaded Engine instead of an obsolete wheel lookup. [The action values](parameter-actions.csv), [historical path comparison](../native-parameter/parameter-path-comparison.csv), and [log](born-parameter-rerun.log) retain the evidence.

## Historical chart replay outcome (#762)

**No accepted native reproduction method exists in the inspected Engine.**
`probe.py` and `verify_global_scan.py` are retained as historical fixed-coordinate
replays, with their original row comparisons intact. The original chart rows,
coordinates, plots and diagnostic comparisons keep their historical identity.
They are not current-main roots or current-main chart evidence. The chart-row
criterion remains `1e-10`; the separate retained-point comparison in `probe.py`
remains `2e-10`. Neither limit was changed.

The original [chart producer at `cbc8e0c1c`](https://github.com/tannerpolley/ePC-SAFT/blob/cbc8e0c1c/analyses/2026-born-unavailable-states/scripts/global_scan.py)
compiled the full equilibrium request, then performed the chart correction in
Python. Its `keep_rows` excludes zero-based row 16, the conserved carbon
balance. Its `correct` function fixes liquid log mole fraction 7 (hydronium) or
5 (bicarbonate), solves for the other 19 coordinates, and checks the remaining
19 native equations at `1e-10`. The Engine supplied row values and exact
Jacobians; it did not supply a native chart solve. There was no separate
carbon-omitted Engine request to recover.

Source inspection at `9541d9f7e57984c175df81b8ca91d5a827a90278` establishes the
missing native operation:

- `PhaseSpec` and `EquilibriumProblem` in `engine/native/equilibrium/problem.hpp`
  declare composition guesses and feed inventories, but no fixed liquid
  species fraction or selected omitted balance.
- `compile_problem` in `engine/native/equilibrium/request_compile.cpp` appends
  every conserved-moiety row for this reacting, nine-species declaration.
  Its special azeotrope declaration requires a neutral, nonreacting binary
  feed and cannot represent either case.
- `CompiledSystem.rows` is read-only in
  `engine/native/bindings/equilibrium_bindings.cpp`; `evaluate` and `jacobian`
  calculate the declared rows without correcting a point.
- `solve_equilibrium_starts` passes the full declaration to the one native
  driver. `continue_equilibrium` also retains those balances; neither accepts
  the omitted row and held species coordinate. The native inventory bounds
  used by current `global_scan.py` solve only the conserved-inventory linear
  problem, not the thermodynamic chart equations.

Thus direct native evaluation cannot correct the old rows after #385, while a
full-equilibrium restart imposes the carbon equation and releases the chart
coordinate. Changing the feed, chemistry, or reference convention would change
the declared problem. Restoring the Python corrector would violate the native
calculation rule. A new native held-composition corrector is not an existing
method authorized by this issue.

### Retained wheel identities

Both replays retain the [original producer provenance](../provenance.json) and
the [#750 wheel comparison](wheel-comparison.json). The comparison's
`026b3031` wheel is the pre-#385 reference, SHA-256
`94b55dfcf72f21b43010c7125d71f103f41fd56f7a6b6fe0774903eef71cc72a`.
The `e9b1aed5` wheel is the post-#385 comparison, SHA-256
`e691b5a9a2fe41e9115d83c914c912a928ed671f1b82ea001c90e6edf8e3d97b`.
The latter did not produce the pre-correction chart rows. The two historical
turning-point checks in that retained comparison give pre-correction row
maxima `2.4868995751603507e-13` and `2.984279490192421e-13`, versus
`1.857945903793734e-7` and `1.8579481775304885e-7` after the correction.
These wheel measurements are reused from #750, not new historical-wheel runs.
No historical Engine was built or installed for #762.

### Current-main reruns

The worktree's editable Engine at the inspected main commit, with
`EPCSAFT_HSL_LIBRARY` unset, was used at the original temperature, 313.15 K,
feeds, supports, reactions and physical parameters. The bounded
[comparison record](chart-replay/summary.json) identifies the loaded Engine,
input hashes, selection and the [one-time inspection command](chart-replay/inspection-command.txt);
[row differences](chart-replay/row-differences.csv)
retain every compared point. Native row evaluations are compared with stored
values; no correction, density re-solve or new row reference enters that table.
All 69 retained `probe.py` point vectors differ at more than `1e-10`, with
maximum difference `1.8579567040433176e-7`. All 1,233 retained criterion-meeting
global-scan chart vectors also differ, with maximum difference
`1.857956135609129e-7`. The largest current noncarbon row among those chart
vectors is `1.85894805326825e-7`; they fail the unchanged chart criterion.

The unchanged commands still fail their comparisons:

```bash
nice -n 10 env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 .venv/bin/python analyses/2026-born-unavailable-states/scripts/probe.py --output-dir build/issue-762/probe
nice -n 10 env -u EPCSAFT_HSL_LIBRARY EPCSAFT_BUILD_JOBS=6 .venv/bin/python analyses/2026-born-unavailable-states/scripts/verify_global_scan.py --output-dir build/issue-762/global-scan
```

The [probe log](chart-replay/probe-rerun.log) retains its old-row assertion;
the [global-scan log](chart-replay/global-rerun.log) identifies case 2,
attempt 222, difference `1.8579464722279226e-7` against `1e-10`.
The [bounded native restart](chart-replay/native-restart.json) tries that
same chart point as a guess for the full declaration, using the existing
ordered-start driver and its default controls. It returns no accepted
coordinates after 44 iterations, with reported residual
`0.6065943598651167`; the retained message includes the failed reference
preparation. This failed attempt is not evidence of nonexistence, nor the
reason for concluding that the chart declaration is absent. The source
inspection above establishes that distinction. The old rows were not
replaced with this attempt or with direct non-root evaluations.

Global existence at the two original target states remains undetermined.
The historical disposition completes the two replay sub-diagnostics under
#762. Corrected turning points and parameter derivatives do not establish a
root at either target, global nonexistence, or physical validation. This is
numerical verification of the reference change and an inspection of native
method availability. No physical uncertainty estimate is added; the evidence
has not received an independent delivered review or manuscript promotion.

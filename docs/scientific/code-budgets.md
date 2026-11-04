# Engine code budgets

A hard line ceiling per runtime layer, and one for the total, constrains
duplication and unnecessary machinery. The single numeric table is
`ARCHITECTURE.yaml:code_budgets`. Run `python3 tools/diffstat.py --budgets` to
display each count, ceiling and band. Both normal local check commands enforce
it. `check-affected.sh` provides the documentation-only route without a wheel
build; `check-all.sh` always builds and installs, before readiness checks.
Numerical evidence is a separate acceptance condition.

## Counting and scope

Count the lines of each runtime file as one pinned formatter lays it out,
without changing the file: C++ with `clang-format-18` and the style
`tools/line-count.clang-format`, Python with the `.venv` ruff (`uv.lock`
version) at 120 columns, where a trailing comma does not choose the layout. The
counter's `formatted_lines` holds the exact commands. Format original source
before excluding documentation or comments; comments are not reflowed.
For Python, mask only first-statement module, class, sync-function and
async-function string docstrings using their exact AST UTF-8 byte spans, then
count the union of lines occupied by retained tokens. Comments, whitespace,
layout and blank lines do not count. Other string expressions, assigned strings
(including blank-looking multiline data), annotations and declarations count;
code sharing a docstring line remains counted. Attribute prose belongs in its
class's NumPy Attributes section.
For C++, scan literals and preprocessing tokens after backslash-newline
continuation: exclude comments outside literals and whitespace-only lines.
Preserve ordinary/character/raw strings, raw-string data, header names and digit
separators. Declarations, braces, directives, macros and inactive branches count;
no preprocessing or compilation changes the count. Unterminated or unrecognised
lexical input fails visibly. This is formatted source size, not executed
instructions. The original-source formatter-disable refusal still runs first.
Count the
complete working `engine/native/` and `engine/src/epcsaft/` trees, except the
type stub `engine/src/epcsaft/_core.pyi`. The stub declares the compiled
extension's types and holds no Python logic; `tools/check-engine-architecture.py`
checks that its top-level names and class members match the installed
extension. Count tracked and untracked source alike, regardless of diff
attributes. Report native layers and Python separately, plus their total.
Reject unknown native layers, runtime file types and symlinks rather than
silently dropping them from the count. Tests, documentation, source data, build
outputs and third-party dependencies are outside this runtime budget. They
remain subject to scope/ownership review; none may hide runtime logic to evade
a ceiling.

The total ceiling may be tighter than the sum of the layer ceilings, never
looser: a layer ceiling bounds that layer, and the total independently bounds
the whole.

Formatted lines are a reproducible size proxy, not a measure of scientific
complexity. The formatter splits statements packed onto one line, and the
counter refuses a runtime file holding a formatter-disable directive
(`clang-format off`/`on`, `# fmt: off`/`on`/`skip`, `# yapf: disable`/`enable`),
naming its file and line. Packing still gains lines where the formatter keeps
the packed form: comma-joined C++ declarations such as `int a = 1, b = 2;`, and
a C++ braced list, which shares lines without a trailing comma and takes one
line per element with it. Removing useful comments, hiding code in generated
strings, relocating runtime code, or splitting/renaming files to
evade a budget is unacceptable. Spare lines never justify a second
thermodynamic or solver owner.

## Ceilings and the band

The band of a layer, or of the total, is 1% of its current count, rounded up.

1. **Raise to count plus band.** Above a hard ceiling, readiness fails.
   Refactor or reduce, or raise the ceiling with owner approval to exactly the
   count plus the band, recording the raise in the issue body in the same
   change. No override flag or automatic upward ratchet exists.
2. **Lower after consolidation.** When a ceiling sits more than one band above
   the count, readiness fails with the ceiling to lower to (count plus band;
   for the total, no more than the sum of the layer ceilings).
   Lower it in the same change; lowering needs no approval. A ceiling up to one
   band above the count needs no table change.
3. **Approve once per issue.** The owner approves a line allowance per layer
   when authorizing an issue, and the issue body records it. A merge whose
   raises stay within that allowance needs no second approval. Growth beyond
   it needs a new owner decision with the scientific requirement and
   alternatives.

A temporary over-budget working state may be simplified and repaired. Do not
freeze uncompilable code because an intermediate write exceeded a number. The
completed candidate must satisfy both band limits before acceptance.

Never satisfy a budget by dropping required physics, domain checks, derivative
terms, failure diagnostics or numerical falsifiers. If the smallest correct
method does not fit, report that evidence and revise the design/budget openly.
Budget changes must not silently change scientific scope.

## History

Past ceiling changes and their allowances are recorded in the owning issue
bodies and in `git log -p -- ARCHITECTURE.yaml`.

Issue #969 retains association conditioning observations across publication,
selected-endpoint refinement and reference central reuse. The approved allowance
is 180 gross / +60 net formatted runtime lines (Core 50/+10, Properties 65/+28,
Equilibrium 65/+22), with at most 220 raw additions. Against the refreshed
base `053b0ddff`, the implementation uses Core 19/+6, Properties 51/+24 and
Equilibrium 57/+21: 127 gross / +51 net. Properties requires its ceiling
to be raised to the measured count plus its rounded-up one-percent band;
`ARCHITECTURE.yaml` owns that ceiling. The other ceilings remain within band.

Issue #975 integration against `4c035248f` uses Properties +1, Bindings +25
and Python +15, total +41 net formatted runtime lines. Bindings increases from
2,589 to 2,618, the measured 2,592 plus band 26. The total ceiling remains
33,636, below the layer sum 33,667. The approved allowance is 185 gross /
+145 net; numerical qualification is separate.

## Duplication constraints

One common Helmholtz definition supplies properties and their derivatives.
Association elimination, density closure, reference transformations, equilibrium
row assembly and physical-input resolution each retain their accountable owner.
Python maps inputs/outputs; applications supply chemistry and specifications as
data. All equilibrium specifications use one driver; association, density,
fitting and discovery retain their distinct numerical operations using the
shared owners. Existing architecture checks cover known
include/constant/call patterns; size and regex checks cannot establish absence
of semantic duplication. Scientific/code review must inspect the actual calls.

This budget is a local readiness rule, not CI/CD and not a numerical accuracy or
runtime-speed threshold.

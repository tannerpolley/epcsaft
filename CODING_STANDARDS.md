# Coding standards

Reviewers apply these judgement calls to every change, and builders meet them
before asking for review. Mechanical rules are owned by their checks:
`tools/check-engine-architecture.py` (layers, trace payloads, policy entries),
`tools/check-analysis-layout.py` and `tools/check-analysis-declarations.py`
(analysis structure and Engine declarations), `tools/check-notebook-numbers.py`
(hand-typed notebook numbers) and `python3 tools/diffstat.py --budgets` (code
budgets). The [native calculation rule](AGENTS.md#native-calculation-rule) is
owned by `AGENTS.md`.

## Numerical verification

Use Numerical Verification-Driven Development:

- A persistent test falsifies a claim through a value, identity, derivative,
  limit or justified reference case, at a tolerance at least as tight as the
  claim.
- A performance or refactoring change shows equivalence within 1e-10 relative,
  with absolute floors for near-zero quantities (owner, 2026-10-06). Bit-for-bit
  comparison is reserved for unchanged default paths; a tighter tolerance needs
  a specific claim that requires it.
- Internal layouts, solver labels, call sequences and iteration counts are not
  correctness criteria; evaluation counts and timings are performance evidence
  for equivalent accepted work.
- A refused, unavailable or tolerance-gated result is a diagnosis to pursue,
  never an end state ([ALG-TERMINATION](engine/docs/science/algorithms.md)).
- A scientific change updates or explicitly revalidates its mathematics,
  algorithms, code and numerical evidence together; code behavior cannot
  redefine accepted doctrine.

## Evidence

- Homogeneous numerical rows are retained once in CSV/TSV; JSON stays compact
  and carries diagnostics, provenance, failures and claim limits.
- Retained evidence applies to its stated inputs, domain and inspected snapshot.
- Expensive installed-wheel replays belong in `analyses/`.

## Engine runtime shape

- `engine/` is the sole `epcsaft` runtime: one Python build, one Meson graph,
  one version, one Python-free native shared library and one extension. ASL alone
  uses a private archive from the same PIC objects. Follow `ARCHITECTURE.yaml`
  dependency direction.
- Compatibility imports, alternate package identities, hidden physical defaults
  and a public C++ ABI each need an approved scope change and a real external
  consumer.
- The repository has no nested repositories or submodules; `downstream/` clones
  stay untracked. Runtime and tests resolve inputs explicitly. Copied Data fixtures are bound to
  their packet lock, and each materialization retains the packet fingerprint and
  exact file hashes. Analysis reproductions exercise the installed wheel through
  public interfaces.

## Budgets

A change stays within its issue's line allowance. `ARCHITECTURE.yaml:code_budgets`
owns the numbers and [code budgets](docs/scientific/code-budgets.md) the counting,
the band and when a ceiling is raised or lowered.

## Documentation

A change updates the authoritative docs it affects in the same commit and
removes the instructions it supersedes.

# Codex environment

Codex runs `setup.sh` automatically for a newly created worktree. The script
derives the repository root from that worktree, creates its local `.venv`, and
builds and installs one non-editable wheel from the active Unified Engine
checkout with the repository's test, lint, numerical-test, and
wheel-build tools plus `mypy` for explicit static-type audits. The setup checks
that `epcsaft`, `epcsaft.equilibrium`, and
`epcsaft.regression` all resolve from that installed worktree wheel and removes
any additional ePC-SAFT distributions from a reused environment.

The setup does not consume a sibling source tree or mutable wheel from another
checkout. The Test action rebuilds and reinstalls the worktree wheel before
pytest, so edits cannot be tested against stale installed code. Wheel output is
isolated under the ignored build tree and carries no release authority.

The Check action runs Ruff's core syntax, import, and undefined-name rules plus
strict mypy. The Check native action configures an isolated Clang 18 compilation
database and runs the repository's high-signal clang-tidy policy over only the
Engine-owned C++ sources. Check native sanitizers builds a separate wheel with
AddressSanitizer and UndefinedBehaviorSanitizer and runs a focused EOS,
equilibrium, and end-to-end scientific smoke suite. Native checks require the
system `clang++-18`, `clang-tidy-18`, and `run-clang-tidy-18` commands; on Debian
or Ubuntu they are provided by the `clang-18`, `clang-tidy-18`, and
`clang-tools-18` packages.

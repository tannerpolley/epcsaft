# Codex environment

Codex runs `setup.sh` automatically for a newly created worktree. The script
derives the repository root from that worktree, creates its local `.venv`, and
builds and installs one non-editable wheel from the active Unified Engine
checkout with the repository's test, lint, numerical-test, and
wheel-build tools. The setup checks that `epcsaft`, `epcsaft.equilibrium`, and
`epcsaft.regression` all resolve from that installed worktree wheel and removes
any additional ePC-SAFT distributions from a reused environment.

The setup does not consume a sibling source tree or mutable wheel from another
checkout. The Test action rebuilds and reinstalls the worktree wheel before
pytest, so edits cannot be tested against stale installed code. Wheel output is
isolated under the ignored build tree and carries no release authority.

The Check action runs Ruff's core syntax, import, and undefined-name rules. The
broader style and typing backlogs are not misrepresented as passing worktree
gates; they require separate, bounded cleanup campaigns.

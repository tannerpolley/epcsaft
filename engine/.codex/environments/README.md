# Codex environment

Codex runs `setup.sh` automatically for a newly created worktree. The script
derives the repository root from that worktree, creates its local `.venv`, and
installs this Unified Engine checkout in editable mode with the repository's
test, lint, type-check, numerical-test, and wheel-build tools. The setup checks
that `epcsaft`, `epcsaft.equilibrium`, and `epcsaft.regression` all resolve from
the active worktree and removes any additional ePC-SAFT distributions from a
reused environment.

The setup does not consume a sibling source tree or mutable wheel from another
checkout. Actions in `environment.toml` always run against the active
worktree's `.venv`; wheel output is isolated under the ignored build tree.

# Codex environment

Codex runs `setup.sh` automatically for a newly created worktree. The script
derives the repository root from that worktree, creates its local `.venv`, and
installs this EOS checkout in editable mode with the repository's test, lint,
type-check, and wheel-build tools.

The setup does not consume a sibling source tree or mutable wheel from another
checkout. Actions in `environment.toml` always run against the active
worktree's `.venv`.

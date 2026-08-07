# Root Codex environment

This is the sole project environment. It builds the current `engine/` tree as
one non-editable wheel into the root `.venv` and supplies every directory's
readiness and scientific-evidence dependencies. It is safe in linked worktrees
because paths are derived from this script location.

Each checkout installs from its own verified wheel path. Identical worktrees
reuse content-addressed wheel bytes from the machine-local cache, while CMake
keeps incremental Ninja state checkout-local and shares compiler results
through `ccache`.

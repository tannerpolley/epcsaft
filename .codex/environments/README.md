# Root Codex environment

[`../../README.md`](../../README.md) owns the setup and validation commands.
This is the sole project environment. It performs an exact synchronization of
all workspace package dependencies selected by the committed root `uv.lock`,
without installing the workspace packages themselves. `uv.lock` is the only
dependency lock; dependency declarations and lock changes land together. It
then builds the checked-out `engine/` tree as one non-editable wheel into the
root `.venv` and supplies every directory's readiness and scientific-evidence
dependencies. A stale lock fails with `--locked` instead of being repaired.

The setup state follows the checked-out branch: paths, Engine metadata, source,
native dependencies, and build tools determine whether the checkout-local
wheel is current. Identical worktrees can reuse matching wheel bytes from the
machine-local content-addressed cache, while Meson keeps incremental Ninja
state checkout-local. If `ccache` is installed, the native compiler results are
also reused; it is optional.

Each checkout installs from its own verified wheel path, with dependency
resolution disabled because the exact root lock already supplied the Engine
runtime dependencies.

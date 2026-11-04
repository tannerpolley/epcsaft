# Root Codex environment

[`../../engine/README.md`](../../engine/README.md) "Build and installation" and "Checks" own the setup and validation commands.
This is the sole project environment. It performs an inexact synchronization (extra packages are kept) of
all workspace package dependencies selected by the committed root `uv.lock`,
without installing the workspace packages themselves. `uv.lock` is the only
dependency lock; dependency declarations and lock changes land together. It
then installs the checked-out `engine/` tree editable into the root `.venv`
through `tools/build-engine.sh` and supplies every directory's readiness and
scientific-evidence dependencies. A stale lock fails with `--locked` instead of
being repaired.

The editable install builds into the checkout's `build/engine-meson/shared` Meson
tree, which `import epcsaft`, native tests, `tools/native-run` and wheel builds
share; Ninja rebuilds only what changed, one Ninja at a time (`tools/locked-ninja`).
Re-running setup on an editable install only compiles what changed. If `ccache` is installed,
compiler results are reused across checkouts; it is optional.

Engine, native-test and `tools/native-run` compilation default to one build
job per 1.5 GiB of available memory (free RAM or, if tighter, the nearest
cgroup headroom), at least one and at most the CPU count
(`tools/build-engine.sh`). `EPCSAFT_BUILD_JOBS` at setup overrides it with any
positive integer, without a CPU cap; for example,
`EPCSAFT_BUILD_JOBS=1 bash tools/check-all.sh` uses one job. Setup records the
resolved count in `build/environment-wheel/build-jobs`, which `tools/native-run`
reads.
The setting changes build scheduling without changing compiler optimization.
Import-time rebuilds keep the job count fixed when setup installed the editable build;
`EPCSAFT_BUILD_JOBS` does not change them.

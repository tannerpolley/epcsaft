#!/usr/bin/env bash
# run_autotools_clang_conda_build sources this; the subshell keeps the
# builder's strict shell options and exit trap out of the wrapper.
set -eo pipefail
(source packaging/conda-forge/build-ipopt-win.sh)

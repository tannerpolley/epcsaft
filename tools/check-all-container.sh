#!/usr/bin/env bash
# Run tools/check-all.sh on this checkout's HEAD commit from a cold build in a digest-pinned ubuntu:24.04
# container, the hosted runner's system, at no GitHub minutes (#576). Noble's apt Ipopt 3.11.9, Ceres 2.2
# and MUMPS 5.6 are this host's; uv is the host binary; CppAD and HiGHS come from their pinned wraps.
# The container ID is in build/check-all-container.cid. EPCSAFT_BUILD_JOBS passes through.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
cid=build/check-all-container.cid
mkdir -p build && rm -f "$cid"
docker run --rm --cidfile "$cid" -e EPCSAFT_BUILD_JOBS -v "$(git rev-parse --path-format=absolute --git-common-dir):/repo.git:ro" \
    -v "$(command -v uv):/usr/local/bin/uv:ro" \
    ubuntu:24.04@sha256:534baea6a22c03a63003dbc8dbe78fe34bc0d7e595d9a9dc9834884ff530eb55 bash -euo pipefail -c '
apt-get update -qq && DEBIAN_FRONTEND=noninteractive apt-get install -qq --no-install-recommends build-essential \
    ca-certificates clang-format-18 cmake coinor-libipopt-dev file git libceres-dev libeigen3-dev ninja-build \
    pkg-config time >/dev/null
git config --global --add safe.directory "*"
git clone -q file:///repo.git /src && cd /src && git checkout -q "$1"
nice -n 10 bash tools/check-all.sh' _ "$(git rev-parse HEAD)"

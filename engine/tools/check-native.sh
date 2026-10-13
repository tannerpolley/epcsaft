#!/usr/bin/env bash
set -euo pipefail

engine_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workspace_root="$(dirname "$engine_root")"
cd "$engine_root"
export PATH="$workspace_root/.venv/bin:$PATH"

# --changed [base] checks only the native translation units changed since base (default: the merge base
# with origin/main) or including a changed file, by the lint build's own dependency records. Every unit is
# checked with no argument, after a change to any .clang-tidy, a Meson file, the subprojects, the Python
# environment files or this script, and while no full sweep has passed with the current toolchain.
all_units='engine/(native|tests/native)/.*\.cpp$'
selective=false
if [[ "${1:-}" == --changed ]]; then
    base="${2:-$(git merge-base origin/main HEAD)}" || { echo "--changed: no merge base with origin/main" >&2; exit 2; }
    git rev-parse --verify --quiet "$base^{commit}" >/dev/null || { echo "--changed: unknown base $base" >&2; exit 2; }
    # One checked command substitution per listing (not mapfile < <(...), which hides a failure): a failed
    # listing must stop the gate, never read as "nothing changed".
    changed="$(git diff --name-only "$base" -- :/)" || { echo "--changed: listing changes since $base failed" >&2; exit 1; }
    untracked="$(git ls-files --full-name --others --exclude-standard -- :/)" \
        || { echo "--changed: listing untracked files failed" >&2; exit 1; }
    changed="$changed"$'\n'"$untracked"
    if grep -Eq '(^|/)\.clang-tidy$|^engine/((.*/)?meson(\.build|\.options|_options\.txt)|build/subprojects/.*|tools/check-native\.sh)$|^((engine/)?pyproject\.toml|uv\.lock|\.codex/environments/.*)$' <<<"$changed"; then
        echo "lint configuration, build or environment files changed since $base: full sweep"
    else
        selective=true
    fi
fi

for required_command in meson cmake ninja clang-18 clang++-18 clang-tidy-18 run-clang-tidy-18; do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        echo "missing required command: $required_command" >&2
        exit 1
    fi
done
if [[ ! -x "$workspace_root/.venv/bin/python" ]]; then
    echo "missing root environment; run bash .codex/environments/setup.sh" >&2
    exit 1
fi

build_root="$workspace_root/build/native-quality"
# The checkout's build-job count (or EPCSAFT_BUILD_JOBS) sets the compile and clang-tidy parallelism.
jobs="${EPCSAFT_BUILD_JOBS:-$(cat "$workspace_root/build/environment-wheel/build-jobs" 2>/dev/null || echo 2)}"
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || { echo "build-job count must be a positive integer (got $jobs)" >&2; exit 2; }
extension_suffix="$($workspace_root/.venv/bin/python -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')"
# The tree is always the Clang debug build, so a warm tree reconfigures instead of rebuilding.
setup_args=()
if [[ -f "$build_root/meson-private/coredata.dat" ]]; then
    setup_args+=(--reconfigure)
elif [[ -d "$build_root" ]]; then
    rm -rf -- "$build_root"
fi
# Clang ships no quadmath.h; the binary128 test units take GCC's after Clang's own headers.
setup_args+=("-Dcpp_args=-idirafter$(gcc -print-file-name=include)")
# ccache as tools/locked-ninja sets it up; debug info names the build directory, so a new worktree
# reuses another's objects only when ccache leaves the directory out of its hash.
export CCACHE_BASEDIR="$workspace_root" CCACHE_NOHASHDIR=1
CC="ccache clang-18" CXX="ccache clang++-18" meson setup \
    "${setup_args[@]}" \
    --buildtype=debug \
    "$build_root" \
    "$engine_root"
meson compile -C "$build_root" -j "$jobs" \
    "_core$extension_suffix"
# Building every native unit's object leaves Ninja holding each unit's current dependency record.
object_list="$("$workspace_root/.venv/bin/python" -c '
import json, re, sys
for entry in json.load(open(sys.argv[1])):
    if re.search(sys.argv[2], entry["file"]): print(entry["output"])' "$build_root/compile_commands.json" "$all_units")" \
    || { echo "listing native units from compile_commands.json failed" >&2; exit 1; }
[[ -n "$object_list" ]] || { echo "no native unit in $build_root/compile_commands.json" >&2; exit 1; }
mapfile -t objects <<<"$object_list"
ninja -C "$build_root" -j "$jobs" "${objects[@]}"
# Each record lists the object's source first, then every file it includes, relative to the build directory.
deps="$(ninja -C "$build_root" -t deps "${objects[@]}")" || { echo "reading Ninja dependency records failed" >&2; exit 1; }
current="$(grep -c ' (VALID)$' <<<"$deps" || true)"
((current == ${#objects[@]})) || { echo "$current current dependency records for ${#objects[@]} native units" >&2; exit 1; }
# The toolchain is the Clang binaries and the content of every header outside the repository that a unit
# includes (system, Python and pybind11 headers); a marker records that a full sweep passed with it.
toolchain="$(cd "$build_root" && {
    readlink -f "$(command -v clang-tidy-18)" "$(command -v clang++-18)" | xargs -d '\n' sha256sum &&
    awk '/^    / && $1 ~ /^(\/|\.\.\/)/ && $1 !~ /^\.\.\/\.\.\/engine\// { print $1 }' <<<"$deps" | sort -u | xargs -d '\n' sha256sum
} | cut -d ' ' -f 1 | sha256sum | cut -c 1-16)" || { echo "hashing the toolchain failed" >&2; exit 1; }
swept="${XDG_CACHE_HOME:-$HOME/.cache}/epcsaft/clang-tidy-swept/$toolchain"
if [[ "$selective" == true && ! -f "$swept" ]]; then
    echo "no full sweep has passed with toolchain $toolchain: full sweep"
    selective=false
fi
files=("$all_units")
if [[ "$selective" == true ]]; then
    selected="$(CHANGED="$changed" awk '
        BEGIN { n = split(ENVIRON["CHANGED"], c, "\n"); for (i = 1; i <= n; i++) changed[c[i]] = 1 }
        /^[^ ]/ { source = "" }
        /^    / { path = $1; sub(/^(\.\.\/)+/, "", path); if (source == "") source = path; if (path in changed) units[source] = 1 }
        END { for (unit in units) print unit }' <<<"$deps" | sort | sed 's/\./\\./g; s/^/\//; s/$/$/')" \
        || { echo "selecting units from the dependency records failed" >&2; exit 1; }
    [[ -n "$selected" ]] || { echo "no native unit is or includes a file changed since $base"; exit 0; }
    mapfile -t files <<<"$selected"
    echo "linting ${#files[@]} native units changed or including a changed file"
fi

run-clang-tidy-18 \
    -clang-tidy-binary clang-tidy-18 \
    -p "$build_root" \
    -j "$jobs" \
    -quiet \
    "${files[@]}"
[[ "$selective" == true ]] || { mkdir -p "$(dirname "$swept")"; touch "$swept"; }

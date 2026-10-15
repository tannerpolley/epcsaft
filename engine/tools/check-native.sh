#!/usr/bin/env bash
set -euo pipefail

[[ "${1:-}" == --prepared && -d "${2:-}" ]] || {
    echo 'use tools/check-pr-native --pr N --head <full-sha> --base <full-sha> to prepare the shared tree' >&2; exit 2;
}
run_root="$(realpath "$2")"
shift 2
results="${EPCSAFT_GATE_RESULTS:-$run_root}"
workspace_root="$run_root/source"
engine_root="$workspace_root/engine"
build_root="$run_root/meson"
cd "$engine_root"
export PATH="$run_root/venv/bin:$PATH"

# --changed [base] checks only the native translation units changed since base (default: the merge base
# with origin/main) or including a changed file, by the lint build's own dependency records. Every unit is
# checked with no argument, after a change to any .clang-tidy, a Meson file, the subprojects, the Python
# environment files, and while no full sweep has passed with the current toolchain.
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
    if grep -Eq '(^|/)\.clang-tidy$|^engine/((.*/)?meson(\.build|\.options|_options\.txt)|build/subprojects/.*)$|^((engine/)?pyproject\.toml|uv\.lock|\.codex/environments/.*)$' <<<"$changed"; then
        echo "lint configuration, build or environment files changed since $base: full sweep"
    else
        selective=true
    fi
fi

jobs="${EPCSAFT_BUILD_JOBS:-$(nproc)}"
[[ "$jobs" =~ ^[1-9][0-9]*$ ]] || { echo "build-job count must be positive" >&2; exit 2; }
# Preparation compiled every native/test object before either check entered this tree.
object_list="$(python -c '
import json, re, sys
for entry in json.load(open(sys.argv[1])):
    if re.search(sys.argv[2], entry["file"]): print(entry["output"])' "$build_root/compile_commands.json" "$all_units")"
[[ -n "$object_list" ]] || { echo "no native unit in compile_commands.json" >&2; exit 1; }
mapfile -t objects <<<"$object_list"
# Each record identifies an object and lists its dependencies relative to the build directory.
deps="$(ninja -C "$build_root" -t deps "${objects[@]}")" || { echo "reading Ninja dependency records failed" >&2; exit 1; }
current="$(grep -c ' (VALID)$' <<<"$deps" || true)"
((current == ${#objects[@]})) || { echo "$current current dependency records for ${#objects[@]} native units" >&2; exit 1; }
# Sweep identity excludes tracked inputs and the generated HiGHS header embedding the Git revision.
# The complete compile-input identity still owns prepared-pass invalidation.
toolchain="$(python - "$run_root" <<'IDENTITY'
import hashlib, json, pathlib, runpy, subprocess, sys
root = pathlib.Path(sys.argv[1])
source = root / 'source'
tracked = {(source / path).resolve() for path in subprocess.check_output(
    ['git', 'ls-files', '-z'], cwd=source, text=True).split('\0') if path}
revision_header = root / 'meson/build/subprojects/HiGHS-1.15.1/highs/HConfig.h'
external = [pathlib.Path(path) for path in json.loads((root / 'external-headers.json').read_text())]
digest = runpy.run_path(str(source / 'tools/check-pr-native'))['digest']
identity = json.loads((root / 'prepared.json').read_text())
identity.pop('head')
identity['external'] = digest(path for path in external if path not in tracked and path != revision_header)
print(hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest())
IDENTITY
)"
swept="${XDG_CACHE_HOME:-$HOME/.cache}/epcsaft/clang-tidy-swept/$toolchain"
if [[ "$selective" == true && ! -f "$swept" ]]; then
    echo "no full sweep has passed with toolchain $toolchain: full sweep"
    selective=false
fi
files=("$all_units")
printf '%s\n' "${objects[@]}" >"$results/lint-objects.txt"
if [[ "$selective" == true ]]; then
    selected="$(CHANGED="$changed" python - "$build_root" "$workspace_root" "$all_units" <<'SELECT'
import json, os, pathlib, re, subprocess, sys
build, source = map(pathlib.Path, sys.argv[1:3])
units = {entry['output']: (pathlib.Path(entry['directory']) / entry['file']).resolve()
         for entry in json.loads((build / 'compile_commands.json').read_text())
         if re.search(sys.argv[3], entry['file'])}
changed = set(os.environ['CHANGED'].splitlines())
selected, unit = set(), None
for line in subprocess.check_output(['ninja', '-C', str(build), '-t', 'deps'], text=True).splitlines():
    if not line.startswith('    '):
        unit = units.get(line.split(': #deps ', 1)[0])
    elif unit is not None:
        path = (build / line.strip()).resolve()
        if path.is_relative_to(source) and path.relative_to(source).as_posix() in changed:
            selected.add(str(unit))
print('\n'.join('^' + re.escape(path) + '$' for path in sorted(selected)))
SELECT
)"
    : >"$results/lint-selection.txt"
    [[ -n "$selected" ]] || { echo "no native unit is or includes a file changed since $base"; exit 0; }
    mapfile -t files <<<"$selected"
    echo "linting ${#files[@]} native units changed or including a changed file"
fi

printf '%s\n' "${files[@]}" >"$results/lint-selection.txt"
python - "$build_root/compile_commands.json" "${files[@]}" <<'MATCH'
import json, os, re, sys
selection = re.compile('|'.join(sys.argv[2:]))
if not any(selection.search(os.path.normpath(os.path.join(entry['directory'], entry['file'])))
           for entry in json.load(open(sys.argv[1]))):
    sys.exit('lint selection matches no compile entry in compile_commands.json')
MATCH
run-clang-tidy-18 \
    -clang-tidy-binary clang-tidy-18 \
    -p "$build_root" \
    -j "$jobs" \
    -quiet \
    "${files[@]}"
[[ "$selective" == true ]] || { mkdir -p "$(dirname "$swept")"; touch "$swept"; }

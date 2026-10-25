#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
stamp="$repo_root/build/environment-wheel/meson-build-dir"
link="$repo_root/compile_commands.json"

if [[ ! -r "$stamp" ]]; then
    echo "Serena C++ setup requires a prepared Engine build: missing $stamp" >&2
    exit 1
fi

build_dir="$(<"$stamp")"
if [[ "$build_dir" != /* ]]; then
    build_dir="$repo_root/$build_dir"
fi
build_dir="$(readlink -f -- "$build_dir")"
compile_db="$build_dir/compile_commands.json"
if [[ ! -f "$compile_db" ]]; then
    echo "Serena C++ setup requires compile_commands.json in the prepared Engine build: $compile_db" >&2
    exit 1
fi
if [[ -e "$link" && ! -L "$link" ]]; then
    echo "Serena C++ setup will not replace an existing non-symlink: $link" >&2
    exit 1
fi
if [[ -L "$link" && "$(readlink -f -- "$link")" == "$compile_db" ]]; then
    echo "Serena C++ compile database already current: $compile_db"
    exit 0
fi

tmp_link="$(mktemp "$repo_root/.compile_commands.json.serena.XXXXXX")"
trap 'rm -f -- "$tmp_link"' EXIT
rm -f -- "$tmp_link"
ln -s -- "$compile_db" "$tmp_link"
mv -Tf -- "$tmp_link" "$link"
echo "Serena C++ compile database: $compile_db"

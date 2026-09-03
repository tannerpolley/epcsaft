#!/usr/bin/env bash
# Additions and deletions per category, for any commit range.
#
#   tools/git-categories.sh                  # since the branch point
#   tools/git-categories.sh HEAD~5..HEAD     # a range
#   tools/git-categories.sh --code           # only the categories that are code
#
# Categories are the ones in engine/docs/greenfield/CATEGORIES.md. The point is that
# "green additions" should mean code someone wrote, not data or prose.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
only_code=0
range=""
for a in "$@"; do case "$a" in --code) only_code=1 ;; *) range="$a" ;; esac; done
[ -n "$range" ] || range="$(git merge-base main HEAD)..HEAD"

git diff --numstat "$range" | awk -v only_code="$only_code" '
function cat(p) {
  # the superseded tree, deleted by the rebuild: never code we wrote
  if (p ~ /^engine\/cpp\//)                                return "removed"
  if (p ~ /^engine\/src\/epcsaft\/_native\//)              return "removed"
  if (p ~ /^engine\/(tests|src)\/.*\/(numerical|behavior)\//) return "removed"
  if (p ~ /^engine\/native\/core\//)                      return "core"
  if (p ~ /^engine\/native\/(surface|bindings)\//)        return "surface"
  if (p ~ /^engine\/src\//)                               return "surface"
  if (p ~ /^engine\/native\//)                            return "engine"
  if (p ~ /^(tools|data\/tools)\//)                       return "canon:enforcer"
  if (p ~ /^engine\/scripts\//)                           return "canon:enforcer"
  if (p ~ /^engine\/tests\/.*\.json$/)                    return "canon:data"
  if (p ~ /^engine\/tests\//)                             return "canon:enforcer"
  if (p ~ /briefs\//)                                     return "canon:brief"
  if (p ~ /\.md$|\.tex$|\.bib$/)                          return "canon:record"
  if (p ~ /^data\//)                                      return "canon:data"
  if (p ~ /^analyses\//)                                  return "analyses"
  return "build"
}
{ c = cat($3); add[c] += $1; del[c] += $2; files[c]++ }
END {
  order[1]="core"; order[2]="surface"; order[3]="engine"; order[4]="build"
  order[5]="canon:enforcer"; order[6]="canon:brief"; order[7]="canon:record"
  order[8]="canon:data"; order[9]="analyses"; order[10]="removed"
  printf "  %-18s %9s %9s %7s\n", "category", "added", "deleted", "files"
  code_add = 0
  for (i = 1; i <= 10; i++) {
    c = order[i]
    is_code = (c=="core" || c=="surface" || c=="engine" || c=="canon:enforcer" || c=="build")
    if (c=="removed") is_code = 0
    if (only_code && !is_code) continue
    if (!(c in add)) continue
    if (is_code) code_add += add[c]
    printf "  %-18s %+9d %9d %7d\n", c, add[c], -del[c], files[c]
  }
  printf "  %-18s %+9d\n", "CODE TOTAL", code_add
}'

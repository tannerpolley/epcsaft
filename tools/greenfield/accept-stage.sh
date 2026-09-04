#!/bin/bash
# Coordinator's acceptance for a stage brief: wipe and rebuild, every gate in brief order, then
# every mutation applied, rebuilt, rerun, expected red, restored. Usage:
#   tools/greenfield/accept-stage.sh <worktree> <brief.json relative to worktree> <log> [first-mutations...]
# A mutant that does not build, or stays green, is reported INVALID; the coordinator sorts those
# into tooling artifact, stale sed, and real finding before reporting (STAGE-LOG.md, Stage 5).
set -u
W="$1"; B="$2"; LOG="$3"; shift 3; FIRST="$*"
cd "$W"; export PATH="$W/.venv/bin:$PATH"; : > "$LOG"
say(){ echo "[$(date +%H:%M:%S)] $*" >> "$LOG"; }
say "clean rebuild: meson setup --wipe"
if ! timeout 900 meson setup --wipe engine/build/greenfield engine >/tmp/accept-setup.log 2>&1; then say "SETUP FAILED"; tail -20 /tmp/accept-setup.log >> "$LOG"; exit 1; fi
if ! timeout 1200 meson compile -C engine/build/greenfield >/tmp/accept-compile.log 2>&1; then say "CLEAN COMPILE FAILED"; tail -30 /tmp/accept-compile.log >> "$LOG"; exit 1; fi
say "clean build ok"
python3 - "$B" > /tmp/accept-gates.tsv <<'PY'
import json,sys;b=json.load(open(sys.argv[1]))
for g in b['gates']:
    m=g['mutation']; print("\t".join([g['id'],g['needs'],g['command'],m['file'],m['sed'],m['expect']]))
PY
fails=0
while IFS=$'\t' read -r id needs cmd mfile sed expect; do
  if timeout 900 bash -c "$cmd" >/tmp/accept-gate-$id.log 2>&1; then say "GATE $id base=green"; else say "GATE $id base=RED"; fails=$((fails+1)); tail -8 /tmp/accept-gate-$id.log | sed 's/^/    /' >> "$LOG"; fi
done < /tmp/accept-gates.tsv
say "base gates done, red=$fails"
order=$(python3 -c "
import csv,sys;rows=list(csv.reader(open('/tmp/accept-gates.tsv'),delimiter='\t'))
first='$FIRST'.split(); ids=[r[0] for r in rows]; print(' '.join([i for i in first if i in ids]+[i for i in ids if i not in first]))")
for id in $order; do
  IFS=$'\t' read -r _ needs cmd mfile sed expect < <(grep -P "^$id\t" /tmp/accept-gates.tsv)
  if [ ! -f "$mfile" ]; then say "MUT $id INVALID (file missing: $mfile)"; continue; fi
  before=$(md5sum "$mfile" | cut -c1-8); sed -i "$sed" "$mfile"; after=$(md5sum "$mfile" | cut -c1-8)
  if [ "$before" = "$after" ]; then say "MUT $id INVALID (sed matched nothing)"; git checkout -q -- "$mfile"; continue; fi
  built=1
  case "$mfile" in *.cpp|*.hpp)
    if [ "$needs" != "none" ] || [ "$id" = "build" ]; then
      if ! timeout 1200 meson compile -C engine/build/greenfield >/tmp/accept-mut-$id-build.log 2>&1; then built=0; fi
    fi;;
  esac
  if [ $built = 0 ]; then
    if [ "$id" = "build" ]; then say "MUT $id RED (does not compile, as the gate requires)"; else say "MUT $id INVALID (mutant does not build)"; tail -5 /tmp/accept-mut-$id-build.log | sed 's/^/    /' >> "$LOG"; fi
  else
    if timeout 900 bash -c "$cmd" >/tmp/accept-mut-$id.log 2>&1; then say "MUT $id STAYS GREEN -> INVALID"; else say "MUT $id RED"; fi
  fi
  git checkout -q -- "$mfile"
  case "$mfile" in *.cpp|*.hpp) timeout 1200 meson compile -C engine/build/greenfield >/dev/null 2>&1 || say "WARN restore compile failed after $id";; esac
done
say "ACCEPTANCE RUN COMPLETE"

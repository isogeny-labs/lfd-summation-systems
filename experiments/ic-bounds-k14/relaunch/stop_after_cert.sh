#!/bin/bash
# Stop the GPU pod once the checker record of the control instance has been copied (or at DEADLINE, epoch seconds):
# the second instance cannot finish before the pod's lifetime cap. Final copy first, then terminate.
set -u
NAME=icb-k14b; DEADLINE=$1
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../../.." && pwd); DEST="$HERE/fetched"
T="python3 $ROOT/scripts/ic_bounds_big_pod.py"; LOG="$DEST/stop-$NAME.log"
say() { echo "$(date -u +%FT%TZ) $*" >> "$LOG"; }
say "stopper start, deadline $(date -u -r "$DEADLINE" +%FT%TZ)"
while true; do
  $T get $NAME out/certcheck-proper-p114689-k14-D6.json "$DEST/certcheck-proper-p114689-k14-D6.json" > /dev/null 2>&1
  [ -s "$DEST/certcheck-proper-p114689-k14-D6.json" ] && { say "checker record copied"; break; }
  [ "$(date +%s)" -ge "$DEADLINE" ] && { say "deadline reached without the checker record"; break; }
  sleep 60
done
pkill -f "mirror.sh $NAME" && say "mirror stopped"
for f in out/chain-status.json out/proper-p114689-k14-D6.json out/lift-proper-p114689-k14-D6to7.json out/certcheck-proper-p114689-k14-D6.status.json \
         out/unit-p163841-k14-D6.status.json logs/chain-proper-p114689-k14-D6.log logs/chain-lift-proper-p114689-k14-try1.log \
         logs/chain-certcheck-proper-p114689-k14.log logs/chain-selftest-engine.log logs/chain-selftest-lift.log; do
  $T get $NAME "$f" "$DEST/$(basename "$f")" --max-bytes 60000000 > /dev/null 2>&1 || say "not copied: $f"
done
$T jobs $NAME > "$DEST/jobs-$NAME.txt" 2>/dev/null
id=$(python3 -c "import json,os;print(json.load(open(os.path.join(os.environ['POD_STATE_DIR'],'$NAME.json')))['id'])")
say "terminating $id: $($T terminate "$id" --yes 2>&1 | tr '\n' ' ')"
say "pods still running: $($T ls 2>&1 | grep -v EXITED | grep -v gone | tr '\n' ';')"

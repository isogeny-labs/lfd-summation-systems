#!/bin/bash
# Copy a pod's records off it every EVERY seconds, from a process outside any agent session, and terminate the
# pod when its work is done.  usage: mirror.sh NAME JOB DEST DONE EVERY FILE...
#   NAME  pod name of scripts/ic_bounds_big_pod.py (token and id under $POD_STATE_DIR)
#   JOB   the pod job to watch;  DEST local directory;  FILE... paths under /workspace to copy (missing ones are skipped)
#   DONE  "status:<file>" = done when DEST/<file> says "finished": true;  "file:<file>" = done when DEST/<file> exists
# Also done when JOB has left the state "running".  Then: one more copy, terminate, exit.
set -u
NAME=$1; JOB=$2; DEST=$3; DONE=$4; EVERY=$5; shift 5
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../../.." && pwd)
T="python3 $ROOT/scripts/ic_bounds_big_pod.py"
mkdir -p "$DEST"; LOG="$DEST/mirror-$NAME.log"
say() { echo "$(date -u +%FT%TZ) $*" >> "$LOG"; }
copy() { for f in "$@"; do $T get "$NAME" "$f" "$DEST/$(basename "$f")" > /dev/null 2>&1; done
         $T tail "$NAME" "$JOB" --bytes 6000 > "$DEST/$JOB.tail.tmp" 2>/dev/null && mv "$DEST/$JOB.tail.tmp" "$DEST/$JOB.tail"; }
say "mirror start name=$NAME job=$JOB every=${EVERY}s pid=$$"
miss=0; seen=0
while true; do
  jobs=$($T jobs "$NAME" 2>/dev/null)
  if [ -n "$jobs" ]; then miss=0; seen=1; echo "$jobs" > "$DEST/jobs-$NAME.txt"; else miss=$((miss+1)); fi
  copy "$@"
  state=$(echo "$jobs" | awk -v j="$JOB" '$1==j{print $2}')
  done=0
  case "$DONE" in
    status:*) grep -q '"finished": true' "$DEST/${DONE#status:}" 2>/dev/null && done=1;;
    file:*)   [ -s "$DEST/${DONE#file:}" ] && done=1;;
  esac
  [ -n "$state" ] && [ "$state" != running ] && done=1
  say "job=$JOB state=${state:-?} done=$done miss=$miss $(ls "$DEST" | wc -l | tr -d ' ') files"
  if [ "$done" = 1 ]; then
    sleep 20; copy "$@"
    id=$(python3 -c "import json,os;print(json.load(open(os.path.join(os.environ['POD_STATE_DIR'],'$NAME.json')))['id'])")
    say "done; terminating pod $id: $($T terminate "$id" --yes 2>&1 | tr '\n' ' ')"
    exit 0
  fi
  # unreachable for an hour after having answered once (or for three hours from the start): the pod is gone
  if { [ "$seen" = 1 ] && [ $((miss*EVERY)) -ge 3600 ]; } || [ $((miss*EVERY)) -ge 10800 ]; then say "pod unreachable; exiting"; exit 1; fi
  sleep "$EVERY"
done

#!/bin/bash
# Launch the binary deep runs on ONE H200 pod (experiments/ic-bounds-binary-deep/PLAN.md).  Run by the main session only.
# Needs RUNPOD_API_KEY and POD_STATE_DIR (a directory that survives the session; tokens are written there).
# usage: launch.sh MAX_SECONDS [DISK_GB] [chain args...]
#   MAX_SECONDS  lifetime cap of the pod (Runpod stops it then); the chain gets --deadline MAX_SECONDS-1800
#   DISK_GB      container disk, default 60 (records and logs only: the chain writes no checkpoints; with checkpoints
#                set it to 20 + 1.1 x basis_worst_gb of the largest closure, see plan-numbers.json)
# Lessons of ic-bounds-k14/RELAUNCH.md built in: the copy loop (mirror.sh) is started FIRST, detached, outside the
# session, and copies the logs too; the device-only pool wrapper is uploaded as nopool_run.py; one checkpoint at most;
# the disk is sized explicitly; the lifetime cap is an argument; data centres are polled once a minute.
set -eu
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../../.." && pwd); cd "$ROOT"
MAXS=${1:?usage: launch.sh MAX_SECONDS [DISK_GB] [chain args]}; DISK=${2:-60}; shift $(( $# >= 2 ? 2 : 1 )) || true
T="python3 scripts/ic_bounds_big_pod.py"; DEST="$ROOT/experiments/ic-bounds-binary-deep/fetched"
EXP=experiments/ic-bounds-binary-deep
: "${POD_STATE_DIR:?set POD_STATE_DIR}"; mkdir -p "$POD_STATE_DIR" "$DEST"; chmod 700 "$POD_STATE_DIR"
N=${POD_NAME:-icb-bin-deep}; CUR=""
cleanup() { [ -n "$CUR" ] && { id=$(python3 -c "import json;print(json.load(open('$POD_STATE_DIR/$CUR.json'))['id'])" 2>/dev/null) && [ -n "$id" ] && $T terminate "$id" --yes; echo "launch of $CUR failed; pod terminated"; }; }
trap cleanup ERR   # a pod whose job did not start must not idle until its lifetime cap
wait_ping() { for i in $(seq 1 60); do $T ping "$1" > /dev/null 2>&1 && return 0; sleep 20; done; echo "pod $1 did not answer in 20 min"; return 1; }
# the copy list: every record, status file and log the chain can write (missing ones are skipped by mirror.sh)
FILES="out/chain-status.json out/selftest-gpu.json out/results.jsonl logs/chain-selftest.log"
for s in $(python3 - <<E
import json
c=json.load(open('$ROOT/$EXP/cells.json'))
for cell in c['order']:
    for s in cell['systems']:
        n=s['name']
        print('out/%s.json logs/chain-%s.log' % (n,n), ' '.join('out/%s.json.D%d.json' % (n,d) for d in range(0,s['max_degree']+1)))
E
); do FILES="$FILES $s"; done
CREATE="--gpu \"NVIDIA H200\" --min-ram-per-gpu ${MIN_RAM:-64} --min-vcpu-per-gpu ${MIN_VCPU:-8} --max-seconds $MAXS --yes-create"
CUR=$N
# poll the data centres once a minute until a pod is granted (the fourth try in the k=14 relaunch)
for try in $(seq 1 180); do
  if eval $T create --name $N $CREATE --datacenters US-NC-1,EU-FR-1 --pip '"numpy threadpoolctl cupy-cuda12x[ctk]"' --volume-gb 0 --disk-gb $DISK; then break; fi
  echo "no capacity (try $try), waiting 60 s"; sleep 60
  [ "$try" = 180 ] && { CUR=""; echo "no pod after 3 hours"; exit 1; }
done
# copy loop first: it terminates the pod when the chain says finished (or the job leaves "running")
# mirror.sh lives with the k = 14 relaunch scripts (fixed 2026-11-12: this line pointed at $HERE/mirror.sh, which does not exist)
nohup setsid "$ROOT/experiments/ic-bounds-k14/relaunch/mirror.sh" $N chain "$DEST" status:chain-status.json 300 $FILES > /dev/null 2>&1 < /dev/null &
wait_ping $N; $T sys $N
# every upload is retried: the first launch (2026-11-12 01:45 UTC) lost its pod to one HTTP 502 of the proxy after 50 uploads
push() { for i in 1 2 3 4 5 6; do $T push "$@" > /dev/null && return 0; echo "push retry $i: $*"; sleep 15; done; echo "push failed: $*"; return 1; }
for f in csd/ic_bounds_big.py csd/ic_bounds_bin_big.py scripts/ic_bounds_bin_big_run.py $EXP/pod/chain_bin.py; do push $N "$f"; done
push $N experiments/ic-bounds-k14/relaunch/nopool_dev_only.py nopool_run.py     # device pool off, pinned pool untouched
push $N $EXP/cells.json
for f in $EXP/instances/*.json; do push $N "$f"; done
echo "uploads done"
started=0
for i in 1 2 3 4 5 6; do
  if $T run $N chain -- python3 chain_bin.py --deadline $(( MAXS - 1800 )) "$@"; then started=1; break; fi
  # a failed answer does not prove the job did not start: look before trying again
  sleep 15; if $T jobs $N 2>/dev/null | grep -q '^chain '; then started=1; break; fi
done
[ "$started" = 1 ] || { echo "chain did not start"; false; }
CUR=""
echo "chain started on $N; records are copied to $DEST by mirror.sh (log: $DEST/mirror-$N.log)"

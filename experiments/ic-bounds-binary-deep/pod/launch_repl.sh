#!/bin/bash
# Launch the replication of ../REPLICATION.md on ONE H200 pod.  Run by the main session only.
# Needs RUNPOD_API_KEY and POD_STATE_DIR.  usage: launch_repl.sh MAX_SECONDS [DISK_GB]
set -eu
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../../.." && pwd); cd "$ROOT"
MAXS=${1:?usage: launch_repl.sh MAX_SECONDS [DISK_GB]}; DISK=${2:-60}
T="python3 scripts/ic_bounds_big_pod.py"; EXP=experiments/ic-bounds-binary-deep; DEST="$ROOT/$EXP/fetched-repl"
: "${POD_STATE_DIR:?set POD_STATE_DIR}"; mkdir -p "$POD_STATE_DIR" "$DEST"; chmod 700 "$POD_STATE_DIR"
N=${POD_NAME:-icb-bin-repl}; CUR=""
cleanup() { [ -n "$CUR" ] && { id=$(python3 -c "import json;print(json.load(open('$POD_STATE_DIR/$CUR.json'))['id'])" 2>/dev/null) && [ -n "$id" ] && $T terminate "$id" --yes; echo "launch of $CUR failed; pod terminated"; }; }
trap cleanup ERR
wait_ping() { for i in $(seq 1 60); do $T ping "$1" > /dev/null 2>&1 && return 0; sleep 20; done; echo "pod $1 did not answer in 20 min"; return 1; }
CUR=$N
for try in $(seq 1 180); do
  if $T create --name $N --gpu "NVIDIA H200" --min-ram-per-gpu ${MIN_RAM:-64} --min-vcpu-per-gpu ${MIN_VCPU:-16} --max-seconds $MAXS --yes-create \
       --datacenters US-NC-1,EU-FR-1 --pip "numpy threadpoolctl cupy-cuda12x[ctk]" --volume-gb 0 --disk-gb $DISK; then break; fi
  echo "no capacity (try $try), waiting 60 s"; sleep 60
  [ "$try" = 180 ] && { CUR=""; echo "no pod after 3 hours"; exit 1; }
done
nohup setsid "$ROOT/experiments/ic-bounds-k14/relaunch/mirror.sh" $N chain "$DEST" status:chain-status.json 180 \
  out/chain-status.json out/repl-r2.tgz out/repl-r3.tgz > /dev/null 2>&1 < /dev/null &
wait_ping $N; $T sys $N
push() { for i in 1 2 3 4 5 6; do $T push "$@" > /dev/null && return 0; echo "push retry $i: $*"; sleep 15; done; echo "push failed: $*"; return 1; }
for f in csd/ic_bounds_big.py csd/ic_bounds_bin_big.py scripts/ic_bounds_bin_big_run.py $EXP/pod/chain_bin.py $EXP/pod/repl_bin.py; do push $N "$f"; done
push $N experiments/ic-bounds-k14/relaunch/nopool_dev_only.py nopool_run.py
push $N $EXP/cells.json; push $N $EXP/cells_cpu.json
for f in $EXP/instances/*.json; do push $N "$f"; done
echo "uploads done"
started=0
for i in 1 2 3 4 5 6; do
  if $T run $N chain -- python3 repl_bin.py --deadline $(( MAXS - 1800 )); then started=1; break; fi
  sleep 15; if $T jobs $N 2>/dev/null | grep -q '^chain '; then started=1; break; fi
done
[ "$started" = 1 ] || { echo "chain did not start"; false; }
CUR=""
echo "replication started on $N; archives are copied to $DEST by mirror.sh (log: $DEST/mirror-$N.log)"

#!/bin/bash
# Launch the runs of ../RELAUNCH.md: a GPU pod for the k = 14 chain and a small CPU pod for the reference engine.
# Needs RUNPOD_API_KEY and POD_STATE_DIR (a directory that survives the session; tokens are written there).
# usage: launch.sh gpu|ref|both
set -eu
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/../../.." && pwd); cd "$ROOT"
T="python3 scripts/ic_bounds_big_pod.py"; DEST="$ROOT/experiments/ic-bounds-k14/relaunch/fetched"
: "${POD_STATE_DIR:?set POD_STATE_DIR}"; mkdir -p "$POD_STATE_DIR" "$DEST"; chmod 700 "$POD_STATE_DIR"
CUR=""
cleanup() { [ -n "$CUR" ] && { id=$(python3 -c "import json;print(json.load(open('$POD_STATE_DIR/$CUR.json'))['id'])" 2>/dev/null) && [ -n "$id" ] && $T terminate "$id" --yes; echo "launch of $CUR failed; pod terminated"; }; }
trap cleanup ERR   # a pod whose job did not start must not idle until its lifetime cap
wait_ping() { for i in $(seq 1 60); do $T ping "$1" > /dev/null 2>&1 && return 0; sleep 20; done; echo "pod $1 did not answer in 20 min"; return 1; }
gpu() {
  N=icb-k14b; CUR=$N
  # a volume keeps /workspace if the pod stops at its lifetime cap; fall back to the container disk if it is refused
  G='--gpu "NVIDIA H200" --min-ram-per-gpu 200 --min-vcpu-per-gpu 16 --max-seconds 27000 --yes-create'
  eval $T create --name $N $G --pip '"numpy threadpoolctl cupy-cuda12x[ctk]"' --volume-gb 250 --disk-gb 40 \
    || eval $T create --name $N $G --pip '"numpy threadpoolctl cupy-cuda12x[ctk]"' --volume-gb 0 --disk-gb 250
  nohup "$HERE/mirror.sh" $N chain "$DEST" status:chain-status.json 300 \
     out/chain-status.json out/selftest-gpu.json out/lift-selftest-gpu.json \
     out/proper-p114689-k14-D6.json out/proper-p114689-k14-D6.status.json out/lift-proper-p114689-k14-D6to7.json \
     out/lift-proper-p114689-k14-D6to7.status.json out/certcheck-proper-p114689-k14-D6.json \
     out/unit-p163841-k14-D6.json out/unit-p163841-k14-D6.status.json out/lift-unit-p163841-k14-D6to7.json \
     out/lift-unit-p163841-k14-D6to7.status.json out/certcheck-unit-p163841-k14-D6.json \
     logs/chain-selftest-engine.log logs/chain-selftest-lift.log logs/chain-lift-proper-p114689-k14-try1.log \
     logs/chain-lift-proper-p114689-k14-try2.log logs/chain-lift-unit-p163841-k14-try1.log logs/chain-lift-unit-p163841-k14-try2.log \
     > /dev/null 2>&1 &
  wait_ping $N; $T sys $N
  for f in csd/ic_bounds_big.py csd/ic_bounds_lift.py scripts/ic_bounds_big_run.py scripts/ic_bounds_lift_run.py scripts/ic_bounds_cert_check.py \
           experiments/ic-bounds-k14/relaunch/nopool_run.py experiments/ic-bounds-k14/relaunch/chain_k14.py \
           experiments/ic-bounds-k14/instances/proper-p114689-s22-k14.json experiments/ic-bounds-k14/instances/unit-p163841-a3-s21-k14.json; do
    $T push $N "$f"; done
  $T run $N chain -- python3 chain_k14.py; CUR=""
}
ref() {
  N=icb-ref11; CUR=$N
  $T create --name $N --flavors cpu3m,cpu5m --vcpu 8 --pip "numpy threadpoolctl" --volume-gb 0 --disk-gb 20 --max-seconds 36000 --yes-create
  nohup "$HERE/mirror.sh" $N ref "$DEST" file:reference-engine-unit-p65537-k11-D6.json 600 \
     out/reference-engine-unit-p65537-k11-D6.json > /dev/null 2>&1 &
  wait_ping $N; $T sys $N
  $T push $N csd/ic_bounds_d5_verify.py refeng_d5_verify.py
  $T push $N experiments/ic-bounds-k14/relaunch/ref_k11.py
  $T run $N ref -- python3 ref_k11.py 11 6; CUR=""
}
case "${1:-}" in gpu) gpu;; ref) ref;; both) ref; gpu;; *) echo "usage: launch.sh gpu|ref|both"; exit 2;; esac

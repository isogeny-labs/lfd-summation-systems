#!/usr/bin/env bash
# k = 14, D = 6 -> 7 reproduction of the recorded lifts (run on a pod; NOT run locally).
# Input: the engine's hky.basis/v1 export of the closed V_6(14) (+ its .json sidecar, which names the
# instance; the lift checks g and random m g against the basis and refuses otherwise).
# Usage: pod_k14.sh <hky_lift binary> <V_6 basis file> <recorded lift JSON> <out dir> [extra lift flags]
#   recorded lift JSON (lfd-summation-systems experiments/ic-bounds-k14/runs/):
#     lift-unit-p65537-k14-D6to7.json     expect 1 in V_7, dim V_7 = 1,683,218 (full space)
#     lift-unit-p163841-k14-D6to7.json    expect 1 in V_7, dim V_7 = 1,683,218
#     lift-proper-p114689-k14-D6to7.json  expect 1 not in V_7, dim V_7 = 1,680,848
# Exit: 0 all recorded fields equal; 3 a field differs; 2 the input was refused.
set -uo pipefail
BIN=${1:?binary}; BASIS=${2:?basis}; REC=${3:?recorded lift json}; OUT=${4:?out dir}; shift 4
mkdir -p "$OUT"
name=$(basename "$REC" .json)
sha256sum "$BIN" "$BASIS" "$REC" | tee "$OUT/$name.inputs.sha256"
nproc > "$OUT/$name.nproc"; free -g > "$OUT/$name.free"
# --cache-gb: follower cache (one follower <= #N x dim T x 4 bytes, about 12 GB at k = 14);
# memory floor without it: B (rank x #N x 4 = 90 GB) + K (<= dim T^2 x 4 = 15 GB) + one batch.
TIMER=""; [ -x /usr/bin/time ] && TIMER="/usr/bin/time -v"
$TIMER "$BIN" lift --basis "$BASIS" --expect "$REC" --out "$OUT/$name.rust.json" \
   --cache-gb "${CACHE_GB:-120}" --batch-tops "${BATCH_TOPS:-8192}" "$@" \
   > "$OUT/$name.stdout" 2> "$OUT/$name.time"
rc=$?
echo "exit $rc"
jq '.verify' "$OUT/$name.rust.json"
# the same comparison, standalone:
"$BIN" verify --got "$OUT/$name.rust.json" --expected "$REC" > "$OUT/$name.verify.json"; echo "verify exit $?"
exit $rc

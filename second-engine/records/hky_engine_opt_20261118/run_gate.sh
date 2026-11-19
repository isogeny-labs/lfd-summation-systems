#!/bin/bash
# Equality gate of the optimisation round (copy of ../hky_engine_speed2_20261117/run_gate.sh):
# the definitional reference (--engine hky of bin/pilot_base, built from the base commit) against the
# candidate (--engine hky3) on the same 110 cells; faults 1-5 are the negative controls (4, 5 attack the
# epoch coordinates and the sweep; they run with --hky3-sweep-ratio 0 so the sweep runs on these small cells).
# usage: run_gate.sh OUTDIR [extra hky3 flags...]   (exit 0 only if every set is EQUAL and every fault is detected)
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
REF=${REF:-$HERE/bin/pilot_base}
CAND=${CAND:-$HERE/bin/pilot_opt}
EQ=${EQ:-$HERE/bin/hky_engine_equality}
OUT=$1; shift
mkdir -p "$OUT"
P1=2013265921; P0=786433
common=(--ladder-t none --budget 1800 --hky-max-cols 2000000 --hky-max-entries 300000000)
EXTRA=("$@"); fail=0
run_set() {
  local name=$1; shift
  if [ ! -s "$OUT/$name.ref.jsonl" ]; then
    "$REF" "${common[@]}" "$@" --engine hky --out "$OUT/$name.ref.jsonl" 2> "$OUT/$name.ref.log" || { echo "$name: reference failed"; fail=1; }
  fi
  "$CAND" "${common[@]}" "$@" --engine hky3 "${EXTRA[@]}" --out "$OUT/$name.cand.jsonl" 2> "$OUT/$name.cand.log" || { echo "$name: candidate failed"; fail=1; }
  printf '%-22s ' "$name"
  "$EQ" compare "$OUT/$name.ref.jsonl" "$OUT/$name.cand.jsonl" | tail -1
  [ "${PIPESTATUS[0]}" -eq 0 ] || fail=1
}
run_set m2-tower-t1-8      --p $P1 --planted 0 --random 2 --hky-cmax 7 --kinds kummer --m 2 --controls tower --t-min 1 --max-t 8
run_set m2-tower-planted   --p $P1 --planted 2 --random 0 --hky-cmax 5 --kinds kummer --m 2 --controls tower --t-min 2 --max-t 7
run_set m2-null-t2-7       --p $P1 --planted 1 --random 2 --hky-cmax 6 --kinds kummer --m 2 --controls null --t-min 2 --max-t 7
run_set m2-naive-t2-4      --p $P1 --planted 1 --random 1 --hky-cmax 6 --kinds kummer --m 2 --controls naive --t-min 2 --max-t 4
run_set m3-tower-t2-3      --p $P1 --planted 1 --random 2 --hky-cmax 5 --kinds kummer --m 3 --controls tower --t-min 2 --max-t 3
run_set m4-tower-t2-3      --p $P1 --planted 0 --random 1 --hky-cmax 7 --kinds kummer --m 4 --controls tower --t-min 2 --max-t 3
run_set m5-tower-t2        --p $P1 --planted 0 --random 1 --hky-cmax 6 --kinds kummer --m 5 --controls tower --t-min 2 --max-t 2
run_set m6-tower-t2        --p $P1 --planted 0 --random 1 --hky-cmax 6 --kinds kummer --m 6 --controls tower --t-min 2 --max-t 2
run_set m7-tower-t2        --p $P1 --planted 0 --random 1 --hky-cmax 6 --kinds kummer --m 7 --controls tower --t-min 2 --max-t 2
run_set p0-kummer-isogeny  --p $P0 --planted 1 --random 2 --hky-cmax 6 --kinds kummer,isogeny --m 2 --controls tower --t-min 1 --max-t 6
common=(--ladder-t none --budget 1800 --hky-max-cols 1500 --hky-max-entries 300000000)
run_set censor-max-cols    --p $P1 --planted 0 --random 1 --hky-cmax 7 --kinds kummer --m 2 --controls tower --t-min 5 --max-t 7
common=(--ladder-t none --budget 1800 --hky-max-cols 2000000 --hky-max-entries 300000000)
for fault in 1 2 3 4 5; do
  "$CAND" "${common[@]}" --p $P1 --planted 2 --random 0 --hky-cmax 5 --kinds kummer --m 2 --controls tower --t-min 2 --max-t 7 \
     --engine hky3 "${EXTRA[@]}" --hky3-fault $fault --hky3-interreduce 0.05 --hky3-interreduce-min 1 --hky3-sweep-ratio 0 --out "$OUT/neg-fault$fault.cand.jsonl" 2> "$OUT/neg-fault$fault.cand.log"
  printf '%-22s ' "neg-fault$fault"
  "$EQ" compare --expect-mismatch "$OUT/m2-tower-planted.ref.jsonl" "$OUT/neg-fault$fault.cand.jsonl" | tail -1
  [ "${PIPESTATUS[0]}" -eq 0 ] || fail=1
done
if [ $fail -eq 0 ]; then echo "GATE PASS"; else echo "GATE FAIL"; fi
exit $fail

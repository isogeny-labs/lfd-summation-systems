#!/bin/bash
# Equality gate of the skip round: research/hky_engine_opt_20261118/run_gate.sh (110 reference cells, the
# definitional reference `--engine hky` of bin/pilot_base against `--engine hky3` of the candidate), with the
# reference records of that round reused (raw/gate-opt/*.ref.jsonl.gz), faults 1-5 as before, and the skip's
# negative controls: faults 6, 7, 8 (skips the proof does not cover), each also checked by the closedness check.
# usage: run_gate.sh OUTDIR [extra hky3 flags...]   e.g. run_gate.sh raw/gate-skip --hky3-skip on --hky3-verify-closed on
# Exit 0 only if every set is EQUAL, every fault 1-5 is detected, and the closedness check of the candidate
# runs found 0 non-zero reductions. Faults 6-8 are reported (detected by equality and/or by the check).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
REFDIR=${REFDIR:-$HERE/../hky_engine_opt_20261118/raw/gate-opt}
REF=${REF:-$HERE/bin/pilot_base}
CAND=${CAND:-$HERE/bin/pilot_skip}
EQ=${EQ:-$HERE/bin/hky_engine_equality}
OUT=$1; shift
mkdir -p "$OUT"
P1=2013265921; P0=786433
common=(--ladder-t none --budget 1800 --hky-max-cols 2000000 --hky-max-entries 300000000)
EXTRA=("$@"); fail=0
verify_sum() {  # closures without 1 checked, reductions, non-zero
  jq -s '[.[] | select(.hky3 != null) | .hky3.stats[]? | select(.verify_products > 0)] |
         "\(length) \(map(.verify_products) | add // 0) \(map(.verify_nonzero) | add // 0) \(map(.skipped_items) | add // 0) \(map(.queued_items) | add // 0)"' -r "$1"
}
run_set() {
  local name=$1; shift
  if [ ! -s "$OUT/$name.ref.jsonl" ]; then
    if [ -s "$REFDIR/$name.ref.jsonl.gz" ]; then gunzip -c "$REFDIR/$name.ref.jsonl.gz" > "$OUT/$name.ref.jsonl"
    else "$REF" "${common[@]}" "$@" --engine hky --out "$OUT/$name.ref.jsonl" 2> "$OUT/$name.ref.log" || { echo "$name: reference failed"; fail=1; }; fi
  fi
  "$CAND" "${common[@]}" "$@" --engine hky3 "${EXTRA[@]}" --out "$OUT/$name.cand.jsonl" 2> "$OUT/$name.cand.log" || { echo "$name: candidate failed"; fail=1; }
  printf '%-22s ' "$name"
  "$EQ" compare "$OUT/$name.ref.jsonl" "$OUT/$name.cand.jsonl" | tail -1
  [ "${PIPESTATUS[0]}" -eq 0 ] || fail=1
  read -r n v bad sk q <<<"$(verify_sum "$OUT/$name.cand.jsonl")"
  printf '%-22s closedness: %s closures, %s reductions, %s non-zero; skipped %s of %s items\n' "" "$n" "$v" "$bad" "$sk" "$q"
  echo "$name $n $v $bad $sk $q" >> "$OUT/verify.txt"
  [ "$bad" = 0 ] || fail=1
}
rm -f "$OUT/verify.txt"
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
for fault in 1 2 3 4 5 6 7 8; do
  "$CAND" "${common[@]}" --p $P1 --planted 2 --random 0 --hky-cmax 5 --kinds kummer --m 2 --controls tower --t-min 2 --max-t 7 \
     --engine hky3 "${EXTRA[@]}" --hky3-fault $fault --hky3-interreduce 0.05 --hky3-interreduce-min 1 --hky3-sweep-ratio 0 --out "$OUT/neg-fault$fault.cand.jsonl" 2> "$OUT/neg-fault$fault.cand.log"
  printf '%-22s ' "neg-fault$fault"
  "$EQ" compare --expect-mismatch "$OUT/m2-tower-planted.ref.jsonl" "$OUT/neg-fault$fault.cand.jsonl" | tail -1
  ok=${PIPESTATUS[0]}
  read -r n v bad sk q <<<"$(verify_sum "$OUT/neg-fault$fault.cand.jsonl")"
  printf '%-22s closedness: %s closures, %s reductions, %s non-zero; skipped %s of %s items\n' "" "$n" "$v" "$bad" "$sk" "$q"
  echo "neg-fault$fault $n $v $bad $sk $q" >> "$OUT/verify.txt"
  if [ "$fault" -le 5 ] && [ "$ok" -ne 0 ]; then fail=1; fi
done
if [ $fail -eq 0 ]; then echo "GATE PASS"; else echo "GATE FAIL"; fi
exit $fail

#!/bin/bash
# One no-1 proxy closure: cell1.sh BIN T C OUTPREFIX [extra flags...]
# m = 2, Kummer tower, draw 0 of the pilot (random 1, planted 0), only the closure c = C.
# Writes OUTPREFIX.jsonl, OUTPREFIX.log, OUTPREFIX.rss (peak RSS kB from /proc via wait) and
# prints one JSON line with wall, cpu, peak RSS, load and the deterministic counters.
set -u
BIN=$1; T=$2; C=$3; OUT=$4; shift 4
rm -f "$OUT.jsonl"
l0=$(cut -d' ' -f1 /proc/loadavg); s0=$(date +%s.%N)
"$BIN" --p 2013265921 --ladder-t none --planted ${PLANTED:-0} --random ${RANDOM_N:-1} --hky-cmin $C --hky-cmax $C --budget ${BUDGET:-200000} \
  --hky-max-cols 3000000 --hky-max-entries 140000000000 --kinds kummer --m 2 --controls tower \
  --t-min $T --max-t $T --engine hky3 "$@" --out "$OUT.jsonl" 2> "$OUT.log" &
pid=$!; peak=0; cpu=0
while kill -0 $pid 2>/dev/null; do
  r=$(awk '/VmHWM/{print $2}' /proc/$pid/status 2>/dev/null); [ -n "$r" ] && peak=$r
  st=$(awk '{print $14+$15}' /proc/$pid/stat 2>/dev/null); [ -n "$st" ] && cpu=$st
  sleep 0.5
done
wait $pid; rc=$?; s1=$(date +%s.%N)
tck=$(getconf CLK_TCK)
jq -c --arg l0 "$l0" --arg l1 "$(cut -d' ' -f1 /proc/loadavg)" --argjson rc $rc --argjson wall "$(echo "$s1 - $s0" | bc)" \
  --argjson peak "$peak" --argjson cpu "$(echo "scale=2; $cpu / $tck" | bc)" --arg bin "$(basename $BIN)" --arg thr "${RAYON_NUM_THREADS:-all}" \
  'select(has("closures")) | {bin: $bin, threads: $thr, t, rc: $rc, wall_s: $wall, cpu_s_approx: $cpu, peak_rss_kb: $peak, load: [$l0, $l1],
    cl: [.closures[] | {c, cols, one, dims, cen: .censored, muladds, rr: .rows_reduced, entries, ms: (.ms|round)}],
    st: [.hky3.stats[]? | {parallel_s, second_s, serial_s, interreduce_s, batches, rounds, interreductions, interreduce_muladds, second_muladds, products, zero, reducers, stored_words, peak_words, dense_rows, sparse_rows}]}' "$OUT.jsonl"

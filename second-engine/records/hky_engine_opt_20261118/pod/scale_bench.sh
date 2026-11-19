#!/bin/bash
# Thread-scaling benchmark of the hky3 closure engine, old (base) against new (this round), on a pod.
# Prepared in research/hky_engine_opt_20261118 (not run there: 4 shared vCPU). bash + coreutils only
# (no python, no jq). One JSON line per run in $OUT/runs.jsonl; the pilot's own record line is embedded
# verbatim under "pilot" (closures, dims, muladds, rows_reduced, entries, and the hky3 stats).
#
# usage: scale_bench.sh OLD_BIN NEW_BIN [OUT_DIR]
#   OLD_BIN  frozen base pilot (pkm_tower_pilot of the base commit, sha256 in ../README.md)
#   NEW_BIN  frozen optimised pilot (this round)
# environment (defaults in brackets):
#   LADDER_CELL  cell run at every thread count                    ["10 6 planted"]
#   MAIN_CELL    the 5-15 min cell, run at the large thread counts  ["12 6 random"]
#   LADDER_THREADS ["1 2 4 8 16 32", plus 64 if nproc >= 64]
#   MAIN_THREADS   ["16 32", plus 64 if nproc >= 64]
#   REPS [1]   CAP_S (per run, timeout) [3600]   OLD_FLAGS [""]   NEW_FLAGS [""]
#   DRY=1 prints the plan and exits.
# A cell is "T C KIND": m = 2, Kummer tower, p = 2013265921, the pilot's draw 0; KIND = planted (the
# consistent system: 1 is in no closure, the closure runs to the end) or random (the study draw).
# Only the closure c = C is computed (--hky-cmin C --hky-cmax C).
# Default cells: the ladder is planted t = 10 c = 6 (60,460 columns, no 1, rank n - 2: a full closure; 3.64e11
# muladds, 228 s base / 134 s new on 4 shared vCPU, so ~8 and ~4 min at one thread on a pod); the main cell is the
# study draw t = 12 c = 6 (190,051 columns; 1 enters only as the last pivot, so ~all of the full closure is built;
# ~2..4e12 muladds, an estimated 5..15 min at 32 threads). Planted t = 11 and 12 c = 6 were not chosen: the planted
# work grew x72 from t = 9 to t = 10, so they would take hours. Estimated total with the defaults: ~1 h.
#
# Order: for each rep, ladder then main; thread counts ascending; the two arms alternate and their
# order flips every rep and every thread count (ABBA), so drift in the host's load hits both arms.
#
# What to read: the counters (muladds, rows_reduced, entries, dims) must be identical between the
# arms and across thread counts (determinism); wall_s gives the scaling; util = cpu/(wall*threads).
set -u
OLD=${1:?OLD_BIN}; NEW=${2:?NEW_BIN}; OUT=${3:-scale_out_$(date -u +%Y%m%dT%H%M%SZ)}
NPROC=$(nproc)
extra=""; [ "$NPROC" -ge 64 ] && extra=" 64"
LADDER_CELL=${LADDER_CELL:-"10 6 planted"}
MAIN_CELL=${MAIN_CELL:-"12 6 random"}
LADDER_THREADS=${LADDER_THREADS:-"1 2 4 8 16 32$extra"}
MAIN_THREADS=${MAIN_THREADS:-"16 32$extra"}
REPS=${REPS:-1}; CAP_S=${CAP_S:-3600}
OLD_FLAGS=${OLD_FLAGS:-}; NEW_FLAGS=${NEW_FLAGS:-}
mkdir -p "$OUT/tmp"
RUNS="$OUT/runs.jsonl"

sha() { sha256sum "$1" | cut -d' ' -f1; }
SHA_OLD=$(sha "$OLD"); SHA_NEW=$(sha "$NEW")
CPU_MODEL=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ *//; s/"/'"'"'/g')
MEM_KB=$(awk '/MemTotal/{print $2}' /proc/meminfo)
HOST_JSON="{\"nproc\":$NPROC,\"cpu\":\"$CPU_MODEL\",\"mem_kb\":$MEM_KB,\"kernel\":\"$(uname -r)\",\"hostname\":\"$(hostname)\"}"
printf '{"bench":"hky3-scale/v1","event":"start","utc":"%s","host":%s,"old":{"path":"%s","sha256":"%s","flags":"%s"},"new":{"path":"%s","sha256":"%s","flags":"%s"},"ladder_cell":"%s","main_cell":"%s","ladder_threads":"%s","main_threads":"%s","reps":%s,"cap_s":%s}\n' \
  "$(date -u +%FT%TZ)" "$HOST_JSON" "$OLD" "$SHA_OLD" "$OLD_FLAGS" "$NEW" "$SHA_NEW" "$NEW_FLAGS" \
  "$LADDER_CELL" "$MAIN_CELL" "$LADDER_THREADS" "$MAIN_THREADS" "$REPS" "$CAP_S" | tee -a "$RUNS"

# children's cumulative user and system CPU seconds of this shell (bash `times`, second line; `times`
# must run in this shell, not in a $( ) subshell, which has no children of its own)
child_cpu() { times > "$OUT/tmp/.times"; tail -1 "$OUT/tmp/.times" | sed 's/m/ /g; s/s//g' | awk '{printf "%.3f %.3f", $1*60+$2, $3*60+$4}'; }

run_one() { # run_one ARM BIN FLAGS CELL THREADS REP
  local arm=$1 bin=$2 flags=$3 cell=$4 thr=$5 rep=$6
  local t c kind; read -r t c kind <<< "$cell"
  local planted=0 random=1; [ "$kind" = planted ] && { planted=1; random=0; }
  local tag="${arm}_t${t}c${c}${kind}_T${thr}_r${rep}"
  local out="$OUT/tmp/$tag.jsonl"; rm -f "$out"
  local l0; l0=$(cut -d' ' -f1-3 /proc/loadavg)
  child_cpu > "$OUT/tmp/.cpu0"; local cpu0; cpu0=$(cat "$OUT/tmp/.cpu0")
  local s0; s0=$(date +%s.%N)
  # shellcheck disable=SC2086
  RAYON_NUM_THREADS=$thr "$bin" --p 2013265921 --ladder-t none --planted $planted --random $random \
    --hky-cmin "$c" --hky-cmax "$c" --budget 400000 --hky-max-cols 3000000 --hky-max-entries 140000000000 \
    --kinds kummer --m 2 --controls tower --t-min "$t" --max-t "$t" --engine hky3 $flags \
    --out "$out" 2> "$OUT/tmp/$tag.log" &
  local pid=$! peak=0 r capped=0
  while kill -0 $pid 2>/dev/null; do
    if [ "$(awk -v a="$s0" -v b="$(date +%s.%N)" -v c="$CAP_S" 'BEGIN{print (b-a > c)}')" = 1 ]; then kill $pid; capped=1; fi
    r=$(awk '/VmHWM/{print $2}' /proc/$pid/status 2>/dev/null) && [ -n "$r" ] && peak=$r
    sleep 0.2
  done
  wait $pid; local rc=$?
  local s1; s1=$(date +%s.%N)
  child_cpu > "$OUT/tmp/.cpu1"; local cpu1; cpu1=$(cat "$OUT/tmp/.cpu1")
  local l1; l1=$(cut -d' ' -f1-3 /proc/loadavg)
  local wall; wall=$(awk -v a="$s0" -v b="$s1" 'BEGIN{printf "%.3f", b-a}')
  local cu cs
  read -r cu cs <<< "$(awk -v x="$cpu0" -v y="$cpu1" 'BEGIN{split(x,a," "); split(y,b," "); printf "%.3f %.3f", b[1]-a[1], b[2]-a[2]}')"
  local util; util=$(awk -v u="$cu" -v s="$cs" -v w="$wall" -v t="$thr" 'BEGIN{ if (w>0) printf "%.4f", (u+s)/(w*t); else print 0 }')
  local pil; pil=$(grep -m1 '"closures"' "$out" 2>/dev/null); [ -z "$pil" ] && pil=null
  printf '{"bench":"hky3-scale/v1","event":"run","utc":"%s","arm":"%s","sha256":"%s","flags":"%s","cell":{"m":2,"t":%s,"c":%s,"kind":"%s"},"threads":%s,"rep":%s,"rc":%s,"timed_out":%s,"wall_s":%s,"cpu_user_s":%s,"cpu_sys_s":%s,"util":%s,"peak_rss_kb":%s,"load_before":"%s","load_after":"%s","pilot":%s}\n' \
    "$(date -u +%FT%TZ)" "$arm" "$([ "$arm" = old ] && echo "$SHA_OLD" || echo "$SHA_NEW")" "$flags" "$t" "$c" "$kind" \
    "$thr" "$rep" "$rc" "$([ $capped -eq 1 ] && echo true || echo false)" "$wall" "$cu" "$cs" "$util" "$peak" "$l0" "$l1" "$pil" >> "$RUNS"
  printf '%s T=%-3s rep %s: wall %8.1f s  cpu %9.1f s  util %s  rss %s kB  rc %s\n' "$arm" "$thr" "$rep" "$wall" "$(awk -v u="$cu" -v s="$cs" 'BEGIN{print u+s}')" "$util" "$peak" "$rc"
}

plan=()
for rep in $(seq 1 "$REPS"); do
  i=0
  for which in ladder main; do
    if [ $which = ladder ]; then cell=$LADDER_CELL; ths=$LADDER_THREADS; else cell=$MAIN_CELL; ths=$MAIN_THREADS; fi
    for thr in $ths; do
      if [ $(( (i + rep) % 2 )) -eq 0 ]; then plan+=("old|$cell|$thr|$rep" "new|$cell|$thr|$rep")
      else plan+=("new|$cell|$thr|$rep" "old|$cell|$thr|$rep"); fi
      i=$((i + 1))
    done
  done
done
printf 'plan: %s runs\n' "${#plan[@]}"; for x in "${plan[@]}"; do echo "  $x"; done
[ "${DRY:-0}" = 1 ] && exit 0
for x in "${plan[@]}"; do
  IFS='|' read -r arm cell thr rep <<< "$x"
  if [ "$arm" = old ]; then run_one old "$OLD" "$OLD_FLAGS" "$cell" "$thr" "$rep"
  else run_one new "$NEW" "$NEW_FLAGS" "$cell" "$thr" "$rep"; fi
done
printf '{"bench":"hky3-scale/v1","event":"end","utc":"%s"}\n' "$(date -u +%FT%TZ)" >> "$RUNS"
echo "done: $RUNS"

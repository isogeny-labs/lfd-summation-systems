#!/bin/bash
# One paper closure job on a pod. Needs bash + coreutils only (no python, no jq).
#
#   run_job.sh <job name> <instance JSON> <D>
#
# Environment (defaults): BIN=/workspace/bin/pkm_tower_pilot, SHA_FILE=<this dir>/BINARY.sha256,
# OUT_ROOT=/workspace/out, CKPT_MIN=30, MEM_SOFT_PCT=85 (SIGTERM: checkpoint + clean censored exit),
# MEM_HARD_PCT=97 (SIGKILL), SAMPLE_S=60, EXPORT=1, ENGINE_FLAGS (extra engine flags), VERIFY=0 (1: engine re-checks closedness at the end,
# doubling peak memory: J1-J3 on 256 GB, not J4), MEM_LIMIT_KB (test override of the detected limit).
#
# Writes /workspace/out/<job>/: record.json (the engine's record), job.json (this script's record:
# host facts, status, peak RSS, wall, rc), basis.bin + basis.bin.json (hky.basis/v1 and sidecar),
# rss.log (one sample per SAMPLE_S), engine.err, ckpt/ (state.ckpt, state.json). Re-running the same
# command after a container restart resumes from ckpt/state.ckpt automatically.
set -u
JOB=$1; INST=$2; D=$3
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=${BIN:-/workspace/bin/pkm_tower_pilot}
SHA_FILE=${SHA_FILE:-$HERE/BINARY.sha256}
OUT_ROOT=${OUT_ROOT:-/workspace/out}
CKPT_MIN=${CKPT_MIN:-30}
MEM_SOFT_PCT=${MEM_SOFT_PCT:-85}
MEM_HARD_PCT=${MEM_HARD_PCT:-97}
SAMPLE_S=${SAMPLE_S:-60}
EXPORT=${EXPORT:-1}
OUT="$OUT_ROOT/$JOB"; CK="$OUT/ckpt"
mkdir -p "$OUT" "$CK"
LOG="$OUT/job.log"
say() { echo "[$(date -u +%FT%TZ)] $*" | tee -a "$LOG"; }
jstr() { printf '"%s"' "$(printf '%s' "$1" | tr -d '\n' | sed 's/\\/\\\\/g; s/"/\\"/g')"; }

# ── 1. binary identity ────────────────────────────────────────────────
want=$(cut -d' ' -f1 "$SHA_FILE")
got=$(sha256sum "$BIN" | cut -d' ' -f1)
if [ "$want" != "$got" ]; then
  say "BINARY SHA256 MISMATCH: want $want got $got"; exit 2
fi
inst_sha=$(sha256sum "$INST" | cut -d' ' -f1)
say "job=$JOB instance=$INST ($inst_sha) D=$D binary=$got"

# ── 2. host facts ─────────────────────────────────────────────────────
cpu_model=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ *//')
nproc_v=$(nproc)
mem_total_kb=$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)
cg_max=""
if [ -r /sys/fs/cgroup/memory.max ]; then cg_max=$(cat /sys/fs/cgroup/memory.max)
elif [ -r /sys/fs/cgroup/memory/memory.limit_in_bytes ]; then cg_max=$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes); fi
limit_kb=$mem_total_kb
if [ -n "$cg_max" ] && [ "$cg_max" != "max" ] && [ "$cg_max" -lt $((mem_total_kb * 1024)) ] 2>/dev/null; then
  limit_kb=$((cg_max / 1024)); fi
limit_kb=${MEM_LIMIT_KB:-$limit_kb}   # test override
disk_free_kb=$(df -Pk "$OUT" | awk 'NR==2 {print $4}')
glibc=$(ldd --version 2>&1 | head -1)
kernel=$(uname -r)
flags=$(grep -m1 '^flags' /proc/cpuinfo | grep -o 'avx512f\|avx2' | sort -u | tr '\n' ' ')
say "host: cpu='$cpu_model' nproc=$nproc_v mem_total_kb=$mem_total_kb cgroup_max=$cg_max limit_kb=$limit_kb disk_free_kb=$disk_free_kb kernel=$kernel simd='$flags' glibc='$glibc'"
ldd "$BIN" >> "$LOG" 2>&1

# ── 3. run (resume if a checkpoint exists) ────────────────────────────
if [ -f "$OUT/job.json" ] && grep -q '"status": "ok"' "$OUT/job.json"; then
  say "already complete (job.json status ok); nothing to do"; exit 0
fi
args=(--instance "$INST" --D "$D" --checkpoint-every "$CKPT_MIN" --out "$OUT/record.json")
# ENGINE_FLAGS (engine-skip round): extra engine flags, e.g. "--hky3-skip on" (J4.sh); recorded in job.log
# shellcheck disable=SC2206
[ -n "${ENGINE_FLAGS:-}" ] && args+=($ENGINE_FLAGS)
if [ -s "$CK/state.ckpt" ]; then
  say "resuming from $CK/state.ckpt: $(cat "$CK/state.json" 2>/dev/null)"
  args+=(--resume "$CK")
else
  args+=(--checkpoint-dir "$CK")
fi
[ "$EXPORT" = 1 ] && args+=(--export-basis "$OUT/basis.bin")
[ "${VERIFY:-0}" = 1 ] && args+=(--hky3-verify-closed on)   # doubles peak memory
t0=$(date +%s)
"$BIN" "${args[@]}" > "$OUT/record.stdout" 2>> "$OUT/engine.err" &
pid=$!
say "pid=$pid args=${args[*]}"
peak_kb=0; guard=""
soft_kb=$((limit_kb * MEM_SOFT_PCT / 100)); hard_kb=$((limit_kb * MEM_HARD_PCT / 100))
last=0
while kill -0 $pid 2>/dev/null; do
  rss=$(awk '/^VmRSS:/ {print $2}' /proc/$pid/status 2>/dev/null || echo 0)
  hwm=$(awk '/^VmHWM:/ {print $2}' /proc/$pid/status 2>/dev/null || echo 0)
  rss=${rss:-0}; hwm=${hwm:-0}
  [ "$hwm" -gt "$peak_kb" ] && peak_kb=$hwm
  [ "$rss" -gt "$peak_kb" ] && peak_kb=$rss
  now=$(date +%s)
  if [ $((now - last)) -ge "$SAMPLE_S" ]; then
    echo "$now $rss $hwm $(cut -d' ' -f1-3 /proc/loadavg) $(df -Pk "$OUT" | awk 'NR==2 {print $4}')" >> "$OUT/rss.log"
    last=$now
  fi
  if [ -z "$guard" ] && [ "$rss" -gt "$soft_kb" ]; then
    guard="memory_guard_soft"; say "RSS ${rss} kB > ${MEM_SOFT_PCT}% of ${limit_kb} kB: SIGTERM (checkpoint, clean stop)"
    kill -TERM $pid
  elif [ "$rss" -gt "$hard_kb" ]; then
    guard="memory_guard_hard"; say "RSS ${rss} kB > ${MEM_HARD_PCT}% of ${limit_kb} kB: SIGKILL"
    kill -KILL $pid
  fi
  sleep 2
done
wait $pid; rc=$?
t1=$(date +%s)
# ── 4. status and record ──────────────────────────────────────────────
status=ok
if [ $rc -ne 0 ]; then
  if [ -n "$guard" ]; then status="censored_$guard"
  elif [ $rc -eq 137 ]; then status="censored_sigkill_likely_oom"
  elif [ $rc -eq 3 ]; then status="censored_engine"
  else status="failed_rc_$rc"; fi
fi
grep -q '"closed": true' "$OUT/record.json" 2>/dev/null || [ "$status" != ok ] || status="failed_not_closed"
bsha=""; [ -f "$OUT/basis.bin" ] && bsha=$(sha256sum "$OUT/basis.bin" | cut -d' ' -f1)
{
  printf '{\n "schema": "hky_paper_runs.pod_job/v1",\n "job": %s,\n "status": %s,\n "rc": %d,\n' "$(jstr "$JOB")" "$(jstr "$status")" "$rc"
  printf ' "instance": %s,\n "instance_sha256": %s,\n "D": %s,\n' "$(jstr "$INST")" "$(jstr "$inst_sha")" "$D"
  printf ' "binary_sha256": %s,\n "wall_s_this_attempt": %d,\n "peak_rss_kb_this_attempt": %d,\n' "$(jstr "$got")" $((t1 - t0)) "$peak_kb"
  printf ' "memory_limit_kb": %d,\n "guard": %s,\n "basis_sha256": %s,\n' "$limit_kb" "$(jstr "$guard")" "$(jstr "$bsha")"
  printf ' "host": {"cpu": %s, "nproc": %d, "mem_total_kb": %d, "cgroup_max": %s, "kernel": %s, "simd": %s},\n' \
    "$(jstr "$cpu_model")" "$nproc_v" "$mem_total_kb" "$(jstr "$cg_max")" "$(jstr "$kernel")" "$(jstr "$flags")"
  printf ' "checkpoint_meta": %s,\n "engine_record": "record.json"\n}\n' "$(cat "$CK/state.json" 2>/dev/null | tr -d '\n' || echo null)"
} > "$OUT/job.json.tmp"
# an empty checkpoint_meta would not be JSON
sed -i 's/"checkpoint_meta": ,/"checkpoint_meta": null,/' "$OUT/job.json.tmp"
mv "$OUT/job.json.tmp" "$OUT/job.json"
say "done status=$status rc=$rc wall=$((t1 - t0))s peak_rss_kb=$peak_kb basis_sha256=$bsha"
[ "$status" = ok ] && exit 0 || exit 1

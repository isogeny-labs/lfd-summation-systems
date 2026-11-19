#!/bin/bash
# Exact-instance gate: every recorded (instance, k, D) closure of the paper that fits locally is recomputed
# with the pilot's --instance loader and compared entry by entry (columns, profile, contains_one, closed).
# Usage: run_gate.sh <pilot binary> <out dir> [cell ...]   cell = <instance stem>:<D>:<runs stem>
# Writes <out>/rows.jsonl (one line per cell: expected, got, match, wall, muladds, binary sha256).
set -u
BIN=$1; OUT=$2; shift 2
PAPER=${PAPER:-/home/user/lfd-summation-systems/experiments/ic-bounds-k14}
mkdir -p "$OUT"
SHA=$(sha256sum "$BIN" | cut -d' ' -f1)
CELLS=("$@")
if [ ${#CELLS[@]} -eq 0 ]; then
  CELLS=(
    unit-p65537-a1-s10-k10:5:unit-p65537-k10-D5
    unit-p65537-a1-s10-k11:5:unit-p65537-k11-D5
    unit-p65537-a1-s10-k12:5:unit-p65537-k12-D5
    unit-p65537-a1-s10-k10:6:unit-p65537-k10-D6
    unit-p65537-a1-s10-k13:5:unit-p65537-k13-D5
    unit-p163841-a3-s21-k13:5:unit-p163841-k13-D5
    unit-p65537-a1-s10-k14:5:unit-p65537-k14-D5
    unit-p163841-a3-s21-k14:5:unit-p163841-k14-D5
    unit-p65537-a1-s10-k11:6:unit-p65537-k11-D6
    unit-p65537-a1-s10-k15:5:unit-p65537-k15-D5
    unit-p65537-a1-s10-k12:6:unit-p65537-k12-D6
    unit-p65537-a1-s10-k13:6:unit-p65537-k13-D6
    unit-p163841-a3-s21-k13:6:unit-p163841-k13-D6
  )
fi
for cell in "${CELLS[@]}"; do
  IFS=: read -r stem D run <<<"$cell"
  inst="$PAPER/instances/$stem.json"; rec="$PAPER/runs/$run.json"
  [ -f "$rec" ] || { echo "skip $cell: no record"; continue; }
  name="$stem-D$D"
  t0=$(date +%s.%N)
  # CKPT_ROOT set: checkpoint every CKPT_MIN minutes and resume after a restart (re-run the same command)
  ck=(); [ -n "${CKPT_ROOT:-}" ] && ck=(--resume "$CKPT_ROOT/$name" --checkpoint-every "${CKPT_MIN:-10}")
  # EXTRA_FLAGS (engine-skip round): extra engine flags, e.g. "--hky3-skip on --hky3-verify-closed on"
  # shellcheck disable=SC2086
  "$BIN" --instance "$inst" --D "$D" "${ck[@]}" ${EXTRA_FLAGS:-} --out "$OUT/$name.json" 2>> "$OUT/$name.err" > /dev/null
  rc=$?
  t1=$(date +%s.%N)
  exp=$(jq -c '{columns, profile, contains_one}' "$rec")
  got=$(jq -c '{columns, profile, contains_one}' "$OUT/$name.json" 2>/dev/null || echo null)
  match=false; [ "$exp" = "$got" ] && [ "$(jq -r .closed "$OUT/$name.json")" = true ] && match=true
  jq -nc --arg cell "$cell" --argjson exp "$exp" --argjson got "$got" --argjson match $match \
     --arg sha "$SHA" --argjson rc $rc --arg wall "$(echo "$t1 - $t0" | bc)" \
     --slurpfile r "$OUT/$name.json" --arg load "$(cut -d' ' -f1-3 /proc/loadavg)" \
     '{cell:$cell, expected:$exp, got:$got, match:$match, rc:$rc, wall_s:($wall|tonumber),
       muladds:$r[0].muladds, maxrss_mb:$r[0].maxrss_mb, skipped_items:$r[0].engine.skipped_items,
       queued_items:$r[0].engine.queued_items, verify_products:$r[0].engine.verify_products,
       verify_nonzero:$r[0].engine.verify_nonzero, extra_flags:"'"${EXTRA_FLAGS:-}"'", threads:$r[0].rayon_threads,
       engine:$r[0].engine.name, instance_sha256:$r[0].instance.sha256,
       generator_matches_file:$r[0].system.generator_matches_file, binary_sha256:$sha, load:$load}' \
     | tee -a "$OUT/rows.jsonl"
done

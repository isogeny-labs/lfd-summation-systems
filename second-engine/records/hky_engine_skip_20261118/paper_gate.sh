#!/bin/bash
# Paper gate of the skip round: the 11 recorded (instance, k, D) closures of
# research/hky_paper_runs_20261118/gate/final_gate.sh, recomputed by the frozen skip binary with
# EXTRA_FLAGS (default: --hky3-skip on --hky3-verify-closed on), compared with the paper's records
# (columns, profile, contains_one, closed). usage: paper_gate.sh BIN OUTDIR
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
BIN=$(readlink -f "$1"); OUT=$(readlink -m "$2")
export EXTRA_FLAGS=${EXTRA_FLAGS:---hky3-skip on --hky3-verify-closed on}
cd "$HERE/../hky_paper_runs_20261118" || exit 2
./gate/run_gate.sh "$BIN" "$OUT" unit-p65537-a1-s10-k10:5:unit-p65537-k10-D5 unit-p65537-a1-s10-k11:5:unit-p65537-k11-D5 \
  unit-p65537-a1-s10-k12:5:unit-p65537-k12-D5 unit-p65537-a1-s10-k10:6:unit-p65537-k10-D6 unit-p65537-a1-s10-k13:5:unit-p65537-k13-D5 \
  unit-p163841-a3-s21-k13:5:unit-p163841-k13-D5 unit-p65537-a1-s10-k14:5:unit-p65537-k14-D5 unit-p163841-a3-s21-k14:5:unit-p163841-k14-D5 \
  unit-p65537-a1-s10-k15:5:unit-p65537-k15-D5 unit-p65537-a1-s10-k11:6:unit-p65537-k11-D6 unit-p65537-a1-s10-k12:6:unit-p65537-k12-D6
jq -s -r '"PAPER GATE: \(length) cells, \(map(select(.match)) | length) match; closedness \(map(.verify_products) | add) reductions, \(map(.verify_nonzero) | add) non-zero; skipped \(map(.skipped_items) | add) of \(map(.queued_items) | add) items"' "$OUT/rows.jsonl"

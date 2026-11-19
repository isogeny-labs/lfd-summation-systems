#!/bin/bash
# Gate of the standalone engine extract (second-engine/engine) against the paper's records.
#
#   1. closures: the 11 recorded unit-instance closures that fit in 16 GB (k = 10..15 at D = 5,
#      k = 10..12 at D = 6), recomputed with --hky3-verify-closed on (closedness re-checked by the
#      engine), compared with experiments/ic-bounds-k14/runs/*.json (columns, profile,
#      contains_one, closed);
#   2. the same cells with --hky3-skip on --hky3-verify-closed on (the configuration of run J4s);
#   3. export + the independent checker hky_basis_check --closedness at k = 10, D = 5 and D = 6;
#   4. checkpoint/resume (SIGKILL and SIGTERM, then resume) at k = 10, D = 5;
#   5. negative controls at k = 10, D = 6, compared with the record's
#      gate/negative_controls_d6.log: the planted target X = 11932 (36 scalar solutions) must change
#      the result; a perturbed X and a wrong a give a unit instance again (same profile, as recorded).
#
# STEPS (default "1 2 3 4 5") selects the steps; the output directory is shared.
#
# The closure and checkpoint steps call the record's own scripts
# (records/hky_paper_runs_20261118/gate/run_gate.sh, ckpt_test.sh) unchanged.
# Usage: [STEPS="…"] extract_gate.sh OUTDIR [SCRATCH]   (needs cargo, jq, bc; CARGO_TARGET_DIR is honoured)
set -u
STEPS=" ${STEPS:-1 2 3 4 5} "
step() { [[ "$STEPS" == *" $1 "* ]]; }
HERE=$(cd "$(dirname "$0")" && pwd)
ENGINE=$(cd "$HERE/.." && pwd)
ROOT=$(cd "$ENGINE/../.." && pwd)
OUT=$(readlink -m "${1:?OUTDIR}"); SCR=$(readlink -m "${2:-${TMPDIR:-/tmp}/hky3-extract-gate}")
mkdir -p "$OUT" "$SCR"
export PAPER="$ROOT/experiments/ic-bounds-k14"
GATE="$ROOT/second-engine/records/hky_paper_runs_20261118/gate"
(cd "$ENGINE" && cargo build --release --examples) || exit 2
TD=${CARGO_TARGET_DIR:-$ENGINE/target}
BIN="$TD/release/examples/pkm_tower_pilot"; CHK="$TD/release/examples/hky_basis_check"
sha256sum "$BIN" "$CHK" | sed "s|$TD/||" > "$OUT/binaries.sha256.$$"
mv "$OUT/binaries.sha256.$$" "$OUT/binaries.sha256"
{ cargo --version; rustc --version; nproc; grep MemTotal /proc/meminfo; } > "$OUT/host.txt"
CELLS=(unit-p65537-a1-s10-k10:5:unit-p65537-k10-D5 unit-p65537-a1-s10-k11:5:unit-p65537-k11-D5
  unit-p65537-a1-s10-k12:5:unit-p65537-k12-D5 unit-p65537-a1-s10-k13:5:unit-p65537-k13-D5
  unit-p163841-a3-s21-k13:5:unit-p163841-k13-D5 unit-p65537-a1-s10-k14:5:unit-p65537-k14-D5
  unit-p163841-a3-s21-k14:5:unit-p163841-k14-D5 unit-p65537-a1-s10-k15:5:unit-p65537-k15-D5
  unit-p65537-a1-s10-k10:6:unit-p65537-k10-D6 unit-p65537-a1-s10-k11:6:unit-p65537-k11-D6
  unit-p65537-a1-s10-k12:6:unit-p65537-k12-D6)
summary() { jq -s -r --arg n "$1" '"\($n): \(length) cells, \(map(select(.match)) | length) match; closedness re-check \(map(.verify_products) | add) products, \(map(.verify_nonzero) | add) non-zero; skipped \(map(.skipped_items // 0) | add) of \(map(.queued_items // 0) | add) items"' "$2"; }
# 1. closures with the closedness re-check
step 1 && EXTRA_FLAGS="--hky3-verify-closed on" "$GATE/run_gate.sh" "$BIN" "$OUT/closures" "${CELLS[@]}" > /dev/null
step 1 && summary closures "$OUT/closures/rows.jsonl" | tee "$OUT/SUMMARY.txt"
# 2. the skip configuration
step 2 && EXTRA_FLAGS="--hky3-skip on --hky3-verify-closed on" "$GATE/run_gate.sh" "$BIN" "$OUT/closures-skip" "${CELLS[@]}" > /dev/null
step 2 && summary closures-skip "$OUT/closures-skip/rows.jsonl" | tee -a "$OUT/SUMMARY.txt"
# 3. export and the independent checker
step 3 && mkdir -p "$OUT/basis-check"
step 3 && for D in 5 6; do
  b="$SCR/k10-D$D.basis"
  "$BIN" --instance "$PAPER/instances/unit-p65537-a1-s10-k10.json" --D $D --export-basis "$b" \
    --out "$OUT/basis-check/k10-D$D.json" > /dev/null 2> "$OUT/basis-check/k10-D$D.err"
  cp "$b.json" "$OUT/basis-check/k10-D$D.basis.sidecar.json"
  "$CHK" --basis "$b" --sidecar "$b.json" --closedness > "$OUT/basis-check/k10-D$D.check.json" 2> "$OUT/basis-check/k10-D$D.check.err"
  echo "basis-check k10 D$D: rc=$? $(jq -c '{ok, n_cols, rank, contains_one, profile, closedness}' "$OUT/basis-check/k10-D$D.check.json" 2>/dev/null | cut -c1-400)" | tee -a "$OUT/SUMMARY.txt"
  rm -f "$b" "$b.json"
done
# 4. checkpoint / resume
step 4 && "$GATE/ckpt_test.sh" "$BIN" "$PAPER/instances/unit-p65537-a1-s10-k10.json" 5 "$OUT/ckpt-k10D5" 3 2>&1 \
  | grep -v "Killed  " > "$OUT/ckpt-k10D5.log"
step 4 && grep -E "identical|byte-identical" "$OUT/ckpt-k10D5.log" | sed 's/^/ckpt k10 D5: /' | tee -a "$OUT/SUMMARY.txt"
# 5. negative controls at the discriminating cell (recorded: 1 in V_6, profile = all of M_<=6)
if step 5; then
  mkdir -p "$OUT/negctl"
  for ov in "--override-x 11932" "--override-x 4272" "--override-a 2"; do
    n=$(echo "$ov" | tr -d ' -')
    # shellcheck disable=SC2086
    "$BIN" --instance "$PAPER/instances/unit-p65537-a1-s10-k10.json" --D 6 $ov --out "$OUT/negctl/k10-D6-$n.json" \
      2> "$OUT/negctl/k10-D6-$n.err" > /dev/null
    got=$(jq -c '{profile, contains_one, closed, muladds}' "$OUT/negctl/k10-D6-$n.json")
    rec=$(grep -F "\"control\":\"$ov" "$GATE/negative_controls_d6.log" | jq -c '{profile, contains_one, closed, muladds}')
    same=false; [ "$got" = "$rec" ] && same=true
    echo "negctl k10 D6 $ov: $got; equal to the record (gate/negative_controls_d6.log): $same" | tee -a "$OUT/SUMMARY.txt"
  done
fi
echo EXTRACT_GATE_DONE | tee -a "$OUT/SUMMARY.txt"

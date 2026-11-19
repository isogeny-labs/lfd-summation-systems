#!/bin/bash
# Checkpoint/resume test: an uninterrupted closure against one killed with SIGKILL midway and resumed
# (and one stopped with SIGTERM, which checkpoints and exits cleanly, then resumed). Compares dims,
# profile, muladds, rows_reduced and the exported basis SHA-256.
# EXTRA_FLAGS: extra engine flags for every run (engine-skip round: "--hky3-skip on").
# Usage: ckpt_test.sh <pilot> <instance.json> <D> <out dir> <kill after N checkpoints>
set -u
BIN=$1; INST=$2; D=$3; OUT=$4; NK=${5:-3}
rm -rf "$OUT"; mkdir -p "$OUT"
summ() { jq -c '{columns, profile, contains_one, closed, muladds, rows_reduced, sha: .io.export.content_sha256, file_sha: .io.export.file_sha256, ckpts: .io.checkpoints, resumes: .io.resumes, resumed_at_batch: .io.resumed_at_batch}' "$1"; }
"$BIN" --instance "$INST" --D "$D" ${EXTRA_FLAGS:-} --export-basis "$OUT/ref.basis" --out "$OUT/ref.json" > /dev/null 2> "$OUT/ref.err"
echo "reference:  $(summ "$OUT/ref.json")"
for sig in KILL TERM; do
  dir="$OUT/ckpt-$sig"
  "$BIN" --instance "$INST" --D "$D" ${EXTRA_FLAGS:-} --checkpoint-dir "$dir" --checkpoint-every-batches 25 \
      --export-basis "$OUT/$sig.basis" --out "$OUT/$sig-first.json" > /dev/null 2> "$OUT/$sig-first.err" &
  pid=$!
  if [ "${NK#t}" != "$NK" ]; then
    # time trigger: kill after ${NK#t} seconds of wall
    end=$((SECONDS + ${NK#t})); while [ $SECONDS -lt $end ] && kill -0 $pid 2>/dev/null; do sleep 0.2; done
  else
    until [ "$(grep -c "hky3 checkpoint" "$OUT/$sig-first.err" 2>/dev/null)" -ge "$NK" ] || ! kill -0 $pid 2>/dev/null; do sleep 0.2; done
  fi
  kill -$sig $pid; wait $pid; rc=$?
  echo "$sig: first process rc=$rc after $(grep -c 'hky3 checkpoint' "$OUT/$sig-first.err") checkpoints; state: $(cat "$dir/state.json")"
  "$BIN" --instance "$INST" --D "$D" ${EXTRA_FLAGS:-} --resume "$dir" --checkpoint-every-batches 25 \
      --export-basis "$OUT/$sig.basis" --out "$OUT/$sig.json" > /dev/null 2> "$OUT/$sig.err"
  echo "$sig resumed: $(summ "$OUT/$sig.json")"
  same=$(jq -n --slurpfile a "$OUT/ref.json" --slurpfile b "$OUT/$sig.json" \
    '[$a[0], $b[0]] | map({columns, profile, contains_one, closed, muladds, rows_reduced, sha: .io.export.content_sha256}) | .[0] == .[1]')
  echo "$sig identical to reference: $same"
  cmp -s "$OUT/ref.basis" "$OUT/$sig.basis" && echo "$sig basis files byte-identical: true" || echo "$sig basis files byte-identical: false"
done
# disk is tight: keep the records, drop the bases and checkpoint states
ls -la "$OUT"/*.basis 2>/dev/null | awk '{print "basis bytes:", $5, $9}'
rm -f "$OUT"/*.basis "$OUT"/ckpt-*/state.ckpt "$OUT"/ckpt-*/state.ckpt.tmp

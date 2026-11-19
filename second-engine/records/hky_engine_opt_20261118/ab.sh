#!/bin/bash
# Interleaved A/B: ab.sh OUTDIR REPS "T C" "label:BIN:flags" ...   (RAYON_NUM_THREADS from the environment)
# One JSON line per run in OUTDIR/runs.jsonl (wall, cpu, peak RSS, load, counters).
set -u
OUT=$1; REPS=$2; CELLS=$3; shift 3
mkdir -p "$OUT/tmp"
for rep in $(seq 1 $REPS); do
  for cell in $CELLS; do
    T=${cell%:*}; C=${cell#*:}
    for arm in "$@"; do
      label=${arm%%:*}; rest=${arm#*:}; bin=${rest%%:*}; flags=${rest#*:}
      ./cell1.sh "$bin" $T $C "$OUT/tmp/${label}_t${T}c${C}_r$rep" $flags \
        | jq -c --arg label "$label" --argjson rep $rep --arg flags "$flags" '. + {label: $label, rep: $rep, flags: $flags}' >> "$OUT/runs.jsonl"
      tail -1 "$OUT/runs.jsonl" | jq -c '{label, rep, t, wall_s, cpu_s_approx, peak_rss_kb, load, m: .cl[0].muladds, dims: .cl[0].dims}'
    done
  done
done

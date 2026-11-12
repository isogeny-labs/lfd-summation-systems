#!/usr/bin/env bash
# Run scripts/ic_bounds_locus.py for a list of curves, one process per curve (each < 6 min, peak RSS ~0.2 GB).
# usage: [EXTEND=6] scripts/ic_bounds_locus_batch.sh PRIME K DMAX OUTDIR a,b [a,b ...]
set -euo pipefail
ulimit -v 1500000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
P=$1; K=$2; DMAX=$3; OUT=$4; shift 4
ROOT=$(cd "$(dirname "$0")/.." && pwd)
for c in "$@"; do
  f="$OUT/p${P}-k${K}-a${c%,*}-b${c#*,}.jsonl"
  rm -f "$f"
  timeout 1400 python3 "$ROOT/scripts/ic_bounds_locus.py" run --prime "$P" --k "$K" --curve "$c" --dmax "$DMAX" ${EXTEND:+--extend-to $EXTEND} --out "$f" > "$f.log" 2>&1 \
    && echo "done $c $(wc -l < "$f")" || echo "FAILED $c"
done

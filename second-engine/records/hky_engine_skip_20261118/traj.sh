#!/bin/bash
# Cumulative multiply-adds (products + second pass + interreductions) and batch wall at matched ranks,
# from --hky3-trace files: traj.sh STEP trace.jsonl ...   (one line per file and rank milestone)
set -u
STEP=$1; shift
for f in "$@"; do
  jq -s -r --arg f "$(basename "$f" .trace.jsonl)" --argjson step "$STEP" '
    reduce .[] as $e ({m: 0, s: 0, next: $step, out: []};
      if $e.ev == "batch" then
        .m += ($e.ops + $e.ops2) | .s += ($e.par_s + $e.sec_s)
        | if $e.rank >= .next then .out += [[.next, .m, .s]] | .next += $step else . end
      else .m += $e.ops | .s += $e.s end)
    | .out[] | "\($f) rank>=\(.[0]) muladds \(.[1]) engine_s \(.[2] | floor)"' "$f"
done

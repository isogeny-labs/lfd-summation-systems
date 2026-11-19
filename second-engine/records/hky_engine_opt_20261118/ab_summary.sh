#!/bin/bash
# ab_summary.sh RUNS.jsonl: per cell and arm, n, wall median/min/max, cpu median, peak RSS median (MB),
# and whether the deterministic counters (dims, muladds, rows_reduced, entries) agree across every run of the cell.
jq -s -r '
  def med: sort | if length == 0 then null else .[(length - 1) / 2 | floor] end;
  group_by(.t) | .[] | . as $cell
  | ($cell | map(.cl[0] | [.dims, .muladds, .rr, .entries]) | unique | length) as $nc
  | (group_by(.label) | .[] |
      [ "t=\(.[0].t) c=\(.[0].cl[0].c)", .[0].label, length,
        (map(.wall_s) | med | . * 10 | round / 10), (map(.wall_s) | min | . * 10 | round / 10), (map(.wall_s) | max | . * 10 | round / 10),
        (map(.cpu_s_approx) | med), (map(.peak_rss_kb) | med / 1024 | round),
        (map(.load[0] | tonumber) | med),
        (if $nc == 1 then "counters identical" else "COUNTERS DIFFER" end) ] | @tsv)
' "$1" | (printf 'cell\tarm\tn\twall_med\twall_min\twall_max\tcpu_med\trss_MB\tload1\tcounters\n'; cat)

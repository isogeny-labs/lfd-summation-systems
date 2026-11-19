#!/bin/bash
# Small-depth check of the new prime's recipe (p = 786433, a = 1, seed 10): 1 must enter V_D exactly at
# the expected last fall degree, floor((k+14)/4) (paper: 4 for k <= 5, 5 for k = 6..9, 6 for k = 10..13),
# i.e. 1 notin V_{L-1} and 1 in V_L. Every closure without 1 also exports its basis and runs the
# independent checker (exact reduction for k <= 9, annihilator test above).
BIN=$1; CHK=$2; OUT=$3; shift 3
mkdir -p "$OUT"
for k in "$@"; do
  L=$(( (k + 14) / 4 ))
  for D in $((L - 1)) $L; do
    n=unit-p786433-k$k-D$D
    "$BIN" --instance instances/unit-p786433-a1-s10-k$k.json --D $D --export-basis "$OUT/$n.basis" \
        --out "$OUT/$n.json" 2> "$OUT/$n.err" > /dev/null
    mode=--exact; [ $k -ge 10 ] && mode=""
    chk=$("$CHK" --basis "$OUT/$n.basis" --closedness $mode 2>/dev/null | jq -c '{ok, closedness, content_sha256}')
    jq -c --arg L $L --argjson chk "$chk" '{k: .instance.k, D, expected_lfd: ($L|tonumber), columns, profile, contains_one, closed, wall_s: .wall_s_this_process, check: $chk}' "$OUT/$n.json"
    rm -f "$OUT/$n.basis"
  done
done

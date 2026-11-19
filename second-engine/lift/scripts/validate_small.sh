#!/usr/bin/env bash
# Validation suite: lift vs direct closure (both computed by this crate; V_D goes through an
# hky.basis/v1 file, scrambled to a non-reduced echelon) plus the recorded oracle values of
# lfd-summation-systems experiments/ic-bounds-k14-independent/validation.json (dim, 1-membership,
# profile of the direct closure V_{D+1}).  Usage: validate_small.sh OUT.jsonl [SCRATCH]
set -euo pipefail
OUT=${1:?out.jsonl}
S=${2:-/tmp}
B=${HKY_LIFT:-$(dirname "$0")/../../../target/release/hky_lift}
: > "$OUT"
run() { # name p a b X k D expect_dim expect_one expect_profile [extra flags]
  local name=$1 p=$2 a=$3 b=$4 X=$5 k=$6 D=$7 ed=$8 eo=$9 ep=${10}; shift 10
  "$B" validate --p $p --a $a --b $b --X $X --k $k --D $D --dir "$S" \
     --expect-dim $ed --expect-one $eo --expect-profile $ep "$@" \
   | jq -c --arg name "$name" --arg flags "$*" '. + {name: $name, flags: $flags}' >> "$OUT" || true
  tail -1 "$OUT" | jq -c '{name,flags,all_equal,dim:.checks[0].lift,one:.checks[2].lift,s:.seconds.lift_incl_read}'
}
# name                   p      a    b    X    k D  dim   one profile
run unit-65537-k3-D3to4  65537  1    0    4271 3 3  57    1 1,7,22,42,57          --parts
run unit-65537-k4-D3to4  65537  1    0    4271 4 3  163   1 1,9,37,93,163         --parts
run unit-65537-k5-D3to4  65537  1    0    4271 5 3  386   1 1,11,56,176,386       --parts
run unit-65537-k6-D4to5  65537  1    0    4271 6 4  1586  1 1,13,79,299,794,1586  --parts
run unit-65537-k7-D4to5  65537  1    0    4271 7 4  3473  1 1,15,106,470,1471,3473 --parts
run unit-65537-k8-D4to5  65537  1    0    4271 8 4  6885  1 1,17,137,697,2517,6885 --parts
run unit-12289-k4-D3to4  12289  1    0    9362 4 3  163   1 1,9,37,93,163
run unit-12289-k6-D4to5  12289  1    0    9362 6 4  1586  1 1,13,79,299,794,1586
run unit-12289-k7-D4to5  12289  1    0    9362 7 4  3473  1 1,15,106,470,1471,3473
run unit-163841-k5-D3to4 163841 3    0   73730 5 3  386   1 1,11,56,176,386
run unit-163841-k6-D4to5 163841 3    0   73730 6 4  1586  1 1,13,79,299,794,1586
run unit-163841-k7-D4to5 163841 3    0   73730 7 4  3473  1 1,15,106,470,1471,3473
run unit-65537-k2-D2to3  65537  1    0    4271 2 2  5     0 0,0,1,5               --parts
run unit-65537-k3-D2to3  65537  1    0    4271 3 2  7     0 0,0,1,7
run unit-65537-k4-D2to3  65537  1    0    4271 4 2  9     0 0,0,1,9
run unit-65537-k6-D3to4  65537  1    0    4271 6 3  605   0 0,0,2,111,605         --parts
run unit-65537-k7-D3to4  65537  1    0    4271 7 3  829   0 0,0,2,115,829
run unit-65537-k8-D3to4  65537  1    0    4271 8 3  1061  0 0,0,2,119,1061
run unit-12289-k6-D3to4  12289  1    0    9362 6 3  605   0 0,0,2,111,605
run unit-163841-k7-D3to4 163841 3    0   73730 7 3  829   0 0,0,2,115,829
run proper-65537-X2-k6-D3to4  65537 4021 9917 2 6 3 605  0 0,0,2,111,605         --parts
run proper-65537-X2-k6-D4to5  65537 4021 9917 2 6 4 1584 0 0,11,77,297,792,1584  --parts
run proper-65537-X2-k7-D3to4  65537 4021 9917 2 7 3 829  0 0,0,2,115,829
run proper-65537-X2-k7-D4to5  65537 4021 9917 2 7 4 3471 0 0,13,104,468,1469,3471 --parts
run proper-65537-X2-k8-D3to4  65537 4021 9917 2 8 3 1061 0 0,0,2,119,1061
run proper-65537-X2-k8-D4to5  65537 4021 9917 2 8 4 6881 0 0,13,133,693,2513,6881 --parts
run proper-15361-k6-D4to5     15361 4021 9917 4491 6 4 1584 0 0,11,77,297,792,1584
run proper-15361-k7-D4to5     15361 4021 9917 7776 7 4 3469 0 0,11,102,466,1467,3469
run proper-15361-k7-D3to4     15361 4021 9917 7776 7 3 829  0 0,0,2,115,829
# variants: sigma-grouped column order 1 (what the engine writes with the split), include-low,
# falls only at the end, tiny batches, unscrambled input
run order1-unit-65537-k7-D4to5    65537 1 0 4271 7 4 3473 1 1,15,106,470,1471,3473 --order 1 --parts
run order1-proper-65537-X2-k7-D4to5 65537 4021 9917 2 7 4 3471 0 0,13,104,468,1469,3471 --order 1
run order1-unit-65537-k7-D3to4    65537 1 0 4271 7 3 829 0 0,0,2,115,829 --order 1
run low-unit-65537-k7-D4to5       65537 1 0 4271 7 4 3473 1 1,15,106,470,1471,3473 --include-low
run low-proper-65537-X2-k8-D4to5  65537 4021 9917 2 8 4 6881 0 0,13,133,693,2513,6881 --include-low
run endfalls-unit-65537-k8-D3to4  65537 1 0 4271 8 3 1061 0 0,0,2,119,1061 --fall-every 0
run endfalls-unit-65537-k7-D4to5  65537 1 0 4271 7 4 3473 1 1,15,106,470,1471,3473 --fall-every 0
run tiny-unit-65537-k7-D4to5      65537 1 0 4271 7 4 3473 1 1,15,106,470,1471,3473 --batch-tops 7 --cache-gb 0.0001
run tiny-proper-15361-k7-D4to5    15361 4021 9917 7776 7 4 3469 0 0,11,102,466,1467,3469 --batch-tops 5 --cache-gb 0.0001
run noscramble-unit-65537-k6-D4to5 65537 1 0 4271 6 4 1586 1 1,13,79,299,794,1586 --scramble 0
# larger: k = 9, D = 4 -> 5 (direct V_5(9) has 12,616 columns)
run unit-65537-k9-D4to5           65537 1 0 4271 9 4 12616 1 1,19,172,988,4048,12616
run proper-65537-X1424-k9-D4to5   65537 4021 9917 1424 9 4 12614 0 0,17,170,986,4046,12614
echo "cases: $(wc -l < "$OUT")  all_equal: $(jq -s 'map(select(.all_equal)) | length' "$OUT")"

#!/bin/sh
# Records of experiments/ic-rootfind/n131.jsonl and rust-scaling.jsonl.  Run from the repository root.
set -eu
cd "$(dirname "$0")"
cargo build --release --quiet
B=target/release/ic_rootfind
OUT=../../experiments/ic-rootfind
$B selftest
: > $OUT/n131.jsonl
for seed in 1 2 3; do
  $B bench --n 131 --k 44 --log2-range 22 --seed $seed --brute-log2 24 >> $OUT/n131.jsonl
  $B bench --n 131 --k 44 --log2-range 22 --seed $seed --columns >> $OUT/n131.jsonl
  $B bench --n 131 --k 44 --log2-range 22 --seed $seed --threads 4 >> $OUT/n131.jsonl
done
: > $OUT/rust-scaling.jsonl
for nk in "34 12" "46 16" "61 21" "97 33" "113 38" "131 44" "163 55" "191 64"; do
  set -- $nk
  [ "$2" -gt 62 ] && continue
  L=$2; [ "$L" -gt 21 ] && L=21
  $B bench --n $1 --k $2 --log2-range $L --seed 1 >> $OUT/rust-scaling.jsonl
done

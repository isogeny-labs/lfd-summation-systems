#!/usr/bin/env bash
# Negative controls on an exported V_D (unit k = 7, D = 4 -> 5, where 1 enters V_5; and proper k = 7):
# a corrupted basis must be rejected (trailer / closedness / generator check) or change the result;
# the test hooks (no fall shifts, a dropped shift) must change the result.  Usage: negative_controls.sh OUT.jsonl DIR
set -uo pipefail
OUT=${1:?}; S=${2:-/tmp}
B=${HKY_LIFT:-$(dirname "$0")/../../../target/release/hky_lift}
: > "$OUT"
one() { # label file [lift flags]
  local label=$1 f=$2; shift 2
  "$B" lift --basis "$f" "$@" > "$S/neg.json"; local rc=$?
  jq -c --arg label "$label" --argjson rc $rc '{label: $label, rc: $rc, status, reason, dim_VD: .dim_VD, dim_VD1, contains_one, profile, low_basis_sha256, input_checks}' "$S/neg.json" | tee -a "$OUT"
}
for inst in "unit 65537 1 0 4271" "proper 65537 4021 9917 2"; do
  set -- $inst; name=$1; p=$2; a=$3; b=$4; X=$5
  f="$S/neg-$name-k7-D4.hkybasis"
  "$B" closure --p $p --a $a --b $b --X $X --k 7 --D 4 --export "$f" --scramble 3 > /dev/null
  one "$name clean" "$f"
  one "$name clean, include-low" "$f" --include-low
  # trailer check: one byte of a coefficient flipped, trailer not updated
  "$B" corrupt --in "$f" --out "$f.flip" --mode flip-byte:2000 > /dev/null; one "$name flip-byte:2000 (raw)" "$f.flip"
  # consistent files (header, profile and trailer rewritten) with a changed space
  for i in 0 100 400 700 828; do
    "$B" corrupt --in "$f" --out "$f.drop$i" --mode drop-row:$i > /dev/null; one "$name drop-row:$i" "$f.drop$i"
    "$B" corrupt --in "$f" --out "$f.drop$i" --mode drop-row:$i > /dev/null; one "$name drop-row:$i --no-refuse" "$f.drop$i" --no-refuse
    "$B" corrupt --in "$f" --out "$f.pert$i" --mode perturb:$i > /dev/null; one "$name perturb:$i" "$f.pert$i"
    "$B" corrupt --in "$f" --out "$f.pert$i" --mode perturb:$i > /dev/null; one "$name perturb:$i --no-refuse" "$f.pert$i" --no-refuse
  done
  one "$name hook: no fall shifts" "$f" --test-no-fall-shifts
  one "$name hook: drop shift 0" "$f" --test-drop-shift 0
  one "$name hook: drop shift 13" "$f" --test-drop-shift 13
  rm -f "$f".*
done

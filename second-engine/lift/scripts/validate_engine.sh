#!/usr/bin/env bash
# Lift on the ENGINE's hky.basis/v1 exports (pkm_tower_pilot --instance --export-basis, column order 1),
# compared with this crate's direct closure V_{D+1} (dim, profile, 1, canonical hash of V_{D+1} n M_<=D
# in the file's column order) and the engine's V_D with our own closure V_D (same canonical hash).
# Usage: validate_engine.sh <pilot binary> OUT.jsonl DIR  "k D" ...
set -uo pipefail
PILOT=$1; OUT=$2; S=$3; shift 3
B=${HKY_LIFT:-$(dirname "$0")/../../../target/release/hky_lift}
I=${INST:-/home/user/crypto-wt/paper-engine/research/hky_paper_runs_20261118/instances}
: > "$OUT"
for kd in "$@"; do
  set -- $kd; k=$1; D=$2
  inst=$I/unit-p786433-a1-s10-k$k.json
  p=$(jq .prime $inst); a=$(jq .a $inst); b=$(jq .b $inst); X=$(jq '.target[0]' $inst)
  f=$S/eng-k$k-D$D.basis
  "$PILOT" --instance $inst --D $D --export-basis $f --out $f.run.json > /dev/null 2> $f.err
  "$B" validate --p $p --a $a --b $b --X $X --k $k --D $D --order 1 --dir $S --vd-basis $f --parts \
    | jq -c --arg k "$k" --arg D "$D" --arg pilot_sha "$(sha256sum "$PILOT" | cut -c1-64)" \
        '. + {name: ("engine-unit-p786433-k" + $k + "-D" + $D + "to" + (($D|tonumber)+1|tostring)), pilot_sha256: $pilot_sha}' >> "$OUT"
  tail -1 "$OUT" | jq -c '{name, all_equal, vd_equal: .vd_basis_vs_our_closure.equal, dim: .checks[0].lift, one: .checks[2].lift, s: .seconds.lift_incl_read}'
  rm -f $f
done

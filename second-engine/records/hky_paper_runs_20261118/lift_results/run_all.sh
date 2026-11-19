#!/usr/bin/env bash
# Every local validation of the Rust lift, with the frozen binary. Usage: run_all.sh <scratch dir>
set -uo pipefail
S=${1:?scratch}; mkdir -p "$S"
cd "$(dirname "$0")"
export HKY_LIFT=$PWD/bin/hky_lift
T=../../../tools/hky_lift/scripts
PILOT=${PILOT:-/tmp/claude-0/-home-user/daf172a2-9908-530c-9c36-29f5d702bf30/scratchpad/eng/pilot}
P=/home/user/lfd-summation-systems/experiments/ic-bounds-k14
sha256sum "$HKY_LIFT" "$PILOT" > binaries.sha256
$T/validate_small.sh validate_small.jsonl "$S" > validate_small.log 2>&1
$T/validate_engine.sh "$PILOT" validate_engine.jsonl "$S" "4 3" "5 3" "6 3" "6 4" "6 5" "7 3" "7 4" "7 5" "8 3" "8 4" "9 4" > validate_engine.log 2>&1
$T/negative_controls.sh negative_controls.jsonl "$S" > negative_controls.log 2>&1
# k = 10, D = 5 -> 6 on the engine's export, against the recorded lift records of the paper repository
for inst in unit-p65537-a1-s10-k10 proper-p65537-s12-k10; do
  f=$S/eng-$inst-D5.basis
  "$PILOT" --instance $P/instances/$inst.json --D 5 --export-basis $f --out $S/eng-$inst-D5.run.json > /dev/null 2> $S/eng-$inst-D5.err
  cp $f.json eng-$inst-D5.basis.sidecar.json
  jq '{columns, profile, contains_one, closed, muladds, wall_s: .wall_s}' $S/eng-$inst-D5.run.json > eng-$inst-D5.run.summary.json 2>/dev/null
  case $inst in
    unit*) exp=$P/runs/lift-unit-p65537-k10-D5to6-gpu.json ;;
    proper*) exp=$PWD/expected-proper-p65537-s12-k10-D5to6.json ;;
  esac
  "$HKY_LIFT" lift --basis $f --expect $exp --out lift-$inst-D5to6.json > /dev/null; echo "$inst rc=$?"
  rm -f $f
done
echo done

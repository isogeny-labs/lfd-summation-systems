#!/bin/bash
# Negative controls at the discriminating cell k=10, D=6 (recorded: 1 in V_6, profile = all of M_<=6).
# A perturbed target abscissa or a wrong curve coefficient must change the result.
BIN=$1; OUT=$2
I=/home/user/lfd-summation-systems/experiments/ic-bounds-k14/instances/unit-p65537-a1-s10-k10.json
for ov in "--override-x 4272" "--override-a 2"; do
  n=$(echo $ov | tr -d ' -')
  "$BIN" --instance $I --D 6 $ov --out $OUT/negctl-k10-D6-$n.json 2>$OUT/negctl-k10-D6-$n.err >/dev/null
  jq -c --arg ov "$ov" '{control:$ov, D, columns, profile, contains_one, closed, muladds, wall_s_this_process, gen_file_match: .system.generator_matches_file}' $OUT/negctl-k10-D6-$n.json
done

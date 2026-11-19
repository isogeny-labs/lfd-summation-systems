#!/usr/bin/env bash
# Markdown tables from the result files of run_all.sh
cd "$(dirname "$0")"
echo "| case | flags | dim V_D | direct V_{D+1}: dim, 1 in, profile | lift: dim, 1 in | profile, hash of V_{D+1} n M_<=D, sigma parts | recorded oracle (paper repo) | dim T / dim K / dim F / rounds | lift s / direct s | all equal |"
echo "|---|---|---|---|---|---|---|---|---|---|"
for f in validate_small.jsonl validate_engine.jsonl; do
jq -r --slurpfile V /home/user/lfd-summation-systems/experiments/ic-bounds-k14-independent/validation.json '
  def ck(n): (.checks[] | select(.field==n));
  def rec: (.name as $n | [$V[0].cases[] | select(.case==$n) | .oracle_reference_direct.dim_low_parts][0]);
  "| \(.name) | \(.flags // "" | gsub("--";"")) \(if .vd_basis_vs_our_closure then "engine export (order 1); engine V_D = our V_D: \(.vd_basis_vs_our_closure.equal)" else "" end) | \(.dim_VD) | \(ck("dim_VD1").direct), \(ck("contains_one").direct), \(ck("profile").direct|map(tostring)|join(" ")) | \(ck("dim_VD1").lift), \(ck("contains_one").lift) | \(ck("profile").equal), \(ck("low_basis_sha256").equal)\(if (.checks|map(.field)|index("sigma_parts")) then ", \(ck("sigma_parts").equal) \(ck("sigma_parts").lift.part_profiles_low|map(last)|map(tostring)|join("+"))" else "" end) | \(if (.checks|map(.field)|index("dim_VD1 vs recorded")) then "\(ck("dim_VD1 vs recorded").equal and ck("contains_one vs recorded").equal and ck("profile vs recorded").equal)" else "-" end)\(if rec then " (recorded parts \(rec|map(tostring)|join("+")))" else "" end) | \(.lift.dim_T) / \(.lift.dim_K) / \(.lift.falls) / \(.lift.rounds) | \(.seconds.lift_incl_read*100|round/100) / \(.seconds.direct_VD1*100|round/100) | **\(.all_equal)** |"' $f
done

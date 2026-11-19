#!/bin/bash
# Markdown gate table from rows.jsonl files (expected = lfd-summation-systems runs/*.json).
echo '| cell (instance : D) | columns | expected profile | got profile | 1 ∈ V_D exp / got | match | muladds | wall (s, shared 4 vCPU) | binary |'
echo '|:--|--:|:--|:--|:-:|:-:|--:|--:|:--|'
cat "$@" | jq -r '"| \(.cell | split(":")[0:2] | join(" : ")) | \(.got.columns) | \(.expected.profile | map(tostring) | join(",")) | \(.got.profile | map(tostring) | join(",")) | \(.expected.contains_one) / \(.got.contains_one) | **\(.match)** | \(.muladds) | \(.wall_s | floor) | \(.binary_sha256[0:8]) |"'

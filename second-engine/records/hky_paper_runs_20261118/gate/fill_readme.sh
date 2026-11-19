#!/bin/bash
# Regenerates README.md §3 between the GATE markers from the rows.jsonl files.
cd "$(dirname "$0")/.."
T=$(./gate/table.sh gate/raw/rows.jsonl gate/raw-frozen/rows.jsonl gate/raw-merged/rows.jsonl gate/raw-final/rows.jsonl)
BODY="<!-- GATE:begin -->
Expected = \`lfd-summation-systems/experiments/ic-bounds-k14/runs/<cell>.json\` (\`columns\`, \`profile\` = dim(V_D ∩ M_≤w) for
w = 0..D, \`contains_one\`); got = this engine's \`--instance\` record; \`match\` requires all three equal and the closure
closed. Binaries (§8): \`d31a012b\` = \`bin/pilot_gate\` (pre-merge engine); \`52b47153\` = first provisional build (same
sources, commit-stamped); \`0dfbffad\` = the merged engine (wt/engine-opt c51685b6a9 merged at 9fb60e9c13); **\`39acdeac\` =
the final pod binary** (9fb60e9c13 + non-fatal checkpoint/export failures, stamped 450d5b295b). Muladds are equal across all
binaries on every shared cell (the optimizer changes no counter). Scripts: \`gate/run_gate.sh\`, \`gate/final_gate.sh\`,
\`gate/table.sh\`; records under \`gate/raw*/\`. Walls: 4 shared vCPU; the pre-merge rows ran at load 15–27, the merged/final
rows at load ≈ 4.

$T

The profile definition is the same on both sides: the records' \`profile\` ends at \`rank\` and is cumulative by degree, as
the engine's \`dims\` (pivot degrees, cumulative). **Every row matches exactly; no discrepancy anywhere.**

Not reproduced locally: unit-p65537 and unit-p163841 k = 13, D = 6 (1 enters at full rank 313,912). One attempt with the merged
engine through \`pod/run_job.sh\` and a 14.5 GB memory limit (this box has 16 GB) reached 12.4 GB RSS after 37 min, the memory
guard sent SIGTERM, and the clean checkpoint then filled the disk (3.4 GB free) and the engine panicked; \`job.json\` recorded
\`censored_memory_guard_soft\` (\`gate/raw-merged/k13/\`). That exposed the disk-full panic, fixed in 450d5b295b (a failed
checkpoint or export is now reported and the closure continues). k = 13, D = 6 needs > 13 GB (the hky2 pod run peaked at
10.4 GB with sparser rows), so it belongs on a pod, as do the k = 14 D = 6 cells (J1/J2). An earlier pre-merge k = 12, D = 6
attempt died the same way (a 3.4 GB checkpoint at 10 min filled the disk); the merged and final binaries completed it without
checkpoints.
<!-- GATE:end -->"
awk -v b="$BODY" '/<!-- GATE:begin -->/{print b; skip=1; next} /<!-- GATE:end -->/{skip=0; next} !skip{print}' README.md > README.md.new
mv README.md.new README.md

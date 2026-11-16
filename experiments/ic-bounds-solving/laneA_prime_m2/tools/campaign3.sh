#!/bin/bash
# reduced part C (replaces the tail of campaign.sh, which was stopped during its k=12 p=12289 repeat)
cd "$(dirname "$0")/.."
ulimit -v 6000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
R="python3 tools/run_scan.py"
for k in 3 5 7 9; do $R $k 40961 20 200 0 none 0 11; $R $k 65537 20 200 0 none 0 11; done
$R 12 2147389441 20 60 0 none 0 11
echo "# end3"

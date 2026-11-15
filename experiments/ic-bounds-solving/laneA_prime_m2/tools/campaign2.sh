#!/bin/bash
# after campaign.sh: exception-frequency statistics, d_ff at k=7..10, d_lf quotient at k=7(,8), k=13
cd "$(dirname "$0")/.."
ulimit -v 6000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
R="python3 tools/run_scan.py"
echo "# start2 $(date -u +%FT%TZ)"
# exception frequency vs 1/p (brute force on all at k=3 to keep P1 honest)
for p in 73 97 137 193 257 401; do $R 3 $p 3000 20000 1 none 0 5; done
for p in 113 257 401 593; do $R 4 $p 1500 20000 0 none 0 5; done
for p in 3137 12289; do $R 5 $p 1000 3000 0 none 0 5; done
for k in 7 8 9 10; do python3 tools/dff_only.py $k 3; done
$R 7 $(python3 -c "
import sys; sys.path.insert(0,'tools')
from common import *
print(prime_above(3*4**7,7))") 2 20 0 quot 2 3
echo "# end2 $(date -u +%FT%TZ)"

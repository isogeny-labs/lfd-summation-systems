#!/bin/bash
# usage: run from laneB_prime_m3. Single thread, 6 GB virtual cap.
ulimit -v 6000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
rm -f data/scan_summaries.jsonl
for k in 2 3 4 5 6 7 8; do
  N=$((1<<k))
  PS=$(python3 - <<PY
import sys; sys.path.insert(0,'tools')
from build_m3 import primes_1mod
k=$k; N=1<<k
fac={2:0.25,3:0.25,4:0.25,5:0.25,6:0.5,7:4,8:4}[k]
print(primes_1mod(N,max(int(8**k*fac),12))[0])
PY
)
  for p in $PS 65537; do
    if [ $k -ge 7 ] && [ $p -eq 65537 ]; then NU=2; MT=30; else NU=20; MT=20000; fi
    python3 tools/scan.py $k $p $NU $MT 700
  done
  # locus
  python3 tools/scan.py $k 65537 6 5000 200 locus
done
echo ALLDONE

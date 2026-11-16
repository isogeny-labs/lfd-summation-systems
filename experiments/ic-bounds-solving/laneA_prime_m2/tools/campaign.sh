#!/bin/bash
# Sequential campaign (2 threads max).  Run from lane dir: bash tools/campaign.sh > data/campaign.log
cd "$(dirname "$0")/.."
ulimit -v 6000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
R="python3 tools/run_scan.py"
P() { python3 -c "
import sys; sys.path.insert(0,'tools')
from common import *
print($1)"; }
echo "# start $(date -u +%FT%TZ)"
# --- Part A: P1 brute force vs inverse, small k, many small primes (exception hunt) ---
for p in 41 73 89 97 113 137 193 257; do $R 3 $p 60 2000 1 none 0 1; done
$R 3 97 12 400 1 ord 12 2
for p in 97 113 193 241 257 337 401 1009 2017; do $R 4 $p 30 1000 1 none 0 1; done
$R 4 1009 12 400 1 ord 12 2
for p in 193 3137 12289 40961; do $R 5 $p 12 600 1 none 0 1; done
$R 5 3137 3 100 0 ord 3 2
$R 5 3137 12 100 0 quot 12 3
# k=6 brute force (one process per prime: 90 s enumeration each)
$R 6 12289 6 300 1 none 0 1
$R 6 12289 12 300 0 quot 12 2
# --- Part B: scan k=3..12 ---
for k in 3 4 5 6 7 8 9 10 11 12; do
  pm=$(P "prime_above(3*4**$k,$k)"); pb=$(P "prime_below(2**31-1,$k)"); p0=$(P "prime_above(2**$k+1,$k)"); ps=$(P "prime_above(2*4**$k,$k)")
  echo "## k=$k p_mid=$pm p_big=$pb p_smallest=$p0 p_2x4k=$ps"
  n=20; [ $k -le 6 ] && n=200
  $R $k $pm $n $((3*n)) 0 none 0 11
  [ $k -le 11 ] && $R $k $pb 20 60 0 none 0 11
  $R $k $p0 $n 300 0 none 0 11
  $R $k $ps $n $((4*n)) 0 none 0 11
  [ $k -le 12 ] && $R $k 12289 $n 300 0 none 0 11
done
for k in 3 5 7 9 12; do $R $k 40961 20 300 0 none 0 11; $R $k 65537 20 300 0 none 0 11; done
$R 12 $(P "prime_below(2**31-1,12)") 20 60 0 none 0 11
echo "# end $(date -u +%FT%TZ)"

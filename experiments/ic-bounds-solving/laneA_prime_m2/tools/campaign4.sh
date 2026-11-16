#!/bin/bash
cd "$(dirname "$0")/.."
ulimit -v 6000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 tools/task5_dfull.py 4 97 5
python3 tools/task5_dfull.py 5 193 5
python3 tools/task5_dfull.py 6 769 5
P=$(python3 -c "
import sys; sys.path.insert(0,'tools')
from common import *
print(prime_above(3*4**13,13))")
echo "k=13 p=$P"
python3 tools/run_scan.py 13 $P 5 12 0 none 0 11
echo "# end4"

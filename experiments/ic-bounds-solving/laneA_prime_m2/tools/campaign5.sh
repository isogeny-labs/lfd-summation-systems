#!/bin/bash
cd "$(dirname "$0")/.."
ulimit -v 6000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 tools/task5_dfull.py 3 41 5
python3 tools/task5_dfull.py 4 97 5
python3 tools/task5_dfull.py 5 193 5
python3 tools/task5_dfull.py 6 769 3
echo "# end5"

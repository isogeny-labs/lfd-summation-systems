#!/bin/bash
ulimit -v 6000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 tools/scan.py 7 8390273 6 5000 300 locus
python3 tools/scan.py 8 67109633 3 5000 300 locus
python3 tools/scan.py 2 13 300 400000 120
python3 tools/scan.py 2 17 300 400000 120
python3 tools/scan.py 3 17 300 400000 150
python3 tools/scan.py 3 41 300 400000 150
python3 tools/scan.py 4 97 300 400000 150
python3 tools/scan.py 4 113 300 400000 150
echo EXTRADONE

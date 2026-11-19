#!/bin/bash
cd "$(dirname "$0")/.."
until [ $(wc -l < raw/ab_main_thr4/runs.jsonl) -ge 27 ]; do sleep 20; done
./ab.sh raw/ab_dense 2 "10:5 14:5" "opt:bin/pilot_opt2:" "dense12:bin/pilot_opt2:--hky3-dense 1/2"
./ab.sh raw/ab_t11c6 1 "11:6" "base:bin/pilot_base:" "opt:bin/pilot_opt2:"
PLANTED=1 RANDOM_N=0 ./ab.sh raw/ab_planted 1 "10:6" "base:bin/pilot_base:" "opt:bin/pilot_opt2:"
echo QUEUE2 DONE

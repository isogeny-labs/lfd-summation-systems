#!/bin/bash
# Gate + checkpoint test for the final pod binary (bin/pkm_tower_pilot_final).
B=bin/pkm_tower_pilot_final
P=/home/user/lfd-summation-systems/experiments/ic-bounds-k14/instances
./gate/run_gate.sh $B gate/raw-final unit-p65537-a1-s10-k10:5:unit-p65537-k10-D5 unit-p65537-a1-s10-k11:5:unit-p65537-k11-D5 unit-p65537-a1-s10-k12:5:unit-p65537-k12-D5 unit-p65537-a1-s10-k10:6:unit-p65537-k10-D6 unit-p65537-a1-s10-k13:5:unit-p65537-k13-D5 unit-p163841-a3-s21-k13:5:unit-p163841-k13-D5 unit-p65537-a1-s10-k14:5:unit-p65537-k14-D5 unit-p163841-a3-s21-k14:5:unit-p163841-k14-D5 unit-p65537-a1-s10-k15:5:unit-p65537-k15-D5 unit-p65537-a1-s10-k11:6:unit-p65537-k11-D6 unit-p65537-a1-s10-k12:6:unit-p65537-k12-D6
./gate/ckpt_test.sh $B $P/unit-p65537-a1-s10-k10.json 5 gate/raw-final/ckpt-k10D5-t3 t3 2>&1 | grep -v "Killed  " > gate/ckpt-final-k10D5-t3.log
echo FINAL_GATE_DONE
